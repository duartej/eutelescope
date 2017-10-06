/*
 *  Created by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *
 */

// alibava includes ".h"
#include "AlibavaCrosstalkCorrection.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaCluster.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <EVENT/LCFloatVec.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>

// ROOT includes
#include "TFile.h"
#include "TGraphErrors.h"


// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <algorithm>
#include <functional>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaCrosstalkCorrection::AlibavaCrosstalkCorrection () :
    AlibavaBaseProcessor("AlibavaCrosstalkCorrection"),
    _xtInitialized(false),
    _f(),
    _nmax(2),
    _coefficients_name("")
{
    // modify processor description
    _description ="AlibavaCrosstalkCorrection corrects the effect of an assymetric "\
                   " (inductive) cross-talk between channels";
    // first of register the input /output collection
    registerInputCollection(LCIO::TRACKERDATA, "InputCollectionName",
            "Input collection name, it should be pedestal subtracted",
            _inputCollectionName, std::string("recodata") );
    registerOutputCollection(LCIO::TRACKERDATA, "OutputCollectionName",
            "Output data collection name", _outputCollectionName, std::string("alibava_clusters_corrected") );
	
    // if needed one can change these to optional parameters
    registerProcessorParameter ("NoiseInputFile",
            "The filename where the pedestal and noise values stored",
            _pedestalFile , std::string("pedestal.slcio"));
    registerProcessorParameter ("NoiseCollectionName",
            "Noise collection name, better not to change",
            _noiseCollectionName, std::string ("noise"));
    
    registerProcessorParameter ("CollectionName",
            "Name of the collection coefficients, if is not found in the event, the"\
            " content of the 'f' property is used",
            _coefficients_name, std::string("xtfactors"));
    
    registerProcessorParameter ("f",
            "The neighbourg coefficients of the assymetric cross-talk, "\
            "ordered by neighbour distance to the seed",
            _f, std::vector<float>({0.006}));
   
    //int nmax_prov(0);
    registerProcessorParameter ("MaximumNeighbourgs",
            "Forces the number of maximum neighbourgs to be considered in the correction",
            _nmax, int(2));
}


void AlibavaCrosstalkCorrection::init () 
{
    streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    // this method is called only once even when the rewind is active
    /* To set of channels to be used
         ex.The format should be like $ChipNumber:StartChannel-EndChannel$
         ex. $0:5-20$ $0:30-100$ $1:50-70$ means from chip 0 channels
         between 5-20 and 30-100, from chip 1 channels between 50-70 will 
         be used (all numbers included). the rest will be masked and not used
         Note that the numbers should be in ascending order and there should
         be no space between two $ character
         */
    if(Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
    {
        Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "
            << ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << std::endl;
    }
    
    /* To choose if processor should skip masked events
     ex. Set the value to 0 for false, to 1 for true
     */
    if(Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
    {
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "
            << ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << std::endl;
    }
    
    // usually a good idea to
    printParameters ();	
}

void AlibavaCrosstalkCorrection::processRunHeader (LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
    
    // Add processor name to the runheader
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
    
    // get and set selected chips
    this->setChipSelection(arunHeader->getChipSelection());

    // set channels to be used (if it is defined)
    this->setChannelsToBeUsed();	
	
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}


void AlibavaCrosstalkCorrection::processEvent (LCEvent * anEvent) 
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        _numberOfSkippedEvents++;
        return;
    }

    
    // As an input collection we get pedestal subtracted alibava data
    // This collection contains AlibavaEvents
    LCCollectionVec * inputColVec=nullptr;
    try
    {
        inputColVec=dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(getInputCollectionName()));
    }
    catch( lcio::DataNotAvailableException ) 
    {
        // do nothing again
        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
        return;
    }
    	
    // Cross-check
    if( inputColVec->getNumberOfElements() != static_cast<int>(this->getChipSelection().size()) )
    {
        streamlog_out(ERROR) << "Inconsistency in the Input Collection '" << 
            getInputCollectionName() << "'. Must contain the same number of elements than " 
            << " chips evaluated (currently " << this->getChipSelection().size() << std::endl;
        return;
    }
    
    // The cross-talk factors, if were included in the DB, 
    // therefore obtain them, if were not done before
    if(!_xtInitialized)
    {
        LCCollectionVec * xtColVec=nullptr;
        try
        {
            xtColVec=dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(_coefficients_name));
            // Get the vector of coefficients
            if(xtColVec->getNumberOfElements() != 1)
            {
                streamlog_out(ERROR) << "Unexpected format for the cross-talk " 
                    << "coeficients '" << _coefficients_name << "'" << std::endl;
                throw lcio::DataNotAvailableException("");
            }
            const EVENT::LCFloatVec * xtfactors = dynamic_cast<EVENT::LCFloatVec*>(xtColVec->getElementAt(0));
            _f.clear();
            for(unsigned int i=0; i < xtfactors->size() && static_cast<int>(i) < _nmax ; ++i)
            {
                // Remember the cross talk factors are obtained directly
                // from the difference bewteen left-right, and we are
                // applying the filter by subtracting to the right channel 
                // (to compensate back)
                // Check that they are positive?? XXX: IT seems to behave 
                // better (convergence of the factorrs) if we don't do 
                // that check (but needs a minimum of N=5
                /*if( (*xtfactors)[i] > 0)
                {
                    // Assuming this is the last acceptable 
                    // left to right cross-talk
                    break;
                }*/
                _f.push_back((-1.0)*(*xtfactors)[i]);
            }

        }
        catch(lcio::DataNotAvailableException) 
        {
            // There is no factors introduced by DB-file, let's
            // use the property
            
            // First time let's if sizes are coherent
            if(static_cast<int>(_f.size()) < this->_nmax)
            {
                streamlog_out(ERROR) << "Requested coefficients: " << this->_nmax 
                    << " BUT only provided " << _f.size() << std::endl;
                return;
            }
        }
        // Initialize and dump message first time
        streamlog_out(MESSAGE) << "Neighbour coefficients for correction: " <<std::endl;
        for(unsigned int i = 0; i < _f.size(); ++i)
        {
            streamlog_out(MESSAGE) << "  [" << i << "-coefficient]:" <<_f[i] << std::endl;
        }
        _xtInitialized=true;
    }

    // The ouput data pedestal and common mode corrected
    LCCollectionVec * newColVec = new LCCollectionVec(LCIO::TRACKERDATA);    
    CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,newColVec);

    for(int i = 0; i < inputColVec->getNumberOfElements(); ++i)
    {
        // Note that the i-index is referring to the input col element (which is related with 
        // the chip number but is not exactly that. The collection elements must follow the same 
        // chip order stored in the getChipSelection
        const int chipnum = this->getChipSelection()[i];
        
        // The new LCIO TrackerData corrected
    	TrackerDataImpl * newdataImpl = new TrackerDataImpl();
        // set chip number for newdataImpl
        chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM]=chipnum;
        chipIDEncoder.setCellID(newdataImpl);
        
        // get the data from the collection 
        TrackerDataImpl * trkdata = dynamic_cast<TrackerDataImpl*>(inputColVec->getElementAt(i));

        EVENT::FloatVec newdatavec(trkdata->getChargeValues());
        // FIR application: y[n] = x[n]+f1*x[n-1]*+f2*x[n-2]+....
        for(unsigned int  n= 0; n < ALIBAVA::NOOFCHANNELS; ++n)
        {
            newdatavec[n] = trkdata->getChargeValues()[n];
            // Note j=Neighbourg order, j-1: Elements filter factor
            for(unsigned int j =1; j <= _f.size() &&  j-1 < n ; ++j)
            {
                // Not using masked channels to correct 
                if(this->isMasked(chipnum,n-j))
                {
                    continue;
                }
                newdatavec[n] -= _f[j-1]*trkdata->getChargeValues()[n-j];
            }
        }
        // Add it to the LCIO class
        newdataImpl->setChargeValues(newdatavec);
        newColVec->push_back(newdataImpl);

    } // end of loop ever detectors
    alibavaEvent->addCollection(newColVec, getOutputCollectionName());
}

void AlibavaCrosstalkCorrection::check (LCEvent * /* evt */ ) 
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaCrosstalkCorrection::end() 
{
    if (_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents
            <<" events skipped since they are masked" << std::endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << std::endl;	
}


