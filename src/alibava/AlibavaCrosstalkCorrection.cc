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
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>


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
    _b1(0.017)
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
            _pedestalFile , string("pedestal.slcio"));
    registerProcessorParameter ("NoiseCollectionName",
            "Noise collection name, better not to change",
            _noiseCollectionName, std::string ("noise"));
    
    registerProcessorParameter ("b1",
            "The first neighbourg coefficient of the assymetric cross-talk",
            _b1, float (0.017));
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
    // Whether to activate or not the noisy channel auto masking
    streamlog_out(MESSAGE4) << "Noisy channel auto-masking is ";
    if(Global::parameters->isParameterSet(ALIBAVA::AUTOMASKINGACTIVE))
    {
        _isAutoMaskingActive = Global::parameters->getIntVal(ALIBAVA::AUTOMASKINGACTIVE);
    }
    streamlog_out(MESSAGE4) << _isAutoMaskingActive;

    if(Global::parameters->isParameterSet(ALIBAVA::AUTOMASKINGCRITERIA))
    {
        _autoMaskingCriterium = Global::parameters->getFloatVal(ALIBAVA::AUTOMASKINGCRITERIA);
    }
    streamlog_out(MESSAGE4) << " ( if ON, using " << _autoMaskingCriterium 
        << " sigmas )" << std::endl;
    
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
    
    // set pedestal and noise values
    this->setPedestals();
    
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
        // The first channel doesn't change
        newdatavec[0] = trkdata->getChargeValues()[0];
        for(int i=1; i < ALIBAVA::NOOFCHANNELS; ++i)
        {
            newdatavec[i] = trkdata->getChargeValues()[i]-_b1*trkdata->getChargeValues()[i-1];
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


