/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 *
 */


// alibava includes ".h"
#include "AlibavaPedestalSubtraction.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"


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


// ROOT includes ".h"

// system includes <>
#include <string>
#include <iostream>
#include <memory>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaPedestalSubtraction::AlibavaPedestalSubtraction():
    AlibavaBaseProcessor("AlibavaPedestalSubtraction")
{
    // modify processor description
    _description = "AlibavaPedestalSubtraction subtracts the "\
                    "provided pedestal values from the input raw data. ";
    // first of register the input collection (the ADC raw data)
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
            "Input raw data collection name",_inputCollectionName, string("rawdata") );
    
    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
            "Output data collection name",_outputCollectionName, string("recodata") );
    // if needed one can change these to optional parameters
    registerProcessorParameter ("PedestalInputFile",
            "The filename where the pedestal and noise values stored",
            _pedestalFile , string("pedestal.slcio"));
    // now the optional parameters
    registerProcessorParameter ("PedestalCollectionName",
            "Pedestal collection name, better not to change",
            _pedestalCollectionName, string ("pedestal"));
    registerProcessorParameter ("NoiseCollectionName",
            "Noise collection name, better not to change",
            _noiseCollectionName, string ("noise"));
}

void AlibavaPedestalSubtraction::init() 
{
    streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    /* To set of channels to be used
     * ex.The format should be like $ChipNumber:StartChannel-EndChannel$
     * ex. $0:5-20$ $0:30-100$ $1:50-70$
     * means from chip 0 channels between 5-20 and 30-100, from chip 1 
     * channels between 50-70 will be used (all numbers included). the
     * rest will be masked and not used. Note that the numbers should 
     * be in ascending order and there should be no space between two $ character
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
     * ex. Set the value to 0 for false, to 1 for true
     * */
    if(Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
    {
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "
            << ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << std::endl;
    }
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
}

void AlibavaPedestalSubtraction::processRunHeader(LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
    // Add processor name to the runheader
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
    
    // get and set selected chips
    setChipSelection(arunHeader->getChipSelection());
    
    // set pedestal and noise values
    setPedestals();
    
    // set channels to be used (if it is defined)
    // Note that it is needed to be placed after setPedestal
    // in order to use the automatic masking
    setChannelsToBeUsed();

    // if you want
    bookHistos();
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}


void AlibavaPedestalSubtraction::processEvent (LCEvent * anEvent) 
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        ++_numberOfSkippedEvents;
        return;
    }
	
    // retrieve the ADCs raw data
    LCCollectionVec * collectionVec=nullptr;
    try
    {
        collectionVec = dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection( getInputCollectionName()) );
    }
    catch( lcio::DataNotAvailableException ) 
    {
        // do nothing again
        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
    }
    
    // Preparing the corrected vector of ADCs data (a vector of channels per chip)
    LCCollectionVec* newDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    // The encoder to be used when storing the data
    CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,newDataCollection);
    
    // The loop for the correction
    const int noOfChip = collectionVec->getNumberOfElements();
    for( int i = 0; i < noOfChip; ++i )
    {
        // get the raw ADCs from the collection
	TrackerDataImpl * trkdata = dynamic_cast<TrackerDataImpl*>(collectionVec->getElementAt(i));
	const int chipnum = getChipNum(trkdata);
        FloatVec adc_vec = trkdata->getChargeValues();
        // Defining the corrected ADCs
        FloatVec corrected_adc_vec;
        corrected_adc_vec.clear();
        
        FloatVec pedVec = getPedestalOfChip(chipnum);
        // now subtract pedestal values from all channels
        for(size_t ichan=0; ichan< adc_vec.size();++ichan) 
        {
            if(isMasked(chipnum,ichan))
            {
                corrected_adc_vec.push_back(0);
                continue;
            }
            // subtract pedestal and store it
            const float newdata = adc_vec[ichan]-pedVec[ichan];
            corrected_adc_vec.push_back(newdata);
        }
        // Store the corrected ADCs in the collection
        TrackerDataImpl * newDataImpl = new TrackerDataImpl();
        newDataImpl->setChargeValues(corrected_adc_vec);
        chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
        chipIDEncoder.setCellID(newDataImpl);
        newDataCollection->push_back(newDataImpl);
    }
    alibavaEvent->addCollection(newDataCollection, getOutputCollectionName());
    // Note memory is freed by the framework
} 

void AlibavaPedestalSubtraction::check (LCEvent * /* evt */ ) 
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaPedestalSubtraction::end() 
{
    if(_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents 
            <<" events skipped since they are masked" << std::endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << std::endl;
}

void AlibavaPedestalSubtraction::fillHistos(TrackerDataImpl * /* trkdata */)
{
    // nothing is done here, if you want to plot histograms use AlibavaHistogramMaker processor
}

void AlibavaPedestalSubtraction::bookHistos()
{
    // nothing is done here, if you want to plot histograms use AlibavaHistogramMaker processor
}


