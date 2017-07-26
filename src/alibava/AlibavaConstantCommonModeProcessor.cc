/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 *
 *  modified by: Eda Yildirim eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 *    - Change algorithm to evaluate common noise, including a 
 *      convergence loop to avoid signal contamination in the 
 *      calculation
 *    - Major changes in the classes logic, remove unnecessary 
 *      methods.
 *    - Remove 1d histo, and add 2d for error
 *
 */

// alibava includes ".h"
#include "AlibavaConstantCommonModeProcessor.h"
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
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes ".h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"

// system includes <>
#include <string>
#include <iostream>
#include <sstream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;

// TODO::put masking option

AlibavaConstantCommonModeProcessor::AlibavaConstantCommonModeProcessor () :
    AlibavaBaseProcessor("AlibavaConstantCommonModeProcessor"),
    _commonmodeCollectionName(ALIBAVA::NOTSET),
    _commonmodeerrorCollectionName(ALIBAVA::NOTSET),
    _NoiseDeviation(2.5),
    _commonModeShiftHistoName("commonModeShiftOverEvents"),
    _commonModeShiftErrorHistoName ("commonModeShiftErrorOverEvents"),
    _adcmax(200),
    _adcmax_err(12),
    _adcmin_err(0),
    _commonmode(),
    _commonmodeerror()
{
    // modify processor description
    _description ="AlibavaConstantCommonModeProcessor computes the common "\
                   "mode values of each chip and their errors";
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
            "Input data collection name (should be pedestal subtracted!)",
            _inputCollectionName, string("recodata") );
    // now the optional parameters
    registerOptionalParameter("CommonModeCollectionName",
            "Common mode collection name, better not to change",
            _commonmodeCollectionName, string ("commonmode"));
    
    registerOptionalParameter("CommonModeErrorCollectionName",
            "Common mode error collection name, better not to change",
            _commonmodeerrorCollectionName, string ("commonmodeerror"));
    
    registerOptionalParameter ("NoiseDeviation",
            "The limit to the deviation of noise. The data exceeds this deviation "\
            "will be considered as signal and not be included in common mode error calculation",
            _NoiseDeviation, float(2.5) );
    
    registerOptionalParameter ("MaxADCsCountsForHistograms",
            "The max ADCs counts for the common mode which define the ranges of the histograms",
            _adcmax, float(200.0) );
    
    registerOptionalParameter("MaxADCsCountsErrorForHistograms",
            "The max ADCs counts for the common mode error which define the ranges of the histograms",
            _adcmax_err, float(10.0) );
    
    registerOptionalParameter("MinADCsCountsErrorForHistograms",
            "The min ADCs counts for the common mode error which define the ranges of the histograms",
            _adcmin_err, float(0.0) );
}

void AlibavaConstantCommonModeProcessor::init()
{
    streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    /* To set of channels to be used
     * ex.The format should be like $ChipNumber:StartChannel-EndChannel$
     * ex. $0:5-20$ $0:30-100$ $1:50-70$ means from chip 0 channels between
     * 5-20 and 30-100, from chip 1 channels between 50-70 will be used 
     * (all numbers included). the rest will be masked and not used. Note 
     * that the numbers should be in ascending order and there should be no 
     * space between two $ character
     */
    if(Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
    {
        Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED 
            <<" is not set! All channels will be used!" << std::endl;
    }
    
    /* To choose if processor should skip masked events
     * ex. Set the value to 0 for false, to 1 for true
     */
    if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
    {
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS 
            <<" is not set! Masked events will be used!" << endl;
    }
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
}

void AlibavaConstantCommonModeProcessor::processRunHeader (LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
    
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
    
    setChipSelection( arunHeader->getChipSelection() );
    setChannelsToBeUsed();
    
    this->bookHistos(arunHeader->getNoOfEvents());
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}

void AlibavaConstantCommonModeProcessor::processEvent(LCEvent * anEvent)
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        _numberOfSkippedEvents++;
        return;
    }
    
    // Get the collection of pedestal corrected ADCs signal data
    LCCollectionVec * collectionVec = nullptr;
    try
    {
        collectionVec = dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(getInputCollectionName()));
    }
    catch(lcio::DataNotAvailableException)
    {
        // do nothing again
    	streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
    }

    // The output collections: the Common Noise  and its error
    // And the encoders (when store the new collections)
    // Common noise
    LCCollectionVec* commonCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    CellIDEncoder<TrackerDataImpl> commonCol_CellIDEncode (ALIBAVA::ALIBAVADATA_ENCODE,commonCollection);
    // common noise error
    LCCollectionVec* commerrCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    CellIDEncoder<TrackerDataImpl> commerrCol_CellIDEncode (ALIBAVA::ALIBAVADATA_ENCODE,commerrCollection);

    const unsigned int noOfDetector =collectionVec->getNumberOfElements();
    for(unsigned int i = 0; i < noOfDetector; ++i)    
    {
        // Get the data
        TrackerDataImpl * trkdata = dynamic_cast<TrackerDataImpl*>(collectionVec->getElementAt(i));
        const int chipnum = getChipNum(trkdata);
        // Find Commonmode and its error (i.e. the common noise)
        auto commonModePair = this->calculateConstantCommonMode(trkdata->getChargeValues(),getChipNum(trkdata));
        // And store it
    	TrackerDataImpl * commonData = new TrackerDataImpl();
    	commonData->setChargeValues(commonModePair.first);
    	commonCol_CellIDEncode[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
        commonCol_CellIDEncode.setCellID(commonData);
        commonCollection->push_back(commonData);
        
        TrackerDataImpl * commerrData = new TrackerDataImpl();
    	commerrData->setChargeValues(commonModePair.second);
    	commerrCol_CellIDEncode[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
    	commerrCol_CellIDEncode.setCellID(commerrData);
    	commerrCollection->push_back(commerrData);
    	
    	// Fill Histos
    	this->fillHistos(commonModePair,anEvent->getEventNumber(),chipnum);
    }
    alibavaEvent->addCollection(commonCollection, getCommonModeCollectionName());
    alibavaEvent->addCollection(commerrCollection, getCommonModeErrorCollectionName());
}

void AlibavaConstantCommonModeProcessor::check (LCEvent * /* evt */ ) 
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}

void AlibavaConstantCommonModeProcessor::end() 
{
    if (_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << std::endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

// The workhorse of this class
std::pair<EVENT::FloatVec,EVENT::FloatVec> AlibavaConstantCommonModeProcessor::calculateConstantCommonMode(const EVENT::FloatVec & datavec, const int & chipnum)
{
    streamlog_out( DEBUG ) << "Calculate Constant common mode:" << std::endl;
    streamlog_out( DEBUG ) << " ++ Chip " << chipnum << " of " << getNumberOfChips() << ", now iterating...";

    // Just a cross-check
    if(datavec.size() == 0)
    {
        streamlog_out(ERROR) << " Number of not masked and valid ADCs" 
            << " (pedestal substracted) channels is ZERO." << std::endl;
    }
    // Map to keep only non-signal ADCs 
    // ichannel: ADC
    std::map<int,float> non_signal_map;
    // initialize the map
    for(int i=0; i < static_cast<int>(datavec.size()); ++i)
    {
        // Non including in the calculation the 
        // masked channels if any
        if(this->isMasked(chipnum,i))
        {
            continue;
        }
        non_signal_map.emplace(i,datavec[i]);
    }

    // 1. Obtain the mean and the standard deviation for the signal ADCs 
    //    (pedestal substracted) (i.e. the common noise)
    float mean_signal  = this->getMean(non_signal_map);
    float stddev_signal= this->getStdDev(non_signal_map,mean_signal);

    // 2. Use the mean and standard deviation to exclude real
    //    signal presence. Criteria, anything above 2.5 sigmas from the mean
    //    is considered signal
    unsigned int last_vec_size = non_signal_map.size();
    do
    {
        // after the check from the 'while' statement, update
        // the last know vector size
        last_vec_size = non_signal_map.size();
        for(auto it = non_signal_map.cbegin(); it != non_signal_map.cend(); /* no increment*/)
        {
            if(isMasked(chipnum,it->first))
            {
                ++it;
                continue;
            }
            // Remove the signal channels 
            if( std::abs(it->second-mean_signal)/stddev_signal > _NoiseDeviation )
            {
                non_signal_map.erase(it++);
                // [XXX: SHOULD I MASKED THIS CHANNEL?]
            }
            else
            {
                // or keep it
                ++it;
            }
        }
        // 3. Recalculate the mean with the excluded signal channels
        mean_signal  = this->getMean(non_signal_map);
        stddev_signal= this->getStdDev(non_signal_map,mean_signal);
    } while( non_signal_map.size() != last_vec_size );

    EVENT::FloatVec commonmodeVec(ALIBAVA::NOOFCHANNELS);
    EVENT::FloatVec commonmodeerrorVec(ALIBAVA::NOOFCHANNELS);
    // Now reloop over all channels and create the output vector,
    // the same value for all the channels
    for(int ichan = 0; ichan < ALIBAVA::NOOFCHANNELS ; ++ichan)
    {
    	commonmodeVec[ichan] = mean_signal;
    	commonmodeerrorVec[ichan] = stddev_signal;
    }
    streamlog_out( DEBUG0 ) << " Noise common mode: " << mean_signal << " Error:"
        << stddev_signal << std::endl;

    return std::make_pair(commonmodeVec,commonmodeerrorVec);
}


void AlibavaConstantCommonModeProcessor::fillHistos(const std::pair<EVENT::FloatVec,EVENT::FloatVec> & commonNoise, const int & nevt, const int & chipnum)
{
    if(commonNoise.first.size() == 0)
    {
        return;
    }
    // Fill the histograms with the corrected data
    // Remember just to use the first, all the values
    // are exactly the same
    if(TH2F * h = dynamic_cast<TH2F*>(_rootObjectMap[this->getHistoNameCMS(chipnum)]) )
    {
        h->Fill(nevt,commonNoise.first[0]);
    }
    if(TH2F * h = dynamic_cast<TH2F*>(_rootObjectMap[this->getHistoNameCMSError(chipnum)]) )
    {
        h->Fill(nevt,commonNoise.second[0]);
    }
}


// getter and setter for _commonmodeCollectionName
void AlibavaConstantCommonModeProcessor::setCommonModeCollectionName(std::string commonmodeCollectionName){
	_commonmodeCollectionName = commonmodeCollectionName;
}
std::string AlibavaConstantCommonModeProcessor::getCommonModeCollectionName(){
	return _commonmodeCollectionName;
}

// getter and setter for _commonmodeerrorCollectionName
void AlibavaConstantCommonModeProcessor::setCommonModeErrorCollectionName(std::string commonmodeerrorCollectionName){
	_commonmodeerrorCollectionName = commonmodeerrorCollectionName;
}
std::string AlibavaConstantCommonModeProcessor::getCommonModeErrorCollectionName(){
	return _commonmodeerrorCollectionName;
}


std::string AlibavaConstantCommonModeProcessor::getHistoNameCMS(const int & ichip)
{
    return _commonModeShiftHistoName+"_"+std::to_string(ichip);
}

std::string AlibavaConstantCommonModeProcessor::getHistoNameCMSError(const int & ichip)
{
    return _commonModeShiftErrorHistoName+"_"+std::to_string(ichip);
}

void AlibavaConstantCommonModeProcessor::bookHistos(const int & _NtotalEvt)
{
    // Check the total number of events, and put the whole number
    // whenever the relevant global is not present
    int NtotalEvt=Global::parameters->getIntVal("MaxRecordNumber");
    if(NtotalEvt == 0)
    {
        NtotalEvt=_NtotalEvt;
    }
    // Calculate the number of bins for the adc related
    // and set the ranges
    int b_adcmax = _adcmax;
    int b_adcmin = -1.0*_adcmax;
    if(_adcmax < 0)
    {
        b_adcmax=b_adcmin;
        b_adcmin=-1.0*b_adcmax;
    }
    int nbins_adc=2.0*_adcmax;

    EVENT::IntVec chipVec = this->getChipSelection();
    for(const auto & ichip: chipVec)
    {
        // be careful not overloading the memory
        unsigned int nbins=NtotalEvt;
        if(NtotalEvt > 20000)
        {
            nbins=20000;
        }
        // a histogram showing the common mode shift (mean over all channels per event) over events
        std::string tempHistoTitle2("Common Mode Shift over Events;Event Number;Common Mode Shift [ADC counts]");
        TH2F * histo2 = new TH2F(this->getHistoNameCMS(ichip).c_str(),tempHistoTitle2.c_str(),nbins,0,NtotalEvt,nbins_adc,b_adcmin,b_adcmax);
        _rootObjectMap.insert(make_pair(this->getHistoNameCMS(ichip), histo2));

        // shows the common mode shift error (std-dev of the common mode shift mean)
        std::string tempHistoTitle3("Common Mode Shift Error (std.dev) ;Event Number;#sigma_{Common Mode Shift} [ADC counts]");
        TH2F * histo3 = new TH2F(this->getHistoNameCMSError(ichip).c_str(),tempHistoTitle3.c_str(),
                nbins,0,NtotalEvt,100,_adcmin_err,_adcmax_err);
        _rootObjectMap.insert(make_pair(this->getHistoNameCMSError(ichip), histo3));
    }
    
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl; 
}
