/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
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
 */


// alibava includes ".h"
#include "AlibavaCommonModeSubtraction.h"
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
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include "TH1F.h"
#include "TH2F.h"

// system includes <>
#include <string>
#include <iostream>
#include <memory>
#include <algorithm>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaCommonModeSubtraction::AlibavaCommonModeSubtraction () :
    AlibavaBaseProcessor("AlibavaCommonModeSubtraction"),
    _commonmodeCollectionName(ALIBAVA::NOTSET),
    _commonmodeerrorCollectionName(ALIBAVA::NOTSET),
    _correctedSignalHistoName("totalcorrectedSignalHisto"),
    _adcmax(200),
    _chanDataHistoName ("Common_and_Pedestal_subtracted_data_channel")
{
    // modify processor description
    _description="AlibavaCommonModeSubtraction subtracts the provided"\
                  " common mode values from the input reco (pedestal subtracted) data. ";   
    
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
            "Input reco data collection name",_inputCollectionName, string("recodata") );

    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
            "Output data collection name",_outputCollectionName, string("recodata_cmmd") );   
    
    // if needed one can change these to optional parameters    
    // now the optional parameters
    registerProcessorParameter("CommonModeCollectionName",
            "Common mode collection name, (see AlibavaConstantCommonModeProcessor)",
            _commonmodeCollectionName, string ("commonmode"));
    
    registerProcessorParameter ("CommonModeErrorCollectionName",
            "Common mode error collection name, (see AlibavaConstantCommonModeProcessor",
            _commonmodeerrorCollectionName, string ("commonmodeerror"));
    
    registerOptionalParameter ("MaxADCsCountsForHistograms",
            "The max ADCs counts for the common mode which define the ranges of the histograms",
            _adcmax, float(200.0) );
}


void AlibavaCommonModeSubtraction::init () 
{
    streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    /* To set of channels to be used 
     ex.The format should be like $ChipNumber:StartChannel-EndChannel$ 
     ex. $0:5-20$ $0:30-100$ $1:50-70$ 
     means from chip 0 channels between 5-20 and 30-100, from chip 1 
     channels between 50-70 will be used (all numbers included). the rest 
     will be masked and not used. Note that the numbers should be in ascending
     order and there should be no space between two $ character
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
        _skipMaskedEvents = bool(Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS));
    }
    else
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS 
            <<" is not set! Masked events will be used!" << std::endl;
    }
    
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
}

void AlibavaCommonModeSubtraction::processRunHeader (LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    // Add processor name to the runheader
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
    
    
    setChipSelection(arunHeader->getChipSelection());
    // set channels to be used (if it is defined)
    setChannelsToBeUsed();
    
    // if you want
    this->bookHistos(arunHeader->getNoOfEvents());
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}


void AlibavaCommonModeSubtraction::processEvent (LCEvent * anEvent) 
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
	
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        _numberOfSkippedEvents++;
        return;
    }
    
    // Get the needed inputs, ADCs (pedestal substracted data), and common mode shift
    LCCollectionVec * dataColVec = nullptr;
    LCCollectionVec * cmmdColVec = nullptr;
    try
    {
        dataColVec = dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(getInputCollectionName()));
    } 
    catch(lcio::DataNotAvailableException )
    {
        // do nothing again
    	streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
    }
    try
    {
        cmmdColVec =dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(_commonmodeCollectionName));
    }
    catch(lcio::DataNotAvailableException)
    {
        // do nothing again
    	streamlog_out( ERROR5 ) << "Collection ("<<_commonmodeCollectionName<<") not found! " << std::endl;
    }
    // Sanity check: these collections have same number of elements
    if( (dataColVec->getNumberOfElements()) != (cmmdColVec->getNumberOfElements()) ) 
    {
        streamlog_out( ERROR5 ) << "Number of elements in collections are not equal!" <<std::endl;
        streamlog_out( ERROR5 ) << getInputCollectionName() << " has " 
            << dataColVec->getNumberOfElements() << " elements while "<< _commonmodeCollectionName <<" has "
            << cmmdColVec->getNumberOfElements() << std::endl;
    }
    
    // The ouput data pedestal and common mode corrected
    LCCollectionVec * newColVec = new LCCollectionVec(LCIO::TRACKERDATA);    
    CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,newColVec);
    
    // Start the processing (per each bettle)
    const unsigned int noOfChips=dataColVec->getNumberOfElements();    	
    for(unsigned int i = 0; i < noOfChips; ++i )
    {
        // get data from the collection: ped-substrated, and common mode
        TrackerDataImpl *dataImpl = dynamic_cast<TrackerDataImpl*>(dataColVec->getElementAt(i));
        TrackerDataImpl *cmmdImpl = dynamic_cast<TrackerDataImpl*>(cmmdColVec->getElementAt(i));
    	TrackerDataImpl * newdataImpl = new TrackerDataImpl();
        // check that they belong to same chip
    	if( (getChipNum(dataImpl)) != (getChipNum(cmmdImpl)) ) 
        {
            streamlog_out( ERROR5 ) << "The chip numbers in the collections is not same! " << std::endl;
        }
        const int chipnum = getChipNum(dataImpl);
    
    	// set chip number for newdataImpl
        chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM]=chipnum;
        chipIDEncoder.setCellID(newdataImpl);
        
        const FloatVec datavec= dataImpl->getChargeValues();
    	const FloatVec cmmdvec= cmmdImpl->getChargeValues();
    	
        // More sanity checks: the size of data sets are equal to ALIBAVA::NOOFCHANNELS
    	if( int(datavec.size()) != ALIBAVA::NOOFCHANNELS )
        {
            streamlog_out( ERROR5 ) << "Number of channels in input data is" 
                << " not equal to ALIBAVA::NOOFCHANNELS! "<< std::endl;
        }
        if( int(cmmdvec.size()) != ALIBAVA::NOOFCHANNELS )
        {
            streamlog_out( ERROR5 ) << "Number of channels in common mode data"
                << " is not equal to ALIBAVA::NOOFCHANNELS! " << std::endl;
        }
    
    		
    	// now subtract common mode values from all channels
        EVENT::FloatVec newdatavec(datavec.size());
        std::transform(datavec.begin(),datavec.end(),cmmdvec.begin(),newdatavec.begin(),
                [](const float & dataPedSub,const float & cmmd) { return dataPedSub-cmmd; } );
        // Just put to zero those masked channels (presumely that is already done
        // in every step, pedestal substraction, common mode, etc... but again)
    	for(size_t ichan=0; ichan<datavec.size();ichan++) 
        {
            if(isMasked(chipnum, ichan)) 
            {
                newdatavec[ichan] = 0.0;
            }
        }
        // fill histograms
    	this->fillHistos(newdatavec,anEvent->getEventNumber(),chipnum);
        
        // And add it to the lcio class
    	newdataImpl->setChargeValues(newdatavec);
    	newColVec->push_back(newdataImpl);
    }
    alibavaEvent->addCollection(newColVec, getOutputCollectionName());
}

void AlibavaCommonModeSubtraction::check (LCEvent * /* evt */ ) 
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaCommonModeSubtraction::end() 
{
    if(_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << std::endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << std::endl;
	
}

void AlibavaCommonModeSubtraction::fillHistos(const EVENT::FloatVec & datavec,const int & evt, const int & chipnum)
{
    // Fill the histograms with the corrected data
    for(unsigned int ichan = 0; ichan < datavec.size(); ++ichan)
    {
        if(this->isMasked(chipnum, ichan))
        {
            continue;
        }

        std::string tempHistoName = this->getChanDataHistoName(chipnum, ichan);
        if( TH1F * histo = dynamic_cast<TH1F*>(_rootObjectMap[tempHistoName]) )
        {
            histo->Fill(datavec[ichan]);
        }
        std::string tempHistoName1 = this->getHistoNameCS(chipnum);
        if( TH2F * histo1 = dynamic_cast<TH2F*>(_rootObjectMap[tempHistoName1]) )
        {
            histo1->Fill(evt,datavec[ichan]);
        }
    }
}

std::string AlibavaCommonModeSubtraction::getChanDataHistoName(int chipnum, int ichan)
{
    return std::string(_chanDataHistoName+"_chip_"+
            std::to_string(chipnum)+"_chan_"+std::to_string(ichan));
}

std::string AlibavaCommonModeSubtraction::getHistoNameCS(const int & ichip)
{
    return std::string(_correctedSignalHistoName+"_"+std::to_string(ichip));
}

void AlibavaCommonModeSubtraction::bookHistos(const int & _NtotalEvt)
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
    for(auto ichip: chipVec)
    {
        // be careful not overloading the memory
        unsigned int nbins=NtotalEvt;
        if(NtotalEvt > 20000)
        {
            nbins=20000;
        }
        // a histogram showing the corrected signal over events
        std::string tempHistoTitle2("Corrected Signal per Event;Event Number;ADC corrected [counts]");
        TH2F * histo2 = new TH2F(this->getHistoNameCS(ichip).c_str(),tempHistoTitle2.c_str(),
                nbins,0,NtotalEvt,nbins_adc,b_adcmin,b_adcmax);
        _rootObjectMap.insert(make_pair(this->getHistoNameCS(ichip), histo2));

        // shows the common mode shift error (std-dev of the common mode shift mean)
        /*std::string tempHistoTitle3("Common Mode Shift Error (std.dev) ;Event Number;#sigma_{Common Mode Shift} [ADC counts]");
        TH2F * histo3 = new TH2F(this->getHistoNameCMSError(ichip).c_str(),tempHistoTitle3.c_str(),
                nbins,0,NtotalEvt,100,_adcmin_err,_adcmax_err);
        _rootObjectMap.insert(make_pair(this->getHistoNameCMSError(ichip), histo3));*/
    }
    // The histograms of corrected signal per channel 
    AIDAProcessor::tree(this)->cd(this->name());
    for(const auto & chipnum: chipVec)
    {
        for(int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ichan++) 
        {
            if(isMasked(chipnum, ichan)) 
            {
                continue;
            }
            std::string tempHistoName = this->getChanDataHistoName(chipnum,ichan);
            std::string tempHistoTitle(tempHistoName+";ADCs;NumberofEntries");
            TH1F * chanDataHisto=new TH1F(tempHistoName.c_str(),tempHistoTitle.c_str(),
                    nbins_adc,b_adcmin,b_adcmax);
            _rootObjectMap.insert(std::make_pair(tempHistoName, chanDataHisto));
        }
    }   
    streamlog_out(MESSAGE1)<< "End of Booking histograms." << std::endl; 
}

