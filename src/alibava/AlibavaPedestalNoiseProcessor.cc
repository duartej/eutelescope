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
#include "AlibavaPedestalNoiseProcessor.h"
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

// ROOT includes ".h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"

// system includes <>
#include <string>
#include <iostream>
#include <sstream>
#include <memory>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaPedestalNoiseProcessor::AlibavaPedestalNoiseProcessor() :
    AlibavaBaseProcessor("AlibavaPedestalNoiseProcessor"),
    _pedestalHistoName ("hpedestal"),
    _noiseHistoName ("hnoise"),
    _temperatureHistoName("htemperature"),
    _chanDataHistoName ("Data_chan"),
    _chanDataFitName ("Fit_chan"),
    _adcmax(1000.0),
    _adcmin(0.0),
    _nbins(1000)
{
    // modify processor description
    _description ="AlibavaPedestalNoiseProcessor computes the pedestal"\
                   " and noise values of each channel";
	
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
            "Input raw data collection name", _inputCollectionName, string("rawdata") );
    registerProcessorParameter ("PedestalOutputFile",
            "The filename to store the pedestal and noise values",_pedestalFile , string("outputped.slcio"));
    // now the optional parameters
    registerOptionalParameter ("PedestalCollectionName",
            "Pedestal collection name, better not to change",_pedestalCollectionName, string ("pedestal"));
    registerOptionalParameter ("NoiseCollectionName",
            "Noise collection name, better not to change",_noiseCollectionName, string ("noise"));
    registerOptionalParameter ("MaxADCsCountsForHistograms",
            "The max ADCs counts which define the ranges of the histograms",
            _adcmax, float(1000.0) );
    registerOptionalParameter ("MinADCsCountsForHistograms",
            "The min ADCs counts which define the ranges of the histograms",
            _adcmax, float(0.0) );
    registerOptionalParameter ("NbinsForHistograms",
            "The number of bins for the histograms",
            _adcmin, float(1000) );

}


void AlibavaPedestalNoiseProcessor::init () 
{
    //streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    /* To set of channels to be used
     * ex.The format should be like $ChipNumber:StartChannel-EndChannel$
     * ex. $0:5-20$ $0:30-100$ $1:50-70$
     * means from chip 0 channels between 5-20 and 30-100, from chip 1 
     * channels between 50-70 will be used (all numbers included). the 
     * rest will be masked and not used. Note that the numbers should be 
     * in ascending order and there should be no space between two $ character
     */
    if(Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
    {
        Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << std::endl;
    }
    /* To choose if processor should skip masked events
     * ex. Set the value to 0 for false, to 1 for true
     */
    if(Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
    {
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
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

void AlibavaPedestalNoiseProcessor::processRunHeader(LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
    
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
    setChipSelection( arunHeader->getChipSelection() );
    // Masking channels (user input).  XXX:  Maybe decide a criteria to do it automaticaly?
    setChannelsToBeUsed();
    
    // An instance of the input output manager (Pedestal and noise)
    AlibavaPedNoiCalIOManager man;
    // Prepare the lcio file with the pedestal and noise, using the
    // run Header of the input lcio (i.e. the uncorrected pedestal lcio file?)
    man.createFile(_pedestalFile, arunHeader->lcRunHeader());
    // and book the needed histograms
    bookHistos();

    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}


void AlibavaPedestalNoiseProcessor::processEvent(LCEvent * anEvent) 
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        ++_numberOfSkippedEvents;
        return;
    }
    
    LCCollectionVec * collectionVec = nullptr;
    
    try
    {
        collectionVec = dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(getInputCollectionName()));
    } 
    catch( lcio::DataNotAvailableException ) 
    {
        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
    }
    
    // fill temperature histogram
    if(TH1F * temperatureHisto = dynamic_cast<TH1F*> (_rootObjectMap[_temperatureHistoName]))
    {
        temperatureHisto->Fill(alibavaEvent->getEventTemp());
    }
        
    // the number of the chips
    const unsigned int noOfDetector = collectionVec->getNumberOfElements();
    for ( size_t i = 0; i < noOfDetector; ++i )
    {
        TrackerDataImpl * trkdata = dynamic_cast<TrackerDataImpl*>( collectionVec->getElementAt(i) ) ;
        // fill the raw ADCs histos
        this->fillHistos(trkdata);
    }	
}

void AlibavaPedestalNoiseProcessor::check (LCEvent * /* evt */ ) 
{
    // nothing to check here - could be used to fill 
    // check plots in reconstruction processor
}


void AlibavaPedestalNoiseProcessor::end() 
{
    // Used the processEvent method to fill the ADCs histograms per channel/chip
    // and calculate the pedestal and noise per channel
    calculatePedestalNoise();
    if(_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents 
            <<" events skipped since they are masked" << std::endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << std::endl;
}

void AlibavaPedestalNoiseProcessor::calculatePedestalNoise()
{
    // Just defined here to avoid overload in the loops
    std::string tempHistoName,tempFitName;
    // Just to avoid the print of the canvas creation when
    // fitting (XXX: Another way?)
    TCanvas * cc = new TCanvas("cc","cc",600,800);
    
    EVENT::IntVec chipSelection = getChipSelection();
    // for each chip 
    for(auto ichip: chipSelection) 
    {
        // Get the pedestal and noise histos
        TH1F * hped = dynamic_cast<TH1F*> (_rootObjectMap[getPedestalHistoName(ichip)]);
        if(hped == nullptr)
        {
            streamlog_out ( ERROR5 ) << "Pedestal histogram '" 
                << getPedestalHistoName(ichip) << "' not booked." << std::endl;
        }
        TH1F * hnoi = dynamic_cast<TH1F*> (_rootObjectMap[getNoiseHistoName(ichip)]);
        if(hnoi == nullptr)
        {
            streamlog_out ( ERROR5 ) << "Pedestal histogram '" 
                << getNoiseHistoName(ichip) << "' not booked." << std::endl;
        }
        
        EVENT::FloatVec pedestalVec,noiseVec;
        // for each channel evaluate the mean pedestal using all the events
        for(int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ++ichan) 
        {
            double ped=0.0, noi=0.0;
            // if channel is masked, set pedestal and noise to 0,
            // otherwise, evaluate the
            if(!isMasked(ichip,ichan))
            {
                // First, get the histograms of the channel and the name
                // of the gaussian to fit
                tempFitName = getChanDataFitName(ichip, ichan);
                tempHistoName = getChanDataHistoName(ichip, ichan);
                TH1F * histo = dynamic_cast<TH1F*> (_rootObjectMap[tempHistoName]);
                TF1 * tempfit = dynamic_cast<TF1*> (_rootObjectMap[tempFitName]);
                // Fit to a gaussian
                histo->Fit(tempfit,"Q");
                ped = tempfit->GetParameter(1);
                noi = tempfit->GetParameter(2);
                hped->SetBinContent(ichan+1,ped);
                hnoi->SetBinContent(ichan+1,noi);
            }
            pedestalVec.push_back(ped);
            noiseVec.push_back(noi);
        }
        AlibavaPedNoiCalIOManager man;
        man.addToFile(_pedestalFile,_pedestalCollectionName, ichip, pedestalVec);
        man.addToFile(_pedestalFile,_noiseCollectionName, ichip, noiseVec);
    }
    delete cc;
}

// --- Some methods to keep consistency for names
std::string AlibavaPedestalNoiseProcessor::getChanDataHistoName(unsigned int ichip, unsigned int ichan)
{
    std::stringstream s;
    s<< _chanDataHistoName<<"_chip"<<ichip<<"_chan" << ichan;
    return s.str();
}

string AlibavaPedestalNoiseProcessor::getChanDataFitName(unsigned int ichip, unsigned int ichan)
{
    std::stringstream s;
    s<< _chanDataFitName<<"_chip"<<ichip<<"_chan"  << ichan;
    return s.str();
}

string AlibavaPedestalNoiseProcessor::getPedestalHistoName(unsigned int ichip)
{
    std::stringstream s;
    s<< _pedestalHistoName<<"_chip" << ichip;
    return s.str();
}

string AlibavaPedestalNoiseProcessor::getNoiseHistoName(unsigned int ichip)
{
    std::stringstream s;
    s<< _noiseHistoName<<"_chip" << ichip;
    return s.str();
}

// Fill the raw ADCs histograms
void AlibavaPedestalNoiseProcessor::fillHistos(TrackerDataImpl * trkdata)
{
    // Get the raw ADCs
    const lcio::FloatVec datavec = trkdata->getChargeValues();	
    int chipnum = getChipNum(trkdata);
    for(unsigned int ichan=0; ichan<datavec.size();++ichan) 
    {
        if(isMasked(chipnum, ichan))
        {
            continue;
        }
        
        std::string tempHistoName = getChanDataHistoName(chipnum, ichan);
        if( TH1F * histo = dynamic_cast<TH1F*> (_rootObjectMap[tempHistoName]) )
        {
            histo->Fill(datavec[ichan]);
        }
    }
}


void AlibavaPedestalNoiseProcessor::bookHistos()
{
    // Calculate the number of bins for the adc related
    // and set the ranges
    if(_adcmax < _adcmin)
    {
        streamlog_out( ERROR ) << "Inconsistent definition of ranges in the histograms: "
            << "[XXX TO FiNISH THE MESSAGE]"  << std::endl;
    }
    AIDAProcessor::tree(this)->cd(this->name());
    EVENT::IntVec chipSelection = getChipSelection();
    
    //the chipSelection should be in ascending order!
    //this is guaranteed with AlibavaConverter::checkIfChipSelectionIsValid()
    
    // temperature of event
    TH1F * temperatureHisto = new TH1F(_temperatureHistoName.c_str(),"Temperature",1000,-50,50);
    _rootObjectMap.insert(std::make_pair(_temperatureHistoName,temperatureHisto));
    
    // Create the histograms for the final pedestals/noise per channel
    for(unsigned int i=0; i<chipSelection.size(); ++i) 
    {
        unsigned int ichip=chipSelection[i];
        // Pedestal: 
        TH1F * pedestalHisto = new TH1F(getPedestalHistoName(ichip).c_str(),"",
                ALIBAVA::NOOFCHANNELS,-0.5,ALIBAVA::NOOFCHANNELS-0.5);
        _rootObjectMap.insert(std::make_pair(getPedestalHistoName(ichip), pedestalHisto)); 
        //title string for pedestal histogram
        const std::string sp("Pedestal (chip "+std::to_string(ichip)+");Channel Number;Pedestal (ADCs)");
        pedestalHisto->SetTitle(sp.c_str());
        
        // Noise: 
        TH1F * noiseHisto = new TH1F(getNoiseHistoName(ichip).c_str(),"",
                ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);
        _rootObjectMap.insert(std::make_pair(getNoiseHistoName(ichip),noiseHisto));
        //title string for noise histogram
        const std::string sn("Noise (chip "+std::to_string(ichip)+");Channel Number;Noise (ADCs)");
        noiseHisto->SetTitle(sn.c_str());
    }
    
    // create the directory to store the partial histograms/functions used
    // to evaluate the final pedestals and noise
    AIDAProcessor::tree(this)->mkdir(getInputCollectionName().c_str());
    AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
    
    // Histograms for the corrected pedestal and noise for each channel
    // (and chip)
    for(unsigned int i=0; i<chipSelection.size(); ++i) 
    {
        const unsigned int ichip=chipSelection[i];
        for( int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ++ichan) 
        {
            if(isMasked(ichip,ichan))
            {
                continue;
            }
            const std::string tempHistoName = getChanDataHistoName(ichip,ichan);
            const std::string tempFitName = getChanDataFitName(ichip,ichan);
            const std::string tempHistoTitle(tempHistoName+";ADCs;Entries");
            
            TH1F * chanDataHisto = new TH1F (tempHistoName.c_str(),"",_nbins,_adcmin,_adcmax);
            _rootObjectMap.insert(std::make_pair(tempHistoName, chanDataHisto));
            chanDataHisto->SetTitle(tempHistoTitle.c_str());
            
            // The gaussian fit associated: 
            //  + mean: corrected pedestal per channel
            //  + sigma: total noise per channel
            TF1 *chanDataFit = new TF1(tempFitName.c_str(),"gaus");
            _rootObjectMap.insert(std::make_pair(tempFitName, chanDataFit));
        }
    }
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms" << std::endl;
}


