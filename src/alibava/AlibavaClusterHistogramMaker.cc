/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 */


// alibava includes ".h"
#include "AlibavaClusterHistogramMaker.h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaCluster.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/tinyxml.h"

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
#include "TH1F.h"
#include "TH2F.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaClusterHistogramMaker::AlibavaClusterHistogramMaker():
    AlibavaBaseHistogramMaker("AlibavaClusterHistogramMaker"),
    // List of Histogram names, initialized here.
    _clusterSizeVsHitAmplitudeHistoName("hClusterSizeVsHitAmplitude"),
    _timeVsHitAmplitudeHistoName("hTimeVsHitAmplitude"),
    _timeVsSNRHistoName("hTimeVsSNR"),
    _SNRVsHitAmplitudeHistoName("hSNRVsHitAmplitude"),
    _etaVsSNRHistoName("hEtaVsSNR"),
    _etaFromClusterVsSNRHistoName("hEtaFromClusterVsSNR"),
    _timeVsSeedHistoName("hTimeVsSeed"),
    _signalVsSeedHistoName("hSignalVsSeed"),
    _etaVSTime("hEta_vs_Time"),
    _etaVSCoG("hEta_vs_CoG"),
    _etaVSUCPFA("hEta_vs_UCPFA"),
    _etaFromClusterVsTime("hEtaFromCluster_vs_Time"),
    _neighboursCharge("hNeighboursCharge"),
    _etaVSClusterSize("hEta_vs_ClusterSize"),
    _etaFromClusterVSClusterSize("hEtaFromCluster_vs_ClusterSize"),
    _clusterSizeVsCoG("hClusterSize_vs_CoG"),
    _clusterSizeVsUCPFA("hClusterSize_vs_UCPFA"),
    _timeVsClusters("hTimeVsClusters")
{
    // modify processor description
    _description ="AlibavaClusterHistogramMaker takes some "\
                   "type of Alibava data and produces histograms ";   
    
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
            "Input raw data collection name",_inputCollectionName, string("rawdata") );
    // Details about HistoXMLFile
    registerProcessorParameter("HistoXMLFile",
            "The path of XML file where the histograms are defined",
            _histoXMLFileName , string("AlibavaHistoList.xml"));
    
    registerProcessorParameter ("TopTagInXMLFile",
            "The top tag in HistoXMLFile",
            _topTag , string("AlibavaHistoList"));
    
    registerProcessorParameter ("TagToProcess",
            "The tag in TopTagInXMLFile. This processor will only consider "\
            "the histogram definitions inside this tag. This tag should be "\
            "inside <TopTagInXMLFile> ... <TopTagInXMLFile/>",
            _tagToProcess , string("myAlibavaClusterHistogramMaker"));
    
    // optional parameters, these parameters are defined in AlibavaBaseHistogramMaker, but they are optional
    registerOptionalParameter ("NoiseInputFile",
            "The filename where the pedestal and noise values stored",
            _pedestalFile , string(ALIBAVA::NOTSET));
    
    registerOptionalParameter ("CalibrationInputFile",
            "The filename where the calibration values stored",
            _calibrationFile , string(ALIBAVA::NOTSET));
    
    registerOptionalParameter ("NoiseCollectionName",
            "Noise collection name, better not to change",
            _noiseCollectionName, string (ALIBAVA::NOTSET));
    
    registerOptionalParameter ("ChargeCalibrationCollectionName",
            "Charge calibration collection name, better not to change",
            _chargeCalCollectionName, string (ALIBAVA::NOTSET));
    
    registerOptionalParameter ("PlotNoise",
            "Choose if noise should be plotted. If you want to plot "\
            "only noise set this variable true and set the noise collection name.",
            _plotPedestalAndNoise, bool(false));
    
    registerOptionalParameter ("EventsToPlot",
            "Choose if pedestal and noise should be plotted. If you want"\
            "to plot only noise or pedestal set this variable true and only"\
            "set the noise or pedeatal collection name you want to be plotted.",
            _plotEvents, IntVec());
    
    registerOptionalParameter ("MultiplySignalBy",
            "In case this variable is set, all signals will be multipled by this value.",
            _multiplySignalby, float(1.0));
    
    registerOptionalParameter ("PlotSomePercentOfEvents",
            "In case this variable is set (say x), x percent of total events "\
            "will be plotted randomly. The number should be between 0-100",
            _plotXPercentOfEvents, float(0.0));
}

void AlibavaClusterHistogramMaker::init () 
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;
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

    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
    // if you want to change any variable defined in HistoXMLFile use this function
    // see example function below
    this->createRulesToChangeXMLValues();
}

void AlibavaClusterHistogramMaker::createRulesToChangeXMLValues()
{
    // here we define which variables in HistoXMLFile we want to change or modify
    // You can only change or modify these variables:
    // xBin, xMin, xMax, yBin, yMin, yMax, title, labelX and labelY
    // others cannot be changed!
	
    this->changeXMLVariable(_maskedEventsHistoName, "xMax", float(_totalNumberOfEvents));
    
    // if MultiplySignalBy is set, it is a good idea to add it to the 
    // label of the "signal" axis
    if(_multiplySignalby != 1.0)
    {
        // For this first create the string you want add
        std::string signalMultipliedby(" ("+std::to_string(_multiplySignalby)+") ");
        addToXMLTitle(_clusterSizeVsHitAmplitudeHistoName, "labelY", "left", signalMultipliedby);
        addToXMLTitle(_timeVsHitAmplitudeHistoName, "labelY", "left", signalMultipliedby);
    }	
}


void AlibavaClusterHistogramMaker::processRunHeader (LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
    
    // Add processor name to the runheader
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
    
    // set total number of events in this run
    _totalNumberOfEvents = arunHeader->getNoOfEvents();
    streamlog_out ( DEBUG1 ) << "N events "<<_totalNumberOfEvents << endl;
    
    // get and set selected chips
    this->setChipSelection(arunHeader->getChipSelection());
    
    // set pedestal and noise values (if it is defined)
    this->setPedestals();
    
    // set channels to be used (if it is defined)
    // just after set pedestals and noise (if automasking
    // is enabled)
    this->setChannelsToBeUsed();
    
    if(_plotPedestalAndNoise)
    {
        plotPedestalAndNoise();
    }
    
    // set calibration if available
    this->setCalibration();
    
    // reads and creates the histograms defined in HistoXMLFile
    this->bookHistos();
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}

void AlibavaClusterHistogramMaker::fillListOfHistos()
{
    // Checks if all the histograms needed by this processor is defined in _histoXMLFileName
    // Unfortunately histogram names are hard coded, so we are checking if all histo names
    // exists in the _rootObjectMap
	
    // We need these histograms
	
    //////////////////////
    // One per each chip
    addToHistoCheckList_PerChip(_clusterSizeVsHitAmplitudeHistoName);
    addToHistoCheckList_PerChip(_timeVsHitAmplitudeHistoName);
    addToHistoCheckList_PerChip(_timeVsSNRHistoName);
    addToHistoCheckList_PerChip(_SNRVsHitAmplitudeHistoName);
    addToHistoCheckList_PerChip(_etaVsSNRHistoName);
    addToHistoCheckList_PerChip(_etaFromClusterVsSNRHistoName);
    addToHistoCheckList_PerChip(_etaFromClusterVsTime);
    addToHistoCheckList_PerChip(_neighboursCharge);
    addToHistoCheckList_PerChip(_timeVsSeedHistoName);
    addToHistoCheckList_PerChip(_signalVsSeedHistoName);
    addToHistoCheckList_PerChip(_etaVSTime);
    addToHistoCheckList_PerChip(_etaVSCoG);
    addToHistoCheckList_PerChip(_etaVSUCPFA);
    addToHistoCheckList_PerChip(_etaVSClusterSize);
    addToHistoCheckList_PerChip(_etaFromClusterVSClusterSize);
    addToHistoCheckList_PerChip(_clusterSizeVsCoG);
    addToHistoCheckList_PerChip(_clusterSizeVsUCPFA);
    addToHistoCheckList_PerChip(_calChargeHistoName);
    addToHistoCheckList_PerChip(_timeVsClusters);
    
    //////////////////////
    // One for all chips
    addToHistoCheckList(_maskedEventsHistoName);
    
    // here it checks
    checkListOfHistosCreatedByXMLFile();
}

void AlibavaClusterHistogramMaker::bookHistos()
{
    // create histograms defined in HistoXMLFile
    this->processHistoXMLFile();
	
    // If you set _plotEvents or _plotXPercentOfEvents
    // events
    if (_plotEvents.size() != 0 || _plotXPercentOfEvents != 0 )
    {
        // create the directory for event histograms
        AIDAProcessor::tree(this)->cd(this->name());
        AIDAProcessor::tree(this)->mkdir(_eventHistoDir.c_str());
    }
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << std::endl;
}



void AlibavaClusterHistogramMaker::processEvent (LCEvent * anEvent) 
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    int eventnum = alibavaEvent->getEventNumber();

    bool printClusters = false;
    /*if(eventnum == 7830) 
    {
        printClusters = true;
    }
    if(eventnum % 1000 == 0 )
    {
        streamlog_out ( DEBUG1 ) << "Looping events "<<eventnum << std::endl;
        printClusters = true;
    }*/
	
    // if _skipMaskedEvents is set
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        _numberOfSkippedEvents++;
        
        // Fill number of masked events histogram
        TH1F * histo = dynamic_cast<TH1F*> (_rootObjectMap[_maskedEventsHistoName]);
        histo->Fill(eventnum);
        
        return;
    }
	
    // if you set PlotSomePercentOfEvents or EventsToPlot paramaters
    bool plotThisEvent = false;
    /*if(isEventToBePlotted(eventnum))
    {
        // bookEventHisto will book a histogram for each chip
        this->bookEventHisto(eventnum);
        plotThisEvent = true;
    }*/
	
	
    /////////////////////////////
    // Now loop over detectors //
    LCCollectionVec * collectionVec=nullptr;
    try
    {
    	collectionVec = dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(getInputCollectionName()));
    }
    catch( lcio::DataNotAvailableException ) 
    {
        // do nothing again
        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
        return;
    }

    // For event-like quantities histograms
    std::map<int,int> clusterSize;
    for(auto ichip: getChipSelection())
    {
        clusterSize[ichip] = 0;
    }

    for(const auto & ptkdata : *collectionVec)
    {
        TrackerDataImpl * trkdata = dynamic_cast<TrackerDataImpl*>(ptkdata);
        AlibavaCluster anAlibavaCluster(trkdata);
    	
        // if this event selected to be plotted fill the histogram
        if(plotThisEvent)
        {
            // here we fill the event histogram for all clusters in this event
            // how it will fill is defined in fillEventHisto
    	    fillEventHisto(eventnum,trkdata);
        }
        // Event-like quantities

    	// for everything else
    	fillHistos(&anAlibavaCluster,alibavaEvent->getEventTime());
    	if(printClusters)
        {
            printAlibavaCluster(trkdata);
        }
        // What chip?
        int ichip = anAlibavaCluster.getChipNum();
    
        ++clusterSize[ichip];
    }
    // Histograms event-like quantities
    for(auto ichip: getChipSelection())
    {
        TH2F * histo = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_timeVsClusters,ichip)]);
        histo->Fill(alibavaEvent->getEventTime(), clusterSize[ichip]);
    }
}

void AlibavaClusterHistogramMaker::printAlibavaCluster(TrackerDataImpl* trkdata)
{
    AlibavaCluster anAlibavaCluster(trkdata);
    anAlibavaCluster.print();
}

void AlibavaClusterHistogramMaker::fillHistos(AlibavaCluster * anAlibavaCluster, float const & time )
{
    // We will first convert his trkdata to AlibavaCluster
    // see AlibavaCluster::AlibavaCluster(TrackerDataImpl * trkdata)
	
    //AlibavaCluster anAlibavaCluster(trkdata);
    int ichip = anAlibavaCluster->getChipNum();
    
    
    // Lets fill Cluster size histogram
    int clusterSize = anAlibavaCluster->getClusterSize();
    
    // And the hit amplitude depending on the cluster size (the last two could be just absorved 
    // by this one)
    TH2F *histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_clusterSizeVsHitAmplitudeHistoName,ichip)]);
    histo2->Fill( clusterSize,_multiplySignalby * anAlibavaCluster->getTotalSignal());

    // Hit amplitude depending on the TDC time
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_timeVsHitAmplitudeHistoName,ichip)]);
    histo2->Fill( time,_multiplySignalby * anAlibavaCluster->getTotalSignal());
    
    // Cluster signal to noise depending on the TDC time
    const float SNR_cluster = anAlibavaCluster->getTotalSNR(getNoiseOfChip(ichip));
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_timeVsSNRHistoName,ichip)]);
    histo2->Fill( time,SNR_cluster);
    
    // The SNR and the signal
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_SNRVsHitAmplitudeHistoName,ichip)]);
    histo2->Fill( SNR_cluster, _multiplySignalby*anAlibavaCluster->getTotalSignal() );
    
    // Cluster Seed channel depending on the TDC time
    const int seedChanNumber = anAlibavaCluster->getSeedChanNum();
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_timeVsSeedHistoName,ichip)]);
    histo2->Fill( time,seedChanNumber );
    // Cluster seed channel vs. the signal in the channel
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_signalVsSeedHistoName,ichip)]);
    histo2->Fill( _multiplySignalby*anAlibavaCluster->getTotalSignal(),seedChanNumber );

    // if calibration was defined
    if( _chargeCalMap.size() > 0)
    {
        histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_calChargeHistoName,ichip)]);
        histo2->Fill(time,_multiplySignalby*anAlibavaCluster->getTotalSignal(_chargeCalMap[ichip])*1e-3);
    }
    
    
    // Then fill eta histograms
    const float eta            = anAlibavaCluster->getEta();
    const float etaFromCluster = anAlibavaCluster->getEtaFromCluster();
    
    // the Eta (seed) vs. time
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_etaVSTime,ichip)]);
    histo2->Fill(time,eta);
    
    // the Eta (cluster, size 2 only) vs. time
    if(clusterSize == 2)
    {
        histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_etaFromClusterVsTime,ichip)]);
        histo2->Fill(time,etaFromCluster);
    }

    // Let's see only for clusters > 2 the distribution of the charge in its 3-first (if have) neighbours
    for(int kindex = 0; kindex < clusterSize; ++kindex)
    {
        // Obtain the channel number (normalized to the seed)
        const int channel = anAlibavaCluster->getChanNum(kindex)-seedChanNumber;
        // Not using the seed
        if(channel == 0)
        {
            continue;
        }
        histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_neighboursCharge,ichip)]);
        // Remeber seed is the first element of the cluster
        histo2->Fill(channel,anAlibavaCluster->getSignal(kindex)/anAlibavaCluster->getSignal(0));
    }

    
    // the Eta (seed) vs. SNR
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_etaVsSNRHistoName,ichip)]);
    histo2->Fill( eta,SNR_cluster );
    
    // the Eta (cluster) vs. SNR
    histo2 = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_etaFromClusterVsSNRHistoName,ichip)]);
    histo2->Fill( etaFromCluster, SNR_cluster );
    
    // center of gravity
    const float CoG = anAlibavaCluster->getCenterOfGravity();
    histo2 = dynamic_cast<TH2F*> (_rootObjectMap[getHistoNameForChip(_etaVSCoG,ichip)]);
    histo2->Fill(CoG, eta);
    
    // Using the Unbiased center position
    const float UCPFA = anAlibavaCluster->getUnbiasedCenterPosition();
    histo2 = dynamic_cast<TH2F*> (_rootObjectMap[getHistoNameForChip(_etaVSUCPFA,ichip)]);
    histo2->Fill(UCPFA, eta);
    
    histo2 = dynamic_cast<TH2F*> (_rootObjectMap[getHistoNameForChip(_clusterSizeVsUCPFA,ichip)]);
    histo2->Fill(UCPFA, clusterSize);
    
    // cluster size
    histo2 = dynamic_cast<TH2F*> (_rootObjectMap[getHistoNameForChip(_etaVSClusterSize,ichip)]);
    histo2->Fill(clusterSize, eta);
    
    histo2 = dynamic_cast<TH2F*> (_rootObjectMap[getHistoNameForChip(_etaFromClusterVSClusterSize,ichip)]);
    histo2->Fill(clusterSize, etaFromCluster);
    
    histo2 = dynamic_cast<TH2F*> (_rootObjectMap[getHistoNameForChip(_clusterSizeVsCoG,ichip)]);
    histo2->Fill(CoG,clusterSize);
}

void AlibavaClusterHistogramMaker::check (LCEvent * /* evt */ ) {
	// nothing to check here
}


void AlibavaClusterHistogramMaker::end() 
{
    if (_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents
            <<" events skipped since they are masked" << std::endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << std::endl;
}

// You have to define what bookEventHisto will do
void AlibavaClusterHistogramMaker::bookEventHisto(int eventnum)
{
    AIDAProcessor::tree(this)->cd(this->name());
    AIDAProcessor::tree(this)->cd(_eventHistoDir.c_str());
    
    EVENT::IntVec chipSelection = getChipSelection();
    
    for(unsigned int i=0; i<chipSelection.size(); i++) 
    {
        int ichip=chipSelection[i];
        
        std::string histoName = getEventHistoName(eventnum,ichip);
        // to be sure that we won't create more than one histogram per event per chip
	if (!doesRootObjectExists(histoName))
        {
            TH1F * eventHisto = new TH1F (histoName.c_str(),"",ALIBAVA::NOOFCHANNELS, -0.5, ALIBAVA::NOOFCHANNELS-0.5);
            
            std::string sp("Event "+std::to_string(eventnum)+" (chip "+std::to_string(ichip)+");Channel Number;");
            // if signal multiplied by some number put it to the label
            if(_multiplySignalby != 1)
            {
                sp += " ("+std::to_string( _multiplySignalby)+") ";
            }
            // in any case here is labelY
            sp+="Signal (ADCs)";
            eventHisto->SetTitle(sp.c_str());
            _rootObjectMap.insert(std::make_pair(histoName, eventHisto));
        }
    } // end of loop over selected chips
}

// You have to define what fillEventHisto will do
void AlibavaClusterHistogramMaker::fillEventHisto(int eventnum, TrackerDataImpl * trkdata)
{
    //noiseHisto->SetBinContent(ichan+1,noiVec[ichan]);
	
    //	streamlog_out (DEBUG1) << "Plotting Event "<< eventnum<<endl;
    AlibavaCluster anAlibavaCluster(trkdata);
    const int ichip = anAlibavaCluster.getChipNum();
	
    std::string histoName =getEventHistoName(eventnum,ichip);
    if (!doesRootObjectExists(histoName)) 
    {
        streamlog_out(ERROR3)<<"Root object for Event "<<eventnum<<" doesn't exists!"<<endl;
        return;
    }
    TH1F * histo = dynamic_cast<TH1F*> (_rootObjectMap[histoName]);
	
    // convert the trkdata to AlibavaCluster
    const int Nmembers = anAlibavaCluster.getClusterSize();
    for (int imember=0; imember<Nmembers; imember++) 
    {
        int channelNum = anAlibavaCluster.getChanNum(imember);
        float signal = _multiplySignalby * anAlibavaCluster.getSignal(imember);
        histo->SetBinContent(channelNum+1, signal);
    }
}
