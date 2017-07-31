/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 *    - Improve algorithms
 */


// alibava includes ".h"
#include "AlibavaCorrelator.h"
#include "AlibavaBaseHistogramMaker.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"

// eutelescope includes ".h"
#include "EUTELESCOPE.h"

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
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerHitImpl.h>


// ROOT includes ".h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;

// Just a global 
static CellIDDecoder<TrackerHitImpl> hitDecoder(eutelescope::EUTELESCOPE::HITENCODING);

AlibavaCorrelator::AlibavaCorrelator () :
    AlibavaBaseHistogramMaker("AlibavaCorrelator"),
    // List of Histogram names, initialized here. As an example we put only 2
    _hHitPos("hHitPos"),
    _hCorX("hCorX"),
    _hCorY("hCorY"),
    _hSyncX("hSyncX"),
    _hSyncY("hSyncY"),
    _detectorIDs()
{
    // modify processor description
    _description = "AlibavaCorrelator gets hit collection and plots correlation histograms ";
    
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERHIT, "InputCollectionName",
            "Input raw data collection name",
            _inputCollectionName, string("rawdata") );
    
    // Details about HistoXMLFile
    registerProcessorParameter ("HistoXMLFile",
            "The path of XML file where the histograms are defined",
            _histoXMLFileName , string("AlibavaHistoList.xml"));

    registerProcessorParameter ("TopTagInXMLFile",
            "The top tag in HistoXMLFile",_topTag , string("AlibavaHistoList"));
    
    registerProcessorParameter ("TagToProcess", "The tag in TopTagInXMLFile. This"\
            " processor will only consider the histogram definitions inside this tag."\
            " This tag should be inside <TopTagInXMLFile> ... <TopTagInXMLFile/>",
            _tagToProcess , string("myAlibavaCorrelator"));
    
    registerProcessorParameter ("DetectorIDs",
            "The list of detector IDs", _detectorIDs , IntVec());
}


void AlibavaCorrelator::init () 
{
    streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    
    /* To choose if processor should skip masked events
     * ex. Set the value to 0 for false, to 1 for true
     */
    if(Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
    {
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "
            << ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" 
            << std::endl;
    }
	
    // here sort _detectorIDs
    std::sort(_detectorIDs.begin(),_detectorIDs.end());
	
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
	
    // if you want to change any variable defined in HistoXMLFile use this function
    // see example function below
    createRulesToChangeXMLValues();
}

void AlibavaCorrelator::createRulesToChangeXMLValues()
{
    // here we define which variables in HistoXMLFile we want to change or modify
    // You can only change or modify these variables:
    // xBin, xMin, xMax, yBin, yMin, yMax, title, labelX and labelY
    // others cannot be changed!
	
    /*! If you want to replace one of these variables: xBin, xMin, xMax, yBin, yMin, yMax
     * We will use the function
     * void changeXMLVariable(string histoName, string variableName, float newValue);
     *    defined in AlibavaBaseHistogramMaker
     * Note that variableNames are case sensitive!
     */
    // For example usually it is good idea to replace maximum event number for example
    // Lets say _someOtherHistoName histograms X axis should be total number of events
    //	changeXMLVariable(_someOtherHistoName, "xMax", float(_totalNumberOfEvents));
    // And _someOtherHistoName histograms Y axis should be 1000 since it is the max value it can get
    //changeXMLVariable(_someHistoName, "yMax", 1000.0);
    // It is really bad idea, but for some reason if you want
    // to hard code the binning of X axis, here how you can do.
    // Note that float number will be changed to integer in 
    // AlibavaBaseHistogramMaker::updateXMLVariables() method but here 
    // we should give it as float
    // changeXMLVariable(_someOtherHistoName, "xBin", float(100));
	
    /*! If you want to replace or modify one of these variables: title, labelX and labelY
     *  We will use the function
     *      		void addToXMLTitle(string histoName, string titleName, string whichSide, string stringToAdd);
     *          defined in AlibavaBaseHistogramMaker
     *          Note that titleName and whichSide are case sensitive!
     */
	
    // Say that you have the signal in X axis of _someHistoName.
    // Assuming that "Signal" is written in HistoXMLFile as this
    // histograms labelY, to get "(_multiplySignalby) Signal" you
    // should add signalMultipliedby string to the left
    //		addToXMLTitle(_someHistoName, "labelY", "left", signalMultipliedby);
    // Again it is usually bad idea but if you want to replace the X label here how you can do
    //		addToXMLTitle(_someOtherHistoName, "labelX", "all", string("Event Number"));
    // And here how you add string to right of title
    //		addToXMLTitle(_someHistoName, "title", "right", string("Some thing"));
}


void AlibavaCorrelator::processRunHeader (LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
	
    // Add processor name to the runheader
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
	
    // set total number of events in this run
    _totalNumberOfEvents = arunHeader->getNoOfEvents();
    streamlog_out ( DEBUG1 ) << "N events "<<_totalNumberOfEvents << endl;

    // reads and creates the histograms defined in HistoXMLFile
    bookHistos();
    
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}

std::string AlibavaCorrelator::getHistoNameForDetector(std::string name, int detID)
{
    return std::string(name+"_d"+std::to_string(detID));
}

std::string AlibavaCorrelator::getHistoNameForDetector(std::string name, int detID1, int detID2)
{
    return std::string(name+"_d"+std::to_string(detID1)+"_d"+std::to_string(detID2));
}

void AlibavaCorrelator::fillListOfHistos()
{
    addToHistoCheckList(_hHitPos);
    addToHistoCheckList(_hCorX);
    addToHistoCheckList(_hCorY);
    addToHistoCheckList(_hSyncX);
    addToHistoCheckList(_hSyncY);
    
    checkListOfHistosCreatedByXMLFile();
}

// Create TH2F histograms per detector
void AlibavaCorrelator::createTH2F(const std::string & histoName)
{
    AIDAProcessor::tree(this)->cd(this->name());
    AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
    TH2F * h = dynamic_cast<TH2F*> (_rootObjectMap[histoName]);
    
    streamlog_out ( MESSAGE1 )  << "hist "<< histoName<<" "<<h << std::endl;
    AIDAProcessor::tree(this)->mkdir(histoName.c_str());
    AIDAProcessor::tree(this)->cd(histoName.c_str());
    
    for(unsigned int idet=0; idet<_detectorIDs.size(); idet++) 
    {
        int detID = _detectorIDs[idet];
    	std::string newHistoName = getHistoNameForDetector(histoName, detID);
	TH2F * hnew = (TH2F*)h->Clone( newHistoName.c_str() );

        std::string title = hnew->GetTitle();
	title = title + std::string(" (det ") + to_string(detID) + std::string(")");
	_rootObjectMap.insert(std::make_pair(newHistoName, hnew));
    }
}

// XXX Change to createTH1F
void AlibavaCorrelator::createTH1F(const std::string & histoName)
{
    AIDAProcessor::tree(this)->cd(this->name());
    AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
    TH1F * h = dynamic_cast<TH1F*> (_rootObjectMap[histoName]);
    
    streamlog_out ( MESSAGE1 )  << "hist "<< histoName<<" "<<h << std::endl;
    AIDAProcessor::tree(this)->mkdir(histoName.c_str());
    AIDAProcessor::tree(this)->cd(histoName.c_str());
    
    for(unsigned int idet=0; idet<_detectorIDs.size(); idet++) 
    {
        int detID = _detectorIDs[idet];
    	std::string newHistoName = getHistoNameForDetector(histoName, detID);
	TH1F * hnew = (TH1F*)h->Clone( newHistoName.c_str() );

        std::string title = hnew->GetTitle();
	title = title + string(" (det ") + to_string(detID) + string(")");
        //hnew->SetTitle(title.c_str());

	_rootObjectMap.insert(std::make_pair(newHistoName, hnew));
    }
}

// XXX Change to createTH2F
void AlibavaCorrelator::createTH2F_comb(const std::string & histoName,bool isTitleSimple)
{
    AIDAProcessor::tree(this)->cd(this->name());
    AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
    TH2F * h = dynamic_cast<TH2F*> (_rootObjectMap[histoName]);
    
    AIDAProcessor::tree(this)->mkdir(histoName.c_str());
    AIDAProcessor::tree(this)->cd(histoName.c_str());
	
    // the _detectorIDs vector has to be sorted
    for (unsigned int idet=0; idet<_detectorIDs.size(); idet++) 
    {
        int detID = _detectorIDs[idet];
        for (unsigned int iCorDet=idet+1; iCorDet<_detectorIDs.size(); iCorDet++) 
        {
            int corDetID = _detectorIDs[iCorDet]; // correlated detector id
	    std::string newHistoName = getHistoNameForDetector(histoName, detID, corDetID);
	    TH2F * hnew = (TH2F*)h->Clone( newHistoName.c_str() );
	    std::string title = hnew->GetXaxis()->GetTitle();
            if(isTitleSimple)
            {
                title = title + string(" det") + to_string(detID);
            }
            else
            {
                title = title + string(" d") + to_string(detID) + string(" - d")+to_string(corDetID);
            }
            //hnew->SetTitle(title.c_str());
            
            _rootObjectMap.insert(std::make_pair(newHistoName, hnew));
        }
    }
}


void AlibavaCorrelator::bookHistos()
{
    // create histograms defined in HistoXMLFile
    processHistoXMLFile();

    // hit plot
    createTH2F(_hHitPos);
    // hCorX
    createTH2F_comb(_hCorX,true);
    // hCorY
    createTH2F_comb(_hCorY,true);
    // hSyncX
    createTH2F_comb(_hSyncX,false);
    // hSyncY
    createTH2F_comb(_hSyncY,false);
	
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << std::endl;
}

bool AlibavaCorrelator::isInDetectorIDsList(int detID)
{
    if(std::find(_detectorIDs.begin(),_detectorIDs.end(),detID) != std::end(_detectorIDs))
    {
        return true;
    }
    else
    {
        return false;
    }
}

void AlibavaCorrelator::processEvent (LCEvent * anEvent)
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    const int eventnum = alibavaEvent->getEventNumber();
    
    LCCollectionVec * collectionVec = nullptr;
    try
    {
        collectionVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
    } 
    catch(DataNotAvailableException& e  ) 
    {
        streamlog_out(MESSAGE2) <<  "No input collection " 
            << getInputCollectionName() << " found on event " 
            << alibavaEvent->getEventNumber()
            << " in run " << alibavaEvent->getRunNumber() << std::endl;
        return;
    }
    
    /////////////////////////////
    // Now loop ever detectors //
    const unsigned int noOfHits = collectionVec->getNumberOfElements();
    
    for(unsigned int ihit=0; ihit < noOfHits; ++ihit)
    {
        TrackerHitImpl * ahit = dynamic_cast<TrackerHitImpl*>(collectionVec->getElementAt(ihit));
	const int detID = hitDecoder(ahit)["sensorID"];
        if( !isInDetectorIDsList(detID) ) 
        {
            continue;
        }
        
        const double* pos = ahit->getPosition();
        // the hit position
        std::string histoName = getHistoNameForDetector(_hHitPos, detID);
        TH2F * h  = dynamic_cast<TH2F*>(_rootObjectMap[histoName]);
        h->Fill(pos[0],pos[1]);
        
        // correlation plots: (only check hits not previously checked)
	for(unsigned int i = ihit+1; i < noOfHits; ++i ) 
        {
            TrackerHitImpl * anotherHit = dynamic_cast< TrackerHitImpl * > ( collectionVec->getElementAt( i ) );
            int anotherDetID = hitDecoder( anotherHit )["sensorID"];
            // only consider hits from other detectors (higher than 
            // the current one) for correlation and synch
            if(detID >= anotherDetID)
            {
                continue;
            }
            const double* anotherPos = anotherHit->getPosition();
            histoName = getHistoNameForDetector(_hCorX, detID, anotherDetID);
            TH2F * hCorX = dynamic_cast<TH2F*> (_rootObjectMap[ histoName ]);
            hCorX->Fill(pos[0], anotherPos[0]);
            
            histoName = getHistoNameForDetector(_hCorY, detID, anotherDetID);
            TH2F * hCorY = dynamic_cast<TH2F*> (_rootObjectMap[ histoName ]);
            hCorY->Fill(pos[1], anotherPos[1]);
            
            histoName = getHistoNameForDetector(_hSyncX, detID, anotherDetID);
            TH2F * hSyncX = dynamic_cast<TH2F*> (_rootObjectMap[ histoName ]);
            hSyncX->Fill(eventnum, pos[0]-anotherPos[0]);
            
            histoName = getHistoNameForDetector(_hSyncY, detID, anotherDetID);
            TH2F * hSyncY = dynamic_cast<TH2F*> (_rootObjectMap[ histoName ]);
            hSyncY->Fill(eventnum, pos[1]-anotherPos[1]);
        }
    }
}


void AlibavaCorrelator::fillHistos(TrackerDataImpl * /* trkdata */ ){
	// nothing here
}

void AlibavaCorrelator::check (LCEvent * /* evt */ ) {
	// nothing to check here
}


void AlibavaCorrelator::end() {
	
	if (_numberOfSkippedEvents > 0)
		streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
	streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
	
}

// You have to define what bookEventHisto will do
void AlibavaCorrelator::bookEventHisto(int ){
	// does nothing
}

// You have to define what fillEventHisto will do
void AlibavaCorrelator::fillEventHisto(int , TrackerDataImpl * ){
	// does nothing
}
