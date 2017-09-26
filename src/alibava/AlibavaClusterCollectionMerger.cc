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



// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaClusterCollectionMerger.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"

// marlin includes
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"


// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

// system includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <memory>
#include <stdlib.h>
#include <algorithm>

using namespace std;
using namespace marlin;
using namespace lcio;
using namespace alibava;
using namespace eutelescope;

AlibavaClusterCollectionMerger::AlibavaClusterCollectionMerger ():
    DataSourceProcessor("AlibavaClusterCollectionMerger"),
    //telescope
    _telescopeFileName(ALIBAVA::NOTSET),
    _telescopePulseCollectionName(ALIBAVA::NOTSET),
    _telescopeSparseCollectionName(ALIBAVA::NOTSET),
    // alibava
    _alibavaFileName(ALIBAVA::NOTSET),
    _alibavaReferenceFileName(ALIBAVA::NOTSET),
    _alibavaPulseCollectionName(ALIBAVA::NOTSET),
    _alibavaSparseCollectionName(ALIBAVA::NOTSET),
    // output
    _outputPulseCollectionName(ALIBAVA::NOTSET),
    _outputSparseCollectionName(ALIBAVA::NOTSET),
    _eventIDDiff(-1),
    _usedSensorIDs(),
    _refPresent(false),
    _refSensorID(-1),
    getReferenceSensorID(nullptr),
    getSensorID(nullptr),
    _checkAlibavaReferenceEvent(nullptr)

{
    // initialize few variables
	
    _description = "Merges alibava (both DUT and Reference if present) and telescope cluster collections";
    ///////////////
    // Telescope //
    ///////////////
    
    // file name
    registerProcessorParameter("InputTelescopeFileName", 
            "This is the input file name that telescope cluster collections stored",
            _telescopeFileName, string("runXXXXXX.slcio") );

    // pulse collection
    registerInputCollection (LCIO::TRACKERPULSE, "TelescopeClusterPulseCollectionName",
            "Name of the cluster pulse collection of telescope data",
            _telescopePulseCollectionName, string("telescope_cluster_pulse") );   
    // sparse collection
    registerInputCollection (LCIO::TRACKERDATA, "TelescopeSparseClusterCollectionName",
            "Name of the sparse cluster collection of telescope data",
            _telescopeSparseCollectionName, string("telescope_sparse_cluster") );

    /////////////
    // Alibava //
    /////////////
    // file name
    registerProcessorParameter("InputAlibavaFileName", 
            "This is the input file name that alibava cluster collections stored",
            _alibavaFileName, std::string("runXXXXXX.slcio") );
    
    registerProcessorParameter("InputAlibavaReferenceFileName", 
            "This is the input file name for the alibava used as reference sensor and the "\
            "cluster collections are stored. If the string is empty, will act as there is no reference sensor",
            _alibavaReferenceFileName, std::string("") );
    
    // pulse collection
    registerInputCollection(LCIO::TRACKERPULSE, "AlibavaClusterPulseCollectionName",
            "Name of the cluster pulse collection of alibava data",
            _alibavaPulseCollectionName, string("alibava_cluster_pulse") );
    
    // sparse collection
    registerInputCollection(LCIO::TRACKERDATA, "AlibavaSparseClusterCollectionName",
            "Name of the sparse cluster collection of alibava data",
            _alibavaSparseCollectionName, string("alibava_sparse_cluster") );
    ////////////
    // Output //
    ////////////    
    // pulse collection
    registerOutputCollection (LCIO::TRACKERPULSE, "OutputClusterPulseCollectionName",
            "Name of the merged/output cluster pulse collection",
            _outputPulseCollectionName, string("merged_cluster_pulse") );
    // sparse collection
    registerOutputCollection(LCIO::TRACKERDATA, "OutputSparseClusterCollectionName",
            "Name of the merged/output sparse cluster collection. DO NOT Change this. This is hard coded in other  ",
            _outputSparseCollectionName, string("original_zsdata") );
 
    registerProcessorParameter ("EventIDDifference",
            "AlibavaEventNumber - TelescopeEventNumber",_eventIDDiff , int(0));
}

AlibavaClusterCollectionMerger * AlibavaClusterCollectionMerger::newProcessor () 
{
    return new AlibavaClusterCollectionMerger;
}



void AlibavaClusterCollectionMerger::init () 
{
    printParameters ();

    // Whether a reference alibava should be used
    if(_alibavaReferenceFileName == "")
    {
        _checkAlibavaReferenceEvent = [] (LCReader* /* noop */, AlibavaEventImpl* /* noop*/) -> bool { return nullptr; };
        _refPresent = false;
    }
    else
    {
        _checkAlibavaReferenceEvent = [] (LCReader* reader, AlibavaEventImpl* & updatedEvt) -> bool 
                { 
                   updatedEvt = static_cast<AlibavaEventImpl*>(reader->readNextEvent()); 
                   if(updatedEvt == nullptr ) return false;
                   else return true;
                };
        _refPresent = true;
        // List the available sensor IDs for the telescope and DUT alibava
        getSensorID = [this] (const int & id) -> int
            {
                this->_usedSensorIDs.insert(id);
                return id;
            };
        // To be used initialy, afterwards, it will call just the datamember
        this->getReferenceSensorID = [this] () -> int 
            {
                if(this->_refSensorID==-1)
                {
                    _refSensorID = *(_usedSensorIDs.rbegin())+1;
                    // Update function to return the update value
                    this->getSensorID = [](const int & id) { return id; };
                }
                return _refSensorID;
            };
    }    
}

void AlibavaClusterCollectionMerger::readDataSource(int /* numEvents */) 
{
	// open telescope file
	LCReader* telescope_lcReader = LCFactory::getInstance()->createLCReader();
        try 
        {
            telescope_lcReader->open(_telescopeFileName );
        }
        catch( IOException& e )
        {
            streamlog_out ( ERROR1 ) << "Can't open the telescope file: " << e.what() << std::endl ;
	}
	
	// open alibava file
	LCReader* alibava_lcReader = LCFactory::getInstance()->createLCReader();
	try
        {
            alibava_lcReader->open( _alibavaFileName );
        }
        catch( IOException& e )
        {
            streamlog_out ( ERROR1 ) << "Can't open the alibava file: " << e.what() << std::endl ;
        }
	// open alibava reference file, if provided
	LCReader* alibavaRef_lcReader = LCFactory::getInstance()->createLCReader();
        if(_refPresent)
        {
            try
            {
                alibavaRef_lcReader->open( _alibavaReferenceFileName );
            }
            catch( IOException& e )
            {
                streamlog_out ( ERROR1 ) << "Can't open the alibava file: " << e.what() << std::endl ;
            }
        }

	// we will copy alibava run header to as the header of output file (use the DUT).
	try
        {
            LCRunHeader* alibava_runHeader = alibava_lcReader->readNextRunHeader();
            ProcessorMgr::instance()->processRunHeader( alibava_runHeader ) ;
	}
        catch( IOException& e )
        {
            streamlog_out ( ERROR1 ) << "Can't access run header of the alibava file: " << e.what() << std::endl ;
            return;
	}
	
	int eventCounter=0;
	EUTelEventImpl*  telescopeEvent=nullptr;
	AlibavaEventImpl* alibavaEvent=nullptr;
	AlibavaEventImpl* alibavaRefEvent=nullptr;
	if(_eventIDDiff<0 )
        {
            for (int i=0; i<abs(_eventIDDiff); i++)
            {
                telescopeEvent = static_cast<EUTelEventImpl*> ( telescope_lcReader->readNextEvent() );
            }
	}
	else if(_eventIDDiff>0)
        {
            for(int i=0; i<_eventIDDiff; i++)
            {
                alibavaEvent    = static_cast<AlibavaEventImpl*> ( alibava_lcReader->readNextEvent() );
                if(_refPresent)
                {
                    alibavaRefEvent = static_cast<AlibavaEventImpl*>(alibavaRef_lcReader->readNextEvent());
                }
            }
	}
	
        // [XXX: WARNING; if there is some problem with any of the input files regarding the expected
        //       collections, the threw exception is not going to be seen...
	bool noCollectionFound = false;
	while( ((telescopeEvent = static_cast<EUTelEventImpl*> ( telescope_lcReader->readNextEvent())) != 0 ) 
                && ((alibavaEvent = static_cast<AlibavaEventImpl*> (alibava_lcReader->readNextEvent())) != 0 )
                && _checkAlibavaReferenceEvent(alibavaRef_lcReader,alibavaRefEvent) )
	{
            noCollectionFound = false;
            if(telescopeEvent->getEventType() == kEORE)
            { 
                streamlog_out ( MESSAGE5 ) << "Reached EORE of telescope data"<< std::endl;
                break;		
            }
	    if( eventCounter % 1000 == 0 )
            {
                streamlog_out ( MESSAGE4 ) << "Looping events "<<alibavaEvent->getEventNumber() << std::endl;
            }
	    
            LCCollectionVec * alibavaPulseColVec = nullptr;
            LCCollectionVec * alibavaSparseColVec = nullptr;
            try
            {
                // get alibava collections
		alibavaPulseColVec = dynamic_cast<LCCollectionVec*>( alibavaEvent->getCollection( _alibavaPulseCollectionName ) ) ;
		alibavaSparseColVec = dynamic_cast<LCCollectionVec*>( alibavaEvent->getCollection( _alibavaSparseCollectionName ) ) ;
            }
            catch ( DataNotAvailableException& e ) 
            {
                noCollectionFound = true;
                streamlog_out( WARNING5 ) <<"[DUT Alibava]: No input collection " 
                    << _alibavaPulseCollectionName << " or "<< _alibavaSparseCollectionName 
                    << " found on alibava event " << alibavaEvent->getEventNumber() << std::endl;
            }
            
            LCCollectionVec * alibavaRefPulseColVec = nullptr;
            LCCollectionVec * alibavaRefSparseColVec = nullptr;
            if(_refPresent)
            {
                try
                {
                    // get alibava collections
                    alibavaRefPulseColVec = dynamic_cast<LCCollectionVec*>( alibavaRefEvent->getCollection( _alibavaPulseCollectionName ) ) ;
                    alibavaRefSparseColVec = dynamic_cast<LCCollectionVec*>( alibavaRefEvent->getCollection( _alibavaSparseCollectionName ) ) ;
                }
                catch (DataNotAvailableException& e ) 
                {
                    noCollectionFound = true;
                    streamlog_out( WARNING5 ) <<"[Reference Alibava]: No input collection " 
                        << _alibavaPulseCollectionName << " or "<< _alibavaSparseCollectionName 
                        << " found on alibava event " << alibavaEvent->getEventNumber() << std::endl;
                }
            }
			
            LCCollectionVec * telescopePulseColVec = nullptr;
            LCCollectionVec * telescopeSparseColVec = nullptr;
	    try
            {
                // get telescope collections
		telescopePulseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopePulseCollectionName ) ) ;
		telescopeSparseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopeSparseCollectionName ) ) ;
            }
            catch( DataNotAvailableException& e ) 
            {
                noCollectionFound = true;
                streamlog_out( WARNING5 ) <<"No input collection " 
                    << _telescopePulseCollectionName << " or "<< _telescopeSparseCollectionName 
                    << " found on telescope event " << telescopeEvent->getEventNumber() << std::endl;
            }
		
	    // create output collections
	    LCCollectionVec * outputPulseColVec = new LCCollectionVec(LCIO::TRACKERPULSE);
	    LCCollectionVec * outputSparseColVec = new LCCollectionVec(LCIO::TRACKERDATA);
		
	    if( !noCollectionFound ) 
            {
                // copy telescope clusters
		copyClustersInCollection(outputPulseColVec, outputSparseColVec, telescopePulseColVec, telescopeSparseColVec);
                // copy alibava cluster
                copyClustersInCollection(outputPulseColVec, outputSparseColVec, alibavaPulseColVec, alibavaSparseColVec);
                // copy alibava Ref cluster
                if(_refPresent)
                {
                    copyClustersInCollection(outputPulseColVec, outputSparseColVec, alibavaRefPulseColVec, alibavaRefSparseColVec,true);
                }
            }
            
            try
            {
                AlibavaEventImpl* outputEvent = new AlibavaEventImpl();
                outputEvent->setRunNumber( alibavaEvent->getRunNumber() );
                outputEvent->setEventNumber(eventCounter);
                outputEvent->setEventType( alibavaEvent->getEventType() );
                outputEvent->setEventSize( alibavaEvent->getEventSize() );
                outputEvent->setEventValue( alibavaEvent->getEventValue() );
                outputEvent->setEventTime( alibavaEvent->getEventTime() );
                outputEvent->setEventTemp( alibavaEvent->getEventTemp() );
                outputEvent->setCalCharge( alibavaEvent->getCalCharge() );
                outputEvent->setCalDelay( alibavaEvent->getCalDelay() );
                if(alibavaEvent->isEventMasked())
                {
                    outputEvent->maskEvent();
                }
		else
                {
                    outputEvent->unmaskEvent();
                }

	        if( !noCollectionFound ) 
                {
                    outputEvent->addCollection(outputPulseColVec, _outputPulseCollectionName);
                    outputEvent->addCollection(outputSparseColVec, _outputSparseCollectionName);
                }
                ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( outputEvent ) ) ;
                // Free memory
                delete outputEvent;
                outputEvent = nullptr;
                //	streamlog_out ( MESSAGE1 ) << "Successfully copied Alibava collections to output event" << endl ;
            }
            catch ( IOException& e) 
            {
                // do nothing again
                streamlog_out( ERROR5 ) << e.what() << std::endl;
            }
            
            eventCounter++;
        }// end of loop over events	
}


void AlibavaClusterCollectionMerger::copyClustersInCollection(LCCollectionVec * outputPulseColVec, 
        LCCollectionVec * outputSparseColVec, LCCollectionVec * inputPulseColVec, 
        LCCollectionVec * inputSparseColVec,
        bool isReferenceSensor)
{
    // Here is the Cell ID Encodes for pulseFrame and sparseFrame
    // CellID Encodes are introduced in eutelescope::EUTELESCOPE
	
    // for sparseFrame (usually called cluster collection)
    CellIDEncoder<TrackerDataImpl> outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputSparseColVec );
    // for pulseFrame
    CellIDEncoder<TrackerPulseImpl> outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputPulseColVec );
	
    // for input
    CellIDDecoder<TrackerPulseImpl> inputPulseColDecoder(inputPulseColVec);
    CellIDDecoder<TrackerDataImpl> inputSparseColDecoder(inputSparseColVec);
	
    // go through input clusters and copy them to output cluster collection
    unsigned int noOfClusters = inputPulseColVec->getNumberOfElements();
    for(unsigned int i = 0; i < noOfClusters; ++i )
    {
        TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl();
        TrackerDataImpl * outputSparseFrame = new TrackerDataImpl();
        
        TrackerPulseImpl* inputPulseFrame = dynamic_cast<TrackerPulseImpl*>(inputPulseColVec->getElementAt(i));
        TrackerDataImpl* inputSparseFrame = dynamic_cast<TrackerDataImpl*>(inputPulseFrame->getTrackerData());
        
        // Note the reference sensor has to be calling this function the last one, 
        // in that case, all the sensors IDs have been used, so first time is 
        // going to obtain the lowest ID number available to be used
        int sensorID(-1);
        if(isReferenceSensor)
        {
            sensorID = this->getReferenceSensorID();
        }
        sensorID = this->getSensorID(static_cast<int>(inputSparseColDecoder(inputSparseFrame)["sensorID"]));
        // set Cell ID for sparse collection
        outputSparseColEncoder["sensorID"] = sensorID;
	outputSparseColEncoder["sparsePixelType"] =static_cast<int>(inputSparseColDecoder(inputSparseFrame)["sparsePixelType"]);
	outputSparseColEncoder["quality"] = static_cast<int>(inputSparseColDecoder(inputSparseFrame)["quality"]);
	outputSparseColEncoder.setCellID( outputSparseFrame );
		
	// copy tracker data
	outputSparseFrame->setChargeValues(inputSparseFrame->getChargeValues());
	// add it to the cluster collection
	outputSparseColVec->push_back( outputSparseFrame );
		
	// prepare a pulse for this cluster 
        // (should be same sensor ID, right?)
        outputPulseColEncoder["sensorID"] = sensorID; //static_cast<int> (inputPulseColDecoder(inputPulseFrame) ["sensorID"]);
	outputPulseColEncoder["type"] = static_cast<int>(inputPulseColDecoder(inputPulseFrame) ["type"]);
	outputPulseColEncoder.setCellID( outputPulseFrame );
		
	outputPulseFrame->setCharge( inputPulseFrame->getCharge() );
	outputPulseFrame->setTrackerData( outputSparseFrame);
        outputPulseColVec->push_back( outputPulseFrame );
    } // end of loop over input clusters	
}


void AlibavaClusterCollectionMerger::end () 
{
    streamlog_out ( MESSAGE5 )  << "AlibavaClusterCollectionMerger Successfully finished" << endl;
}

