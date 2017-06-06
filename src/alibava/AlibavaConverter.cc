/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 *    - Fix data binary types (time_t -> unsigned int) for the
 *      time of start of the run
 *
 */



// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaConverter.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"

// marlin includes
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"


// lcio includes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

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
using namespace alibava;

AlibavaConverter::AlibavaConverter ():
    DataSourceProcessor("AlibavaConverter"),
    _fileName(ALIBAVA::NOTSET),
    _geoID(0),
    _runNumber(-1),
    _rawDataCollectionName("rawdata"),
    _rawChipHeaderCollectionName("chipheader"),
    _chipSelection(),
    _startEventNum(-1),
    _stopEventNum(-1),
    _storeHeaderPedestalNoise(false)
{
    // initialize few variables
    _description = "Reads data streams produced by Alibava and produces the corresponding LCIO output";
    registerProcessorParameter("InputFileName", "This is the input file name",
            _fileName, string("runXXXXXX.dat") );
    registerProcessorParameter("GeoID", "The geometry identification number", _geoID, static_cast<int> ( 0 ));
    registerProcessorParameter("RunNumber", "Run number of file (formatted)",_runNumber, int(-1) );
    registerOutputCollection (LCIO::TRACKERDATA, "RawDataCollectionName",
            "Name of the collection",_rawDataCollectionName, std::string("rawdata") );
    registerOutputCollection (LCIO::TRACKERDATA, "RawChipHeaderCollectionName",
            "Name of the collection",_rawChipHeaderCollectionName, std::string("chipheader") );
    
    // now optional parameters
    registerOptionalParameter("ChipSelection", 
            "Selection of chip that you want to store data from. Chip numbers "\
            "start from 0. If not set, all data (i.e. chip 0 and 1) will be stored",
            _chipSelection, EVENT::IntVec() );
    registerOptionalParameter("StartEventNum", "The event number that AlibavaConverter "\
            "should start storing. Default value is -1, in this case it will store every event",
            _startEventNum, int(-1) );
    registerOptionalParameter("StopEventNum", "The event number that AlibavaConverter "\
            "should stop storing. Default value is -1, in this case it will store every event",
            _stopEventNum, int(-1) );
    registerOptionalParameter("StoreHeaderPedestalNoise", "Alibava stores a pedestal "\
            "and noise set in the run header. These values are not used in te rest of "\
            "the analysis, so it is optional to store it. By default it will not be "\
            "stored, but it you want you can set this variable to true to store it "\
            "in the header of slcio file",
            _storeHeaderPedestalNoise, bool(false) );
}

AlibavaConverter * AlibavaConverter::newProcessor() 
{
    return new AlibavaConverter;
}


void AlibavaConverter::init() 
{
    checkIfChipSelectionIsValid();
    if(_startEventNum!=-1 && _stopEventNum==-1)
    {
        streamlog_out( WARNING5 )<< "First "<< _startEventNum 
            << " will be skipped!" <<endl;
    }
    else if(_startEventNum!=-1 && _stopEventNum!=-1)
    {
        streamlog_out( WARNING5 )<< "Only events from " 
            <<_startEventNum <<" to "<<_stopEventNum <<" will be stored! "<<endl;
    }
    else if (_startEventNum==-1 && _stopEventNum!=-1)
    {
        streamlog_out( WARNING5 )<< "Events after event number: "<<
            _stopEventNum<<" will be skipped! "<<endl;	
    }
    printParameters();
}

void AlibavaConverter::readDataSource(int /* numEvents */) 
{
    // this event counter is used to stop the processing when it is
    // greater than numEvents.
    int eventCounter = 0;
    
    // this is to make the output messages nicer
    streamlog::logscope scope(streamlog::out);
    scope.setName(name());
    
    // print out a debug message
    streamlog_out( MESSAGE5 ) << "Reading " << _fileName 
        << " with AlibavaConverter " << _runNumber << std::endl;
    
    ////////////////
    // Open File  //
    ////////////////
    ifstream infile;
    infile.open(_fileName.c_str());
    if(!infile.is_open()) 
    {
        streamlog_out( ERROR5 ) << "AlibavaConverter could not read the file "
            <<_fileName<<" correctly. Please check the path and file names that" 
            << " have been input" << std::endl;
        exit(-1);
    }
    else
    {
        streamlog_out( MESSAGE4 )<<"Input file "<<_fileName<<" is opened!"<<std::endl;
    }
    /////////////////
    // Read Header //
    /////////////////
    // --- Time of start of the run 
    unsigned int date;
    infile.read(reinterpret_cast< char *> (&date), sizeof(unsigned int));
    // --- Run type: 1:Calibration; 2:Laser Syn.; 3:Laser; 4: Rad. source; 5: Pedestal
    int type;
    infile.read(reinterpret_cast< char *> (&type), sizeof(int));
    // --- Header lenght
    unsigned int lheader; 
    infile.read(reinterpret_cast< char *> (&lheader), sizeof(unsigned int));
    // read the next field until reach the lenght of the header
    std::string header;
    header.clear();
    for(unsigned int ic=0; ic<lheader; ic++)
    {
        char tmp_c;
        infile.read(&tmp_c, sizeof(char));
        header.append(1, tmp_c);
    }
    header = trim_str(header);
    
    // Firmware version
    int version(-1);
    if (header[0]!='V' && header[0]!='v')
    {
        version = 0;
    }
    else
    {
        version = int(header[1]-'0');
        header = header.substr(5);
    }
    
    ////////////////////////////////////
    // Read header pedestal and noise //
    // /////////////////////////////////
	
    // Alibava stores a pedestal and noise set in the run header. 
    // These values are not used in te rest of the analysis, so it is 
    // optional to store it. By default it will not be stored, but it you
    // want you can set _storeHeaderPedestalNoise variable to true.
    float tmp_float;
    FloatVec headerPedestal;
    // first pedestal
    for(int ichan=0; ichan<ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS; ichan++) 
    {
        infile.read(reinterpret_cast< char *> (&tmp_float), sizeof(double));
        headerPedestal.push_back(tmp_float);
    }	
    // now noise
    FloatVec headerNoise;
    for(int ichan=0; ichan<ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS; ichan++) 
    {
        infile.read(reinterpret_cast< char *> (&tmp_float), sizeof(double));
        headerNoise.push_back(tmp_float);
    }
    ////////////////////
    // Process Header //
    ////////////////////
    LCRunHeaderImpl * arunHeader = new LCRunHeaderImpl();
    AlibavaRunHeaderImpl* runHeader = new AlibavaRunHeaderImpl(arunHeader);
    runHeader->setDetectorName(Global::GEAR->getDetectorName());
    runHeader->setHeader(header);
    runHeader->setHeaderVersion(version);
    runHeader->setDataType(type);
    runHeader->setDateTime(string(ctime(reinterpret_cast<time_t*>(&date))));
    if(_storeHeaderPedestalNoise) 
    {
        runHeader->setHeaderPedestal(headerPedestal);
        runHeader->setHeaderNoise(headerNoise);
    }
    runHeader->setRunNumber(_runNumber);
    runHeader->setChipSelection(_chipSelection);
	
    // get number of events from header (only valid if the event is
    // different from calibration type)
    std::string tmpstring = getSubStringUpToChar(header,";",0);
    int noofevents = atoi(tmpstring.c_str());
    runHeader->setNoOfEvents(noofevents);
    //runHeader->addProcessor(type());
	
    ProcessorMgr::instance()->processRunHeader( runHeader->lcRunHeader() ) ;
    
    delete arunHeader;
    delete runHeader;

    ////////////////
    // Read Event //
    ////////////////
    // Expected different streams depending the firmware version
    // used to store the data
    if(version<2) 
    {
        // this code is not written for version<=1.
	streamlog_out( ERROR5 )<<" Unexpected data version found (version="<<version<<"<2). Data is not saved"<<std::endl;
	return;
    }
	
    do
    {
        if( eventCounter % 1000 == 0 )
        {
            streamlog_out ( MESSAGE4 ) << "Processing event "<< eventCounter << " in run " << _runNumber<<std::endl;
        }
	
        unsigned int headerCode, eventTypeCode=0;
        bool breaktheloop=false;
        do
        {
            infile.read(reinterpret_cast< char *> (&headerCode), sizeof(unsigned int));
            if(infile.bad() || infile.eof())
            {
                // [JDC] End of file reached, so breaking the loop
                // in a clean way
                breaktheloop=true;
                break;
            }
	    eventTypeCode = (headerCode>>16) & 0xFFFF;
        } while( eventTypeCode != 0xcafe );
        if(breaktheloop)
        {
            // [JDC] Breaking the loop in a clean way
            break;
        }
	
        eventTypeCode = headerCode & 0x0fff;

        unsigned int userEventTypeCode = headerCode & 0x1000;
        if(userEventTypeCode)
        {	
            streamlog_out( ERROR5 )<<" Unexpected data type found (type= User type). Data is not saved"<<std::endl;
            return;
        }
        unsigned int eventSize;
        infile.read(reinterpret_cast< char *> (&eventSize), sizeof(unsigned int));
        
        double value;
        infile.read(reinterpret_cast< char *> (&value), sizeof(double));
		
	//see AlibavaGUI.cc
        double charge = int(value) & 0xff;
        double delay = int(value) >> 16;
        charge = charge * 1024;
	
        // The timestamp
        unsigned int clock;
        // Thomas 13.05.2015: Firmware 3 introduces the clock to the header!
        // for now this is not stored...
        if (version==3)
        {
            infile.read(reinterpret_cast< char *> (&clock), sizeof(unsigned int));
        }
        // Time digital converter 
        unsigned int tdcTime;
        infile.read(reinterpret_cast< char *> (&tdcTime), sizeof(unsigned int));
        // Temperature measured on daughter board
	unsigned short temp;  // temperature measured on Daughter board        
        infile.read(reinterpret_cast< char *> (&temp), sizeof(unsigned short));
        
        // vector for data
        FloatVec all_data;
        all_data.reserve(ALIBAVA::NOOFCHIPS*ALIBAVA::NOOFCHANNELS);
        // vector for chip header
        FloatVec all_chipheaders;
        all_chipheaders.reserve(ALIBAVA::NOOFCHIPS*ALIBAVA::CHIPHEADERLENGTH);
        
        // An auxiliary variable (contains the data before pushing in the 
        // previous vector
        short tmp_short;
        // Chip header
        unsigned short chipHeader[2][ALIBAVA::CHIPHEADERLENGTH];
        // iterate over number of chips and store the chip headers 
        // and data
	FloatVec::iterator it;
        for(int ichip=0; ichip < ALIBAVA::NOOFCHIPS; ++ichip)
        {
            infile.read(reinterpret_cast< char *> (chipHeader[ichip]), ALIBAVA::CHIPHEADERLENGTH*sizeof(unsigned short));            
            // store chip header in all_chipheaders vector
            streamlog_out (DEBUG0) << "chip " << ichip << " header: " ;
            for(int j = 0; j < ALIBAVA::CHIPHEADERLENGTH; ++j)
            {
                streamlog_out (DEBUG0) << " " << chipHeader[ichip][j];
                all_chipheaders.push_back(float(chipHeader[ichip][j]));
            }
            streamlog_out (DEBUG0) << std::endl;            
            // store data in all_data vector
	    for(int ichan=0; ichan < ALIBAVA::NOOFCHANNELS; ichan++) 
            {
                infile.read(reinterpret_cast< char *> (&tmp_short), sizeof(unsigned short));
                all_data.push_back(float(tmp_short));
            }
        }
        ///////////////////
        // Process Event //
        ///////////////////
	// now write these to AlibavaEvent
        AlibavaEventImpl* anEvent = new AlibavaEventImpl();
        anEvent->setRunNumber(_runNumber);
        anEvent->setEventNumber(eventCounter);
        anEvent->setEventType(eventTypeCode);
        anEvent->setEventSize(eventSize);
        anEvent->setEventValue(value);
        if(version==3)
        {
            anEvent->setEventClock(clock);
        }
        anEvent->setEventTime(tdc_time(tdcTime));
	anEvent->setEventTemp(get_temperature(temp));
	anEvent->setCalCharge(charge);
	anEvent->setCalDelay(delay);
	anEvent->unmaskEvent();
        
        // creating LCCollection for raw data
        LCCollectionVec* rawDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
        CellIDEncoder<TrackerDataImpl> chipIDEncoder(ALIBAVA::ALIBAVADATA_ENCODE,rawDataCollection);
 
        // creating LCCollection for raw chip header
        LCCollectionVec* rawChipHeaderCollection = new LCCollectionVec(LCIO::TRACKERDATA);
        CellIDEncoder<TrackerDataImpl> chipIDEncoder2(ALIBAVA::ALIBAVADATA_ENCODE,rawChipHeaderCollection);
        
	// for this to work the _chipselection has to be sorted in ascending order!!!
	for(unsigned int ichip=0; ichip<_chipSelection.size(); ichip++) 
        {
            // store raw data
	    FloatVec chipdata;
            chipdata.clear();
            // separate data for each chip
            chipdata.insert(chipdata.end(), 
                    all_data.begin()+_chipSelection[ichip]*ALIBAVA::NOOFCHANNELS, 
                    all_data.begin()+(_chipSelection[ichip]+1)*ALIBAVA::NOOFCHANNELS);
            TrackerDataImpl * arawdata = new TrackerDataImpl();
            arawdata->setChargeValues(chipdata);
            chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = _chipSelection[ichip];
            chipIDEncoder.setCellID(arawdata);
            rawDataCollection->push_back(arawdata);

            // store chip header
            FloatVec chipHeader_vec;
            chipHeader_vec.clear();
            //separate chip headers for each chip
            chipHeader_vec.insert(chipHeader_vec.end(), 
                    all_chipheaders.begin()+_chipSelection[ichip]*ALIBAVA::CHIPHEADERLENGTH, 
                    all_chipheaders.begin()+(_chipSelection[ichip]+1)*ALIBAVA::CHIPHEADERLENGTH);
            TrackerDataImpl * achipheader = new TrackerDataImpl();
            achipheader->setChargeValues(chipHeader_vec);
            chipIDEncoder2[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = _chipSelection[ichip];
            chipIDEncoder2.setCellID(achipheader);
            rawChipHeaderCollection->push_back(achipheader);
        }
        anEvent->addCollection(rawDataCollection, _rawDataCollectionName);
        anEvent->addCollection(rawChipHeaderCollection,_rawChipHeaderCollectionName);
        if(_startEventNum!=-1 && eventCounter<_startEventNum) 
        {
            streamlog_out( MESSAGE5 )<<" Skipping event "<<eventCounter<<". StartEventNum is set to "
                <<_startEventNum<<std::endl;
            ++eventCounter;
            continue;
        }
        if(_stopEventNum!=-1 && eventCounter>_stopEventNum) 
        {
            streamlog_out( MESSAGE5 )<<" Reached StopEventNum: "<<_stopEventNum 
                <<". Last saved event number is "<<eventCounter<<std::endl;
            break;
        }
        ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( anEvent ) ) ;
        ++eventCounter;
        // Free memory
        delete anEvent;
    }while( !(infile.bad() || infile.eof()) );
    
    infile.close();
    if(_stopEventNum!=-1 && eventCounter<_stopEventNum)
    {
        streamlog_out( MESSAGE5 )<<" Stopped before reaching StopEventNum: "
            <<_stopEventNum<<". The file has "<<eventCounter<<" events."<<std::endl;
    }
}


void AlibavaConverter::end() 
{
    streamlog_out ( MESSAGE5 )  << "AlibavaConverter Successfully finished" << std::endl;
}

double AlibavaConverter::tdc_time(unsigned int tdcTime)
{
    unsigned short fpart = tdcTime & 0xffff;
    short ipart = (tdcTime & 0xffff0000)>>16;
    if (ipart<0)
    {
        fpart *= -1;
    }
    //double tt = 100.*(1. -(ipart + (fpart/65535.)));
    double tt = 100.0*(ipart + (fpart/65535.));
    return tt;
}

double AlibavaConverter::get_temperature(unsigned short temp)
{
    if(temp==0)
    {
        return 9999.;
    }
    else
    {
        return 0.12*temp - 39.8;
    }
}

void AlibavaConverter::checkIfChipSelectionIsValid()
{
    bool resetChipSelection = false;
    
    // check if there is chip selection or if there is chip selection but not valid
    if(_chipSelection.size()==0)
    {
        streamlog_out( WARNING5 )<< "You didn't select any chip"<<std::endl;
        resetChipSelection = true;
    }
    else if( int(_chipSelection.size())> ALIBAVA::NOOFCHIPS) 
    {
        streamlog_out( WARNING5 )<< "The number of chips you selected ( "
            <<_chipSelection.size()<<" ) is more than an alibava daughter board can have ( "
            <<ALIBAVA::NOOFCHIPS<<" )"<<std::endl;
        resetChipSelection= true;
    }
    else
    {
        // first sort the chip numbers in ascending order
        std::sort(_chipSelection.begin(),_chipSelection.end());
        // check if the selected chips make sense
        for(int ichip=0; ichip<int(_chipSelection.size()); ++ichip) 
        {
            bool del_this_chip = false;
            if(_chipSelection[ichip]<0)
            {
                streamlog_out( ERROR5 )<< "Selected chip cannot have negative value. "<<std::endl;
                del_this_chip =true;
            }
            else if(_chipSelection[ichip]>=ALIBAVA::NOOFCHIPS)
            {
                streamlog_out( ERROR5 )<< "Chip numbering has to start from zero \"0\" "
                    << "and cannot be greater than " << ALIBAVA::NOOFCHIPS-1 << std::endl;
                del_this_chip = true;
            }
            
            // if this chip selection is not valid, delete it.
            if(del_this_chip) 
            {
                streamlog_out( ERROR5 )<< "Chip "<<_chipSelection[ichip]
                    <<" is deleted from the chip selection list"<< std::endl;
                _chipSelection.erase(_chipSelection.begin()+ichip);
                ichip = ichip-1;
            }
        }// end of ichip loop
	
        // check again if there is any selected chip left
        if(_chipSelection.size()==0)
        {
            resetChipSelection= true;
        }
    }
    
    if (resetChipSelection) 
    {
        streamlog_out( WARNING5 )<< "I will save data from all chips"<< std::endl;
        _chipSelection.clear();
        for (int ichip=0; ichip<ALIBAVA::NOOFCHIPS; ichip++) 
        {
            _chipSelection.push_back(ichip);
        }
    }

    // now there is valid chip selection!!!
    streamlog_out( WARNING5 )<< "Final applied chip selection: ";
    for(int ichip=0; ichip<int(_chipSelection.size()); ++ichip)
    {
        streamlog_out( WARNING5 )<< _chipSelection[ichip] << " ";
    }
    streamlog_out( WARNING5 )<<std::endl;
    streamlog_out( WARNING5 )<< "Only data coming from these chips will be stored"<<std::endl;	
}

