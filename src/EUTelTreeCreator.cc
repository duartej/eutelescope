/*
 * Created by Eda Yildirim
 *  (2015 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// ROOT includes:

// eutelescope includes ".h"
#include "EUTelTreeCreator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelGeometryTelescopeGeoDescription.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// marlinutil includes
#include "marlinutil/SimpleLine.h"

// gear includes <.h>


// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// ROOT
#include "TFile.h"
#include "TTree.h"

// system includes <>
#include <cmath>
//#include <memory>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace marlin;
using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
EUTelTreeCreator::EUTelTreeCreator():
    Processor("EUTelTreeCreator"),
    _inputHitCollectionName(),
    _inputTrackCollectionName(),
    _telescopePlanes(),
    _dutPlanes(),
    _file(nullptr),
    _filename(""),
    _tree(nullptr),
    _branches(),
    _branches_I(),
    _branches_int(),
    _iRun(0),
    _iEvt(0)
{
    // modify processor description
    _description =  "EUTelTreeCreator Creates a simple tree from a merged collection"\
                     " of hits (where they are included the DUT hits) and tracks";

    registerInputCollection(LCIO::TRACKERHIT,"InputHitCollectionName", "Input hit collection name."\
            " Hits should be in global coordinates and pre-aligned", _inputHitCollectionName, std::string(""));
    
    registerInputCollection(LCIO::TRACK,"InputTrackCollectionName", "Input track collection name to be used for the "\
            "DUT missing coordinate search", _inputTrackCollectionName, std::string(""));
    
    registerProcessorParameter("OutputFilename", "The ntuple ROOT filename", _filename, std::string("track_duts_ntuple.root"));
    
    registerProcessorParameter("TelescopePlanes", "The id of the telescope planes",
            _telescopePlanes, std::vector<int>({0,1,2,3,4}));
    
    registerProcessorParameter("DUTPlanes", "The id of the DUT planes",
            _dutPlanes, std::vector<int>({5,6}));
}

/*EUTelTreeCreator::~EUTelTreeCreator()
{
    // Free memory
    for(auto & el: _branches)
    {
        if(el.second != nullptr)
        {
            delete el.second;
            el.second = nullptr;
        }
    }
    
    for(auto & el: _branches_I)
    {
        if(el.second != nullptr)
        {
            delete el.second;
            el.second = nullptr;
        }
    }
}*/

void EUTelTreeCreator::clear_branches()
{
    for(auto & el: _branches)
    {
        if(el.second != nullptr)
        {
            el.second->clear();
        }
    }
    
    for(auto & el: _branches_I)
    {
        if(el.second != nullptr)
        {
            el.second->clear();
        }
    }
}

void EUTelTreeCreator::init() 
{
    // Create the files and the tree
    _file = new TFile(_filename.c_str(),"RECREATE");

    _tree = new TTree("events","Ntuple events");

    // Event related
    _branches_int["EventNumber"] = _iEvt;
    _branches_int["RunNumber"]   = _iRun;
    _branches_int["ntracks"]     = 0;
    // DUT hits position
    // the sensor plane ID
    _branches_I["dut_hit_plane"] = new std::vector<int>;
    // DUT hits position
    _branches["hit_X_dut"] = new std::vector<float>;
    _branches["hit_Y_dut"] = new std::vector<float>;
    _branches["hit_Z_dut"] = new std::vector<float>;
    _branches["hit_XLocal_dut"] = new std::vector<float>;
    _branches["hit_YLocal_dut"] = new std::vector<float>;
    _branches["hit_ZLocal_dut"] = new std::vector<float>;

    // Create the elements to fill the branches
    // track related
    _branches["trk_chi2"] = new std::vector<float>;
    _branches["trk_Ndf"]  = new std::vector<float>;
    _branches["trk_dxdz"] = new std::vector<float>;
    _branches["trk_dydz"] = new std::vector<float>;
    
    // Hits position associated to the tracks
    // The corresponding track element
    _branches_I["trk_hit_index"]= new std::vector<int>;
    // The sensor plane of the hit
    _branches_I["trk_hit_plane"]= new std::vector<int>;
    // the position
    _branches["trk_hit_measX"] = new std::vector<float>;
    _branches["trk_hit_measY"] = new std::vector<float>;
    _branches["trk_hit_measZ"] = new std::vector<float>;
    _branches["trk_hit_fitX"]  = new std::vector<float>;
    _branches["trk_hit_fitY"]  = new std::vector<float>;
    _branches["trk_hit_fitZ"]  = new std::vector<float>;
    // The equivalent in the local coordiante system (at the sensor plane)
    _branches["trk_hit_measXLocal"] = new std::vector<float>;
    _branches["trk_hit_measYLocal"] = new std::vector<float>;
    _branches["trk_hit_measZLocal"] = new std::vector<float>;
    _branches["trk_hit_fitXLocal"]  = new std::vector<float>;
    _branches["trk_hit_fitYLocal"]  = new std::vector<float>;
    _branches["trk_hit_fitZLocal"]  = new std::vector<float>;
    

    // And attach them to the tree
    for(auto & el: _branches)
    {
        _tree->Branch(el.first.c_str(),&el.second);
    }

    for(auto & el: _branches_I)
    {
        _tree->Branch(el.first.c_str(),&el.second);
    }

    for(auto & el: _branches_int)
    {
        _tree->Branch(el.first.c_str(),&el.second);
    }

    
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
    
    // set to zero the run and event counters
    _iRun = 0;
    _iEvt = 0;
}


void EUTelTreeCreator::processRunHeader (LCRunHeader * rdr) 
{
    std::unique_ptr<EUTelRunHeaderImpl> header = std::make_unique<EUTelRunHeaderImpl>(rdr);
    header->addProcessor(type());
    
    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, instead of barely
    // quitting ask the user what to do.
    if( header->getGeoID() == 0 )
    {
        streamlog_out ( WARNING0 ) <<  "The geometry ID in the run header is set to zero." << std::endl
            <<  "This may mean that the GeoID parameter was not set" << std::endl;
    }
    // increment the run counter
    ++_iRun;
}

void EUTelTreeCreator::processEvent (LCEvent * event) 
{
    ++_iEvt;
    
    EUTelEventImpl * evt = static_cast<EUTelEventImpl*>(event);
    
    if( evt->getEventType() == kEORE ) 
    {
        streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << std::endl;
        return;
    } 
    else if ( evt->getEventType() == kUNKNOWN ) 
    {
        streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
            << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
    }

    LCCollectionVec * inputHitCollection = nullptr;
    LCCollectionVec * inputTrackCollection = nullptr;
    try
    {
        inputHitCollection = static_cast<LCCollectionVec*>(event->getCollection( _inputHitCollectionName ));
    }
    catch(DataNotAvailableException& e)
    {
        streamlog_out  ( MESSAGE2 ) <<  "No input collection [" << _inputHitCollectionName 
            << "] found on event " << event->getEventNumber()
            << " in run " << event->getRunNumber() << std::endl;
        return ;
    }
    try
    {
        inputTrackCollection = static_cast<LCCollectionVec*>(event->getCollection(_inputTrackCollectionName));
    }
    catch (DataNotAvailableException& e  )
    {
        streamlog_out  ( MESSAGE2 ) <<  "No input collection " << _inputTrackCollectionName 
            << " found on event " << event->getEventNumber()
            << " in run " << event->getRunNumber() << std::endl;
        return ;
    }

    // Clear the elements before the filling
    clear_branches();    
    
    _branches_int["EventNumber"] = _iEvt;
    _branches_int["RunNumber"]   = _iRun;
    
    // prepare an encoder for the hit collection
    CellIDDecoder<TrackerHitImpl> inputCellIDDecoder( inputHitCollection );
    
    std::vector<TrackerHitImpl*> dutHits;

    // Here identify which hits come from the DUTs
    for(int iInputHits = 0; iInputHits < inputHitCollection->getNumberOfElements(); ++iInputHits)
    {
        TrackerHitImpl * inputHit = dynamic_cast<TrackerHitImpl*>(inputHitCollection->getElementAt(iInputHits));
        int sensorID    = inputCellIDDecoder(inputHit)["sensorID"];
        
        // check if is DUT hit
        if(std::find(_dutPlanes.begin(),_dutPlanes.end(),sensorID) != _dutPlanes.end())
        {
            // and then keep track of it
            dutHits.push_back(inputHit);
        }
    }
 
    CellIDDecoder<TrackerHit> trkHitDecoder(EUTELESCOPE::HITENCODING);
    // Loop over the tracks
    for(int itrk = 0; itrk < inputTrackCollection->getNumberOfElements(); ++itrk)
    {
        // Use the two first points at a different plane as P1 and P2
        // extracted from the tracks
        TrackImpl * track = dynamic_cast<TrackImpl*>(inputTrackCollection->getElementAt(itrk));
        // Fill track related
        _branches["trk_chi2"]->push_back(track->getChi2());
        _branches["trk_Ndf"]->push_back(track->getNdf());
        _branches["trk_dxdz"]->push_back(track->getOmega());
        _branches["trk_dydz"]->push_back(track->getPhi());
        
        // Now the hits
        TrackerHitVec hitvec = track->getTrackerHits();
        for(unsigned int k = 0 ; k < hitvec.size(); ++k)
        {
            const int sensorID = trkHitDecoder(hitvec[k])["sensorID"];
            _branches_I["trk_hit_index"]->push_back(itrk);
            _branches_I["trk_hit_plane"]->push_back(sensorID);
            std::string xname("trk_hit_");
            std::string yname("trk_hit_");
            std::string zname("trk_hit_");
            std::string xnameLocal("trk_hit_");
            std::string ynameLocal("trk_hit_");
            std::string znameLocal("trk_hit_");

            if((trkHitDecoder(hitvec[k])["properties"] & kFittedHit) > 0)
            {
                // Fitted hits
                xname += "fitX";
                yname += "fitY";
                zname += "fitZ";
                xnameLocal += "fitXLocal";
                ynameLocal += "fitYLocal";
                znameLocal += "fitZLocal";
            }
            else
            {
                // Or measured
                xname += "measX";
                yname += "measY";
                zname += "measZ";
                xnameLocal += "measXLocal";
                ynameLocal += "measYLocal";
                znameLocal += "measZLocal";
            }
      	    const double* pos = hitvec[k]->getPosition();
            _branches[xname]->push_back(static_cast<float>(pos[0]));
            _branches[yname]->push_back(static_cast<float>(pos[1]));
            _branches[zname]->push_back(static_cast<float>(pos[2]));
            // The local fit
            double pos_local[3];
            geo::gGeometry().master2Local(sensorID,pos,pos_local);
            _branches[xnameLocal]->push_back(static_cast<float>(pos_local[0]));
            _branches[ynameLocal]->push_back(static_cast<float>(pos_local[1]));
            _branches[znameLocal]->push_back(static_cast<float>(pos_local[2]));
        }
    }
    // Number of tracks
    _branches_int["ntracks"] = inputTrackCollection->getNumberOfElements();
    
    // keep track of how many duts are 
    // Loop over the DUTs
    for(auto & dutHit: dutHits)
    {
        // Extract the sensorid
        const int sensorID_dut = trkHitDecoder(dutHit)["sensorID"];
        _branches_I["dut_hit_plane"]->push_back(sensorID_dut);

      	const double* posdut = dutHit->getPosition();
        _branches["hit_X_dut"]->push_back(posdut[0]);
        _branches["hit_Y_dut"]->push_back(posdut[1]);
        _branches["hit_Z_dut"]->push_back(posdut[2]);
        // The local fit
        double posdut_local[3];
        geo::gGeometry().master2Local(sensorID_dut,posdut,posdut_local);
        _branches["hit_XLocal_dut"]->push_back(static_cast<float>(posdut_local[0]));
        _branches["hit_YLocal_dut"]->push_back(static_cast<float>(posdut_local[1]));
        _branches["hit_ZLocal_dut"]->push_back(static_cast<float>(posdut_local[2]));
    }
    // Fill the tree
    _tree->Fill();
}



void EUTelTreeCreator::end()
{
    _file->Write();
    _file->Close();
}

