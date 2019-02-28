// eutelescope inlcudes
#include "EUTelAPIXTbTrackTuple.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelUtility.h"
#include "EUTelVirtualCluster.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"

// eutelescope geometry
#include "EUTelGenericPixGeoDescr.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCObject.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

#include <set>
#include <algorithm>
#include <numeric>
//PROV --DELETE
#include <cassert>

using namespace eutelescope;

EUTelAPIXTbTrackTuple::EUTelAPIXTbTrackTuple()
    : Processor("EUTelAPIXTbTrackTuple"), _inputTrackColName(""),
      _inputTrackerHitColName(""), _inputMeasHitColName(""),
      _inputTelPulseCollectionName(""),
      _inputDutPulseCollectionName(""), _telZsColName(""), _dutZsColName(""),
      _path2file(""), _DUTIDs(std::vector<int>()), _nRun(0), _nEvt(0),
      _runNr(0), _evtNr(0), _isFirstEvent(false), _file(nullptr), _eutracks(nullptr),
      _nTrackParams(0), _xPos(nullptr), _yPos(nullptr), _dxdz(nullptr), _dydz(nullptr),
      _trackIden(nullptr), _trackNum(nullptr), _chi2(nullptr), _ndof(nullptr),
      _hitpattern(nullptr),
      _zstree(nullptr), _nPixHits(0), p_col(nullptr), p_row(nullptr), p_tot(nullptr),
      p_iden(nullptr), p_lv1(nullptr), p_hitTime(nullptr), p_frameTime(nullptr),
      _mhits(nullptr),_nmHits(0),_mhitXpos(nullptr),_mhitYpos(nullptr),
      _mhitEtaX(nullptr),_mhitEtaY(nullptr),_mhitTOT(nullptr),
      _mhitId(nullptr),_mhitBCID(nullptr),_mhitSize(nullptr),
      _mhitsSizeX(nullptr),_mhitsSizeY(nullptr),
      _euhits(nullptr), _nHits(0), _hitXPos(nullptr), _hitYPos(nullptr), _hitZPos(nullptr),
      _hitSensorId(nullptr) {
  // processor description
  _description = "Prepare tbtrack style n-tuple with track fit results";

  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName",
                          "Name of the input Track collection",
                          _inputTrackColName, std::string("fittracks"));

  registerInputCollection(LCIO::TRACKERHIT, "InputTrackerHitCollectionName",
                          "Name of the plane-wide hit-data hit collection",
                          _inputTrackerHitColName, std::string("fitpoints"));
  
  registerInputCollection(LCIO::TRACKERHIT, "InputMeasuredHitCollectionName",
                          "Name of the measured local hit collection",
                          _inputMeasHitColName, std::string("local_hit"));
  
  registerProcessorParameter("DutZsColName",
                             "DUT zero surpressed data colection name",
                             _dutZsColName, std::string("zsdata_apix"));

  registerProcessorParameter("OutputPath",
                             "Path/File where root-file should be stored",
                             _path2file, std::string("NTuple.root"));

  registerProcessorParameter("DUTIDs",
                             "Int std::vector containing the IDs of the DUTs",
                             _DUTIDs, std::vector<int>());
}

void EUTelAPIXTbTrackTuple::init() {
  // usually a good idea to
  printParameters();

  _isFirstEvent = true;

  _nRun = 0;
  _nEvt = 0;

  prepareTree();

  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                             EUTELESCOPE::DUMPGEOROOT);
  
  // XXX Assumes the gear file contain the planes ordered in z, anyway
  // do we need the order?
  int plane=0;
  for(unsigned int sensorID : geo::gGeometry().sensorIDsVec()) 
  {
      _indexIDMap[sensorID] = plane;
      plane++;
  }

  for (auto dutID : _DUTIDs) {
    // Later we need to shift the sensor since in EUTel centre of sensor is 0|0
    // while in TBmon(II) it is in the lower left corner
    geo::EUTelGenericPixGeoDescr *geoDescr =
        geo::gGeometry().getPixGeoDescr(dutID);
    double xSize, ySize;
    geoDescr->getSensitiveSize(xSize, ySize);

    _xSensSize[dutID] = xSize;
    _ySensSize[dutID] = ySize;
  }
}

void EUTelAPIXTbTrackTuple::processRunHeader(LCRunHeader *runHeader) {
  auto eutelHeader = std::make_unique<EUTelRunHeaderImpl>(runHeader);
  eutelHeader->addProcessor(type());
  _nRun++;

  // Decode and print out Run Header information - just a check
  _runNr = runHeader->getRunNumber();
}

void EUTelAPIXTbTrackTuple::processEvent(LCEvent *event) {
  _nEvt++;
  _evtNr = event->getEventNumber();
  EUTelEventImpl *euEvent = static_cast<EUTelEventImpl *>(event);

  if (euEvent->getEventType() == kEORE) {
    streamlog_out(DEBUG5) << "EORE found: nothing else to do." << std::endl;
    return;
  }
  // Clear all event info containers
  clear();

  // try to read in hits (e.g. fitted hits in local frame)
  if (!readHits(_inputTrackerHitColName, event)) {
    return;
  }
  
  // try to read in measured hits in the local frame
  if (!readMeasHits(_inputMeasHitColName, event)) {
    return;
  }

  // read in raw data
  if (!readZsHits(_dutZsColName, event)) {
    return;
  }

  // read in tracks
  if (!readTracks(event)) {
    return;
  }

  // fill the trees
  _zstree->Fill();
  _eutracks->Fill();
  _mhits->Fill();
  _euhits->Fill();


  _isFirstEvent = false;
}

void EUTelAPIXTbTrackTuple::end() {
  // write version number
  _versionNo->push_back(1.3);
  _versionTree->Fill();
  // Maybe some stats output?
  _file->Write();
}

// Read in TrackerHit(Impl) to later dump them
bool EUTelAPIXTbTrackTuple::readHits(std::string hitColName, LCEvent *event) {
  LCCollection *hitCollection = nullptr;

  try {
    hitCollection = event->getCollection(hitColName);
  } catch (lcio::DataNotAvailableException &e) {
    streamlog_out(DEBUG2) << "Hit collection " << hitColName
                          << " not found in event " << event->getEventNumber()
                          << "!" << std::endl;
    return false;
  }

  int nHit = hitCollection->getNumberOfElements();
  _nHits = nHit;

  for (int ihit = 0; ihit < nHit; ihit++) {
    TrackerHitImpl *meshit =
        dynamic_cast<TrackerHitImpl *>(hitCollection->getElementAt(ihit));
    const double *pos = meshit->getPosition();
    LCObjectVec clusterVec = (meshit->getRawHits());

    UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder(EUTELESCOPE::HITENCODING);
    int sensorID = hitDecoder(meshit)["sensorID"];

    // Only dump DUT hits
    if (std::find(_DUTIDs.begin(), _DUTIDs.end(), sensorID) == _DUTIDs.end()) {
      continue;
    }

    double x = pos[0];
    double y = pos[1];
    double z = pos[2];

    // offset by half sensor/sensitive size
    _hitXPos->push_back(x + _xSensSize.at(sensorID) / 2.0);
    _hitYPos->push_back(y + _ySensSize.at(sensorID) / 2.0);
    _hitZPos->push_back(z);
    _hitSensorId->push_back(sensorID);
  }

  return true;
}

// Read in TrackerHit(Impl) to later dump them
bool EUTelAPIXTbTrackTuple::readMeasHits(std::string hitColName, LCEvent *event) 
{
    LCCollection *hitCollection = nullptr;
    try 
    {
        hitCollection = event->getCollection(hitColName);
    } catch (lcio::DataNotAvailableException &e) 
    {
        streamlog_out(DEBUG2) << "Hit collection " << hitColName
            << " not found in event " << event->getEventNumber()
            << "!" << std::endl;
        return false;
    }
    
    const int nmHit = hitCollection->getNumberOfElements();

    _nmHits = 0;    
    for(int ihit = 0; ihit < nmHit; ++ihit) 
    {
        TrackerHitImpl *meshit = dynamic_cast<TrackerHitImpl *>(hitCollection->getElementAt(ihit));
        const double *pos = meshit->getPosition();
        
        UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder(EUTELESCOPE::HITENCODING);
        const int sensorID = hitDecoder(meshit)["sensorID"];
        
        // Only dump DUT hits
        if(std::find(_DUTIDs.begin(), _DUTIDs.end(), sensorID) == _DUTIDs.end()) 
        {
            continue;
        }
        ++_nmHits;
        // the Raw data
        const auto raw = meshit->getRawHits();
        if(raw.size() == 0 || raw[0] == nullptr)
        {
            streamlog_out(ERROR) << "Severe issue: hit does" <<
                "not contain its raw data" << std::endl;
            return false;
        }
        // Note that the TrackDataImpl::getChargeValues() function returns a vector which
        // contains [X-coordinate,Y-coordinate, TOT ,time, ??, ?? , ?? ,?? ] per each element of the cluster
        auto & zsdata = static_cast<TrackerDataImpl*>(raw[0])->getChargeValues();
        std::map<std::pair<int,int>,int> cluster;
        for(unsigned int k=0; k < zsdata.size()-7; k+=8)
        {
            //charge += zsdata[k+2];
            cluster[std::pair<int,int>(zsdata[k],zsdata[k+1])] = zsdata[k+2];
        }
        // Adding up charge:
        const int charge = std::accumulate(
                std::begin(cluster),
                std::end(cluster),
                0,
                [](int value,const std::map<std::pair<int,int>,int>::value_type& p)
                {
                    return value+p.second;
                });
        // Flattening the cluster
        std::map<int,int> cluster_x;
        std::map<int,int> cluster_y;
        std::for_each(cluster.begin(),cluster.end(), 
                [&](std::pair<const std::pair<int,int>,int> & el)
                {
                    if(cluster_x.find(el.first.first) != cluster_x.end())
                    {
                        cluster_x[el.first.first] += el.second;
                    }
                    else
                    {
                        cluster_x.insert(std::make_pair(el.first.first,el.second));
                    }
                    if(cluster_y.find(el.first.second) != cluster_y.end())
                    {
                        cluster_y[el.first.second] += el.second;
                    }
                    else
                    {
                        cluster_y.insert(std::make_pair(el.first.second,el.second));
                    }
                });

        // Fill the tree
        // offset by half sensor/sensitive size
        _mhitXpos->push_back(pos[0]+_xSensSize.at(sensorID)/2.0);
        _mhitYpos->push_back(pos[1]+_ySensSize.at(sensorID)/2.0);
        _mhitTOT->push_back(charge);
        _mhitId->push_back(sensorID);
        _mhitBCID->push_back(zsdata[3]);
        _mhitSize->push_back(zsdata.size()/8);
        _mhitSizeX->push_back( cluster_x.size() );
        _mhitEtaX->push_back( static_cast<float>(cluster_x.cbegin()->second)/static_cast<float>(charge) );
        _mhitSizeY->push_back( cluster_y.size() );
        _mhitEtaY->push_back( static_cast<float>(cluster_y.cbegin()->second)/static_cast<float>(charge) );
    }
    return true;
}


EUTelVirtualCluster * EUTelAPIXTbTrackTuple::get_cluster_from_raw_data(EVENT::LCObject * raw, EVENT::LCEvent * evt )
{
    // XXX: HARDCODED!! In order to avoid the use of the cluster collection from
    // the data, the hardcoded encoding string is included here. This is dangerous,
    // as the decoding could change at some point. Maybe it would be better to 
    // include the cluster collection explicitly.. 
    const auto st_encoding("sensorID:7,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5,quality:5");

    EUTelVirtualCluster * cluster =  nullptr;
    UTIL::CellIDDecoder<TrackerPulseImpl> clDec(st_encoding);
    const int tmp = clDec(nullptr)["type"];
    const ClusterType type(static_cast<ClusterType>(tmp));
    if(type == kEUTelDFFClusterImpl) 
    {
        // digital fixed cluster implementation. Remember it can come from
        // both RAW and ZS data
        cluster = new EUTelDFFClusterImpl(static_cast<TrackerDataImpl *>(raw));
    }
    else if(type == kEUTelFFClusterImpl) 
    {
        // fixed cluster implementation. Remember it can come from
        // both RAW and ZS data
        cluster = new EUTelFFClusterImpl(static_cast<TrackerDataImpl *>(raw));
    } 
    else if(type == kEUTelBrickedClusterImpl)
    {
        // bricked cluster implementation
        // Remember it can come from both RAW and ZS data
        cluster = new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *>(raw));
    }
    else if(type == kEUTelSparseClusterImpl) 
    {
        // ok the cluster is of sparse type, but we also need to know
        // the kind of pixel description used. This information is
        // stored in the corresponding cluster data collection.
        LCCollectionVec *sparseClusterCollectionVec =  
            dynamic_cast<LCCollectionVec *>(evt->getCollection(_dutZsColName));
        TrackerDataImpl *oneCluster = dynamic_cast<TrackerDataImpl *>(sparseClusterCollectionVec->getElementAt(0));
        UTIL::CellIDDecoder<TrackerDataImpl> anotherDecoder(sparseClusterCollectionVec);
        SparsePixelType pixelType = static_cast<SparsePixelType>(static_cast<int>(anotherDecoder(oneCluster)["sparsePixelType"]));

        // now we know the pixel type. So we can properly create a new
        // instance of the sparse cluster
        if(pixelType == kEUTelGenericSparsePixel) 
        {
            cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel>(static_cast<TrackerDataImpl *>(raw));
        }
        else
        {
            streamlog_out(ERROR4) << "Unknown pixel type. Quitting..." << std::endl;
            throw UnknownDataTypeException("Pixel type unknown");
        }
    }
    else
    {
        streamlog_out(ERROR4) << "Unknown cluster type. Quitting..." << std::endl;
        throw UnknownDataTypeException("Cluster type unknown");
    }
    
    return cluster;
}


// Read in TrackerHit to later dump
bool EUTelAPIXTbTrackTuple::readTracks(LCEvent *event) {
  LCCollection *trackCol = nullptr;

  try {
    trackCol = event->getCollection(_inputTrackColName);
  } catch (lcio::DataNotAvailableException &e) {
    streamlog_out(DEBUG2) << "Track collection " << _inputTrackColName
                          << " not found in event " << event->getEventNumber()
                          << "!" << std::endl;
    return false;
  }

  // setup cellIdDecoder to decode the hit properties
  UTIL::CellIDDecoder<TrackerHitImpl> hitCellDecoder(EUTELESCOPE::HITENCODING);

  int nTrackParams = 0;
  for (int itrack = 0; itrack < trackCol->getNumberOfElements(); itrack++) {
    lcio::Track *fittrack =
        dynamic_cast<lcio::Track *>(trackCol->getElementAt(itrack));

    std::vector<EVENT::TrackerHit *> trackhits = fittrack->getTrackerHits();
    double chi2 = fittrack->getChi2();
    double ndof = fittrack->getNdf();
    double dxdz = fittrack->getOmega();
    double dydz = fittrack->getPhi();

    // extract the hit pattern. In bits ordered by z ()
    // [ Pln ... Pl2 Pl1 Pl0 ]
    UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder(EUTELESCOPE::HITENCODING);
    int hitpattern = 0;
    for(unsigned int ih=0; ih < trackhits.size(); ++ih)
    {
        TrackerHitImpl * meashit = dynamic_cast<TrackerHitImpl*>(trackhits[ih]);
        if( (hitCellDecoder(meashit)["properties"] & kFittedHit) != 0)
        {
            // fitted hits, skip it
            continue;
        }
        const int plane =  _indexIDMap[hitDecoder(meashit)["sensorID"]];
        //const int hitmask = (1 << plane);
        // Activate the corresponding bit to that plane
        hitpattern |= (1 << plane); //hitmask;
    }

    /* Get the (fitted) hits belonging to this track,
       they are in global frame when coming from the straight track fitter */
    for (unsigned int ihit = 0; ihit < trackhits.size(); ihit++) {
      TrackerHitImpl *fittedHit =
          dynamic_cast<TrackerHitImpl *>(trackhits.at(ihit));
      const double *pos = fittedHit->getPosition();
      if ((hitCellDecoder(fittedHit)["properties"] & kFittedHit) == 0) {
        continue;
      }

      //UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder(EUTELESCOPE::HITENCODING);
      int sensorID = hitDecoder(fittedHit)["sensorID"];

      // Dump the (fitted) hits for the DUTs
      if (std::find(_DUTIDs.begin(), _DUTIDs.end(), sensorID) ==
          _DUTIDs.end()) {
        continue;
      }

      nTrackParams++;
      /* Transform to local coordinates */
      double pos_loc[3];
      geo::gGeometry().master2Local(sensorID, pos, pos_loc);

      double x = pos_loc[0];
      double y = pos_loc[1];

      // eutrack tree
      _xPos->push_back(x);
      _yPos->push_back(y);
      _dxdz->push_back(dxdz);
      _dydz->push_back(dydz);
      _trackIden->push_back(sensorID);
      _trackNum->push_back(itrack);
      _chi2->push_back(chi2);
      _ndof->push_back(ndof);
      _hitpattern->push_back(hitpattern);
    }
  }

  _nTrackParams = nTrackParams;
  return true;
}

// Read in raw (zs) TrackerData(Impl) to later dump
bool EUTelAPIXTbTrackTuple::readZsHits(std::string colName, LCEvent *event) {
  LCCollectionVec *zsInputCollectionVec = nullptr;

  try {
    zsInputCollectionVec =
        dynamic_cast<LCCollectionVec *>(event->getCollection(colName));
  } catch (DataNotAvailableException &e) {
    streamlog_out(DEBUG2) << "Raw ZS data collection " << colName
                          << " not found in event " << event->getEventNumber()
                          << "!" << std::endl;
    return false;
  }

  UTIL::CellIDDecoder<TrackerDataImpl> cellDecoder(zsInputCollectionVec);
  for (unsigned int plane = 0; plane < zsInputCollectionVec->size(); plane++) {
    TrackerDataImpl *zsData = dynamic_cast<TrackerDataImpl *>(
        zsInputCollectionVec->getElementAt(plane));
    SparsePixelType type = static_cast<SparsePixelType>(
        static_cast<int>(cellDecoder(zsData)["sparsePixelType"]));
    int sensorID = cellDecoder(zsData)["sensorID"];

    if (type == kEUTelGenericSparsePixel) {
      auto sparseData = std::make_unique<
          EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(zsData);

      for (auto &apixPixel : *sparseData) {
        _nPixHits++;
        p_iden->push_back(sensorID);
        p_row->push_back(apixPixel.getYCoord());
        p_col->push_back(apixPixel.getXCoord());
        p_tot->push_back(static_cast<int>(apixPixel.getSignal()));
        p_lv1->push_back(static_cast<int>(apixPixel.getTime()));
      }

    } else if (type == kEUTelMuPixel) {
      auto sparseData =
          std::make_unique<EUTelTrackerDataInterfacerImpl<EUTelMuPixel>>(
              zsData);
      for (auto &binaryPixel : *sparseData) {
        _nPixHits++;
        p_iden->push_back(sensorID);
        p_row->push_back(binaryPixel.getYCoord());
        p_col->push_back(binaryPixel.getXCoord());
        p_hitTime->push_back(binaryPixel.getHitTime());
        p_frameTime->push_back(binaryPixel.getFrameTime());
      }
    } else {
      throw UnknownDataTypeException("Unknown sparsified pixel");
    }
  }
  return true;
}

void EUTelAPIXTbTrackTuple::clear() {
  /* Clear zsdata */
  p_col->clear();
  p_row->clear();
  p_tot->clear();
  p_iden->clear();
  p_lv1->clear();
  p_hitTime->clear();
  p_frameTime->clear();
  _nPixHits = 0;
  
  // Clear measured hits
  _nmHits=0;
  _mhitXpos->clear();
  _mhitYpos->clear(); 
  _mhitEtaX->clear();
  _mhitEtaY->clear();
  _mhitTOT->clear();
  _mhitId->clear();
  _mhitBCID->clear();
  _mhitSize->clear();
  _mhitSizeX->clear();
  _mhitSizeY->clear();

  /* Clear hittrack */
  _xPos->clear();
  _yPos->clear();
  _dxdz->clear();
  _dydz->clear();
  _trackNum->clear();
  _trackIden->clear();
  _chi2->clear();
  _ndof->clear();
  _hitpattern->clear();
  // Clear hits
  _hitXPos->clear();
  _hitYPos->clear();
  _hitZPos->clear();
  _hitSensorId->clear();
}

void EUTelAPIXTbTrackTuple::prepareTree() {
  _file = new TFile(_path2file.c_str(), "RECREATE");

  _xPos = new std::vector<double>();
  _yPos = new std::vector<double>();
  _dxdz = new std::vector<double>();
  _dydz = new std::vector<double>();
  _trackIden = new std::vector<int>();
  _trackNum = new std::vector<int>();
  _chi2 = new std::vector<double>();
  _ndof = new std::vector<double>();
  _hitpattern = new std::vector<int>();
  
  _mhitXpos = new std::vector<double>();
  _mhitYpos = new std::vector<double>();
  _mhitEtaX = new std::vector<double>();
  _mhitEtaY = new std::vector<double>();
  _mhitTOT  = new std::vector<int>();
  _mhitId   = new std::vector<int>();
  _mhitBCID = new std::vector<int>();
  _mhitSize = new std::vector<int>();
  _mhitSizeX= new std::vector<int>();
  _mhitSizeY= new std::vector<int>();

  p_col = new std::vector<int>();
  p_row = new std::vector<int>();
  p_tot = new std::vector<int>();
  p_iden = new std::vector<int>();
  p_lv1 = new std::vector<int>();
  p_hitTime = new std::vector<int>();
  p_frameTime = new std::vector<double>();

  _hitXPos = new std::vector<double>();
  _hitYPos = new std::vector<double>();
  _hitZPos = new std::vector<double>();
  _hitSensorId = new std::vector<int>();

  _versionNo = new std::vector<double>();
  _versionTree = new TTree("version", "version");
  _versionTree->Branch("no", &_versionNo);

  _euhits = new TTree("fitpoints", "fitpoints");
  _euhits->SetAutoSave(1000000000);

  _euhits->Branch("nHits", &_nHits);
  _euhits->Branch("xPos", &_hitXPos);
  _euhits->Branch("yPos", &_hitYPos);
  _euhits->Branch("zPos", &_hitZPos);
  _euhits->Branch("sensorId", &_hitSensorId);

  _zstree = new TTree("rawdata", "rawdata");
  _zstree->SetAutoSave(1000000000);
  _zstree->Branch("nPixHits", &_nPixHits);
  _zstree->Branch("euEvt", &_nEvt);
  _zstree->Branch("col", &p_col);
  _zstree->Branch("row", &p_row);
  _zstree->Branch("tot", &p_tot);
  _zstree->Branch("lv1", &p_lv1);
  _zstree->Branch("iden", &p_iden);
  _zstree->Branch("hitTime", &p_hitTime);
  _zstree->Branch("frameTime", &p_frameTime);

  _mhits = new TTree("meashits", "meashits");
  _mhits->SetAutoSave(1000000000);

  _mhits->Branch("nHits", &_nmHits);
  _mhits->Branch("xPos", &_mhitXpos);
  _mhits->Branch("yPos", &_mhitYpos);
  _mhits->Branch("eta_x", &_mhitEtaX);
  _mhits->Branch("eta_y", &_mhitEtaY);
  _mhits->Branch("tot", &_mhitTOT);
  _mhits->Branch("sensorId", &_mhitId);
  _mhits->Branch("bcid", &_mhitBCID);
  _mhits->Branch("clusterSize", &_mhitSize);
  _mhits->Branch("clusterSize_x", &_mhitSizeX);
  _mhits->Branch("clusterSize_y", &_mhitSizeY);

  // Tree for storing all track param info
  _eutracks = new TTree("tracks", "tracks");
  _eutracks->SetAutoSave(1000000000);
  _eutracks->Branch("nTrackParams", &_nTrackParams);
  _eutracks->Branch("euEvt", &_nEvt);
  _eutracks->Branch("xPos", &_xPos);
  _eutracks->Branch("yPos", &_yPos);
  _eutracks->Branch("dxdz", &_dxdz);
  _eutracks->Branch("dydz", &_dydz);
  _eutracks->Branch("trackNum", &_trackNum);
  _eutracks->Branch("iden", &_trackIden);
  _eutracks->Branch("chi2", &_chi2);
  _eutracks->Branch("ndof", &_ndof);
  _eutracks->Branch("hitPattern",&_hitpattern);

  _euhits->AddFriend(_zstree);
  _euhits->AddFriend(_eutracks);
  _euhits->AddFriend(_mhits);
}
