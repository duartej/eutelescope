/*
 * Created by Eda Yildirim
 *  (2015 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


// ROOT includes:

// eutelescope includes ".h"
#include "EUTelMissingCoordinateEstimatorUseTracks.h"
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

// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>


// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// ROOT
#include "TH2F.h"

// system includes <>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
//#include <memory>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace marlin;
using namespace gear;
using namespace eutelescope;

namespace auxfunc
{
    // Extracted from http://en.cppreference.com/w/cpp/language/constexpr
    constexpr int factorial(int n)
    {
        return n <= 1? 1 : (n * factorial(n - 1));
    }
    
    // Just the mean of a floats
    float get_mean(const std::vector<float> & v)
    {
        if(v.size() == 0)
        {   
            return 0.0;
        }
        return std::accumulate(v.begin(),v.end(),0.0)/float(v.size());
    }

    // And the standard deviation
    float get_std_dev(const std::vector<float> & v,const float & mean)
    {
        if(v.size()==0)
        {
            return 0.0;
        }
        // standard deviation = sqrt( E[(x-E[x])^2] )
        std::vector<float> diff(v.size());
        // Get x-E[x] in `diff`
        std::transform(v.begin(),v.end(),diff.begin(),
                [mean](const float & element) { return element-mean; } );
        // obtain the (x-E[x])^2 by using diff*diff (inner_product)
        const float sq_sum = std::inner_product(diff.begin(),diff.end(),diff.begin(),0.0);
        return std::sqrt(sq_sum/static_cast<float>(diff.size()));
    }
}

// definition of static members, mainly the histonames
std::string EUTelMissingCoordinateEstimatorUseTracks::_hitplotname   = "matchedhits";
std::string EUTelMissingCoordinateEstimatorUseTracks::_missinghits   = "missinghits";
std::string EUTelMissingCoordinateEstimatorUseTracks::_duthitNtracks  = "Nduthits_Ntracks";

// definition of static members mainly used to name histograms
EUTelMissingCoordinateEstimatorUseTracks::EUTelMissingCoordinateEstimatorUseTracks():
    Processor("EUTelMissingCoordinateEstimatorUseTracks"),
    _inputHitCollectionName(),
    _inputTrackCollectionName(),
    _outputHitCollectionName(),
    _dutPlanes(),
    _missingCoordinate(),
    _maxResidual_prov(),
    _maxResidual(),
    _telescopeResolution(0),
    _iRun(0),
    _iEvt(0),
    _missingHitPos(0),
    _knownHitPos(0),
    _nDutHits(0),
    _nDutHitsCreated(0),
    _maxExpectedCreatedHitPerDUTHit(10),
    _numberOfCreatedHitPerDUTHit()
{
    // modify processor description
    _description =  "EUTelMissingCoordinateEstimatorUseTracks As the name suggest this processor"\
                     " finds the position of the missing coordinate on a DUT micro-strip detector "\
                     "using the reconstructed tracks from the telescope. Previsously, the DUTs "\
                     "must be pre-aligned with the telescope";
    registerInputCollection(LCIO::TRACKERHIT,"InputHitCollectionName", "Input hit collection name."\
            " Hits should be in global coordinates and pre-aligned", _inputHitCollectionName, std::string(""));
    
    registerInputCollection(LCIO::TRACK,"InputTrackCollectionName", "Input track collection name to be used for the "\
            "DUT missing coordinate search", _inputTrackCollectionName, std::string(""));
    
    registerOutputCollection(LCIO::TRACKERHIT,"OutputHitCollectionName", "Output hit collection name", _outputHitCollectionName, std::string(""));    
    
    registerProcessorParameter("DUTPlanes","This is the list of sensorIDs that missing "\
            "coordinate of their hits needs to be found. Notice that if the specified "\
            "coordinate already exists it will be overwritten", _dutPlanes, EVENT::IntVec());
    
    registerProcessorParameter("MissingCoordinate","The coordinate axis that needs to be estimated [X or Y]",
            _missingCoordinate, std::string("X") );
    
    registerProcessorParameter("MaxResidual","This processor will look for a closest track "\
            "(in known coordinate) to determine if the hits are correlated. The hits will be "\
            "considered as correlated if the residual is smaller than MaxResidual", 
            _maxResidual_prov, EVENT::FloatVec({0.05}));   
    
    registerProcessorParameter("TelescopeResolution","Assumed telescope resolution. It is used"\
            " to define the binning size (and eventually the uncertainty in the missing coordinate",
            _telescopeResolution, float(0.05) );   
}


void EUTelMissingCoordinateEstimatorUseTracks::init() 
{
    // Check consistency between residuals and duts
    if(_dutPlanes.size() != _maxResidual_prov.size())
    {
        streamlog_out(WARNING2) << "Misleading configuration, number of DUTs planes must " 
            << "be equal to the number of set residuals. Forcing this behaviour" << std::endl;
        const float residual = _maxResidual_prov[0];
        _maxResidual_prov.clear();
        for(unsigned int k=0; k < _dutPlanes.size(); ++k)
        {
            _maxResidual_prov.push_back(residual);
        }
    }

    // Initialize the residual map
    for(unsigned int k=0; k < _dutPlanes.size(); ++k)
    {
        _maxResidual[_dutPlanes[k]] = _maxResidual_prov[k];
    }

    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
    
    // set to zero the run and event counters
    _iRun = 0;
    _iEvt = 0;
    
    // check if _missingCoordinate is valid
    // if it is given as lowercase make them uppercase letters
    // and set missing and known hit position variables
    std::transform(_missingCoordinate.begin(), _missingCoordinate.end(), _missingCoordinate.begin(), ::toupper);

    if(_missingCoordinate == std::string("X"))
    {
        _missingHitPos = 0;
        _knownHitPos = 1;
    }
    if(_missingCoordinate == std::string("Y")) 
    {
        _missingHitPos = 1;
        _knownHitPos = 0;
    }
    
    if(_missingCoordinate == std::string("X") || _missingCoordinate == std::string("Y") ) 
    {
        streamlog_out (DEBUG4) << "MissingCoordinate value set as: "<< _missingCoordinate << std::endl;
    }
    else 
    {
        streamlog_out (ERROR4) << "MissingCoordinate value ("<<_missingCoordinate<<") is not valid!"<< std::endl;
        exit(-1);
    }
    
    // set counters to zero
    _nDutHits = 0;
    _nDutHitsCreated = 0;
    
    for(unsigned int i=0; i < _maxExpectedCreatedHitPerDUTHit+1; i++)
    {
        _numberOfCreatedHitPerDUTHit.push_back(0);
    }

    // And initialize histograms
    bookHistos();
}


void EUTelMissingCoordinateEstimatorUseTracks::processRunHeader (LCRunHeader * rdr) 
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

void EUTelMissingCoordinateEstimatorUseTracks::processEvent (LCEvent * event) 
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
    LCCollectionVec * outputHitCollection= nullptr;
    try
    {
        inputHitCollection = static_cast<LCCollectionVec*>(event->getCollection( _inputHitCollectionName ));
    }
    catch(DataNotAvailableException& e)
    {
        streamlog_out  ( MESSAGE2 ) <<  "No input collection " << _inputHitCollectionName 
            << " found on event " << event->getEventNumber()
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
    
    
    try
    {
        outputHitCollection  = static_cast<LCCollectionVec*> (event->getCollection( _outputHitCollectionName ));
    }
    catch(...)
    {
        outputHitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
    }
    // prepare an encoder for the hit collection
    CellIDEncoder<TrackerHitImpl> outputCellIDEncoder(EUTELESCOPE::HITENCODING, outputHitCollection);    
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
            ++_nDutHits;
        }
        else
        {
            // Store all telescope hits in the new collection, 
            // we will store DUT hits after updating its position
            outputHitCollection->push_back( cloneHit(inputHit) );
        }
    }
 
    // with countCreatedDutHits vector we will count how many hits we create out of one DUT hit
    std::vector<unsigned int> countCreatedDutHits(dutHits.size(),0);
    
    /*
     The line that passes through 2 points can be written as L(t)= P1 + V*t
     where V is the displacement vector and P1 is the starting point
     so L(t) becomes: L(t) = (x1,y1,z1) + (x2-x1, y2-y1, z2-z1)*t
     * x=x1+(x2−x1)t
     * y=y1+(y2−y1)t
     * z=z1+(z2−z1)t
     */

    // NOTE FROM EUTelAPIXTbTrackTuple.cc --> dxdz = track->getOmega()
    //                                        dydz = track->getPhi()
    // See 191 at src/EUTelDafFitter.cc !!!
 
    CellIDDecoder<TrackerHit> trkHitDecoder(EUTELESCOPE::HITENCODING);
    
    // Loop over the tracks to extract the lines 
    // the two points used to define a line are stored here
    // Reference point,direction vector (unitary) 
    std::vector<std::pair<LCVector3D,LCVector3D> > fitHits;
    fitHits.reserve(inputTrackCollection->getNumberOfElements());
    for(int itrk = 0; itrk < inputTrackCollection->getNumberOfElements(); ++itrk)
    {
        // Any point it is obtained as: P(t) = RefPoint + (dxdz,dydz,1.0)*t
        TrackImpl * track = dynamic_cast<TrackImpl*>(inputTrackCollection->getElementAt(itrk));
        // See L.191 at src/EUTelDafFitter.cc 
        const float dxdz = track->getOmega();
        const float dydz = track->getPhi();
        fitHits.push_back(std::make_pair(LCVector3D(track->getReferencePoint()[0],track->getReferencePoint()[1],track->getReferencePoint()[2]),
                    LCVector3D(dxdz,dydz,1.0)));
    }

    // Preparing to fill some histos: the sensor id of the duts and the 
    // number of hits found for that dut
    std::map<int,int> dutSensorId_Nhits;
    // Loop over the DUTs
    int iDutHit = -1;
    for(auto & dutHit: dutHits)
    {
        ++iDutHit;
        const double * dutHitPos = dutHit->getPosition();
        // Extract the sensorid
        const int sensorID_dut = trkHitDecoder(dutHit)["sensorID"];
        dutSensorId_Nhits[sensorID_dut]++;

        for(auto & ref_direction: fitHits)
        {
            // Use the know Z-position of the DUT to find the path-length
            // t = (z-z1)/(z2-z1)
            double t = (dutHitPos[2] - ref_direction.first.z())/ref_direction.second.z();

            // Find the known coordinate coresponding to that z on the line
            double knownPredictedPos = ref_direction.first[_knownHitPos]+ref_direction.second[_knownHitPos]*t;
            
            // Check if is close enough (XXX use track covariance matrix, ...)
            if( std::fabs(knownPredictedPos-dutHitPos[_knownHitPos]) < _maxResidual[sensorID_dut])
            {
                // Extract the estimate position in the unknown coordinate
                double estimatedPredictedPos = ref_direction.first[_missingHitPos]+ref_direction.second[_missingHitPos]*t;
                // And copy the new hits, first copy the old positions
                double newhitpos[3] = { dutHit->getPosition()[0], dutHit->getPosition()[1], dutHit->getPosition()[2] };
                // and update with the just estimated unknown coordinate
                newhitpos[_missingHitPos]=estimatedPredictedPos;

                // Copy the hit and afterwards set the new position
                TrackerHitImpl * newHit = cloneHit(dutHit);
                newHit->setPosition( &(*newhitpos) );
                // -- XXX: It could be possible to check if this one is the closest one, by
                //        creating a map of this hit with the distance (in the if), and use
                //        always the closest
                outputHitCollection->push_back(newHit);
                // Histograming ...
                const std::string histoname(_hitplotname+"_"+std::to_string(sensorID_dut));
                dynamic_cast<TH2F*>(_histoMap[histoname])->Fill(newhitpos[0],newhitpos[1]);
                
                ++_nDutHitsCreated;
                ++countCreatedDutHits[iDutHit];
           }
           else
           {
               // Histograming ...
               const std::string histoname(_missinghits+"_"+std::to_string(sensorID_dut));
               dynamic_cast<TH2F*>(_histoMap[histoname])->Fill(dutHitPos[_knownHitPos],(dutHitPos[_knownHitPos]-knownPredictedPos)/_maxResidual[sensorID_dut]);
           }
        }
    }

    // Histograming number of elements per event
    for(auto & id_n: dutSensorId_Nhits)
    {
        const std::string histoname_ntrk(_duthitNtracks+"_"+std::to_string(id_n.first));
        dynamic_cast<TH2F*>(_histoMap[histoname_ntrk])->Fill(inputTrackCollection->getNumberOfElements(),id_n.second);
    }

    for(unsigned int iDutHit=0; iDutHit< dutHits.size(); iDutHit++)
    {
        if(_maxExpectedCreatedHitPerDUTHit < countCreatedDutHits[iDutHit])
        {
            ++_numberOfCreatedHitPerDUTHit[_maxExpectedCreatedHitPerDUTHit];
        }
	else 
        {
            ++_numberOfCreatedHitPerDUTHit[countCreatedDutHits[iDutHit]];
        }
    }
    
    try
    {
        event->getCollection( _outputHitCollectionName ) ;
    }
    catch(...)
    {
        event->addCollection( outputHitCollection, _outputHitCollectionName );
    }
    
    if( isFirstEvent() )
    {
        _isFirstEvent = false;
    }
}


TrackerHitImpl* EUTelMissingCoordinateEstimatorUseTracks::cloneHit(TrackerHitImpl *inputHit)
{
    TrackerHitImpl * newHit = new TrackerHitImpl;

    // copy hit position
    const double* hitPos = inputHit->getPosition();
    newHit->setPosition( &hitPos[0] );
    
    // copy cov. matrix
    newHit->setCovMatrix( inputHit->getCovMatrix() );
    
    // copy type
    newHit->setType( inputHit->getType() );
    
    // copy rawhits
    LCObjectVec clusterVec = inputHit->getRawHits();
    newHit->rawHits() = clusterVec;
   
    // copy cell IDs
    newHit->setCellID0( inputHit->getCellID0() );
    newHit->setCellID1( inputHit->getCellID1() );
    
    // copy EDep
    newHit->setEDep( inputHit->getEDep() );

    // copy EDepError
    newHit->setEDepError( inputHit->getEDepError() );
    
    // copy Time
    newHit->setTime( inputHit->getTime() );
    
    // copy Quality
    newHit->setQuality( inputHit->getQuality() );

    return newHit;
}


void EUTelMissingCoordinateEstimatorUseTracks::end()
{
    streamlog_out ( MESSAGE4 )  << "Number of hits you had from all DUTs "<< _nDutHits << std::endl;
    streamlog_out ( MESSAGE4 )  << "Number of hits created with the estimated missing coordinate "<< _nDutHitsCreated << std::endl;
    for (unsigned int i=0; i<_numberOfCreatedHitPerDUTHit.size(); i++)
    {
        streamlog_out ( MESSAGE4 )  << "You created "<< i ;
        if(i==_maxExpectedCreatedHitPerDUTHit) 
        {
            streamlog_out ( MESSAGE4 )  << " or more";
        }
	streamlog_out ( MESSAGE4 )  << " hits per DUT hit "<< _numberOfCreatedHitPerDUTHit[i] <<" many times"<<std::endl;
    }
    streamlog_out ( MESSAGE4 )  << "Successfully finished" << std::endl;
    
}

void EUTelMissingCoordinateEstimatorUseTracks::bookHistos() 
{
    // create the directory in the root file
    AIDAProcessor::tree(this)->cd(this->name());

    // And the files
    for(auto & sensorID: _dutPlanes)
    {
        const std::string histoname(_hitplotname+"_"+std::to_string(sensorID));
        const std::string histoname_missing(_missinghits+"_"+std::to_string(sensorID));
        const std::string histoname_NtrkHit(_duthitNtracks+"_"+std::to_string(sensorID));
        // Get the sensor layouta
        //const double constant=1.2;
        const double xMin =  -(geo::gGeometry().siPlaneXSize ( sensorID )/2);
        const double xMax = ( geo::gGeometry().siPlaneXSize ( sensorID )/2);   
        
        const double yMin = -(geo::gGeometry().siPlaneYSize ( sensorID )/2);
        const double yMax = (geo::gGeometry().siPlaneYSize ( sensorID )/2); 
        
        const int xNBin = geo::gGeometry().siPlaneXNpixels(sensorID);
        const int yNBin = static_cast<int>(std::fabs(yMax-yMin)/_telescopeResolution);

        _histoMap.insert(std::make_pair(histoname,new TH2F(histoname.c_str(),";x [mm]; y[mm]",xNBin,xMin,xMax,yNBin,yMin,yMax)));
        _histoMap.insert(std::make_pair(histoname_missing,new TH2F(histoname_missing.c_str(),";x [mm]; y[mm]",xNBin,xMin,xMax,yNBin,xMin,xMax)));
        _histoMap.insert(std::make_pair(histoname_NtrkHit,new TH2F(histoname_NtrkHit.c_str(),";N_{trk}; N_{hit}^{DUT}",50,-0.5,49.5,10,-0.5,9.5)));
        // Note to change the resolution if it's lower than the pitch of the DUT
        if( geo::gGeometry().siPlaneXPitch(sensorID) > _maxResidual[sensorID] )
        {
            streamlog_out(WARNING2) << "Force the residual to be at least the pitch size: " 
                <<  geo::gGeometry().siPlaneXPitch(sensorID) << " [mm]" << std::endl;
            _maxResidual[sensorID] =  geo::gGeometry().siPlaneXPitch(sensorID);
        }
    }
}
