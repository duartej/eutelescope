// eutelescope includes ".h"
#include "EUTelPreAlignment.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// system includes <>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <memory>
#include <cstdio>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

namespace auxfunc
{
    static UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder (EUTELESCOPE::HITENCODING);

    bool isReferenceDetectorHit(const TrackerHitImpl * trkH, const int & refPlaneID)
    {
        return hitDecoder(trkH)["sensorID"] == refPlaneID;
    }
}


EUTelPreAlign::EUTelPreAlign(): Processor("EUTelPreAlign")
{
    _description = "Apply alignment constants to hit collection";

    registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName", 
            "The name of the input hit collection", 
            _inputHitCollectionName, std::string("hit"));
    
    registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedID, 0);

    registerOptionalParameter("AlignmentConstantLCIOFile",
            "Name of LCIO db file where alignment constantds will be stored", 
            _alignmentConstantLCIOFile, std::string("alignment.slcio") );

    registerOptionalParameter("HotPixelCollectionName", 
            "This is the name of the hot pixel collection that clusters should be checked against (optional).", 
            _hotPixelCollectionName, std::string(""));

    registerProcessorParameter ("Events", 
            "How many events should be used for an approximation to the X,Y shifts (pre-alignment)? (default=50000)", 
            _events, 50000 );
 
    registerOptionalParameter("ResidualsXMin",
            "Minimal values of the hit residuals in the X direction for a correlation band. "\
            "Note: these numbers are ordered according to the z position of the sensors and "\
            "NOT according to the sensor id.",
            _residualsXMin, std::vector<float > (6, -10.) );

    registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y "\
            "direction for a correlation band. Note: these numbers are ordered according to "\
            "the z position of the sensors and NOT according to the sensor id.",
            _residualsYMin, std::vector<float > (6, -10.) );

    registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X "\
            "direction for a correlation band. Note: these numbers are ordered according to "\
            "the z position of the sensors and NOT according to the sensor id.",
            _residualsXMax, std::vector<float > (6,  10.) );

    registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y "\
            "direction for a correlation band. Note: these numbers are ordered according to "\
            "the z position of the sensors and NOT according to the sensor id.",
            _residualsYMax, std::vector<float > (6,  10.) );

    registerOptionalParameter ("MinNumberOfCorrelatedHits", 
            "If there are more then this number of correlated hits (planes->track candidate) (default=5)",
            _minNumberOfCorrelatedHits, static_cast <int> (5) );

    registerOptionalParameter("HistogramFilling", 
            "Switch on or off the histogram filling", 
            _fillHistos, bool(true) );
  
    registerOptionalParameter("DumpGEAR", 
            "Dump alignment into GEAR file instead of prealignment database", 
            _dumpGEAR, bool(false) );
  
    registerOptionalParameter("NewGEARSuffix", 
            "Suffix for the new GEAR file, set to empty string (this is not default!) to overwrite old GEAR file", 
            _GEARFileSuffix, std::string("_pre") );

    registerOptionalParameter("ExcludedPlanes", 
            "The list of sensor IDs that shall be excluded.", 
            _ExcludedPlanes, std::vector<int>() );

    registerOptionalParameter("ExcludedPlanesXCoord", 
            "The list of sensor IDs for which the X coordinate shall be excluded.", 
            _ExcludedPlanesXCoord, std::vector<int>() );

    registerOptionalParameter("ExcludedPlanesYCoord", 
            "The list of sensor IDs for which the Y coordinate  shall be excluded.", 
            _ExcludedPlanesYCoord, std::vector<int>() );
}

void EUTelPreAlign::init()
{
    _iRun = 0;  _iEvt = 0;

    _sensorIDVec = geo::gGeometry().sensorIDsVec();
    _sensorIDtoZOrderMap.clear();
    for(size_t index = 0; index < _sensorIDVec.size(); index++) 
    {
        _sensorIDtoZOrderMap.insert( std::make_pair(_sensorIDVec.at(index), (int)index) );
    }
    
    for(std::vector<int>::iterator it = _sensorIDVec.begin(); it != _sensorIDVec.end(); it++) 
    {
        int sensorID = *it;
        if(sensorID == _fixedID) 
        {
            _fixedZ = geo::gGeometry().siPlaneZPosition(sensorID); 
        } 
        else
        {
            _preAligners.push_back( PreAligner(	geo::gGeometry().siPlaneXPitch(sensorID)/10.,
                        geo::gGeometry().siPlaneYPitch(sensorID)/10.,
                        geo::gGeometry().siPlaneZPosition(sensorID),
                        sensorID ) );	
        }
    }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    std::string tempHistoName = "";
    std::string basePath; 

    if( _fillHistos ) 
    {
        // Allow any plane to be the fixed reference:
	for(size_t i = 0; i < _sensorIDVec.size(); i++) 
        {
            int sensorID = _sensorIDVec.at(i);
            
            basePath = "plane_" + to_string( sensorID );
            AIDAProcessor::tree(this)->mkdir(basePath.c_str());
            basePath.append("/");
            
            tempHistoName = "hitXCorr_fixed_to_" + to_string( sensorID ) ;
            AIDA::IHistogram1D * histo1Da = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.);
            _hitXCorr.insert( make_pair( sensorID, histo1Da) );

	    tempHistoName = "hitYCorr_fixed_to_" + to_string( sensorID) ;
	    AIDA::IHistogram1D * histo1Db = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.) ;
	    _hitYCorr.insert( make_pair( sensorID, histo1Db) );
	    
            // 2-dimensional plots using the hits of the reference sensor matched 
            // with the current plane, note that for microstrips...
            const double constant=1.2;
            const double xMin =  -(geo::gGeometry().siPlaneXSize ( sensorID )/2)*constant;
            const double xMax = ( geo::gGeometry().siPlaneXSize ( sensorID )/2)*constant;   

            const double yMin = -(geo::gGeometry().siPlaneYSize ( sensorID )/2)*constant;
            const double yMax = (geo::gGeometry().siPlaneYSize ( sensorID )/2)+constant; 

            const int xNBin =    geo::gGeometry().siPlaneXNpixels ( sensorID );
            const int yNBin =    geo::gGeometry().siPlaneYNpixels ( sensorID );
            tempHistoName = "H2D_hitCorr_fixed_to_" + to_string( sensorID) ;
	    AIDA::IHistogram2D * histo2Db = AIDAProcessor::histogramFactory(this)->createHistogram2D( 
                    (basePath + tempHistoName).c_str(), xNBin, xMin, xMax, yNBin, yMin, yMax ) ;
	    _2dHitHistos.insert( make_pair( sensorID, histo2Db) );
        }
    }
#endif
    
    // this method is called only once even when the rewind is active
    printParameters ();
}

void EUTelPreAlign::processRunHeader(LCRunHeader* rdr) 
{
    std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
    runHeader->addProcessor( type() );
    ++_iRun;
}


void  EUTelPreAlign::FillHotPixelMap(LCEvent *event)
{
    if( _hotPixelCollectionName.empty()) return;
    
    LCCollectionVec* hotPixelCollectionVec = nullptr;
    try 
    {
        hotPixelCollectionVec = static_cast< LCCollectionVec* > ( event->getCollection( _hotPixelCollectionName  ) );
        streamlog_out ( DEBUG5 ) << "Hotpixel database " << _hotPixelCollectionName.c_str() << " found" << endl; 
    }
    catch (...)
    {
        streamlog_out ( WARNING5 ) << "Hotpixel database " << _hotPixelCollectionName.c_str() << " not found" << endl; 
        return;
    }
    
    CellIDDecoder<TrackerDataImpl> cellDecoder( hotPixelCollectionVec );
    for(int i=0; i<  hotPixelCollectionVec->getNumberOfElements(); i++)
    {
        TrackerDataImpl* hotPixelData = dynamic_cast< TrackerDataImpl *> ( hotPixelCollectionVec->getElementAt( i ) );
        SparsePixelType  type         = static_cast<SparsePixelType> (static_cast<int> (cellDecoder( hotPixelData )["sparsePixelType"]));
        
        int sensorID = static_cast<int>( cellDecoder( hotPixelData )["sensorID"] );

        if( type  ==  kEUTelGenericSparsePixel )
        {  
            auto m26Data = std::make_unique<EUTelSparseClusterImpl<EUTelGenericSparsePixel>>(hotPixelData);

	    std::vector<EUTelGenericSparsePixel*> m26PixelVec;
	    EUTelGenericSparsePixel m26Pixel;
	    //Push all single Pixels of one plane in the m26PixelVec

	    for(unsigned int iPixel = 0; iPixel < m26Data->size(); iPixel++ ) 
            {
                std::vector<int> m26ColVec();
                m26Data->getSparsePixelAt( iPixel, &m26Pixel);

                try
                {
                    _hotPixelMap[sensorID].push_back(std::make_pair(m26Pixel.getXCoord(), m26Pixel.getYCoord()));
                }
                catch(...)
                {
                    streamlog_out ( ERROR5 ) << " cannot add pixel to hotpixel map! SensorID: "  
                        << sensorID << ", X:" << m26Pixel.getXCoord() << ", Y:" << m26Pixel.getYCoord() << endl; 
                    abort();
                }
            }
	}          
    }
}

void EUTelPreAlign::processEvent(LCEvent* event)
{
    if( isFirstEvent()) 
    {
        FillHotPixelMap(event);
	_isFirstEvent = false;
    }
    ++_iEvt;

    if(_iEvt > _events) 
    {
        return;
    }
    
    EUTelEventImpl* evt = static_cast<EUTelEventImpl*> (event);
    
    if( evt->getEventType() == kEORE ) 
    {
        streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
        return;
    }
    else if(evt->getEventType() == kUNKNOWN ) 
    {
        streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
            << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }
    
    LCCollectionVec * inputCollectionVec = nullptr;
    try
    {
        inputCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    }
    catch( DataNotAvailableException& e) 
    {
        streamlog_out(WARNING2) <<  "No input collection " 
            << _inputHitCollectionName << " found on event " << event->getEventNumber()
            << " in run " << event->getRunNumber() << endl;
        return;
    }
        
    std::vector<float> residX;
    std::vector<float> residY;
    std::vector<float> currentX;
    std::vector<float> currentY;
    std::vector<PreAligner*> prealign;

        
    // 1. Extract the distance of each hit in the reference detector 
    // against all the hits in any other detector. Only distances
    // lower than the requested by the user are stored
    //
    // -- Loop over  all the hits in the fixed plane:
    for(auto * rawHit : *inputCollectionVec)
    {
        TrackerHitImpl * refHit = dynamic_cast<TrackerHitImpl*>(rawHit);
        if( ! auxfunc::isReferenceDetectorHit(refHit,_fixedID) )
        {
            continue;
        }

        const double* refPos = refHit->getPosition();
        residX.clear();
        residY.clear();
        currentX.clear();
        currentY.clear();
        prealign.clear();
        
        // loop over all the available hits (but the ones coming from
        // the reference detector)
        for(auto * irawHit : *inputCollectionVec)
        {
            TrackerHitImpl * hit = dynamic_cast<TrackerHitImpl*>(irawHit);
            // Must be not from the reference detector
            if(auxfunc::isReferenceDetectorHit(hit,_fixedID))
            {
                continue;
            }
            
            //Hits with a hot pixel are ignored
            if(hitContainsHotPixels(hit))
            {
                continue;
            }
           
            // And its position
            const double * pos = hit->getPosition();
            
            // Get the detector ID of this hit
            int iHitID = auxfunc::hitDecoder(hit)["sensorID"]; 
            // Obtain the PreAligner (i.e. the histogram) related
            // with the detector of the current hit
            auto paIt = std::find_if(_preAligners.begin(),_preAligners.end(), 
                    [&iHitID](const PreAligner & _pa) { return _pa.getIden() == iHitID; });
            if(paIt == _preAligners.end())
            {
                streamlog_out ( ERROR5 ) << "Mismatched hit at " << pos[2] << endl;
            }

            PreAligner & pa = *paIt;
            // Calculate the reference hit with respect this current one
            double correlationX =  refPos[0] - pos[0];
            double correlationY =  refPos[1] - pos[1];
            
            const int idZ = _sensorIDtoZOrderMap[ iHitID ];
            // Just store it only if is inside the allowed correlation
            if( (_residualsXMin[idZ] < correlationX ) 
                        && ( correlationX < _residualsXMax[idZ]) 
                        && (_residualsYMin[idZ] < correlationY ) 
                        && ( correlationY < _residualsYMax[idZ]) )
            {
                residX.push_back( correlationX );
                residY.push_back( correlationY );
                currentX.push_back( pos[0] );
                currentY.push_back( pos[1] );
                prealign.push_back(&pa);
            }
        }
        
        // 2. Accepted hits must be present in at least a minimum number
        //    of detector planes, in that case is considered a serie of 
        //    points correlated 
        if(prealign.size() > static_cast< unsigned int >(_minNumberOfCorrelatedHits) 
                && residX.size() == residY.size() ) 
        {
            for(unsigned int ii = 0; ii < prealign.size(); ++ii ) 
            {
                prealign[ii]->addPoint( residX[ii], residY[ii] );
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
                if( _fillHistos ) 
                {
                    ( dynamic_cast<AIDA::IHistogram1D*> (_hitXCorr[ prealign[ii]->getIden() ] ) )->fill( residX[ii] );
                    ( dynamic_cast<AIDA::IHistogram1D*> (_hitYCorr[ prealign[ii]->getIden() ] ) )->fill( residY[ii] );
                    ( dynamic_cast<AIDA::IHistogram2D*>(_2dHitHistos[prealign[ii]->getIden()]))->fill(currentX[ii],currentY[ii] );
                }
#endif
            }
        }
    }
}

bool EUTelPreAlign::hitContainsHotPixels( TrackerHitImpl   * hit) 
{
    // if no hot pixel map was loaded, just return here
    if( _hotPixelMap.size() == 0) 
    {
        return 0;
    }
    
    try
    {
        LCObjectVec clusterVector = hit->getRawHits();
        if ( hit->getType() == kEUTelSparseClusterImpl ) 
        {
            TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );
            eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel >* cluster =
                new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel >(clusterFrame);
            int sensorID = cluster->getDetectorID();
            
            for(unsigned int iPixel = 0; iPixel < cluster->size(); iPixel++ )
            {
                EUTelGenericSparsePixel m26Pixel;
                cluster->getSparsePixelAt( iPixel, &m26Pixel);
                {
                    try
                    {
                        if( std::find(_hotPixelMap.at(sensorID).begin(), 
                                    _hotPixelMap.at(sensorID).end(),
                                    std::make_pair(m26Pixel.getXCoord(),m26Pixel.getYCoord()))
                                != _hotPixelMap.at(sensorID).end())
                        { 
                            delete cluster;                        			  
                            return true; // if TRUE  this hit will be skipped
                        }
                    }
                    catch(const std::out_of_range& oor)
                    {
                        streamlog_out(DEBUG0) << " Could not find hot pixel map for sensor ID " 
                            << sensorID << ": " << oor.what() << endl;
                    }
                }
            }
            
            delete cluster;
            return false;
        } 
        else if( hit->getType() == kEUTelBrickedClusterImpl )
        {
            // fixed cluster implementation. Remember it
            //  can come from
            //  both RAW and ZS data
            streamlog_out ( WARNING5 ) << " Hit type kEUTelBrickedClusterImpl is not " 
                << "implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;
        } 
        else if( hit->getType() == kEUTelDFFClusterImpl ) 
        {
            // fixed cluster implementation. Remember it can come from
            // both RAW and ZS data
            streamlog_out ( WARNING5 ) << " Hit type kEUTelDFFClusterImpl is not "
                << "implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;
        } 
        else if ( hit->getType() == kEUTelFFClusterImpl ) 
	{
            // fixed cluster implementation. Remember it can come from
	    // both RAW and ZS data
	    streamlog_out ( WARNING5 ) << " Hit type kEUTelFFClusterImpl is not "
                << "implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;
        }
        else
	{
            streamlog_out ( WARNING5 ) << " Hit type is not known and is not "
                << "implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;
        }
    }
    catch (exception& e)
    {
        streamlog_out(ERROR4) << "something went wrong in EUTelPreAlign::hitContainsHotPixels: " << e.what() << endl;
    }
    catch(...)
    { 
        // if anything went wrong in the above return FALSE, meaning do not skip this hit
        streamlog_out(ERROR4) << "something went wrong in EUTelPreAlign::hitContainsHotPixels " << endl;
        return 0;
    }
    
    // if none of the above worked return FALSE, meaning do not skip this hit
    return 0;
}
      
void EUTelPreAlign::end()
{
    LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );
    
    // First check if we did not found any histograms for this sensor,
    // i.e. the reference sensor, then just fill the alignment offset to 0
    for(const auto & sensorID: _sensorIDVec)
    {
        auto pa_it = std::find_if(_preAligners.begin(),_preAligners.end(),
                [&sensorID](const PreAligner & _pa) { return _pa.getIden() == sensorID; });
        if(pa_it == _preAligners.end())
        {
            EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
            constant->setXOffset( 0.0 );
            constant->setYOffset( 0.0 );
            constant->setSensorID( sensorID );
            constantsCollection->push_back( constant );
            streamlog_out ( MESSAGE5 ) << (*constant) << endl;
            continue; 
        }
    }
    
    for(auto & pa : _preAligners)
    {
        const int sensorID = pa.getIden();
        std::vector<int>::iterator it = std::find(_ExcludedPlanes.begin(),_ExcludedPlanes.end(),sensorID);
        std::vector<int>::iterator itXCoord = std::find(_ExcludedPlanesXCoord.begin(),_ExcludedPlanesXCoord.end(),sensorID);
        std::vector<int>::iterator itYCoord = std::find(_ExcludedPlanesYCoord.begin(),_ExcludedPlanesYCoord.end(),sensorID);
        
        // Set the offset of each detector (but those explicitely excluded) by 
        // using the peak of the histogram of distances (refHit-current_det_hit)
        EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();        
        if(it == _ExcludedPlanes.end())
        {
            if( itXCoord == _ExcludedPlanesXCoord.end() && abs( pa.getPeakX() ) < 1000 )
                constant->setXOffset( -1.0*pa.getPeakX() );
            else
                constant->setXOffset( 0.0 );
            
            if(  itYCoord == _ExcludedPlanesYCoord.end() && std::fabs(pa.getPeakY() ) < 1000. )
                constant->setYOffset( -1.0 *pa.getPeakY() );
            else
                constant->setYOffset( 0.0 );
        }
        else
        {
            constant->setXOffset(0.0);
            constant->setYOffset(0.0);
        }
        constant->setSensorID( sensorID );
        constantsCollection->push_back( constant );
        
        //Also update the EUTelGeometry descr.
        double updatedXOff = geo::gGeometry().siPlaneXPosition(sensorID) + pa.getPeakX();
        double updatedYOff = geo::gGeometry().siPlaneYPosition(sensorID) + pa.getPeakY();
        
        geo::gGeometry().setPlaneXPosition(sensorID, updatedXOff);
        geo::gGeometry().setPlaneYPosition(sensorID, updatedYOff);
        
        streamlog_out ( MESSAGE5 ) << (*constant) << endl;
    }
    
    //if we dont dump into the new gear file, we write out the old database
    if(!_dumpGEAR)
    {
        LCWriter* lcWriter = LCFactory::getInstance()->createLCWriter();
        try 
        {
            lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW    );
        } 
        catch ( IOException& e ) 
        {
            streamlog_out ( ERROR4 ) << e.what() << endl;
            exit(-1);
        }
        
        streamlog_out ( MESSAGE5 ) << "Writing to " << _alignmentConstantLCIOFile << endl;
        
        LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
        lcHeader->setRunNumber( 0 );
        lcWriter->writeRunHeader(lcHeader);
        delete lcHeader;
        LCEventImpl * event = new LCEventImpl;
        event->setRunNumber( 0 );
        event->setEventNumber( 0 );
        LCTime * now = new LCTime;
        event->setTimeStamp( now->timeStamp() );
        delete now;
        
        streamlog_out( DEBUG5 ) << " adding Collection " << "alignment " << endl;
        event->addCollection( constantsCollection, "alignment" );
        lcWriter->writeEvent( event );
        delete event;
        lcWriter->close();
    }
    else
    {
        //Write updated GEAR file
        marlin::StringParameters* MarlinStringParams = marlin::Global::parameters;
        std::string outputFilename = (MarlinStringParams->getStringVal("GearXMLFile")).substr(0, (MarlinStringParams->getStringVal("GearXMLFile")).size()-4);
        streamlog_out(MESSAGE5) << "Writing updated GEAR file with filename: " << outputFilename+"_pre.xml" << std::endl;
        geo::gGeometry().writeGEARFile(outputFilename+_GEARFileSuffix+".xml");
        
        //in case we don't write out the collection, we need to delete ourself as we don't pass the collection to LCIO
        delete constantsCollection;
    }
}

