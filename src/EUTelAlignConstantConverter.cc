/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelAlignConstantConverter.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Global.h"
#include "marlin/StringParameters.h"

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <Exceptions.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace marlin;
using namespace eutelescope;

EUTelAlignConstantConverter::EUTelAlignConstantConverter():Processor("EUTelAlignConstantConverter") 
{
    std::vector<std::string> alignmentNames = {"prealign"};
    // modify processor description
    _description = "EUTelAlignConstantConverter creates an updated GEAR file with the loaded alignment constants.";
    registerOptionalParameter("AlignmentCollectionNames", "The name of the alignment constants in the same order as they were created", 
            _collectionNames, alignmentNames );
    registerOptionalParameter("NewGEARSuffix", "Suffix for the new GEAR file, set to empty string "\
            "(this is not default!) to overwrite old GEAR file", _GEARFileSuffix, std::string("_aligned") );
}

void EUTelAlignConstantConverter::init()
{
    //this method is called only once even when the rewind is active usually a good idea to
    printParameters ();
    
    // set to zero the run and event counters
    _iRun = 0;
    _iEvt = 0;
    
    //Getting access to geometry description
    geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);

    //check if the GEAR manager pointer is not null!
    if( Global::GEAR == 0x0 ) 
    {
        streamlog_out(ERROR2) << "The GearMgr is not available, for an unknown reason." << std::endl;
        throw InvalidGeometryException("GEAR manager is not initialised");
    }
 
    // Build the sensorID map to relate the alignement constants. Each
    // element contain the vector of alignment constants to be (ordered) applied 
    // to the given sensor
    const std::vector<int> sensorIDVec = geo::gGeometry().sensorIDsVec();
    for(const auto & id: sensorIDVec)
    {
        _alignmentMap[id] = {};
    }

    // FIXME: if  Global::GEAR->getDetectorName() != "" -> ERROR message saying that, otherwise
    // peta because there is no detector name in the LCIO db files
  
    streamlog_out( MESSAGE4 ) << "End of initialisation" << std::endl;
}

void EUTelAlignConstantConverter::processRunHeader(LCRunHeader* rdr) 
{
    auto header = std::make_unique<EUTelRunHeaderImpl>(rdr);
    header->addProcessor( type() ) ;
    
    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, instead of barely
    // quitting ask the user what to do.
    if( header->getGeoID() != (int)geo::gGeometry().getSiPlanesLayoutID() ) 
    {
        streamlog_out( ERROR2 ) << "Error during the geometry consistency check: " << std::endl;
        streamlog_out( ERROR2 ) << "The run header says the GeoID is " << header->getGeoID() << std::endl;
        streamlog_out( ERROR2 ) << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << std::endl;
    }

    // increment the run counter
    ++_iRun;
}

void EUTelAlignConstantConverter::processEvent(LCEvent* event) 
{
    ++_iEvt;

    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
    
    // Get the alignment collection from the event
    for(const auto & name: _collectionNames)
    {
        LCCollectionVec * alignColVec = nullptr;
        try 
        {
            alignColVec = dynamic_cast <LCCollectionVec*>(evt->getCollection(name));
        } 
        catch(DataNotAvailableException& e) 
        {
            streamlog_out(ERROR) << "No Alignment ["<< name <<"] collection found " << std::endl;
            throw DataNotAvailableException("No alignment collection found. Check input files!");
        }
        // Fill the alignment map
        for(unsigned int i=0; i < alignColVec->size(); ++i)
        {
            EUTelAlignmentConstant * aligncte = static_cast<EUTelAlignmentConstant*>(alignColVec->getElementAt(i));
            _alignmentMap[aligncte->getSensorID()].push_back(aligncte);
        }
    }
}

void EUTelAlignConstantConverter::end() 
{
    // Prepare the elements to perform the gear update
    // -----------------------------------------------

    for(const auto & themap: _alignmentMap)
    {
        const int sensorID = themap.first;
        streamlog_out(DEBUG) << "Alignment on sensor " << sensorID << " calculated in " 
            << themap.second.size() << " steps: " << std::endl ;
        int steps = 0;
        for(const auto * acte: themap.second)
        {
            const double alpha = acte->getAlpha();
            const double beta  = acte->getBeta();
            const double gamma = acte->getGamma();
            const double xOffset = acte->getXOffset();
            const double yOffset = acte->getYOffset();
            const double zOffset = acte->getZOffset();
            
            streamlog_out(DEBUG) << " [" << steps << "]: " 
                << std::setw(8) << "xOffset= " << xOffset << " [mm], " 
                << std::setw(8) << "yOffset= " << yOffset << " [mm], " 
                << std::setw(8) << "zOffset= " << zOffset << " [mm], " 
                << std::setw(6) << "alpha: " << alpha*1e3 << " [mrad], "
                << std::setw(6) << "beta: " << beta*1e3 << " [mrad], "
                << std::setw(6) << "gamma: " << gamma*1e3 << " [mrad], "
                << std::endl; 
            // The old rotation matrix is well defined by GEAR file (and in the previous alignment steps,
            // as it is updated at the end of each step)
            Eigen::Matrix3d rotOld = geo::gGeometry().rotationMatrixFromAngles(sensorID);
            // The position vector (idem)
            Eigen::Vector3d oldOffset;
            oldOffset << geo::gGeometry().siPlaneXPosition(sensorID), geo::gGeometry().siPlaneYPosition(sensorID), geo::gGeometry().siPlaneZPosition(sensorID);
        

            //The new rotation matrix is obtained via the alpha, beta, gamma from MillepedeII
            Eigen::Matrix3d rotAlign = geo::gGeometry().rotationMatrixFromAngles( -alpha, -beta, -gamma);
            //The corrected rotation is given by: rotAlign*rotOld, from this rotation we can extract the
            //updated alpha', beta' and gamma'
            Eigen::Vector3d newCoeff = geo::gGeometry().getRotationAnglesFromMatrix(rotAlign*rotOld);

            //std::cout << "Old rotation matrix:\n " << rotOld << std::endl; 
            //std::cout << "Align rotation matrix:\n " << rotAlign << std::endl; 
            //std::cout << "Updated coefficients:\n " << newCoeff << std::endl; 
            //std::cout << "This results in the updated rotations (alpha', beta', gamma'): " << newCoeff[0] << ", " << newCoeff[1] << ", " << newCoeff[2] << std::endl;
            
            geo::gGeometry().setPlaneXPosition(sensorID, oldOffset[0]-xOffset);
            geo::gGeometry().setPlaneYPosition(sensorID, oldOffset[1]-yOffset);
            geo::gGeometry().setPlaneZPosition(sensorID, oldOffset[2]-zOffset);

	    geo::gGeometry().setPlaneXRotationRadians(sensorID,  newCoeff[0]);
            geo::gGeometry().setPlaneYRotationRadians(sensorID,  newCoeff[1]);
            geo::gGeometry().setPlaneZRotationRadians(sensorID,  newCoeff[2]);

            ++steps;
        }
    }
    marlin::StringParameters* MarlinStringParams = marlin::Global::parameters;
    std::string outputFilename = (MarlinStringParams->getStringVal("GearXMLFile")).substr(0, (MarlinStringParams->getStringVal("GearXMLFile")).size()-4);
    std::cout << "GEAR Filename: " << outputFilename+_GEARFileSuffix+".xml" << std::endl;
    geo::gGeometry().writeGEARFile(outputFilename+_GEARFileSuffix+".xml");
    streamlog_out( MESSAGE2 ) << std::endl << "Successfully finished" << std::endl;
}
