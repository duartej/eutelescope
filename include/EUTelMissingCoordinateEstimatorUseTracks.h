/*
 * Created by Jordi Duarte-Campderros
 *  (2017 CERN/IFCA)
 *
 * Missing Coordinate Estimator using a input track collection
 * 
 * Given a DUT with only with measured coordinate (a micro-strip sensor),
 * and an input reconstructed tracks (coming from a telescope),
 * the mising coordinate will be calculated by using the tracks.
 * If a track matchs the detector plane within a iluminated strip
 * (plus uncertainties: track and strip resolutions), the missing
 * coordinate is assigned by using the track.
 * One needs to used this after pre-alignment in order to align
 * the DUT with the telescope.
 *
 *  email:jorge.duarte.campderros@cern.ch
 */

#ifndef EUTELMISSINGCOORDINATEESTIMATOUSETRACKS_H
#define EUTELMISSINGCOORDINATEESTIMATORSETRACKS_H

// eutelescope includes ".h"
#include "EUTelUtility.h"


// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include "TObject.h"

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <set>

namespace eutelescope 
{
    class EUTelMissingCoordinateEstimatorUseTracks : public marlin::Processor 
    {
        private:
            DISALLOW_COPY_AND_ASSIGN(EUTelMissingCoordinateEstimatorUseTracks)
                
        public:
            virtual Processor * newProcessor() 
            {
                return new EUTelMissingCoordinateEstimatorUseTracks;
            }
            
            //! Default constructor
            EUTelMissingCoordinateEstimatorUseTracks ();
            
            virtual void init ();
            virtual void processRunHeader (LCRunHeader * run);
            virtual void processEvent (LCEvent * evt);
            virtual void end();
            //! Histogram booking
            /*! Does nothing
             */
            void bookHistos();
        
        protected:
            //! Input TrackerHit collection name
            /*! This is the name the user wants to give to the input hit
             *  collection.
             */
            std::string _inputHitCollectionName;
            
            //! Input Tracks collection name
            /*! This is the name the user wants to give to the input track
             *  collection from the telescope
             */
            std::string _inputTrackCollectionName;
        
            //! Output TrackerHit collection name
            /*! This is the name the user wants to give to the output hit
             *  collection.
             */
            std::string _outputHitCollectionName;
            
            //! DUT Planes
            /*! This is the list of sensorIDs that missing coordinate of their 
             *  hits needs to be found. Notice that if the specified coordinate 
             *  already exists it will be overwritten
             */
            EVENT::IntVec _dutPlanes;
        
            //! Missing Coordinate
           /*! The coordinate axis that needs to be estimated. 
             *  You have to set this to either X or Y.
             */
            std::string _missingCoordinate;
        
            //! Max Residual
            /*! This processor will look for a closest hits (in known coordinate)
             *  to determine if the hits are correlated to the track
             *  The hits will be considered as correlated if the residual is smaller
             *  than MaxResidual
             */
            EVENT::FloatVec _maxResidual_prov;
            std::map<int,float> _maxResidual;
        
            //! Telescope resolution
            /*! Use the telescope resolution to obtain the binning size for the unknown
             *  coordinate binning (defining also its error)
             */
            float _telescopeResolution;

            //! Clone Hit
            /*! This method is used to clone TrackerHitImpl object
             */
            TrackerHitImpl* cloneHit(TrackerHitImpl *inputHit);

        private:
            //! Run number
            int _iRun;
        
            //! Event number
            int _iEvt;

            //! Histogram maps: retrieve them back using dynamic_cast
            std::map<std::string,TObject*> _histoMap;

            static std::string _hitplotname;
            static std::string _missinghits;
            static std::string _duthitNtracks;
            
            //! Missing hit position
            unsigned int _missingHitPos;

            //! Known hit position
            unsigned int _knownHitPos;
        
            //! Number of DUT hits
            unsigned int _nDutHits;
        
            //! Number of Dut hits created
            unsigned int _nDutHitsCreated;
	       
 	    //! Number of Expected created hit per DUT hit
 	    unsigned int _maxExpectedCreatedHitPerDUTHit;
	
	    //! Count number of created hit per DUT Hit
	    std::vector<unsigned int> _numberOfCreatedHitPerDUTHit; 
    };
    
    //! A global instance of the processor
    EUTelMissingCoordinateEstimatorUseTracks gEUTelMissingCoordinateEstimatorUseTracks;
    
}
#endif
