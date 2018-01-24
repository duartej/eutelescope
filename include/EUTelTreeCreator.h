/*
 * Created by Jordi Duarte-Campderros
 *  (2017 CERN/IFCA)
 *
 * Creates a simple tree from a merged collection
 * of hits (where they are included the DUT hits)
 * and tracks
 *
 *  email:jorge.duarte.campderros@cern.ch
 */

#ifndef EUTELTREECREATOR_H
#define EUTELTREECREATOR_H

// eutelescope includes ".h"
#include "EUTelUtility.h"


// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
//#include <EVENT/LCCollection.h>
//#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <vector>
#include <map>

// forward declaration
class TTree;
class TFile;

namespace eutelescope 
{
    class EUTelTreeCreator : public marlin::Processor 
    {
        private:
            DISALLOW_COPY_AND_ASSIGN(EUTelTreeCreator)
                
        public:
            virtual Processor * newProcessor() 
            {
                return new EUTelTreeCreator;
            }
            
            //! Default constructor
            EUTelTreeCreator ();
            
            virtual void init ();
            virtual void processRunHeader (LCRunHeader * run);
            virtual void processEvent (LCEvent * evt);
            virtual void end();
        
        private:
            // Auxiliary function
            void clear_branches();

            //! The sensitive axis of the sensors
            int _sensitive_axis;

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

            //! The telescope planes
            std::vector<int> _telescopePlanes;

            //! The DUT planes 
            std::vector<int> _dutPlanes;

            // The File to store
            TFile * _file;
            std::string _filename;

            // And the tree
            TTree * _tree;

            // The auxiliary map for the tree filling
            std::map<std::string,std::vector<int>*> _branches_I;
            std::map<std::string,std::vector<float>*> _branches;
            std::map<std::string,int> _branches_int;
            std::map<std::string,float> _branches_float;

            int _iRun;
            int _iEvt;
    };
    
    //! A global instance of the processor
    EUTelTreeCreator gEUTelTreeCreator;
}
#endif
