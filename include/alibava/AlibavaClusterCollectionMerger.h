/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros (IFCA/CERN)
 *  * 2017-08-18 Included the ability to deal with a reference
 *    alibava sensor
 */


#ifndef ALIBAVACLUSTERCOLLECTIONMERGER_H
#define ALIBAVACLUSTERCOLLECTIONMERGER_H 1

// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaRunHeaderImpl.h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include<string>
#include<vector>
#include<set>
#include<functional>

namespace alibava 
{
    // forward declaration
    class AlibavaEventImpl;

    class AlibavaClusterCollectionMerger : public marlin::DataSourceProcessor    
    {
        public:
            //! Default constructor
	    AlibavaClusterCollectionMerger ();
	    
	    //! New processor
	    /*! Return a new instance of a AlibavaClusterCollectionMerger. It is
	     *  called by the Marlin execution framework and shouldn't be used
	     *  by the final user.
	     */
	    virtual AlibavaClusterCollectionMerger * newProcessor ();
	    
	    //! Creates events from the eudaq software
	    virtual void readDataSource (int numEvents);
	    
	    //! Init method
	    /*! It is called at the beginning of the cycle and it prints out
	     *  the parameters.
	     */
	    virtual void init ();
	    
	    //! End method
	    /*! It prints out a good bye message
	     */
	    virtual void end ();
	    
	protected:
            ///////////////
            // Telescope //
            ///////////////
            //! This is the input file name that telescope cluster collections stored
            std::string _telescopeFileName;
            
            //! Name of the cluster pulse collection of telescope data
            std::string _telescopePulseCollectionName;
            
            //! Name of the sparse cluster collection of telescope data
            std::string _telescopeSparseCollectionName;
            
            
            /////////////
            // Alibava //
            /////////////
            //! This is the input file name that alibava cluster collections stored
            std::string _alibavaFileName;
            //! This is the input file name that alibava (used as the reference sensor)
            // cluster collections stored
            std::string _alibavaReferenceFileName;
            
            //! Name of the cluster pulse collection of alibava data
            std::string _alibavaPulseCollectionName;
            
            //! Name of the sparse cluster collection of alibava data
            std::string _alibavaSparseCollectionName;
            
            
	    ////////////
	    // Output //
	    ////////////
	    
	    //! Name of the merged/output cluster pulse collection
	    std::string _outputPulseCollectionName;
	    
	    //! Name of the merged/output sparse cluster collection.
	    /*! This parameter is hard coded in other EUTelescope
	     *  processors. It has to be "original_zsdata"
	     */
	    std::string _outputSparseCollectionName;
	    
	    //! Event ID difference
	    int _eventIDDiff;

            //! used sensor ids
            std::set<int> _usedSensorIDs;
	    
	private:
            //! Whether or not is present the reference alibava sensor
            bool _refPresent;
            //! The reference sensor ID 
            int _refSensorID;
            //! Functional to evaluate the reference sensor ID, which is defined
            //! as the lowest number not used by any of the other planes
            // XXX TO BE DEPRECATED
            std::function<int()> getReferenceSensorID;
            //! Functional to update the sensorIDs, but once they are fill, just
            //! a dummy function returning the input
            // XXX TO BE DEPRECATED
            std::function<int(const int &)>getSensorID;
            // ID for the dut and reference introduced by the user
            int _dutID;
            int _referenceID;
            //! A function to be define its behaviour depending if the
            //! a reference alibava sensor must be processed or not
            std::function<bool(LCReader *,AlibavaEventImpl*& )> _checkAlibavaReferenceEvent;
            void copyClustersInCollection(LCCollectionVec * outputPulseColVec, LCCollectionVec * outputSparseColVec, 
                    LCCollectionVec * inputPulseColVec, LCCollectionVec * inputSparseColVec,
                    bool isDUTSensor=false,bool isReferenceSensor=false);
	};	
	//! A global instance of the processor
	AlibavaClusterCollectionMerger gAlibavaClusterCollectionMerger;	
} // end of alibava namespace
#endif

