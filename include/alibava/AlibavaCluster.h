/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  - add extra noise data-member 
 *  - some code style
 *
 */

#ifndef ALIBAVACLUSTER_H_
#define ALIBAVACLUSTER_H_ 1

// alibava includes ".h"
#include "ALIBAVA.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <vector>

namespace alibava
{
    class AlibavaCluster 
    {
        private:
            // The channel number of the strip composing the cluster
            std::vector<int> _channums;
            // The evaluated signal for each channel
            std::vector<float> _signals;
            // The noise corresponding to each channel
            std::vector<float> _noise;
            float _eta;
            int _chipNum;
            int _seedChanNum;
            int _clusterID;
            bool _isSensitiveAxisX;
            int _signalPolarity;
            
        public:
            AlibavaCluster();
            AlibavaCluster(lcio::TrackerDataImpl* trkdata);
            ~AlibavaCluster();

            // Get the channel number correponding to the 
            // imember-strip of the cluster
	    int getChanNum(int imember);
            // Get the signal corresponding to the i-strip of the cluster
	    float getSignal(int imember);
            // Overload function to obtain the signal in electrons for a given strip
            float getSignal(const int & imember,const std::vector<float> calibration);
            // Get the total signal of the cluster
	    float getTotalSignal();
            // Overload function to obtain the signal of the cluster in electrons
	    float getTotalSignal(const std::vector<float> calibration);
            // Get the noise of a given cluster
            float getNoise(int imember);
            // Get the total noise of the cluster
            float getNoise();
            // Get the signal to noise ratio 
	    float getTotalSNR(EVENT::FloatVec noiseVec);
            // Get the signal to noise ratio 
            float getTotalSNR();

	    // Add a strip to be part of the cluster
	    void add(int achannum, float asignal, float noise);
	    int getClusterSize();
	    bool has_seed();
	    void print();
	    
	    float getCenterOfGravity();
             
            // Calculate the center position of the cluster
            // using an unbiased center position finder algorithm
            // Head-tail algorithm 
            // (R. Turchetta, “Spatial resolution of silicon microstrip detectors”, 1993)
            float getUnbiasedCenterPosition();
	    
	    void createTrackerData(lcio::TrackerDataImpl * alibavaCluster);
	    
            // The eta is calculated using only the strips
            // belonging to the cluster
            float getEtaFromCluster();

	    // setter / getter for _eta
	    float getEta();
	    void setEta(float eta);
	    
	    // setter / getter for _chipNum
	    int getChipNum();
	    void setChipNum(int chipnum);
	    
	    // setter / getter for _seedChanNum
	    int getSeedChanNum();
	    void setSeedChanNum(int seedChanNum);
	    
	    // setter / getter for _clusterID
	    int getClusterID();
	    void setClusterID(int clusterID);
	    
	    // setter / getter for _isSensitiveAxisX
	    void setIsSensitiveAxisX(bool isSensitiveAxisX);
	    bool getIsSensitiveAxisX();
	    
	    // setter / getter for _signalPolarity
	    int getSignalPolarity();
	    void setSignalPolarity(int signalPolarity);
	};
	
} // end of alibava namespace

#endif /* ALIBAVACLUSTER_H_ */
