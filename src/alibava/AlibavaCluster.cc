/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 */

// alibava includes ".h"
#include "AlibavaCluster.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <iomanip>
#include <map>
#include <iterator>


using namespace std;
using namespace lcio;
using namespace alibava;

AlibavaCluster::AlibavaCluster() :
    _channums(),
    _signals(),
    _eta(-1),
    _chipNum(-1),
    _seedChanNum(-1),
    _clusterID(-1),
    _isSensitiveAxisX(true),
    _signalPolarity(-1)
{
    // does nothing
}

AlibavaCluster::AlibavaCluster(TrackerDataImpl* trkdata) :
    _channums(),
    _signals(),
    _eta(-1),
    _chipNum(-1),
    _seedChanNum(-1),
    _clusterID(-1),
    _isSensitiveAxisX(true),
    _signalPolarity(-1)
{
    CellIDDecoder<TrackerDataImpl> clusterIDDecoder(ALIBAVA::ALIBAVACLUSTER_ENCODE);
    _clusterID = static_cast<int>(clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERID]);
    _seedChanNum = static_cast<int>(clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_SEED]);
    _chipNum = static_cast<int>(clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_CHIPNUM]);
    
    int sensitiveAxisX= static_cast<int>(clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX]);
    if(sensitiveAxisX == 0)
    {
        _isSensitiveAxisX = false;
    }
    else
    {
        _isSensitiveAxisX = true;
    }
		
    int negSignal = static_cast<int>(clusterIDDecoder( trkdata )[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE]);
    if(negSignal == 0)
    {
        _signalPolarity = 1;
    }
    else
    {
        _signalPolarity = -1;
    }
	
    EVENT::FloatVec data = trkdata->getChargeValues();
    // first number is eta,
    // then it goes like channel number, signal, noise // see createTrackerData
    // So data size has to be multiple of 3 (minus the eta)!
    if((data.size()-1) % 3 != 0)
    {
        streamlog_out (ERROR5) << "Size in TrackerData that stores " 
            << "cluster information is 3-multiple! Either channel number, " 
            << " signal or noise information is missing!"<< std::endl;
    }
    // first number is eta, others are channel number and corresponding signal
    this->setEta(data.at(0));
    for(unsigned int imember=1; imember<data.size(); ++imember) 
    {
        _channums.push_back( int(data[imember]) );
        ++imember;
        _signals.push_back( data[imember] );
        ++imember;
        _noise.push_back( data[imember] );
    }
}

AlibavaCluster::~AlibavaCluster(){
	
}

void AlibavaCluster::createTrackerData(lcio::TrackerDataImpl * alibavaCluster)
{
    // put channel information in TrackerData
    // first number we put is channel number then signal
    EVENT::FloatVec dataToStore;
	
    // first store eta
    dataToStore.push_back(getEta());
    // then channel number, signal and noise for each member
    for(int imember=0; imember<getClusterSize(); imember++) 
    {
        dataToStore.push_back( float(getChanNum(imember)) );
        dataToStore.push_back( getSignal(imember) );
        dataToStore.push_back( getNoise(imember) );
    }
    alibavaCluster->setChargeValues(dataToStore);
    // Cell ID encoding will done in the main processor	
}

/*
float AlibavaCluster::getEta(){
	// if clustersize is 1, nothing to calculate, returns -1 
	if (getClusterSize()==1)
		return -1;
	
	float leftNeighs = 0;
	float rightNeighs = 0;
	float centerOfGravity = getCenterOfGravity();
	
	for (int imember=0; imember<getClusterSize(); imember++) {
		// since we will compare this number with float
		float memberChanNum = float ( getChanNum(imember) );
		if ( memberChanNum < centerOfGravity )
			leftNeighs += getSignal(imember);
		else if ( memberChanNum > centerOfGravity )
			rightNeighs += getSignal(imember);
	}
	
	float eta = leftNeighs/(leftNeighs+rightNeighs);
	return eta;
}
*/

float AlibavaCluster::getEta()
{
    return _eta;
}

// [JDC]
float AlibavaCluster::getEtaFromCluster()
{
    // [1] It would define the postion by subtracting the
    // right strip  x = X_R - p f(eta)
    // In case of clusters > 2, it will use the two higher
    // signal channels
    
    // Calculate the charge sharing for clusters > 2
    if(getClusterSize()==1)
    {
        return -1;
    }
    // Find the left cluster and its signal 
    //const int seedStripNum = getSeedChanNum();
    // Note that the seed is always the first element of 
    // the cluster
    //const float seedSignal = getSignalPolarity()*getSignal(0);
    //float leftSignal = 0;
    //float rightSignal = 0;
    std::multimap<float,int> channelsMemberOrderMap;
    // Store the multimap to use the highest signal channels
    for(int imember=0; imember<getClusterSize(); ++imember)
    {
        channelsMemberOrderMap.insert(std::pair<float,int>(static_cast<float>(getSignalPolarity()*getSignal(imember)),getChanNum(imember)));
    }
    // Use only the two last elements (the highest signals) of the map
    auto highestSignal = channelsMemberOrderMap.rbegin();
    auto prevToHighestSignal = std::next(channelsMemberOrderMap.rbegin());
    
    // A guess
    float leftSignal = highestSignal->first;
    float rightSignal= prevToHighestSignal->first;
    // Check if the guess was right, i.e. the highest signal channel
    // is the leftest channel, otherwise, change the order
    if( highestSignal->second > prevToHighestSignal->second)
    {
        leftSignal = prevToHighestSignal->first;
        rightSignal= highestSignal->first;
//std::cout << getClusterSize() << "-- LEFT Signal: " << leftSignal << " (Ch:" << prevToHighestSignal->second << "). RIGHT signal: " << rightSignal << " (Ch:" << highestSignal->second << ")" << std::endl;
//    }
//else{
//std::cout << getClusterSize() << "-- LEFT Signal: " << leftSignal << " (Ch:" << highestSignal->second << "). RIGHT signal: " << rightSignal << " (Ch:" << prevToHighestSignal->second << ")" << std::endl;
}

    /*{
        const int ichan = getChanNum(imember);
        // Seed channel to the left description (see [1])
        if(ichan <= seedStripNum)
        {
            leftSignal += getSignalPolarity()*getSignal(imember);
        }
        else
        {
            rightSignal += getSignalPolarity()*getSignal(imember);
        }
    }*/
    //Get the left (the lowest channel number) and right 
    // (the highest channel number) signals
    //const float leftSignal  = static_cast<float>(getSignalPolarity()*getSignal(channelsMemberOrderMap.begin()->second));
    //const float rightSignal = static_cast<float>(getSignalPolarity()*getSignal(channelsMemberOrderMap.end()->second));
    
    /*float eta = -1;
    if(leftSignal > rightSignal) 
    {
        // then seed channel is on the right
        eta = leftSignal/(leftSignal + seedSignal);
    }
    else 
    {
        // seed channel is on the left
        eta = seedSignal/(seedSignal + rightSignal);
    }
    return eta;*/
    
    return leftSignal/(leftSignal+rightSignal);
}

void AlibavaCluster::setEta(float eta)
{
    _eta = eta;
}


float AlibavaCluster::getCenterOfGravity()
{
    if(getClusterSize()==1)
    {
        return getChanNum(0); // channel number of the only member
    }

    float totalsignal = 0;
    float weightedSum = 0;
    for(int imember=0; imember<getClusterSize(); imember++) 
    {
        totalsignal += getSignal(imember);
        weightedSum += getSignal(imember)*getChanNum(imember);
    }
	
    return weightedSum/totalsignal;
}

// 
float AlibavaCluster::getUnbiasedCenterPosition()
{
    // Note as it is considered number of channels, the
    // pitch is set to 1 unity (the separation between 
    // two channels)
    
    if(getClusterSize()==1)
    {
        return float(getSeedChanNum());
    }

    // The center is the seed
    float xrec = float(getSeedChanNum());
    // find the edge channel using a map order (low, i.e left 
    // to high, i.e right)
    std::map<int,int> channelMember;
    for(int imember=0; imember < getClusterSize(); ++imember)
    {
        channelMember.emplace(getChanNum(imember),imember);
    }
    // The edges 
    const int leftChannel  = channelMember.begin()->second;
    const float signalLeft = getSignalPolarity()*getSignal(leftChannel);

    const int rightChannel  = channelMember.end()->second;
    const float signalRight = getSignalPolarity()*getSignal(rightChannel);
    
    // special 2-clusters
    if(getClusterSize()==2)
    {
        return xrec+(signalRight-signalLeft)/(3.0*std::max(signalLeft,signalRight));
    }

    float partialSignal = 0.0;
    // adding up all the charges but the edges strips
    for(auto it = std::next(channelMember.begin()); it != std::prev(channelMember.end()); ++it)
    {
        partialSignal += getSignalPolarity()*getSignal(it->second);
    }
    partialSignal /= float(getClusterSize()-2);

    return xrec+(signalRight-signalLeft)/(2.0*partialSignal);  
}

int AlibavaCluster::getChanNum(int imember) 
{
    if(imember >= getClusterSize()) 
    {
        streamlog_out(ERROR5)<< "Not enough member in this AlibavaCluster"<<std::endl;
        return -1;
    }
    return _channums[imember];
}

float AlibavaCluster::getSignal(int imember)
{
    if(imember >= getClusterSize()) 
    {
        streamlog_out(ERROR5)<< "Not enough members in this AlibavaCluster"<<std::endl;
        return 0;
    }
    return _signals[imember];
}

// Overload function getting electrons (instead of ADC)
float AlibavaCluster::getSignal(const int & imember, const std::vector<float> calv)
{
    if(imember >= getClusterSize()) 
    {
        streamlog_out(ERROR5)<< "Not enough members in this AlibavaCluster"<<std::endl;
        return 0;
    }
    // Check the validity of calv
    if(calv.size() != ALIBAVA::NOOFCHANNELS)
    {
        streamlog_out(ERROR5)<< "Inconsistency in the calibrated vector"<<std::endl;
        return 0;
    }
    return _signals[imember]*calv[imember];
}

float AlibavaCluster::getNoise(int imember)
{
    if(imember >= getClusterSize())
    {
        streamlog_out(ERROR5)<< "Not enough members in this AlibavaCluster"<<std::endl;
        return 0;
    }
    return _noise[imember];
}

float AlibavaCluster::getNoise()
{
    float totalnoise=0;
    for(int imember=0; imember<getClusterSize(); imember++)
    {
        totalnoise += _noise[imember];
    }
    return totalnoise;
}

float AlibavaCluster::getTotalSNR(FloatVec noiseVec) 
{
    float snr_total = 0;
    for(int imember=0; imember<getClusterSize(); imember++) 
    {
        int channel = getChanNum(imember);
        if(noiseVec[channel]==0) 
        {
            streamlog_out (ERROR5)<<"Channels noise is zero! Check this cluster! Returning zero!"<<endl;
            return 0;
        }
	snr_total += getSignal(imember)/noiseVec[channel];
    }
    return snr_total * getSignalPolarity();
}

float AlibavaCluster::getTotalSNR()
{
    return this->getTotalSNR(this->_noise);
}

float AlibavaCluster::getTotalSignal() 
{
    float totalsignal=0;
    for(int imember=0; imember<getClusterSize(); imember++)
    {
        totalsignal += _signals[imember];
    }
    return totalsignal;
}

// Overload for the electrons
float AlibavaCluster::getTotalSignal(const std::vector<float> calv)
{ 
    // Check the validity of calv
    if(calv.size() != ALIBAVA::NOOFCHANNELS)
    {
        streamlog_out(ERROR5)<< "Inconsistency in the calibrated vector"<<std::endl;
        return 0;
    }
 
    float totalsignal=0;
    for(int imember=0; imember<getClusterSize(); imember++)
    {
        totalsignal += _signals[imember]*calv[imember];
    }
    return totalsignal;
}

void AlibavaCluster::add(int achannum, float asignal,float noise)
{
    _channums.push_back(achannum);
    _signals.push_back(asignal);
    _noise.push_back(noise);
}

int AlibavaCluster::getClusterSize()
{
    if(_channums.size() != _signals.size()) 
    {
        streamlog_out (ERROR5)<<"There is something wrong with "
            <<"this AlibavaCluster, number of channels is not "
            << "consistent with number of signals, registered "
            << "in this cluster!"<<std::endl;
    }
    return int (_channums.size());
}

bool AlibavaCluster::has_seed()
{
    if(_seedChanNum != -1)
    {
        return true;
    }
    return false;
}

///////////////////////
// Setters - Getters //
///////////////////////

// setter / getter for _chipNum
int AlibavaCluster::getChipNum()
{
    return _chipNum;
}

void AlibavaCluster::setChipNum(int chipnum)
{
    _chipNum = chipnum;
}

// setter / getter for _seedChanNum
int AlibavaCluster::getSeedChanNum()
{
    if (!has_seed()) 
    {
        streamlog_out(ERROR3)<<"Seed channel not defined!"<<std::endl;
    }
    return _seedChanNum;
}
void AlibavaCluster::setSeedChanNum(int seedChanNum)
{
    _seedChanNum = seedChanNum;
}

// setter / getter for _clusterID
int AlibavaCluster::getClusterID()
{
    if(_clusterID<0) 
    {
        streamlog_out(ERROR3)<<"Cluster ID is not defined!"<<std::endl;
    }
    return _clusterID;
}

void AlibavaCluster::setClusterID(int clusterID)
{
    _clusterID = clusterID;
}

// setter / getter for _isSensitiveAxisX
void AlibavaCluster::setIsSensitiveAxisX(bool isSensitiveAxisX)
{
    _isSensitiveAxisX = isSensitiveAxisX;
}

bool AlibavaCluster::getIsSensitiveAxisX()
{
    return _isSensitiveAxisX;
}
/*
// setter / getter for _sensorIDOffset
int AlibavaCluster::getSensorIDOffset(){
	return _sensorIDOffset;
}
void AlibavaCluster::setSensorIDOffset(int sensorIDOffset){
	_sensorIDOffset = sensorIDOffset;
}

// setter / getter for _missingCorrdinateValue
int AlibavaCluster::getMissingCorrdinateValue(){
	return _missingCorrdinateValue;
}
void AlibavaCluster::setMissingCorrdinateValue(int missingCorrdinateValue){
	_missingCorrdinateValue = missingCorrdinateValue;
}
*/
// setter / getter for _signalPolarity
int AlibavaCluster::getSignalPolarity()
{
    return _signalPolarity;
}

void AlibavaCluster::setSignalPolarity(int signalPolarity)
{
    _signalPolarity = signalPolarity;
}

void AlibavaCluster::print()
{
    streamlog_out(MESSAGE1)<<"************ Cluster ************"<<endl;
    streamlog_out(MESSAGE1)<<"ClusterID: "<< getClusterID()<< " ChipNum: "<<getChipNum()<<std::endl;
    
    std::string axis;
    if(getIsSensitiveAxisX()) 
    {
        axis = string("X");
    }
    else
    {
        axis = string("Y");
    }
    streamlog_out(MESSAGE1)<<"SensitiveAxis: "<< axis << " SignalPolarity: "<<getSignalPolarity()<<std::endl;
	
    streamlog_out(MESSAGE1)<<"SeedChannel: "<< getSeedChanNum() << " ClusterSize: "
        <<getClusterSize()<< " Eta: " <<getEta() <<std::endl;
    streamlog_out(MESSAGE1)<<"ClusterMembers: channelNum, signal"<< std::endl;
    streamlog_out(MESSAGE1)<<"               ";
    for(int imember=0; imember<getClusterSize(); ++imember)
    {
        streamlog_out(MESSAGE1)<<"["<<getChanNum(imember)<<", "<<getSignal(imember)<< "] ";
    }
    streamlog_out(MESSAGE1)<<std::endl;
    streamlog_out(MESSAGE1)<<"***********************************"<<std::endl;
}









































