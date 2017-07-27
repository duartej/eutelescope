/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - (Partially yet) Clean and coding style (Allman)
 *
 */


// alibava includes ".h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "AlibavaBaseProcessor.h"
#include "AlibavaPedNoiCalIOManager.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// REALLY??
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

// ROOT includes ".h"
#include "TObject.h"

// system includes <>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <map>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaBaseProcessor::AlibavaBaseProcessor (std::string processorName) :
    Processor(processorName),
    _rootObjectMap(),
    _inputCollectionName(ALIBAVA::NOTSET),
    _outputCollectionName(ALIBAVA::NOTSET),
    _pedestalFile(ALIBAVA::NOTSET),
    _calibrationFile(ALIBAVA::NOTSET),
    _pedestalCollectionName(ALIBAVA::NOTSET),
    _noiseCollectionName(ALIBAVA::NOTSET),
    _chargeCalCollectionName(ALIBAVA::NOTSET),
    _channelsToBeUsed(),
    _isAutoMaskingActive(false),
    _autoMaskingCriterium(2.5),
    _skipMaskedEvents(false),
    _numberOfSkippedEvents(0),
    _chipSelection(),
    _pedestalMap(),
    _noiseMap(),
    _chargeCalMap(),
    _isPedestalValid(false),
    _isNoiseValid(false),
    _isCalibrationValid(false)
{
    // modify processor description
    _description = "AlibavaBaseProcessor";
    // reset all the final arrays
    setAllMasksTo(false);
}

// checks if the root object exists in _rootObjectMap
bool AlibavaBaseProcessor::doesRootObjectExists(std::string aHistoName){
	map<string,TObject*>::const_iterator it = _rootObjectMap.find(aHistoName);
	return it!=_rootObjectMap.end();
}


///////////////////////////
// Input/Output Collection
///////////////////////////

// getter and setter for _inputCollectionName
void AlibavaBaseProcessor::setInputCollectionName(std::string inputCollectionName){
	_inputCollectionName = inputCollectionName;
}
std::string AlibavaBaseProcessor::getInputCollectionName(){
	return _inputCollectionName;
}

// getter and setter for _outputCollectionName
void AlibavaBaseProcessor::setOutputCollectionName(std::string outputCollectionName){
	_outputCollectionName = outputCollectionName;
}
std::string AlibavaBaseProcessor::getOutputCollectionName(){
	return _outputCollectionName;
}


///////////////////////////
// Pedestal and Noise
///////////////////////////

// used to set pedestal and noise values
// Note that chipSelection has to be set first!!!
void AlibavaBaseProcessor::setPedestals()
{
    // first get selected chips
    EVENT::IntVec selectedchips = getChipSelection();
    if(selectedchips.size()==0)
    {
        streamlog_out(ERROR5)<< "No selected chips found! Couldn't set the pedestal, noise values!"<<endl;
    }
    AlibavaPedNoiCalIOManager man;
    EVENT::FloatVec vec_float;
    
    // for each selected chip, get and save pedestal and noise values
    for(unsigned int ichip=0; ichip<selectedchips.size(); ichip++) 
    {
        int chipnum = selectedchips[ichip];
        // if pedestalCollectionName set
        if(getPedestalCollectionName()!= string(ALIBAVA::NOTSET)) 
        {
            // get pedestal for this chip
	     vec_float.clear();
             vec_float = man.getPedNoiCalForChip(_pedestalFile,_pedestalCollectionName, chipnum);
             _pedestalMap.insert(make_pair(chipnum, vec_float));
        }
        else
        {
            streamlog_out(DEBUG5)<< "The pedestal values for chip "
                <<chipnum<<" is not set, since pedestalCollectionName is not set!"<<std::endl;
        }
		
	// if noiseCollectionName set
        if(getNoiseCollectionName()!= string(ALIBAVA::NOTSET))
        {
            // get noise for this chip
            vec_float.clear();
            vec_float = man.getPedNoiCalForChip(_pedestalFile,_noiseCollectionName, chipnum);
            _noiseMap.insert(make_pair(chipnum, vec_float));
        }
        else
        {
            streamlog_out(DEBUG5)<< "The noise values for chip "<<chipnum
                <<" is not set, since noiseCollectionName is not set!"<<std::endl;
        }
    }
    checkPedestals();
}

void AlibavaBaseProcessor::checkPedestals()
{
    _isPedestalValid = false;
    _isNoiseValid = false;
    
    // first get selected chips
    EVENT::IntVec selectedchips = getChipSelection();
    if(selectedchips.size()==0)
    {
        streamlog_out(ERROR5)<< "No selected chips found! Couldn't "
            <<"set the pedestal, noise values!"<<std::endl;
    }
    
    // for each selected chip get and save pedestal and noise values
    for(const auto & chipnum: selectedchips)
    {
        if(getPedestalCollectionName()!= std::string(ALIBAVA::NOTSET)) 
        {
            // check pedestal values for this chip
            const EVENT::FloatVec vec_float = _pedestalMap[chipnum];
            if( int(vec_float.size()) != ALIBAVA::NOOFCHANNELS)
            {
                streamlog_out(ERROR5)<< "The pedestal values for chip "
                    <<chipnum<<" is not set properly!"<<std::endl;
            }
            else
            {
                _isPedestalValid = true;
            }
        }
        if(getNoiseCollectionName()!= std::string(ALIBAVA::NOTSET))
        {
            // check noise values for this chip
	    const EVENT::FloatVec vec_float = _noiseMap[chipnum];
            if( int(vec_float.size()) != ALIBAVA::NOOFCHANNELS)
            {
                streamlog_out(ERROR5)<< "The noise values for chip "
                    <<chipnum<<" is not set properly!"<<std::endl;
            }
            else
            {
                _isNoiseValid = true;
            }
        }
    }
	
    if(_isPedestalValid)
    {
        streamlog_out(MESSAGE5)<< "The pedestal values for all selected chips are valid!"<<std::endl;
        this->printOutPedNoiCalCollection(this->_pedestalMap);
    }
    else
    {
        streamlog_out(WARNING5)<< "The pedestal values for all selected chips are not set properly!"<<std::endl;
    }
	
    if(_isNoiseValid)
    {
        streamlog_out(MESSAGE5)<< "The noise values for all selected chips are valid!"<<std::endl;
        this->printOutPedNoiCalCollection(this->_noiseMap);
    }
    else
    {
        streamlog_out(WARNING5)<< "The noise values for all selected chips are not set properly!"<<std::endl;
    }	
}

// [JDC]: disentangle the pedestals from the noise: TO Be deployed
void AlibavaBaseProcessor::setNoise()
{
    // first get selected chips
    EVENT::IntVec selectedchips = getChipSelection();
    if(getChipSelection().size()==0)
    {
        streamlog_out(ERROR5)<< "No selected chips found! Couldn't set the noise values!"<<endl;
    }
    AlibavaPedNoiCalIOManager man;
    
    // for each selected chip, get and save pedestal and noise values
    for(auto chipnum: getChipSelection()) 
    {
	// if noiseCollectionName set
        if(getNoiseCollectionName()!= string(ALIBAVA::NOTSET))
        {
            // get noise for this chip
            _noiseMap.insert(std::make_pair(chipnum, 
                        man.getPedNoiCalForChip(_pedestalFile,_noiseCollectionName, chipnum)));
        }
        else
        {
            streamlog_out(DEBUG5)<< "The noise values for chip "<<chipnum
                <<" is not set, since noiseCollectionName is not set!"<<std::endl;
        }
    }
    checkNoise();
}

void AlibavaBaseProcessor::checkNoise()
{
    _isNoiseValid = true;
    
    if(getChipSelection().size()==0)
    {
        streamlog_out(ERROR5)<< "No selected chips found! Couldn't "
            <<"set the pedestal, noise values!"<<std::endl;
    }
    
    // for each selected chip get and save pedestal and noise values
    for(const auto & chipnum: getChipSelection())
    {
        if(getNoiseCollectionName()!= std::string(ALIBAVA::NOTSET))
        {
            // check noise values for this chip
            if(int(_noiseMap[chipnum].size()) != ALIBAVA::NOOFCHANNELS)
            {
                streamlog_out(ERROR5)<< "The noise values for chip "
                    <<chipnum<<" is not set properly!"<<std::endl;
                _isNoiseValid = false;
            }
        }
    }
    
    if(_isNoiseValid)
    {
        streamlog_out(MESSAGE5)<< "The noise values for all selected chips are valid!"<<std::endl;
        this->printOutPedNoiCalCollection(this->_noiseMap);
    }
    else
    {
        streamlog_out(WARNING5)<< "The noise values for all selected chips are not set properly!"<<std::endl;
    }	
}

///////////////////////////
// Pedestal
///////////////////////////
// getter and setter for _pedestalCollectionName
void AlibavaBaseProcessor::setPedestalCollectionName(std::string pedestalCollectionName){
	_pedestalCollectionName = pedestalCollectionName;
}
std::string AlibavaBaseProcessor::getPedestalCollectionName()
{
    return _pedestalCollectionName;
}


// to access the pedestal values of a chip
EVENT::FloatVec AlibavaBaseProcessor::getPedestalOfChip(int chipnum){
	return _pedestalMap[chipnum];
}

// to access the pedestal value of a channel
float AlibavaBaseProcessor::getPedestalAtChannel(int chipnum, int channum){
	if (isPedestalValid()){
		EVENT::FloatVec vec_float = getPedestalOfChip(chipnum);
		return vec_float[channum];
	}
	else {
		streamlog_out(ERROR5)<< "The pedestal values for chip "<<chipnum<<" is not set properly!"<<endl;
		return 0;
	}
}

bool AlibavaBaseProcessor::isPedestalValid(){
	return _isPedestalValid;
}

///////////////////////////
// Noise
///////////////////////////
// getter and setter for _noiseCollectionName
void AlibavaBaseProcessor::setNoiseCollectionName(std::string noiseCollectionName){
	_noiseCollectionName = noiseCollectionName;
}
std::string AlibavaBaseProcessor::getNoiseCollectionName(){
	return _noiseCollectionName;
}

// to access the noise values of a chip
EVENT::FloatVec AlibavaBaseProcessor::getNoiseOfChip(int chipnum){
	return _noiseMap[chipnum];
}
// to access the noise value of a channel
float AlibavaBaseProcessor::getNoiseAtChannel(int chipnum, int channum){
	if (isNoiseValid()){
		EVENT::FloatVec vec_float = getNoiseOfChip(chipnum);
		return vec_float[channum];
	}
	else {
		//streamlog_out(ERROR5)<< "The noise values for chip "<<chipnum<<" is not set properly!"<<endl;
		return 0;
	}
	
}

bool AlibavaBaseProcessor::isNoiseValid(){
	return _isNoiseValid;
}

// [JDC] Find noisy channels
void AlibavaBaseProcessor::maskNoisyChannels(const int & chip,const float & criteria)
{
    streamlog_out(MESSAGE5) << "Automaticaly masking noisy channel, defined as "
        << "|Noise_i - <Noise>| > " << criteria << " sigma" << std::endl;
    if(not this->isNoiseValid())
    {
        streamlog_out(ERROR5) << " The noise values for chip " << chip 
            << " is not set properly" << std::endl;
    }
    const EVENT::FloatVec noise_vec = this->getNoiseOfChip(chip);
    // -- Map to keep only non-noisy channels
    // ichannel: noise
    std::map<int,float> non_noisy_map;
    // initialize the map
    for(int i=0; i < static_cast<int>(noise_vec.size()); ++i)
    {
        // Not taking into account the masked channels
        if(this->isMasked(chip,i))
        {
            continue;
        }
        non_noisy_map.emplace(i,noise_vec[i]);
    }
    
    // 1. Obtain the mean and the standard deviation for the noise
    float mean_noise = this->getMean(non_noisy_map);
    float stddev_noise= this->getStdDev(non_noisy_map,mean_noise);
    
    // 2. Use the mean and standard deviation to find noisy channels, 
    //    while using criteria defined by the user
    
    // Remove iteratively noisy channels until converge (the
    // map will not loose any member anymore)
    unsigned int last_vec_size = noise_vec.size();
    do
    {
        // after the check from the 'while' statement, update
        // the last know vector size
        last_vec_size = non_noisy_map.size();
        for(auto it = non_noisy_map.cbegin(); it != non_noisy_map.cend(); /* no increment*/)
        {
            if(this->isMasked(chip,it->first))
            {
                ++it;
                continue;
            }
            // Remove the noisy channels and tag it as masked channel
            if( std::abs(it->second-mean_noise)/stddev_noise > criteria )
            {
                streamlog_out(MESSAGE5) << "[AUTO-MASKED] noisy channel: " << it->first 
                    << " with noise: " << it->second << " (Mean noise all channels: " 
                    << mean_noise << ", standard deviation: " << stddev_noise << ")" << std::endl;
                this->_isMasked[chip][it->first] = true;
                non_noisy_map.erase(it++);
            }
            else
            {
                // or keep it
                ++it;
            }
        }
        // 3. Recalculate the mean with the excluded signal channels
        mean_noise  = this->getMean(non_noisy_map);
        stddev_noise= this->getStdDev(non_noisy_map,mean_noise);
    } while( non_noisy_map.size() != last_vec_size );
    
    // The new vector
    EVENT::FloatVec noise_vec_new(ALIBAVA::NOOFCHANNELS,0);
    for(auto ich_noise: non_noisy_map)
    {
        noise_vec_new[ich_noise.first] = ich_noise.second;
    }

    // and update the noise map with the masked channels (zero value)
    _noiseMap[chip] = noise_vec_new;
}

///////////////////////////
// Charge Calibration
///////////////////////////


// used to set pedestal and noise values
// Note that chipSelection has to be set first!!!
void AlibavaBaseProcessor::setCalibration()
{
    // first get selected chips
    EVENT::IntVec selectedchips = getChipSelection();
    if(selectedchips.size()==0)
    {
        streamlog_out(ERROR5)<< "No selected chips found! Couldn't set calibration values!"<<std::endl;
    }
    
    AlibavaPedNoiCalIOManager man;

    // for each selected chip get and save pedestal and noise values
    for(const auto & chipnum: selectedchips)
    {
        // get charge calibration for this chip
        EVENT::FloatVec vec_float= man.getPedNoiCalForChip(_calibrationFile,_chargeCalCollectionName, chipnum);
        _chargeCalMap.insert(std::make_pair(chipnum, vec_float));
    }
    checkCalibration();
}

void AlibavaBaseProcessor::checkCalibration()
{
    _isCalibrationValid = true;
    
    // for each selected chip get and save pedestal and noise values
    for(auto chipnum: getChipSelection())
    {
        if( static_cast<int>(_chargeCalMap[chipnum].size()) != ALIBAVA::NOOFCHANNELS )
        {
            streamlog_out(ERROR5)<< "The charge calibration values for chip "
                << chipnum <<" is not set properly!"<<std::endl;
            _isCalibrationValid = false;
        }
    }
	
    if(_isCalibrationValid)
    {
        streamlog_out(MESSAGE5)<< "The calibration values for all selected chips are valid!"<<std::endl;
        this->printOutPedNoiCalCollection(this->_chargeCalMap);
    }
    else
    {
        streamlog_out(ERROR5)<< "No selected chips found! Couldn't set calibration values!"<<std::endl;
    }
}


// getter and setter for _chargeCalCollectionName
void AlibavaBaseProcessor::setChargeCalCollectionName(std::string chargeCalCollectionName)
{
    _chargeCalCollectionName = chargeCalCollectionName;
}

std::string AlibavaBaseProcessor::getChargeCalCollectionName()
{
    return _chargeCalCollectionName;
}

// to access the charge calibration values of a chip
EVENT::FloatVec AlibavaBaseProcessor::getChargeCalOfChip(int chipnum)
{
    return _chargeCalMap[chipnum];
}

// to access the charge calibration value of a channel
float AlibavaBaseProcessor::getChargeCalAtChannel(int chipnum, int channum)
{
    if(_isCalibrationValid)
    {
        return getChargeCalOfChip(chipnum)[channum];
    }
    else 
    {
        streamlog_out(ERROR5)<< "The noise values for chip "<<chipnum<<" is not set properly!"<<endl;
        return 0;
    }
}


///////////////////////////
// Others
///////////////////////////

// [JDC] XXX TO BE DEPLOYED
// Check the _channelsToBeUsed is properly filled
// bool AlibavaBaseProcessor::isChannelsToBeUsedActive()
// /* To set of channels to be used
// ex.The format should be like $ChipNumber:StartChannel-EndChannel$
//	 ex. $0:5-20$ $0:30-100$ $1:50-70$ means from chip 0 channels
//     between 5-20 and 30-100, from chip 1 channels between 50-70 will 
//     be used (all numbers included). the rest will be masked and not used
//   Note that the numbers should be in ascending order and there should
//     be no space between two $ character
//   */
//if(Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
//{
//    Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
//    return true;
//}
//else 
//{
//    streamlog_out ( MESSAGE4 ) << "The Global Parameter "
//        << ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << std::endl;
//    return false;
//}
//

// Print-out function for noise, pedestals or calibration
void AlibavaBaseProcessor::printOutPedNoiCalCollection(const std::map<int,EVENT::FloatVec> & themap)
{
    int channelprintnum = 8;
    
    std::string wcol("");
    if(themap == this->_pedestalMap)
    {
        wcol = "Pedestals";
    }
    else if(themap == this->_noiseMap)
    {
        wcol = "Noise";
    }
    else if(themap == this->_chargeCalMap)
    {
        wcol = "Calibration";
    }
    
    streamlog_out( MESSAGE5 ) <<"******************************************************"<<endl;
    streamlog_out( MESSAGE5 ) <<"****** Current values for the " 
        << std::setw(12) << wcol << " map: ******"<<endl;
    
    std::cout.precision(3);
    for(auto ichip: _chipSelection) 
    {
        streamlog_out( MESSAGE5 ) <<"***************** Chip "<< ichip <<" ****************";
        for(int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ) 
        {
            if(ichan % channelprintnum ==0)
            {
                streamlog_out( MESSAGE5 ) <<endl;
                streamlog_out( MESSAGE5 ) <<" Channels "<< std::setw(3) 
                    << ichan<<" - "<< std::setw(3) << ichan+channelprintnum-1 
                    << " : ";
            }
            streamlog_out( MESSAGE5 ) << std::setw(8) << themap.at(ichip)[ichan] << " ";
            ++ichan;
        }
        streamlog_out( MESSAGE5 )<< endl;
    }
    streamlog_out( MESSAGE5 ) <<"******************************************************"<<endl;
}
	

// getter and setter for _chipSelection
void AlibavaBaseProcessor::setChipSelection(EVENT::IntVec chipselection){
	_chipSelection = chipselection;
}

EVENT::IntVec AlibavaBaseProcessor::getChipSelection(){
	return _chipSelection;
}

// get number of selected chips
unsigned int AlibavaBaseProcessor::getNumberOfChips(){
	return _chipSelection.size();
}

// returns true if the chip is in the list of selected chip
bool AlibavaBaseProcessor::isChipValid(int ichip){
	for (unsigned int i=0; i<_chipSelection.size(); i++)
		if (_chipSelection[i] == ichip) return true;
	
	// if it reaches here it means that the chip number is not in the selected chip numbers, then return false
	return false;
	
}
// returns true if the channel number is valid
bool AlibavaBaseProcessor::isChannelValid(int ichan){
	if ( ichan < ALIBAVA::NOOFCHANNELS )
		return true;
	else
		return false;
}

// to access the noise value of a channel
bool AlibavaBaseProcessor::isMasked(int ichip, int ichan){
	// if channels to be used not identified use all channels
	if (_channelsToBeUsed.size()==0)
		return false;
	else if (isChipValid(ichip) && isChannelValid(ichan))
		return _isMasked[ichip][ichan];
	else{
		streamlog_out( ERROR5 ) <<"Trying to access mask value of non existing chip/channel. Returning zero."<< endl;
		return 0;
	}
}

// [JDC] XXX
// More useful overloaded??
//bool AlibavaBaseProcessor::isMasked(const EVENT::FloatVec & rawdatavec

void AlibavaBaseProcessor::setChannelsToBeUsed()
{
    // FIXME:: Probably does not need this function. A better parameter could be
    //         defined instead of the _channelsToBeUsed string
    //
    // Let's decode this StringVec.
    /*  The format of _channelsToBeUsed string parameter shoul be like
     *  $0:5-20$ $0:30-100$ $1:50-70$ means from chip 0 channels between 5-20 and 30-100, from chip 1 channels between 50-70
     *  will be used (all numbers included).
     *   the rest will be masked and not used
     *   Note that the numbers should be in ascending order
     */
    streamlog_out (DEBUG5) <<"Setting channels to be used! "<<endl;
	
    // first mask all channels later we will unmask the ones selected
    setAllMasksTo(true);
    streamlog_out (DEBUG5) <<"All channels masked for now!"<<endl;
    
    // now loop over strings to decode it
    for (unsigned int istr=0; istr<_channelsToBeUsed.size(); ++istr) 
    {
        string istring = _channelsToBeUsed[istr];
        int onchip, fromchannel, tochannel;
        decodeMaskingString(istring, &onchip, &fromchannel, &tochannel);
        streamlog_out (DEBUG5) <<"Processing channel selection: "<<istring
            << " ... un-masking channels from "<<fromchannel<<" to "
            <<tochannel<<" on chip "<<onchip<< std::endl;
        if (isMaskingValid(onchip,fromchannel,tochannel)) 
        {
            // unmask the selected channels
            for(int ichan = fromchannel; ichan<tochannel+1; ++ichan)
            {
                _isMasked[onchip][ichan]=false;
            }
    		
    	}
    }
    // And now the automatic masking if requested
    if(_isAutoMaskingActive)
    {
        for(auto ichip: this->getChipSelection())
        {
            this->maskNoisyChannels(ichip,_autoMaskingCriterium);
        }
    }
    printChannelMasking();
}

void AlibavaBaseProcessor::decodeMaskingString(string istring, int *onchip, int *fromchannel, int *tochannel )
{

    // FIXME:: Probably does not need this function. A better parameter could be
    //         defined instead of the _channelsToBeUsed string
    //         See setChannelsToBeUsed function
	
    //lets first set everything to -1 to check later
    *onchip = -1;
    *fromchannel = -1;
    *tochannel = -1;
    
    size_t tmppos=0;
    string tmpstring;
    while(tmppos<istring.size()) 
    {   	
    	//first find the pos of " character
    	tmppos = istring.find("$",0);
    	tmppos = tmppos+1;
    	
    	//find the chip number
    	tmpstring = getSubStringUpToChar(istring,":",tmppos);
    	*onchip = atoi(tmpstring.c_str());
    	tmppos = tmppos + tmpstring.size()+1;
    	
    	
    	// start chan
    	tmpstring = getSubStringUpToChar(istring,"-",tmppos);
    	*fromchannel = atoi(tmpstring.c_str());
    	tmppos = tmppos + tmpstring.size()+1;
    	
    	if (tmppos >= istring.size()) break;
    	//end chan
    	tmpstring = getSubStringUpToChar(istring,"$",tmppos);
    	tmppos = tmppos + tmpstring.size()+1;
    	*tochannel = atoi(tmpstring.c_str());
    	
    }
}

bool AlibavaBaseProcessor::isMaskingValid(int onchip, int fromchannel, int tochannel ){
	
	bool isValid = true;
	
	// check that channel ranges read correctly
	// channels should be in ascending order!
	
	if(!isChipValid(onchip)){
		streamlog_out( ERROR5 ) <<" Chip selection on channel masking is not valid!"<< endl;
		isValid = false;
	}
	if (!isChannelValid(fromchannel) || !isChannelValid(tochannel)){
		streamlog_out( ERROR5 ) <<" Channel selection on channel masking is not valid!"<< endl;
		isValid = false;
	}
	// start channel cannot be bigger than end chan
	if (fromchannel>tochannel) {
		streamlog_out( ERROR5 ) <<" Channel range on channel masking is not valid! Channels should be in ascending order!"<< endl;
		isValid = false;
	}
	return isValid;
}


void AlibavaBaseProcessor::printChannelMasking()
{
    int channelprintnum = 16;
    
    streamlog_out( MESSAGE5 ) <<" ChannelsToBeUsed set as: ";
    for(unsigned int istr=0; istr<_channelsToBeUsed.size(); ++istr)
    {
        streamlog_out( MESSAGE5 ) <<_channelsToBeUsed[istr]<<", ";
    }
    streamlog_out( MESSAGE5 )<< endl;
    
    streamlog_out( MESSAGE5 ) <<"******************************************************"<<endl;
    streamlog_out( MESSAGE5 ) <<"************** Applied Channel Masking: **************"<<endl;
    
    for(auto ichip: _chipSelection) 
    {
        streamlog_out( MESSAGE5 ) <<"***************** Chip "<< ichip <<" ****************";
        for(int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ) 
        {
            if(ichan % channelprintnum ==0)
            {
                streamlog_out( MESSAGE5 ) <<endl;
                streamlog_out( MESSAGE5 ) <<" Channels "<< std::setw(3) 
                    << ichan<<" - "<< std::setw(3) << ichan+channelprintnum-1 
                    << " : ";
            }
            streamlog_out( MESSAGE5 ) << !(isMasked(ichip, ichan)) << " ";
            ++ichan;
        }
        streamlog_out( MESSAGE5 )<< endl;
    }
    streamlog_out( MESSAGE5 ) <<"******************************************************"<<endl;
}

//returns the min chip number
unsigned int AlibavaBaseProcessor::getMinChipNumber(){
	//the chipSelection should be in ascending order!
	//this is guaranteed with AlibavaConverter::checkIfChipSelectionIsValid()
	return _chipSelection.front();
}
//returns the max chip number
unsigned int AlibavaBaseProcessor::getMaxChipNumber(){
	//the chipSelection should be in ascending order!
	//this is guaranteed with AlibavaConverter::checkIfChipSelectionIsValid()
	return _chipSelection.back();
	
}


void AlibavaBaseProcessor::setAllMasksTo(bool abool){
	for (int i=0; i<ALIBAVA::NOOFCHIPS; i++)
		for (int j=0; j<ALIBAVA::NOOFCHANNELS; j++)
			_isMasked[i][j]=abool;
	
}

// only valid for AlibavaData nor for AlibavaClusters
int AlibavaBaseProcessor::getChipNum(TrackerDataImpl * trkdata){
	
	CellIDDecoder<TrackerDataImpl> chipIDDecoder(ALIBAVA::ALIBAVADATA_ENCODE);
	int chipnum = static_cast<int> ( chipIDDecoder( trkdata )[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] );
	return chipnum;
}


// Helper functions to calculate mean and std-deviation
EVENT::FloatVec AlibavaBaseProcessor::convertIntoVec(const std::map<int,float> & m)
{
    EVENT::FloatVec v;
    for(const auto & _f: m)
    {
        v.push_back(_f.second);
    }
    return v;
}

float AlibavaBaseProcessor::getMean(const std::map<int,float> & m)
{
    return this->getMean(this->convertIntoVec(m));
}

float AlibavaBaseProcessor::getMean(const EVENT::FloatVec & v)
{
    if(v.size() == 0)
    {
        return 0.0;
    }
    return std::accumulate(v.begin(),v.end(),0.0)/float(v.size());
}

float AlibavaBaseProcessor::getStdDev(const std::map<int,float> & m,const float & mean)
{
    return this->getStdDev(this->convertIntoVec(m),mean);
}

float AlibavaBaseProcessor::getStdDev(const EVENT::FloatVec & v,const float & mean)
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


