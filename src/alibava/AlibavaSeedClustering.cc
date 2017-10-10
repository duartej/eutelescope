/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  Modified by J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 *    - Modify cluster algorithm (improve clearity and
 *      performance)
 */

// alibava includes ".h"
#include "AlibavaSeedClustering.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaCluster.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <EVENT/LCFloatVec.h>
#include <IMPL/LCCollectionVec.h> 
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>


// ROOT includes ".h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <algorithm>
#include <functional>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaSeedClustering::AlibavaSeedClustering () :
    AlibavaBaseProcessor("AlibavaSeedClustering"),
    _seedCut(3),
    _neighCut(2),
    _isSensitiveAxisX(true),
    _signalPolarity(-1),
    _etaHistoName("hEta"),
    _clusterSizeHistoName("hClusterSize"),
    _clusterSizePerEvtHistoName("hClusterPerEvent"),
    _clusterChargePerTDCTimeHistoName("hClusterChargePerTDCTime"),
    _xtLCIOFile("crosstalkfactors.slcio"),
    _nNeighbourgs(8),
    _neighbourgsHistoName("hSeedNeighbourgs")
{
    // modify processor description
    _description ="AlibavaSeedClustering finds clusters using seed "\
                   "and neighbour cuts ";
    // first of register the input /output collection
    registerInputCollection(LCIO::TRACKERDATA, "InputCollectionName",
            "Input collection name, it should be pedestal subtracted",
            _inputCollectionName, std::string("recodata") );
    registerOutputCollection(LCIO::TRACKERDATA, "OutputCollectionName",
            "Output data collection name",
            _outputCollectionName, std::string("alibava_clusters") );
	
    // if needed one can change these to optional parameters
    registerProcessorParameter ("NoiseInputFile",
            "The filename where the pedestal and noise values stored",
            _pedestalFile , string("pedestal.slcio"));
    registerProcessorParameter ("NoiseCollectionName",
            "Noise collection name, better not to change",
            _noiseCollectionName, std::string ("noise"));
    
    registerProcessorParameter ("SeedSNRCut",
            "The signal/noise ratio that channels have to pass to be considered as seed channel",
            _seedCut, float (3));
    registerProcessorParameter ("NeighbourSNRCut",
            "The signal/noise ratio that neigbour channels have to pass to be added to the cluster",
            _neighCut, float (2));
    registerProcessorParameter ("IsSensitiveAxisX",
            "The default sensitive axis of the strip sensor(s) according to "\
            " telescope is 'X'. If sensitive axis is Y then set this parameter "\
            "to false.",_isSensitiveAxisX, bool(true));
    registerProcessorParameter ("SignalPolarity",
            "Polarity of the signal. Set this parameter to -1 for negative "\
            "signals, any other value will be disregarded and the signal "\
            "will be assumed to be positive ",
            _signalPolarity, int (-1));
    
    registerProcessorParameter ("xtLCIOFile",
            "Name of the LCIO file where is going to be stored the"\
            " cross-talk factors",_xtLCIOFile,std::string("crosstalk_factors.slcio"));
}


void AlibavaSeedClustering::init () 
{
    streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    // this method is called only once even when the rewind is active
    /* To set of channels to be used
	 ex.The format should be like $ChipNumber:StartChannel-EndChannel$
	 ex. $0:5-20$ $0:30-100$ $1:50-70$ means from chip 0 channels
         between 5-20 and 30-100, from chip 1 channels between 50-70 will 
         be used (all numbers included). the rest will be masked and not used
	 Note that the numbers should be in ascending order and there should
         be no space between two $ character
	 */
    if(Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
    {
        Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "
            << ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << std::endl;
    }
    
    /* To choose if processor should skip masked events
     ex. Set the value to 0 for false, to 1 for true
     */
    if(Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
    {
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "
            << ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << std::endl;
    }
    // Whether to activate or not the noisy channel auto masking
    streamlog_out(MESSAGE4) << "Noisy channel auto-masking is ";
    if(Global::parameters->isParameterSet(ALIBAVA::AUTOMASKINGACTIVE))
    {
        _isAutoMaskingActive = Global::parameters->getIntVal(ALIBAVA::AUTOMASKINGACTIVE);
    }
    streamlog_out(MESSAGE4) << _isAutoMaskingActive;

    if(Global::parameters->isParameterSet(ALIBAVA::AUTOMASKINGCRITERIA))
    {
        _autoMaskingCriterium = Global::parameters->getFloatVal(ALIBAVA::AUTOMASKINGCRITERIA);
    }
    streamlog_out(MESSAGE4) << " ( if ON, using " << _autoMaskingCriterium 
        << " sigmas )" << std::endl;
    
    // check signal Polarity
    if(_signalPolarity != -1)
    {
        _signalPolarity = 1;
    }
    // usually a good idea to
    printParameters ();	
}

void AlibavaSeedClustering::processRunHeader (LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
    
    // Add processor name to the runheader
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    arunHeader->addProcessor(type());
    
    // get and set selected chips
    this->setChipSelection(arunHeader->getChipSelection());
    
    // set pedestal and noise values
    this->setPedestals();
    
    // set channels to be used (if it is defined)
    this->setChannelsToBeUsed();	
	
    // if you want
    this->bookHistos();
	
    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}


void AlibavaSeedClustering::processEvent (LCEvent * anEvent) 
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        _numberOfSkippedEvents++;
        return;
    }
    const float tdctime = alibavaEvent->getEventTime();
    
    // As an input collection we get pedestal subtracted alibava data
    // This collection contains AlibavaEvents
    LCCollectionVec * inputColVec=nullptr;
    try
    {
        inputColVec=dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(getInputCollectionName()));
    }
    catch( lcio::DataNotAvailableException ) 
    {
        // do nothing again
        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
        return;
    }
    // The Alibava Cluster collection
    LCCollectionVec * clusterColVec = new LCCollectionVec(LCIO::TRACKERDATA);
    // cell id encode for AlibavaCluster
    CellIDEncoder<TrackerDataImpl> clusterIDEncoder(ALIBAVA::ALIBAVACLUSTER_ENCODE,clusterColVec);
    	
    // Cross-check
    if( inputColVec->getNumberOfElements() != static_cast<int>(this->getChipSelection().size()) )
    {
        streamlog_out(ERROR) << "Inconsistency in the Input Collection '" << 
            getInputCollectionName() << "'. Must contain the same number of elements than " 
            << " chips evaluated (currently " << this->getChipSelection().size() << std::endl;
        return;
    }

    for(int i = 0; i < inputColVec->getNumberOfElements(); ++i)
    {
        // Note that the i-index is referring to the input col element (which is related with 
        // the chip number but is not exactly that. The collection elements must follow the same 
        // chip order stored in the getChipSelection
        const int noOfChip = this->getChipSelection()[i];

        // get the data from the collection 
        TrackerDataImpl * trkdata = dynamic_cast<TrackerDataImpl*>(inputColVec->getElementAt(i));
        // Get the clusters of hits (see findClusters function)
        std::vector<AlibavaCluster> clusters = this->findClusters(trkdata);
        // The number of clusters per event histogram
        TH1I * h = dynamic_cast<TH1I*>(_rootObjectMap[getHistoNameForChip(_clusterSizePerEvtHistoName,noOfChip)]);
        h->Fill(clusters.size());

        // The cluster charge per tdc time
        TH2F * hClusterSignalvsTime = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_clusterChargePerTDCTimeHistoName,noOfChip)]);
        // loop over clusters
        for(auto acluster: clusters) 
        {
            // Fill charge vs. time histogram
            hClusterSignalvsTime->Fill(tdctime,acluster.getSignalPolarity()*acluster.getTotalSignal());

            // create a TrackerDataImpl for each cluster
	    TrackerDataImpl * alibavaCluster = new TrackerDataImpl();
            acluster.createTrackerData(alibavaCluster);
		
            // now store chip number, seed channel, cluster ID, cluster size, sensitive axis and signal polarity in CellIDEncode
            // set cluster ID
            clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERID] = acluster.getClusterID();
            // set sensitive axis
            clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX] = int(acluster.getIsSensitiveAxisX());
	    // set signal polarity
            clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE] = int(acluster.getSignalPolarity() == -1);
	    // set chip number
            clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CHIPNUM] = acluster.getChipNum();
            // set cluster size
            clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERSIZE] = acluster.getClusterSize();
            // set seed channel number
            clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_SEED] = acluster.getSeedChanNum();
            clusterIDEncoder.setCellID(alibavaCluster);
            // And fill it to perisistify
            clusterColVec->push_back(alibavaCluster);
        } // end of loop over clusters
    } // end of loop ever detectors
    alibavaEvent->addCollection(clusterColVec, getOutputCollectionName());		
}

std::vector<AlibavaCluster> AlibavaSeedClustering::findClusters(TrackerDataImpl * trkdata)
{	
    // first get chip number
    const int chipnum = this->getChipNum(trkdata);
  
    // then get the data vector
    const FloatVec dataVec = trkdata->getChargeValues();
    
    // we will need noise vector too
    const FloatVec noiseVec = this->getNoiseOfChip(chipnum);
    
    // then check which channels we can add to a cluster
    // obviously not the ones masked. 
    // Note that all channels are set to true initially
    std::vector<bool> channel_can_be_used(dataVec.size(),true);
    // And mask channels that cannot pass NeighbourSNRCut, 
    // add the *channel numbers* as seed candidate if it pass SeedSNRCut
    std::vector<int> seedCandidates;
    for(int ichan=0; ichan<int(dataVec.size()); ++ichan)
    {
        // Not use this channel if is already masked
        if(this->isMasked(chipnum,ichan))
        {
            channel_can_be_used[ichan]=false;
            continue;
        }
        // if it is here, it is not masked, so ..
        //channel_can_be_used.push_back(true);
	// calculate snr = signal/noise
        const float snr = (_signalPolarity * dataVec[ichan])/noiseVec[ichan];
	// mask channels that cannot pass NeighbourSNRCut
        // therefore, we don't need to check it in the next loop
        if(snr < _neighCut)
        {
            channel_can_be_used[ichan]=false;
        }
        else if(snr > _seedCut) 
        {
            seedCandidates.push_back(ichan);
        }
    }
    // sort seed channels according to their SNR, highest comes first!
    std::sort(seedCandidates.begin(),seedCandidates.end(),
            [&] (const int & ichanLeft,const int & ichanRight) 
            { return ((_signalPolarity*dataVec[ichanLeft])/noiseVec[ichanLeft] > 
                (_signalPolarity*dataVec[ichanRight])/noiseVec[ichanRight]); });
    
    // Define some functors to be used to provide order to the cluster
    // finding loop (first left, then right)
    // -- the pseudo-iterator (with respect the central seed)
    auto low_strip_functor = [] (int & _index) { --_index; };
    auto high_strip_functor = [] (int & _index) { ++_index; };
    std::vector<std::function<void (int&)> > 
        next_strip_functor{ low_strip_functor, high_strip_functor };
    
    // -- the definition of edge of the chip
    auto low_isEdge_functor = [&] (const int & _index) { return (_index < 0); }; 
    auto high_isEdge_functor = [&] (const int & _index) { return (_index >= int(channel_can_be_used.size())); }; 
    std::vector<std::function<bool (const int &)> >
        isEdge_functor{ low_isEdge_functor, high_isEdge_functor };
    // -- and the initial neighbour to the central seed. This 
    //    vector must be used as: seedChan+initial_neighbour[k]
    const std::vector<int> initial_neighbour = { -1, 1 };
        
    // Get the pointer to the histo
    TH2F * hSeedNeighbours = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_neighbourgsHistoName,chipnum)]);

    // now form clusters starting from the seed channel 
    // that has highest SNR
    int clusterID = 0;
    // form clusters and store them in a vector
    std::vector<AlibavaCluster> clusterVector;
    for(const int & seedChan: seedCandidates)
    {
        // if this seed channel was used in another cluster or
        // is masked, skip it
	if(!channel_can_be_used[seedChan])
        {
            continue;
        }
        // Fill the cluster data
        AlibavaCluster acluster;
        acluster.setChipNum(chipnum);
        acluster.setSeedChanNum(seedChan);
        // [JDC] Here there is a problem: calculate eta BEFORE! 
        //       But the algorithm for calculate the cluster is
        //       made afterwards... It seems unlogical
        //       If the calculation is done now, it's going to take into 
        //       account the possibility of charge sharing between
        //       adjacents strips (even if they are not included 
        //       in the cluster)
        acluster.setEta( this->calculateEta(trkdata,seedChan) );
        acluster.setIsSensitiveAxisX(_isSensitiveAxisX);
        acluster.setSignalPolarity(_signalPolarity);
        acluster.setClusterID(clusterID);
        ++clusterID;

	// add seed channel to the cluster
        acluster.add(seedChan, dataVec[seedChan]);
        // mask seed channel so no other cluster can use it!
        channel_can_be_used[seedChan]=false;
        
        // We will check if one of the neighbours not bonded. 
        // If bonded, not use the entire cluster
        bool thereIsNonBondedChan = false;

        // Search for clusters, (left direction, right direction)
        for(unsigned int k=0; k<next_strip_functor.size(); ++k)
        {
            // Find the initial neighbour 
            const int ichanInitial=seedChan+initial_neighbour[k];

            // start the inclusion of neighbours
            for(int ichan=ichanInitial;  ; next_strip_functor[k](ichan))
            {
                // first, check if the strip is in the edge
                if(isEdge_functor[k](ichan))
                {
                    // The neighbour is not bonded [JDC??]
                    thereIsNonBondedChan = true;
                    // break the loop
                    break;
                }
                // or if the channel is masked
                if(isMasked(chipnum,ichan))
                {
                    // The neighbour is not bonded [JDC??]
                    thereIsNonBondedChan=true;
                    break;
                }
                // Fill the neighbour histo (left < 0; right >0)
                //hSeedNeighbours->Fill(ichan-seedChan,dataVec[ichan]/dataVec[seedChan]);

                // And the channel if possible
                if(channel_can_be_used[ichan])
                {
                    acluster.add(ichan, dataVec[ichan]);
                    // and mask it so that it will not be added to any other cluster
                    channel_can_be_used[ichan]=false;
                }
                else
                {
                    // if it is not possible to add it, the cluster is over
                    break;
                }
            }
        }
	
	// now if there is no neighbour not bonded
	if(thereIsNonBondedChan == false)
        {
            // fill the histograms and add them to the cluster vector
            this->fillHistos(acluster);
            clusterVector.push_back(acluster);
        }

        // Cross-talk diagnosis plot
        if(acluster.getClusterSize() == 1)
        {
            for(int ineigh=1; ineigh < _nNeighbourgs && ineigh > 0 && ineigh < ALIBAVA::NOOFCHANNELS; ++ineigh)
            {
                // Be sure the neighbourgs are in equal conditions 
                // (no seed next-to-neighbourg)
                if(std::find_if(seedCandidates.begin(),seedCandidates.end(), 
                            [&seedChan,&ineigh] (const int & ch) { return (ch == seedChan-ineigh-1 || ch == seedChan+ineigh+1); }) != seedCandidates.end()) 
                {
                   // Found another cluster/seed too close
                   break;
                }
                hSeedNeighbours->Fill(ineigh,(dataVec[seedChan-ineigh]-dataVec[seedChan+ineigh])/dataVec[seedChan]);
            }
        }
    }
    
    return clusterVector;
}

float AlibavaSeedClustering::calculateEta(TrackerDataImpl *trkdata, int seedChan)
{
    // first get chip number
    const int chipnum = getChipNum(trkdata);
    // then get the data vector
    EVENT::FloatVec dataVec = trkdata->getChargeValues();

    // we will multiply all signal values by _signalPolarity to work on positive signal always
    float seedSignal = _signalPolarity * dataVec.at(seedChan);
    
    // now we make an unrealistic signal
    float unrealisticSignal = -10000; // you cannot get this signal from alibava
	
    const int leftChan = seedChan - 1;
    float leftSignal = unrealisticSignal;
    // check if the channel on the left is masked
    if( leftChan >= 0 && isMasked(chipnum,leftChan)==false ) 
    {
        leftSignal = _signalPolarity * dataVec.at(leftChan);
    }
    
    const int rightChan = seedChan+1;
    float rightSignal = unrealisticSignal;
    // check if the channel on the right is masked
    if( rightChan < int( dataVec.size() ) && isMasked(chipnum, rightChan) == false ) 
    {
        rightSignal = _signalPolarity * dataVec.at(rightChan);
    }
	
    // now compare right and left channel to see which one is higher

    // if both right anf left channel is masked. Simply return -1
    // this case should not be saved by clustering algorithm anyways
    if(rightSignal == unrealisticSignal && leftSignal == unrealisticSignal ) 
    {
        streamlog_out (DEBUG1) << "Both neighbours are masked!"<<endl;
        return -1;
    }
	
    float eta = -1;
    // compare left and right signals
    // here both signal has to be positive (if not noise)
    // if one of the channel is masked, since it will be 
    // set to unrealisticSignal other signal will always 
    // be higher then the masked one.
    // Eta calculation: chargeOnLeftChannel / (chargeOnLeftChannel + chargeOnRightChannel)
    if( leftSignal > rightSignal) 
    {
        // then seed channel is on the right
	eta = leftSignal / ( leftSignal + seedSignal );
    }
    else 
    {
        // seed channel is on the left
	eta = seedSignal / (seedSignal + rightSignal);
    }
    return eta;	
}


void AlibavaSeedClustering::check (LCEvent * /* evt */ ) 
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}


void AlibavaSeedClustering::end() 
{
    if (_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents
            <<" events skipped since they are masked" << std::endl;
    }

    // Set in the proper directory, we are going to 
    // store some new histograms and Graphs
    AIDAProcessor::tree(this)->cd(this->name());

    // Prepare the averaged of the relative charge difference 
    // between the adjacent channels to a seed
    TCanvas *canvas = new TCanvas("aux_canvas");
    TGraphErrors * gr = new TGraphErrors();
    gr->SetTitle(";Neighbourg distance; <(Q_{seed-1}-Q_{seed+1})/Q_{seed}>");

    // The collection to be filled with the cross-talk factors
    LCCollectionVec * constantsCollection = new LCCollectionVec(LCIO::LCFLOATVEC);
    // The LCIO vector to be stored
    EVENT::LCFloatVec * xtconstant = new EVENT::LCFloatVec();
    // Evaluate the cross-talk factors
    for(const auto & chipnum: this->getChipSelection())
    {        
        // Get the pointer to the histo 
        TH2F * hSeedNeighbours = dynamic_cast<TH2F*>(_rootObjectMap[getHistoNameForChip(_neighbourgsHistoName,chipnum)]);
        int igr = 0;
        for(int xbin = 2; xbin < hSeedNeighbours->GetNbinsX()+1; ++xbin, ++igr)
        {
            const int sd = hSeedNeighbours->GetXaxis()->GetBinCenter(xbin);
            const std::string seed_distance(std::to_string(static_cast<int>(sd)));
            // XXX: Just using the two first neighbourgs
            /*if( seed_distance > 2 )
            {
                break;
            }*/
            // Get the projection histogram
            const std::string projection_name(std::string(hSeedNeighbours->GetName())+"_Neighbourg_"+seed_distance);
            TH1D * hproj = hSeedNeighbours->ProjectionY(projection_name.c_str(),xbin,xbin);
            
            // A first gaussian estimation
            TF1 * gaus = new TF1(std::string(projection_name+"_gaus").c_str(),"gaus",-0.07,0.07);
            hproj->Fit(gaus,"QNR");
            const std::vector<double> gaus_parameters { gaus->GetParameter(0),gaus->GetParameter(1), gaus->GetParameter(2) };

            // A fine model, adding a polynomial for the continous background
            TF1 * mod = new TF1(std::string(projection_name+"_gaus_cheb").c_str(),"gaus(0)+cheb2(3)");
            // - Initialize the gausian parameters, the polynomial is going to be easy to fit
            for(int k=0; k < gaus->GetNumberFreeParameters(); ++k)
            {
                // Note that the gausian in the model is set the first one,
                // therefore can do that
                mod->SetParameter(k,gaus->GetParameter(k));
            }
            // Now fit the full model
            hproj->Fit(mod,"Q");

            // Update the mean graph
            gr->SetPoint(igr,sd,mod->GetParameter(1));
            gr->SetPointError(igr,0.0,mod->GetParameter(2)/std::sqrt(hproj->GetEntries()));

            // Store it in the DB
            xtconstant->push_back(mod->GetParameter(1));
            
            // Delete histos and graphs
            delete gaus;
            gaus = nullptr;
            delete mod;
            mod = nullptr;
        }
    }
    // Store the graph
    gr->Write("gRelSeedChargeMean");
    delete gr;
    gr = nullptr;
    delete canvas;
    canvas =nullptr;
    // Store the DB
    constantsCollection->push_back(xtconstant);
    
    // Creating the DB to store the cross-talk coefficients
    LCWriter* lcWriter = LCFactory::getInstance()->createLCWriter();
    try 
    {
        lcWriter->open( _xtLCIOFile,LCIO::WRITE_NEW);
    } 
    catch(IOException& e) 
    {
        streamlog_out ( ERROR4 ) << e.what() << std::endl;
        exit(-1);
    }
    streamlog_out ( MESSAGE5 ) << "Cross-talk factors being written to " << _xtLCIOFile << endl;
    
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
        
    event->addCollection(constantsCollection, "xtfactors");
    lcWriter->writeEvent(event);
    delete event;
    lcWriter->close();


    streamlog_out ( MESSAGE4 ) << "Successfully finished" << std::endl;	
}

void AlibavaSeedClustering::fillHistos(AlibavaCluster anAlibavaCluster)
{
    const int ichip = anAlibavaCluster.getChipNum();
    
    float eta = anAlibavaCluster.getEta();
    TH1F * hEta = dynamic_cast<TH1F*> (_rootObjectMap[getHistoNameForChip(_etaHistoName,ichip)]);
    hEta->Fill(eta);
    
    int clusterSize = anAlibavaCluster.getClusterSize();
    TH1F * hClusterSize = dynamic_cast<TH1F*> (_rootObjectMap[getHistoNameForChip(_clusterSizeHistoName,ichip)]);
    hClusterSize->Fill( clusterSize );
}

void AlibavaSeedClustering::bookHistos()
{
    AIDAProcessor::tree(this)->cd(this->name());
    
    for(const auto & ichip: this->getChipSelection()) 
    {
        // Clustersize histogram
	std::string histoName(getHistoNameForChip(_clusterSizeHistoName,ichip));
        std::string title("Cluster Size (chip "+to_string(ichip)+string(");Cluster Size;Entries"));
        TH1F * hClusterSize = new TH1F(histoName.c_str(),title.c_str(),10, 0, 10);
        _rootObjectMap.insert(make_pair(histoName, hClusterSize));
        
        // Clustersize per event
	histoName=getHistoNameForChip(_clusterSizePerEvtHistoName,ichip);
        title = "Number of cluster per Event (chip "+to_string(ichip)+string(");Cluster Size;Entries");
        TH1I * hClusterSizePerEvt = new TH1I(histoName.c_str(),title.c_str(),10, 0, 10);
        _rootObjectMap.insert(make_pair(histoName, hClusterSizePerEvt));
        
        // Eta histogram ClusterSize > 1
        histoName = getHistoNameForChip(_etaHistoName,ichip);
	title = std::string("Eta distribution ClusterSize > 1 (chip "+to_string(ichip)+string(");Eta;Number of Entries"));
        TH1F * hEta = new TH1F (histoName.c_str(),title.c_str(),100, -0.1, 1.1);
	_rootObjectMap.insert(make_pair(histoName, hEta));
        
        // Cluster charge per TDC time 
	histoName=getHistoNameForChip(_clusterChargePerTDCTimeHistoName,ichip);
        title = "Charge cluster vs TDC time (chip "+to_string(ichip)+string(");TDC time [ns];Cluster Charge [ADC];");
        TH2F * hClusterChargePerTDCtime = new TH2F(histoName.c_str(),title.c_str(),100, 0, 30,1000,0,1000);
        _rootObjectMap.insert(make_pair(histoName, hClusterChargePerTDCtime));
	
        histoName=getHistoNameForChip(_neighbourgsHistoName,ichip);
        title = "Seed neighbourgs (chip "+to_string(ichip)+string(");Neighbourg distance to Seed; Relative charge to seed;");
        TH2F * hSeedNeighbourg = new TH2F(histoName.c_str(),title.c_str(),_nNeighbourgs, -0.5, _nNeighbourgs-0.5,200,-0.4,0.4);
        _rootObjectMap.insert(make_pair(histoName, hSeedNeighbourg));

    } // end of loop over selected chips
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}

std::string AlibavaSeedClustering::getHistoNameForChip(string histoName, int ichip)
{
    return std::string(histoName+"_chip"+std::to_string(ichip));
}




