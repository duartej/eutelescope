/*
 *  Performs the calibration: converts ADC counts into
 *  charge (number of electrons) 
 *
 *  J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *
 */


// alibava includes ".h"
#include "AlibavaCalibrateProcessor.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#include "EUTelExceptions.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes ".h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"

// system includes <>
#include <string>
#include <iostream>
#include <numeric>


using namespace lcio;
using namespace marlin;
using namespace alibava;


AlibavaCalibrateProcessor::AlibavaCalibrateProcessor() :
    AlibavaBaseProcessor("AlibavaCalibrateProcessor"),
    _calibrateHistoName ("hcalibrate"),
    _temperatureHistoName("htemperature"),
    _chanDataHistoName ("Data_chan"),
    _chanDataFitName ("Fit_chan"),
    _adcmax(1000.0),
    _adcmin(0.0),
    _nbins(1000),
    _nPulses(-1),
    _initialCharge(-1),
    _finalCharge(-1),
    _deltaCharge(-1),
    _nSamplesPerPulse(-1)
{
    // modify processor description
    _description ="AlibavaCalibrateProcessor obtains the corresponding charge"\
                   " to a ADC count";
	
    // first of register the input collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
            "Input raw data collection name", _inputCollectionName, std::string("rawdata_cal") );
    registerProcessorParameter ("CalibrationOutputFile",
            "The filename to store the pedestal and noise values",_calibrationFile , std::string("outputped.slcio"));
    // now the optional parameters
    registerOptionalParameter ("ChargeCalibrationCollectionName",
            "Charge calibration collection name",_chargeCalCollectionName,std::string ("charge_calibration"));
    registerOptionalParameter ("MaxADCsCountsForHistograms",
            "The max ADCs counts which define the ranges of the histograms",
            _adcmax, float(1000.0) );
    registerOptionalParameter ("MinADCsCountsForHistograms",
            "The min ADCs counts which define the ranges of the histograms",
            _adcmin, float(0.0) );
    registerOptionalParameter ("NbinsForHistograms",
            "The number of bins for the histograms",
            _nbins, int(1000) );
}


void AlibavaCalibrateProcessor::init () 
{
    //streamlog_out ( MESSAGE4 ) << "Running init" << std::endl;
    /* To set of channels to be used
     * ex.The format should be like $ChipNumber:StartChannel-EndChannel$
     * ex. $0:5-20$ $0:30-100$ $1:50-70$
     * means from chip 0 channels between 5-20 and 30-100, from chip 1 
     * channels between 50-70 will be used (all numbers included). the 
     * rest will be masked and not used. Note that the numbers should be 
     * in ascending order and there should be no space between two $ character
     */
    if(Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
    {
        Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
    }
    else 
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << std::endl;
    }
    /* To choose if processor should skip masked events
     * ex. Set the value to 0 for false, to 1 for true
     */
    if(Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
    {
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
    }
    else
    {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS 
            <<" is not set! Masked events will be used!" << std::endl;
    }
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();
}

void AlibavaCalibrateProcessor::processRunHeader(LCRunHeader * rdr) 
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << std::endl;
    
    auto arunHeader = std::make_unique<AlibavaRunHeaderImpl>(rdr);
    // Obtain the calibration run information: 
    //    Step; InitialCharge; FinalCharge; StepSize
    const std::string header(arunHeader->getHeader()+";");
    const std::string token(";");
    std::string::size_type n = 0;
    std::string::size_type i = 0;
    std::vector<int> elements;
    while(i != header.size())
    {
        n = header.find(token,i);
        elements.push_back( std::stof(header.substr(i,n)) );
        i = n+1;
    }
    // Very precise notation, otherwise some format problem is 
    // present
    if(elements.size() != 4)
    {
        throw(lcio::Exception("Inconsistency in the Run Header of the input file."\
                    "The HEADER STRING does not have the expected format: "\
                    "\"Step\"; \"InitialCharge\"; \"FinalCharge\"; StepSize  "));
    }
    _nPulses        = elements[0];
    _initialCharge  = elements[1];
    _finalCharge    = elements[2];
    _deltaCharge    = elements[3];
    // The delta charge to be used in each step should be
    if(std::fabs(float(_deltaCharge)-(_finalCharge-_initialCharge)/float(_nPulses)) > 1e-9)
    {
        throw(lcio::Exception("Inconsistency in the Run Header of the calibration"\
                    " input file. The HEADER STRING provides a deltaCharge which does not"\
                    " fulfill deltaCharge=(finalCharge-initialCharge)/numberPulses"));
    }
    // Note that the AlibavaRunHeaderImpl::getNoOfEvents() does not return the
    // number of events, given that this information is not stored in the header
    // (see AlibavaConverter class)
    // Ask to the user to be present in the steering file
    // [JDC]-- FIXME -- POTENTIAL PROBLEM WITH THIS APPROACH
    //         Try to find a way to store the real number of events in the 
    //         AlibavaConverter, otherwise, it relies in the user to 
    //         include in the MaxRecordNumber the valid number. Just
    //         let's asumme that it should be 32*100=3200 and send a warning
    //         otherwise...
    const int totalEvts= Global::parameters->getIntVal("MaxRecordNumber");
    if(totalEvts <= 0)
    {
        throw eutelescope::InvalidParameterException("Invalid 'MaxRecordNumber'. Please check"\
                " the number of events in the calibration file, and put it in the" \
                " steering file");
    }
    if(totalEvts != 3200 )
    {
        streamlog_out(WARNING) << "'MaxRecordNumber' parameter read from the steering"\
            " file set to " << _nSamplesPerPulse << " ! Usually, the calibration runs"\
            " are set to 3200 events" << std::endl;
    }
    _nSamplesPerPulse = totalEvts/_nPulses;
    // Print-out the info
    streamlog_out(MESSAGE) << "******************************************" << std::endl;
    streamlog_out(MESSAGE) << "Parameters for the calibration" << std::endl;
    streamlog_out(MESSAGE) << "Initial charge (IC): " << _initialCharge << " [electrons]" << std::endl;
    streamlog_out(MESSAGE) << "Final   charge (FC): " << _finalCharge   << " [electrons]" << std::endl;
    streamlog_out(MESSAGE) << "Charge   steps (DC): " << _deltaCharge   << " [electrons]" << std::endl;
    streamlog_out(MESSAGE) << "Number of pulses (Np): " << _nPulses << std::endl;
    streamlog_out(MESSAGE) << "Number of Samples per step (Ns): " << _nSamplesPerPulse << std::endl;
    streamlog_out(MESSAGE) << "Total processed events (Np*Ns): "  << _nSamplesPerPulse*_nPulses << std::endl;
    streamlog_out(MESSAGE) << "*** Charge at event n: IC+DC*(Ns+n)/Ns ***" << std::endl;
    streamlog_out(MESSAGE) << "******************************************" << std::endl;
    // Calibration info stored---

    arunHeader->addProcessor(type());
    setChipSelection( arunHeader->getChipSelection() );
    // Masking channels (user input)
    setChannelsToBeUsed();
    
    // An instance of the input output manager (Pedestal and noise)
    AlibavaPedNoiCalIOManager man;
    // Prepare the lcio file with the charge calibration, using the
    // run Header of the input lcio
    man.createFile(_calibrationFile, arunHeader->lcRunHeader());
    // and book the needed histograms
    bookHistos();

    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}


void AlibavaCalibrateProcessor::processEvent(LCEvent * anEvent) 
{
    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);
    if(_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) 
    {
        ++_numberOfSkippedEvents;
        return;
    }
    
    LCCollectionVec * collectionVec = nullptr;
    try
    {
        collectionVec = dynamic_cast<LCCollectionVec*>(alibavaEvent->getCollection(getInputCollectionName()));
    } 
    catch( lcio::DataNotAvailableException ) 
    {
        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << std::endl;
    }
    
    // fill temperature histogram
    if(TH1F * temperatureHisto = dynamic_cast<TH1F*> (_rootObjectMap[_temperatureHistoName]))
    {
        temperatureHisto->Fill(alibavaEvent->getEventTemp());
    }
    
    // the number of the chips
    const unsigned int noOfDetector = collectionVec->getNumberOfElements();
    for ( size_t i = 0; i < noOfDetector; ++i )
    {
        TrackerDataImpl * trkdata = dynamic_cast<TrackerDataImpl*>( collectionVec->getElementAt(i) ) ;
        // fill the raw ADCs histos
        this->fillHistos(trkdata,alibavaEvent->getCalCharge(),anEvent->getEventNumber());
    }
}

void AlibavaCalibrateProcessor::check(LCEvent * /* evt */ ) 
{
    // nothing to check here - could be used to fill 
    // check plots in reconstruction processor
}


void AlibavaCalibrateProcessor::end() 
{
    // Used the processEvent method to fill the ADCs histograms per channel/chip
    // Now obtain the calibration curves
    
    // Just defined here to avoid overload in the loops
    std::string tempHistoName,tempFitName;
    // Just to avoid the print of the canvas creation when
    // fitting (XXX: Another way?)
    TCanvas * cc = new TCanvas("cc","cc",600,800);
    
    EVENT::IntVec chipSelection = getChipSelection();
    // for each chip 
    for(auto ichip: chipSelection) 
    {
        // Get the Graph for that chip
        const std::string graphname(getCalChargeHistoName(ichip)+"_electrons");
        TH1F * grph = dynamic_cast<TH1F*>(_rootObjectMap[graphname]);
        // and for each channel evaluate the mean number of electrons per
        // ADC 
        EVENT::FloatVec calibrationVec;
        for(int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ++ichan) 
        {
            int elec_per_adc = 0;
            // if channel is masked, set pedestal and noise to 0,
            // otherwise, evaluate the
            if(!isMasked(ichip,ichan))
            {
                // First, get the histograms of the channel and the name
                // of the gaussian to fit
                tempFitName = getChanDataFitName(ichip, ichan);
                tempHistoName = getChanDataHistoName(ichip, ichan);
                TProfile * histo = dynamic_cast<TProfile*>(_rootObjectMap[tempHistoName]);
                TF1 * tempfit = dynamic_cast<TF1*>(_rootObjectMap[tempFitName]);
                // Fit to a line 
                histo->Fit(tempfit,"Q");
                // Slope: [ADC counts/number of electrons*1e-3] --> 1.0/slope
                elec_per_adc = std::round(1.0/(tempfit->GetParameter(1)*1e-3));
                const float err = tempfit->GetParError(1)*elec_per_adc;
                const int bin = grph->FindBin(ichan);
                grph->SetBinContent(bin,elec_per_adc);
                grph->SetBinError(bin,err);
            }
            calibrationVec.push_back(elec_per_adc);
        }
        AlibavaPedNoiCalIOManager man;
        man.addToFile(_calibrationFile,_chargeCalCollectionName, ichip, calibrationVec);
    }
    delete cc;
    if(_numberOfSkippedEvents > 0)
    {
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents 
            <<" events skipped since they are masked" << std::endl;
    }
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << std::endl;
}

// --- Some methods to keep consistency for names
std::string AlibavaCalibrateProcessor::getChanDataHistoName(unsigned int ichip, unsigned int ichan)
{
    std::stringstream s;
    s<< _chanDataHistoName<<"_chip"<<ichip<<"_chan" << ichan;
    return s.str();
}

std::string AlibavaCalibrateProcessor::getChanDataFitName(unsigned int ichip, unsigned int ichan)
{
    std::stringstream s;
    s<< _chanDataFitName<<"_chip"<<ichip<<"_chan"  << ichan;
    return s.str();
}

std::string AlibavaCalibrateProcessor::getCalChargeHistoName(unsigned int ichip)
{
    return _calibrateHistoName+"_chip"+std::to_string(ichip);
}

// Fill the raw ADCs histograms
void AlibavaCalibrateProcessor::fillHistos(TrackerDataImpl * trkdata, const int & injected_pulse, const int & evt)
{
    // NOTE the charge injection pattern follows:
    //
    //     EVENT   |   STRIP    | CHARGE SIGN
    // ------------+------------+-------------
    //  evt%2 == 0 |  i%2 == 0  |  NEGATIVE
    //  evt%2 == 0 |  i%2 == 1  |  POSITIVE
    //  ======================================
    //  evt%2 == 1 |  i%2 == 0  |  POSITIVE
    //  evt%2 == 1 |  i%2 == 1  |  NEGATIVE
    //  
    //  That conditions can be summarize with
    //  (-1)**[ (EVENT%2) + ((STRIP+1)%2) ]
    // Get the injected charge (electrons) --> a function
    // Note that charge will depend on the parity of the event
    //const float injected_pulse = getInjectedCharge(evt);

    // Get the raw ADCs
    const lcio::FloatVec datavec = trkdata->getChargeValues();
    const int chipnum = getChipNum(trkdata);

    TH2F * hCh = dynamic_cast<TH2F*>(_rootObjectMap[ getCalChargeHistoName(chipnum)]);
    if( hCh == nullptr )
    {
        streamlog_out(ERROR) << "Histogram not properly booked: '" 
            << getCalChargeHistoName(chipnum) 
            << "'. The program is going to crash" << std::endl;
    }

    for(unsigned int ichan=0; ichan<datavec.size();++ichan) 
    {
        if(isMasked(chipnum, ichan))
        {
            continue;
        }
        //Per channel histogram
        hCh->Fill(ichan,datavec[ichan]);
        
        // Pulse charge vs. ADC counts
        if( TProfile * histo = dynamic_cast<TProfile*>(_rootObjectMap[getChanDataHistoName(chipnum,ichan)]) )
        {
            // Note that the injected charge will depend on event parity
            // and channel parity (see above)
            //const int sign = std::pow(-1,0.((evt+1)%2))*std::pow(-1,(ichan%2));
            const int sign = std::pow(-1,((evt%2)+(ichan+1)%2));
            histo->Fill(sign*injected_pulse*1e-3,datavec[ichan]);
        }
    }
}


void AlibavaCalibrateProcessor::bookHistos()
{
    // Calculate the number of bins for the adc related
    // and set the ranges
    if(_adcmax < _adcmin)
    {
        streamlog_out( ERROR ) << "Inconsistent definition of ranges in the histograms: "
            << "[XXX TO FiNISH THE MESSAGE]"  << std::endl;
    }
    AIDAProcessor::tree(this)->cd(this->name());
    EVENT::IntVec chipSelection = getChipSelection();
    
    //the chipSelection should be in ascending order!
    //this is guaranteed with AlibavaConverter::checkIfChipSelectionIsValid()
    
    // temperature of event
    TH1F * temperatureHisto = new TH1F(_temperatureHistoName.c_str(),"Temperature",1000,-50,50);
    _rootObjectMap.insert(std::make_pair(_temperatureHistoName,temperatureHisto));
    
    // Create the histograms
    for(unsigned int i=0; i<chipSelection.size(); ++i) 
    {
        unsigned int ichip=chipSelection[i];
        TH2F * calHisto = new TH2F(getCalChargeHistoName(ichip).c_str(),"",
                ALIBAVA::NOOFCHANNELS,-0.5,ALIBAVA::NOOFCHANNELS-0.5,_nbins,_adcmin,_adcmax);
        _rootObjectMap.insert(std::make_pair(getCalChargeHistoName(ichip), calHisto)); 
        //title string for pedestal histogram
        const std::string sp("ADC average per channel (chip "+std::to_string(ichip)+");Channel Number;Signal [ADCs]");
        calHisto->SetTitle(sp.c_str());
        // And the graph of n-electrons per ADC per each channel
        const std::string graphname(getCalChargeHistoName(ichip)+"_electrons");
        const std::string gtitle = std::string("Electrons per ADC in average per channel (chip "+\
                std::to_string(ichip)+");Channel Number; electrons per ADC count");
        TH1F * grph = new TH1F(graphname.c_str(),gtitle.c_str(),ALIBAVA::NOOFCHANNELS,-0.5,ALIBAVA::NOOFCHANNELS-0.5);
        _rootObjectMap.insert(std::make_pair(graphname,grph));
    }
    
    // create the directory to store the partial histograms/functions used
    // to evaluate the final pedestals and noise
    AIDAProcessor::tree(this)->mkdir(getInputCollectionName().c_str());
    AIDAProcessor::tree(this)->cd(getInputCollectionName().c_str());
    
    // Histograms for calibration for each channel
    // (and chip)
    for(unsigned int i=0; i<chipSelection.size(); ++i) 
    {
        const unsigned int ichip=chipSelection[i];
        for( int ichan=0; ichan<ALIBAVA::NOOFCHANNELS; ++ichan) 
        {
            if(isMasked(ichip,ichan))
            {
                continue;
            }
            const std::string tempHistoName = getChanDataHistoName(ichip,ichan);
            const std::string tempHistoTitle(tempHistoName+";Injected pulse [e x10^{3}] ;Signal [ADC]");
            
            TProfile * chanDataHisto = new TProfile(tempHistoName.c_str(),"",2*_nPulses+1,
                    (-1)*_finalCharge*1e-3,_finalCharge*1e-3,_adcmin,_adcmax);
            chanDataHisto->SetTitle(tempHistoTitle.c_str());
            _rootObjectMap.insert(std::make_pair(tempHistoName, chanDataHisto));
            
            // The gaussian fit associated: 
            //  + mean: corrected pedestal per channel
            //  + sigma: total noise per channel
            const std::string tempFitName = getChanDataFitName(ichip,ichan);
            TF1 *chanDataFit = new TF1(tempFitName.c_str(),"pol1");
            _rootObjectMap.insert(std::make_pair(tempFitName, chanDataFit));
        }
    }
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms" << std::endl;
}

// Get the injected pulse (in electrons)
float AlibavaCalibrateProcessor::getInjectedCharge(const int & eventNumber)
{
    // WARNING: XXX ---> maybe use std::modf function??
    return _initialCharge+_deltaCharge*((eventNumber)/(_nSamplesPerPulse));
}

