/* 
 * Calibrate alibava ADCs counts
 *
 *  J. Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *    - Clean and coding style (Allman)
 */

#ifndef ALIBAVACALIBRATEPROCESSOR_H
#define ALIBAVACALIBRATEPROCESSOR_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>


namespace alibava 
{
    //! Pedestal and noise  processor for Marlin.
    class AlibavaCalibrateProcessor: public alibava::AlibavaBaseProcessor   
    {
        public:
            //! Returns a new instance of AlibavaCalibrateProcessor
	    /*! This method returns an new instance of the this processor.  It
	     *  is called by Marlin execution framework and it shouldn't be
	     *  called/used by the final user.
	     *
	     *  @return a new AlibavaCalibrateProcessor.
	     */
            virtual Processor * newProcessor() 
            {
                return new AlibavaCalibrateProcessor;
            }
            //! Default constructor
            AlibavaCalibrateProcessor();
            
            //! Called at the job beginning.
	    /*! This is executed only once in the whole execution. It prints
	     *  out the processor parameters and reset all needed data
	     *  members. In principle this could also be the right place to
	     *  book histograms, but since those are also grouped according to
	     *  the detector numbers we need to have this parameter available.
	     */
	    virtual void init();
		
	    //! Called for every run.
	    /*! At the beginning of every run the run header is read and
	     *  processed by this method. As a first thing, the input run
	     *  header is dynamically re-casted to a EUTelRunHeaderImpl and
	     *  then important things like the number of detectors and the
	     *  pixel detector boundaries are dumped from the file. After that
	     *  the EUTelPedestalNoiseProcess::bookHistos() is called.
	     *
	     *  @param run the LCRunHeader of the this current run
	     */
	    virtual void processRunHeader(LCRunHeader * run);
		
	    //! Called every event
	    /*! Since the behavior of the PedestalNoise processor is different
	     *  if this is the first or one of the following loop, this method
	     *  is just calling
	     *  AlibavaCalibrateProcessor::firstLoop(LCEvent*) or
	     *  AlibavaCalibrateProcessor::otherLoop(LCEvent*)
	     *
	     *  @param evt the current LCEvent event as passed by the
	     *  ProcessMgr
	     */
	    virtual void processEvent(LCEvent * evt);
	    
	    //! Check event method
	    /*! This method is called by the Marlin execution framework as
	     *  soon as the processEvent is over. It can be used to fill check
	     *  plots. For the time being there is nothing to check and do in
	     *  this slot.
	     *
	     *  @param evt The LCEvent event as passed by the ProcessMgr
	     *
	     */
	    virtual void check(LCEvent * evt);
	    
	    //! Book histograms
	    /*! This method is used to prepare the needed directory structure
	     *  within the current ITree folder and books all required
	     *  histograms. Histogram pointers are stored into
	     *  EUTelPedestalNoiseProcess::_rootObjectMap so that they can be
	     *  recalled and filled from anywhere in the code.  Apart from the
	     *  histograms listed in AlibavaCalibrateProcessor::fillHistos()
	     *  there is also a common mode histo described here below:
	     *
	     *  \li commonModeHisto: 1D histogram to store the calculated
	     *  common mode value for each event. This histogram is booked and
	     *  filled only if the loop counter is greater-equal than 1,
	     *  because for _iLoop == 0 there is no common mode suppression.
	     *  This histo is not filled with the other because it needs to be
	     *  updated every event.
	     *
	     *  @see AlibavaCalibrateProcessor::fillHistos() for the todos
	     */
	    void bookHistos();
	    
	    //! Fill histograms
	    /*! This method is used to fill in histograms for each channel. 
	     */
	    void fillHistos(TrackerDataImpl * trkdata, const int & evt);
	    

	    //! Called after data processing.
	    /*! This method is called when the loop on events is finished. It
	     *  is checking whether the calculation is properly finished or
	     *  not.
	     *  A very common error occurs when the file finished without a
	     *  EORE or when the MaxRecordNumber was set to low to loop over
	     *  all the needed record. To check this is very easy because we
	     *  just have to crosscheck if _iLoop is equal to noOfCMIterations.
	     */
	    virtual void end();

            /* Calculate the injected charge (in electrons) in the given pulse
             * defined by the number of events
             */
            float getInjectedCharge(const int & evtNumber);
	
        protected:
            //! Name of the Pedestal histogram 
	    std::string _calibrateHistoName;
	    
            //! Name of the Temperature histogram
	    std::string _temperatureHistoName;
	    
	    //! The name of the histogram used to calculate calibration curves
	    /*! For every channel a histogram will be created
	     *  and filled with the readings of that channel
	     *  
	     *  The actual name of the histogram will be 
	     *  _chanRawDataHistoName+chanNum
	     */
	    std::string _chanDataHistoName;
	    
	    //! The name of the fits used to calculate pedestal and noise
	    /*! For every channel a histogram will be created
	     *  and filled with the readings of that channel
	     *
	     *  The actual name of the histogram will be
	     *  _chanRawDataFitName+"_chip"+chipnum+"_chan"+chanNum
	     */
	    std::string _chanDataFitName;

            //! Histogram tune: x-range and Nbins
            float _adcmax;
            float _adcmin;
            int _nbins;

            // Provisional
            int _nPulses;
            int _initialCharge;
            int _finalCharge;
            int _deltaCharge;
            int _nSamplesPerPulse;

	    //! The function that returns name of the histogram for each channel
	    std::string getChanDataHistoName(unsigned int ichip, unsigned int ichan);
	    
	    //! The function that returns name of the fit for each channel
	    std::string getChanDataFitName(unsigned int ichip, unsigned int ichan);

	    //! The function that returns name of the calibration histogram for each chip
	    std::string getCalChargeHistoName(unsigned int ichip);
    };
    
    //! A global instance of the processor
    AlibavaCalibrateProcessor gAlibavaCalibrateProcessor;	
}

#endif
