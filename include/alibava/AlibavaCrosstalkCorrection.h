/*
 * Created by Jordi Duarte-Campderros
 *  (2017 IFCA-CERN) jorge.duarte.campderros@cern.ch
 *
 * Processor to correct the effect of the assymetric 
 * cross-talking
 *
 */

#ifndef ALIBAVACROSSTALKCORRECTION_H
#define ALIBAVACROSSTALKCORRECTION_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"


// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>


namespace alibava 
{
    class AlibavaCrosstalkCorrection:public alibava::AlibavaBaseProcessor   
    {
        public:
            //! Returns a new instance of AlibavaCrosstalkCorrection
	    /*! This method returns an new instance of the this processor.  It
	     *  is called by Marlin execution framework and it shouldn't be
	     *  called/used by the final user.
	     *
	     *  @return a new AlibavaCrosstalkCorrection.
	     */
	    virtual Processor * newProcessor () 
            {
                return new AlibavaCrosstalkCorrection;
            }
	    
	    //! Default constructor
	    AlibavaCrosstalkCorrection ();
	    
	    //! Called at the job beginning.
	    /*! This is executed only once in the whole execution. It prints
	     *  out the processor parameters and reset all needed data
	     *  members. In principle this could also be the right place to
	     *  book histograms, but since those are also grouped according to
	     *  the detector numbers we need to have this parameter available.
	     */
	    virtual void init ();
	    
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
	    virtual void processRunHeader (LCRunHeader * run);
	    
	    //! Called every event
	    /*! Since the behavior of the PedestalNoise processor is different
	     *  if this is the first or one of the following loop, this method
	     *  is just calling
	     *  AlibavaCrosstalkCorrection::firstLoop(LCEvent*) or
	     *  AlibavaCrosstalkCorrection::otherLoop(LCEvent*)
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
	    virtual void check (LCEvent * evt);
	    
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
	    
	    //////////////////////////
	    // Processor Parameters //
	    //////////////////////////
	    
        private:
	    // The first neighbour coefficient
	    float _b1;
    };
    //! A global instance of the processor
    AlibavaCrosstalkCorrection gAlibavaCrosstalkCorrection;
}

#endif
