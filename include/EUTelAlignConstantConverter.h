// Version: $Id$
/*
 *   Copy the alignment constants into a gear file
 *
 *   author: Jordi Duarte-Campderros (Aug. 2017)
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELALIGNCONSTANTCONVERTER_H
#define EUTELALIGNCONSTANTCONVERTER_H

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
//#include <EVENT/LCRunHeader.h>
//#include <EVENT/LCEvent.h>

// system includes <>
#include<string>
#include<vector>
#include<map>

//forward declarations

namespace eutelescope 
{
    // forward declaration
    class EUTelAlignmentConstant;

    class EUTelAlignConstantConverter : public marlin::Processor 
    {
        public:
            virtual Processor * newProcessor() { return new EUTelAlignConstantConverter; }
    
            EUTelAlignConstantConverter();
            virtual void init() ;
            virtual void processRunHeader( LCRunHeader* run );
            virtual void processEvent( LCEvent * evt );
            virtual void check( LCEvent * /*evt*/ ){;};
            virtual void end();
        
        private:
            int _iRun;
            int _iEvt;

            std::string _GEARFileSuffix;
            std::vector<std::string> _collectionNames;
            std::map<int,std::vector<EUTelAlignmentConstant*> > _alignmentMap;
    };
    
    //! A global instance of the processor
    EUTelAlignConstantConverter gEUTelAlignConstantConverter;
}
#endif

