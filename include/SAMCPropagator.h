#ifndef SAMCPropagator_HH
#define SAMCPropagator_HH 1

#include "TNamed.h"

/** Propator class which sets the various event data members
 *  as it propagates forward and backwards through the detector.
 *
 */
class SAMCPropagator : public TNamed {

   protected:

   public:
      SAMCPropagator(){}
      virtual ~SAMCPropagator(){}


   ClassDef(SAMCPropagator,0)
};

#endif

