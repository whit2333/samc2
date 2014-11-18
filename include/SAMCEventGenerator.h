#ifndef SAMCEventGenerator_HH
#define SAMCEventGenerator_HH

#include "TNamed.h"

/** Event generator.
 *  Generates the inputs for an event.
 */
class SAMCEventGenerator : public TNamed {

   public:
      SAMCEventGenerator(){}
      virtual ~SAMCEventGenerator(){}


   ClassDef(SAMCEventGenerator,1)
};

#endif

