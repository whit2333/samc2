#ifndef SAMCEventGenerator_HH
#define SAMCEventGenerator_HH

#include "TNamed.h"
#include "TRandom.h"
#include "SAMCManager.h"
#include "SAMCEvent.h"

/** Event generator.
 *  Generates the inputs for an event.
 */
class SAMCEventGenerator : public TNamed {

   protected:

      Int_t     fBeamProfileType;
      Double_t  fRasterXSize ;
      Double_t  fRasterYSize ;
      Double_t  fGausXSigma  ;
      Double_t  fGausYSigma  ;
      Double_t  fBeamXCenter ;
      Double_t  fBeamYCenter ;
      Double_t  fDelta_th;
      Double_t  fDelta_ph;
      Double_t  fDelta_dp;
      Double_t  fz0      ;
      Double_t  fT_L     ;

   public:

      SAMCEventGenerator();
      virtual ~SAMCEventGenerator();

      void Print();
      void Init();

      virtual Int_t GenerateEvent(SAMCEvent& event);

      void  SetBeamProfileType(Int_t t){fBeamProfileType = t;}
      Int_t GetBeamProfileType(){return fBeamProfileType;}

   ClassDef(SAMCEventGenerator,1)
};

#endif

