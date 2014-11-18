#ifndef SAMCMaterial_HH
#define SAMCMaterial_HH

#include "TNamed.h"

class SAMCMaterial : public TNamed {

   public:

      SAMCMaterial();
      virtual ~SAMCMaterial();

      SAMCMaterial(const SAMCMaterial& v) ;
      const SAMCMaterial& operator=(const SAMCMaterial& v) ;

      int    fZ;         // Z
      double fA;         // A
      double fM;         // Mass
      double fT;         // Thickness g/cm^2
      double fTR;        // Thickness in Rad_Len
      double fDensity;   // density g/cm^3
      double fbt;        // bt
      double fX0;        // Radiation Length g/cm^2
      double fL;         // Length:distance along Z axix in TCS

      ClassDef(SAMCMaterial,1)
};

#endif

