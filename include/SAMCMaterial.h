#ifndef SAMCMaterial_HH
#define SAMCMaterial_HH

#include "TNamed.h"

class SAMCMaterial : public TNamed {

   public:

      SAMCMaterial(){
         // These numbers are made up
         Z   = 1;
         A   = 1;
         M   = 0.938;
         rho = 1.0;
         T   = 1.0;
         TR  = 1.0;
         bt  = 0.01;
         X0  = 1.0;
         L   = 1.0;
      }

      virtual ~SAMCMaterial(){
      }

      std::string Name; // Name
      int    Z;         // Z
      double A;         // A
      double M;         // Mass
      double T;         // Thickness g/cm^2
      double TR;        // Thickness in Rad_Len
      double rho;       // density g/cm^3
      double bt;        // bt
      double X0;        // Radiation Length g/cm^2
      double L;         // Length:distance along Z axix in TCS

      ClassDef(SAMCMaterial,1)
};

#endif

