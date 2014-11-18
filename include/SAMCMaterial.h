#ifndef SAMCMaterial_HH
#define SAMCMaterial_HH

#include "TNamed.h"

class SAMCMaterial : public TNamed {

   public:

      SAMCMaterial();
      virtual ~SAMCMaterial();

      SAMCMaterial(const SAMCMaterial& v) ;
      const SAMCMaterial& operator=(const SAMCMaterial& v) ;

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

