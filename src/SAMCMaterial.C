#include "SAMCMaterial.h"

SAMCMaterial::SAMCMaterial(){
   // These numbers are made up
   fZ       = 1;
   fA       = 1;
   fM       = 0.938;
   fDensity = 1.0;
   fT       = 1.0;
   fTR      = 1.0;
   fbt      = 0.01;
   fX0      = 1.0;
   fL       = 1.0;
}
//______________________________________________________________________________
SAMCMaterial::~SAMCMaterial(){
}
//______________________________________________________________________________
SAMCMaterial::SAMCMaterial(const SAMCMaterial& v) : TNamed(v) {
   fZ       = v.fZ ;
   fA       = v.fA ;
   fM       = v.fM ;
   fDensity = v.fDensity;
   fT       = v.fT ;
   fTR      = v.fTR;
   fbt      = v.fbt;
   fX0      = v.fX0;
   fL       = v.fL ;

}
//______________________________________________________________________________
const SAMCMaterial& SAMCMaterial::operator=(const SAMCMaterial& v) {
   if (this != &v) {
      TNamed::operator=(v);
      fZ       = v.fZ ;
      fA       = v.fA ;
      fM       = v.fM ;
      fDensity = v.fDensity;
      fT       = v.fT ;
      fTR      = v.fTR;
      fbt      = v.fbt;
      fX0      = v.fX0;
      fL       = v.fL ;
   }
   return *this;
}
//______________________________________________________________________________
void SAMCMaterial::Print(Option_t * ) {
   std::cout << " - Material : " << GetName() << std::endl;
   std::cout << "   fZ       : " << fZ       << std::endl;
   std::cout << "   fA       : " << fA       << std::endl;
   std::cout << "   fM       : " << fM       << std::endl;
   std::cout << "   fDensity : " << fDensity << std::endl;
   std::cout << "   fT       : " << fT       << std::endl;
   std::cout << "   fTR      : " << fTR      << std::endl;
   std::cout << "   fbt      : " << fbt      << std::endl;
   std::cout << "   fX0      : " << fX0      << std::endl;
   std::cout << "   fL       : " << fL       << std::endl;

}
//______________________________________________________________________________

