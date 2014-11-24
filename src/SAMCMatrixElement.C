#include "SAMCMatrixElement.h"

//______________________________________________________________________________
SAMCMatrixElement::SAMCMatrixElement(){
   Init();
}
//____________________________________________________________
SAMCMatrixElement::~SAMCMatrixElement(){
   Clear();
}
//____________________________________________________________
void SAMCMatrixElement::Init(){
   fDefaultPolyOrder=0; 
   fPoly.resize(fDefaultPolyOrder);
   fPW.resize(0);
   fIsZero = true;  // if the WHOLE LINE is zero, we skip it
   fOrder  = 0; 
   fV      = 0; 
}
//____________________________________________________________
void SAMCMatrixElement::Clear(){
   fPoly.clear();
   fPW.clear();
   Init();
}
//____________________________________________________________
void SAMCMatrixElement::Print(){

   using namespace std;
   int N = fPW.size();
   int M = fPoly.size();

   std::cout << "-----------------------------------" << std::endl;
   std::cout << "Exponents: " << std::endl;
   for(int i=0;i<N;i++){
      std::cout << i << "\t" << fPW[i] << std::endl;
   }
   std::cout << "Polynomials: " << std::endl;
   std::cout << "Order: " << fOrder << std::endl;
   for(int i=0;i<M;i++){
      std::cout << i          << "\t" 
         << scientific << setprecision(4) 
         << fPoly[i] << std::endl;
   }

}
//______________________________________________________________________________

