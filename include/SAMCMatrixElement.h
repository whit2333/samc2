#ifndef SAMCMatrixElement_H
#define SAMCMatrixElement_H  1

#include <iostream>
#include <cstdlib> 
#include <iomanip> 
#include <vector> 
#include "TObject.h"


class SAMCMatrixElement {

   private: 
      bool                 fIsZero;           // whether or not the element is zero
      int                  fOrder;            // order of the polynomial
      int                  fDefaultPolyOrder; // default order of the polynomial
      double               fV;                // its computed value
      std::vector<int>     fPW;               // exponents of matrix element (e.g., D100 = { 1, 0, 0 })
      std::vector<double>  fPoly;             // the associated polynomial

   public:
      SAMCMatrixElement();
      ~SAMCMatrixElement();

      void Init();
      void Print(); 
      void Clear();

      void SetValue(double v){fV = v;} 
      void SetBoolean(bool v){fIsZero = v;} 
      void SetPolyOrder(int order){fOrder = order;} 
      void SetExpoSize(int i){fPW.resize(i);} 
      void SetPolySize(int i){fPoly.resize(i);} 
      void SetExpo(int i,int v){fPW[i] = v;} 
      void SetPoly(int i,double v){fPoly[i] = v;}

      int GetPolyOrder(){return fOrder;} 
      int GetExpo(int i){return fPW[i];} 
      int GetExpoSize(){return fPW.size();} 
      double GetValue(){return fV;} 
      double GetPoly(int i){return fPoly[i];}
      bool IsElementZero(){return fIsZero;} 
      //bool Match( const SAMCMatrixElement& rhs ) const;


   ClassDef(SAMCMatrixElement,1)
};

#endif 

