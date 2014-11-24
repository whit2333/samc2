#ifndef SAMCAnalyzerOptics_H
#define SAMCAnalyzerOptics_H 1

#include <cstdlib> 
#include <iostream>
#include <string>
#include <iomanip>
#include <vector> 
#include "SAMCMatrixElement.h"
#include "TString.h"
#include <fstream>

class SAMCAnalyzerOptics{

   private:
      int fP0;                            // Momentum setting.  Needed for CSR to chose the right DB file.  
      int fOpticsPackage;                 // switch to choose the analyzer Hall A optics package 
      int fScatAngle;                     // scattering angle (used to switch from LHRS to RHRS, based on if th > 0 or not) 
      std::vector<SAMCMatrixElement> fFP;          // transport to focal plane: fFP stores t,y,p in that order (size 3)  
      std::vector<SAMCMatrixElement> fY,fD,fP,fT;  // focal plane to target: Y = y, D = delta p, P = phi, T = theta 
      std::vector<SAMCMatrixElement> fYTA,fPTA;    // focal plane to target: YTA and PTA utilize abs(th_fp)
      // FIXME: Where are these variables?! 
      std::vector<SAMCMatrixElement> fL;

      enum EFPMatrixElemTags {T000,Y000,P000}; 	

   public: 
      SAMCAnalyzerOptics();
      ~SAMCAnalyzerOptics();

      void Init();
      void ClearVectors();
      void LoadData(int); 
      void PrintMatrixElement(std::string); 
      void SetOpticsPackage(int);
      void SetMomentumSetting(int p){fP0 = p;}  
      void SetScatteringAngle(int th){fScatAngle = th;} 
      void CalculateFocalPlaneCoords(std::vector<double>,std::vector<double> &);
      void CalculateTargetCoords(std::vector<double>,std::vector<double> &);
      void CalculateMatrix(const double,std::vector<SAMCMatrixElement> &); 

      int GetOpticsPackage(){return fOpticsPackage;} 
      double CalculateTargetVar(std::vector<SAMCMatrixElement> &,const double [][5]); 

};

#endif
