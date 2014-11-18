#ifndef SAMCManager_HH
#define SAMCManager_HH 1

#include "SAMCConfig.h"
#include "SAMCMaterial.h"

#include "TNamed.h"
#include "TObject.h"
#include "TROOT.h"
#include <iostream>
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectory.h"
#include <set>
#include <string>
#include "TBrowser.h"
#include "TObjString.h"
#include "TMultiGraph.h"

/** Manager singleton class.
 *  
 */
class SAMCManager : public TObject {

   private:

      Int_t   fVerbosity;

   public:
      // Todo: rename all the variables below

      //default unit
      //MeV cm g rad
      //cross section output unit = nbarn/sr
      bool   IsDebug            ;//= false; // output to debugfile
      bool   IsMultiScat        ;//= true;  // enable multi-scattering
      bool   IsEnergyLoss       ;//= false; // enable Energy Loss
      int    Which_Beam_Profile ;//= 1;     // 0=Norm, 1=Gaus, 2=gaus+norm= triangle
      int    Which_Kin          ;//= 0;     // 0= phase space, 1 = elastic, 2= quasi-elastic
      double HRS_L              ;//= 116;   // Length Of HRS ( Left or Right,here _L         = Length) (cm)
      double Q1_Exit_Eff_L      ;//= 40;    // Q1 Effective Length  (cm)
      double D_Entrance_Eff_L   ;//= 400;   // D Entrance Effective Length  (cm)
      double D_Exit_Eff_L       ;//= 700;   // D Exit Effective Length  (cm)
      double Q3_Entrance_Eff_L  ;//= 775;   // Q3 Entrance Effective Length  (cm)
      double Q3_Exit_Eff_L      ;//= 865;   // Q3 Exit Effective Length  (cm)
      double FP_Eff_L           ;//= 1130;  // FP Effective Length  (cm)
      double Q1_Radius          ;//= 14.92; // Q1 Radius  (cm)
      double D_X_Radius         ;//= 40;    // Dipole X dir Radius  (cm)
      double D_Y_L              ;//= 12.5;  // Dipole Y Length when x                       = 0  (cm)
      double Q3_Entrance_Radius ;//= 30;    // Q3 Entrance Radius  (cm)
      double Q3_Exit_Radius     ;//= 30;    // Q3 Exit Radius  (cm)
      double VDC_Res_x          ;//= 0.013; // resolution of Wire chamer x  (cm)
      double VDC_Res_y          ;//= 0.013; // resolution of Wire chamer y  (cm)
      double VDC_Res_th         ;//= 0.3;   // resolution of Wire chamer th (mr)
      double VDC_Res_ph         ;//= 0.3;   // resolution of Wire chamer ph (mr)
      double delta_dp           ;//= 0.14;  // d(dp) dp full width for generator
      double delta_th           ;//= 0.18;  // d(th) th full width for generator(tan(th))
      double delta_ph           ;//= 0.09;  // d(ph) th full width for generator(tan(ph))
      double gaus_x_sigma       ;//= 0;     // if Which_Beam_Profile == 1/2 sigma of x beam for generator(cm)
      double gaus_y_sigma       ;//= 0;     // if Which_Beam_Profile == 1/2 sigma of y beam for generator(cm)
      double raster_x_size      ;//= 0;     // raster x full size for generator(cm)
      double raster_y_size      ;//= 0;     // raster y full size for generator(cm)
      double beam_x_center      ;//= 0;     // beam x center for generator(cm)
      double beam_y_center      ;//= 0;     // beam y center for generator(cm)
      double z0                 ;//= 0;     // target center for generator(cm)
      double T_L                ;//= 0;     // target length for generator(cm)
      double T_H                ;//= 0;     // target height for generator(cm)
      double D_x                ;//= 0;     // X offset(in TCS) between TCS and HCS
      double D_y                ;//= 0;     // Y offset(in TCS) between TCS and HCS

      double E0                 ;//= 0.0;   // =incident beam energy for generator (MeV)
      double P0                 ;//= 0.0;   // =HRS Setting Momentum for generator(MeV)
      double fTheta             ;// HRS angle


      SAMCMaterial fTargetMaterial;
      SAMCMaterial fMat0;
      SAMCMaterial fMat1;
      SAMCMaterial fMat2;
      SAMCMaterial fMat3;
      SAMCMaterial fMat4;
      SAMCMaterial fMat5;
      SAMCMaterial fMat6;

      double fTheta_Target; // target rotation (viewed from above)


      std::string   fFile_Name;       //input file name
      std::string   RfunDB_FileName; //Rfun DB input file name

      int       fNCuts        ;//= -1;
      int**     fXY           ;//= NULL;
      double**  fLineProperty ;//= NULL;

      int    fNumberOfEvents;// Number of events simulated

      TFile       * fOutputFile;
      std::string   fOutputFileName;
      TTree       * fOutputTree;
      std::string   fOutputTreeName;

      std::string   fUserOutputGenFileName;


   protected:

      SAMCManager();
      static SAMCManager * fgSAMCManager;

   public :

      static SAMCManager * GetManager(){
         if(!fgSAMCManager) fgSAMCManager = new SAMCManager();
         return(fgSAMCManager);
      }
      static SAMCManager * Instance(){
         if(!fgSAMCManager) fgSAMCManager = new SAMCManager();
         return(fgSAMCManager);
      }

      virtual ~SAMCManager();

      Int_t   GetVerbosity() const {return fVerbosity;}
      void    SetVerbosity(Int_t v) {fVerbosity = v;}

      Int_t LoadConfig(const char * filename); // *MENU*
      void  PrintConfig(); // *MENU*

      ClassDef(SAMCManager,0)
};

#endif 
