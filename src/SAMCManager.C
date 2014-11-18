#include "SAMCManager.h"
#include <iomanip>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream

//ClassImp(SAMCManager)
//_____________________________________________________________
// A description of the class starts with the line above, and
// will take place here !
//

SAMCManager * SAMCManager::fgSAMCManager = 0;
//_____________________________________________________________________________

SAMCManager::SAMCManager() {
   fVerbosity = 2;
   IsDebug            = false; // output to debugfile
   IsMultiScat        = true;  // enable multi-scattering
   IsEnergyLoss       = false; // enable Energy Loss
   Which_Beam_Profile = 1;     // 0=Norm, 1=Gaus, 2=gaus+norm= triangle
   Which_Kin          = 0;     // 0= phase space, 1 = elastic, 2= quasi-elastic
   HRS_L              = 116;   // Length Of HRS ( Left or Right,here _L         = Length) (cm)
   Q1_Exit_Eff_L      = 40;    // Q1 Effective Length  (cm)
   D_Entrance_Eff_L   = 400;   // D Entrance Effective Length  (cm)
   D_Exit_Eff_L       = 700;   // D Exit Effective Length  (cm)
   Q3_Entrance_Eff_L  = 775;   // Q3 Entrance Effective Length  (cm)
   Q3_Exit_Eff_L      = 865;   // Q3 Exit Effective Length  (cm)
   FP_Eff_L           = 1130;  // FP Effective Length  (cm)
   Q1_Radius          = 14.92; // Q1 Radius  (cm)
   D_X_Radius         = 40;    // Dipole X dir Radius  (cm)
   D_Y_L              = 12.5;  // Dipole Y Length when x                       = 0  (cm)
   Q3_Entrance_Radius = 30;    // Q3 Entrance Radius  (cm)
   Q3_Exit_Radius     = 30;    // Q3 Exit Radius  (cm)
   VDC_Res_x          = 0.013; // resolution of Wire chamer x  (cm)
   VDC_Res_y          = 0.013; // resolution of Wire chamer y  (cm)
   VDC_Res_th         = 0.3;   // resolution of Wire chamer th (mr)
   VDC_Res_ph         = 0.3;   // resolution of Wire chamer ph (mr)
   delta_dp           = 0.14;  // d(dp) dp full width for generator
   delta_th           = 0.18;  // d(th) th full width for generator(tan(th))
   delta_ph           = 0.09;  // d(ph) th full width for generator(tan(ph))
   gaus_x_sigma       = 0;     // if Which_Beam_Profile == 1/2 sigma of x beam for generator(cm)
   gaus_y_sigma       = 0;     // if Which_Beam_Profile == 1/2 sigma of y beam for generator(cm)
   raster_x_size      = 0;     // raster x full size for generator(cm)
   raster_y_size      = 0;     // raster y full size for generator(cm)
   beam_x_center      = 0;     // beam x center for generator(cm)
   beam_y_center      = 0;     // beam y center for generator(cm)
   z0                 = 0;     // target center for generator(cm)
   T_L                = 0;     // target length for generator(cm)
   T_H                = 0;     // target height for generator(cm)

   D_x                = 0;     // X offset(in TCS) between TCS and HCS
   D_y                = 0;     // Y offset(in TCS) between TCS and HCS
   E0                 = 0.0;   // =incident beam energy for generator (MeV)
   P0                 = 0.0;   // =HRS Setting Momentum for generator(MeV)

   File_Name       = ""; //input file name
   RfunDB_FileName = ""; //Rfun DB input file name

   fNCuts        = -1;
   fXY           = NULL;
   fLineProperty = NULL;

   fNumberOfEvents = 0;

}
//_____________________________________________________________________________

SAMCManager::~SAMCManager(){
      fgSAMCManager=0;
}
//_____________________________________________________________________________

