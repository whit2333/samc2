#include "SAMCManager.h"
#include <iomanip>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include "yaml-cpp/yaml.h"

//ClassImp(SAMCManager)

SAMCManager * SAMCManager::fgSAMCManager = 0;

//_____________________________________________________________________________
SAMCManager::SAMCManager() {

   fVerbosity      = 2;
   fOutputFile     = 0;
   fOutputTree     = 0;
   fNumberOfEvents = 0;

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

   fFile_Name       = ""; //input file name
   RfunDB_FileName = ""; //Rfun DB input file name

   fNCuts        = -1;
   fXY           = NULL;
   fLineProperty = NULL;


}
//_____________________________________________________________________________
SAMCManager::~SAMCManager(){
      fgSAMCManager=0;
}
//_____________________________________________________________________________
Int_t SAMCManager::LoadConfig(const char * filename){

   // yaml-cpp library format
   YAML::Node aNode = YAML::LoadFile("filename");

   fNumberOfEvents    = aNode["Nevents"].as<int>();
   IsMultiScat        = aNode["debug"]["MS"].as<int>();
   IsEnergyLoss       = aNode["debug"]["EL"].as<int>();
   Which_Kin          = aNode["debug"]["PS"].as<int>();

   HRS_L              = aNode["HRS_parameters"]["HRS_length"].as<double>();
   Q1_Exit_Eff_L      = aNode["HRS_parameters"]["Q1_length"].as<double>();
   D_Entrance_Eff_L   = aNode["HRS_parameters"]["D_entrance_length"].as<double>();
   D_Exit_Eff_L       = aNode["HRS_parameters"]["D_exit_legnth"].as<double>();
   Q3_Entrance_Eff_L  = aNode["HRS_parameters"]["Q3_entrance_length"].as<double>();
   Q3_Exit_Eff_L      = aNode["HRS_parameters"]["Q3_exit_length"].as<double>();
   FP_Eff_L           = aNode["HRS_parameters"]["FP_length"].as<double>();
   Q1_Radius          = aNode["HRS_parameters"]["Q1_radius"].as<double>();
   D_X_Radius         = aNode["HRS_parameters"]["D_x_radius"].as<double>();
   D_Y_L              = aNode["HRS_parameters"]["D_y_half_length"].as<double>();
   Q3_Entrance_Radius = aNode["HRS_parameters"]["Q3_entrance_radius"].as<double>();
   Q3_Exit_Radius     = aNode["HRS_parameters"]["Q3_exit_radius"].as<double>();
   VDC_Res_x          = aNode["HRS_parameters"]["WC_resolution_x"].as<double>();
   VDC_Res_y          = aNode["HRS_parameters"]["WC_resolution_y"].as<double>();
   VDC_Res_th         = aNode["HRS_parameters"]["WC_resolution_th"].as<double>();
   VDC_Res_ph         = aNode["HRS_parameters"]["WC_resolution_ph"].as<double>();
   D_x                = aNode["HRS_parameters"]["D_x"].as<double>();
   D_y                = aNode["HRS_parameters"]["D_y"].as<double>();

   delta_dp           = aNode["monte_carlo_parameters"]["delta_dp"].as<double>();
   delta_th           = aNode["monte_carlo_parameters"]["delta_th"].as<double>();
   delta_ph           = aNode["monte_carlo_parameters"]["delta_ph"].as<double>();
   Which_Beam_Profile = aNode["monte_carlo_parameters"]["beam_profile"].as<int>();
   raster_x_size      = aNode["monte_carlo_parameters"]["raster_x"].as<double>();
   raster_y_size      = aNode["monte_carlo_parameters"]["raster_y"].as<double>();
   gaus_x_sigma       = aNode["monte_carlo_parameters"]["beam_profile_x"].as<double>();
   gaus_y_sigma       = aNode["monte_carlo_parameters"]["beam_profile_y"].as<double>();
   beam_x_center      = aNode["monte_carlo_parameters"]["beam_center_x"].as<double>();
   beam_y_center      = aNode["monte_carlo_parameters"]["beam_center_y"].as<double>();
   z0                 = aNode["monte_carlo_parameters"]["target_center_z0"].as<double>();
   T_L                = aNode["monte_carlo_parameters"]["target_length"].as<double>();
   T_H                = aNode["monte_carlo_parameters"]["target_height"].as<double>();
   RfunDB_FileName    = aNode["monte_caralo_parameters"]["run_db_filename"].as<std::string>();

   E0                 = aNode["HRS_settings"]["Es"].as<double>();
   P0                 = aNode["HRS_settings"]["P0"].as<double>();

   return 0;
}
//_____________________________________________________________________________


