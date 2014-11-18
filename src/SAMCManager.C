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

   fVerbosity             = 2;
   fOutputFile            = 0;
   fOutputFileName        = "";
   fOutputTree            = 0;
   fOutputTreeName        = "";
   fNumberOfEvents        = 0;
   fUserOutputGenFileName = "";

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
   fOutputFileName    = aNode["output_file"].as<std::string>();

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

   fTargetMaterial.Name = aNode["target"]["target_name"].as<string>();
   fTargetMaterial.Z    = aNode["target"]["target_Z"].as<double>();
   fTargetMaterial.A    = aNode["target"]["target_A"].as<double>();
   fTargetMaterial.T    = aNode["target"]["target_thickness"].as<double>();
   fTargetMaterial.rho  = aNode["target"]["target_density"].as<double>();
   fTargetMaterial.L    = T_L;//in HCS

   fMat0.Name = aNode["target"]["target_window_initial_name"].as<string>();
   fMat0.Z    = aNode["target"]["target_window_initial_Z"].as<double>();
   fMat0.A    = aNode["target"]["target_window_initial_A"].as<double>();
   fMat0.T    = aNode["target"]["target_window_initial_thickness"].as<double>();
   fMat0.rho  = aNode["target"]["target_window_initial_density"].as<double>();
   //fMat0.L    = T_L;//in HCS
   fMat1.Name = aNode["target"]["target_window_final_name"].as<string>();
   fMat1.Z    = aNode["target"]["target_window_final_Z"].as<double>();
   fMat1.A    = aNode["target"]["target_window_final_A"].as<double>();
   fMat1.T    = aNode["target"]["target_window_final_thickness"].as<double>();
   fMat1.rho  = aNode["target"]["target_window_final_density"].as<double>();
   //fMat1.L    = T_L;//in HCS

   fTheta_Target  = aNode["target"]["target_rotation"].as<double>();

   return 0;
}
//_____________________________________________________________________________
void SAMCManager::PrintConfig(){
   // this gross use of printf needs fixed -whit
   printf("%-*s=%*s %-*s %-*s\n",  15,"samc_rootfilename",10,fOutputFileName.c_str(),8,"",  40,"(Save Results to Root File Name)");
   if ( !fUserOutputGenFileName.empty() ) {
      printf("%-*s=%*s %-*s %-*s\n",  15,"fUserOutputGenFileName", 10, fUserOutputGenFileName.c_str(),8,"",  40,"(User-Defined Generator File Name)");
   }
   printf("%-*s=%*d %-*s %-*s\n"   , 15 , "Num_Of_Events"     , 10 , fNumberOfEvents     , 8 , ""    , 40 , "(Number Of Events)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "HRS_L"             , 10 , HRS_L             , 8 , "cm"  , 40 , "(Length of HRS)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "Q1_Exit_Eff_L"     , 10 , Q1_Exit_Eff_L     , 8 , "cm"  , 40 , "(Length of Magnetic)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "D_Entrance_Eff_L"  , 10 , D_Entrance_Eff_L  , 8 , "cm"  , 40 , "(Length of Magnetic)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "D_Exit_Eff_L"      , 10 , D_Exit_Eff_L      , 8 , "cm"  , 40 , "(Length of Magnetic)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "Q3_Entrance_Eff_L" , 10 , Q3_Entrance_Eff_L , 8 , "cm"  , 40 , "(Length of Magnetic)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "Q3_Exit_Eff_L"     , 10 , Q3_Exit_Eff_L     , 8 , "cm"  , 40 , "(Length of Magnetic)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "FP_Eff_L"          , 10 , FP_Eff_L          , 8 , "cm"  , 40 , "(Length of Magnetic)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "VDC_Res_x"         , 10 , VDC_Res_x         , 8 , "cm"  , 40 , "(Resolution of VDC x)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "VDC_Res_y"         , 10 , VDC_Res_y         , 8 , "cm"  , 40 , "(Resolution of VDC y)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "VDC_Res_th"        , 10 , VDC_Res_th        , 8 , "mr"  , 40 , "(Resolution of VDC th)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "VDC_Res_ph"        , 10 , VDC_Res_ph        , 8 , "mr"  , 40 , "(Resolution of VDC ph)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "delta_dp"          , 10 , delta_dp*100      , 8 , "%"   , 40 , "(full width of \xce\x94(dp) for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "delta_th"          , 10 , delta_th          , 8 , ""    , 40 , "(full width of \xce\x94(tan(\xce\xb8tg)) for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "delta_ph"          , 10 , delta_ph          , 8 , ""    , 40 , "(full width of \xce\x94(tan(\xcf\x86tg)) for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "raster_x_size"     , 10 , raster_x_size     , 8 , "cm"  , 40 , "(raster x full size for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "raster_y_size"     , 10 , raster_y_size     , 8 , "cm"  , 40 , "(raster y full size for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "beam_x_center"     , 10 , beam_x_center     , 8 , "cm"  , 40 , "(beam x center for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "beam_y_center"     , 10 , beam_y_center     , 8 , "cm"  , 40 , "(beam y center for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "z0"                , 10 , z0                , 8 , "cm"  , 40 , "(target center for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "T_L"               , 10 , T_L               , 8 , "cm"  , 40 , "(target length for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "T_H"               , 10 , T_H               , 8 , "cm"  , 40 , "(target height for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "D_x"               , 10 , D_x               , 8 , "cm"  , 40 , "(X offset(in TCS) between TCS and HCS)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "D_y"               , 10 , D_y               , 8 , "cm"  , 40 , "(Y offset(in TCS) between TCS and HCS)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "E0"                , 10 , E0                , 8 , "MeV" , 40 , "(Incident Beam Energy for generator)");
   printf("%-*s=%*.2f %-*s %-*s\n" , 15 , "P0"                , 10 , P0                , 8 , "MeV" , 40 , "(HRS Setting Momentum for generator)");
}
//______________________________________________________________________________


