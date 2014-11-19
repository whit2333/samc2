#include "SAMC.h"
#include "SAMCEvent.h"
#include "SAMCPropagator.h"
#include "SAMCEventGenerator.h"
#include "SAMCManager.h"
#include "SAMCFortran.h"

#include "getopt.h"

//______________________________________________________________________________
void usage(const char* aCommand) {
   printf("--------------------\n");
   printf("%-*s%s -dh filename\n",10,"Usage:",aCommand);
   printf("%-*s -d debug\n",10," ");
   printf("%-*s -h help\n",10," ");
   printf("%-*s filename: absolute file location\n",10," ");
   exit(8);
}
//______________________________________________________________________________
void getargs(int argc,char** argv) {

   const char* cmd = *argv;
   int  run_set    = 0;
   int  seed_set   = 10000;
   const struct option longopts[] =
   {
      {"run",       required_argument,  0, 'r'},
      {"seed",      required_argument,  0, 's'},
      {"verbose",   no_argument,        0, 'v'},
      {"debug",     no_argument,        0, 'd'},
      {"help",      no_argument,        0, 'h'},
      {0,0,0,0}
   };

   int index = 0;
   int iarg  = 0;
   //opterr    = 1; //turn off getopt error message

   while(iarg != -1) {
      iarg = getopt_long(argc, argv, "vr:s:hd", longopts, &index);
      switch (iarg)
      {
         case 'r':
            std::cout << "run should be set to " << optarg << std::endl;
            run_set = atoi(optarg);
            break;

         case 's':
            seed_set = seed_set + atoi(optarg);
            std::cout << "seed :  " << seed_set << std::endl;
            break;

         case 'h':
            //print_usage();
            usage(cmd);
            break;

         case 'v':
            SAMCManager::Instance()->SetVerbosity(2);
            break;

         case 'd':
            SAMCManager::Instance()->IsDebug=true;
            break;

      }
   }
   std::string theRest  = "";
   for (int i = optind; i < argc; i++) {
      theRest        += argv[i];
   }
   //std::cout << theRest << std::endl;

   if ( SAMCManager::Instance()->fFile_Name.empty() ) {
      SAMCManager::Instance()->fFile_Name = theRest;
   }

   if ( SAMCManager::Instance()->fFile_Name.empty() )
   {
      printf("[Error %s: Line %d] No File Input.\n",__FILE__,__LINE__);
      usage(cmd);
   }
}
//______________________________________________________________________________
//bool IsANumber(std::string aString) {
//   std::istringstream iss(aString);
//   std::ostringstream oss;
//   int x;
//   iss >> x;
//   oss << x;
//   if ( aString == oss.str() )
//      return true;
//   else
//      return false;
//}
//______________________________________________________________________________
int ReadDatabase( const std::string& aFileName ) {

   const int LEN = 200;

   char buf[LEN];

   // Read data from database

   FILE* fi = fopen( aFileName.c_str(),"r" );
   if( !fi ) return -1;

   int i,j,k;
   j=0;
   while ( fgets(buf,LEN,fi) )
   {
      i=0;
      while ( buf[i]==' ' )
      {
         ++i;
      }
      if ( buf[i]!='#' )
      {
         ++j;
      }
      //else it's comment, skipped
   }
   fclose(fi);

   SAMCManager::Instance()->fNCuts=j;

   if ( SAMCManager::Instance()->fXY ) {
      for ( i = 0; i < NELEMENTS; ++i ) {
         delete [] SAMCManager::Instance()->fXY[i];
      }
      delete [] SAMCManager::Instance()->fXY;
   }
   if ( SAMCManager::Instance()->fLineProperty ) {
      for ( i = 0; i < 3; ++i ) {
         delete [] SAMCManager::Instance()->fLineProperty[i];
      }
      delete [] SAMCManager::Instance()->fLineProperty;
   }

   SAMCManager::Instance()->fXY=new int*[SAMCManager::Instance()->fNCuts];
   SAMCManager::Instance()->fLineProperty=new double*[SAMCManager::Instance()->fNCuts];
   for ( i = 0; i < SAMCManager::Instance()->fNCuts; ++i ) {
      SAMCManager::Instance()->fXY[i]=new int[NELEMENTS];
      SAMCManager::Instance()->fLineProperty[i]=new double[3];
   }
   fi = fopen(aFileName.c_str(),"r");
   j=0;
   while ( fgets(buf,LEN,fi) )
   {
      i=0;
      while ( buf[i]==' ' )
      {
         ++i;
      }
      if ( buf[i]!='#' )
      {
         k=0;
         while ( k<NELEMENTS ) {
            if ( buf[i]=='0' || buf[i]=='1' ) {
               //printf("#%c\n",buf[i]);
               SAMCManager::Instance()->fXY[j][k++]=atoi(&buf[i]);
               buf[i]=' ';
            }
            ++i;
         }
         //printf("####%s\n",buf);
         sscanf ( buf, "%lf %lf %lf", &(SAMCManager::Instance()->fLineProperty)[j][LINE_SLOPE],&(SAMCManager::Instance()->fLineProperty)[j][LINE_INTERSECTION],&(SAMCManager::Instance()->fLineProperty)[j][LINE_SIGN] );  
         ++j;
      }
      //else it's comment, skipped
   }
   fclose(fi);
   //printf("fNCuts=%d\n",fNCuts);
   //for ( i = 0; i < fNCuts; ++i ) {
   //	for ( j = 0; j < NELEMENTS; ++j ) {
   //		printf("%d ",fXY[i][j]);
   //	}
   //	for ( j = 0; j < 3; ++j ) {
   //		printf("%g ",fLineProperty[i][j]);
   //	}
   //	printf("\n");
   //}
   //printf("---------------------\n");
   return 0;
}
//______________________________________________________________________________
int main(int argc, char** argv) {

   srand(time(NULL));

   getargs(argc,argv);

   //int i,j,k,
   //k = 0;
   int fail_events = 0;

   SAMCManager * man = SAMCManager::Instance();
   man->LoadConfig(man->fFile_Name.c_str());

   std::string samc_rootfilename   = man->fOutputFileName;
   std::string userdefgen_filename = ""; // Not sure what this is really used for -whit
   int Num_Of_Events = man->fNumberOfEvents;


   if ( man->IsDebug ) {
      if ( man->IsMultiScat ){
         printf("Multi-Scattering Enabled.\n");
      } else{
         printf("Multi-Scattering Disabled.\n");
      }

      if ( man->IsEnergyLoss ){
         printf("Energy Loss Enabled.\n");
      } else {
         printf("Energy Loss Disabled.\n");
      }

      switch ( man->Which_Kin ) {
         case 1:
            printf("Elastic.\n");
            break;
         case 2:
            printf("Quasi-Elastic.\n");
            break;
         case 0:
         default:
            printf("Phase Space.\n");
            break;
      }

      man->PrintConfig();
   }

   //man->RfunDB_FileName = "";
   std::cout << man->RfunDB_FileName << std::endl;

   if ( man->RfunDB_FileName.empty() || man->RfunDB_FileName=="NONE" ) {
      man->fNCuts=0;
   } else {
      if( ReadDatabase(man->RfunDB_FileName) ) {
         printf("[Error %s: Line %d] %s Rfun DB File Not Found.\n",__FILE__,__LINE__,man->RfunDB_FileName.c_str());
         return 11;
      }
   }

   int err = 0;
   ifstream userdefgen_file;
   std::string comment;
   double usrgen[6];

   if ( !man->fUserOutputGenFileName.empty() ) {
      userdefgen_file.open(man->fUserOutputGenFileName.c_str());
      if ( !userdefgen_file.good() ) {
         printf("[Error %s: Line %d] %s File Not Found.\n",__FILE__,__LINE__,man->fUserOutputGenFileName.c_str());
         return 10;
      }
      getline(userdefgen_file,comment);//read comment;
   }

   // --------------------------------------------------------
   //
   TFile* f = new TFile(samc_rootfilename.c_str(),"recreate");
   TTree* T = new TTree("SAMC","Tree with Acceptance Simulation");

   SAMCEvent* Event   = new SAMCEvent();
   T->Branch("samcEvent","SAMCEvent",&Event);


   TStopwatch timer;
   timer.Start();

   //Set the Seed using CPU Time, but not fix the seed,
   //to avoid same random values in different runs.
   gRandom->SetSeed(0);

   //only init rfunction once
   left_init_r_function_();
   right_init_r_function_();

   fail_events=0;

   bool   IsDebug             = man->IsDebug            ;//= false; // output to debugfile
   double E0                  = man->E0;
   int    Which_Beam_Profile  = man->Which_Beam_Profile;
   double delta_dp            = man->delta_dp           ;// d(dp) dp full width for generator
   double delta_th            = man->delta_th           ;// d(th) th full width for generator(tan(th))
   double delta_ph            = man->delta_ph           ;// d(ph) th full width for generator(tan(ph))
   double gaus_x_sigma        = man->gaus_x_sigma       ;// if Which_Beam_Profile == 1/2 sigma of x beam for generator(cm)
   double gaus_y_sigma        = man->gaus_y_sigma       ;// if Which_Beam_Profile == 1/2 sigma of y beam for generator(cm)
   double raster_x_size       = man->raster_x_size      ;// raster x full size for generator(cm)
   double raster_y_size       = man->raster_y_size      ;// raster y full size for generator(cm)
   double beam_x_center       = man->beam_x_center      ;// beam x center for generator(cm)
   double beam_y_center       = man->beam_y_center      ;// beam y center for generator(cm)
   double z0                  = man->z0                 ;// target center for generator(cm)
   double T_L                 = man->T_L                ;// target length for generator(cm)
   double T_H                 = man->T_H                ;// target height for generator(cm)

   // Create the event generator and track propagator
   SAMCEventGenerator event_generator = SAMCEventGenerator();
   SAMCPropagator     propagator      = SAMCPropagator();
   event_generator.Print();

   // -------------------------------
   // Event Loop
   for(int iEvent = 0; iEvent < Num_Of_Events; iEvent++ ) {

      // TODO: need to clean up event class
      Event->Clear();
      Event->Win_Before_Mag.clear();
      Event->Win_After_Mag.clear();

      Event->Id          = iEvent;
      Event->E_s         = gRandom->Gaus(E0,E0*3e-5);//beam dispersion is 3e-5
      Event->theta       = man->fTheta;//atof(inputdata[j++].c_str());

      // Why is this being set every event? -whit
      Event->Target     = man->fTargetMaterial;
      Event->Win_i      = man->fMat0;
      Event->Win_f      = man->fMat1;
      Event->T_theta    = man->fTheta_Target;

      Event->AddOneSAMCMaterial( Event->Win_Before_Mag, man->fMat2 );
      Event->AddOneSAMCMaterial( Event->Win_Before_Mag, man->fMat3 );
      Event->AddOneSAMCMaterial( Event->Win_Before_Mag, man->fMat4 );

      Event->AddOneSAMCMaterial( Event->Win_After_Mag, man->fMat5 );
      Event->AddOneSAMCMaterial( Event->Win_After_Mag, man->fMat6 );

      // Here there are two event generators
      // 1. Predefined user events read from a file
      // 2. Random events thrown here.
      if ( !man->fUserOutputGenFileName.empty() ) {
         int buf_i = 0;
         for ( buf_i = 0; buf_i < 6; ++buf_i ) {
            userdefgen_file>>usrgen[buf_i];
         }
         Event->beam_x     = usrgen[0];
         Event->beam_y     = usrgen[1];
         Event->reactz_gen = usrgen[2]+z0;//Add offset from the seed, for HRS-R only, Z. Ye 07/26/2013
         Event->th_tg_gen  = usrgen[3];
         Event->ph_tg_gen  = usrgen[4];
         Event->dp_gen     = usrgen[5];
      } else {

         event_generator.GenerateEvent(*Event);

      }
 
      propagator.InitTrack(*Event);

      err = propagator.PropagateTrack(*Event);

      if ( err ) { exit(err); }
      if ( IsDebug && err==0 ) {
         Event->Print();
      }

      T->Fill();

      //i-=err;//if err>0 This event is bad, simulate again. To make sure total passed number of events == one user inputs.
      if ( Event->IsPassed==0 && man->fUserOutputGenFileName.empty() ) {
         ++fail_events; //but also record the bad one to calculate acceptance
         //				--i;
      }
      if ( ( iEvent+1 )%1000==0 || iEvent==0 || iEvent==(Num_Of_Events-1) ) {
         if ( man->fUserOutputGenFileName.empty() ) {
            printf("\t%d Good Events Simulated. And %d Bad Events Saved.\n",iEvent+1,fail_events);
         }
         else {
            printf("\t%d User-Defined Events Simulated.\n",iEvent+1);
         }
      }

   } // end of event loop

   if ( !man->fUserOutputGenFileName.empty() ) {
      userdefgen_file.close();
   }

   T->Write();
   delete T;
   f->Close();
   delete f;

   printf("Results saved at %s\n",samc_rootfilename.c_str());
   printf("Time at the end of job = %f seconds\n",timer.CpuTime());

   return 0;
}
//______________________________________________________________________________

