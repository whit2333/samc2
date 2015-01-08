#include "SAMCPropagator.h"
#include "SAMCConstants.h"
#include "SAMCFortran.h"


//___________________________________________________________________
SAMCPropagator::SAMCPropagator(){
}
//___________________________________________________________________
SAMCPropagator::~SAMCPropagator(){
}
//___________________________________________________________________
void SAMCPropagator::Print(Option_t * opt){
}
//______________________________________________________________________________
void SAMCPropagator::AddOneSAMCMaterial(std::vector<SAMCMaterial>& aWin,const SAMCMaterial& mat) {
   // A temporary work around function to get things to compile.
   // Ideally this would not be int the event class!!!!
   AddOneSAMCMaterial(aWin, mat.fX0, mat.fDensity, mat.fL, mat.fA, mat.fZ, std::string(mat.GetName()) );
}
//______________________________________________________________________________
void SAMCPropagator::AddOneSAMCMaterial(std::vector<SAMCMaterial>& aWin,const double& aX0,const double& arho,const double& aL,const double& aA,const int& aZ, std::string aName) {

   //To make sure read correctly, so I reverse the order compared with input file.
   SAMCMaterial a;
   a.SetName( aName.c_str() );
   a.fZ    = aZ;
   a.fA    = aA;
   a.fL    = aL;
   a.fDensity  = arho;
   a.fT    = aL*arho;
   a.fX0   = aX0;
   if ( fabs(aX0)<1e-10 )
   {
      a.fTR=0;
      a.fbt=0;
   }
   else
   {
      a.fTR = a.fT/aX0;
      a.fbt = b(aZ)*aL*arho/aX0;
   }
   aWin.push_back(a);
}
//___________________________________________________________________
Int_t SAMCPropagator::PropagateTrack(SAMCEvent& event) {
   int Num_Event_Add = 0;

   SetSAMCMaterial(event.Target);
   SetSAMCMaterial(event.Win_i);
   SetSAMCMaterial(event.Win_f);

   //Generator();                // generate original target variables

   Num_Event_Add = RefineTg(event); // transport to front of magnetic, then back to get refined target variables

   if ( Num_Event_Add > 0 ) {
      event.IsPassed    = 0;
      event.IsQualified = 0;
   } else {
      //transfer to focal plane using John.LeRose matrix
      Num_Event_Add = ToFp(event,event.x_tg_ref,event.y_tg_ref,event.th_tg_ref,event.ph_tg_ref,event.dp_ref);
   }

   if ( Num_Event_Add>0 ) {
      event.IsPassed    = 0;
      event.IsQualified = 0;
   } else {
      Num_Event_Add = ReconstructTg(event,event.x_fp,event.y_fp,event.th_fp,event.ph_fp,event.x_tg_ref);
   }

   return 0;
}

//___________________________________________________________________
void SAMCPropagator::InitTrack(SAMCEvent& event) {

   //set value for Member Data derived from variables from file
   //File provides E_s,theta,Target.(Z,A,T,rho) Win_i.(Z,A,T,rho)
   //Win_f.(Z,A,T,rho) T_theta
   //Win_Before_Mag(Name,Z,A,L,rho,X0)
   //Win_After_Mag(Name,Z,A,L,rho,X0)
   //beam_x,beam_y,reactz_gen,th_tg_gen,ph_tg_gen,dp_gen
   //z0,HRS_L,VDC_Res(x,y,th,ph),D_(x,y),T_L,P0
   //IsMultiScat,IsEnergyLoss,Which_Kin,FP_Eff_L

   /*Set SAMCMaterial.(Z,A,M,X0,T,TR,bt){{{*/
   //SetSAMCMaterial(Target);
   //SetSAMCMaterial(Win_i);
   //SetSAMCMaterial(Win_f);
   //Target.Print();
   //Win_i.Print();
   //Win_f.Print();

   /*}}}*/

   SAMCManager * man = SAMCManager::Instance();

   /*Set Beam Info(HCS and TCS){{{*/
   //know beam_x,beam_y,reactz_gen,E_s,theta,HRS_L
   //Set s,s_TCS,x_tg_gen,y_tg_gen,p_TCS,p_P,p_P_TCS

   event.theta_rad                = event.theta*DegToRad();                         // rad
   event.s(0)                     = event.beam_x;                                   // cm
   event.s(1)                     = event.beam_y;                                   // cm
   event.reactz_gen              += -(event.beam_x)*TMath::Tan(man->fTheta_Target*DegToRad()); //
   event.s(2)                     = event.reactz_gen;                               // cm
   event.target_edgepoint_TRCS(0) = event.theta/TMath::Abs(event.theta)*man->T_H/2;             // if theta>0,T_H/2, if<0, -T_H/2 in Target Rotation Coordinate System(T_theta = 0) not TCS, check Coordinate.svg
   event.target_edgepoint_TRCS(1) = 0;                                        //
   event.target_edgepoint_TRCS(2) = man->T_L/2;                               // and z0= 0

   TLorentzVector lp_TRCS;//the interaction point in TRCS at T_theta=0
   lp_TRCS     = event.s;
   lp_TRCS(2) -= man->z0;
   //Printf("s(%g,%g,%g),target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",s(0),s(1),s(2),target_edgepoint_TRCS(0),target_edgepoint_TRCS(1),target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
   lp_TRCS.RotateY(-man->fTheta_Target*DegToRad());

   event.s_TCS = event.s;
   event.s_TCS.RotateZ(PI/2);//passive ratation around HCS, so -(-PI/2)
   event.s_TCS.RotateX(event.theta_rad);//passive ratation around HCS, so -(-theta_rad)
   event.s_TCS(0)-=man->D_x;
   event.s_TCS(1)-=man->D_y;

   event.p_inter_point_TCS = event.s_TCS;
   event.reactz_TCS        = event.s_TCS.Z();

   event.s_TCS(0) -= event.s_TCS.Z()*event.th_tg_gen;
   event.s_TCS(1) -= event.s_TCS.Z()*event.ph_tg_gen;
   event.s_TCS(2) -= event.s_TCS.Z();

   event.x_tg_gen = event.s_TCS(0);
   event.y_tg_gen = event.s_TCS(1);

   //Fix me: I don't think we need this
   //SAMCMaterial halftarget=Target;
   //halftarget.T /= 2;
   //E_s-=Ion_Loss(E_s,Win_i);
   //E_s-=Bremss_Loss(E_s,Win_i.bt+btr);
   //E_s-=Ion_Loss(E_s,halftarget);
   //E_s-=Bremss_Loss(E_s,halftarget.bt+btr);
   event.s(3)     = event.E_s;//MeV
   event.s_TCS(3) = event.E_s;//MeV

   event.p_TCS.SetVect(event.s_TCS.Vect());
   TLorentzVector lz(0,0,1,0);//in HCS
   event.p_P = lz;
   event.p_P.RotateZ(PI/2);
   event.p_P.RotateX(event.theta_rad);
   event.p_P.RotateY(-atan(event.th_tg_gen));
   event.p_P.RotateX(atan(event.ph_tg_gen));
   event.Angle       = event.p_P.Angle(lz.Vect());
   event.Angle_Deg   = event.Angle*RadToDeg();
   event.sinsq_Angle = sin(event.Angle/2)*sin(event.Angle/2);
   event.sinsq       = sin(event.Angle/2)*sin(event.Angle/2);//sin(Angle/2)^2
   //p_P_TCS=lz;//Now think it's in TCS
   //p_P_TCS.RotateY(atan(th_tg_gen));//p_P_TCS.X()/p_P_TCS.Z()=th_tg
   //p_P_TCS.RotateX(-atan(ph_tg_gen));//p_P_TCS.Y()/p_P_TCS.Z()=ph_tg
   event.p_P_TCS.SetX(event.th_tg_gen);
   event.p_P_TCS.SetY(event.ph_tg_gen);
   event.p_P_TCS.SetZ(1);
   event.p_P_TCS.SetVect(event.p_P_TCS.Vect().Unit());
   /*}}}*/

   size_t i;
   size_t imax;
   SAMCMaterial m0;
   /*Set Windows{{{*/
   //know Win_Before_Mag,Win_f,HRS_L,theta
   //Add Target,Win_f,Air if necessary to Win_Before_Mag
   //Win_After_Mag don't need to be changed since it's vacuum.
   //for Ion_Loss Z,A,T,rho
   //for Bremss_Loss bt
   //for MultiScattering TR
   //for Transport L
   //Thickness definition diagram
   //Win_i |<--T-->|

   //Target|<--T-->|

   //Win_f ---------
   //      |<--T
   //      ---------

   //Target diagram,z0,TCS,reactz_gen,T
   //|<----T--->|
   //|  !  !   !|
   //   TCS!   !
   //      z0  !
   //          reactz_gen
   //          <-lL
   std::vector<SAMCMaterial>::iterator it = event.Win_Before_Mag.begin();//here only have the SAMCMaterial before mag defined in input file.
   double lHL                   = 0; //if lHL<HRS_Lcm or <3.57m(See Hall_A_NIM), i assume it's air
   int FirstInWinBeforeMagBlock = 0;
   imax                         = event.Win_Before_Mag.size();

   for ( i=0; i<imax; i++ ) {
      lHL += event.Win_Before_Mag[i].fL;
   }

   double lphrad = event.theta_rad - man->fTheta_Target*DegToRad();//scattering angle in y axis in TCS or x axis in TRCS
   //it doesn't I add atan(ph_tg_gen) since it's so small for target and windows

   //Win_Before_Mag add Win_f
   if( (event.Win_f.fZ != 0) && (event.Win_f.fA!= 0) && (event.Win_f.fDensity != 0) && (event.Win_f.fTR != 0) )
   {
      m0 = event.Win_f;//Win_f added to Win_Before_Mag
      m0.SetName("Win_f");
      event.Win_Before_Mag.insert(it,m0);
      FirstInWinBeforeMagBlock++;
   }
   //Win_Before_Mag add Target
   it = event.Win_Before_Mag.begin();

   m0 = event.Target;//target after interaction point added to Win_Before_Mag
   double lL; //assume target is rentangle lL=distance between reactz_gen and edge of target
   //in TRCS x plane<--> y plane in TCS
   lL = event.target_edgepoint_TRCS(2) - lp_TRCS(2);//verticle distance between interaction point and edge plane of target
   if ( lL > man->T_L ) {
      Printf("target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",event.target_edgepoint_TRCS(0),event.target_edgepoint_TRCS(1),event.target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
   }
   lL*=tan(lphrad);
   lL+=lp_TRCS(0);//x on edge

   if ( fabs(lL)<fabs(event.target_edgepoint_TRCS(0)) ) {//HCS 0=x top view
      //in the target
      lL=fabs((event.target_edgepoint_TRCS(2)-lp_TRCS(2))/cos(lphrad));
   }
   else {
      //out of target before hiting the edge of target, the edge means the downstream face.
      lL=fabs((event.target_edgepoint_TRCS(0)-lp_TRCS(0))/sin(lphrad));
   }

   m0.fL=lL;
   event.Win_Before_Mag.insert(it,m0);
   FirstInWinBeforeMagBlock++;

   //Win_Before_Mag add Air before last SAMCMaterial in the input file
   imax = event.Win_Before_Mag.size();
   if ( (lHL)<(man->HRS_L) )
   {
      //http://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/104.html
      double total = 0.000124+0.755267+0.231781+0.012827;//total mass of component of air
      m0.SetName("Air");
      m0.fZ = int(6*0.000124/total+7*0.755267/total+8*0.231781/total+18*0.012827/total);
      m0.fA = m0.fZ/0.499;
      m0.fL = man->HRS_L-lHL;

      m0.fDensity  = 1.205e-03;
      m0.fT        = m0.fL*m0.fDensity;
      m0.fX0       = 36.66;
      m0.fbt       = event.b(m0.fZ)*m0.fT/m0.fX0;
      it           = event.Win_Before_Mag.end();
      it--;
      event.Win_Before_Mag.insert(it,m0);

   } else if ( (lHL)>man->HRS_L ) {
      printf("[Error %s: Line %d] Total windows length=%f>HRS_L=%f.\n",__FILE__,__LINE__,(lHL+event.Win_Before_Mag[imax-1].fL),man->HRS_L);
      exit(-3);
   }

   //Correct First SAMCMaterial.L in Win_Before_Mag Block of inputfile
   it = event.Win_Before_Mag.begin();
   (it+FirstInWinBeforeMagBlock)->fL += event.reactz_TCS;
   for ( i=0; i<FirstInWinBeforeMagBlock; i++ )
      (it+FirstInWinBeforeMagBlock)->fL-=(it+i)->fL;

   for ( i = 0; i < event.Win_Before_Mag.size(); ++i ) {

      event.Win_Before_Mag[i].fT = event.Win_Before_Mag[i].fL*event.Win_Before_Mag[i].fDensity;

      if ( fabs(event.Win_Before_Mag[i].fX0)<1e-10 ) {
         event.Win_Before_Mag[i].fTR=0;
      }
      else {
         event.Win_Before_Mag[i].fTR = event.Win_Before_Mag[i].fT/event.Win_Before_Mag[i].fX0;
      }
      event.Win_Before_Mag[i].fbt =  event.b(event.Win_Before_Mag[i].fZ)*event.Win_Before_Mag[i].fTR;
   }
   /*}}}*/

   /*Set Target Info(TCS){{{*/
   //know P0,dp_gen,Which_Kin,IsEnergyLoss,E_s,E_p,sinsq,Target
   //Set dp_gen,s_TCS,x_tg_gen,y_tg_gen,p_TCS,p_P,p_P_TCS
   //Set Q2,q2,btr,Win_Before_Mag[0].bt
   //x_tg_gen,y_tg_gen,see Set Beam Info(HCS and TCS)
   switch ( man->Which_Kin )
   {
      case 1: //elastic
         event.E_p   = gRandom->Gaus(0,man->P0*3e-4);
         event.E_p  += event.E_s/(1+2*event.E_s*event.sinsq/event.Target.fM);
         break;
      case 2: //quasi-elastic
      case 0: //phase default
      default:
         event.E_p=man->P0*(1+event.dp_gen);
   }
   event.Q2  = 4*event.E_s*event.E_p*event.sinsq;
   event.q2  = -event.Q2;
   event.btr = AP*(log(event.Q2/(ELECTRON_MASS*ELECTRON_MASS))-1);//b*t_r
   ////Correct Win_Before_Mag[0](Target) bt
   //Win_Before_Mag[0].bt+=btr;

   event.p_TCS(3)   = event.E_p;//MeV
   event.p_P(3)     = event.E_p;
   event.p_P_TCS(3) = event.E_p;
   event.dp_gen     = (event.E_p-man->P0)/man->P0;
   /*}}}*/

}
//______________________________________________________________________________
int SAMCPropagator::RefineTg(SAMCEvent& event) {

   //return Number of Event needs to be added. just 1 ^_^
   //refine target variables from Generator because John.LeRose matrix only works for vacuum
   //know Win_Before_Mag,p_TCS,p_P_TCS
   event.p_TCS_ref   = event.p_inter_point_TCS;
   event.p_P_TCS_ref = event.p_P_TCS;

   SAMCManager * man = SAMCManager::Instance();
   size_t i;
   size_t imax;
   std::vector<SAMCMaterial> Win_Empty;//tmp use

   SAMCMaterial mixture;

   double offset = event.p_TCS(2)-event.p_TCS_ref(2);//L=along Z in TCS
   mixture.fL     = offset;

   //Printf("p_TCS(%g,%g,%g),(th=%g,ph=%g)",p_TCS(0),p_TCS(1),p_TCS(2),th_tg_gen,ph_tg_gen);
   //FIXME: I want to add energy loss in GetRef_Plane
   //but the SAMC will run forever and no good event.
   if ( event.E_p<ELECTRON_MASS ) {
      printf("[Warning %s: Line %d] E_p=%g<ELECTRON_MASS=%g\n",__FILE__,__LINE__,event.E_p,ELECTRON_MASS);
      return 1;
   }

   GetRef_Plane(event.p_TCS_ref,event.p_P_TCS_ref,event.Win_Before_Mag,Win_Empty,0,offset);

   event.x_tg_ref  = event.p_TCS_ref(0);
   event.y_tg_ref  = event.p_TCS_ref(1);
   event.th_tg_ref = event.p_P_TCS_ref(0)/event.p_P_TCS_ref(2);
   event.ph_tg_ref = event.p_P_TCS_ref(1)/event.p_P_TCS_ref(2);
   //Printf("p_TCS_ref(%g,%g,%g),(th=%g,ph=%g)",p_TCS_ref(0),p_TCS_ref(1),p_TCS_ref(2),th_tg_ref,ph_tg_ref);

   event.dp_ref = event.dp_gen;
   event.E_p    = event.p_P_TCS_ref(3);

   if ( man->IsEnergyLoss )
   {
      //Fix me: Do I need to count the energy loss due to windows after target?
      //I think so. So I add them.
      imax = event.Win_Before_Mag.size();
      for ( i=0; i<imax; i++ )
      {
         event.E_p -= Ion_Loss(event.E_p,event.Win_Before_Mag[i]);
         event.E_p -= Bremss_Loss(event.E_p,event.Win_Before_Mag[i].fbt);
      }
      imax=event.Win_After_Mag.size();
      for ( i=0; i<imax; i++ )
      {
         event.E_p-=Ion_Loss(event.E_p,event.Win_After_Mag[i]);
         event.E_p-=Bremss_Loss(event.E_p,event.Win_After_Mag[i].fbt);
      }
      if ( event.E_p<0 ) {
         printf("[Warning %s: Line %d] E_p after Energy Loss=%f<0.\n",__FILE__,__LINE__,event.E_p);
         return 1;
      }
      event.Q2             = 4*event.E_s*event.E_p*event.sinsq;
      event.q2             = -event.Q2;
      event.p_TCS_ref(3)   = event.E_p;//MeV
      event.p_P_TCS_ref(3) = event.E_p;
      event.dp_ref         = (event.E_p-man->P0)/man->P0;
   }


   event.cs_M     = sigma_M(event.E_s, event.Angle_Deg);
   event.cs_Final = event.cs_M;//Just use the Mott cross section so far
   //////////////////////////////////////////////////////////////////
   // This is where you add your own cross section.
   // e.g.:
   // cs_Final = Your_Cross_Section(E_s, E_p, Angle_Deg, A, Z);
   //////////////////////////////////////////////////////////////////
   return 0;
}
//______________________________________________________________________________
int SAMCPropagator::ToFp(SAMCEvent& event, const double& ax,const double& ay,const double& ath,const double& aph,const double& adp)
{
   SAMCManager * man = SAMCManager::Instance();
   //return Number of Event needs to be added. just 1 ^_^
   //cm->meter, ath,aph,adp are correct
   //must be float, otherwise cannot pass to fortran correctly. float<->dimension
   float  matrix[MSIZE] = {ax/100.,ath,ay/100,aph,adp};
   int    msize         = MSIZE;
   double xtest = 0;
   double ytest = 0;
   unsigned int i = 0;
   for ( i = 0; i < 2; ++i ) {
      event.q1ex[i]=0;
      event.dent[i]=0;
      event.dext[i]=0;
      event.q3en[i]=0;
      event.q3ex[i]=0;
   }
   event.IsPassedQ1Ex     = false;
   event.IsPassedDipoleEn = false;
   event.IsPassedDipoleEx = false;
   event.IsPassedQ3En     = false;
   event.IsPassedQ3Ex     = false;
   if ( event.theta>0 )
   {
      //for left arm
      /*Q1 Exit{{{*/
      xtest   = x_e_q1ex_(matrix,msize)*100; //cm
      ytest   = y_e_q1ex_(matrix,msize)*100; //cm
      event.q1ex[0] = xtest/100;
      event.q1ex[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>man->Q1_Radius*man->Q1_Radius )
      {
         //(xtets!=xtest)==true to avoid nan(Not a number)
         //printf("Blocked by Q1 Exit.\n");
         return 1;
      }
      event.IsPassedQ1Ex=true;
      /*}}}*/
      /*Dipole Entrance{{{*/
      // Transport electron to dipole entrance, trapezoid define by (jjl)
      // -40cm<x<40cm (+x is down in HRS frame)
      // y=+-(12.5*(1-(1.25*x/840)) (smallest gap at positive x, x in cm)
      xtest   = x_e_dent_(matrix,msize)*100; //cm
      ytest   = y_e_dent_(matrix,msize)*100; //cm
      event.dent[0] = xtest/100;
      event.dent[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)> man->D_X_Radius || fabs(ytest)>man->D_Y_L*(1-1.25*xtest/840) ) //nan
      {
         //printf("Blocked by Dipole Entrance.\n");
         return 1;
      }
      event.IsPassedDipoleEn=true;
      /*}}}*/

      /*Dipole Exit{{{*/
      xtest   = x_e_dext_(matrix,msize)*100; //cm
      ytest   = y_e_dext_(matrix,msize)*100; //cm
      event.dext[0] = xtest/100;
      event.dext[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>man->D_X_Radius || fabs(ytest)>man->D_Y_L*(1-1.25*xtest/840) ) //nan
      {
         //printf("Blocked by Dipole Exit.\n");
         return 1;
      }
      event.IsPassedDipoleEx=true;
      /*}}}*/

      /*Q3 Entrance{{{*/
      xtest   = x_e_q3en_(matrix,msize)*100; //cm
      ytest   = y_e_q3en_(matrix,msize)*100; //cm
      event.q3en[0] = xtest/100;
      event.q3en[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>man->Q3_Entrance_Radius*man->Q3_Entrance_Radius )
      {
         //printf("Blocked by Q3 Entrance.\n");
         return 1;
      }
      event.IsPassedQ3En=true;

      /*}}}*/

      /*Q3 Exit{{{*/
      xtest   = x_e_q3ex_(matrix,msize)*100; //cm
      ytest   = y_e_q3ex_(matrix,msize)*100; //cm
      event.q3ex[0] = xtest/100;
      event.q3ex[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>man->Q3_Exit_Radius*man->Q3_Exit_Radius )
      {
         //printf("Blocked by Q3 Exit.\n");
         return 1;
      }
      event.IsPassedQ3Ex=true;
      /*}}}*/

      event.x_fp  = x_e_fp_(matrix,msize)*100.; //cm
      event.y_fp  = y_e_fp_(matrix,msize)*100.; //cm
      event.th_fp = t_e_fp_(matrix,msize); //(tantheta)
      event.ph_fp = p_e_fp_(matrix,msize); //(tanphi)
   }
   else
   {
      //for right arm
      /*Q1 Exit{{{*/
      xtest   = x_h_q1ex_(matrix,msize)*100; //cm
      ytest   = y_h_q1ex_(matrix,msize)*100; //cm
      event.q1ex[0] = xtest/100;
      event.q1ex[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>man->Q1_Radius*man->Q1_Radius )
      {
         //(xtets!=xtest)==true to avoid nan(Not a number)
         //printf("Blocked by Q1 Exit.\n");
         return 1;
      }
      event.IsPassedQ1Ex=true;
      /*}}}*/

      /*Dipole Entrance{{{*/
      // Transport electron to dipole entrance, trapezoid define by (jjl)
      // -40cm<x<40cm (+x is down in HRS frame)
      // y=+-(D_Y_L*(1-(1.25*x/840)) (smallest gap at positive x, x in cm)
      xtest   = x_h_dent_(matrix,msize)*100; //cm
      ytest   = y_h_dent_(matrix,msize)*100; //cm
      event.dent[0] = xtest/100;
      event.dent[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>man->D_X_Radius || fabs(ytest)>man->D_Y_L*(1-1.25*xtest/840) ) //nan
      {
         //printf("Blocked by Dipole Entrance.\n");
         return 1;
      }
      event.IsPassedDipoleEn=true;
      /*}}}*/

      /*Dipole Exit{{{*/
      xtest   = x_h_dext_(matrix,msize)*100; //cm
      ytest   = y_h_dext_(matrix,msize)*100; //cm
      event.dext[0] = xtest/100;
      event.dext[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>man->D_X_Radius || fabs(ytest)>man->D_Y_L*(1-1.25*xtest/840) ) //nan
      {
         //printf("Blocked by Dipole Exit.\n");
         return 1;
      }
      event.IsPassedDipoleEx=true;
      /*}}}*/

      /*Q3 Entrance{{{*/
      xtest   = x_h_q3en_(matrix,msize)*100; //cm
      ytest   = y_h_q3en_(matrix,msize)*100; //cm
      event.q3en[0] = xtest/100;
      event.q3en[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>man->Q3_Entrance_Radius*man->Q3_Entrance_Radius )
      {
         //printf("Blocked by Q3 Entrance.\n");
         return 1;
      }
      event.IsPassedQ3En=true;
      /*}}}*/

      /*Q3 Exit{{{*/
      xtest   = x_h_q3ex_(matrix,msize)*100; //cm
      ytest   = y_h_q3ex_(matrix,msize)*100; //cm
      event.q3ex[0] = xtest/100;
      event.q3ex[1] = ytest/100;
      if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>man->Q3_Exit_Radius*man->Q3_Exit_Radius )
      {
         //printf("Blocked by Q3 Exit.\n");
         return 1;
      }
      event.IsPassedQ3Ex=true;
      /*}}}*/

      event.x_fp  = x_h_fp_(matrix,msize)*100.; //cm
      event.y_fp  = y_h_fp_(matrix,msize)*100.; //cm
      event.th_fp = t_h_fp_(matrix,msize); //(tantheta)
      event.ph_fp = p_h_fp_(matrix,msize); //(tanphi)
   }

   if ( man->IsMultiScat) {

      //SAMCMaterial mixture=GetMixture(Win_After_Mag);
      //mixture.TR*=sqrt(ph_fp*ph_fp+th_fp*th_fp);
      //mixture.L=0;
      //p_FP.SetXYZT(x_fp,y_fp,0,p_TCS_ref(3));
      //p_P_FP.SetXYZT(th_fp,ph_fp,1,p_TCS_ref(3));
      //Transport(p_FP,p_P_FP,mixture,IsMultiScat);//Just change th_fp,ph_fp
      //th_fp=p_P_FP(0)/p_P_FP(2);
      //ph_fp=p_P_FP(1)/p_P_FP(2);

      //SAMCMaterial mixture=GetMixture(Win_After_Mag);
      event.p_FP.SetXYZT(  event.x_fp ,event.y_fp ,0,event.p_TCS_ref(3));
      event.p_P_FP.SetXYZT(event.th_fp,event.ph_fp,1,event.p_TCS_ref(3));
      SAMCMaterial Vacuum;
      Vacuum.fL=0;
      for ( i = 0; i < event.Win_After_Mag.size(); ++i ) {
         Vacuum.fL-=   event.Win_After_Mag[i].fL;
      }
      Vacuum.SetName("Vacuum");
      Vacuum.fTR=0;
      Transport(event.p_FP,event.p_P_FP,Vacuum,man->IsMultiScat);//Just change th_fp,ph_fp
      for ( i = 0; i < event.Win_After_Mag.size(); ++i ) {
         Transport(event.p_FP,event.p_P_FP,event.Win_After_Mag[i],man->IsMultiScat);//Just change th_fp,ph_fp
      }
      //Transport(p_FP,p_P_FP,mixture,IsMultiScat);//Just change th_fp,ph_fp
      event.th_fp = event.p_P_FP(0)/event.p_P_FP(2);
      event.ph_fp = event.p_P_FP(1)/event.p_P_FP(2);
      event.x_fp  = event.p_FP(0);
      event.y_fp  = event.p_FP(1);
   }

   //Smearing
   TRandom* tr  = new TRandom();
   event.x_fp         = tr->Gaus(event.x_fp,man->VDC_Res_x);
   event.y_fp         = tr->Gaus(event.y_fp,man->VDC_Res_y);
   event.th_fp        = tr->Gaus(event.th_fp,man->VDC_Res_th/1000.);
   event.ph_fp        = tr->Gaus(event.ph_fp,man->VDC_Res_ph/1000.);
   event.th_fp_no_ort = event.th_fp;
   delete tr;

   matrix[0] = event.x_fp/100.;
   msize     = 1;

   //orthogonalize theta see John.LeRose Webpage
   //http://hallaweb.jlab.org/news/minutes/tranferfuncs.html
   if ( event.theta>0 )
   {
      event.th_fp-=ltxfit_(matrix,msize);
   } else {
      event.th_fp-=rtxfit_(matrix,msize);
   }
   return 0;
}
//______________________________________________________________________________
int SAMCPropagator::ReconstructTg(SAMCEvent& event,const double& ax,const double& ay,const double& ath,const double& aph,const double& axtg)
{
   //return Number of Event needs to be added. just 1 ^_^
   float matrix[MSIZE]={ax/100.,ath,ay/100,aph,axtg/100};
   int msize=MSIZE;

   SAMCManager * man = SAMCManager::Instance();

   event.x_tg_rec = axtg;
   float rf_y,rf_d,rf_th,rf_ph;
   if ( event.theta>0 )
   {
      //printf("%f %f %f %f %f\n",matrix[0],matrix[1],matrix[2],matrix[3],matrix[4]);
      event.y_tg_rec  = ly00_(matrix,msize)*100;
      event.th_tg_rec = ltheta_(matrix,msize);
      event.ph_tg_rec = lphi_(matrix,msize);
      event.dp_rec    = ldelta_(matrix,msize);
      rf_y      = event.y_tg_rec/100;
      rf_d      = event.dp_rec;
      rf_th     = event.th_tg_rec;
      rf_ph     = event.ph_tg_rec;
      if ( man->fNCuts <= 0 ) {
         event.rvalue = left_rfunction_(&rf_y,&rf_d,&rf_th,&rf_ph);
      }
   } else {
      event.y_tg_rec  = ry00_(matrix,msize)*100;
      event.th_tg_rec = rtheta_(matrix,msize);
      event.ph_tg_rec = rphi_(matrix,msize);
      event.dp_rec    = rdelta_(matrix,msize);
      rf_y      = event.y_tg_rec/100;
      rf_d      = event.dp_rec;
      rf_th     = event.th_tg_rec;
      rf_ph     = event.ph_tg_rec;
      if ( man->fNCuts<=0 ) {
         event.rvalue = right_rfunction_(&rf_y,&rf_d,&rf_th,&rf_ph);
      }
   }
   if ( man->fNCuts>0 ) {
      event.rvalue = CalcRValue(rf_th,rf_ph,rf_y,rf_d);
   }
   TLorentzVector rec(event.x_tg_rec,event.y_tg_rec,0,man->P0*(1+event.dp_rec));

   rec(0)  +=  man->D_x;
   rec(1)  +=  man->D_y;
   rec(0)  +=  event.reactz_TCS*event.th_tg_rec;
   rec(1)  +=  event.reactz_TCS*event.ph_tg_rec;
   rec(2)  +=  event.reactz_TCS;
   rec.RotateX(-event.theta_rad);
   rec.RotateZ(-PI/2);
   event.reactz_rec=rec(2);

   event.Angle_rec = acos( (cos(event.theta_rad)-event.ph_tg_rec*sin(event.theta_rad)) 
         / sqrt(1.0+pow(event.th_tg_rec,2)+pow(event.ph_tg_rec,2)) );
   // cerr <<"New Angle is " << Angle_rec*RadToDeg() <<endl; 
   event.Qsq = 4.0*event.E_s*event.E_p*sin(event.Angle_rec/2.0)*sin(event.Angle_rec/2.0);
   //	cerr <<"Qsq is " << Qsq <<endl;
   event.Xbj = event.Qsq/2.0/(event.E_s-event.E_p)/PROTON_MASS;
   //	cerr <<"Xbj is " << Xbj <<endl;

   if (  //fabs(reactz_rec/100)<0.1 &&
         //fabs(y_tg_rec/100)<=0.01 &&
         fabs(event.dp_rec)    < man->delta_dp/2 &&
         fabs(event.th_tg_rec) < man->delta_th/2 &&
         fabs(event.ph_tg_rec) < man->delta_ph/2
         //fabs(ph_tg_rec)<delta_ph/2
      ){
      event.IsQualified = 1;
   } else {
      event.IsQualified = 0;
   }
   return 0;
}
//______________________________________________________________________________
//___________________________________________________________________
//______________________________________________________________________________
double SAMCPropagator::Ion_Loss(const double& aE0,const SAMCMaterial& aSAMCMaterial)
{
   //aT: g/cm^2, arho: g/cm^3
   //Particle Booklet Equ(27.9)
   //printf("Z=%d,A=%f,T=%f,rho=%f\n",aSAMCMaterial.Z,aSAMCMaterial.A,aSAMCMaterial.T,aSAMCMaterial.rho);
   double lK       = 0.307075;// cm^2/g for A                                                                                = 1 g/mol
   double lbetasq  = 1-ELECTRON_MASS*ELECTRON_MASS/(aE0*aE0);
   double lxi      = lK/2*aSAMCMaterial.fZ/aSAMCMaterial.fA*aSAMCMaterial.fT/lbetasq;//aT: g/cm^2
   double lhbarwsq = 28.816*28.816*aSAMCMaterial.fDensity*aSAMCMaterial.fZ/aSAMCMaterial.fA*1e-12;//MeV arho is density of absorber
   double j        = 0.200;
   double Delta_p  = lxi*(log(2*ELECTRON_MASS*lxi/lhbarwsq)+j);
   double lw       = 4*lxi;
   double result   = 0;
   if ( aSAMCMaterial.fZ!=0 && aSAMCMaterial.fA!=0 && aSAMCMaterial.fT!=0 && aSAMCMaterial.fDensity !=0 )
      result=gRandom->Landau(Delta_p,lw);
   if ( result>(aE0-ELECTRON_MASS) )
      result=aE0-ELECTRON_MASS;
   if ( result<0 )
      result=0;
   return result;
}
//______________________________________________________________________________
double SAMCPropagator::Bremss_Loss(const double& aE0,const double& abt)
{
   //Bremsstrahlung Energy Loss for external and internal(equivalent radiator)
   //Xiaodong Jiang, PhD.thesis Equ 5.15
   //http://filburt.mit.edu/oops/Html/Pub/theses/xjiang.ps
   //*0.999 to avoid lose all energy
   double result=0;
   if ( abt!=0 )
      result=aE0*pow(gRandom->Rndm()*0.999,1./abt);
   if ( result>(aE0-ELECTRON_MASS) )
      result=aE0-ELECTRON_MASS;
   if ( result<0 )
      result=0;
   return result;
}
//______________________________________________________________________________
double SAMCPropagator::eta(const int& aZ)
{
   //Phys.Rev.D 12,1884 A46
   return log(1440*pow(aZ,-2/3.))/log(183*pow(aZ,-1/3.));
}
//______________________________________________________________________________
double SAMCPropagator::b(const int& aZ)
{
   //Phys.Rev.D 12,1884 A45
   if ( aZ!=0 )
      return 4./3.*( 1+1./9.*(aZ+1)/(aZ+eta(aZ))/log(183*pow(aZ,-1/3.)) );
   else
      return 0;
}
//______________________________________________________________________________
double SAMCPropagator::Rad_Len(const int& aZ,const double& aA)
{
   //particle book equation 27.20
   //Lrad=elastic form factor F_el, scattering on the nucleus
   //Lrad_prime=inelastic form factor F_inel, scattering on the shell electrons
   //f(Z)=Coulomb correction
   double Lrad,Lrad_prime,f_Z;
   if ( aZ==1 )
   {
      Lrad=5.31;
      Lrad_prime=6.144;
   }
   else if ( aZ==2 )
   {
      Lrad=4.79;
      Lrad_prime=5.621;
   }
   else if ( aZ==3 )
   {
      Lrad=4.74;
      Lrad_prime=5.805;
   }
   else if ( aZ==4 )
   {
      Lrad=4.71;
      Lrad_prime=5.924;
   }
   else
   {
      Lrad=log(184.15*pow(aZ,-1./3.));
      Lrad_prime=log(1194.*pow(aZ,-2./3.));
   }
   double a=ALPHA*aZ;
   a=a*a;
   f_Z=a*(1./(1+a)+0.20206-0.0369*a+0.0083*a*a-0.002*a*a*a);
   if ( aZ!=0 && aA!=0 )
      return 716.408*aA/(aZ*aZ*(Lrad-f_Z)+aZ*Lrad_prime);
   else
      return 0;
}
//______________________________________________________________________________
void SAMCPropagator::Transport(TLorentzVector& aPos,
                          TLorentzVector& aMom,
                          const SAMCMaterial& aSAMCMaterial,
                          const bool& aIsMultiScatt)
{
   //need aSAMCMaterial.TR and aSAMCMaterial.L
   //aPos and aMom are inputs, also outputs
   //aPos: position aMom: momentum
   double lE = aMom(3);
   double ms_phi,ms_theta;
   if ( aIsMultiScatt && fabs(aSAMCMaterial.fTR)>1e-15 )
   {
      ms_phi   = MultiScattering(lE,aSAMCMaterial.fTR);//rad
      ms_theta = MultiScattering(lE,aSAMCMaterial.fTR);//rad
   } else {
      ms_phi   = 0;
      ms_theta = 0;
   }
   //ms_theta=2*PI*gRandom->Rndm();
   double lth = aMom.X()/aMom.Z();
   double lph = aMom.Y()/aMom.Z();

   //pass L/2
   aPos.SetX(aPos.X()+aSAMCMaterial.fL/2*lth);
   aPos.SetY(aPos.Y()+aSAMCMaterial.fL/2*lph);
   aPos.SetZ(aPos.Z()+aSAMCMaterial.fL/2);

   //change the angle and pass the rest L/2
   lth = (tan(ms_theta)+lth)/(1-tan(ms_theta)*lth);
   lph = (tan(ms_phi)+lph)/(1-tan(ms_phi)*lph);
   aMom.SetX(lth);
   aMom.SetY(lph);
   aMom.SetZ(1);
   aPos.SetX(aPos.X()+aSAMCMaterial.fL/2*lth);
   aPos.SetY(aPos.Y()+aSAMCMaterial.fL/2*lph);
   aPos.SetZ(aPos.Z()+aSAMCMaterial.fL/2);
}
//______________________________________________________________________________
double SAMCPropagator::MultiScattering(const double& aE,const double& aTR) {
   //only for electron
   double lPsq    = aE*aE-ELECTRON_MASS*ELECTRON_MASS;
   double bcp     = lPsq/aE;
   double ltheta0 = 13.6/bcp*sqrt(aTR)*(1+0.038*log(aTR));
   if ( aTR!=0 )
   {
      //return gRandom->Gaus(0,ltheta0/2.3548);//rad sigma=width/(2*sqrt(2)
      //ltheta/w, sg=sigma of y_tg in data
      //w=1.461e6*sg^2-6976*sg+9.316
      return gRandom->Gaus(0,ltheta0/1.3548);//rad sigma=width because of the
      //calculation of energy loss is behind the multiscattering
      //I add more multiscattering
   }
   return 0;
}
//______________________________________________________________________________
void SAMCPropagator::SetSAMCMaterial(SAMCMaterial& aSAMCMaterial)
{
   aSAMCMaterial.fM = aSAMCMaterial.fA*AMU;//MeV
   if( aSAMCMaterial.fL == 0 && aSAMCMaterial.fDensity != 0 ) {
      aSAMCMaterial.fL = aSAMCMaterial.fT/aSAMCMaterial.fDensity;
   }
   aSAMCMaterial.fX0 = Rad_Len(aSAMCMaterial.fZ,aSAMCMaterial.fA);
   if( aSAMCMaterial.fX0 != 0 ) {
      aSAMCMaterial.fTR = aSAMCMaterial.fT/aSAMCMaterial.fX0;
   } else {
      aSAMCMaterial.fTR = 0;
   }
   aSAMCMaterial.fbt = b(aSAMCMaterial.fZ)*aSAMCMaterial.fTR;
}
//______________________________________________________________________________
SAMCMaterial SAMCPropagator::GetMixture(const std::vector<SAMCMaterial>& aWin)
{
   size_t i;
   size_t imax=aWin.size();
   SAMCMaterial mixture;
   mixture.SetName("mixture");
   mixture.fL=0;
   mixture.fA=0;
   //get mixture TR
   for ( i=0; i<imax; i++ )
   {
      mixture.fL+=aWin[i].fL;
      mixture.fA+=aWin[i].fA;
   }
   mixture.fTR=0;
   for ( i=0; i<imax; i++ )
   {
      mixture.fTR+=aWin[i].fA/mixture.fA*aWin[i].fTR;
   }
   return mixture;
}
//______________________________________________________________________________
void SAMCPropagator::GetRef_Plane(TLorentzVector& aoPos,
                             TLorentzVector& aoMom,
                             const std::vector<SAMCMaterial>& aWinBefore,
                             const std::vector<SAMCMaterial>& aWinAfter,
                             const double& aL,
                             const double& aOffset)
{
   //Get Refinement aoPos and aoMom for each plane on q1,q2,d,q3,fp
   //ao means input and output
   //aL=Vacuum Length, aOffset=distance between interaction point and z=0 in TCS
   std::vector<SAMCMaterial> AllWins;
   SAMCMaterial Vacuum;
   SAMCManager * man = SAMCManager::Instance();

   AllWins.clear();
   AllWins=aWinBefore;
   Vacuum.SetName("Vacuum");
   Vacuum.fL=aL;
   Vacuum.fA=0;
   Vacuum.fTR=0;
   AllWins.push_back(Vacuum);
   size_t i;
   for ( i=0; i<aWinAfter.size(); i++ )
      AllWins.push_back(aWinAfter[i]);
   Vacuum.fL=0;
   //Printf("Pos(%g,%g,%g),(th=%g,ph=%g)",aoPos(0),aoPos(1),aoPos(2),aoMom(0)/aoMom(2),aoMom(1)/aoMom(2));
   for ( i = 0; i < AllWins.size(); ++i ) {
      Transport(aoPos,aoMom,AllWins[i],man->IsMultiScat);//move to mag
      //Printf("Pos(%g,%g,%g),(th=%g,ph=%g)",aoPos(0),aoPos(1),aoPos(2),aoMom(0)/aoMom(2),aoMom(1)/aoMom(2));
      Vacuum.fL-=AllWins[i].fL;
   }
   Vacuum.fL+=aOffset;
   Transport(aoPos,aoMom,Vacuum,man->IsMultiScat);//back to ztg=0 in TCS
   //Printf("Pos(%g,%g,%g),(th=%g,ph=%g)",aoPos(0),aoPos(1),aoPos(2),aoMom(0)/aoMom(2),aoMom(1)/aoMom(2));
}
//______________________________________________________________________________
double SAMCPropagator::sigma_M(const double& aE,const double& aTheta)
{
   //Mott Cross Section
   //aE=MeV,aTheta=deg
   double ltheta=aTheta/2.*PI/180;
   double mott;
   if ( ltheta!=0 )
      mott=pow(ALPHA*cos(ltheta)/(2*aE*pow(sin(ltheta),2)),2);
   else
      mott=0;
   return mott*MEV2SR_TO_NBARNSR; //nbarn
   //return mott; //nbarn
}
//______________________________________________________________________________
double SAMCPropagator::CalcRValue(const double& ath,const double& aph,const double& ay,const double& adp)
{
   int i,j,k;
   SAMCManager * man = SAMCManager::Instance();

   i=0;
   double* lxy=new double[NELEMENTS];
   lxy[i++]=ath;
   lxy[i++]=aph;
   lxy[i++]=ay;
   lxy[i++]=adp;
   double prod=-1000;
   if ( man->fNCuts>0 ) {
      double lnormal;
      double lomega;
      for ( i = 0; i < man->fNCuts; ++i ) {
         //slope*x+intersection-y=0

         lnormal = sqrt((man->fLineProperty[i][LINE_SLOPE]) * (man->fLineProperty[i][LINE_SLOPE]+1));//Get normalized factor (sqrt(a*a+b*b)) a = slope b = -1
         lomega  = 0;
         k       = 0;

         for ( j = 0; j < NELEMENTS; ++j ) {
            lomega += man->fXY[i][j]*pow(-1.0,man->fXY[i][j]+k)*pow(man->fLineProperty[i][LINE_SLOPE],k)*lxy[j];
            //if ( i==1 ) {
            //	Printf("i=%d,j=%d,k=%d,lomega=%g",i,j,k,lomega);
            //	Printf("fXY[%d][%d]=%d,k=%d,fLineProperty[%d][%d]=%g,lxy[%d]=%g,lomega=%d*pow(-1,%d)*pow(%g,%d)*%g=%g",i,j,fXY[i][j],k,i,LINE_SLOPE,fLineProperty[i][LINE_SLOPE],j,lxy[j],fXY[i][j],fXY[i][j]+k,fLineProperty[i][LINE_SLOPE],k,lxy[j],fXY[i][j]*pow(-1.0,fXY[i][j]+k)*pow(fLineProperty[i][LINE_SLOPE],k)*lxy[j]);

            //}
            k += man->fXY[i][j];
         }
         lomega += man->fLineProperty[i][LINE_INTERSECTION];
         lomega /= lnormal*man->fLineProperty[i][LINE_SIGN];
         //Printf("lomega[%d]=%g",i,lomega);
         if ( i==0 ) {
            prod = lomega;
         }
         else {
            prod=PROD_AND(prod,lomega);
         }
      }
   }

   delete [] lxy;
   return prod;
}
//______________________________________________________________________________
double SAMCPropagator::PROD_AND(const double& ax,const double& ay)
{
   double prod;
   prod=TMath::Min(ax,ay);
   //prod=ax+ay-sqrt(ax*ax+ay*ay);
   return prod;
}
//______________________________________________________________________________



