#include "SAMCAnalyzerOptics.h"

//______________________________________________________________________________
SAMCAnalyzerOptics::SAMCAnalyzerOptics(){
   fOpticsPackage = 0;
   fP0            = -1;
   fScatAngle     = 0;
}
//_________________________________________________________________________________
SAMCAnalyzerOptics::~SAMCAnalyzerOptics(){
   ClearVectors();
}
//_________________________________________________________________________________
void SAMCAnalyzerOptics::ClearVectors(){
   fFP.clear();
   fY.clear();
   fD.clear();
   fP.clear();
   fT.clear();
   fYTA.clear();
   fPTA.clear();
   fL.clear();
}
//_________________________________________________________________________________
void SAMCAnalyzerOptics::SetOpticsPackage(int choice){

   fOpticsPackage = choice;

   TString OpticsName;

   switch(fOpticsPackage){
      case 0: // J LeRose 
         OpticsName = Form("J. LeRose's");
         break;
      case 1: // E06-014 
         OpticsName = Form("E06-014");
         break;
      case 2: // CSR
         OpticsName = Form("CSR");
         break;
      default: 
         std::cout << "[SAMC]: Invalid optics package!  Did you forget to set it?  Exiting..." << std::endl;
         exit(1);
   }

   std::cout << "[SAMCAnalyzerOptics]: Using " << OpticsName << " optics matrix for target variables" << std::endl;

   if(fOpticsPackage==2){
      // for CSR only 
      LoadData(1);
   }else{
      for(int i=1;i<4;i++) LoadData(i); 
   }
}
//_________________________________________________________________________________
void SAMCAnalyzerOptics::LoadData(int FileIndex){

   // import the general matrix elements
   // y,t,p       => transport to focal plane
   // Y,T,P,D     => focal plane to target  
   // L           => pathlength from z=0 (target) to focal plane in meters
   // XF,YF,TF,PF => forward: target to focal plane (?) 

   SAMCMatrixElement ME; // matrix element object for storing data  
   ME.SetPolyOrder(0);

   std::string iVar;
   int         SIZE;
   int         iOpt,iDEG,iI,iJ,iK;
   double      iVal1,iVal2,iVal3,iVal4,iVal5,iVal6,iVal7; 
   std::vector<int> E;
   std::vector<double> V;

   // for CSR optics
   // choose DB file based on momentum  
   int MyP0; 

   if(fP0>=450){
      MyP0 = 1102;
   }else if(fP0>0){
      MyP0 = 399;
   }else{
      std::cout << "[SAMCAnalyzerOptics::LoadData]: Invalid momentum setting!  Is it set properly?  Exiting..." << std::endl;
      exit(1);
   }

   std::string fn;
   std::string prefix = "./input/optics/";
   char p0[2048]; 
   sprintf(p0,"%d",MyP0);
   std::string p0_str = (std::string)p0; 

   switch(fOpticsPackage){
      case 0: // J. LeRose 
         std::cout << "[SAMCAnalyzerOptics]: Loading E06-014 optics anyway..." << std::endl;
         prefix += "E06014/";
         break; 
      case 1: // E06-014 
         prefix += "E06014/";
         break; 
      case 2: // CSR 
         prefix += "CSR/configs/" + p0_str + "/";
         break;
      default:
         std::cout << "[SAMCAnalyzerOptics]: WARNING: Invalid optics package!  Exiting..." << std::endl;
         exit(1); 
   }

   std::string fprefix = "NONE"; 
   if(fScatAngle>0) fprefix = "db_L"; 
   if(fScatAngle<0) fprefix = "db_R";
   if(fScatAngle==0){
      std::cout <<"[SAMCAnalyzerOptics::LoadData]: Invalid scattering angle!  Is it set?  Exiting..." << std::endl;
      exit(1);
   }

   if(FileIndex==1){
      // y,t,p,Y,T,P,D
      fn = prefix + fprefix + ".vdc_1.dat"; 
      ME.SetExpoSize(3); 
   }else if(FileIndex==2){
      // L
      fn = prefix + fprefix + ".vdc_2.dat"; 
      ME.SetExpoSize(4); 
   }else if(FileIndex==3){
      // XF,YF,TF,PF
      fn = prefix + fprefix + ".vdc_3.dat"; 
      ME.SetExpoSize(5); 
   }

   std::ifstream infile; 
   infile.open(fn.c_str()); 
   if(infile.fail()){
      std::cout << "[SAMCAnalyzerOptics]: Cannot open the file: " << fn << std::endl;
      exit(1);
   }else{
      std::cout << "[SAMCAnalyzerOptics]: Loading data for file: " << fn << std::endl;
      while(!infile.eof()){
         if(FileIndex==1){
            infile >> iVar  >> iI    >> iJ    >> iK 
               >> iVal1 >> iVal2 >> iVal3 >> iVal4 >> iVal5 >> iVal6 >> iVal7 
               >> iDEG;  
            E.push_back(iI);
            E.push_back(iJ);
            E.push_back(iK);
            V.push_back(iVal1);
            V.push_back(iVal2);
            V.push_back(iVal3);
            V.push_back(iVal4);
            V.push_back(iVal5);
            V.push_back(iVal6);
            V.push_back(iVal7);
         }else if(FileIndex==2){
            infile >> iVar >> iOpt >> iI >> iJ >> iK >> iVal1;   
            E.push_back(iOpt);
            E.push_back(iI);
            E.push_back(iJ);
            E.push_back(iK);
            V.push_back(iVal1);
         }else if(FileIndex==3){
            infile >> iVar >> iOpt >> iI >> iJ >> iK >> iDEG >> iVal1;   
            E.push_back(iOpt);
            E.push_back(iI);
            E.push_back(iJ);
            E.push_back(iK);
            E.push_back(iDEG);
            V.push_back(iVal1);
         }
         SIZE = E.size();
         ME.SetExpoSize(SIZE); 
         for(int i=0;i<SIZE;i++){
            ME.SetExpo(i,E[i]);
         }
         SIZE = V.size();  
         ME.SetPolySize(SIZE); 
         for(int i=0;i<SIZE;i++){
            ME.SetPoly(i,V[i]);
            if(V[i]!=0){
               ME.SetBoolean(false);
               ME.SetPolyOrder(i+1);
            } 
         }
         if( ME.IsElementZero() ){
            SIZE = E.size();
            std::cout << "[SAMCAnalyzerOptics]: Skipping element: " << iVar << "\t";
            for(int i=0;i<SIZE;i++){
               std::cout << E[i] << "\t";
            }
            std::cout << std::endl;
            ME.Clear();
            V.clear();
            E.clear();
            continue;
         }
         // push back on the right matrix element vector 
         if(iVar=="t" || iVar=="y" || iVar=="p") fFP.push_back(ME); 
         if(iVar=="Y")   fY.push_back(ME); 
         if(iVar=="D")   fD.push_back(ME); 
         if(iVar=="T")   fT.push_back(ME); 
         if(iVar=="P")   fP.push_back(ME); 
         if(iVar=="YTA") fYTA.push_back(ME); 
         if(iVar=="PTA") fPTA.push_back(ME); 
         if(iVar=="L")   fL.push_back(ME); 
         ME.Clear();
         V.clear();
         E.clear();
      }
      infile.close();
      std::cout << "[SAMCAnalyzerOptics]: Data load successful." << std::endl;
   }

}
//_________________________________________________________________________________
void SAMCAnalyzerOptics::CalculateTargetCoords(std::vector<double> FPVar,std::vector<double> &TgVar){
   // calculates the target coordinates from the focal plane coordinates 

   const int kNUM_PRECOMP_POW = 10;         // wtf? 
   double x_fp=0,y_fp=0,th_fp=0,ph_fp=0;    // focal plane variables 
   double y_tg=0,th_tg=0,ph_tg=0,dp=0;      // target variables
   double powers[kNUM_PRECOMP_POW][5];

   // clear the target vector (to be safe) 
   TgVar.clear();

   // fill focal plane variables.  I think the standard units are meters
   x_fp  = FPVar[0]/100.; // convert to meters  
   y_fp  = FPVar[1]/100.; // convert to meters
   th_fp = FPVar[2]; 
   ph_fp = FPVar[3];

   // calculate the powers we need
   for(int i=0;i<kNUM_PRECOMP_POW;i++){
      powers[i][0] = pow(x_fp,i);
      powers[i][1] = pow(th_fp,i);
      powers[i][2] = pow(y_fp,i);
      powers[i][3] = pow(ph_fp,i);
      powers[i][4] = pow(fabs(th_fp),i);
   }

   // calculate the matrices we need
   CalculateMatrix(x_fp,fD);
   CalculateMatrix(x_fp,fT);
   CalculateMatrix(x_fp,fY);
   CalculateMatrix(x_fp,fP);
   // CalculateMatrix(x_fp,fYTA);
   // CalculateMatrix(x_fp,fPTA);

   // calculate the coordinates at the target
   // FIXME: what about the PTA, YTA matrix elements?? They are EMPTY! 
   th_tg = CalculateTargetVar(fT,powers);
   ph_tg = CalculateTargetVar(fP,powers);
   // ph_tg = CalculateTargetVar(fP,powers) + CalculateTargetVar(fPTA,powers);
   // y_tg  = CalculateTargetVar(fY,powers) + CalculateTargetVar(fYTA,powers);
   y_tg  = CalculateTargetVar(fY,powers);  
   dp    = CalculateTargetVar(fD,powers);         

   // fill the vector of target variables.  note the order!  
   y_tg *= 100.;           // convert back to cm 
   TgVar.push_back(y_tg);  
   TgVar.push_back(th_tg);
   TgVar.push_back(ph_tg);
   TgVar.push_back(dp); 

}
//_____________________________________________________________________________
void SAMCAnalyzerOptics::CalculateFocalPlaneCoords(std::vector<double> DCSVar,std::vector<double> &FPVar){
   // calculates focal plane coordinates from detector coordinates

   double tan_rho, cos_rho, tan_rho_loc, cos_rho_loc;
   // TRANSPORT coordinates (projected to z=0)
   double x,y,th,ph;
   // Rotating TRANSPORT coordinates
   double r_x,r_y,r_th,r_ph;

   // clear the focal plane vector (to be safe) 
   FPVar.clear(); 

   // NOTE: - Optics matrix elements are in meters.  We convert
   //         x and y values to meters first.  At the end, we
   //         convert back to cm.  

   // FIXME: I think the input from SAMC is in fact TCS variables!  
   //        Need to translate to DCS! 

   // FIXME: I'm fairly certain the input vector should be the DCS variables. 
   double d_x  = DCSVar[0]/100.;  // convert to meters
   double d_y  = DCSVar[1]/100.;  // convert to meters
   double d_th = DCSVar[2];
   double d_ph = DCSVar[3];

   // tan rho (for the central ray) is stored as a matrix element 
   tan_rho = fFP[T000].GetPoly(0);
   cos_rho = 1.0/sqrt(1.0+tan_rho*tan_rho);

   // Analyzer's THaTrack methods: 
   // GetDX     => x position in DCS 
   // GetDY     => y position in DCS 
   // GetDTheta => tangent of DCS theta
   // GetDPhi   => tangent of DCS phi 

   // first calculate the transport frame coordinates
   // th = (track->GetDTheta()+tan_rho) /
   // 	(1.0-track->GetDTheta()*tan_rho);
   // x  = track->GetDX() * cos_rho * (1.0 + theta * tan_rho);
   // ph = track->GetDPhi() / (1.0-track->GetDTheta()*tan_rho) / cos_rho;
   // y  = track->GetDY() + tan_rho*phi*track->GetDX()*cos_rho;

   // first calculate the transport frame coordinates
   th = (d_th + tan_rho)/(1.0 - d_th*tan_rho);
   x  = d_x*cos_rho*(1.0 + th*tan_rho);
   ph = d_ph/(1.0 - d_th*tan_rho)/cos_rho;
   y  = d_y + d_x*ph*tan_rho*cos_rho;

   if(x==0) std::cout << "[SAMCAnalyzerOptics::CalculateFocalPlaneCoords]: x_fp = 0!" << std::endl;

   // then calculate the rotating transport frame coordinates
   r_x = x;

   // calculate the focal-plane matrix elements
   // if(mode == kTransport){
   // 	CalculateMatrix( x, fFPMatrixElems );
   // }else if (mode == kRotatingTransport){
   // 	CalculateMatrix( r_x, fFPMatrixElems );
   // }

   // FIXME: I believe that the focal plane is in the rotated coordinate system 
   CalculateMatrix(r_x,fFP);

   r_y = y - fFP[Y000].GetValue();       // Y000

   // Calculate now the tan(rho) and cos(rho) of the local rotation angle.
   tan_rho_loc = fFP[T000].GetValue();   // T000
   cos_rho_loc = 1.0/sqrt(1.0 + tan_rho_loc*tan_rho_loc);

   r_ph = (d_ph - fFP[P000].GetValue()/* P000 */)/(1.0 - d_th*tan_rho_loc)/cos_rho_loc;
   r_th = (d_th + tan_rho_loc)/(1.0 - d_th*tan_rho_loc);

   // set the values we calculated
   // track->Set(x, y, theta, phi);
   // track->SetR(r_x, r_y, r_theta, r_phi);

   // convert x and y back to cm
   r_x *= 100.;
   r_y *= 100.;
   FPVar.push_back(r_x); 
   FPVar.push_back(r_y); 
   FPVar.push_back(r_th); 
   FPVar.push_back(r_ph); 

}
//_________________________________________________________________________________
void SAMCAnalyzerOptics::CalculateMatrix(const double x,std::vector<SAMCMatrixElement>& matrix){
   // calculates the values of the matrix elements for a given location
   // by evaluating a polynomial in x of order it->order with 
   // coefficients given by it->poly

   double arg=0; 
   int size = matrix.size(); 

   if(size>0){
      for( std::vector<SAMCMatrixElement>::iterator it=matrix.begin();
            it!=matrix.end(); it++ ){
         it->SetValue(0.0);
         if( it->GetPolyOrder() > 0){
            for(int i=it->GetPolyOrder()-1; i>=1; i--){
               arg = x*( it->GetValue() + it->GetPoly(i) ); 
               it->SetValue(arg);
            }
            arg = it->GetValue() + it->GetPoly(0);
            it->SetValue(arg);
         }
      }
   }

}
//_____________________________________________________________________________
double SAMCAnalyzerOptics::CalculateTargetVar(std::vector<SAMCMatrixElement>& matrix,
      const double powers[][5]){
   // calculates the value of a variable at the target
   // the x-dependence is already in the matrix, so only 1-3 (or np) used
   double retval= 0.0;
   double v     = 0.0;
   int size = matrix.size(); 

   if(size>0){
      for( std::vector<SAMCMatrixElement>::iterator it=matrix.begin();
            it!=matrix.end(); it++ ){ 
         if(it->GetValue() != 0.0){
            v = it->GetValue();
            unsigned int np = it->GetExpoSize(); // generalize for extra matrix elems.
            for(unsigned int i=0; i<np; i++){
               v *= powers[it->GetExpo(i)][i+1];
            }
            retval += v;
            //      retval += it->v * powers[it->pw[0]][1] 
            //	              * powers[it->pw[1]][2]
            //	              * powers[it->pw[2]][3];
         }
      }
   }

   return retval;

}
//_________________________________________________________________________________
void SAMCAnalyzerOptics::PrintMatrixElement(std::string type){

   int SIZE = 0;
   std::vector<SAMCMatrixElement> MyMatrixElement; 
   if(type=="FP"){
      SIZE = fFP.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fFP[i]);
   }else if(type=="Y"){
      SIZE = fY.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fY[i]);
   }else if(type=="D"){
      SIZE = fD.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fD[i]);
   }else if(type=="T"){
      SIZE = fT.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fT[i]);
   }else if(type=="P"){
      SIZE = fP.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fP[i]);
   }else if(type=="YTA"){
      SIZE = fYTA.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fYTA[i]);
   }else if(type=="PTA"){
      SIZE = fPTA.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fPTA[i]);
   }else if(type=="L"){
      SIZE = fL.size();
      for(int i=0;i<SIZE;i++) MyMatrixElement.push_back(fL[i]);
   }

   std::cout << type << " matrix elements: " << std::endl;
   for(int i=0;i<SIZE;i++) MyMatrixElement[i].Print();


}
