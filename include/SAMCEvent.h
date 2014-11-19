#ifndef SAMCEvent_HH 
#define SAMCEvent_HH 1

#include <string>
#include <vector>
#include "SAMCMaterial.h"
#include "SAMCFortran.h"
#include "TLorentzVector.h"
#include "TVector3.h"

//#include "Your_Cross_Section.h" //The file you store your own cross section models, call it when caluclating cs_Final

#define MSIZE 5

/** Event class.
 *  This class needs to shed a bunch of member functions.
 * 
 *  HCS: Hall Coordinate System
 *  TCS: Target Coordinate System
 *
 */
class SAMCEvent {

   public:
      SAMCEvent();
      virtual ~SAMCEvent() ;
      SAMCEvent(const SAMCEvent &);
      SAMCEvent& operator=(const SAMCEvent &);

      int Process() ;

      void AddOneSAMCMaterial(std::vector<SAMCMaterial>& aWin, const SAMCMaterial& );
      void AddOneSAMCMaterial(std::vector<SAMCMaterial>& aWin,
                              const double& aX0,
                              const double& arho,
                              const double& aL,
                              const double& aA,
                              const int& aZ,
                              std::string aName) ;

   private:
 
      //generate original target variables
      void Generator();

      //transport to front of magnetic, then back to get refined target variables
      int RefineTg() ;

      //transfer to focal plane using John.LeRose matrix
      int ToFp(const double& ax,const double& ay,const double& ath,const double& aph,const double& adp) ;

      int ReconstructTg(const double& ax,const double& ay,const double& ath,const double& aph,const double& axtg) ;

      /*Bool_t IntersectPlaneWithRay( const TVector3& xax,{{{*/
      Bool_t IntersectPlaneWithRay( const TVector3& xax,
                                    const TVector3& yax,
                                    const TVector3& org,
                                    const TVector3& ray_start,
                                    const TVector3& ray_vect,
                                    Double_t& length,
                                    TVector3& intersect ) ;

   private:

      double Ion_Loss(const double& aE0,const SAMCMaterial& aSAMCMaterial) ;

      double Bremss_Loss(const double& aE0,const double& abt) ;

      double eta(const int& aZ); 
      double b(const int& aZ) ;
      double Rad_Len(const int& aZ,const double& aA) ;
      void   Transport(TLorentzVector& aPos,TLorentzVector& aMom,const SAMCMaterial& aSAMCMaterial,const bool& aIsMultiScatt=false) ;
      double MultiScattering(const double& aE,const double& aTR) ;
      void   SetSAMCMaterial(SAMCMaterial& aSAMCMaterial) ;
      SAMCMaterial GetMixture(const std::vector<SAMCMaterial>& aWin) ;
      void GetRef_Plane(TLorentzVector& aoPos,TLorentzVector& aoMom,const std::vector<SAMCMaterial>& aWinBefore,const std::vector<SAMCMaterial>& aWinAfter,const double& aL,const double& aOffset) ;
      double sigma_M(const double& aE,const double& aTheta) ;

      double CalcRValue(const double& ath,const double& aph,const double& ay,const double& adp) ;

      double PROD_AND(const double& ax,const double& ay) ;


   public:

      void Clear(){
         IsPassed    = 1;
         IsQualified = 1;
      }

      /*Member Data from file{{{*/
      int                       Id;             // Event Id
      int                       IsPassed;       // if pass through the magnet, 1=true, 0=false(not pass)
      int                       IsQualified;    // if rec var is in the range of gen, 1=true, 0=false(not qualified)
      double                    theta;          // scattering angle(deg)
      SAMCMaterial              Target;         // target SAMCMaterial
      SAMCMaterial              Win_i;          // front[inital] window, stick with target
      SAMCMaterial              Win_f;          // back[final] window, stick with target
      std::vector<SAMCMaterial> Win_Before_Mag; // windows after target and before magnetic
      std::vector<SAMCMaterial> Win_After_Mag;  // windows after target and after magnetic
      double                    T_theta;        // target angle(deg), the central line of target is along beam
      //----
      //|  |
      //|  |--->central line, so the T_theta is defined like phi_tg
      //|  |
      //----

      //TLorentzVector ( vec, t ) = ( k,E ) in ROOT
      //HCS: Hall Coordinate System
      //TCS: Target Coordinate System
      TLorentzVector s;                     // s(s3,E_s) 4-momentum position of the incident electron in HCS
      TLorentzVector target_edgepoint_TRCS; // 4-momentum position of the edge point in HCS (if theta>0,T_H/2,if <0,-T_H/2)
      TLorentzVector s_TCS;                 // s_TCS(s3,E_s) 4-momentum position of the incident electron in TCS, and ztg=0
      TLorentzVector p_TCS;                 // p_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=0
      TLorentzVector p_inter_point_TCS;     // p_inter_point_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=reactz_TCS
      TLorentzVector p_P;                   // p_P 4-momentum momentum of the outgoing electron in HCS
      TLorentzVector p_P_TCS;               // p_P_TCS 4-momentum momentum of the outgoing electron in TCS

      double Angle;       // real scattering angel, not from file,rad
      double Angle_Deg;   // real scattering angel, not from file,Deg
      double sinsq_Angle; // sin(Angle/2)^2
      double theta_rad;   // scattering angle in SAMC
      double sinsq;       // sin(theta/2)^2
      double reactz_TCS;  //
      double btr;         // b*tr (equivalent radiator)
      double Q2;          // MeV*MeV
      double q2;          // q2=-Q2 MeV*MeV

      /*Member Data from random{{{*/
      double beam_x;     // cm
      double beam_y;     // cm
      double reactz_gen; // interaction z position cm
      double E_s;        // =incident beam energy (MeV)
      double E_p;        // =incident beam energy (MeV)
      double th_tg_gen; // theta target
      double ph_tg_gen; // phi target
      double dp_gen;    // dp at target
      /*}}}*/

      /*Member Data derived from variables above{{{*/
      double x_tg_gen; //
      double y_tg_gen; //
      /*}}}*/

      /*Member Data for refinement{{{*/
      double dp_ref; //dp for John.LeRose matrix
      double x_tg_ref; //x tg for John.LeRose matrix
      double y_tg_ref; //y tg for John.LeRose matrix
      double th_tg_ref; //th tg for John.LeRose matrix
      double ph_tg_ref; //ph tg for John.LeRose matrix
      TLorentzVector p_TCS_ref;//
      TLorentzVector p_P_TCS_ref;//
      /*}}}*/

      /*Member Data for focal plane {{{*/
      TLorentzVector p_FP;             // p_FP 4-momentum momentum of the outgoing electron in focal plane
      TLorentzVector p_P_FP;           // p_P_FP 4-momentum momentum of the outgoing electron in focal plane
      double         x_fp;             //
      double         y_fp;             //
      double         th_fp;            //
      double         ph_fp;            //
      double         th_fp_no_ort;     // th_fp without TXFIT orthogonalization
      double         q1ex[2];          // 0:x 1:y
      double         dent[2];
      double         dext[2];
      double         q3en[2];
      double         q3ex[2];
      bool           IsPassedQ1Ex;
      bool           IsPassedDipoleEn;
      bool           IsPassedDipoleEx;
      bool           IsPassedQ3En;
      bool           IsPassedQ3Ex;
      /*}}}*/

      /*Member Data for reconstruct target varibles {{{*/
      double x_tg_rec; //cm
      double y_tg_rec; //cm
      double th_tg_rec; //
      double ph_tg_rec; //
      double dp_rec;
      double reactz_rec;
      double rvalue;
      double cs_M;//mott cross section
      double cs_Final;//final cross section; 
      double Angle_rec;//Calculated from reconstructed variables,
      double Qsq; //Calculated from reconstructed variables, Q2 is from the generated variables 
      double Xbj;//Calculated from reconstructed variables,
      /*}}}*/

      void Print() ; // *MENU* 
      void PrintMaterial(const SAMCMaterial& aMaterial) ;

   ClassDef(SAMCEvent,1)
};
#endif

