#ifndef SAMCPropagator_HH
#define SAMCPropagator_HH 1

#include <vector>
#include "TNamed.h"
#include "SAMCEvent.h"
#include "SAMCMaterial.h"
#include "SAMCManager.h"

/** Propator class which sets the various event data members
 *  as it propagates forward and backwards through the detector.
 *
 */
class SAMCPropagator : public TNamed {

   public:

      SAMCMaterial              fTarget_mat;        // target SAMCMaterial
      SAMCMaterial              fTargetWin_i_mat;   // front[inital] window, stick with target
      SAMCMaterial              fTargetWin_f_mat;   // back[final] window, stick with target
      std::vector<SAMCMaterial> fMaterials_before;  // windows after target and before magnetic
      std::vector<SAMCMaterial> fMaterials_after;   // windows after target and after magnetic
      double                    fT_theta;           // target angle(deg), the central line of target is along beam

   public:
      SAMCPropagator();
      virtual ~SAMCPropagator();

      void InitTrack(SAMCEvent& event);

      Int_t PropagateTrack(SAMCEvent& event);

      //transport to front of magnetic, then back to get refined target variables
      int RefineTg(SAMCEvent& event) ;

      //transfer to focal plane using John.LeRose matrix
      int ToFp(SAMCEvent& event,const double& ax,const double& ay,const double& ath,const double& aph,const double& adp) ;

      int ReconstructTg(SAMCEvent& event,const double& ax,const double& ay,const double& ath,const double& aph,const double& axtg) ;

      void Print(Option_t * opt =""); // *MENU*

      double Ion_Loss(const double& aE0,const SAMCMaterial& aSAMCMaterial) ;
      double Bremss_Loss(const double& aE0,const double& abt) ;
      double eta(const int& aZ); 
      double b(const int& aZ) ;
      double Rad_Len(const int& aZ,const double& aA) ;
      void   Transport(TLorentzVector& aPos,TLorentzVector& aMom,const SAMCMaterial& aSAMCMaterial,const bool& aIsMultiScatt=false) ;
      double MultiScattering(const double& aE,const double& aTR) ;
      void   SetSAMCMaterial(SAMCMaterial& aSAMCMaterial) ;
      SAMCMaterial GetMixture(const std::vector<SAMCMaterial>& aWin) ;
      void   GetRef_Plane(TLorentzVector& aoPos,TLorentzVector& aoMom,const std::vector<SAMCMaterial>& aWinBefore,const std::vector<SAMCMaterial>& aWinAfter,const double& aL,const double& aOffset) ;
      double sigma_M(const double& aE,const double& aTheta) ;
      double CalcRValue(const double& ath,const double& aph,const double& ay,const double& adp) ;
      double PROD_AND(const double& ax,const double& ay) ;

      void AddOneSAMCMaterial(std::vector<SAMCMaterial>& aWin, const SAMCMaterial& );
      void AddOneSAMCMaterial(std::vector<SAMCMaterial>& aWin,
                              const double& aX0,
                              const double& arho,
                              const double& aL,
                              const double& aA,
                              const int& aZ,
                              std::string aName) ;


   ClassDef(SAMCPropagator,0)
};

#endif

