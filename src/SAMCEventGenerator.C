#include "SAMCEventGenerator.h"

//______________________________________________________________________________
SAMCEventGenerator::SAMCEventGenerator(){
   SetNameTitle("SAMCEventGenerator","standard samc event generator");
   SAMCManager * man = SAMCManager::Instance();
   fBeamProfileType  = man->Which_Beam_Profile;
   fRasterXSize      = man->raster_x_size;
   fRasterYSize      = man->raster_y_size;
   fGausXSigma       = man->gaus_x_sigma;
   fGausYSigma       = man->gaus_y_sigma;
   fBeamXCenter      = man->beam_x_center;
   fBeamYCenter      = man->beam_y_center;
   fDelta_th         = man->delta_th;
   fDelta_ph         = man->delta_ph;
   fDelta_dp         = man->delta_dp;
   fz0               = man->z0;
   fT_L              = man->T_L;
}
//______________________________________________________________________________
SAMCEventGenerator::~SAMCEventGenerator(){
}
//______________________________________________________________________________
void SAMCEventGenerator::Print() {
   std::cout << GetName() << std::endl;
   std::cout << GetTitle() << std::endl;
   std::cout << "  fBeamProfileType : " << fBeamProfileType << std::endl;
   std::cout << "      fRasterXSize : " << fRasterXSize << std::endl;
   std::cout << "      fRasterYSize : " << fRasterYSize << std::endl;
   std::cout << "       fGausXSigma : " << fGausXSigma  << std::endl;
   std::cout << "       fGausYSigma : " << fGausYSigma  << std::endl;
   std::cout << "      fBeamXCenter : " << fBeamXCenter << std::endl;
   std::cout << "      fBeamYCenter : " << fBeamYCenter << std::endl;
   std::cout << "         fDelta_th : " << fDelta_th << std::endl;
   std::cout << "         fDelta_ph : " << fDelta_ph << std::endl;
   std::cout << "         fDelta_dp : " << fDelta_dp << std::endl;
   std::cout << "         fz0       : " << fz0       << std::endl;
   std::cout << "         fT_L      : " << fT_L      << std::endl;
}
//______________________________________________________________________________
void SAMCEventGenerator::Init(){
   SAMCManager * man = SAMCManager::Instance();
   fBeamProfileType  = man->Which_Beam_Profile;
   fRasterXSize      = man->raster_x_size;
   fRasterYSize      = man->raster_y_size;
   fGausXSigma       = man->gaus_x_sigma;
   fGausYSigma       = man->gaus_y_sigma;
   fBeamXCenter      = man->beam_x_center;
   fBeamYCenter      = man->beam_y_center;
   fDelta_th         = man->delta_th;
   fDelta_ph         = man->delta_ph;
   fDelta_dp         = man->delta_dp;
   fz0               = man->z0;
   fT_L              = man->T_L;
}
//______________________________________________________________________________
Int_t SAMCEventGenerator::GenerateEvent(SAMCEvent& event){

   if ( fBeamProfileType == 0 ) {
      event.beam_x = fBeamXCenter + (gRandom->Rndm()-0.5)*fRasterXSize;
      event.beam_y = fBeamYCenter + (gRandom->Rndm()-0.5)*fRasterYSize;
   }
   else if ( fBeamProfileType == 1 ) {
      event.beam_x = gRandom->Gaus(fBeamXCenter,fGausXSigma);
      event.beam_y = gRandom->Gaus(fBeamYCenter,fGausYSigma);
   }
   else if ( fBeamProfileType == 2 ) {
      event.beam_x = fBeamXCenter + (gRandom->Rndm()-0.5)*fRasterXSize;
      event.beam_y = fBeamYCenter + (gRandom->Rndm()-0.5)*fRasterYSize;
      event.beam_x = gRandom->Gaus(event.beam_x,fGausXSigma);
      event.beam_y = gRandom->Gaus(event.beam_y,fGausYSigma);
   }
   event.reactz_gen = fz0 + (gRandom->Rndm()-0.5)*fT_L;
   //cout<<"---DEBUG---"<<endl;
   //cout<<"z0="<<z0<<endl;
   //cout<<"T_L="<<T_L<<endl;
   //cout<<"beam_x="<<Event->beam_x<<endl;
   //cout<<"reactz_gen="<<Event->reactz_gen<<endl;
   //cout<<"T_theta="<<Event->T_theta<<endl;
   //cout<<"-----------"<<endl;
   event.th_tg_gen = (gRandom->Rndm()-0.5)*fDelta_th;
   event.ph_tg_gen = (gRandom->Rndm()-0.5)*fDelta_ph;
   event.dp_gen    = (gRandom->Rndm()-0.5)*fDelta_dp;
   return 0;
}
//______________________________________________________________________________
