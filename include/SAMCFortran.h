#ifndef SAMCFortran_HH
#define SAMCFortran_HH


extern "C" {

   //must be like the following
   //change l(r)txfit_ in fortran to real function ...
   //otherwise no result is passed.

   //for left arm
   float ltxfit_(float*,int&);
   float ldelta_(float*,int&);
   float ltheta_(float*,int&);
   float lphi_(float*,int&);
   float ly00_(float*,int&);
   float x_e_q1ex_(float*,int&);
   float y_e_q1ex_(float*,int&);
   float x_e_dent_(float*,int&);
   float y_e_dent_(float*,int&);
   float x_e_dext_(float*,int&);
   float y_e_dext_(float*,int&);
   float x_e_q3en_(float*,int&);
   float y_e_q3en_(float*,int&);
   float x_e_q3ex_(float*,int&);
   float y_e_q3ex_(float*,int&);
   float x_e_fp_(float*,int&);
   float y_e_fp_(float*,int&);
   float p_e_fp_(float*,int&);
   float t_e_fp_(float*,int&);
   void left_init_r_function_();
   void left_getindex_(float*,float*,int*);
   float left_rfunction_(float*,float*,float*,float*);

   //for right arm
   float rtxfit_(float*,int&);
   float rdelta_(float*,int&);
   float rtheta_(float*,int&);
   float rphi_(float*,int&);
   float ry00_(float*,int&);
   float x_h_q1ex_(float*,int&);
   float y_h_q1ex_(float*,int&);
   float x_h_dent_(float*,int&);
   float y_h_dent_(float*,int&);
   float x_h_dext_(float*,int&);
   float y_h_dext_(float*,int&);
   float x_h_q3en_(float*,int&);
   float y_h_q3en_(float*,int&);
   float x_h_q3ex_(float*,int&);
   float y_h_q3ex_(float*,int&);
   float x_h_fp_(float*,int&);
   float y_h_fp_(float*,int&);
   float p_h_fp_(float*,int&);
   float t_h_fp_(float*,int&);
   void right_init_r_function_();
   void right_getindex_(float*,float*,int*);
   float right_rfunction_(float*,float*,float*,float*);
}

#endif

