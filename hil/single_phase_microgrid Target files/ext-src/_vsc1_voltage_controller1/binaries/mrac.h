#ifndef MRAC_H
#define MRAC_H

// MRAC Controller Struct definition
typedef struct MRAC
{
  // general
  float sampling_time;
  float gamma;
  float ts_times_gamma;

  // reference model W_{m}(z)
  float pole_0;
  float pole_1;
  float pole_2;
  // phase and gain correction for W_{m}(z)
  float rho;
  float phi;
  // recurrence equation gains
  float alpha0;
  float alpha1;
  float alpha2;
  float alpha3;
  float beta1;
  float beta2;
  float beta3;

  // augmented error
  float a_err_kminus1;

  // recurrence equation variables for y_{m} ==> Y(z)/R(z) = W_{m}(z)
  float ym_kminus1;
  float ym_kminus2;
  float ym_kminus3;

  // recurrence equation variables for v ==> V(z)/U_{bounded}(z) = W_{m}(z)
  float v_kminus1;
  float v_kminus2;
  float v_kminus3;
  float ub_kminus1;
  float ub_kminus2;
  float ub_kminus3;
  float ub_kminus4;

  // output boundary
  float u_boundary;

  // X[k-1]
  float x1_kminus1;
  float x2_kminus1;
  float x3_kminus1;
  float r_kminus1;
  float ds_kminus1;
  float dc_kminus1;
  // X[k-2]
  float x1_kminus2;
  float x2_kminus2;
  float x3_kminus2;
  float r_kminus2;
  float ds_kminus2;
  float dc_kminus2;
  // X[k-3]
  float x1_kminus3;
  float x2_kminus3;
  float x3_kminus3;
  float r_kminus3;
  float ds_kminus3;
  float dc_kminus3;

  // zeta[k]
  float zeta_x1;
  float zeta_x2;
  float zeta_x3;
  float zeta_u;
  float zeta_r;
  float zeta_ds;
  float zeta_dc;
  // zeta[k-1]
  float zeta_x1_kminus1;
  float zeta_x2_kminus1;
  float zeta_x3_kminus1;
  float zeta_u_kminus1;
  float zeta_r_kminus1;
  float zeta_ds_kminus1;
  float zeta_dc_kminus1;
  // zeta[k-2]
  float zeta_x1_kminus2;
  float zeta_x2_kminus2;
  float zeta_x3_kminus2;
  float zeta_u_kminus2;
  float zeta_r_kminus2;
  float zeta_ds_kminus2;
  float zeta_dc_kminus2;
  // zeta[k-3]
  float zeta_x1_kminus3;
  float zeta_x2_kminus3;
  float zeta_x3_kminus3;
  float zeta_u_kminus3;
  float zeta_r_kminus3;
  float zeta_ds_kminus3;
  float zeta_dc_kminus3;

  // theta[k]
  float theta_x1;
  float theta_x2;
  float theta_x3;
  float theta_u;
  float theta_r;
  float theta_ds;
  float theta_dc;
  // zeta[k-1]
  float theta_x1_kminus1;
  float theta_x2_kminus1;
  float theta_x3_kminus1;
  float theta_u_kminus1;
  float theta_r_kminus1;
  float theta_ds_kminus1;
  float theta_dc_kminus1;

  // normalization signal
  float inv_m2_kminus1;

  // teste
  float csi;
} MRAC;

void initMRAC(MRAC *mrac, float sampling_time);
void computeMRAC(MRAC *mrac, float r_al, float r_be, float x1, float x2, float x3, float ds, float dc, float *u_ctrl, float *r_c, float *y_m);

extern inline float computeReferenceCorrection(float ral, float rbe, float gain, float phase);
extern inline void computeWmCoeffs(MRAC *mrac);
extern inline void computeZeta(MRAC *mrac, float rc, float x1, float x2, float x3, float ds, float dc);
extern inline void computeTheta(MRAC *mrac);
extern inline float computeU(MRAC *mrac, float rc, float x1, float x2, float x3, float ds, float dc);
extern inline float computeBoundary(MRAC *mrac, float u);
extern inline void computeAugmentedError(MRAC *mrac, float t_err, float v);
extern inline void computeM2(MRAC *mrac);
extern inline float computeWm(MRAC *mrac, float u, float u_kminus1, float u_kminus2, float u_kminus3, float y_kminus1, float y_kminus2, float y_kminus3);
extern inline void updateMRAC(MRAC *mrac, float ym, float v, float u_bounded, float x1, float x2, float x3, float rc, float ds, float dc);

void setTheta(MRAC *mrac, float theta_x1, float theta_x2, float theta_x3, float theta_r, float theta_ds, float theta_dc);
void setBoundary(MRAC *mrac, float boundary);
void setGain(MRAC *mrac, float rho);
void setPhaseShift(MRAC *mrac, float phi);
void setPole0(MRAC *mrac, float pole_0);
void setPole1(MRAC *mrac, float pole_1);
void setPole2(MRAC *mrac, float pole_2);
void setGamma(MRAC *mrac, float gamma);
void setWmCoeffs(MRAC *mrac, float alpha0, float alpha1, float alpha2, float alpha3, float beta1, float beta2, float beta3);

void getTheta(MRAC *mrac, float *theta_x1, float *theta_x2, float *theta_x3, float *theta_r, float *theta_ds, float *theta_dc);

#endif