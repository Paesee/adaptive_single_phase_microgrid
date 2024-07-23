// include header file
#include "mrac.h"

// include math.h (it must be included here for HIL compatibility, can't be called in the library header file)
#include <math.h>

/*
 * MRAC Controller Main Functions Implementation
 */

/// @brief initialize MRAC controller
/// @param mrac 
/// @param sampling_time 
/// @param gamma 
/// @param settling_time 
/// @param gain_correction 
/// @param phase_shift 
void initMRAC(MRAC *mrac, float sampling_time)
{
  // General initialization
  mrac->sampling_time = sampling_time;
  mrac->gamma = 1.0;
  mrac->ts_times_gamma = sampling_time * mrac->gamma;

  // Reference model initialization
  mrac->pole_0 = 0.0;
  mrac->pole_1 = 0.0;
  mrac->pole_2 = 0.0;

  mrac->rho = 1.0;
  mrac->phi = 0.0;

  // Set alpha and beta to initial values (these should be set properly according to the system requirements)
  mrac->alpha0 = 0.0;
  mrac->alpha1 = 0.0;
  mrac->alpha2 = 0.0;
  mrac->alpha3 = 0.0;
  mrac->beta1 = 0.0;
  mrac->beta2 = 0.0;
  mrac->beta3 = 0.0;

  // Augmented error
  mrac->a_err_kminus1 = 0.0;

  // Recurrence equation variables for y_{m}
  mrac->ym_kminus1 = 0.0;
  mrac->ym_kminus2 = 0.0;
  mrac->ym_kminus3 = 0.0;

  // Recurrence equation variables for v
  mrac->v_kminus1 = 0.0;
  mrac->v_kminus2 = 0.0;
  mrac->v_kminus3 = 0.0;
  mrac->ub_kminus1 = 0.0;
  mrac->ub_kminus2 = 0.0;
  mrac->ub_kminus3 = 0.0;
  mrac->ub_kminus4 = 0.0;

  // Output boundary
  mrac->u_boundary = 99999.0; // Default value, should be set properly

  // Initialize previous states to zero
  mrac->x1_kminus1 = 0.0;
  mrac->x2_kminus1 = 0.0;
  mrac->x3_kminus1 = 0.0;
  mrac->r_kminus1 = 0.0;
  mrac->ds_kminus1 = 0.0;
  mrac->dc_kminus1 = 0.0;

  mrac->x1_kminus2 = 0.0;
  mrac->x2_kminus2 = 0.0;
  mrac->x3_kminus2 = 0.0;
  mrac->r_kminus2 = 0.0;
  mrac->ds_kminus2 = 0.0;
  mrac->dc_kminus2 = 0.0;

  mrac->x1_kminus3 = 0.0;
  mrac->x2_kminus3 = 0.0;
  mrac->x3_kminus3 = 0.0;
  mrac->r_kminus3 = 0.0;
  mrac->ds_kminus3 = 0.0;
  mrac->dc_kminus3 = 0.0;

  // Initialize zeta values to zero
  mrac->zeta_x1 = 0.0;
  mrac->zeta_x2 = 0.0;
  mrac->zeta_x3 = 0.0;
  mrac->zeta_r = 0.0;
  mrac->zeta_ds = 0.0;
  mrac->zeta_dc = 0.0;

  mrac->zeta_x1_kminus1 = 0.0;
  mrac->zeta_x2_kminus1 = 0.0;
  mrac->zeta_x3_kminus1 = 0.0;
  mrac->zeta_r_kminus1 = 0.0;
  mrac->zeta_ds_kminus1 = 0.0;
  mrac->zeta_dc_kminus1 = 0.0;

  mrac->zeta_x1_kminus2 = 0.0;
  mrac->zeta_x2_kminus2 = 0.0;
  mrac->zeta_x3_kminus2 = 0.0;
  mrac->zeta_r_kminus2 = 0.0;
  mrac->zeta_ds_kminus2 = 0.0;
  mrac->zeta_dc_kminus2 = 0.0;

  mrac->zeta_x1_kminus3 = 0.0;
  mrac->zeta_x2_kminus3 = 0.0;
  mrac->zeta_x3_kminus3 = 0.0;
  mrac->zeta_r_kminus3 = 0.0;
  mrac->zeta_ds_kminus3 = 0.0;
  mrac->zeta_dc_kminus3 = 0.0;

  // Initialize theta values to zero
  mrac->theta_x1 = 0.0;
  mrac->theta_x2 = 0.0;
  mrac->theta_x3 = 0.0;
  mrac->theta_u = 0.0;
  mrac->theta_r = 0.0;
  mrac->theta_ds = 0.0;
  mrac->theta_dc = 0.0;

  mrac->theta_x1_kminus1 = 0.0;
  mrac->theta_x2_kminus1 = 0.0;
  mrac->theta_x3_kminus1 = 0.0;
  mrac->theta_r_kminus1 = 0.0;
  mrac->theta_ds_kminus1 = 0.0;
  mrac->theta_dc_kminus1 = 0.0;

  // Normalization signal
  mrac->inv_m2_kminus1 = 1.0; // Start with 1 to avoid division by zero
}

/// @brief compute MRAC control law
/// @param mrac 
/// @param r_al 
/// @param r_be 
/// @param x1 
/// @param x2 
/// @param x3 
/// @param ds 
/// @param dc 
/// @param u_ctrl 
/// @param r_c 
/// @param y_m 
void computeMRAC(MRAC *mrac, float r_al, float r_be, float x1, float x2, float x3, float ds, float dc, float *u_ctrl, float *y_m)
{
  // correct reference magnitude and phase shift
  float rc = computeReferenceCorrection(r_al, r_be, mrac->rho, mrac->phi);

  // compute reference model
  float ym = computeWm(mrac, rc, mrac->r_kminus1, mrac->r_kminus2, mrac->r_kminus3, mrac->ym_kminus1, mrac->ym_kminus2, mrac->ym_kminus3);

  // compute zeta
  computeZeta(mrac, rc, x1, x2, x3, ds, dc);

  // compute theta
  computeTheta(mrac);

  // compute control action
  float u = computeU(mrac, rc, x1, x2, x3, ds, dc);

  // compute output boundary
  float u_bounded = computeBoundary(mrac, u);

  // control action filtering
  float v = computeWm(mrac, u_bounded, mrac->ub_kminus1, mrac->ub_kminus2, mrac->ub_kminus3, mrac->v_kminus1, mrac->v_kminus2, mrac->v_kminus3);

  // tracking error
  float t_err = x2 - ym; // Y = [0 1 0] * [I_{Lf} V_{C} I_{0}]'

  // augmented error
  computeAugmentedError(mrac, t_err, v);

  // update normalization signal;
  computeM2(mrac);

  // update mrac
  updateMRAC(mrac, ym, v, u_bounded, x1, x2, x3, rc, ds, dc);

  // output control action
  *u_ctrl = u_bounded;
  *y_m = ym;
}

/*
 * MRAC Controller Auxiliary Functions Implementation
 */

/// @brief compute gain and phase correction for reference
/// @param ral 
/// @param rbe 
/// @param gain 
/// @param phase 
/// @return 
extern inline float computeReferenceCorrection(float ral, float rbe, float gain, float phase)
{
  return (gain * (cos(phase) * ral + sin(phase) * rbe));
}

/// @brief compute W_{m} reccurence equations gains according to poles selection
/// @param mrac 
extern inline void computeWmCoeffs(MRAC *mrac)
{
  //TODO
}

/// @brief compute zeta
/// @param mrac 
/// @param rc 
/// @param x1 
/// @param x2 
/// @param x3 
/// @param ds 
/// @param dc 
extern inline void computeZeta(MRAC *mrac, float rc, float x1, float x2, float x3, float ds, float dc)
{
  // zeta = [x1, x2, x3, r, ds, dc]
  mrac->zeta_x1 = mrac->alpha0 * x1               + mrac->alpha1 * mrac->x1_kminus1 + mrac->alpha2 * mrac->x1_kminus2 + mrac->alpha3 * mrac->x1_kminus3 - mrac->beta1 * mrac->zeta_x1_kminus1 - mrac->beta2 * mrac->zeta_x1_kminus2 - mrac->beta3 * mrac->zeta_x1_kminus3;
  mrac->zeta_x2 = mrac->alpha0 * x2               + mrac->alpha1 * mrac->x2_kminus1 + mrac->alpha2 * mrac->x2_kminus2 + mrac->alpha3 * mrac->x2_kminus3 - mrac->beta1 * mrac->zeta_x2_kminus1 - mrac->beta2 * mrac->zeta_x2_kminus2 - mrac->beta3 * mrac->zeta_x2_kminus3;
  mrac->zeta_x3 = mrac->alpha0 * x3               + mrac->alpha1 * mrac->x3_kminus1 + mrac->alpha2 * mrac->x3_kminus2 + mrac->alpha3 * mrac->x3_kminus3 - mrac->beta1 * mrac->zeta_x3_kminus1 - mrac->beta2 * mrac->zeta_x3_kminus2 - mrac->beta3 * mrac->zeta_x3_kminus3;
  mrac->zeta_u  = mrac->alpha0 * mrac->ub_kminus1 + mrac->alpha1 * mrac->ub_kminus2 + mrac->alpha2 * mrac->ub_kminus3 + mrac->alpha3 * mrac->ub_kminus4 - mrac->beta1 * mrac->zeta_u_kminus1  - mrac->beta2 * mrac->zeta_u_kminus2  - mrac->beta3 * mrac->zeta_u_kminus3;
  mrac->zeta_r  = mrac->alpha0 * rc               + mrac->alpha1 * mrac->r_kminus1  + mrac->alpha2 * mrac->r_kminus2  + mrac->alpha3 * mrac->r_kminus3  - mrac->beta1 * mrac->zeta_r_kminus1  - mrac->beta2 * mrac->zeta_r_kminus2  - mrac->beta3 * mrac->zeta_r_kminus3;
  mrac->zeta_ds = mrac->alpha0 * ds               + mrac->alpha1 * mrac->ds_kminus1 + mrac->alpha2 * mrac->ds_kminus2 + mrac->alpha3 * mrac->ds_kminus3 - mrac->beta1 * mrac->zeta_ds_kminus1 - mrac->beta2 * mrac->zeta_ds_kminus2 - mrac->beta3 * mrac->zeta_ds_kminus3;
  mrac->zeta_dc = mrac->alpha0 * dc               + mrac->alpha1 * mrac->dc_kminus1 + mrac->alpha2 * mrac->dc_kminus2 + mrac->alpha3 * mrac->dc_kminus3 - mrac->beta1 * mrac->zeta_dc_kminus1 - mrac->beta2 * mrac->zeta_dc_kminus2 - mrac->beta3 * mrac->zeta_dc_kminus3;
}

/// @brief compute theta
/// @param mrac 
extern inline void computeTheta(MRAC *mrac)
{
  // zeta vector common factor only once
  float common_factor = mrac->ts_times_gamma * mrac->a_err_kminus1 * mrac->inv_m2_kminus1;
  // compute theta
  mrac->theta_x1 = mrac->theta_x1_kminus1 - common_factor * mrac->zeta_x1_kminus1; // theta_{x1} => X[1] = I_{Lf}
  mrac->theta_x2 = mrac->theta_x2_kminus1 - common_factor * mrac->zeta_x2_kminus1; // theta_{x2} => X[2] = V_{C}
  mrac->theta_x3 = mrac->theta_x3_kminus1 - common_factor * mrac->zeta_x3_kminus1; // theta_{x3} => X[3] = I_{0}
  mrac->theta_u  = mrac->theta_u_kminus1  - common_factor * mrac->zeta_u_kminus1; // theta_{x3} => X[3] = U
  mrac->theta_r  = mrac->theta_r_kminus1  - common_factor * mrac->zeta_r_kminus1;  // theta_{r} => R_{C} = corrected reference
  mrac->theta_ds = mrac->theta_ds_kminus1 - common_factor * mrac->zeta_ds_kminus1; // theta_{ds} => D_{sin} = sin disturbance rejection
  mrac->theta_dc = mrac->theta_dc_kminus1 - common_factor * mrac->zeta_dc_kminus1; // theta_{dc} => D_{cos} = cos disturbance rejection
}

/// @brief compute control action u
/// @param mrac 
/// @param x1 
/// @param x2 
/// @param x3 
/// @param rc 
/// @param ds 
/// @param dc 
/// @return 
extern inline float computeU(MRAC *mrac, float rc, float x1, float x2, float x3, float ds, float dc)
{
  return (mrac->theta_x1 * x1 + mrac->theta_x2 * x2 + mrac->theta_x3 * x3 + mrac->theta_u * mrac->ub_kminus1 + mrac->theta_r * rc + mrac->theta_ds * ds + mrac->theta_dc * dc);
}

/// @brief compute output boundary
/// @param mrac 
/// @param u 
/// @return 
extern inline float computeBoundary(MRAC *mrac, float u)
{
  if(u != 0)
  {
    float u_bounded = u;
    float uu = u * u;
    float inv_sqrt_uu = 1.0 / sqrt(uu);
    if ((uu * inv_sqrt_uu) >= mrac->u_boundary)
      u_bounded = (mrac->u_boundary * u) * inv_sqrt_uu;
    return u_bounded;
  }else{
    return 0;
  }
}

/// @brief compute augmented error
/// @param mrac
/// @param t_err 
/// @param v 
/// @return 
extern inline void computeAugmentedError(MRAC *mrac, float t_err, float v)
{
  mrac->csi = - v + mrac->theta_x1 * mrac->zeta_x1 + mrac->theta_x2 * mrac->zeta_x2 + mrac->theta_x3 * mrac->zeta_x3 + mrac->theta_u * mrac->zeta_u + mrac->theta_r * mrac->zeta_r + mrac->theta_ds * mrac->zeta_ds + mrac->theta_dc * mrac->zeta_dc;
  mrac->a_err_kminus1 = t_err + mrac->csi;
}

/// @brief compute normalization signal
/// @param mrac 
/// @return 
extern inline void computeM2(MRAC *mrac)
{
  mrac->inv_m2_kminus1 = 1.0 / (1 + mrac->zeta_x1 * mrac->zeta_x1 + mrac->zeta_x2 * mrac->zeta_x2 + mrac->zeta_x3 * mrac->zeta_x3 + mrac->zeta_u * mrac->zeta_u + mrac->zeta_r * mrac->zeta_r + mrac->zeta_ds * mrac->zeta_ds + mrac->zeta_dc * mrac->zeta_dc + mrac->csi * mrac->csi);
}

extern inline float computeWm(MRAC *mrac, float u, float u_kminus1, float u_kminus2, float u_kminus3, float y_kminus1, float y_kminus2, float y_kminus3)
{
  return (mrac->alpha0 * u + mrac->alpha1 * u_kminus1 + mrac->alpha2 * u_kminus2 + mrac->alpha3 * u_kminus3 - mrac->beta1 * y_kminus1 - mrac->beta2 * y_kminus2 - mrac->beta3 * y_kminus3);
}

extern inline void updateMRAC(MRAC *mrac, float ym, float v, float u_bounded, float x1, float x2, float x3, float rc, float ds, float dc)
{
  // Recurrence equation variables for y_m
  mrac->ym_kminus3 = mrac->ym_kminus2;
  mrac->ym_kminus2 = mrac->ym_kminus1;
  mrac->ym_kminus1 = ym;

  // Recurrence equation variables for v
  mrac->v_kminus3 = mrac->v_kminus2;
  mrac->v_kminus2 = mrac->v_kminus1;
  mrac->v_kminus1 = v;
  mrac->ub_kminus4 = mrac->ub_kminus3;
  mrac->ub_kminus3 = mrac->ub_kminus2;
  mrac->ub_kminus2 = mrac->ub_kminus1;
  mrac->ub_kminus1 = u_bounded;

  // X
  mrac->x1_kminus3 = mrac->x1_kminus2;
  mrac->x1_kminus2 = mrac->x1_kminus1;
  mrac->x1_kminus1 = x1;
  mrac->x2_kminus3 = mrac->x2_kminus2;
  mrac->x2_kminus2 = mrac->x2_kminus1;
  mrac->x2_kminus1 = x2;
  mrac->x3_kminus3 = mrac->x3_kminus2;
  mrac->x3_kminus2 = mrac->x3_kminus1;
  mrac->x3_kminus1 = x3;
  mrac->r_kminus3 = mrac->r_kminus2;
  mrac->r_kminus2 = mrac->r_kminus1;
  mrac->r_kminus1 = rc;
  mrac->ds_kminus3 = mrac->ds_kminus2;
  mrac->ds_kminus2 = mrac->ds_kminus1;
  mrac->ds_kminus1 = ds;
  mrac->dc_kminus3 = mrac->dc_kminus2;
  mrac->dc_kminus2 = mrac->dc_kminus1;
  mrac->dc_kminus1 = dc;

  // Zeta values
  mrac->zeta_x1_kminus3 = mrac->zeta_x1_kminus2;
  mrac->zeta_x1_kminus2 = mrac->zeta_x1_kminus1;
  mrac->zeta_x1_kminus1 = mrac->zeta_x1;
  mrac->zeta_x2_kminus3 = mrac->zeta_x2_kminus2;
  mrac->zeta_x2_kminus2 = mrac->zeta_x2_kminus1;
  mrac->zeta_x2_kminus1 = mrac->zeta_x2;
  mrac->zeta_x3_kminus3 = mrac->zeta_x3_kminus2;
  mrac->zeta_x3_kminus2 = mrac->zeta_x3_kminus1;
  mrac->zeta_x3_kminus1 = mrac->zeta_x3;
  mrac->zeta_u_kminus3 = mrac->zeta_u_kminus2;
  mrac->zeta_u_kminus2 = mrac->zeta_u_kminus1;
  mrac->zeta_u_kminus1 = mrac->zeta_u;
  mrac->zeta_r_kminus3 = mrac->zeta_r_kminus2;
  mrac->zeta_r_kminus2 = mrac->zeta_r_kminus1;
  mrac->zeta_r_kminus1 = mrac->zeta_r;
  mrac->zeta_ds_kminus3 = mrac->zeta_ds_kminus2;
  mrac->zeta_ds_kminus2 = mrac->zeta_ds_kminus1;
  mrac->zeta_ds_kminus1 = mrac->zeta_ds;
  mrac->zeta_dc_kminus3 = mrac->zeta_dc_kminus2;
  mrac->zeta_dc_kminus2 = mrac->zeta_dc_kminus1;
  mrac->zeta_dc_kminus1 = mrac->zeta_dc;

  // Theta values
  mrac->theta_x1_kminus1 = mrac->theta_x1;
  mrac->theta_x2_kminus1 = mrac->theta_x2;
  mrac->theta_x3_kminus1 = mrac->theta_x3;
  mrac->theta_u_kminus1 = mrac->theta_u;
  mrac->theta_r_kminus1 = mrac->theta_r;
  mrac->theta_ds_kminus1 = mrac->theta_ds;
  mrac->theta_dc_kminus1 = mrac->theta_dc;
}

/*
 * Setters
 */

/// @brief set theta for better initilization
/// @param mrac 
/// @param theta_x1 
/// @param theta_x2 
/// @param theta_x3 
/// @param theta_r 
/// @param theta_ds 
/// @param theta_dc 
void setTheta(MRAC *mrac, float theta_x1, float theta_x2, float theta_x3, float theta_u, float theta_r, float theta_ds, float theta_dc)
{
  mrac->theta_x1 = theta_x1;
  mrac->theta_x2 = theta_x2;
  mrac->theta_x3 = theta_x3;
  mrac->theta_u  = theta_u;
  mrac->theta_r  = theta_r;
  mrac->theta_ds = theta_ds;
  mrac->theta_dc = theta_dc;
}

/// @brief set output boundary
/// @param mrac 
/// @param boundary 
void setBoundary(MRAC *mrac, float boundary)
{
  mrac->u_boundary = boundary;
}

/// @brief set reference model W_{m} gain correction factor rho
/// @param mrac 
/// @param rho 
void setGain(MRAC *mrac, float rho)
{
  mrac->rho = rho;
}

/// @brief set reference model W_{m} phase correction factor phi
/// @param mrac 
/// @param phi 
void setPhaseShift(MRAC *mrac, float phi)
{
  mrac->phi = phi;
}

/// @brief set one of the poles of the transfer function W_{m}
/// @param mrac 
/// @param pole_0 
void setPole0(MRAC *mrac, float pole_0)
{
  mrac->pole_0 = pole_0;
  computeWmCoeffs(mrac);
}

/// @brief set one of the poles of the transfer function W_{m}
/// @param mrac 
/// @param pole_1 
void setPole1(MRAC *mrac, float pole_1)
{
  mrac->pole_1 = pole_1;
  computeWmCoeffs(mrac);
}

/// @brief set one of the poles of the transfer function W_{m}
/// @param mrac 
/// @param pole_2 
void setPole2(MRAC *mrac, float pole_2)
{
  mrac->pole_2 = pole_2;
  computeWmCoeffs(mrac);
}

/// @brief set adaptive gain
/// @param mrac 
/// @param gamma 
void setGamma(MRAC *mrac, float gamma)
{
  mrac->gamma = gamma;
  mrac->ts_times_gamma = mrac->sampling_time * gamma;
}

/// @brief set transfer function W_{m} directly (may be used instead of setting poles position)
/// @param mrac 
/// @param alpha0 
/// @param alpha1 
/// @param alpha2 
/// @param alpha3 
/// @param beta0 
/// @param beta1 
/// @param beta2 
void setWmCoeffs(MRAC *mrac, float alpha0, float alpha1, float alpha2, float alpha3, float beta1, float beta2, float beta3)
{
  mrac->alpha0 = alpha0;
  mrac->alpha1 = alpha1;
  mrac->alpha2 = alpha2;
  mrac->alpha3 = alpha3;
  mrac->beta1 = beta1;
  mrac->beta2 = beta2;
  mrac->beta3 = beta3;
}

/*
 * Getters
 */

/// @brief get theta vector
/// @param mrac 
/// @param theta_x1 
/// @param theta_y 
/// @param theta_r 
/// @param theta_sin 
/// @param theta_cos 
void getTheta(MRAC *mrac, float *theta_x1, float *theta_x2, float *theta_x3, float *theta_u, float *theta_r, float *theta_ds, float *theta_dc)
{
  *theta_x1 = mrac->theta_x1;
  *theta_x2 = mrac->theta_x2;
  *theta_x3 = mrac->theta_x3;
  *theta_u = mrac->theta_u;
  *theta_r = mrac->theta_r;
  *theta_ds = mrac->theta_ds;
  *theta_dc = mrac->theta_dc;
}