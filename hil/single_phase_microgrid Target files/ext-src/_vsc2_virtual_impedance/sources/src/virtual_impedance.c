#include "virtual_impedance.h"

// Virtual Impedance functions implementation
void initVirtualImpedance(VirtualImpedance *vz, float sampling_time, float resistance, float inductance, float zeta, float pole)
{
  // general
  vz->sampling_time = sampling_time;
  vz->resistance = resistance;
  vz->inductance = inductance;
  vz->zeta = zeta;
  vz->omega_pole = pole;

  // recurrence equation gains
  calculateVZcoefficients(vz);

  // previous interest variables
  vz->current_kminus1 = 0;
  vz->current_kminus2 = 0;
  vz->voltage_kminus1 = 0;
  vz->voltage_kminus2 = 0;
}

void executeVirtualImpedance(VirtualImpedance *vz, float current, float *voltage_drop)
{
  // calculate voltage drop
  float v_drop = vz->alpha0 * current + vz->alpha1 * vz->current_kminus1 + vz->alpha2 * vz->current_kminus2 - vz->beta1 * vz->voltage_kminus1 - vz->beta2 * vz->voltage_kminus2;
  // store previous variables
  vz->current_kminus2 = vz->current_kminus1;
  vz->current_kminus1 = current;
  vz->voltage_kminus2 = vz->voltage_kminus1;
  vz->voltage_kminus1 = v_drop;
  // output voltage drop
  *voltage_drop = v_drop;
}

void vz_setSamplingTime(VirtualImpedance *vz, float sampling_time)
{
  vz->sampling_time = sampling_time;
  calculateVZcoefficients(vz);
}

void setResistance(VirtualImpedance *vz, float resistance)
{
  vz->resistance = resistance;
  calculateVZcoefficients(vz);
}

void setInductance(VirtualImpedance *vz, float inductance)
{
  vz->inductance = inductance;
  calculateVZcoefficients(vz);
}

void setDamping(VirtualImpedance *vz, float zeta)
{
  vz->zeta = zeta;
  calculateVZcoefficients(vz);
}

void setComplexPole(VirtualImpedance *vz, float pole)
{
  vz->omega_pole = pole;
  calculateVZcoefficients(vz);
}

void calculateVZcoefficients(VirtualImpedance *vz)
{

  float omega_pole_2 = vz->omega_pole*vz->omega_pole;
  float sampling_time_2 = vz->sampling_time*vz->sampling_time;
  float inv_divider = 1.0 / (omega_pole_2 * sampling_time_2 + 4 * vz->zeta * vz->omega_pole * vz->sampling_time + 4);
  vz->alpha0 = (sampling_time_2 * omega_pole_2 * vz->resistance + 2 * vz->inductance * vz->sampling_time * omega_pole_2) * inv_divider;
  vz->alpha1 = (2 * sampling_time_2 * omega_pole_2 * vz->resistance) * inv_divider;
  vz->alpha2 = (sampling_time_2 * omega_pole_2 * vz->resistance - 2 * vz->inductance * vz->sampling_time * omega_pole_2) * inv_divider;
  vz->beta1 = (2 * sampling_time_2 * omega_pole_2 - 8) * inv_divider;
  vz->beta2 = (sampling_time_2 * omega_pole_2 - 4 * vz->zeta * vz->sampling_time * vz->omega_pole + 4) * inv_divider;
}