#include <math.h>
#include <mathfunctions.h>
#include <physicalmodels/numa2d.h>

int Numa2DCalculateStandardAtmosphere_1(void *a_p,double a_z,double *a_ExnerP,double *a_P,double *a_rho,double *a_T)
{
  Numa2D *physics = (Numa2D*) a_p;

  double R      = physics->m_R;
  double gamma  = physics->m_gamma;
  double g      = physics->m_g;

  /* reference quantities at zero altitude */
  double P0, T0;
  P0   = physics->Pref;
  T0   = physics->Tref;

  double inv_gamma_m1 = 1.0/(gamma-1.0);
  double Cp = gamma * inv_gamma_m1 * R;

  double theta  = T0;
  *a_ExnerP = 1.0 - (g/(Cp*theta))*a_z;
  *a_P      = P0 * raiseto((*a_ExnerP),gamma*inv_gamma_m1);
  *a_rho    = (P0/(R*theta)) * raiseto((*a_ExnerP),inv_gamma_m1);
  *a_T      = (*a_rho) * theta;

  return(0);
}

int Numa2DCalculateStandardAtmosphere_2(void *a_p,double a_z,double *a_ExnerP,double *a_P,double *a_rho,double *a_T)
{
  Numa2D *physics = (Numa2D*) a_p;

  double R      = physics->m_R;
  double gamma  = physics->m_gamma;
  double g      = physics->m_g;

  /* reference quantities at zero altitude */
  double P0, T0;
  P0 = physics->Pref;
  T0 = physics->Tref;

  double BV = 0.01; /* Brunt-Vaisala frequency */
  double inv_gamma_m1 = 1.0/(gamma-1.0);
  double Cp = gamma * inv_gamma_m1 * R;

  double term   = BV*BV*a_z/g;
  double theta  = T0 * exp(term);
  *a_ExnerP = 1.0 + (g*g/(Cp*T0*BV*BV)) * (exp(-term) - 1.0);
  *a_P      = P0 * raiseto((*a_ExnerP),gamma*inv_gamma_m1);
  *a_rho    = (P0/(R*theta)) * raiseto((*a_ExnerP),inv_gamma_m1);
  *a_T      = (*a_rho) * theta;

  return(0);
}
