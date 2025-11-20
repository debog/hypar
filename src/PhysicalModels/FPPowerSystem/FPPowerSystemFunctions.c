#include <math.h>
#include <physicalmodels/fppowersystem.h>

double FPPowerSystemDriftFunction(int a_dir,void *a_p,double a_x,double a_y, double a_t)
{
  FPPowerSystem *params = (FPPowerSystem*) a_p;

  double drift = 0;
  if (a_dir == 0)   drift = params->m_O_s * (a_y - 1.0);
  else if (a_dir == 1) {
    if      (a_t < params->m_tf)    params->Pmax = params->m_E*params->m_V/params->m_g1;
    else if (a_t < params->tcl)   params->Pmax = 0.0;
    else                        params->Pmax = params->m_E*params->m_V/params->m_g2;
    drift = (1.0/(2.0*params->H)) * (params->Pm - params->Pmax*sin(a_x) - params->D*(a_y-1.0));
  }

  return drift;
}

double FPPowerSystemDissipationFunction(int a_dir,void *a_p,double a_t)
{
  FPPowerSystem *params = (FPPowerSystem*) a_p;

  double dissp = 0;
  if (a_dir == 1) {
    double gamma = params->D / (2.0*params->H);
    dissp =   (1.0/(2.0*params->H)) * (1.0/(2.0*params->H))
            * ((params->m_l*params->m_q*params->m_q)/(params->m_l*gamma+1.0))
            * (1.0 - exp(-(gamma+1.0/params->m_l)*a_t));
  }

  return(dissp);
}
