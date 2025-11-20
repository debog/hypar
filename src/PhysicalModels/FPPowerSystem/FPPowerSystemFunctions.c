#include <math.h>
#include <physicalmodels/fppowersystem.h>

double FPPowerSystemDriftFunction(int dir,void *p,double x,double y, double t)
{
  FPPowerSystem *params = (FPPowerSystem*) p;

  double drift = 0;
  if (dir == 0)   drift = params->m_O_s * (y - 1.0);
  else if (dir == 1) {
    if      (t < params->m_tf)    params->Pmax = params->m_E*params->m_V/params->m_g1;
    else if (t < params->tcl)   params->Pmax = 0.0;
    else                        params->Pmax = params->m_E*params->m_V/params->m_g2;
    drift = (1.0/(2.0*params->H)) * (params->Pm - params->Pmax*sin(x) - params->D*(y-1.0));
  }

  return drift;
}

double FPPowerSystemDissipationFunction(int dir,void *p,double t)
{
  FPPowerSystem *params = (FPPowerSystem*) p;

  double dissp = 0;
  if (dir == 1) {
    double gamma = params->D / (2.0*params->H);
    dissp =   (1.0/(2.0*params->H)) * (1.0/(2.0*params->H))
            * ((params->m_l*params->m_q*params->m_q)/(params->m_l*gamma+1.0))
            * (1.0 - exp(-(gamma+1.0/params->m_l)*t));
  }

  return(dissp);
}
