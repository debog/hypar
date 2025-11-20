#include <math.h>
#include <physicalmodels/fppowersystem1bus.h>

double FPPowerSystem1BusDriftFunction(int a_dir,void *a_p,double a_x,double a_y, double a_t)
{
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) a_p;

  double drift = 0;
  if (a_dir == _XDIR_) {
    drift = params->omegaB * (a_y - params->omegaS);
  } else if (a_dir == _YDIR_) {
    drift = (params->omegaS/(2*params->H))
          * (params->Pm_avg - params->Pmax*sin(a_x) - params->D*(a_y-params->omegaS));
  }

  return drift;
}

double FPPowerSystem1BusDissipationFunction(int a_dir1,int a_dir2,void *a_p,double a_t)
{
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) a_p;

  double sigma  = params->sigma;
  double omegaS = params->omegaS;
  double omegaB = params->omegaB;
  double lambda = params->lambda;
  double H      = params->H;
  double D      = params->D;

  double dissp = 0;
  if (a_dir1 == _YDIR_) {

    double term = (sigma*sigma*omegaS*omegaS) / (4.0*H*H);
    double expterm = exp(-a_t/lambda);

    if (a_dir2 == _XDIR_) {
      dissp = term * (lambda*omegaB) * (lambda*(1-expterm) - a_t*expterm);
      /* dissp = term * lambda*omegaB*lambda; */
    } else if (a_dir2 == _YDIR_) {
      double gamma = D*omegaS / (2.0*H);
      dissp = term * (lambda*(1-expterm) + (gamma*lambda*(a_t*expterm-lambda*(1-expterm))) );
      /* dissp = term * params->lambda * (1 - params->lambda*gamma); */
    }

  }

  return(dissp);
}
