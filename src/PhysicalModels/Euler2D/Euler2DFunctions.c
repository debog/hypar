#include <math.h>
#include <basic.h>
#include <physicalmodels/euler2d.h>

int Euler2DRoeAverage(double *a_uavg,double *a_uL,double *a_uR,void *a_p)
{
  Euler2D *param  = (Euler2D*) a_p;
  _Euler2DRoeAverage_(a_uavg,a_uL,a_uR,param);
  return(0);
}
