#include <physicalmodels/euler1d.h>

#define _Euler1DGetFlowVar_GPU(u,rho,v,e,P,gamma) \
  { \
    rho = u[0]; \
    v   = u[1] / rho; \
    e   = u[2]; \
    P   = (e - 0.5*rho*v*v) * (gamma-1.0); \
  }

#define _gpuEuler1DGetFlowVar2_(u,np,rho,v,e,P,gamma) \
  { \
    rho = u[0]; \
    v   = u[np] / rho; \
    e   = u[2*np]; \
    P   = (e - 0.5*rho*v*v) * (gamma-1.0); \
  }


