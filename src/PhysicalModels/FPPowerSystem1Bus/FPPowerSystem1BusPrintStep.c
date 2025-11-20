#include <stdio.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <hypar.h>

int FPPowerSystem1BusPrintStep(void* a_s,void *a_m,double a_t)
{
  HyPar             *solver = (HyPar*)              a_s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*)  solver->m_physics;
  printf("Domain integral of the probability density function: %1.16E\n",params->m_pdf_integral);
  return(0);
}
