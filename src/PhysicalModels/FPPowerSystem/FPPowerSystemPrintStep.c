#include <stdio.h>
#include <physicalmodels/fppowersystem.h>
#include <hypar.h>

int FPPowerSystemPrintStep(void* s,void *m,double t)
{
  HyPar          *solver = (HyPar*)         s;
  FPPowerSystem  *params = (FPPowerSystem*) solver->m_physics;
  printf("Domain integral of the probability density function: %1.16E\n",params->m_pdf_integral);
  return(0);
}
