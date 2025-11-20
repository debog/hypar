#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

/*

  Rohde, A., "Eigenvalues and Eigenvectors of the Euler Equations
  in General Geometries", AIAA 2001-2609

*/

int Euler2DLeftEigenvectors(double *a_u,double *a_L,void *a_p,int a_dir)
{
  Euler2D *param  = (Euler2D*)  a_p;
  _Euler2DLeftEigenvectors_(a_u,a_L,param,a_dir);
  return(0);
}

int Euler2DRightEigenvectors(double *a_u,double *a_R,void *a_p,int a_dir)
{
  Euler2D *param  = (Euler2D*)  a_p;
  _Euler2DRightEigenvectors_(a_u,a_R,param,a_dir);
  return(0);
}
