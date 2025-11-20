/*! @file ShallowWater1DJacobian.c
    @author Debojyoti Ghosh
    @brief Contains the functions compute flux Jacobians for the 1D shallow water system
*/

#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/shallowwater1d.h>

/*! Function to compute the flux Jacobian of the 1D shallow water equations, given the
    solution at a grid point. The Jacobian is square matrix of size nvar=2, and
    is returned as a 1D array (double) of 4 elements in row-major format.
*/
int ShallowWater1DJacobian(
                    double  *a_Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 4 */
                    double  *a_u,   /*!< solution at a grid point (array of size nvar = 2) */
                    void    *a_p,   /*!< object containing the physics-related parameters */
                    int     a_dir,  /*!< dimension (x/y/z) (not used, since this is 1D system) */
                    int     a_nvars,/*!< number of vector components */
                    int     a_upw   /*!< 0 -> send back complete Jacobian,
                                       1 -> send back Jacobian of right(+)-moving flux,
                                      -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  ShallowWater1D  *param = (ShallowWater1D*) a_p;
  static double   R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                  L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _ShallowWater1DEigenvalues_      (a_u,D,param,0);
  _ShallowWater1DLeftEigenvectors_ (a_u,L,param,0);
  _ShallowWater1DRightEigenvectors_(a_u,R,param,0);

  int aupw = absolute(a_upw), k;
  k = 0; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+a_upw)*max(0,D[k]) + 0.5*aupw*(1-a_upw)*min(0,D[k]) );
  k = 3; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+a_upw)*max(0,D[k]) + 0.5*aupw*(1-a_upw)*min(0,D[k]) );

  MatMult2(2,DL,D,L);
  MatMult2(2,a_Jac,R,DL);

  return(0);
}
