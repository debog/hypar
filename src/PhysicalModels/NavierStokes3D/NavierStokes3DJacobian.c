/*! @file NavierStokes3DJacobian.c
    @author Debojyoti Ghosh
    @brief Contains the functions compute flux Jacobians for the 3D Navier-Stokes system
*/

#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/navierstokes3d.h>

/*! Function to compute the flux Jacobian of the 3D Navier-Stokes equations, given the
    solution at a grid point. The Jacobian is square matrix of size nvar=5, and
    is returned as a 1D array (double) of 25 elements in row-major format.
*/
int NavierStokes3DJacobian(
                    double  *Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 25 */
                    double  *u,   /*!< solution at a grid point (array of size nvar = 5) */
                    void    *p,   /*!< object containing the physics-related parameters */
                    int     dir,  /*!< dimension (0 -> x, 1 -> y, 2 -> z) */
                    int     nvars,/*!< number of vector components */
                    int     upw   /*!< 0 -> send back complete Jacobian,
                                       1 -> send back Jacobian of right(+)-moving flux,
                                      -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  NavierStokes3D *param = (NavierStokes3D*) p;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _NavierStokes3DEigenvalues_      (u,_NavierStokes3D_stride_,D,param->gamma,dir);
  _NavierStokes3DLeftEigenvectors_ (u,_NavierStokes3D_stride_,L,param->gamma,dir);
  _NavierStokes3DRightEigenvectors_(u,_NavierStokes3D_stride_,R,param->gamma,dir);

  int aupw = absolute(upw), k;
  k = 0;  D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 6;  D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 12; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 18; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  k = 24; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );

  MatMult5(_MODEL_NVARS_,DL,D,L);
  MatMult5(_MODEL_NVARS_,Jac,R,DL);

  return(0);
}

/*! Function to compute the Jacobian of the fast flux (representing the acoustic waves)
    of the 3D Navier-Stokes equations, given the solution at a grid point. The Jacobian
    is square matrix of size nvar=5, and is returned as a 1D array (double) of 25 elements
    in row-major format.
*/
int NavierStokes3DStiffJacobian(
                          double  *Jac, /*!< Jacobian matrix: 1D array of size nvar^2 = 25 */
                          double  *u,   /*!< solution at a grid point (array of size nvar = 5) */
                          void    *p,   /*!< object containing the physics-related parameters */
                          int     dir,  /*!< dimension (0 -> x, 1 -> y, 2 -> z) */
                          int     nvars,/*!< number of vector components */
                          int     upw   /*!< 0 -> send back complete Jacobian,
                                             1 -> send back Jacobian of right(+)-moving flux,
                                            -1 -> send back Jacobian of left(-)-moving flux */
                   )
{
  NavierStokes3D *param = (NavierStokes3D*) p;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* get the eigenvalues and left,right eigenvectors */
  _NavierStokes3DEigenvalues_      (u,_NavierStokes3D_stride_,D,param->gamma,dir);
  _NavierStokes3DLeftEigenvectors_ (u,_NavierStokes3D_stride_,L,param->gamma,dir);
  _NavierStokes3DRightEigenvectors_(u,_NavierStokes3D_stride_,R,param->gamma,dir);

  int aupw = absolute(upw), k;
  if (dir == _XDIR_) {
    k = 0;  D[k] = 0.0;
    k = 6;  D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
    k = 12; D[k] = 0.0;
    k = 18; D[k] = 0.0;
    k = 24; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  } else if (dir == _YDIR_) {
    k = 0;  D[k] = 0.0;
    k = 6;  D[k] = 0.0;
    k = 12; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
    k = 18; D[k] = 0.0;
    k = 24; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  } else if (dir == _ZDIR_) {
    k = 0;  D[k] = 0.0;
    k = 6;  D[k] = 0.0;
    k = 12; D[k] = 0.0;
    k = 18; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
    k = 24; D[k] = absolute( (1-aupw)*D[k] + 0.5*aupw*(1+upw)*max(0,D[k]) + 0.5*aupw*(1-upw)*min(0,D[k]) );
  }

  MatMult5(_MODEL_NVARS_,DL,D,L);
  MatMult5(_MODEL_NVARS_,Jac,R,DL);

  return(0);
}
