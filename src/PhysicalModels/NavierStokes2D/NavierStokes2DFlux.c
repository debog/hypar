/*! @file NavierStokes2DFlux.c
    @author Debojyoti Ghosh
    @brief Functions to compute the hyperbolic flux for 2D Navier-Stokes equations
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

/*!
  Compute the hyperbolic flux function for the 2D Navier-Stokes equations:
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ (e+p)v \end{array}\right]
  \f}
  Note: the flux function needs to be computed at the ghost points as well.
*/
int NavierStokes2DFlux(
                        double  *a_f, /*!< Array to hold the computed flux vector (same layout as a_u) */
                        double  *a_u, /*!< Array with the solution vector */
                        int     a_dir,/*!< Spatial dimension (x or y) for which to compute the flux */
                        void    *a_s, /*!< Solver object of type #HyPar */
                        double  a_t   /*!< Current simulation time */
                      )
{
  HyPar           *solver = (HyPar*)   a_s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->m_physics;
  int             *dim    = solver->m_dim_local;
  int             ghosts  = solver->m_ghosts;
  static int      index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,_MODEL_NDIMS_);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    double rho, vx, vy, e, P;
    _NavierStokes2DGetFlowVar_((a_u+_MODEL_NVARS_*p),rho,vx,vy,e,P,param->m_gamma);
    _NavierStokes2DSetFlux_((a_f+_MODEL_NVARS_*p),rho,vx,vy,e,P,a_dir);
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}

/*! Compute the stiff flux, given the solution vector. The stiff flux is the component of
    the total flux that represents the acoustic modes (see #_NavierStokes2DSetStiffFlux_).
    Here, the linearized approximation to the stiff flux is computed as:
    \f{equation}{
      {\bf f}\left({\bf u}\right) = A_f{\bf u}
    \f}
    where \f$A_f = A_f\left({\bf u}_{ref}\right)\f$ is the fast Jacobian (#NavierStokes2D::fast_jac)
    evaluated for the solution at the beginning of each time step (\f${\bf u}_{ref}\f$ is
    #NavierStokes2D::solution). This is done in NavierStokes2DPreStep().\n\n
  Note: the flux function needs to be computed at the ghost points as well.
*/
int NavierStokes2DStiffFlux(
                              double  *a_f, /*!< Array to hold the computed flux vector (same layout as a_u) */
                              double  *a_u, /*!< Array with the solution vector */
                              int     a_dir,/*!< Spatial dimension (x or y) for which to compute the flux */
                              void    *a_s, /*!< Solver object of type #HyPar */
                              double  a_t   /*!< Current simulation time */
                           )
{
  HyPar             *solver = (HyPar*)   a_s;
  NavierStokes2D    *param  = (NavierStokes2D*) solver->m_physics;
  int               *dim    = solver->m_dim_local;
  int               ghosts  = solver->m_ghosts;
  static const int  JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,_MODEL_NDIMS_);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    double *Af = param->m_fast_jac+(2*p+a_dir)*JacSize;
    MatVecMult4(_MODEL_NVARS_,(a_f+_MODEL_NVARS_*p),Af,(a_u+_MODEL_NVARS_*p));
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}

/*! Compute the non-stiff flux, given the solution vector. The non-stiff flux is the component of
    the total flux that represents the entropy modes (see #_NavierStokes2DSetNonStiffFlux_).
    Here, the linearized approximation to the non-stiff flux is computed as:
    \f{equation}{
      {\bf f}\left({\bf u}\right) - A_f{\bf u}
    \f}
    where \f${\bf f}\left({\bf u}\right)\f$ is the total flux computed in NavierStokes2DFlux(),
    and \f$A_f{\bf u}\f$ is the linearized stiff flux computed in NavierStokes2DStiffFlux().\n\n
  Note: the flux function needs to be computed at the ghost points as well.
*/
int NavierStokes2DNonStiffFlux(
                              double  *a_f, /*!< Array to hold the computed flux vector (same layout as a_u) */
                              double  *a_u, /*!< Array with the solution vector */
                              int     a_dir,/*!< Spatial dimension (x or y) for which to compute the flux */
                              void    *a_s, /*!< Solver object of type #HyPar */
                              double  a_t   /*!< Current simulation time */
                           )
{
  HyPar             *solver = (HyPar*)   a_s;
  NavierStokes2D    *param  = (NavierStokes2D*) solver->m_physics;
  int               *dim    = solver->m_dim_local;
  int               ghosts  = solver->m_ghosts;
  static const int  JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];
  static double     ftot[_MODEL_NVARS_], fstiff[_MODEL_NVARS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,_MODEL_NDIMS_);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    /* compute total flux */
    double rho, vx, vy, e, P;
    _NavierStokes2DGetFlowVar_((a_u+_MODEL_NVARS_*p),rho,vx,vy,e,P,param->m_gamma);
    _NavierStokes2DSetFlux_(ftot,rho,vx,vy,e,P,a_dir);
    /* compute stiff stuff */
    double *Af = param->m_fast_jac+(2*p+a_dir)*JacSize;
    MatVecMult4(_MODEL_NVARS_,fstiff,Af,(a_u+_MODEL_NVARS_*p));
    /* subtract stiff flux from total flux */
    _ArraySubtract1D_((a_f+_MODEL_NVARS_*p),ftot,fstiff,_MODEL_NVARS_);
    /* Done */
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}
