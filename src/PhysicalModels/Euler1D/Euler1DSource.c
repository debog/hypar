/*! @file Euler1DSource.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute the gravitation source terms for the 1D Euler equations.
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <mpivars.h>
#include <hypar.h>

static int Euler1DSourceFunction(double*,double*,double*,void*,void*,double);

/*! Compute the gravitational source terms for the 1D Euler equations. The source term
    is computed according to the balanced formulation introduced in the reference below.
    The source term is reformulated and "discretized" in a similar fashion as the hyperbolic
    flux to ensure that the hydrostatic balance is maintained to machine precision.
    + Xing, Shu, "High Order Well-Balanced WENO Scheme for the Gas Dynamics Equations
                  Under Gravitational Fields", J. Sci. Comput., 54, 2013, pp. 645--662,
                  http://dx.doi.org/10.1007/s10915-012-9585-8.
*/
int Euler1DSource(
                  double  *a_source, /*!< Computed source terms (array size & layout same as a_u) */
                  double  *a_u,      /*!< Solution (conserved variables) */
                  void    *a_s,      /*!< Solver object of type #HyPar */
                  void    *a_m,      /*!< MPI object of type #MPIVariables */
                  double  a_t   /*!< Current solution time */
                 )
{
  HyPar         *solver = (HyPar* ) a_s;
  MPIVariables  *mpi = (MPIVariables*) a_m;
  Euler1D       *param  = (Euler1D*) solver->m_physics;

  if (param->m_grav == 0.0)  return(0); /* no gravitational forces */

  int     v, done, p, p1, p2;
  double  *SourceI = solver->m_flux_i; /* interace a_source term       */
  double  *SourceC = solver->m_flux_c; /* cell-centered a_source term  */
  double  *SourceL = solver->m_f_l;
  double  *SourceR = solver->m_f_r;

  int     ndims   = solver->m_ndims;
  int ghosts = solver->m_ghosts;
  int     *dim    = solver->m_dim_local;
  double  *x      = solver->m_x;
  double  *dxinv  = solver->m_dxinv;
  int     index[ndims],index1[ndims],index2[ndims],dim_interface[ndims];

  /* set interface dimensions */
  _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_XDIR_]++;
  /* calculate the split a_source function exp(-phi/RT) */
  IERR Euler1DSourceFunction(SourceC,a_u,x,solver,mpi,a_t); CHECKERR(ierr);
  /* calculate the left and right interface a_source terms */
  IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,a_u,x, 1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,a_u,x,-1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
  /* calculate the final interface a_source term */
  IERR param->SourceUpwind(SourceI,SourceL,SourceR,a_u,_XDIR_,solver,a_t);
  /* calculate the final cell-centered a_source term */
  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_XDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
    double dx_inverse;   _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dx_inverse);
    double rho, vel, e, P; _Euler1DGetFlowVar_((a_u+_MODEL_NVARS_*p),rho,vel,e,P,param);
    double term[_MODEL_NVARS_] = {0.0, rho, rho*vel};
    for (v=0; v<_MODEL_NVARS_; v++) {
      a_source[_MODEL_NVARS_*p+v] += (  (term[v]*(1.0/param->m_grav_field[p]))
                                    * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse );
    }
    vel = P; /* useless statement to avoid compiler warning */
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(0);
}

/*! Compute the gravitational source function that is then "discretized" in a way similar to
    the hyperbolic flux function for the balanced formulation introduced in the reference below.
    The source term is reformulated and "discretized" in a similar fashion as the hyperbolic
    flux to ensure that the hydrostatic balance is maintained to machine precision.
    + Xing, Shu, "High Order Well-Balanced WENO Scheme for the Gas Dynamics Equations
                  Under Gravitational Fields", J. Sci. Comput., 54, 2013, pp. 645--662,
                  http://dx.doi.org/10.1007/s10915-012-9585-8.
*/
int Euler1DSourceFunction(
                          double  *a_f, /*!< Computed source function (array size and layout same as a_u) */
                          double  *a_u, /*!< Solution (conserved variables) */
                          double  *a_x, /*!< Spatial coordinates */
                          void    *a_s, /*!< Solver object of type #HyPar */
                          void    *a_m, /*!< MPI object of type #MPIVariables */
                          double  a_t   /*!< Current solution time */
                         )
{
  HyPar         *solver = (HyPar* )       a_s;
  Euler1D       *param  = (Euler1D*)      solver->m_physics;

  int ghosts = solver->m_ghosts;
  int     *dim    = solver->m_dim_local;
  int     ndims   = solver->m_ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i; for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    (a_f+_MODEL_NVARS_*p)[0] = 0.0;
    (a_f+_MODEL_NVARS_*p)[1] = param->m_grav_field[p];
    (a_f+_MODEL_NVARS_*p)[2] = param->m_grav_field[p];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
