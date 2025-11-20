/*! @file NonLinearInterpolation.c
    @author Debojyoti Ghosh
    @brief Compute the coefficients of a non-linear interpolation method (such as WENO, CRWENO)
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the interpolation coefficients of a nonlinear interpolation method, based on
    the smoothness of the flux function of the solution passed as an argument. This function
    is provided separately so that these coefficients can be pre-computed and stored for
    future use.

    In the implementation of non-linear methods (such as WENO5), the calculation
    of the non-linear coefficients, and the actual evaluation of the interpolant are separated
    into different functions. This provides flexibility on when the nonlinear coefficients are
    evaluated. Some scenarios are as follows:
    + For explicit time integration, it will be computed every time the hyperbolic flux term
      is being computed.
    + For implicit time integration, consistency or linearization  may require that the
      coefficients be computed and "frozen" at the beginning of each stage. Thus, this function
      can be called at the beginning of each time integration stage, while HyperbolicFunction()
      is called with the argument \b LimFlag = 0.
*/
int NonLinearInterpolation(
                            double  *a_u, /*!< Solution array */
                            void    *a_s, /*!< Solver object of type #HyPar */
                            void    *a_m, /*!< MPI object of type #MPIVariables */
                            double  a_t,  /*!< Current solution time */
                            /*! The flux function \f${\bf f}_d\left({\bf u}\a_right)\f$, whose
                            properties (typically smoothness) is used to evaluate the nonlinear
                            coefficients */
                            int     (*FluxFunction)(double*,double*,int,void*,double)
                          )
{
  HyPar         *solver = (HyPar*) a_s;
  MPIVariables  *mpi    = (MPIVariables*)   a_m;
  double        *FluxC  = solver->m_flux_c; /* cell centered flux */
  _DECLARE_IERR_;

  int flag = (FluxFunction && solver->m_flag_nonlinearinterp && solver->SetInterpLimiterVar);
  if (flag) {;
    int     ndims  = solver->m_ndims;
    int     nvars  = solver->m_nvars;
    int     ghosts = solver->m_ghosts;
    int     *dim   = solver->m_dim_local;
    double  *x     = solver->m_x;

    int size = 1, d;
    for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

    /* apply boundary conditions and exchange data over MPI interfaces */
    IERR solver->ApplyBoundaryConditions(solver,mpi,a_u,NULL,a_t);      CHECKERR(ierr);
    IERR solver->ApplyIBConditions(solver,mpi,a_u,a_t);                 CHECKERR(ierr);
    IERR MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,a_u);     CHECKERR(ierr);

    int offset = 0;
    for (d = 0; d < ndims; d++) {
      /* evaluate cell-centered flux */
      IERR FluxFunction(FluxC,a_u,d,solver,a_t);                          CHECKERR(ierr);
      /* calculate non-linear interpolation coefficients */
      IERR solver->SetInterpLimiterVar(FluxC,a_u,x+offset,d,solver,mpi);CHECKERR(ierr);
      offset += dim[d] + 2*ghosts;
    }
  }

  return(0);
}

