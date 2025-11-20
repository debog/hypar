/*! @file ParabolicFunctionNC1Stage.c
    @author Debojyoti Ghosh
    @brief Evaluate the parabolic term using a 1-stage finite-difference discretization.
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Evaluate the parabolic term using a conservative finite-difference spatial discretization:
    The parabolic term is assumed to be of the form:
    \f{equation}{
      {\bf P}\left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {\partial^2 {\bf g}_d\left(\bf u\right)} {\partial x_d^2},
    \f}
    which is discretized at a grid point as:
    \f{equation}{
      \left.{\bf P}\left({\bf u}\right)\right|_j = \sum_{d=0}^{D-1} \frac { \mathcal{L}_d \left[ {\bf g}_d \right]  } {\Delta x_d^2},
    \f}
    where \f$d\f$ is the spatial dimension index, \f$D\f$ is the total number of spatial dimensions (#HyPar::m_ndims), and \f$j\f$ is
    the grid index along \f$d\f$. \f$\mathcal{L}_d\f$ represents the finite-difference approximation to the Laplacian along the \f$d\f$
    spatial dimension, computed using #HyPar::SecondDerivativePar.

    \b Note: this form of the parabolic term \b does \b not allow for cross-derivatives.

    To use this form of the parabolic term:
    + specify \b "par_space_type" in solver.inp as \b "nonconservative-1stage" (#HyPar::m_spatial_type_par).
    + the physical model must specify \f${\bf g}_d\left({\bf u}\right)\f$ through #HyPar::GFunction.
*/
int ParabolicFunctionNC1Stage(
                                double  *a_par, /*!< array to hold the computed parabolic term */
                                double  *a_u,   /*!< solution */
                                void    *a_s,   /*!< Solver object of type #HyPar */
                                void    *a_m,   /*!< MPI object of type #MPIVariables */
                                double  a_t     /*!< Current simulation time */
                             )
{
  HyPar         *solver = (HyPar*)        a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  double        *Func   = solver->m_flux_c;
  double        *Deriv2 = solver->m_deriv2;
  int           d, v, i, done;
  _DECLARE_IERR_;

  int     ndims  = solver->m_ndims;
  int     nvars  = solver->m_nvars;
  int     ghosts = solver->m_ghosts;
  int     *dim   = solver->m_dim_local;
  double  *dxinv = solver->m_dxinv;
  int     size   = solver->m_npoints_local_wghosts;

  if (!solver->GFunction) return(0); /* zero parabolic terms */
  solver->m_count_par++;

  int index[ndims];
  _ArraySetValue_(a_par,size*nvars,0.0);

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    int size_deriv = 1; for (i=0; i<ndims; i++) size_deriv *= dim[i];

    /* calculate the diffusion function */
    IERR solver->GFunction(Func,a_u,d,solver,a_t); CHECKERR(ierr);
    IERR solver->SecondDerivativePar(Deriv2,Func,d,solver,mpi); CHECKERR(ierr);

    /* calculate the final term - second derivative of the diffusion function */
    done = 0; _ArraySetValue_(index,ndims,0);
    int p;
    while (!done) {
      _ArrayIndex1D_(ndims,dim,index,ghosts,p);
      for (v=0; v<nvars; v++)
        a_par[nvars*p+v] += (   dxinv[offset+ghosts+index[d]]*dxinv[offset+ghosts+index[d]]
                            * Deriv2[nvars*p+v] );
      _ArrayIncrementIndex_(ndims,dim,index,done);
    }

    offset += dim[d] + 2*ghosts;
  }

  if (solver->m_flag_ib) _ArrayBlockMultiply_(a_par,solver->m_iblank,size,nvars);
  return(0);
}
