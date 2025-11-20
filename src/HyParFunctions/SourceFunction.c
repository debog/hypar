/*! @file SourceFunction.c
    @author Debojyoti Ghosh, Youngdae Kim
    @brief Evaluate the source term
*/

#include <stdlib.h>
#include <string.h>
#include <basic.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <boundaryconditions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Evaluate the source term \f${\bf S}\left({\bf u}\right)\f$ in the governing equation,
    if the physical model specifies one. In addition, if the simulation requires a sponge
    boundary treatment, the sponge BC function is called.
*/
int SourceFunction(
                    double  *a_source,  /*!< the computed source term */
                    double  *a_u,       /*!< solution */
                    void    *a_s,       /*!< solver object of type #HyPar */
                    void    *a_m,       /*!< MPI object of type #MPIVariables */
                    double  a_t         /*!< Current simulation time */
                  )
{
  HyPar           *solver   = (HyPar*)        a_s;
  MPIVariables    *mpi      = (MPIVariables*) a_m;

  /* initialize to zero */
  int size = solver->m_ndof_cells_wghosts;
#if defined(HAVE_CUDA)
  if (solver->m_use_gpu) {
    gpuArraySetValue(a_source,size, 0.0);
  } else {
#endif
    _ArraySetValue_(a_source,size,0.0);
#if defined(HAVE_CUDA)
  }
#endif

  /* call the a_source function of the physics model, if available */
  if (solver->SFunction) {
    solver->SFunction(a_source,a_u,solver,mpi,a_t);
    solver->m_count_sou++;
  }

  /* Apart from other a_source terms, implement sponge BC as a a_source */
  DomainBoundary* boundary = (DomainBoundary*) solver->m_boundary;
  int n;
  int nb = solver->m_n_boundary_zones;
#if defined(HAVE_CUDA)
  if (solver->m_use_gpu) {
    for (n = 0; n < nb; n++) {
      if (!strcmp(boundary[n].m_bctype,_SPONGE_)) {
        fprintf(stderr,"ERROR: Sponge BC not yet implemented on GPU.\n");
      }
    }
  } else {
#endif
    for (n = 0; n < nb; n++) {
      if (!strcmp(boundary[n].m_bctype,_SPONGE_)) {
        BCSpongeSource( &boundary[n],
                        solver->m_ndims,
                        solver->m_nvars,
                        solver->m_ghosts,
                        solver->m_dim_local,
                        solver->m_x,
                        a_u,
                        a_source  );
      }
    }
#if defined(HAVE_CUDA)
  }
#endif


  return(0);
}
