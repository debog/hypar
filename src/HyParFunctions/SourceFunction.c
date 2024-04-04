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
                    double  *source,  /*!< the computed source term */
                    double  *u,       /*!< solution */
                    void    *s,       /*!< solver object of type #HyPar */
                    void    *m,       /*!< MPI object of type #MPIVariables */
                    double  t         /*!< Current simulation time */
                  )
{
  HyPar           *solver   = (HyPar*)        s;
  MPIVariables    *mpi      = (MPIVariables*) m;

  /* initialize to zero */
  int size = solver->ndof_cells_wghosts;
#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    gpuArraySetValue(source,size, 0.0);
  } else {
#endif
    _ArraySetValue_(source,size,0.0);
#if defined(HAVE_CUDA)
  }
#endif

  /* call the source function of the physics model, if available */
  if (solver->SFunction) {
    solver->SFunction(source,u,solver,mpi,t);
    solver->count_sou++;
  }

  /* Apart from other source terms, implement sponge BC as a source */
  DomainBoundary* boundary = (DomainBoundary*) solver->boundary;
  int n;
  int nb = solver->nBoundaryZones;
#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    for (n = 0; n < nb; n++) {
      if (!strcmp(boundary[n].bctype,_SPONGE_)) {
        fprintf(stderr,"ERROR: Sponge BC not yet implemented on GPU.\n");
      }
    }
  } else {
#endif
    for (n = 0; n < nb; n++) {
      if (!strcmp(boundary[n].bctype,_SPONGE_)) {
        BCSpongeSource( &boundary[n],
                        solver->ndims,
                        solver->nvars,
                        solver->ghosts,
                        solver->dim_local,
                        solver->x,
                        u,
                        source  );
      }
    }
#if defined(HAVE_CUDA)
  }
#endif


  return(0);
}
