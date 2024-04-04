/*! @file PetscCreatePointList.c
    @brief Create list of computational points.
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <arrayfunctions.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>
#include <petscinterface_struct.h>

/*!
  Create a list of computational points for each simulation
  domain: This is a list of all the
  grid points on which the PDE is solved. Thus, it is the total
  number of grid points minus the ghost points and blanked out
  points.

  Note: this list is local, not global.
*/
int PetscCreatePointList(void *obj /*!< Object of type #PETScContext */)
{
  PETScContext* ctxt = (PETScContext*) obj;
  SimulationObject* sim = (SimulationObject*) ctxt->simobj;

  int nsims = ctxt->nsims;
  int ndims = sim[0].solver.ndims;

  /* count the number of computational points */
  ctxt->npoints = 0;
  ctxt->ndofs = 0;
  ctxt->offsets = (int*) calloc(nsims, sizeof(int));
  for (int ns = 0; ns < nsims; ns++) {
    ctxt->npoints += sim[ns].solver.npoints_local;
    ctxt->offsets[ns] = ctxt->ndofs;
    ctxt->ndofs += (sim[ns].solver.npoints_local*sim[ns].solver.nvars);
  }

  int nv = ndims+1;
  ctxt->points.resize(nsims, nullptr);

  for (int ns = 0; ns < nsims; ns++) {

    int npoints = sim[ns].solver.npoints_local;
    ctxt->points[ns] = (int*) calloc (npoints*nv,sizeof(int));

    const int* const dim( sim[ns].solver.dim_local );
    const int ghosts( sim[ns].solver.ghosts );

    int done = 0, i = 0;
    std::vector<int> index(ndims,0);
    while (!done) {
      int p; _ArrayIndex1D_(ndims, dim, index.data(), ghosts, p);
      _ArrayCopy1D_(index.data(), (ctxt->points[ns]+i*nv), ndims);
      (ctxt->points[ns]+i*nv)[ndims] = p;
      _ArrayIncrementIndex_(ndims,dim,index,done);
      i++;
    }

    if (i != npoints) {
      fprintf(stderr,"Error in PetscCreatePointList() on rank %d:\n", sim[ns].mpi.rank);
      fprintf(stderr,"Inconsistency in point count - %d, %d.\n",
              i, npoints);
      return 1;
    }
  }

  int global_npoints;
  int global_ndofs;
  MPISum_integer(&global_npoints,&(ctxt->npoints),1,&sim[0].mpi.world);
  MPISum_integer(&global_ndofs,&(ctxt->ndofs),1,&sim[0].mpi.world);
  if (!ctxt->rank) {
    printf("PETSc:    total number of computational points is %d.\n",global_npoints);
    printf("PETSc:    total number of computational DOFs is %d.\n",global_ndofs);
  }

  return 0;
}

#endif
