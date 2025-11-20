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
int PetscCreatePointList(void *a_obj /*!< Object of type #PETScContext */)
{
  PETScContext* ctxt = (PETScContext*) a_obj;
  SimulationObject* sim = (SimulationObject*) ctxt->m_simobj;

  int nsims = ctxt->m_nsims;
  int ndims = sim[0].solver.m_ndims;

  /* count the number of computational points */
  ctxt->m_npoints = 0;
  ctxt->m_ndofs = 0;
  ctxt->m_offsets = (int*) calloc(nsims, sizeof(int));
  for (int ns = 0; ns < nsims; ns++) {
    ctxt->m_npoints += sim[ns].solver.m_npoints_local;
    ctxt->m_offsets[ns] = ctxt->m_ndofs;
    ctxt->m_ndofs += (sim[ns].solver.m_npoints_local*sim[ns].solver.m_nvars);
  }

  int nv = ndims+1;
  ctxt->m_points.resize(nsims, nullptr);

  for (int ns = 0; ns < nsims; ns++) {

    int npoints = sim[ns].solver.m_npoints_local;
    ctxt->m_points[ns] = (int*) calloc (npoints*nv,sizeof(int));

    const int* const dim( sim[ns].solver.m_dim_local );
    const int ghosts( sim[ns].solver.m_ghosts );

    int done = 0, i = 0;
    std::vector<int> index(ndims,0);
    while (!done) {
      int p; _ArrayIndex1D_(ndims, dim, index.data(), ghosts, p);
      _ArrayCopy1D_(index.data(), (ctxt->m_points[ns]+i*nv), ndims);
      (ctxt->m_points[ns]+i*nv)[ndims] = p;
      _ArrayIncrementIndex_(ndims,dim,index,done);
      i++;
    }

    if (i != npoints) {
      fprintf(stderr,"Error in PetscCreatePointList() on rank %d:\n", sim[ns].mpi.m_rank);
      fprintf(stderr,"Inconsistency in point count - %d, %d.\n",
              i, npoints);
      return 1;
    }
  }

  int global_npoints;
  int global_ndofs;
  MPISum_integer(&global_npoints,&(ctxt->m_npoints),1,&sim[0].mpi.m_world);
  MPISum_integer(&global_ndofs,&(ctxt->m_ndofs),1,&sim[0].mpi.m_world);
  if (!ctxt->m_rank) {
    printf("PETSc:    total number of computational points is %d.\n",global_npoints);
    printf("PETSc:    total number of computational DOFs is %d.\n",global_ndofs);
  }

  return 0;
}

#endif
