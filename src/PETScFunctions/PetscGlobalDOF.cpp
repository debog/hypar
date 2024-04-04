/*! @file PetscGlobalDOF.cpp
    @author Debojyoti Ghosh
    @brief Compute the global DOF index for all the grid points
*/

#ifdef with_petsc

#include <stdlib.h>
#include <vector>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>
#include <petscinterface.h>

static int ApplyPeriodicity(  int     dir,    /*!< Spatial dimension along which to apply periodicity */
                              int     ndims,  /*!< Number of spatial dimensions */
                              int     *size,  /*!< Integer array with the number of grid points in
                                                   each spatial dimension */
                              int     ghosts, /*!< Number of ghost points */
                              double  *phi    /*!< The array on which to apply the boundary condition */ )
{
  int bounds[ndims], index1[ndims], index2[ndims], offset[ndims],
      done, p1 = 0, p2 = 0;
  _ArrayCopy1D_(size,bounds,ndims); bounds[dir] = ghosts;

  done = 0; _ArraySetValue_(index1,ndims,0);
  while (!done) {
    _ArraySetValue_(offset,ndims,0); offset[dir] = -ghosts;
    _ArrayIndex1DWO_(ndims,size,index1,offset,ghosts,p1);
    _ArrayCopy1D_(index1,index2,ndims);
    index2[dir] = index1[dir] + size[dir]-ghosts;
    _ArrayIndex1D_(ndims,size,index2,ghosts,p2);

    phi[p1] = phi[p2];
    _ArrayIncrementIndex_(ndims,bounds,index1,done);
  }

  done = 0; _ArraySetValue_(index1,ndims,0);
  while (!done) {
    _ArraySetValue_(offset,ndims,0); offset[dir] = size[dir];
    _ArrayIndex1DWO_(ndims,size,index1,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,size,index1,ghosts,p2);

    phi[p1] = phi[p2];
    _ArrayIncrementIndex_(ndims,bounds,index1,done);
  }
  return(0);
}

/*! Compute the global DOF index for all the grid points: The "global DOF index"
    is the component number (or block component number for #HyPar::nvars > 1) of
    a grid point in the global solution vector. It is also the row number (or
    block row number) of the grid point in the global matrix representing, for
    example, the Jacobian of the right-hand-side.

    #PETScContext::globalDOF is an integer array with the same layout as the solution
    array #HyPar::u (but with one component) containing the global DOF index for the
    corresponding grid points. It has the same number of ghost points as #HyPar::u.
    + This array is initialized to -1.
    + The global DOF indices are computed for all non-ghost grid points.
    + If any boundaries are periodic, periodic boundary conditions are applied to fill
      the appropriate ghost points.
    + Ghost points corresponding to internal (MPI) boundaries are filled using
      MPIExchangeBoundariesnD().
    + Thus, ghost points corresponding to physical, non-periodic boundaries retain the
      initial value of -1.
*/
int PetscGlobalDOF(void* c /*!< Object of type #PETScContext*/)
{
  PETScContext* ctxt = (PETScContext*) c;
  SimulationObject* sim = (SimulationObject*) ctxt->simobj;
  int nsims = ctxt->nsims;

  ctxt->globalDOF.resize(nsims, nullptr);

  /* compute MPI offset */
  std::vector<int> local_sizes(ctxt->nproc ,0);
  local_sizes[ctxt->rank] = ctxt->npoints;
  MPIMax_integer(local_sizes.data(),local_sizes.data(),ctxt->nproc,&sim[0].mpi.world);

  int MPIOffset = 0;
  for (int i=0; i<ctxt->rank; i++) MPIOffset += local_sizes[i];

  int simOffset = 0;
  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    int   *dim        = solver->dim_local,
          ndims       = solver->ndims,
          ghosts      = solver->ghosts,
          npoints     = solver->npoints_local,
          npoints_wg  = solver->npoints_local_wghosts,
          nv          = ndims + 1, i;

    ctxt->globalDOF[ns] = (double*) calloc( npoints_wg, sizeof(double) );
    _ArraySetValue_(ctxt->globalDOF[ns], npoints_wg, -1.0);

    for (int i = 0; i < npoints; i++) {
      int p = (ctxt->points[ns]+i*nv)[ndims];
      ctxt->globalDOF[ns][p] = (double) (i + simOffset + MPIOffset);
    }

    for (int i=0; i<ndims; i++) {
      if (solver->isPeriodic[i]) {
        ApplyPeriodicity(i,ndims,dim,ghosts,ctxt->globalDOF[ns]);
      }
    }
    MPIExchangeBoundariesnD(ndims,1,dim,ghosts,mpi,ctxt->globalDOF[ns]);

    simOffset += npoints;
  }

  return 0;
}

#endif
