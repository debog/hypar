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

static int ApplyPeriodicity(  int     a_dir,    /*!< Spatial dimension along which to apply periodicity */
                              int     a_ndims,  /*!< Number of spatial dimensions */
                              int     *a_size,  /*!< Integer array with the number of grid points in
                                                   each spatial dimension */
                              int     a_ghosts, /*!< Number of ghost points */
                              double  *a_phi    /*!< The array on which to apply the boundary condition */ )
{
  int bounds[a_ndims], index1[a_ndims], index2[a_ndims], offset[a_ndims],
      done, p1 = 0, p2 = 0;
  _ArrayCopy1D_(a_size,bounds,a_ndims); bounds[a_dir] = a_ghosts;

  done = 0; _ArraySetValue_(index1,a_ndims,0);
  while (!done) {
    _ArraySetValue_(offset,a_ndims,0); offset[a_dir] = -a_ghosts;
    _ArrayIndex1DWO_(a_ndims,a_size,index1,offset,a_ghosts,p1);
    _ArrayCopy1D_(index1,index2,a_ndims);
    index2[a_dir] = index1[a_dir] + a_size[a_dir]-a_ghosts;
    _ArrayIndex1D_(a_ndims,a_size,index2,a_ghosts,p2);

    a_phi[p1] = a_phi[p2];
    _ArrayIncrementIndex_(a_ndims,bounds,index1,done);
  }

  done = 0; _ArraySetValue_(index1,a_ndims,0);
  while (!done) {
    _ArraySetValue_(offset,a_ndims,0); offset[a_dir] = a_size[a_dir];
    _ArrayIndex1DWO_(a_ndims,a_size,index1,offset,a_ghosts,p1);
    _ArrayIndex1D_(a_ndims,a_size,index1,a_ghosts,p2);

    a_phi[p1] = a_phi[p2];
    _ArrayIncrementIndex_(a_ndims,bounds,index1,done);
  }
  return(0);
}

/*! Compute the global DOF index for all the grid points: The "global DOF index"
    is the component number (or block component number for #HyPar::m_nvars > 1) of
    a grid point in the global solution vector. It is also the row number (or
    block row number) of the grid point in the global matrix representing, for
    example, the Jacobian of the right-hand-side.

    #PETScContext::globalDOF is an integer array with the same layout as the solution
    array #HyPar::m_u (but with one component) containing the global DOF index for the
    corresponding grid points. It has the same number of ghost points as #HyPar::m_u.
    + This array is initialized to -1.
    + The global DOF indices are computed for all non-ghost grid points.
    + If any boundaries are periodic, periodic boundary conditions are applied to fill
      the appropriate ghost points.
    + Ghost points corresponding to internal (MPI) boundaries are filled using
      MPIExchangeBoundariesnD().
    + Thus, ghost points corresponding to physical, non-periodic boundaries retain the
      initial value of -1.
*/
int PetscGlobalDOF(void* a_c /*!< Object of type #PETScContext*/)
{
  PETScContext* ctxt = (PETScContext*) a_c;
  SimulationObject* sim = (SimulationObject*) ctxt->m_simobj;
  int nsims = ctxt->m_nsims;

  ctxt->m_globalDOF.resize(nsims, nullptr);

  /* compute MPI offset */
  std::vector<int> local_sizes(ctxt->m_nproc ,0);
  local_sizes[ctxt->m_rank] = ctxt->m_npoints;
  MPIMax_integer(local_sizes.data(),local_sizes.data(),ctxt->m_nproc,&sim[0].mpi.m_world);

  int MPIOffset = 0;
  for (int i=0; i<ctxt->m_rank; i++) MPIOffset += local_sizes[i];

  int simOffset = 0;
  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    int   *dim        = solver->m_dim_local,
          ndims       = solver->m_ndims,
          ghosts      = solver->m_ghosts,
          npoints     = solver->m_npoints_local,
          npoints_wg  = solver->m_npoints_local_wghosts,
          nv          = ndims + 1, i;

    ctxt->m_globalDOF[ns] = (double*) calloc( npoints_wg, sizeof(double) );
    _ArraySetValue_(ctxt->m_globalDOF[ns], npoints_wg, -1.0);

    for (int i = 0; i < npoints; i++) {
      int p = (ctxt->m_points[ns]+i*nv)[ndims];
      ctxt->m_globalDOF[ns][p] = (double) (i + simOffset + MPIOffset);
    }

    for (int i=0; i<ndims; i++) {
      if (solver->m_is_periodic[i]) {
        ApplyPeriodicity(i,ndims,dim,ghosts,ctxt->m_globalDOF[ns]);
      }
    }
    MPIExchangeBoundariesnD(ndims,1,dim,ghosts,mpi,ctxt->m_globalDOF[ns]);

    simOffset += npoints;
  }

  return 0;
}

#endif
