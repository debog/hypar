/*! @file InitialSolution.c
    @author Debojyoti Ghosh, Youngdae Kim
    @brief Read in initial solution from file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <common.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <io.h>
#include <mpivars.h>
#include <simulation_object.h>

int VolumeIntegral(double*,double*,void*,void*);

/*! Read in initial solution from file, and compute grid spacing
    and volume integral of the initial solution */
int InitialSolution ( void  *s,   /*!< Array of simulation objects of type #SimulationObject */
                      int   nsims /*!< Number of simulation objects */
                    )
{
  SimulationObject* simobj = (SimulationObject*) s;
  int n, flag, d, i, offset, ierr;

  for (n = 0; n < nsims; n++) {

    int ghosts = simobj[n].solver.m_ghosts;

    char fname_root[_MAX_STRING_SIZE_] = "initial";
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(n, index, (int)log10(nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }

    ierr = ReadArray( simobj[n].solver.m_ndims,
                      simobj[n].solver.m_nvars,
                      simobj[n].solver.m_dim_global,
                      simobj[n].solver.m_dim_local,
                      simobj[n].solver.m_ghosts,
                      &(simobj[n].solver),
                      &(simobj[n].mpi),
                      simobj[n].solver.m_x,
                      simobj[n].solver.m_u,
                      fname_root,
                      &flag );
    if (ierr) {
      fprintf(stderr, "Error in InitialSolution() on rank %d.\n",
              simobj[n].mpi.m_rank);
      return ierr;
    }
    if (!flag) {
      fprintf(stderr,"Error: initial solution file not found.\n");
      return(1);
    }
    CHECKERR(ierr);

    /* exchange MPI-boundary values of u between processors */
    MPIExchangeBoundariesnD(  simobj[n].solver.m_ndims,
                              simobj[n].solver.m_nvars,
                              simobj[n].solver.m_dim_local,
                              simobj[n].solver.m_ghosts,
                              &(simobj[n].mpi),
                              simobj[n].solver.m_u  );

    /* calculate dxinv */
    offset = 0;
    for (d = 0; d < simobj[n].solver.m_ndims; d++) {
      for (i = 0; i < simobj[n].solver.m_dim_local[d]; i++) {
        simobj[n].solver.m_dxinv[i+offset+ghosts]
          = 2.0 / (simobj[n].solver.m_x[i+1+offset+ghosts]-simobj[n].solver.m_x[i-1+offset+ghosts]);
      }
      offset += (simobj[n].solver.m_dim_local[d] + 2*ghosts);
    }

    /* exchange MPI-boundary values of dxinv between processors */
    offset = 0;
    for (d = 0; d < simobj[n].solver.m_ndims; d++) {
      ierr = MPIExchangeBoundaries1D( &(simobj[n].mpi),
                                      &(simobj[n].solver.m_dxinv[offset]),
                                      simobj[n].solver.m_dim_local[d],
                                      ghosts,
                                      d,
                                      simobj[n].solver.m_ndims ); CHECKERR(ierr);
      if (ierr) {
        fprintf(stderr, "Error in InitialSolution() on rank %d.\n",
                simobj[n].mpi.m_rank);
        return ierr;
      }
      offset += (simobj[n].solver.m_dim_local[d] + 2*ghosts);
    }

    /* fill in ghost values of dxinv at physical boundaries by extrapolation */
    offset = 0;
    for (d = 0; d < simobj[n].solver.m_ndims; d++) {
      double *dxinv = &(simobj[n].solver.m_dxinv[offset]);
      int    *dim = simobj[n].solver.m_dim_local;
      if (simobj[n].mpi.m_ip[d] == 0) {
        /* fill left boundary along this dimension */
        for (i = 0; i < ghosts; i++) dxinv[i] = dxinv[ghosts];
      }
      if (simobj[n].mpi.m_ip[d] == simobj[n].mpi.m_iproc[d]-1) {
        /* fill right boundary along this dimension */
        for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) dxinv[i] = dxinv[dim[d]+ghosts-1];
      }
      offset  += (dim[d] + 2*ghosts);
    }

    /* calculate volume integral of the initial solution */
    ierr = VolumeIntegral(  simobj[n].solver.m_volume_integral_initial,
                            simobj[n].solver.m_u,
                            &(simobj[n].solver),
                            &(simobj[n].mpi) ); CHECKERR(ierr);
    if (ierr) {
      fprintf(stderr, "Error in InitialSolution() on rank %d.\n",
              simobj[n].mpi.m_rank);
      return ierr;
    }
    if (!simobj[n].mpi.m_rank) {
      if (nsims > 1) printf("Volume integral of the initial solution on domain %d:\n", n);
      else           printf("Volume integral of the initial solution:\n");
      for (d=0; d<simobj[n].solver.m_nvars; d++) {
        printf("%2d:  %1.16E\n",d,simobj[n].solver.m_volume_integral_initial[d]);
      }
    }
    /* Set initial total boundary flux integral to zero */
    _ArraySetValue_(simobj[n].solver.m_total_boundary_integral,simobj[n].solver.m_nvars,0);

  }

#if defined(HAVE_CUDA)
  if (simobj[0].solver.m_use_gpu) {
    for (int n = 0; n < nsims; n++) {
      int npoints_local_wghosts = simobj[n].solver.m_npoints_local_wghosts;
      int nvars                 = simobj[n].solver.m_nvars;
      int size_x                = simobj[n].solver.m_size_x;

      gpuMemcpy(simobj[n].solver.m_gpu_x,
                simobj[n].solver.m_x, simobj[n].solver.m_size_x*sizeof(double),
                gpuMemcpyHostToDevice);
      gpuMemcpy(simobj[n].solver.m_gpu_dxinv, simobj[n].solver.m_dxinv,
                simobj[n].solver.m_size_x*sizeof(double),
                gpuMemcpyHostToDevice);

#ifdef CUDA_VAR_ORDERDING_AOS
      gpuMemcpy(simobj[n].solver.m_gpu_u,
                simobj[n].solver.m_u,
                simobj[n].solver.m_ndof_cells_wghosts*sizeof(double),
                gpuMemcpyHostToDevice);
#else
      double *h_u = (double *) malloc(simobj[n].solver.m_ndof_cells_wghosts*sizeof(double));
      for (int i=0; i<npoints_local_wghosts; i++) {
        for (int v=0; v<nvars; v++) {
          h_u[i+v*npoints_local_wghosts] = simobj[n].solver.m_u[i*nvars+v];
        }
      }
      gpuMemcpy(simobj[n].solver.m_gpu_u,
                h_u,
                simobj[n].solver.m_ndof_cells_wghosts*sizeof(double),
                gpuMemcpyHostToDevice);
      free(h_u);
#endif
    }
  }
#endif

  return 0;
}
