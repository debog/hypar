/*! @file Initialize.c
    @author Debojyoti Ghosh, Youngdae Kim
    @brief Initialization function
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <mpivars.h>
#include <simulation_object.h>

/*! Initialization function called at the beginning of a simulation. This function
    does the following:
    + allocates memory for MPI related arrays
    + initializes the values for MPI variables
    + creates sub-communicators and communication groups
    + allocates memory for arrays to store solution, right-hand-side,
      flux, and other working vectors.
    + initializes function counters to zero
*/
int Initialize( void *s,    /*!< Array of simulation objects of type #SimulationObject */
                int  nsims  /*!< Number of simulation objects */
              )
{
  SimulationObject* simobj = (SimulationObject*) s;
  int i,d,n;

  if (nsims == 0) {
    return 1;
  }

#if defined(HAVE_CUDA)
  if (simobj[0].solver.m_use_gpu && (simobj[0].solver.m_gpu_device_no >= 0)) {
      gpuSetDevice(simobj[0].solver.m_gpu_device_no);
  }
#endif

  if (!simobj[0].mpi.m_rank)  printf("Partitioning domain and allocating data arrays.\n");

  for (n = 0; n < nsims; n++) {

    /* this is a full initialization, not a barebones one */
    simobj[n].is_barebones = 0;

    /* allocations */
    simobj[n].mpi.m_ip           = (int*) calloc (simobj[n].solver.m_ndims,sizeof(int));
    simobj[n].mpi.m_is           = (int*) calloc (simobj[n].solver.m_ndims,sizeof(int));
    simobj[n].mpi.m_ie           = (int*) calloc (simobj[n].solver.m_ndims,sizeof(int));
    simobj[n].mpi.m_bcperiodic   = (int*) calloc (simobj[n].solver.m_ndims,sizeof(int));
    simobj[n].solver.m_dim_local = (int*) calloc (simobj[n].solver.m_ndims,sizeof(int));
    simobj[n].solver.m_is_periodic= (int*) calloc (simobj[n].solver.m_ndims,sizeof(int));

#if defined(HAVE_CUDA)
    simobj[n].mpi.m_wctime = 0;
    simobj[n].mpi.m_wctime_total = 0;
#endif

#ifndef serial
    _DECLARE_IERR_;

    /* Domain partitioning */
    int total_proc = 1;
    for (i=0; i<simobj[n].solver.m_ndims; i++) total_proc *= simobj[n].mpi.m_iproc[i];
    if (simobj[n].mpi.m_nproc != total_proc) {
      fprintf(stderr,"Error on rank %d: total number of processes is not consistent ", simobj[n].mpi.m_rank);
      fprintf(stderr,"with number of processes along each dimension.\n");
      if (nsims > 1) fprintf(stderr,"for domain %d.\n", n);
      fprintf(stderr,"mpiexec was called with %d processes, ",simobj[n].mpi.m_nproc);
      fprintf(stderr,"total number of processes from \"solver.inp\" is %d.\n", total_proc);
      return(1);
    }

    /* calculate ndims-D rank of each process (ip[]) from rank in MPI_COMM_WORLD */
    IERR MPIRanknD( simobj[n].solver.m_ndims,
                    simobj[n].mpi.m_rank,
                    simobj[n].mpi.m_iproc,
                    simobj[n].mpi.m_ip); CHECKERR(ierr);

    /* calculate local domain sizes along each dimension */
    for (i=0; i<simobj[n].solver.m_ndims; i++) {
      simobj[n].solver.m_dim_local[i] = MPIPartition1D( simobj[n].solver.m_dim_global[i],
                                                      simobj[n].mpi.m_iproc[i],
                                                      simobj[n].mpi.m_ip[i] );
    }

    /* calculate local domain limits in terms of global domain */
    IERR MPILocalDomainLimits(  simobj[n].solver.m_ndims,
                                simobj[n].mpi.m_rank,
                                &(simobj[n].mpi),
                                simobj[n].solver.m_dim_global,
                                simobj[n].mpi.m_is,
                                simobj[n].mpi.m_ie  );
    CHECKERR(ierr);

    /* create sub-communicators for parallel computations along grid lines in each dimension */
    IERR MPICreateCommunicators(simobj[n].solver.m_ndims,&(simobj[n].mpi)); CHECKERR(ierr);

    /* initialize periodic BC flags to zero */
    for (i=0; i<simobj[n].solver.m_ndims; i++) simobj[n].mpi.m_bcperiodic[i] = 0;

    /* create communication groups */
    IERR MPICreateIOGroups(&(simobj[n].mpi)); CHECKERR(ierr);

#else

    for (i=0; i<simobj[n].solver.m_ndims; i++) {
      simobj[n].mpi.m_ip[i]            = 0;
      simobj[n].solver.m_dim_local[i]  = simobj[n].solver.m_dim_global[i];
      simobj[n].mpi.m_iproc[i]         = 1;
      simobj[n].mpi.m_is[i]            = 0;
      simobj[n].mpi.m_ie[i]            = simobj[n].solver.m_dim_local[i];
      simobj[n].mpi.m_bcperiodic[i]    = 0;
    }

#endif

    simobj[n].solver.m_npoints_global
      = simobj[n].solver.m_npoints_local
      = simobj[n].solver.m_npoints_local_wghosts
      = 1;
    for (i=0; i<simobj[n].solver.m_ndims; i++) {
      simobj[n].solver.m_npoints_global *= simobj[n].solver.m_dim_global[i];
      simobj[n].solver.m_npoints_local *= simobj[n].solver.m_dim_local [i];
      simobj[n].solver.m_npoints_local_wghosts *= (simobj[n].solver.m_dim_local[i]+2*simobj[n].solver.m_ghosts);
    }

    /* Allocations */
    simobj[n].solver.m_index = (int*) calloc ((short)simobj[n].solver.m_ndims,sizeof(int));
    simobj[n].solver.m_stride_with_ghosts = (int*) calloc ((short)simobj[n].solver.m_ndims,sizeof(int));
    simobj[n].solver.m_stride_without_ghosts = (int*) calloc ((short)simobj[n].solver.m_ndims,sizeof(int));
    int accu1 = 1, accu2 = 1;
    for (i=0; i<simobj[n].solver.m_ndims; i++) {
      simobj[n].solver.m_stride_with_ghosts[i]    = accu1;
      simobj[n].solver.m_stride_without_ghosts[i] = accu2;
      accu1 *= (simobj[n].solver.m_dim_local[i]+2*simobj[n].solver.m_ghosts);
      accu2 *=  simobj[n].solver.m_dim_local[i];
    }

#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuMalloc((void**)&simobj[n].solver.m_gpu_dim_local, simobj[n].solver.m_ndims*sizeof(int));
      gpuMemcpy(  simobj[n].solver.m_gpu_dim_local,
                  simobj[n].solver.m_dim_local,
                  simobj[n].solver.m_ndims*sizeof(int),
                  gpuMemcpyHostToDevice );
    }
#endif

    /* state variables */
    int size = 1;
    for (i=0; i<simobj[n].solver.m_ndims; i++) {
      size *= (simobj[n].solver.m_dim_local[i]+2*simobj[n].solver.m_ghosts);
    }
    simobj[n].solver.m_ndof_cells_wghosts = simobj[n].solver.m_nvars*size;
    simobj[n].solver.m_u = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuMalloc((void**)&simobj[n].solver.m_gpu_u, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_gpu_u, 0, simobj[n].solver.m_nvars*size*sizeof(double));
    }
#endif
#ifdef with_petsc
    if (simobj[n].solver.m_use_petsc_ts) {
      simobj[n].solver.m_u0      = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_uref    = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_rhsref  = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_rhs     = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
    } else simobj[n].solver.m_u0 = simobj[n].solver.m_uref = simobj[n].solver.m_rhsref = simobj[n].solver.m_rhs = NULL;
#endif
#ifdef with_librom
    simobj[n].solver.m_u_rom_predicted = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
#endif
    simobj[n].solver.m_hyp     = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
    simobj[n].solver.m_par     = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
    simobj[n].solver.m_source  = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
    simobj[n].solver.m_iblank  = (double*) calloc (size              ,sizeof(double));

#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuMalloc((void**)&simobj[n].solver.m_hyp, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_par, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_source, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_hyp, 0, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_par, 0, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_source, 0, simobj[n].solver.m_nvars*size*sizeof(double));
    } else {
#endif
      simobj[n].solver.m_hyp     = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_par     = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_source  = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
#if defined(HAVE_CUDA)
    }
#endif

    simobj[n].solver.m_iblank  = (double*) calloc (size,sizeof(double));
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuMalloc((void**)&simobj[n].solver.m_gpu_iblank, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_gpu_iblank, 0, size*sizeof(double));
    }
#endif

    /* grid */
    size = 0;
    for (i=0; i<simobj[n].solver.m_ndims; i++) {
      size += (simobj[n].solver.m_dim_local[i]+2*simobj[n].solver.m_ghosts);
    }
    simobj[n].solver.m_x     = (double*) calloc (size,sizeof(double));
    simobj[n].solver.m_dxinv = (double*) calloc (size,sizeof(double));
    simobj[n].solver.m_size_x = size;
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuMalloc((void**)&simobj[n].solver.m_gpu_x, size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_gpu_dxinv, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_gpu_x, 0, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_gpu_dxinv, 0, size*sizeof(double));
    }
#endif

    /* cell-centered arrays needed to compute fluxes */
    size = 1;
    for (i=0; i<simobj[n].solver.m_ndims; i++) {
      size *= (simobj[n].solver.m_dim_local[i]+2*simobj[n].solver.m_ghosts);
    }
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuMalloc((void**)&simobj[n].solver.m_flux_c, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_u_c, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_deriv1, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_deriv2, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_flux_c, 0, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_u_c, 0, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_deriv1, 0, simobj[n].solver.m_nvars*size*sizeof(double));
      gpuMemset(simobj[n].solver.m_deriv2, 0, simobj[n].solver.m_nvars*size*sizeof(double));
    } else {
#endif
      simobj[n].solver.m_u_c     = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_flux_c  = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_deriv1 = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
      simobj[n].solver.m_deriv2 = (double*) calloc (simobj[n].solver.m_nvars*size,sizeof(double));
#if defined(HAVE_CUDA)
    }
#endif

    /* node-centered arrays needed to compute fluxes */
    size = 1;  for (i=0; i<simobj[n].solver.m_ndims; i++) size *= (simobj[n].solver.m_dim_local[i]+1);
    size *= simobj[n].solver.m_nvars;
    simobj[n].solver.m_ndof_nodes = size;
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuMalloc((void**)&simobj[n].solver.m_flux_i, size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_u_l, size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_u_r, size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_f_l, size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_f_r, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_flux_i, 0, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_u_l, 0, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_u_r, 0, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_f_l, 0, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_f_r, 0, size*sizeof(double));
    } else {
#endif
      simobj[n].solver.m_flux_i = (double*) calloc (size,sizeof(double));
      simobj[n].solver.m_u_l    = (double*) calloc (size,sizeof(double));
      simobj[n].solver.m_u_r    = (double*) calloc (size,sizeof(double));
      simobj[n].solver.m_f_l    = (double*) calloc (size,sizeof(double));
      simobj[n].solver.m_f_r    = (double*) calloc (size,sizeof(double));
#if defined(HAVE_CUDA)
    }
#endif

    /* allocate MPI send/receive buffer arrays */
    int bufdim[simobj[n].solver.m_ndims], maxbuf = 0;
    for (d = 0; d < simobj[n].solver.m_ndims; d++) {
      bufdim[d] = 1;
      for (i = 0; i < simobj[n].solver.m_ndims; i++) {
        if (i == d) bufdim[d] *= simobj[n].solver.m_ghosts;
        else        bufdim[d] *= simobj[n].solver.m_dim_local[i];
      }
      if (bufdim[d] > maxbuf) maxbuf = bufdim[d];
    }
    maxbuf *= (simobj[n].solver.m_nvars*simobj[n].solver.m_ndims);
    simobj[n].mpi.m_maxbuf  = maxbuf;
    simobj[n].mpi.m_sendbuf = (double*) calloc (2*simobj[n].solver.m_ndims*maxbuf,sizeof(double));
    simobj[n].mpi.m_recvbuf = (double*) calloc (2*simobj[n].solver.m_ndims*maxbuf,sizeof(double));
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      simobj[n].mpi.m_cpu_dim = (int *) calloc(simobj[n].solver.m_ndims, sizeof(int));
      _ArrayCopy1D_(simobj[n].solver.m_dim_local, simobj[n].mpi.m_cpu_dim, simobj[n].solver.m_ndims);
      gpuMalloc((void**)&simobj[n].mpi.m_gpu_sendbuf, 2*simobj[n].solver.m_ndims*simobj[n].mpi.m_maxbuf*sizeof(double));
      gpuMalloc((void**)&simobj[n].mpi.m_gpu_recvbuf, 2*simobj[n].solver.m_ndims*simobj[n].mpi.m_maxbuf*sizeof(double));
      gpuMemset(simobj[n].mpi.m_gpu_sendbuf, 0, 2*simobj[n].solver.m_ndims*simobj[n].mpi.m_maxbuf*sizeof(double));
      gpuMemset(simobj[n].mpi.m_gpu_recvbuf, 0, 2*simobj[n].solver.m_ndims*simobj[n].mpi.m_maxbuf*sizeof(double));
    }
#endif

    /* allocate the volume and boundary integral arrays */
    simobj[n].solver.m_volume_integral        = (double*) calloc (simobj[n].solver.m_nvars  ,sizeof(double));
    simobj[n].solver.m_volume_integral_initial = (double*) calloc (simobj[n].solver.m_nvars  ,sizeof(double));
    simobj[n].solver.m_total_boundary_integral = (double*) calloc (simobj[n].solver.m_nvars,sizeof(double));
    simobj[n].solver.m_conservation_error     = (double*) calloc (simobj[n].solver.m_nvars,sizeof(double));
    for (i=0; i<simobj[n].solver.m_nvars; i++) simobj[n].solver.m_conservation_error[i] = -1;
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      int total_offset = 0;
      for (d=0; d<simobj[n].solver.m_ndims; d++) {
          simobj[n].solver.m_gpu_npoints_boundary_offset[d] = total_offset;
          simobj[n].solver.m_gpu_npoints_boundary[d] = 1;

          for (i=0; i<simobj[n].solver.m_ndims; i++) {
              if (i != d) simobj[n].solver.m_gpu_npoints_boundary[d] *= simobj[n].solver.m_dim_local[i];
          }
          total_offset += 2*simobj[n].solver.m_gpu_npoints_boundary[d];
      }
      simobj[n].solver.m_stage_boundary_buffer_size = (total_offset*simobj[n].solver.m_nvars);
      gpuMalloc((void**)&simobj[n].solver.m_stage_boundary_buffer, simobj[n].solver.m_stage_boundary_buffer_size*sizeof(double));
      gpuMemset(simobj[n].solver.m_stage_boundary_buffer, 0, simobj[n].solver.m_stage_boundary_buffer_size*sizeof(double));

      size = 2*simobj[n].solver.m_ndims*simobj[n].solver.m_nvars;
      gpuMalloc((void**)&simobj[n].solver.m_stage_boundary_integral, size*sizeof(double));
      gpuMalloc((void**)&simobj[n].solver.m_step_boundary_integral, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_stage_boundary_integral, 0, size*sizeof(double));
      gpuMemset(simobj[n].solver.m_step_boundary_integral, 0, size*sizeof(double));
    } else {
#endif
      simobj[n].solver.m_stage_boundary_integral = (double*) calloc (2*simobj[n].solver.m_ndims*simobj[n].solver.m_nvars,sizeof(double));
      simobj[n].solver.m_step_boundary_integral  = (double*) calloc (2*simobj[n].solver.m_ndims*simobj[n].solver.m_nvars,sizeof(double));
#if defined(HAVE_CUDA)
    }
#endif

    /* initialize function call counts to zero */
    simobj[n].solver.m_count_hyp
      = simobj[n].solver.m_count_par
      = simobj[n].solver.m_count_sou
      = 0;
#ifdef with_petsc
    simobj[n].solver.m_count_rhs_function
      = simobj[n].solver.m_count_i_function
      = simobj[n].solver.m_count_i_jacobian
      = simobj[n].solver.m_count_i_jac_function
      = 0;
#endif

    /* Initialize iblank to 1*/
    _ArraySetValue_(simobj[n].solver.m_iblank,simobj[n].solver.m_npoints_local_wghosts,1);
#if defined(HAVE_CUDA)
    if (simobj[n].solver.m_use_gpu) {
      gpuArraySetValue(simobj[n].solver.m_gpu_iblank, simobj[n].solver.m_npoints_local_wghosts, 1.0);
    }
#endif

  }

  return 0;
}
