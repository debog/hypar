/*! @file Initialize.c
    @author Debojyoti Ghosh
    @brief Initialization function
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
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

  for (n = 0; n < nsims; n++) {

    /* this is a full initialization, not a barebones one */
    simobj[n].is_barebones = 0;

    /* allocations */
    simobj[n].mpi.ip           = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    simobj[n].mpi.is           = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    simobj[n].mpi.ie           = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    simobj[n].mpi.bcperiodic   = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    simobj[n].solver.dim_local = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    simobj[n].solver.isPeriodic= (int*) calloc (simobj[n].solver.ndims,sizeof(int));

#ifndef serial
    _DECLARE_IERR_;
  
    /* Domain partitioning */
    if (!simobj[n].mpi.rank) {
      if (nsims == 1) printf("Partitioning domain.\n");
      else            printf("Partitioning domain %d.\n", n);
    }
    int total_proc = 1;
    for (i=0; i<simobj[n].solver.ndims; i++) total_proc *= simobj[n].mpi.iproc[i];
    if (simobj[n].mpi.nproc != total_proc) {
      fprintf(stderr,"Error on rank %d: total number of processes is not consistent ", simobj[n].mpi.rank);
      fprintf(stderr,"with number of processes along each dimension.\n");
      if (nsims > 1) fprintf(stderr,"for domain %d.\n", n);
      fprintf(stderr,"mpiexec was called with %d processes, ",simobj[n].mpi.nproc);
      fprintf(stderr,"total number of processes from \"solver.inp\" is %d.\n", total_proc);
      return(1);
    }
  
    /* calculate ndims-D rank of each process (ip[]) from rank in MPI_COMM_WORLD */
    IERR MPIRanknD( simobj[n].solver.ndims,
                    simobj[n].mpi.rank,
                    simobj[n].mpi.iproc,
                    simobj[n].mpi.ip); CHECKERR(ierr);
  
    /* calculate local domain sizes along each dimension */
    for (i=0; i<simobj[n].solver.ndims; i++) {
      simobj[n].solver.dim_local[i] = MPIPartition1D( simobj[n].solver.dim_global[i],
                                                      simobj[n].mpi.iproc[i],
                                                      simobj[n].mpi.ip[i] );
    }
  
    /* calculate local domain limits in terms of global domain */
    IERR MPILocalDomainLimits(  simobj[n].solver.ndims,
                                simobj[n].mpi.rank,
                                &(simobj[n].mpi),
                                simobj[n].solver.dim_global,
                                simobj[n].mpi.is,
                                simobj[n].mpi.ie  );
    CHECKERR(ierr);
  
    /* create sub-communicators for parallel computations along grid lines in each dimension */
    IERR MPICreateCommunicators(simobj[n].solver.ndims,&(simobj[n].mpi)); CHECKERR(ierr);
  
    /* initialize periodic BC flags to zero */
    for (i=0; i<simobj[n].solver.ndims; i++) simobj[n].mpi.bcperiodic[i] = 0;
  
    /* create communication groups */
    IERR MPICreateIOGroups(&(simobj[n].mpi)); CHECKERR(ierr);

#else

    for (i=0; i<simobj[n].solver.ndims; i++) {
      simobj[n].mpi.ip[i]            = 0;
      simobj[n].solver.dim_local[i]  = simobj[n].solver.dim_global[i];
      simobj[n].mpi.iproc[i]         = 1;
      simobj[n].mpi.is[i]            = 0;
      simobj[n].mpi.ie[i]            = simobj[n].solver.dim_local[i];
      simobj[n].mpi.bcperiodic[i]    = 0;
    }

#endif

    simobj[n].solver.npoints_global 
      = simobj[n].solver.npoints_local 
      = simobj[n].solver.npoints_local_wghosts 
      = 1;
    for (i=0; i<simobj[n].solver.ndims; i++) {
      simobj[n].solver.npoints_global *= simobj[n].solver.dim_global[i];
      simobj[n].solver.npoints_local *= simobj[n].solver.dim_local [i];
      simobj[n].solver.npoints_local_wghosts *= (simobj[n].solver.dim_local[i]+2*simobj[n].solver.ghosts);
    }
  
    /* Allocations */
    if (!simobj[n].mpi.rank) printf("Allocating data arrays.\n");
    simobj[n].solver.index = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    simobj[n].solver.stride_with_ghosts = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    simobj[n].solver.stride_without_ghosts = (int*) calloc (simobj[n].solver.ndims,sizeof(int));
    int accu1 = 1, accu2 = 1;
    for (i=0; i<simobj[n].solver.ndims; i++) {
      simobj[n].solver.stride_with_ghosts[i]    = accu1;
      simobj[n].solver.stride_without_ghosts[i] = accu2;
      accu1 *= (simobj[n].solver.dim_local[i]+2*simobj[n].solver.ghosts);
      accu2 *=  simobj[n].solver.dim_local[i];
    }

    int size;

    /* state variables */
    size = 1;
    for (i=0; i<simobj[n].solver.ndims; i++) {
      size *= (simobj[n].solver.dim_local[i]+2*simobj[n].solver.ghosts);
    }
    simobj[n].solver.u = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
#ifdef with_petsc
    if (simobj[n].solver.use_petscTS) {
      simobj[n].solver.u0      = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
      simobj[n].solver.uref    = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
      simobj[n].solver.rhsref  = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
      simobj[n].solver.rhs     = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    } else simobj[n].solver.u0 = simobj[n].solver.uref = simobj[n].solver.rhsref = simobj[n].solver.rhs = NULL;
#endif
    simobj[n].solver.hyp     = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.par     = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.source  = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.iblank  = (double*) calloc (size              ,sizeof(double));
    /* grid */
    size = 0;
    for (i=0; i<simobj[n].solver.ndims; i++) size += (simobj[n].solver.dim_local[i]+2*simobj[n].solver.ghosts);
    simobj[n].solver.x     = (double*) calloc (size,sizeof(double));
    simobj[n].solver.dxinv = (double*) calloc (size,sizeof(double));
    /* arrays needed to compute fluxes */
    size = 1;  for (i=0; i<simobj[n].solver.ndims; i++) size *= (simobj[n].solver.dim_local[i]+2*simobj[n].solver.ghosts);
    simobj[n].solver.uC     = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.fluxC  = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.Deriv1 = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.Deriv2 = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    size = 1;  for (i=0; i<simobj[n].solver.ndims; i++) size *= (simobj[n].solver.dim_local[i]+1);
    simobj[n].solver.fluxI = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.uL    = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.uR    = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.fL    = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    simobj[n].solver.fR    = (double*) calloc (simobj[n].solver.nvars*size,sizeof(double));
    /* allocate MPI send/receive buffer arrays */
    int bufdim[simobj[n].solver.ndims], maxbuf = 0;
    for (d = 0; d < simobj[n].solver.ndims; d++) {
      bufdim[d] = 1;
      for (i = 0; i < simobj[n].solver.ndims; i++) {
        if (i == d) bufdim[d] *= simobj[n].solver.ghosts;
        else        bufdim[d] *= simobj[n].solver.dim_local[i];
      }
      if (bufdim[d] > maxbuf) maxbuf = bufdim[d];
    }
    maxbuf *= simobj[n].solver.nvars;
    simobj[n].mpi.maxbuf  = maxbuf;
    simobj[n].mpi.sendbuf = (double*) calloc (2*simobj[n].solver.ndims*maxbuf,sizeof(double));
    simobj[n].mpi.recvbuf = (double*) calloc (2*simobj[n].solver.ndims*maxbuf,sizeof(double));
    /* allocate the volume and boundary integral arrays */
    simobj[n].solver.VolumeIntegral        = (double*) calloc (simobj[n].solver.nvars  ,sizeof(double));
    simobj[n].solver.VolumeIntegralInitial = (double*) calloc (simobj[n].solver.nvars  ,sizeof(double));
    simobj[n].solver.StageBoundaryIntegral = (double*) calloc (2*simobj[n].solver.ndims*simobj[n].solver.nvars,sizeof(double));
    simobj[n].solver.StepBoundaryIntegral  = (double*) calloc (2*simobj[n].solver.ndims*simobj[n].solver.nvars,sizeof(double));
    simobj[n].solver.TotalBoundaryIntegral = (double*) calloc (simobj[n].solver.nvars,sizeof(double));
    simobj[n].solver.ConservationError     = (double*) calloc (simobj[n].solver.nvars,sizeof(double));

    /* initialize function call counts to zero */
    simobj[n].solver.count_hyp 
      = simobj[n].solver.count_par 
      = simobj[n].solver.count_sou 
      = 0;
#ifdef with_petsc
    simobj[n].solver.count_RHSFunction 
      = simobj[n].solver.count_IFunction
      = simobj[n].solver.count_IJacobian 
      = simobj[n].solver.count_IJacFunction 
      = 0;
#endif

    /* Initialize iblank to 1*/
    _ArraySetValue_(simobj[n].solver.iblank,simobj[n].solver.npoints_local_wghosts,1);

  }

  return(0);
}
