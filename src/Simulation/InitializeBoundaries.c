/*! @file InitializeBoundaries.c
    @author Debojyoti Ghosh, Youngdae Kim
    @brief Initialize the boundary implementation
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#if defined(HAVE_CUDA)
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#endif
#include <mathfunctions.h>
#include <boundaryconditions.h>
#include <mpivars.h>
#include <simulation_object.h>

static int CalculateLocalExtent(void*,void*);

/*! This function initializes the variables and functions related to implementing
    the boundary conditions.
    + Rank 0 reads in the boundary conditions file and broadcasts the information
      to all processors.
    + Depending on the type of boundary, additional information is read in. For
      example, for Dirichlet boundary, the Dirichlet value is read in.
    + Allocate and initialize arrays and variables related to implementing the
      boundary conditions.
    + Each rank finds out if the subdomain it owns abuts any of the boundaries
      specified.

    Note that boundary conditions are implemented as boundary objects of the
    type #DomainBoundary.
*/
int InitializeBoundaries( void  *s,   /*!< Array of simulation objects of type #SimulationObject */
                          int   nsims /*!< number of simulation objects */
                        )
{
  SimulationObject *sim = (SimulationObject*) s;
  int ns;
  _DECLARE_IERR_;

  for (ns = 0; ns < nsims; ns++) {

    DomainBoundary  *boundary = NULL;
    HyPar           *solver   = &(sim[ns].solver);
    MPIVariables    *mpi      = &(sim[ns].mpi);
    int nb, ferr;

    /* root process reads boundary condition file */
    if (!mpi->rank) {

      char filename[_MAX_STRING_SIZE_] = "boundary";
      char filename_backup[_MAX_STRING_SIZE_] = "boundary";
      if (nsims > 1) {
        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(ns, index, (int)log10(nsims)+1);
        strcat(filename, "_");
        strcat(filename, index);
      }
      strcat(filename, ".inp");
      strcat(filename_backup, ".inp");

      FILE *in;
      in = fopen(filename,"r");
      if (!in) {
        in = fopen(filename_backup, "r");
        if (!in) {
          fprintf(stderr,"Error: boundary condition file %s or %s not found.\n",
                  filename, filename_backup );
          return(1);
        } else {
          if (nsims > 1) printf("Domain %d: ", ns);
          printf("Reading boundary conditions from %s.\n", filename_backup);
        }
      } else {
        if (nsims > 1) printf("Domain %d: ", ns);
        printf("Reading boundary conditions from %s.\n", filename);
      }

      /* read number of boundary conditions and allocate */
      ferr = fscanf(in,"%d",&solver->nBoundaryZones); if (ferr != 1) return(1);
      boundary = (DomainBoundary*) calloc (solver->nBoundaryZones,sizeof(DomainBoundary));
      for (nb = 0; nb < solver->nBoundaryZones; nb++) {
        boundary[nb].DirichletValue = boundary[nb].SpongeValue
                                   = boundary[nb].FlowVelocity
                                   = boundary[nb].UnsteadyDirichletData
                                   = NULL;
        boundary[nb].UnsteadyDirichletSize = NULL;
      }

      /* read each boundary condition */
      for (nb = 0; nb < solver->nBoundaryZones; nb++) {
        int d, v;
        boundary[nb].xmin = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
        boundary[nb].xmax = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */

        ferr = fscanf(in,"%s",boundary[nb].bctype); if (ferr != 1) return(1);
        ferr = fscanf(in,"%d",&boundary[nb].dim  ); if (ferr != 1) return(1);
        ferr = fscanf(in,"%d",&boundary[nb].face ); if (ferr != 1) return(1);
        for (d=0; d < solver->ndims; d++) {
          ferr = fscanf(in,"%lf %lf", &boundary[nb].xmin[d], &boundary[nb].xmax[d]);
          if (ferr != 2) return(1);
        }

        /* read in boundary type-specific additional data if required */

        if (!strcmp(boundary[nb].bctype,_DIRICHLET_)) {
          boundary[nb].DirichletValue = (double*) calloc (solver->nvars,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the Dirichlet value for each variable on this boundary */
          for (v = 0; v < solver->nvars; v++) ferr = fscanf(in,"%lf",&boundary[nb].DirichletValue[v]);
        }

        if (!strcmp(boundary[nb].bctype,_SPONGE_)) {
          boundary[nb].SpongeValue = (double*) calloc (solver->nvars,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the sponge value for each variable on this boundary */
          for (v = 0; v < solver->nvars; v++) ferr = fscanf(in,"%lf",&boundary[nb].SpongeValue[v]);
        }

        if (    (!strcmp(boundary[nb].bctype,_SLIP_WALL_))
            ||  (!strcmp(boundary[nb].bctype,_NOSLIP_WALL_)) ) {
          boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the wall velocity */
          for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].FlowVelocity[v]);
        }

        if (    (!strcmp(boundary[nb].bctype,_SW_SLIP_WALL_))
            ||  (!strcmp(boundary[nb].bctype,_SW_NOSLIP_WALL_)) ) {
          boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the wall velocity */
          for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].FlowVelocity[v]);
        }

        if (!strcmp(boundary[nb].bctype,_SUBSONIC_INFLOW_)) {
          boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density and velocity */
          ferr = fscanf(in,"%lf",&boundary[nb].FlowDensity);
          for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].FlowVelocity[v]);
        }

        if (!strcmp(boundary[nb].bctype,_SUBSONIC_OUTFLOW_)) {
          /* read in the outflow pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].FlowPressure);
        }

        if (!strcmp(boundary[nb].bctype,_SUBSONIC_AMBIVALENT_)) {
          boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density, velocity, and pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].FlowDensity);
          for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].FlowVelocity[v]);
          ferr = fscanf(in,"%lf",&boundary[nb].FlowPressure);
        }

        if (!strcmp(boundary[nb].bctype,_SUPERSONIC_INFLOW_)) {
          boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density, velocity and pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].FlowDensity);
          for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].FlowVelocity[v]);
          ferr = fscanf(in,"%lf",&boundary[nb].FlowPressure);
        }

        if (!strcmp(boundary[nb].bctype,_TURBULENT_SUPERSONIC_INFLOW_)) {
          boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density, velocity and pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].FlowDensity);
          for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].FlowVelocity[v]);
          ferr = fscanf(in,"%lf",&boundary[nb].FlowPressure);
          ferr = fscanf(in,"%s" , boundary[nb].UnsteadyDirichletFilename);
        }

        if (    (!strcmp(boundary[nb].bctype,_THERMAL_SLIP_WALL_))
            ||  (!strcmp(boundary[nb].bctype,_THERMAL_NOSLIP_WALL_)) ){
          boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the wall velocity */
          for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].FlowVelocity[v]);
          /* read in the filename where temperature data is available */
          ferr = fscanf(in,"%s" , boundary[nb].UnsteadyTemperatureFilename);
        }

        /* if boundary is periodic, let the MPI and HyPar know */
        if (!strcmp(boundary[nb].bctype,_PERIODIC_)) {
          solver->isPeriodic[boundary[nb].dim] = 1;
        }
        /*
          The MPI function to exchange internal (MPI) boundary information will handle
          periodic boundaries ONLY IF number of process along that dimension is more
          than 1.
        */
        if ((!strcmp(boundary[nb].bctype,_PERIODIC_)) && (mpi->iproc[boundary[nb].dim] > 1)) {
          mpi->bcperiodic[boundary[nb].dim] = 1;
        }

        /* some checks */
        if (boundary[nb].dim >= solver->ndims) {
          fprintf(stderr,"Error in reading boundary condition %d: dim %d is invalid (ndims = %d).\n",
                  nb,boundary[nb].dim,solver->ndims);
          return(1);
        }
        printf("  Boundary %30s:  Along dimension %2d and face %+1d\n",
                  boundary[nb].bctype,boundary[nb].dim,boundary[nb].face);
      }

      fclose(in);
      printf("%d boundary condition(s) read.\n",solver->nBoundaryZones);
    }

    /* tell other processes how many BCs are there and let them allocate */
    IERR MPIBroadcast_integer(&solver->nBoundaryZones,1,0,&mpi->world); CHECKERR(ierr);
    if (mpi->rank) {
      boundary = (DomainBoundary*) calloc (solver->nBoundaryZones,sizeof(DomainBoundary));
      for (nb = 0; nb < solver->nBoundaryZones; nb++) {
        boundary[nb].xmin = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
        boundary[nb].xmax = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
        boundary[nb].DirichletValue = boundary[nb].SpongeValue
                                   = boundary[nb].FlowVelocity
                                   = boundary[nb].UnsteadyDirichletData
                                   = NULL;
        boundary[nb].UnsteadyDirichletSize = NULL;
      }
    }

    /* communicate BC data to other processes */
    for (nb = 0; nb < solver->nBoundaryZones; nb++) {
      IERR MPIBroadcast_character(boundary[nb].bctype,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_integer  (&boundary[nb].dim  ,1                ,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_integer  (&boundary[nb].face ,1                ,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double   (boundary[nb].xmin  ,solver->ndims    ,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double   (boundary[nb].xmax  ,solver->ndims    ,0,&mpi->world); CHECKERR(ierr);
    }
    IERR MPIBroadcast_integer(solver->isPeriodic,solver->ndims,0,&mpi->world);CHECKERR(ierr);

    /* broadcast periodic boundary info for MPI to all processes */
    IERR MPIBroadcast_integer(mpi->bcperiodic,solver->ndims,0,&mpi->world);CHECKERR(ierr);

    /* On other processes, if necessary, allocate and receive boundary-type-specific data */
    for (nb = 0; nb < solver->nBoundaryZones; nb++) {
      if (!strcmp(boundary[nb].bctype,_DIRICHLET_)) {
        if (mpi->rank)  boundary[nb].DirichletValue = (double*) calloc (solver->nvars,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].DirichletValue,solver->nvars,0,&mpi->world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].bctype,_SPONGE_)) {
        if (mpi->rank)  boundary[nb].SpongeValue = (double*) calloc (solver->nvars,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].SpongeValue,solver->nvars,0,&mpi->world); CHECKERR(ierr);
      }

      if (    (!strcmp(boundary[nb].bctype,_SLIP_WALL_))
          ||  (!strcmp(boundary[nb].bctype,_NOSLIP_WALL_)) ) {
        if (mpi->rank) boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
      }

      if (    (!strcmp(boundary[nb].bctype,_SW_SLIP_WALL_))
          ||  (!strcmp(boundary[nb].bctype,_SW_NOSLIP_WALL_)) ) {
        if (mpi->rank) boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].bctype,_SUBSONIC_INFLOW_)) {
        if (mpi->rank) boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].FlowDensity,1            ,0,&mpi->world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].bctype,_SUBSONIC_OUTFLOW_)) {
        IERR MPIBroadcast_double(&boundary[nb].FlowPressure,1,0,&mpi->world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].bctype,_SUBSONIC_AMBIVALENT_)) {
        if (mpi->rank) boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].FlowDensity,1            ,0,&mpi->world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
        IERR MPIBroadcast_double(&boundary[nb].FlowPressure,1,0,&mpi->world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].bctype,_SUPERSONIC_INFLOW_)) {
        if (mpi->rank) boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].FlowDensity ,1            ,0,&mpi->world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].FlowVelocity ,solver->ndims,0,&mpi->world); CHECKERR(ierr);
        IERR MPIBroadcast_double(&boundary[nb].FlowPressure,1            ,0,&mpi->world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].bctype,_TURBULENT_SUPERSONIC_INFLOW_)) {
        if (mpi->rank) boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].FlowDensity ,1            ,0,&mpi->world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].FlowVelocity ,solver->ndims,0,&mpi->world); CHECKERR(ierr);
        IERR MPIBroadcast_double(&boundary[nb].FlowPressure,1            ,0,&mpi->world); CHECKERR(ierr);
        /* allocate arrays and read in unsteady boundary data */
        IERR BCReadTurbulentInflowData(&boundary[nb],mpi,solver->ndims,solver->nvars,solver->dim_local); CHECKERR(ierr);
      }

      if (    (!strcmp(boundary[nb].bctype,_THERMAL_SLIP_WALL_))
          ||  (!strcmp(boundary[nb].bctype,_THERMAL_NOSLIP_WALL_)) ) {
        if (mpi->rank) boundary[nb].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
        /* allocate arrays and read in boundary temperature data */
        IERR BCReadTemperatureData(&boundary[nb],mpi,solver->ndims,solver->nvars,solver->dim_local); CHECKERR(ierr);
      }

    }

    solver->boundary = boundary;

    /* each process calculates its own part of these boundaries */
    IERR CalculateLocalExtent(solver,mpi); CHECKERR(ierr);

#if defined(HAVE_CUDA)
    int bounds[GPU_MAX_NDIMS];
    if (sim[0].solver.use_gpu) {
      for (nb = 0; nb < solver->nBoundaryZones; nb++) {
        _ArraySubtract1D_(bounds,boundary[nb].ie,boundary[nb].is,solver->ndims);

        _ArrayProduct1D_(bounds,solver->ndims,boundary[nb].gpu_npoints_bounds);
        boundary[nb].gpu_npoints_local_wghosts = solver->npoints_local_wghosts;

        _ArrayProduct1D_(bounds,solver->ndims,boundary[nb].gpu_npoints_bounds);
        boundary[nb].gpu_npoints_local_wghosts = solver->npoints_local_wghosts;

        gpuMalloc((void**)&boundary[nb].gpu_is, solver->ndims*sizeof(int));
        gpuMalloc((void**)&boundary[nb].gpu_ie, solver->ndims*sizeof(int));
        gpuMalloc((void**)&boundary[nb].gpu_bounds, solver->ndims*sizeof(int));
        gpuMemcpy(boundary[nb].gpu_is, boundary[nb].is, solver->ndims*sizeof(int), gpuMemcpyHostToDevice);
        gpuMemcpy(boundary[nb].gpu_ie, boundary[nb].ie, solver->ndims*sizeof(int), gpuMemcpyHostToDevice);
        gpuMemcpy(boundary[nb].gpu_bounds, bounds, solver->ndims*sizeof(int), gpuMemcpyHostToDevice);
        if (   (!strcmp(boundary[nb].bctype,_SLIP_WALL_))
            || (!strcmp(boundary[nb].bctype,_NOSLIP_WALL_)) ) {
            gpuMalloc((void**)&boundary[nb].gpu_FlowVelocity, solver->ndims*sizeof(double));
            gpuMemcpy(  boundary[nb].gpu_FlowVelocity,
                        boundary[nb].FlowVelocity,
                        solver->ndims*sizeof(double),
                        gpuMemcpyHostToDevice);
        }
      }
    }
#endif

    /* initialize function pointers for each boundary condition */
    for (nb = 0; nb < solver->nBoundaryZones; nb++) {
#if defined(HAVE_CUDA)
      BCInitialize(&boundary[nb], solver->use_gpu);
#else
      BCInitialize(&boundary[nb], 0);
#endif
    }

  }

  return 0;
}

/*! For each of the boundary conditions, compute its extent on the
    local sub-domain of each rank (if at all this subdomain abuts
    that boundary), and accordingly set the bounding grid indices.
*/
int CalculateLocalExtent(
                          void *s, /*!< Solver object of type #HyPar */
                          void *m  /*!< MPI object of type #MPIVariables */
                        )
{
  HyPar           *solver   = (HyPar*)        s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  DomainBoundary  *boundary = solver->boundary;

  int n;
  for (n = 0; n < solver->nBoundaryZones; n++) {
    /* allocate */
    boundary[n].is   = (int*) calloc (solver->ndims,sizeof(int)); /* deallocated in BCCleanup.c */
    boundary[n].ie   = (int*) calloc (solver->ndims,sizeof(int)); /* deallocated in BCCleanup.c */

    int d,dim = boundary[n].dim;

    if (!strcmp(boundary[n].bctype,_SPONGE_)) {
      /* Sponge boundary condition */
      boundary[n].on_this_proc = 1;
      int offset = 0;
      for (d=0; d<solver->ndims; d++) {
        int is, ie;
        FindInterval(boundary[n].xmin[d],boundary[n].xmax[d],
                     &solver->x[offset+solver->ghosts],
                     solver->dim_local[d],&is,&ie);
        boundary[n].is[d] = is;
        boundary[n].ie[d] = ie;
        if ((ie-is) <= 0) boundary[n].on_this_proc = 0;
        offset += solver->dim_local[d] + 2*solver->ghosts;
      }
    } else {
      /* other boundary conditions */
      if (boundary[n].face == 1) {

        if (mpi->ip[dim] == 0) {
          boundary[n].on_this_proc = 1;
          int offset = 0;
          for (d=0; d<solver->ndims; d++) {
            if (d == dim) {
              boundary[n].is[d] = -solver->ghosts;
              boundary[n].ie[d] = 0;
            } else {
              int is, ie;
              FindInterval(boundary[n].xmin[d],boundary[n].xmax[d],
                           &solver->x[offset+solver->ghosts],
                           solver->dim_local[d],&is,&ie);
              boundary[n].is[d] = is;
              boundary[n].ie[d] = ie;
              if ((ie-is) <= 0) boundary[n].on_this_proc = 0;
            }
            offset += solver->dim_local[d] + 2*solver->ghosts;
          }
        } else  boundary[n].on_this_proc = 0;

      } else if (boundary[n].face == -1) {

        if (mpi->ip[dim] == mpi->iproc[dim]-1) {
          boundary[n].on_this_proc = 1;
          int offset = 0;
          for (d=0; d<solver->ndims; d++) {
            if (d == dim) {
              boundary[n].is[d] = solver->dim_local[dim];
              boundary[n].ie[d] = solver->dim_local[dim] + solver->ghosts;
            } else {
              int is, ie;
              FindInterval(boundary[n].xmin[d],boundary[n].xmax[d],
                           &solver->x[offset+solver->ghosts],
                           solver->dim_local[d],&is,&ie);
              boundary[n].is[d] = is;
              boundary[n].ie[d] = ie;
              if ((ie-is) <= 0) boundary[n].on_this_proc = 0;
            }
            offset += solver->dim_local[d] + 2*solver->ghosts;
          }
        } else  boundary[n].on_this_proc = 0;
      }
    }

  }

  return(0);
}
