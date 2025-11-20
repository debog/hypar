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
    if (!mpi->m_rank) {

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
      ferr = fscanf(in,"%d",&solver->m_n_boundary_zones); if (ferr != 1) return(1);
      boundary = (DomainBoundary*) calloc (solver->m_n_boundary_zones,sizeof(DomainBoundary));
      for (nb = 0; nb < solver->m_n_boundary_zones; nb++) {
        boundary[nb].m_DirichletValue = boundary[nb].m_SpongeValue
                                   = boundary[nb].m_FlowVelocity
                                   = boundary[nb].m_UnsteadyDirichletData
                                   = NULL;
        boundary[nb].m_UnsteadyDirichletSize = NULL;
      }

      /* read each boundary condition */
      for (nb = 0; nb < solver->m_n_boundary_zones; nb++) {
        int d, v;
        boundary[nb].m_xmin = (double*) calloc (solver->m_ndims,sizeof(double)); /* deallocated in BCCleanup.c */
        boundary[nb].m_xmax = (double*) calloc (solver->m_ndims,sizeof(double)); /* deallocated in BCCleanup.c */

        ferr = fscanf(in,"%s",boundary[nb].m_bctype); if (ferr != 1) return(1);
        ferr = fscanf(in,"%d",&boundary[nb].m_dim  ); if (ferr != 1) return(1);
        ferr = fscanf(in,"%d",&boundary[nb].m_face ); if (ferr != 1) return(1);
        for (d=0; d < solver->m_ndims; d++) {
          ferr = fscanf(in,"%lf %lf", &boundary[nb].m_xmin[d], &boundary[nb].m_xmax[d]);
          if (ferr != 2) return(1);
        }

        /* read in boundary type-specific additional data if required */

        if (!strcmp(boundary[nb].m_bctype,_DIRICHLET_)) {
          boundary[nb].m_DirichletValue = (double*) calloc (solver->m_nvars,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the Dirichlet value for each variable on this boundary */
          for (v = 0; v < solver->m_nvars; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_DirichletValue[v]);
        }

        if (!strcmp(boundary[nb].m_bctype,_SPONGE_)) {
          boundary[nb].m_SpongeValue = (double*) calloc (solver->m_nvars,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the sponge value for each variable on this boundary */
          for (v = 0; v < solver->m_nvars; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_SpongeValue[v]);
        }

        if (    (!strcmp(boundary[nb].m_bctype,_SLIP_WALL_))
            ||  (!strcmp(boundary[nb].m_bctype,_NOSLIP_WALL_)) ) {
          boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the wall velocity */
          for (v = 0; v < solver->m_ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_FlowVelocity[v]);
        }

        if (    (!strcmp(boundary[nb].m_bctype,_SW_SLIP_WALL_))
            ||  (!strcmp(boundary[nb].m_bctype,_SW_NOSLIP_WALL_)) ) {
          boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the wall velocity */
          for (v = 0; v < solver->m_ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_FlowVelocity[v]);
        }

        if (!strcmp(boundary[nb].m_bctype,_SUBSONIC_INFLOW_)) {
          boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density and velocity */
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowDensity);
          for (v = 0; v < solver->m_ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_FlowVelocity[v]);
        }

        if (!strcmp(boundary[nb].m_bctype,_SUBSONIC_OUTFLOW_)) {
          /* read in the outflow pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowPressure);
        }

        if (!strcmp(boundary[nb].m_bctype,_SUBSONIC_AMBIVALENT_)) {
          boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density, velocity, and pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowDensity);
          for (v = 0; v < solver->m_ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_FlowVelocity[v]);
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowPressure);
        }

        if (!strcmp(boundary[nb].m_bctype,_SUPERSONIC_INFLOW_)) {
          boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density, velocity and pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowDensity);
          for (v = 0; v < solver->m_ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_FlowVelocity[v]);
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowPressure);
        }

        if (!strcmp(boundary[nb].m_bctype,_TURBULENT_SUPERSONIC_INFLOW_)) {
          boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read in the inflow density, velocity and pressure */
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowDensity);
          for (v = 0; v < solver->m_ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_FlowVelocity[v]);
          ferr = fscanf(in,"%lf",&boundary[nb].m_FlowPressure);
          ferr = fscanf(in,"%s" , boundary[nb].m_UnsteadyDirichletFilename);
        }

        if (    (!strcmp(boundary[nb].m_bctype,_THERMAL_SLIP_WALL_))
            ||  (!strcmp(boundary[nb].m_bctype,_THERMAL_NOSLIP_WALL_)) ){
          boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
                                       /* deallocated in BCCleanup.c */
          /* read the wall velocity */
          for (v = 0; v < solver->m_ndims; v++) ferr = fscanf(in,"%lf",&boundary[nb].m_FlowVelocity[v]);
          /* read in the filename where temperature data is available */
          ferr = fscanf(in,"%s" , boundary[nb].m_UnsteadyTemperatureFilename);
        }

        /* if boundary is periodic, let the MPI and HyPar know */
        if (!strcmp(boundary[nb].m_bctype,_PERIODIC_)) {
          solver->m_is_periodic[boundary[nb].m_dim] = 1;
        }
        /*
          The MPI function to exchange internal (MPI) boundary information will handle
          periodic boundaries ONLY IF number of process along that dimension is more
          than 1.
        */
        if ((!strcmp(boundary[nb].m_bctype,_PERIODIC_)) && (mpi->m_iproc[boundary[nb].m_dim] > 1)) {
          mpi->m_bcperiodic[boundary[nb].m_dim] = 1;
        }

        /* some checks */
        if (boundary[nb].m_dim >= solver->m_ndims) {
          fprintf(stderr,"Error in reading boundary condition %d: dim %d is invalid (ndims = %d).\n",
                  nb,boundary[nb].m_dim,solver->m_ndims);
          return(1);
        }
        printf("  Boundary %30s:  Along dimension %2d and face %+1d\n",
                  boundary[nb].m_bctype,boundary[nb].m_dim,boundary[nb].m_face);
      }

      fclose(in);
      printf("%d boundary condition(s) read.\n",solver->m_n_boundary_zones);
    }

    /* tell other processes how many BCs are there and let them allocate */
    IERR MPIBroadcast_integer(&solver->m_n_boundary_zones,1,0,&mpi->m_world); CHECKERR(ierr);
    if (mpi->m_rank) {
      boundary = (DomainBoundary*) calloc (solver->m_n_boundary_zones,sizeof(DomainBoundary));
      for (nb = 0; nb < solver->m_n_boundary_zones; nb++) {
        boundary[nb].m_xmin = (double*) calloc (solver->m_ndims,sizeof(double)); /* deallocated in BCCleanup.c */
        boundary[nb].m_xmax = (double*) calloc (solver->m_ndims,sizeof(double)); /* deallocated in BCCleanup.c */
        boundary[nb].m_DirichletValue = boundary[nb].m_SpongeValue
                                   = boundary[nb].m_FlowVelocity
                                   = boundary[nb].m_UnsteadyDirichletData
                                   = NULL;
        boundary[nb].m_UnsteadyDirichletSize = NULL;
      }
    }

    /* communicate BC data to other processes */
    for (nb = 0; nb < solver->m_n_boundary_zones; nb++) {
      IERR MPIBroadcast_character(boundary[nb].m_bctype,_MAX_STRING_SIZE_,0,&mpi->m_world); CHECKERR(ierr);
      IERR MPIBroadcast_integer  (&boundary[nb].m_dim  ,1                ,0,&mpi->m_world); CHECKERR(ierr);
      IERR MPIBroadcast_integer  (&boundary[nb].m_face ,1                ,0,&mpi->m_world); CHECKERR(ierr);
      IERR MPIBroadcast_double   (boundary[nb].m_xmin  ,solver->m_ndims    ,0,&mpi->m_world); CHECKERR(ierr);
      IERR MPIBroadcast_double   (boundary[nb].m_xmax  ,solver->m_ndims    ,0,&mpi->m_world); CHECKERR(ierr);
    }
    IERR MPIBroadcast_integer(solver->m_is_periodic,solver->m_ndims,0,&mpi->m_world);CHECKERR(ierr);

    /* broadcast periodic boundary info for MPI to all processes */
    IERR MPIBroadcast_integer(mpi->m_bcperiodic,solver->m_ndims,0,&mpi->m_world);CHECKERR(ierr);

    /* On other processes, if necessary, allocate and receive boundary-type-specific data */
    for (nb = 0; nb < solver->m_n_boundary_zones; nb++) {
      if (!strcmp(boundary[nb].m_bctype,_DIRICHLET_)) {
        if (mpi->m_rank)  boundary[nb].m_DirichletValue = (double*) calloc (solver->m_nvars,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].m_DirichletValue,solver->m_nvars,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].m_bctype,_SPONGE_)) {
        if (mpi->m_rank)  boundary[nb].m_SpongeValue = (double*) calloc (solver->m_nvars,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].m_SpongeValue,solver->m_nvars,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (    (!strcmp(boundary[nb].m_bctype,_SLIP_WALL_))
          ||  (!strcmp(boundary[nb].m_bctype,_NOSLIP_WALL_)) ) {
        if (mpi->m_rank) boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].m_FlowVelocity,solver->m_ndims,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (    (!strcmp(boundary[nb].m_bctype,_SW_SLIP_WALL_))
          ||  (!strcmp(boundary[nb].m_bctype,_SW_NOSLIP_WALL_)) ) {
        if (mpi->m_rank) boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].m_FlowVelocity,solver->m_ndims,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].m_bctype,_SUBSONIC_INFLOW_)) {
        if (mpi->m_rank) boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].m_FlowDensity,1            ,0,&mpi->m_world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].m_FlowVelocity,solver->m_ndims,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].m_bctype,_SUBSONIC_OUTFLOW_)) {
        IERR MPIBroadcast_double(&boundary[nb].m_FlowPressure,1,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].m_bctype,_SUBSONIC_AMBIVALENT_)) {
        if (mpi->m_rank) boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].m_FlowDensity,1            ,0,&mpi->m_world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].m_FlowVelocity,solver->m_ndims,0,&mpi->m_world); CHECKERR(ierr);
        IERR MPIBroadcast_double(&boundary[nb].m_FlowPressure,1,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].m_bctype,_SUPERSONIC_INFLOW_)) {
        if (mpi->m_rank) boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].m_FlowDensity ,1            ,0,&mpi->m_world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].m_FlowVelocity ,solver->m_ndims,0,&mpi->m_world); CHECKERR(ierr);
        IERR MPIBroadcast_double(&boundary[nb].m_FlowPressure,1            ,0,&mpi->m_world); CHECKERR(ierr);
      }

      if (!strcmp(boundary[nb].m_bctype,_TURBULENT_SUPERSONIC_INFLOW_)) {
        if (mpi->m_rank) boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
        IERR MPIBroadcast_double(&boundary[nb].m_FlowDensity ,1            ,0,&mpi->m_world); CHECKERR(ierr);
        IERR MPIBroadcast_double(boundary[nb].m_FlowVelocity ,solver->m_ndims,0,&mpi->m_world); CHECKERR(ierr);
        IERR MPIBroadcast_double(&boundary[nb].m_FlowPressure,1            ,0,&mpi->m_world); CHECKERR(ierr);
        /* allocate arrays and read in unsteady boundary data */
        IERR BCReadTurbulentInflowData(&boundary[nb],mpi,solver->m_ndims,solver->m_nvars,solver->m_dim_local); CHECKERR(ierr);
      }

      if (    (!strcmp(boundary[nb].m_bctype,_THERMAL_SLIP_WALL_))
          ||  (!strcmp(boundary[nb].m_bctype,_THERMAL_NOSLIP_WALL_)) ) {
        if (mpi->m_rank) boundary[nb].m_FlowVelocity = (double*) calloc (solver->m_ndims,sizeof(double));
        IERR MPIBroadcast_double(boundary[nb].m_FlowVelocity,solver->m_ndims,0,&mpi->m_world); CHECKERR(ierr);
        /* allocate arrays and read in boundary temperature data */
        IERR BCReadTemperatureData(&boundary[nb],mpi,solver->m_ndims,solver->m_nvars,solver->m_dim_local); CHECKERR(ierr);
      }

    }

    solver->m_boundary = boundary;

    /* each process calculates its own part of these boundaries */
    IERR CalculateLocalExtent(solver,mpi); CHECKERR(ierr);

#if defined(HAVE_CUDA)
    int bounds[GPU_MAX_NDIMS];
    if (sim[0].solver.m_use_gpu) {
      for (nb = 0; nb < solver->m_n_boundary_zones; nb++) {
        _ArraySubtract1D_(bounds,boundary[nb].m_ie,boundary[nb].m_is,solver->m_ndims);

        _ArrayProduct1D_(bounds,solver->m_ndims,boundary[nb].m_gpu_npoints_bounds);
        boundary[nb].m_gpu_npoints_local_wghosts = solver->m_npoints_local_wghosts;

        _ArrayProduct1D_(bounds,solver->m_ndims,boundary[nb].m_gpu_npoints_bounds);
        boundary[nb].m_gpu_npoints_local_wghosts = solver->m_npoints_local_wghosts;

        gpuMalloc((void**)&boundary[nb].m_gpu_is, solver->m_ndims*sizeof(int));
        gpuMalloc((void**)&boundary[nb].m_gpu_ie, solver->m_ndims*sizeof(int));
        gpuMalloc((void**)&boundary[nb].m_gpu_bounds, solver->m_ndims*sizeof(int));
        gpuMemcpy(boundary[nb].m_gpu_is, boundary[nb].m_is, solver->m_ndims*sizeof(int), gpuMemcpyHostToDevice);
        gpuMemcpy(boundary[nb].m_gpu_ie, boundary[nb].m_ie, solver->m_ndims*sizeof(int), gpuMemcpyHostToDevice);
        gpuMemcpy(boundary[nb].m_gpu_bounds, bounds, solver->m_ndims*sizeof(int), gpuMemcpyHostToDevice);
        if (   (!strcmp(boundary[nb].m_bctype,_SLIP_WALL_))
            || (!strcmp(boundary[nb].m_bctype,_NOSLIP_WALL_)) ) {
            gpuMalloc((void**)&boundary[nb].m_gpu_FlowVelocity, solver->m_ndims*sizeof(double));
            gpuMemcpy(  boundary[nb].m_gpu_FlowVelocity,
                        boundary[nb].m_FlowVelocity,
                        solver->m_ndims*sizeof(double),
                        gpuMemcpyHostToDevice);
        }
      }
    }
#endif

    /* initialize function pointers for each boundary condition */
    for (nb = 0; nb < solver->m_n_boundary_zones; nb++) {
#if defined(HAVE_CUDA)
      BCInitialize(&boundary[nb], solver->m_use_gpu);
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
  DomainBoundary  *boundary = solver->m_boundary;

  int n;
  for (n = 0; n < solver->m_n_boundary_zones; n++) {
    /* allocate */
    boundary[n].m_is   = (int*) calloc (solver->m_ndims,sizeof(int)); /* deallocated in BCCleanup.c */
    boundary[n].m_ie   = (int*) calloc (solver->m_ndims,sizeof(int)); /* deallocated in BCCleanup.c */

    int d,dim = boundary[n].m_dim;

    if (!strcmp(boundary[n].m_bctype,_SPONGE_)) {
      /* Sponge boundary condition */
      boundary[n].m_on_this_proc = 1;
      int offset = 0;
      for (d=0; d<solver->m_ndims; d++) {
        int is, ie;
        FindInterval(boundary[n].m_xmin[d],boundary[n].m_xmax[d],
                     &solver->m_x[offset+solver->m_ghosts],
                     solver->m_dim_local[d],&is,&ie);
        boundary[n].m_is[d] = is;
        boundary[n].m_ie[d] = ie;
        if ((ie-is) <= 0) boundary[n].m_on_this_proc = 0;
        offset += solver->m_dim_local[d] + 2*solver->m_ghosts;
      }
    } else {
      /* other boundary conditions */
      if (boundary[n].m_face == 1) {

        if (mpi->m_ip[dim] == 0) {
          boundary[n].m_on_this_proc = 1;
          int offset = 0;
          for (d=0; d<solver->m_ndims; d++) {
            if (d == dim) {
              boundary[n].m_is[d] = -solver->m_ghosts;
              boundary[n].m_ie[d] = 0;
            } else {
              int is, ie;
              FindInterval(boundary[n].m_xmin[d],boundary[n].m_xmax[d],
                           &solver->m_x[offset+solver->m_ghosts],
                           solver->m_dim_local[d],&is,&ie);
              boundary[n].m_is[d] = is;
              boundary[n].m_ie[d] = ie;
              if ((ie-is) <= 0) boundary[n].m_on_this_proc = 0;
            }
            offset += solver->m_dim_local[d] + 2*solver->m_ghosts;
          }
        } else  boundary[n].m_on_this_proc = 0;

      } else if (boundary[n].m_face == -1) {

        if (mpi->m_ip[dim] == mpi->m_iproc[dim]-1) {
          boundary[n].m_on_this_proc = 1;
          int offset = 0;
          for (d=0; d<solver->m_ndims; d++) {
            if (d == dim) {
              boundary[n].m_is[d] = solver->m_dim_local[dim];
              boundary[n].m_ie[d] = solver->m_dim_local[dim] + solver->m_ghosts;
            } else {
              int is, ie;
              FindInterval(boundary[n].m_xmin[d],boundary[n].m_xmax[d],
                           &solver->m_x[offset+solver->m_ghosts],
                           solver->m_dim_local[d],&is,&ie);
              boundary[n].m_is[d] = is;
              boundary[n].m_ie[d] = ie;
              if ((ie-is) <= 0) boundary[n].m_on_this_proc = 0;
            }
            offset += solver->m_dim_local[d] + 2*solver->m_ghosts;
          }
        } else  boundary[n].m_on_this_proc = 0;
      }
    }

  }

  return(0);
}
