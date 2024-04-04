/*! @file VlasovInitialize.c
    @author John Loffeld
    @brief Initialize the Vlasov module
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <bandedmatrix.h>
#include <physicalmodels/vlasov.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef fftw
#ifdef serial
#include <fftw3.h>
#else
#include <fftw3-mpi.h>
#endif
#endif

/*! Compute the maximum CFL number */
double VlasovComputeCFL (void*,void*,double,double);
/*! Compute the advection term */
int VlasovAdvection (double*,double*,int,void*,double);
/*! Compute the upwind flux at interfaces */
int VlasovUpwind (double*,double*,double*,double*,
                  double*,double*,int,void*,double);
/*! Write E-field and potential to file */
int VlasovWriteEFieldAndPotential(void*, void*,double);

int VlasovPreStep(double*,void*,void*,double);
int VlasovPostStage(double*,void*,void*,double);

int VlasovEField(double*, void*, double);

/*! Initialize the Vlasov physics module -
    allocate and set physics-related parameters, read physics-related inputs
    from file, and set the physics-related function pointers in #HyPar
*/
int VlasovInitialize(void *s, /*!< Solver object of type #HyPar */
                     void *m  /*!< Object of type #MPIVariables containing MPI-related info */
                    )
{
  HyPar        *solver    = (HyPar*)        s;
  MPIVariables *mpi       = (MPIVariables*) m;
  Vlasov       *physics   = (Vlasov*)       solver->physics;

  int *dim_global = solver->dim_global;
  int *dim_local = solver->dim_local;
  int ghosts = solver->ghosts;

  if (solver->nvars != _MODEL_NVARS_) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in VlasovInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    }
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in VlasovInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    }
    return(1);
  }

  /* default is prescribed electric field */
  physics->self_consistent_electric_field = 0;
  physics->ndims_x = 1;
  physics->ndims_v = 1;

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {
    FILE *in;
    in = fopen("physics.inp","r");
    if (in) {
      printf("Reading physical model inputs from file \"physics.inp\".\n");
      char word[_MAX_STRING_SIZE_];
      int ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          int ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "self_consistent_electric_field")) {
            /* read whether electric field is self-consistent or prescribed */
            int ferr = fscanf(in,"%d", &physics->self_consistent_electric_field);
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "x_ndims")) {
            /* read number of spatial dimensions */
            int ferr = fscanf(in,"%d", &physics->ndims_x);
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "v_ndims")) {
            /* read number of velocity dimensions */
            int ferr = fscanf(in,"%d", &physics->ndims_v);
            if (ferr != 1) return(1);
          } else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            int ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",
                   word, useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
      } else {
        fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
        return(1);
      }
    }
    fclose(in);
  }

  if ((physics->ndims_x+physics->ndims_v) != solver->ndims) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in VlasovInitialize:\n");
      fprintf(stderr, "  space + vel dims not equal to ndims!\n");
    }
    return(1);
  }

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in VlasovInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

#ifndef serial
  /* Broadcast parsed problem data */
  MPIBroadcast_integer(&physics->ndims_x,1,0,&mpi->world);
  MPIBroadcast_integer(&physics->ndims_v,1,0,&mpi->world);
  MPIBroadcast_integer((int *) &physics->self_consistent_electric_field,
                       1,0,&mpi->world);
#endif

  /* compute local number of x-space points with ghosts */
  physics->npts_local_x_wghosts = 1;
  physics->npts_local_x = 1;
  physics->npts_global_x_wghosts = 1;
  physics->npts_global_x = 1;
  for (int d=0; d<physics->ndims_x; d++) {
    physics->npts_local_x_wghosts *= (dim_local[d]+2*ghosts);
    physics->npts_local_x *= dim_local[d];
    physics->npts_global_x_wghosts *= (dim_global[d]+2*ghosts);
    physics->npts_global_x *= dim_global[d];
  }

  /* allocate array to hold the electric field (needs to have ghosts);
     note that number of electric field components is the number of
     spatial dimensions */
  physics->e_field = (double*) calloc(  physics->npts_local_x_wghosts
                                      * physics->ndims_x,
                                        sizeof(double)  );
  physics->potential = (double*) calloc(  physics->npts_local_x_wghosts
                                      * physics->ndims_x,
                                        sizeof(double)  );

  /* Put the mpi object in the params for access in other functions */
  physics->m = m;

  if (physics->self_consistent_electric_field) {
#ifdef fftw
    /* If using FFTW, make sure MPI is enabled
       Currently we only support the distrubuted memory version
    */
#ifdef serial
    fprintf(stderr,"Error in VlasovInitialize(): Using FFTW requires MPI to be enabled.\n");
    return(1);
#else

    if (physics->ndims_x > 1) {
      fprintf(stderr,"Error in VlasovInitialize():\n");
      fprintf(stderr,"  Self-consistent electric field is implemented for only 1 space dimension.\n");
      return 1;
    }

    /* Create a scratch buffer for moving between real and complex values */
    physics->sum_buffer = (double*) calloc(dim_local[0], sizeof(double));

    /* Initialize FFTW and set up data buffers used for the transforms */
    fftw_mpi_init();
    physics->alloc_local = fftw_mpi_local_size_1d(dim_global[0], mpi->comm[0],
                                                  FFTW_FORWARD, 0,
                                                  &physics->local_ni,
                                                  &physics->local_i_start,
                                                  &physics->local_no,
                                                  &physics->local_o_start);
    if (dim_local[0] != physics->local_ni) {
      fprintf(stderr,"Error in VlasovInitialize(): The FFTW data distribution is incompatible with the HyPar one.\n");
      fprintf(stderr,"Decompose the spatial dimension so that the degrees of freedom are evenly divided.\n");
      return(1);
    }

    physics->phys_buffer_e = fftw_alloc_complex(physics->alloc_local);
    physics->fourier_buffer_e = fftw_alloc_complex(physics->alloc_local);

    physics->plan_forward_e = fftw_mpi_plan_dft_1d(dim_global[0],
                                                 physics->phys_buffer_e,
                                                 physics->fourier_buffer_e,
                                                 mpi->comm[0],
                                                 FFTW_FORWARD,
                                                 FFTW_ESTIMATE);

    physics->plan_backward_e = fftw_mpi_plan_dft_1d(dim_global[0],
                                                  physics->fourier_buffer_e,
                                                  physics->phys_buffer_e,
                                                  mpi->comm[0],
                                                  FFTW_BACKWARD,
                                                  FFTW_ESTIMATE);

    physics->phys_buffer_phi = fftw_alloc_complex(physics->alloc_local);
    physics->fourier_buffer_phi = fftw_alloc_complex(physics->alloc_local);

    physics->plan_forward_phi = fftw_mpi_plan_dft_1d(dim_global[0],
                                                 physics->phys_buffer_phi,
                                                 physics->fourier_buffer_phi,
                                                 mpi->comm[0],
                                                 FFTW_FORWARD,
                                                 FFTW_ESTIMATE);

    physics->plan_backward_phi = fftw_mpi_plan_dft_1d(dim_global[0],
                                                  physics->fourier_buffer_phi,
                                                  physics->phys_buffer_phi,
                                                  mpi->comm[0],
                                                  FFTW_BACKWARD,
                                                  FFTW_ESTIMATE);
#endif
#else
  fprintf(stderr,"Error in VlasovInitialize():\n");
  fprintf(stderr,"  Self-consistent electric field requires FFTW library.\n");
  return(1);
#endif
  }

  /* initializing physical model-specific functions */
  solver->PreStep       = VlasovPreStep;
  solver->ComputeCFL    = VlasovComputeCFL;
  solver->FFunction     = VlasovAdvection;
  solver->Upwind        = VlasovUpwind;
  solver->PhysicsOutput = VlasovWriteEFieldAndPotential;
  solver->PostStage     = VlasovPostStage;

  int ierr = VlasovEField(solver->u, solver, 0.0);
  if (ierr) return ierr;

  return 0;
}
