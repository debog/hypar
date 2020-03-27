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
#include <fftw3-mpi.h>
#endif

/*! Compute the maximum CFL number */
double VlasovComputeCFL (void*,void*,double,double);
/*! Compute the advection term */
int VlasovAdvection (double*,double*,int,void*,double);
/*! Compute the upwind flux at interfaces */
int VlasovUpwind (double*,double*,double*,double*,
                  double*,double*,int,void*,double);

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
  int          *dim       = solver->dim_global;
  int          *dim_local = solver->dim_local;
  int           ghosts    = solver->ghosts;
  int           i, ferr;

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
  physics->self_consistent_electric_field = false;

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {
    FILE *in;
    in = fopen("physics.inp","r");
    if (in) {
      printf("Reading physical model inputs from file \"physics.inp\".\n");
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "self_consistent_electric_field")) {
            /* read whether electric field is self-consistent or prescribed */
            ferr = fscanf(in,"%d", &physics->self_consistent_electric_field);
            if (ferr != 1) return(1);
          } else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
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

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in VlasovInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

#ifndef serial
  /* Broadcast parsed problem data */
  IERR MPIBroadcast_integer((int *) &physics->self_consistent_electric_field,
                            1,0,&mpi->world); CHECKERR(ierr);
#endif

  if (physics->self_consistent_electric_field) {
#ifdef fftw
    /* If using FFTW, make sure MPI is enabled
       Currently we only support the distrubuted memory version
    */
#ifdef serial
    if (!mpi->rank) {
      fprintf(stderr,"Error in VlasovInitialize: Using FFTW requires MPI to be enabled.\n");
    }
    return(1);
#endif

    /* Put the mpi object in the params for access in other functions */
    physics->m = m;
  
    /* Create a scratch buffer for moving between real and complex values */
    physics->sum_buffer = (double*) calloc(dim_local[0], sizeof(double));
  
    /* Create a buffer to hold the electric field and do halo exchange */
    physics->field = (double*) calloc(dim_local[0] + 2*ghosts, sizeof(double));
  
    /* Initialize FFTW and set up data buffers used for the transforms */
    fftw_mpi_init();
    physics->alloc_local = fftw_mpi_local_size_1d(dim[0], mpi->comm[0],
                                                  FFTW_FORWARD, 0,
                                                  &physics->local_ni,
                                                  &physics->local_i_start,
                                                  &physics->local_no,
                                                  &physics->local_o_start);
    if (dim_local[0] != physics->local_ni) {
      fprintf(stderr,"Error in VlasovInitialize:  The FFTW data distribution is incompatible with the HyPar one.\n");
      fprintf(stderr,"Decompose the spatial dimension so that the degrees of freedom are evenly divided.\n");
      return(1);
    }
  
    physics->phys_buffer = fftw_alloc_complex(physics->alloc_local);
    physics->fourier_buffer = fftw_alloc_complex(physics->alloc_local);
  
    physics->plan_forward = fftw_mpi_plan_dft_1d(dim[0],
                                                 physics->phys_buffer,
                                                 physics->fourier_buffer,
                                                 mpi->comm[0],
                                                 FFTW_FORWARD,
                                                 FFTW_ESTIMATE);
  
    physics->plan_backward = fftw_mpi_plan_dft_1d(dim[0],
                                                  physics->fourier_buffer,
                                                  physics->phys_buffer,
                                                  mpi->comm[0],
                                                  FFTW_BACKWARD,
                                                  FFTW_ESTIMATE);
#else
  if (!mpi->rank) {
    fprintf(stderr,"Error in VlasovInitialize():\n");
    fprintf(stderr,"  Self-consistent electric field requires FFTW library.\n");
  }
  return(1);
#endif
  }

  /* initializing physical model-specific functions */
  solver->ComputeCFL = VlasovComputeCFL;
  solver->FFunction  = VlasovAdvection;
  solver->Upwind     = VlasovUpwind;

  return(0);
}
