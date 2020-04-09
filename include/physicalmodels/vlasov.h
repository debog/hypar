/*! @file vlasov.h
    @brief Vlasov Equation
    @author John Loffeld

  1D-1V Vlasov Equation:
  \f{equation}{
      \frac{\partial f}{\partial t} 
      +
      v \frac{\partial f}{\partial x}
      +
      E \frac{\partial f}{\partial v}
      = 0,
  \f}
  where 
    + \f$f\f$ is the distribution function
    + \f$x\f$ is the spatial coordinate
    + \f$v\f$ is the velocity
    + \f$E\f$ is the electric field

  Reference:
    + Henon, "Vlasov equation?", Astronomy and Astrophysics, 114, 1982
*/

#include <stdbool.h>

#ifdef fftw
#include <fftw3-mpi.h>
#endif

/*! \def _VLASOV_ 
    1D-1V Vlasov equation
*/
#define _VLASOV_  "vlasov"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 2
/*! Number of variables per grid point */
#define _MODEL_NVARS_ 1


/*! \def Vlasov
    \brief Structure containing variables and parameters for the Vlasov equation
 *  This structure contains the physical parameters, variables, and function pointers 
 *  specific to the Vlasov equation.
*/
/* \brief Structure containing variables and parameters for the Vlasov equation
 * This structure contains the physical parameters, variables, and function pointers 
 * specific to the Vlasov equation.
*/
typedef struct vlasov_parameters {

  /*! Use a self-consistent electric field? Requires FFTW library */
  bool self_consistent_electric_field;

  /*! Number of spatial dimensions */
  int ndims_x;

  /*! Number of velocity dimensions */
  int ndims_v;

  /*! Number of spatial grid points (local) */
  int npts_local_x;

  /*! Number of spatial grid points (global) */
  long npts_global_x;

  /*! Number of spatial grid points with ghosts (local) */
  int npts_local_x_wghosts;

  /*! Number of spatial grid points with ghosts (global) */
  long npts_global_x_wghosts;

  /*! electric field */
  double *e_field;

  /*! Pointer to MPI object of type #MPIVariables */
  void  *m;

#ifdef fftw

  /*! Buffer sum */
  double *sum_buffer;
  /*! Forward FFT plan */
  fftw_plan plan_forward;
  /*! Backward FFT plan */
  fftw_plan plan_backward;
  /*! buffer */
  fftw_complex *phys_buffer;
  /*! buffer */
  fftw_complex *fourier_buffer;
  /*! */
  ptrdiff_t alloc_local;
  /*! */
  ptrdiff_t local_ni, local_i_start;
  /*! */
  ptrdiff_t local_no, local_o_start;
#endif

} Vlasov;

/*! Initialize the Vlasov physics module */
int VlasovInitialize (void*,void*);
/*! Clean up the Vlasov physics module */
int VlasovCleanup    (void*);
