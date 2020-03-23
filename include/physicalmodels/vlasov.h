#include <stdbool.h>

/*
  Vlasov Equation

Reference:
+ Henon, "Vlasov equation?", Astronomy and Astrophysics, 114, 1982

  df         df                     df 
  --  +  v ------  +  (E + v x B) ------  =  0
  dt        dx_i                   dv_i

*/

#ifdef fftw
#include <fftw3-mpi.h>
#endif

#define _VLASOV_  "vlasov"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
/*! Number of spatial dimensions */
#define _MODEL_NDIMS_ 2
/*! Number of variables per grid point */
#define _MODEL_NVARS_ 1


typedef struct vlasov_parameters {
  bool self_consistent_electric_field;

#ifdef fftw
  void  *m;
  double *sum_buffer;
  double *field;
  fftw_plan plan_forward;
  fftw_plan plan_backward;
  fftw_complex *phys_buffer;
  fftw_complex *fourier_buffer;
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni, local_i_start;
  ptrdiff_t local_no, local_o_start;
#endif
} Vlasov;

int    VlasovInitialize (void*,void*);
int    VlasovCleanup    (void*);
