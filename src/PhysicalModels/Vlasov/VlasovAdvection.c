/*! @file VlasovAdvection.c
    @author Some Dumb Guy
    @brief Contains the function to compute the hyperbolic flux for the Vlasov equations over the domain.
*/

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/vlasov.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the Fourier transform sample frequency
    for bin.  Note that bin is 0-based.
    Calling this routine one bin at a time is not efficient,
    but we do not care about speed here.
*/
static int FFTFreqNum(int bin,  /*!< The bin number for the Fourier frequency to be returned */
                      int N     /*!< The total number of Fourier frequencies */
                     )
{
  double remainder;
  int last_positive = 0;

  remainder = N % 2;
  // Note that last_positive is an integer
  if (remainder) {
    // N is odd
    last_positive = (N-1) / 2;
  } else {
    // N is even
    last_positive = N/2 - 1;
  }

  if (bin <= last_positive) {
    return bin;
  } else {
    return -(N - bin);
  }
}

/*! Compute the hyperbolic flux over the local domain in the case where the electric field is prescribed.\n
    \f{equation}{
      E = 0.1 cos(x)
    \f}
*/
static int VlasovAdvectionPrescribed( double *f,   /*!< Array to hold the computed flux (same size and layout as u) */
                                      double *u,   /*!< Array containing the conserved solution */
                                      int     dir, /*!< Spatial dimension */
                                      void   *s,   /*!< Solver object of type #HyPar */
                                      double  t    /*!< Current time */
                                    )
{

  HyPar  *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->physics;

  int *dim    = solver->dim_local;
  int  ghosts = solver->ghosts;
  int  ndims  = solver->ndims;

  int index[ndims], bounds[ndims], bounds_noghost[ndims], offset[ndims];
  int i;

  // set bounds for array index to include ghost points
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  // set bounds for array index to NOT include ghost points
  _ArrayCopy1D_(dim,bounds_noghost,ndims);

  // set offset such that index is compatible with ghost point arrangement
  _ArraySetValue_(offset,ndims,-ghosts);

  // determine the velocity based on the direction
  if (dir == 0) {
    // the velocity coordinate is the velocity
    int done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);

      double velocity;
      _GetCoordinate_(1,index[1]-ghosts,dim,ghosts,solver->x,velocity);

      f[p] = velocity * u[p];
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  }
  else {
    // the electric field is the velocity
    int done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);

      double velocity;
      double x; _GetCoordinate_(0,index[0]-ghosts,dim,ghosts,solver->x,x);
      velocity = 0.1 * cos(x);

      f[p] = velocity * u[p];
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  }

  return(0);
}

/*! Compute the hyperbolic flux over the local domain in the case where the electric field is self-consistent.\n
    The field is solved from the solution values using a Poisson solve in Fourier space.
*/
static int VlasovAdvectionSelfConsistent( double *f,   /*!< Array to hold the computed flux (same size and layout as u) */
                                          double *u,   /*!< Array containing the conserved solution */
                                          int     dir, /*!< Spatial dimension */
                                          void   *s,   /*!< Solver object of type #HyPar */
                                          double  t    /*!< Current time */
                                        )
{
  HyPar        *solver = (HyPar*)         s;
  Vlasov       *param  = (Vlasov*)        solver->physics;
  MPIVariables *mpi    = (MPIVariables *) param->m;

#ifndef fftw
  if (!mpi->rank) {
    fprintf(stderr, "Error in VlasovAdvectionSelfConsistent():\n");
    fprintf(stderr, "  Using a self-consistent electric field requires FFTW.\n");
  }
  exit(1);
#endif

  int *dim    = solver->dim_local;
  int  N      = solver->dim_global[0];
  int  ghosts = solver->ghosts;
  int  ndims  = solver->ndims;

  double       *sum_buffer     = param->sum_buffer;
  double       *field          = param->field;
  fftw_complex *phys_buffer    = param->phys_buffer;
  fftw_complex *fourier_buffer = param->fourier_buffer;
  fftw_plan     plan_forward   = param->plan_forward;
  fftw_plan     plan_backward  = param->plan_backward;
  ptrdiff_t     local_ni       = param->local_ni;
  ptrdiff_t     local_i_start  = param->local_i_start;
  ptrdiff_t     local_no       = param->local_no;
  ptrdiff_t     local_o_start  = param->local_o_start;

  int index[ndims], bounds[ndims], bounds_noghost[ndims], offset[ndims];

  // set bounds for array index to include ghost points
  _ArrayCopy1D_(dim,bounds,ndims);
  for (int i = 0; i < ndims; i++) bounds[i] += 2*ghosts;

  // set bounds for array index to NOT include ghost points
  _ArrayCopy1D_(dim,bounds_noghost,ndims);

  // set offset such that index is compatible with ghost point arrangement
  _ArraySetValue_(offset,ndims,-ghosts);

  // Determine the velocity based on the direction
  if (dir == 0) {
    // The velocity coordinate is the HyPar "velocity" for this component
    int done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);

      double velocity;
      _GetCoordinate_(1,index[1]-ghosts,dim,ghosts,solver->x,velocity);

      f[p] = velocity * u[p];
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  }
  else {
    // The electric field is the HyPar "velocity" for this component.
    // Solve for the electric field

    // First, integrate the particle distribution over velocity.
    // Since the array dimension we want is not unit stride,
    // first manually add up the local part of the array.
    int done = 0; _ArraySetValue_(index,ndims,0);
    _ArraySetValue_(sum_buffer,dim[0],0);
    while (!done) {
      //int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
      int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);

      // accumulate f at this spatial location
      double dvinv; _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dvinv);
      double x; _GetCoordinate_(0,index[0],dim,ghosts,solver->x,x);
      double v; _GetCoordinate_(1,index[1],dim,ghosts,solver->x,v);

      sum_buffer[index[0]] += u[p] / dvinv;

      _ArrayIncrementIndex_(ndims,bounds_noghost,index,done);
    }

    // Now we can add up globally using MPI reduction 
    for (int i = 0; i < dim[0]; i++) {
      MPISum_double(&sum_buffer[i], &sum_buffer[i], 1, &mpi->comm[1]);
    }

    // Find the average velocity over all x
    double average_velocity = 0.0;
    for (int i = 0; i < dim[0]; i++) {
      average_velocity += sum_buffer[i];
    }
    MPISum_double(&average_velocity, &average_velocity, 1, &mpi->comm[0]);
    average_velocity /= (double) N;

    // Copy velocity-integrated values into complex-valued FFTW buffer
    for (int i = 0; i < dim[0]; i++) {
      phys_buffer[i][0] = sum_buffer[i] - average_velocity;
      phys_buffer[i][1] = 0.0;
    }

    // Execute the FFT
    MPI_Barrier(mpi->comm[0]);
    fftw_execute(plan_forward);
    MPI_Barrier(mpi->comm[0]);

    // Simultaneously do a Poisson solve and take derivative in frequency space
    int freq_start = 0;
    if (local_o_start == 0) {
      freq_start = 1;
      fourier_buffer[0][0] = 0.0;
      fourier_buffer[0][1] = 0.0;
    }

    for (int i = freq_start; i < local_no; i++) {
      double freq_num = FFTFreqNum(i + local_o_start, N);
      double thek = freq_num;
      double temp = fourier_buffer[i][0];
      // Swapping values is due to multiplication by i
      fourier_buffer[i][0] = - fourier_buffer[i][1] / thek;
      fourier_buffer[i][1] = temp / thek;
    }

    // Do an inverse Fourier transform to get back physical solved values
    MPI_Barrier(mpi->comm[0]);
    fftw_execute(plan_backward);
    MPI_Barrier(mpi->comm[0]);

    // copy the solved electric field into the velocity buffer
    for (int i = 0; i < dim[0]; i++) {
      field[i + ghosts] = - phys_buffer[i][0] / (double) N;
    }

    // Do halo exchange on the field
    MPIExchangeBoundaries1D(mpi, field, dim[0], ghosts, 0, 2);

    // copy the solved electric field out into the velocity buffer
    done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);

      f[p] = field[index[0]] * u[p];

      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  }

  return(0);
}

/*! Compute the hyperbolic flux over the local domain.\n
    \f{equation}{
      {\bf F}\left({\bf u}(x,v)\right) = c {\bf u}(x,v),
    \f}
    where
     \f{equation}{
      c = v,
    \f}
    if dir = 0, and
     \f{equation}{
      c = E,
    \f}
    the electric field if dir = 1.
*/
int VlasovAdvection( double  *f,   /*!< Array to hold the computed flux (same size and layout as u) */
                     double  *u,   /*!< Array containing the conserved solution */
                     int      dir, /*!< Spatial dimension */
                     void    *s,   /*!< Solver object of type #HyPar */
                     double   t    /*!< Current time */
                   )
{
  HyPar  *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->physics;

  bool self_consistent_electric_field = param->self_consistent_electric_field;

  if (self_consistent_electric_field) {
    return VlasovAdvectionSelfConsistent(f, u, dir, s, t);
  } else {
    return VlasovAdvectionPrescribed(f, u, dir, s, t);
  }
}
