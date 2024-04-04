/*! @file VlasovEField.c
    @author John Loffeld, Ping-Hsuan Tsai
    @brief Contains the function to compute the electric field
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

/*! Sets the prescribed electric field. */
static int SetEFieldPrescribed( double* u,/*!< Conserved solution */
                                void*   s,/*!< Solver object of type #HyPar */
                                double  t /*!< Current time */
                              )
{

  HyPar *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  int* dim_local = solver->dim_local;
  int ghosts = solver->ghosts;

  int ndims_x = param->ndims_x;

  int dim_x[ndims_x];
  _ArrayCopy1D_(dim_local,dim_x,ndims_x);

  int bounds[ndims_x];
  _ArrayCopy1D_(dim_x,bounds,ndims_x);
  for (int d=0; d<ndims_x; d++) bounds[d] += 2*ghosts;

  int offset[ndims_x];
  _ArraySetValue_(offset,ndims_x,-ghosts);

  int done = 0;
  int index[ndims_x]; _ArraySetValue_(index, ndims_x, 0);

  while (!done) {

    double xvec[ndims_x];
    for (int d=0; d<ndims_x; d++) {
      _GetCoordinate_(d,index[d]-ghosts,dim_local,ghosts,solver->x,xvec[d]);
    }

    int p;
    _ArrayIndex1DWO_(ndims_x,dim_x,index,offset,ghosts,p);
    param->e_field[p*ndims_x+0] = 0.1 * cos(xvec[0]);

    _ArrayIncrementIndex_(ndims_x,bounds,index,done);
  }

  return 0;
}

/*! Compute the self-consistent electric field over the local domain: The field
 * is solved from the solution values using a Poisson solve in Fourier space. */
static int SetEFieldSelfConsistent(double* u,/*!< Conserved solution */
                                   void*   s,/*!< Solver object of type #HyPar */
                                   double  t /*!< Current time */
                                        )
{
  HyPar  *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  if (param->ndims_x > 1) {
    fprintf(stderr,"Error in SetEFieldSelfConsistent():\n");
    fprintf(stderr,"  Implemented for 1 spatial dimension only.\n");
    return 1;
  }

#ifndef fftw

  fprintf(stderr,"Error in SetEFieldSelfConsistent():\n");
  fprintf(stderr,"  Using a self-consistent electric field requires FFTW.\n");
  exit(1);

#else

  int *dim    = solver->dim_local;
  int  N      = solver->dim_global[0];
  int  ghosts = solver->ghosts;
  int  ndims  = solver->ndims;

  double       *sum_buffer       = param->sum_buffer;
  double       *field            = param->e_field;
  fftw_complex *phys_buffer_e    = param->phys_buffer_e;
  fftw_complex *fourier_buffer_e = param->fourier_buffer_e;
  fftw_plan     plan_forward_e   = param->plan_forward_e;
  fftw_plan     plan_backward_e  = param->plan_backward_e;
  ptrdiff_t     local_ni         = param->local_ni;
  ptrdiff_t     local_i_start    = param->local_i_start;
  ptrdiff_t     local_no         = param->local_no;
  ptrdiff_t     local_o_start    = param->local_o_start;

  fftw_complex *phys_buffer_phi    = param->phys_buffer_phi;
  fftw_complex *fourier_buffer_phi = param->fourier_buffer_phi;
  fftw_plan     plan_forward_phi   = param->plan_forward_phi;
  fftw_plan     plan_backward_phi  = param->plan_backward_phi;

  int index[ndims], bounds[ndims], bounds_noghost[ndims], offset[ndims];

  // set bounds for array index to include ghost points
  _ArrayCopy1D_(dim,bounds,ndims);
  for (int i = 0; i < ndims; i++) bounds[i] += 2*ghosts;

  // set bounds for array index to NOT include ghost points
  _ArrayCopy1D_(dim,bounds_noghost,ndims);

  // set offset such that index is compatible with ghost point arrangement
  _ArraySetValue_(offset,ndims,-ghosts);

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

  // Find the average density over all x
  double average_velocity = 0.0;
  for (int i = 0; i < dim[0]; i++) {
    average_velocity += sum_buffer[i];
  }
  MPISum_double(&average_velocity, &average_velocity, 1, &mpi->comm[0]);
  average_velocity /= (double) N;

  // Copy velocity-integrated values into complex-valued FFTW buffer
  for (int i = 0; i < dim[0]; i++) {
    phys_buffer_e[i][0] = sum_buffer[i] - average_velocity;
    phys_buffer_e[i][1] = 0.0;
  }

  // Execute the FFT
  MPI_Barrier(mpi->comm[0]);
  fftw_execute(plan_forward_e);
  MPI_Barrier(mpi->comm[0]);

  // Simultaneously do a Poisson solve and take derivative in frequency space
  int freq_start = 0;
  if (local_o_start == 0) {
    freq_start = 1;
    fourier_buffer_e[0][0] = 0.0;
    fourier_buffer_e[0][1] = 0.0;
  }

  // Simultaneously do a Poisson solve in frequency space
  int freq_start_phi = 0;
  if (local_o_start == 0) {
    freq_start_phi = 1;
    fourier_buffer_phi[0][0] = 0.0;
    fourier_buffer_phi[0][1] = 0.0;
  }

  for (int i = freq_start_phi; i < local_no; i++) {
    double freq_num = FFTFreqNum(i + local_o_start, N);
    double thek = freq_num;
    // Swapping values is due to multiplication by i
    fourier_buffer_phi[i][0] = fourier_buffer_e[i][0] / (thek*thek);
    fourier_buffer_phi[i][1] = fourier_buffer_e[i][1]/ (thek*thek);
  }


  for (int i = freq_start; i < local_no; i++) {
    double freq_num = FFTFreqNum(i + local_o_start, N);
    double thek = freq_num;
    double temp = fourier_buffer_e[i][0];
    // Swapping values is due to multiplication by i
    fourier_buffer_e[i][0] = - fourier_buffer_e[i][1] / thek;
    fourier_buffer_e[i][1] = temp / thek;
  }

  // Do an inverse Fourier transform to get back physical solved values
  MPI_Barrier(mpi->comm[0]);
  fftw_execute(plan_backward_e);
  MPI_Barrier(mpi->comm[0]);

  // Do an inverse Fourier transform to get back physical solved values
  MPI_Barrier(mpi->comm[0]);
  fftw_execute(plan_backward_phi);
  MPI_Barrier(mpi->comm[0]);

  // copy the solved electric field into the e buffer
  for (int i = 0; i < dim[0]; i++) {
    param->e_field[i + ghosts] = - phys_buffer_e[i][0] / (double) N;
  }

  // copy the solved potential field into the potential buffer
  for (int i = 0; i < dim[0]; i++) {
    param->potential[i + ghosts] = phys_buffer_phi[i][0] / (double) N;
  }

  // Do halo exchange on the e
  MPIExchangeBoundaries1D(mpi, param->e_field, dim[0], ghosts, 0, ndims);

  // Do halo exchange on the potential
  MPIExchangeBoundaries1D(mpi, param->potential, dim[0], ghosts, 0, ndims);

#endif

  return 0;
}

/*! Compute electric field.\n
*/
int VlasovEField( double* u, /*!< Conserved solution */
                  void*   s, /*!< Solver object of type #HyPar */
                  double  t  /*!< Current time */
                )
{
  HyPar  *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->physics;

  int ierr;

  if (param->self_consistent_electric_field) {
    ierr = SetEFieldSelfConsistent(u, s, t);
  } else {
    ierr = SetEFieldPrescribed(u, s, t);
  }

  return ierr;
}

