/*

  This code generates the initial solution for the
  isotropic turbulence decay problem, using the technique
  described in the 1981 NASA report by Rogallo.

  It needs the FFTW3 library installed. To compile it,
  make sure the fftw3.h and libfftw3 are available in
  the include and linking paths.

  The initial solution is written for the local (MPI) sub-domains,
  and the number of files written depend on the number of I/O ranks.
  Although this code is serial, it will compute the decomposition
  of the global domain and generate the initial solution accordingly.

  + "input_mode" in solver.inp must be set to "parallel n" where n
    is the number of I/O ranks.
  + "ip_file_type" in solver.inp must be set to "binary".

  Note: this code does *not* involve the allocation of the global domain,

*/

#define _MAX_STRING_SIZE_ 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define _ArrayIncrementIndex_(N,imax,i,done) \
  { \
    int arraycounter = 0; \
    while (arraycounter < (N)) { \
      if (i[arraycounter] == imax[arraycounter]-1) { \
        i[arraycounter] = 0; \
        arraycounter++; \
      } else { \
        i[arraycounter]++; \
        break; \
      } \
    } \
    if (arraycounter == (N)) done = 1; \
    else          done = 0; \
  }

#define _ArrayIndex1D_(N,imax,i,ghost,index)  \
  { \
    index = i[N-1]+(ghost); \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost))); \
    } \
  }

static const double PI = 4 * atan(1.0);

double raiseto(double x, double a) {
  return exp(a*log(x));
}

double magnitude(double a, double b) {
  return sqrt(a*a + b*b);
}

void velocityComponent( const int N,
                        const int dir,
                        const double kp,
                        const double u0,
                        double* const uvel )
{
  long N3 = N*N*N;
  long i,j,k;

  double kk = sqrt(3 * (N/2)*(N/2));
  int kkmax = (int) kk;

  fftw_complex *uhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
  fftw_complex *u = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));

  fftw_plan inv_trans_u;
  inv_trans_u = fftw_plan_dft_3d(N, N, N, uhat, u, 1, FFTW_MEASURE);

  /* Specifying the velocities in Fourier space */
  for (i=0; i<N3; i++) uhat[i][0] = uhat[i][1] = 0.0;

  for (i = 1; i < N/2; i++){
    for (j = 0; j < N/2; j++){
      for (k = 0; k < N/2; k++){
        double kk   = sqrt(i*i + j*j + k*k);
        double th1  = 2*PI * ((double)rand())/((double)RAND_MAX);
        double th2  = 2*PI * ((double)rand())/((double)RAND_MAX);
        double phi1 = 2*PI * ((double)rand())/((double)RAND_MAX);

        double E = 16.0 * sqrt(2.0/PI) * (u0*u0/kp) * raiseto(kk/kp, 4.0)
                        * exp(-2.0*(kk/kp)*(kk/kp));
        double alfa_real = sqrt(E/(4*PI*kk*kk))*cos(th1)*cos(phi1);
        double alfa_imag = sqrt(E/(4*PI*kk*kk))*sin(th1)*cos(phi1);
        double beta_real = sqrt(E/(4*PI*kk*kk))*cos(th2)*sin(phi1);
        double beta_imag = sqrt(E/(4*PI*kk*kk))*sin(th2)*sin(phi1);

        if (dir == 0) {
          uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
          uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));
        } else if (dir == 1) {
          uhat[k+N*(j+N*i)][0] = (beta_real*j*k-alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
          uhat[k+N*(j+N*i)][1] = (beta_imag*j*k-alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));
        } else {
          uhat[k+N*(j+N*i)][0] = -(beta_real*sqrt(i*i+j*j))/kk;
          uhat[k+N*(j+N*i)][1] = -(beta_imag*sqrt(i*i+j*j))/kk;
        }

      }
    }
  }
  for (i = 0; i < 1; i++){
    for (k = 0; k < N/2; k++){
      for (j = 1; j < N/2; j++){
        double kk   = sqrt(i*i + j*j + k*k);
        double th1  = 2*PI * ((double)rand())/((double)RAND_MAX);
        double th2  = 2*PI * ((double)rand())/((double)RAND_MAX);
        double phi1 = 2*PI * ((double)rand())/((double)RAND_MAX);

        double E = 16.0 * sqrt(2.0/PI) * (u0*u0/kp) * raiseto(kk/kp, 4.0)
                        * exp(-2.0*(kk/kp)*(kk/kp));
        double alfa_real = sqrt(E/(4*PI*kk*kk))*cos(th1)*cos(phi1);
        double alfa_imag = sqrt(E/(4*PI*kk*kk))*sin(th1)*cos(phi1);
        double beta_real = sqrt(E/(4*PI*kk*kk))*cos(th2)*sin(phi1);
        double beta_imag = sqrt(E/(4*PI*kk*kk))*sin(th2)*sin(phi1);

        if (dir == 0) {
          uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
          uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));
        } else if (dir == 1) {
          uhat[k+N*(j+N*i)][0] = (beta_real*j*k-alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
          uhat[k+N*(j+N*i)][1] = (beta_imag*j*k-alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));
        } else {
          uhat[k+N*(j+N*i)][0] = -(beta_real*sqrt(i*i+j*j))/kk;
          uhat[k+N*(j+N*i)][1] = -(beta_imag*sqrt(i*i+j*j))/kk;
        }

      }
    }
  }
  for (i = 0; i < 1; i++){
    for (j = 0; j < 1; j++){
      for (k = 0; k < N/2; k++){
        uhat[k+N*(j+N*i)][0] = 0;
        uhat[k+N*(j+N*i)][1] = 0;
      }
    }
  }

  /* The following is necessary to ensure that the inverse Fourier
     transform yields a real velocity field with zero imaginary
     components */

  for (i=N/2+1; i<N; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=N/2+1; k<N; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*((N-j)+N*(N-i))][0];
        uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*((N-j)+N*(N-i))][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=0; k<1; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[k+N*((N-j)+N*(N-i))][0];
        uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*(N-i))][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=0; j<1; j++) {
      for (k=N/2+1; k<N; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*(j+N*(N-i))][0];
        uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*(N-i))][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=N/2+1; k<N; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*((N-j)+N*i)][0];
        uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*((N-j)+N*i)][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=0; j<1; j++) {
      for (k=N/2+1; k<N; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*(j+N*i)][0];
        uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*i)][1];
      }
    }
  }
  for (i=0; i<1; i++) {
    for (j=N/2+1; j<N; j++) {
      for (k=0; k<1; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[k+N*((N-j)+N*i)][0];
        uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*i)][1];
      }
    }
  }
  for (i=N/2+1; i<N; i++) {
    for (j=0; j<1; j++) {
      for (k=0; k<1; k++) {
        uhat[k+N*(j+N*i)][0] =  uhat[k+N*(j+N*(N-i))][0];
        uhat[k+N*(j+N*i)][1] = -uhat[k+N*(j+N*(N-i))][1];
      }
    }
  }

  /* Inverse Fourier transform */
  fftw_execute(inv_trans_u);
  fftw_free(uhat);
  fftw_destroy_plan(inv_trans_u);

  double imag_velocity = 0;
  for (i = 0; i < N3; i++){
    double uu = u[i][1];
    imag_velocity += (uu*uu);
  }
  imag_velocity = sqrt(imag_velocity / ((double)N3));
  printf("RMS of imaginary component of computed velocity: %1.6e\n",imag_velocity);

  for (i = 0; i < N3; i++){
    uvel[i] = u[i][0];
  }

  fftw_free(u);

  return;
}

void setVelocityField(const int N,
                      double* const uvel,
                      double* const vvel,
                      double* const wvel )
{
  double kp = 4.0;
  double u0 = 0.3;

  long N3 = N*N*N;
  long i,j,k;
  double dx = 2*PI / ((double)N);

  printf("Computing u-velocity field.\n");
  velocityComponent( N, 0, kp, u0, uvel );
  printf("Computing v-velocity field.\n");
  velocityComponent( N, 1, kp, u0, vvel );
  printf("Computing w-velocity field.\n");
  velocityComponent( N, 2, kp, u0, wvel );

  double rms_velocity = 0;
  for (i = 0; i < N3; i++){
    double uu, vv, ww;
    uu = uvel[i];
    vv = vvel[i];
    ww = wvel[i];
    rms_velocity += (uu*uu + vv*vv + ww*ww);
  }
  rms_velocity = sqrt(rms_velocity / (3*((double)N3)));

  /* scale the velocity components so that rms velocity matches u0 */
  double factor = u0 / rms_velocity;
  printf("Scaling factor = %1.16E\n",factor);
  for (i = 0; i < N3; i++){
    uvel[i] *= factor;
    vvel[i] *= factor;
    wvel[i] *= factor;
  }

  rms_velocity = 0;
  for (i = 0; i < N3; i++){
    double uu, vv, ww;
    uu = uvel[i];
    vv = vvel[i];
    ww = wvel[i];
    rms_velocity += (uu*uu + vv*vv + ww*ww);
  }
  rms_velocity = sqrt(rms_velocity / (3*((double)N3)));
  printf("RMS velocity (component-wise): %1.16E\n",rms_velocity);

  /* calculate the divergence of velocity */
  double DivergenceNorm = 0;
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
        double u1, u2, v1, v2, w1, w2;
        u1 = (i==0   ? uvel[k+N*(j+N*(N-1))] : uvel[k+N*(j+N*(i-1))] );
        u2 = (i==N-1 ? uvel[k+N*(j+N*(0  ))] : uvel[k+N*(j+N*(i+1))] );
        v1 = (j==0   ? vvel[k+N*((N-1)+N*i)] : vvel[k+N*((j-1)+N*i)] );
        v2 = (j==N-1 ? vvel[k+N*((0  )+N*i)] : vvel[k+N*((j+1)+N*i)] );
        w1 = (k==0   ? wvel[(N-1)+N*(j+N*i)] : wvel[(k-1)+N*(j+N*i)] );
        w2 = (k==N-1 ? wvel[(0  )+N*(j+N*i)] : wvel[(k+1)+N*(j+N*i)] );
        double Divergence = ( (u2-u1) + (v2-v1) + (w2-w1) ) / (2.0*dx);
        DivergenceNorm += (Divergence*Divergence);
      }
    }
  }
  DivergenceNorm = sqrt(DivergenceNorm / ((double)N3));
  printf("Velocity divergence: %1.16E\n",DivergenceNorm);

  /* calculate the Taylor microscales */
  double TaylorMicroscale[3];
  double Numerator[3] = {0,0,0};
  double Denominator[3] = {0,0,0};
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
        double u1, u2, uu, v1, v2, vv, w1, w2, ww;
        u1 = (i==0   ? uvel[k+N*(j+N*(N-1))] : uvel[k+N*(j+N*(i-1))] );
        u2 = (i==N-1 ? uvel[k+N*(j+N*(0  ))] : uvel[k+N*(j+N*(i+1))] );
        v1 = (j==0   ? vvel[k+N*((N-1)+N*i)] : vvel[k+N*((j-1)+N*i)] );
        v2 = (j==N-1 ? vvel[k+N*((0  )+N*i)] : vvel[k+N*((j+1)+N*i)] );
        w1 = (k==0   ? wvel[(N-1)+N*(j+N*i)] : wvel[(k-1)+N*(j+N*i)] );
        w2 = (k==N-1 ? wvel[(0  )+N*(j+N*i)] : wvel[(k+1)+N*(j+N*i)] );
        uu  = uvel[k+N*(j+N*i)];
        vv  = vvel[k+N*(j+N*i)];
        ww  = wvel[k+N*(j+N*i)];

        double du, dv, dw;
        du = (u2 - u1) / (2.0*dx);
        dv = (v2 - v1) / (2.0*dx);
        dw = (w2 - w1) / (2.0*dx);

        Numerator[0] += (uu*uu);
        Numerator[1] += (vv*vv);
        Numerator[2] += (ww*ww);

        Denominator[0] += (du*du);
        Denominator[1] += (dv*dv);
        Denominator[2] += (dw*dw);
      }
    }
  }
  Numerator[0] /= (N*N*N); Denominator[0] /= (N*N*N);
  Numerator[1] /= (N*N*N); Denominator[1] /= (N*N*N);
  Numerator[2] /= (N*N*N); Denominator[2] /= (N*N*N);

  TaylorMicroscale[0] = sqrt(Numerator[0]/Denominator[0]);
  TaylorMicroscale[1] = sqrt(Numerator[1]/Denominator[1]);
  TaylorMicroscale[2] = sqrt(Numerator[2]/Denominator[2]);

  printf("Taylor microscales: %1.16E, %1.16E, %1.16E\n",
         TaylorMicroscale[0],TaylorMicroscale[1],TaylorMicroscale[2]);

  return;
}

void GetStringFromInteger(int a,char *A,int width)
{
  int i;
  for (i=0; i<width; i++) {
    char digit = (char) (a%10 + '0');
    a /= 10;
    A[width-1-i] = digit;
  }
  return;
}

void MPIGetFilename(char *root,int rank,char *filename)
{
  char  tail[_MAX_STRING_SIZE_] = "";

  GetStringFromInteger(rank,tail,4);
  strcpy(filename,"");
  strcat(filename,root);
  strcat(filename,"." );
  strcat(filename,tail);

  return;
}

int MPIRanknD(int ndims,int rank,int* iproc,int *ip)
{
  int i,term    = 1;
  for (i=0; i<ndims; i++) term *= iproc[i];
  for (i=ndims-1; i>=0; i--) {
    term /= iproc[i];
    ip[i] = rank/term;
    rank -= ip[i]*term;
  }
  return(0);
}

int MPIPartition1D(int nglobal,int nproc,int rank)
{
  int nlocal;
  if (nglobal%nproc == 0) nlocal = nglobal/nproc;
  else {
    if (rank == nproc-1)  nlocal = nglobal/nproc + nglobal%nproc;
    else                  nlocal = nglobal/nproc;
  }
  return(nlocal);
}

int MPILocalDomainLimits(int ndims,int p,int *iproc,int *dim_global,int *is, int *ie)
{
  int i;
  int ip[ndims];
  MPIRanknD(ndims,p,iproc,ip);

  for (i=0; i<ndims; i++) {
    int imax_local, isize, root = 0;
    imax_local = MPIPartition1D(dim_global[i],iproc[i],root );
    isize      = MPIPartition1D(dim_global[i],iproc[i],ip[i]);
    if (is)  is[i] = ip[i]*imax_local;
    if (ie)  ie[i] = ip[i]*imax_local + isize;
  }
  return(0);
}

int main()
{
  FILE *out, *in;
  int NI, NJ, NK, ndims, nvars, size, bytes, N_IORanks, i, j, k;
  int *dim_global,*dim_local,*iproc;
  char ip_file_type[_MAX_STRING_SIZE_], input_mode[_MAX_STRING_SIZE_], fnameout[_MAX_STRING_SIZE_];
  strcpy(ip_file_type,"ascii");
  strcpy(fnameout,"initial_par.inp");

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {

    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return 0;

  } else {

    char word[_MAX_STRING_SIZE_];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) {
          fscanf(in,"%d",&ndims);
          dim_global = (int*) calloc (ndims,sizeof(int));
          dim_local  = (int*) calloc (ndims,sizeof(int));
          iproc      = (int*) calloc (ndims,sizeof(int));
        }  else if (!strcmp(word, "nvars")) {
          fscanf(in,"%d",&nvars);
        } else if (!strcmp(word, "size")) {
          int i;
          if (!dim_global) {
            fprintf(stderr,"Error in ReadInputs(): dim_global not allocated.\n");
            fprintf(stderr,"Please specify ndims before dimensions.\n"         );
            return 0;
          } else {
            for (i=0; i<ndims; i++) fscanf(in,"%d",&dim_global[i]);
          }
        } else if (!strcmp(word, "iproc")) {
          int i;
          if (!iproc) {
            fprintf(stderr,"Error in ReadInputs(): iproc not allocated.\n");
            fprintf(stderr,"Please specify ndims before iproc.\n"         );
            return 0;
          } else {
            for (i=0; i<ndims; i++) fscanf(in,"%d",&iproc[i]);
          }
        } else if (!strcmp(word, "ip_file_type" )) {
          fscanf(in,"%s",ip_file_type);
        } else if (!strcmp(word, "input_mode")) {
          fscanf(in,"%s",input_mode);
          if (strcmp(input_mode,"serial")) fscanf(in,"%d",&N_IORanks);
        }
      }
    } else {
      fprintf(stderr,"Error: Illegal format in file \"solver.inp\".\n");
      return 0;
    }
    fclose(in);

    /* Print to screen the inputs read */
    printf("\tNo. of dimensions                          : %d\n",ndims);
    printf("\tNo. of variables                           : %d\n",nvars);
    printf("\tDomain size                                : ");
    for (i=0; i<ndims; i++) printf ("%d ",dim_global[i]);
    printf("\n");
    printf("\tProcesses along each dimension             : ");
    for (i=0; i<ndims; i++) printf ("%d ",iproc[i]);
    printf("\n");
    printf("\tInitial solution file type                 : %s\n",ip_file_type);
    printf("\tInitial solution read mode                 : %s\n",input_mode  );
    printf("\tNumber of IO ranks                         : %d\n",N_IORanks   );
  }

  if (ndims != 3) {
    printf("ndims is not 3 in solver.inp. this code is to generate 3D exact conditions\n");
    return(0);
  }
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Error: ip_file_type *must* be specified and set to \"binary\" in solver.inp.\n");
    return(0);
  }
  if (strcmp(input_mode,"parallel")) {
    printf("Error: input_mode is not \"parallel\".\n");
    return 0;
  }

  NI = dim_global[0];
  NJ = dim_global[1];
  NK = dim_global[2];
  printf("Grid:  %d x %d x %d\n",NI,NJ,NK);

  if ((NI != NJ) || (NI != NK) || (NJ != NK)) {
    printf("Error: NI,NJ,NK not equal. Bye!\n");
    return(0);
  }
  long N = NI;
  long N3 = (long) N * (long) N * (long) N;
  double dx = 2*PI / ((double)N);

  /* Calculating the velocity components through a Fourier transform */
  double* u = (double*) calloc (N3, sizeof(double));
  double* v = (double*) calloc (N3, sizeof(double));
  double* w = (double*) calloc (N3, sizeof(double));
  setVelocityField(N, u, v, w);

  printf("Generating grid.\n");
  double *Xg = (double*) calloc (NI+NJ+NK, sizeof(double));
  double *x = Xg, *y = Xg+NI, *z = Xg+NI+NJ;
  for (i = 0; i < NI; i++){
    for (j = 0; j < NJ; j++){
      for (k = 0; k < NK; k++){
        x[i] = i*dx;
        y[j] = j*dx;
        z[k] = k*dx;
      }
    }
  }

  int nproc = 1;
  for (i=0; i<ndims; i++) nproc *= iproc[i];
  if (nproc%N_IORanks != 0) N_IORanks = 1;
  printf("Splitting data into %d processes. Will generate %d files (one for each file IO rank).\n",
         nproc,N_IORanks);

  int proc,IORank;
  int GroupSize = nproc / N_IORanks;
  for (IORank = 0; IORank < N_IORanks; IORank++) {
    printf("Generating and writing local solutions for IORank %d.\n",IORank);
    char out_filename[_MAX_STRING_SIZE_];
    MPIGetFilename(fnameout,IORank,out_filename);

    int Start = IORank      * GroupSize;
    int End   = (IORank+1)  * GroupSize;

    out = fopen(out_filename,"wb");
    for (proc=Start; proc < End; proc++) {

      int ip[ndims],is[ndims],ie[ndims];
      double *Xl, *Ul;
      MPIRanknD(ndims,proc,iproc,ip);
      MPILocalDomainLimits(ndims,proc,iproc,dim_global,is,ie);
      for (i=0; i<ndims; i++) dim_local[i] = ie[i]-is[i];

      size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
      Xl = (double*) calloc (size, sizeof(double));
      int offsetl=0, offsetg=0;
      for (i=0; i<ndims; i++) {
        int p; for (p=0; p<dim_local[i]; p++) Xl[p+offsetl] = Xg[p+is[i]+offsetg];
        offsetl += dim_local[i];
        offsetg += dim_global[i];
      }
      x = Xl;
      y = Xl + dim_local[0];
      z = Xl + dim_local[0] + dim_local[1];

      size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
      Ul = (double*) calloc (size, sizeof(double));
      int done = 0; int index[ndims]; for(i=0; i<ndims; i++) index[i]=0;
      while (!done) {
        int p; _ArrayIndex1D_(ndims,dim_local ,index,0,p);
        i = index[0]; j = index[1]; k = index[2];

        int ig = i + is[0];
        int jg = j + is[1];
        int kg = k + is[2];

        double rho, uvel, vvel, wvel, P, E;
        rho = 1.0;
        uvel = u[kg+N*(jg+N*ig)];
        vvel = v[kg+N*(jg+N*ig)];
        wvel = w[kg+N*(jg+N*ig)];
        P = 1.0/1.4;
        E = P/0.4 + 0.5*rho*(uvel*uvel+vvel*vvel+wvel*wvel);

        Ul[5*p+0] = rho;
        Ul[5*p+1] = rho*uvel;
        Ul[5*p+2] = rho*vvel;
        Ul[5*p+3] = rho*wvel;
        Ul[5*p+4] = E;
        _ArrayIncrementIndex_(ndims,dim_local,index,done);
      }

      size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
      bytes = fwrite(Xl,sizeof(double),size,out);
      if (bytes != size) printf("Error: Unable to write grid data to file %s.\n",fnameout);
      size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
      bytes = fwrite(Ul,sizeof(double),size,out);
      if (bytes != size) printf("Error: Unable to write solution data to file %s.\n",fnameout);

      free(Xl);
      free(Ul);
    }
    fclose(out);

  }

  free(dim_global);
  free(dim_local);
  free(iproc);
  free(Xg);
  free(u);
  free(v);
  free(w);

  return 0;
}
