/*
  This code computes the energy spectrum of a 3D compressible flow
  solution.

  + The grid dimensions must be the same in all three dimensions.

  + The solution must be in HyPar's binary solution format written out
    through the parallel mode (op.bin.xxxx).

  For unsteady output, the code will write out text files containing
  energy spectrum at each output time step (one file for each time step).
  spectrum_xxxxx.dat

  For steady output, the code will write out one ASCII text containing
  the energy spectrum, spectrum.dat.

  This code requires the FFTW3 library (https://www.fftw.org/). Make sure
  the fftw3.h and libfftw3 are available in the include and linking paths.
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

#define _ArrayIndex1D_(N,imax,i,index)  \
  { \
    index = i[N-1]; \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter])) + (i[arraycounter])); \
    } \
  }

#define _ArrayIndex1DWO_(N,imax,i,offset,index) \
  { \
    index = 0;\
    int arraycounter; \
    for (arraycounter = 0; arraycounter < (N); arraycounter++) { \
      index = ((index*(imax[arraycounter])) + (i[arraycounter]+offset[arraycounter]));\
    } \
  }

void fourierTransform(  const long N,
                        double* const uu,
                        fftw_complex* const uhat )
{
  long i,j,k;
  long N3 = N*N*N;

  fftw_complex *u = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
  fftw_plan transform_u;
  transform_u = fftw_plan_dft_3d(N, N, N, u, uhat, -1, FFTW_MEASURE);

  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        u[k+N*(j+N*i)][0] = uu[k+N*(j+N*i)];
        u[k+N*(j+N*i)][1] = 0;
      }
    }
  }

  fftw_execute(transform_u);

  fftw_free(u);

  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        uhat[k+N*(j+N*i)][0] /= ((double)N3);
        uhat[k+N*(j+N*i)][1] /= ((double)N3);
      }
    }
  }

  fftw_destroy_plan(transform_u);
  return;
}

void energySpectrum(const long N,
                    const fftw_complex* const uhat,
                    const fftw_complex* const vhat,
                    const fftw_complex* const what,
                    const char* const fname )
{
  long i,j,k;
  long N3 = N*N*N;
  double kk = sqrt(3 * (N/2)*(N/2));
  int kkmax = (int) kk;

  double *freq = (double*) calloc(kkmax+1, sizeof(double));
  double *Eng = (double*) calloc(kkmax+1, sizeof(double));

  double total_energy = 0.0;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        long p = k+N*(j+N*i);
        long isq, jsq, ksq;
        if (i > N/2)  isq = (i-N) * (i-N);
        else    isq = i*i;
        if (j > N/2)  jsq = (j-N) * (j-N);
        else    jsq = j*j;
        if (k > N/2)  ksq = (k-N) * (k-N);
        else    ksq = k*k;
        double kk = sqrt(isq + jsq + ksq);
        freq[(int)kk] = kk;
        Eng[(int)kk] = Eng[(int)kk]
          + 0.5 * ( (uhat[p][0]*uhat[p][0] + uhat[p][1]*uhat[p][1])
                  + (vhat[p][0]*vhat[p][0] + vhat[p][1]*vhat[p][1])
                  + (what[p][0]*what[p][0] + what[p][1]*what[p][1]) );
        total_energy = total_energy
          + 0.5 * ( (uhat[p][0]*uhat[p][0] + uhat[p][1]*uhat[p][1])
                  + (vhat[p][0]*vhat[p][0] + vhat[p][1]*vhat[p][1])
                  + (what[p][0]*what[p][0] + what[p][1]*what[p][1]) );
      }
    }
  }
  printf("  Total Energy: %1.16E\n",total_energy);

  printf("Writing energy spectrum file %s.\n", fname);
  FILE *out;
  out = fopen(fname,"w");
  for (i = 1; i < kkmax; i++) fprintf(out,"%1.16E\t%1.16E\n",freq[i],Eng[i]/total_energy);
  fclose(out);
  free(freq);
  free(Eng);

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

void IncrementFilename(char *f)
{
  if (f[13] == '9') {
    f[13] = '0';
    if (f[12] == '9') {
      f[12] = '0';
      if (f[11] == '9') {
        f[11] = '0';
        if (f[10] == '9') {
          f[10] = '0';
          if (f[9] == '9') {
            f[9] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[9]++;
          }
        } else {
          f[10]++;
        }
      } else {
        f[11]++;
      }
    } else {
      f[12]++;
    }
  } else {
    f[13]++;
  }
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

void computeEnergySpectrum( const int NI,
                            double* const u,
                            double* const v,
                            double* const w,
                            const char* const fname )
{
  long i,j,k;
  long N = NI;
  long N3 = N*N*N;

  double rms_velocity = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++) {
        long p = k+N*(j+N*i);
        rms_velocity +=  (u[p]*u[p] + v[p]*v[p] + w[p]*w[p]);
      }
    }
  }
  rms_velocity = sqrt(rms_velocity / (3*N3));
  printf("  RMS velocity (component-wise): %1.16E\n",rms_velocity);

  /* calculate the divergence of velocity */
  double pi = 4.0*atan(1.0);
  double dx = 2*pi / ((double)N);
  double DivergenceNorm = 0;
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
        double u1, u2, v1, v2, w1, w2;
        u1 = (i==0   ? u[k+N*(j+N*(N-1))] : u[k+N*(j+N*(i-1))] );
        u2 = (i==N-1 ? u[k+N*(j+N*(0  ))] : u[k+N*(j+N*(i+1))] );
        v1 = (j==0   ? v[k+N*((N-1)+N*i)] : v[k+N*((j-1)+N*i)] );
        v2 = (j==N-1 ? v[k+N*((0  )+N*i)] : v[k+N*((j+1)+N*i)] );
        w1 = (k==0   ? w[(N-1)+N*(j+N*i)] : w[(k-1)+N*(j+N*i)] );
        w2 = (k==N-1 ? w[(0  )+N*(j+N*i)] : w[(k+1)+N*(j+N*i)] );
        double Divergence = ( (u2-u1) + (v2-v1) + (w2-w1) ) / (2.0*dx);
        DivergenceNorm += (Divergence*Divergence);
      }
    }
  }
  DivergenceNorm = sqrt(DivergenceNorm / ((double)N3));
  printf("  Velocity divergence: %1.16E\n",DivergenceNorm);

  /* calculate the Taylor microscales */
  double TaylorMicroscale[3];
  double Numerator[3] = {0,0,0};
  double Denominator[3] = {0,0,0};
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
        double u1, u2, uc, v1, v2, vc, w1, w2, wc;
        u1 = (i==0   ? u[k+N*(j+N*(N-1))] : u[k+N*(j+N*(i-1))] );
        u2 = (i==N-1 ? u[k+N*(j+N*(0  ))] : u[k+N*(j+N*(i+1))] );
        v1 = (j==0   ? v[k+N*((N-1)+N*i)] : v[k+N*((j-1)+N*i)] );
        v2 = (j==N-1 ? v[k+N*((0  )+N*i)] : v[k+N*((j+1)+N*i)] );
        w1 = (k==0   ? w[(N-1)+N*(j+N*i)] : w[(k-1)+N*(j+N*i)] );
        w2 = (k==N-1 ? w[(0  )+N*(j+N*i)] : w[(k+1)+N*(j+N*i)] );
        uc  = u[k+N*(j+N*i)];
        vc  = v[k+N*(j+N*i)];
        wc  = w[k+N*(j+N*i)];

        double du, dv, dw;
        du = (u2 - u1) / (2.0*dx);
        dv = (v2 - v1) / (2.0*dx);
        dw = (w2 - w1) / (2.0*dx);

        Numerator[0] += (uc*uc);
        Numerator[1] += (vc*vc);
        Numerator[2] += (wc*wc);

        Denominator[0] += (du*du);
        Denominator[1] += (dv*dv);
        Denominator[2] += (dw*dw);
      }
    }
  }
  Numerator[0] /= ((double)N3); Denominator[0] /= ((double)N3);
  Numerator[1] /= ((double)N3); Denominator[1] /= ((double)N3);
  Numerator[2] /= ((double)N3); Denominator[2] /= ((double)N3);

  TaylorMicroscale[0] = sqrt(Numerator[0]/Denominator[0]);
  TaylorMicroscale[1] = sqrt(Numerator[1]/Denominator[1]);
  TaylorMicroscale[2] = sqrt(Numerator[2]/Denominator[2]);

  printf("  Taylor microscales: %1.16E, %1.16E, %1.16E\n",
         TaylorMicroscale[0],TaylorMicroscale[1],TaylorMicroscale[2]);

  fftw_complex *uhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
  printf("  Computing Fourier transform (u).\n");
  fourierTransform( N, u, uhat );
  free(u);

  fftw_complex *vhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
  printf("  Computing Fourier transform (v).\n");
  fourierTransform( N, v, vhat );
  free(v);

  fftw_complex *what = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
  printf("  Computing Fourier transform (w).\n");
  fourierTransform( N, w, what );
  free(w);

  energySpectrum(N, uhat, vhat, what, fname );

  fftw_free(uhat);
  fftw_free(vhat);
  fftw_free(what);

  return;
}

int main()
{
  FILE    *in;
  int     ndims, nvars, i, N_IORanks, proc, GroupSize;
  long    size, bytes;
  int     *dim_global,*dim_local,*iproc,IORank; proc;
  char    output_mode[_MAX_STRING_SIZE_], op_overwrite[_MAX_STRING_SIZE_];

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
            return(0);
          } else for (i=0; i<ndims; i++) fscanf(in,"%d",&dim_global[i]);
        } else if (!strcmp(word, "iproc")) {
          int i;
          if (!iproc) {
            fprintf(stderr,"Error in ReadInputs(): iproc not allocated.\n");
            fprintf(stderr,"Please specify ndims before iproc.\n"         );
            return(0);
          } else for (i=0; i<ndims; i++) fscanf(in,"%d",&iproc[i]);
        } else if (!strcmp(word, "output_mode")) {
          fscanf(in,"%s",output_mode);
          if (strcmp(output_mode,"serial")) fscanf(in,"%d",&N_IORanks);
        } else if (!strcmp(word, "op_overwrite")) {
          fscanf(in,"%s",op_overwrite);
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
    printf("\tSolution output  mode                      : %s\n",output_mode  );
    printf("\tNumber of IO ranks                         : %d\n",N_IORanks   );
  }

  if (strcmp(output_mode,"parallel")) {
    printf("Error: output_mode is not \"parallel\". Why are you using this code?\n");
    return 0;
  }

  FILE *inps[N_IORanks];

  /* open the files for each IO group */
  printf("Opening op.bin.xxxx files.\n");
  for (IORank = 0; IORank < N_IORanks; IORank++) {
    char filename[_MAX_STRING_SIZE_];
    MPIGetFilename("op.bin",IORank,filename);
    inps[IORank] = fopen(filename,"rb");
    if (!inps[IORank]) {
      printf("Error: Could not open %s for reading.\n",filename);
      return(0);
    }
  }

  int nproc = 1;
  for (i=0; i<ndims; i++) nproc *= iproc[i];
  if (nproc%N_IORanks != 0) {
    printf("Error: nproc is not a multiple of N_IORanks. HyPar could not have run in this mode. Something wrong/fishy!\n");
    return 0;
  }

  if (!strcmp(op_overwrite,"no")) {
    /* for unsteady solution output, the output file for each IO group
       contains the solutions for all output time steps */
    char out_filename[_MAX_STRING_SIZE_];
    strcpy(out_filename,"spectrum_00000.dat");
    int flag = 1;
    int count = 1;
    while (flag) {
      printf("Count = %d.\n", count);
      count++;
      /* allocate global velocity array */
      size = 1;
      for (i=0; i<ndims; i++) size *= (long) dim_global[i];
      double* Ug = (double*) calloc (size, sizeof(double));
      double* Vg = (double*) calloc (size, sizeof(double));
      double* Wg = (double*) calloc (size, sizeof(double));

      for (IORank=0; IORank < N_IORanks; IORank++) {
        /* for each IO group, calculate its range of processes */
        GroupSize = nproc / N_IORanks;
        int Start = IORank      * GroupSize;
        int End   = (IORank+1)  * GroupSize;
        /* for each process in this IO group, read the solution and
           put it in the global arrays */
        for (proc=Start; proc < End; proc++) {
          /* calculate local dimensions */
          int ip[ndims],is[ndims],ie[ndims];
          MPIRanknD(ndims,proc,iproc,ip);
          MPILocalDomainLimits(ndims,proc,iproc,dim_global,is,ie);
          for (i=0; i<ndims; i++) dim_local[i] = ie[i]-is[i];

          /* allocate local arrays for grid and solution */
          double *Xl, *Ul;
          size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
          Xl = (double*) calloc (size, sizeof(double));
          size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
          Ul = (double*) calloc (size, sizeof(double));

          /* read the local solution */
          size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
          bytes = fread(Xl,sizeof(double),size,inps[IORank]);
          if (bytes != size) printf("Error: Unable to read grid.\n");
          size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
          bytes = fread(Ul,sizeof(double),size,inps[IORank]);
          if (bytes != size) printf("Error: Unable to read solution.\n");

          /* copy to global arrays */
          int done = 0; int index[ndims]; for(i=0; i<ndims; i++) index[i]=0;
          while (!done) {
            long p1; _ArrayIndex1DWO_(ndims,dim_global,index,is,p1);
            int p2; _ArrayIndex1D_(ndims,dim_local,index,p2);
            double rho = Ul[nvars*p2+0];
            Ug[p1] = Ul[nvars*p2+1] / rho;
            Vg[p1] = Ul[nvars*p2+2] / rho;
            Wg[p1] = Ul[nvars*p2+3] / rho;
            _ArrayIncrementIndex_(ndims,dim_local,index,done);
          }

          /* free local grid and solution arrays */
          free(Xl);
          free(Ul);
        }
      }

      /* compute kinetic energy spectrum */
      computeEnergySpectrum( dim_global[0], Ug, Vg, Wg, out_filename);

      IncrementFilename(out_filename);

      /* check if the local solution files have reached their end */
      for (IORank=0; IORank<N_IORanks; IORank++) {
        int check = fgetc(inps[IORank]);
        if (check == EOF) flag = 0;
        else fseek(inps[IORank],-sizeof(char),SEEK_CUR);
      }
    }

  } else {

    /* for steady solution */
    char out_filename[_MAX_STRING_SIZE_];
    strcpy(out_filename,"spectrum.dat");
    /* allocate global grid and solution arrays */
    size = 1;
    for (i=0; i<ndims; i++) size *= (long) dim_global[i];
    double* Ug = (double*) calloc (size, sizeof(double));
    double* Vg = (double*) calloc (size, sizeof(double));
    double* Wg = (double*) calloc (size, sizeof(double));

    for (IORank=0; IORank < N_IORanks; IORank++) {
      /* for each IO group, calculate its range of processes */
      GroupSize = nproc / N_IORanks;
      int Start = IORank      * GroupSize;
      int End   = (IORank+1)  * GroupSize;
      /* for each process in this IO group, read the solution and
         put it in the global arrays */
      for (proc=Start; proc < End; proc++) {
        /* calculate local dimensions */
        int ip[ndims],is[ndims],ie[ndims];
        MPIRanknD(ndims,proc,iproc,ip);
        MPILocalDomainLimits(ndims,proc,iproc,dim_global,is,ie);
        for (i=0; i<ndims; i++) dim_local[i] = ie[i]-is[i];

        /* allocate local arrays for grid and solution */
        double *Xl, *Ul;
        size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
        Xl = (double*) calloc (size, sizeof(double));
        size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
        Ul = (double*) calloc (size, sizeof(double));

        /* read the local solution */
        size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
        bytes = fread(Xl,sizeof(double),size,inps[IORank]);
        if (bytes != size) printf("Error: Unable to read grid.\n");
        size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
        bytes = fread(Ul,sizeof(double),size,inps[IORank]);
        if (bytes != size) printf("Error: Unable to read solution.\n");

        /* copy to global arrays */
        int done = 0; int index[ndims]; for(i=0; i<ndims; i++) index[i]=0;
        while (!done) {
          long p1; _ArrayIndex1DWO_(ndims,dim_global,index,is,p1);
          int p2; _ArrayIndex1D_(ndims,dim_local,index,p2);
          double rho = Ul[nvars*p2+0];
          Ug[p1] = Ul[nvars*p2+1] / rho;
          Vg[p1] = Ul[nvars*p2+2] / rho;
          Wg[p1] = Ul[nvars*p2+3] / rho;
          _ArrayIncrementIndex_(ndims,dim_local,index,done);
        }

        /* free local grid and solution arrays */
        free(Xl);
        free(Ul);
      }
    }

    /* compute kinetic energy spectrum */
    computeEnergySpectrum( dim_global[0], Ug, Vg, Wg, out_filename);

  }

  /* close the files for each IO group */
  for (IORank = 0; IORank < N_IORanks; IORank++) fclose(inps[IORank]);

  free(dim_local);
  free(dim_global);
  free(iproc);
  return(0);
}
