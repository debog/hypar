/*
  This code computes the energy spectrum of a 3D compressible flow
  solution.

  + The grid dimensions must be the same in all three dimensions.
  + The solution must be in HyPar's binary solution format.
  + This code looks for "op.bin"; if the simulation was unsteady with
    multiple "op_<idx>.bin" files, create a symbolic link named "op.bin"
    pointing to the specific "op_<idx>.bin" solution file for which the
    spectrum should be computed.

  This code requires the FFTW3 library (https://www.fftw.org/). Make sure
  the fftw3.h and libfftw3 are available in the include and linking paths.
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

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
                    const fftw_complex* const what )
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
  printf("Total Energy: %1.16E\n",total_energy);

  FILE *out;
  out = fopen("spectrum.dat","w");
  for (i = 1; i < kkmax; i++) fprintf(out,"%1.16E\t%1.16E\n",freq[i],Eng[i]/total_energy);
  fclose(out);
  free(freq);
  free(Eng);

  return;
}

int main()
{
  long i,j,k;
  int NI,NJ,NK;
  char op_file_format[50];
  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
          fscanf(in,"%d",&NK);
        } else if (!strcmp(word, "op_file_format")) fscanf(in,"%s",op_file_format);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  printf("Grid size: %d x %d x %d.\n",NI,NJ,NK);
  long N = NI;
  long N3 = N*N*N;

  if ((!strcmp(op_file_format,"binary")) || (!strcmp(op_file_format,"bin"))) {

    FILE *inp;
    double *x = (double*) calloc (N,sizeof(double));
    double *y = (double*) calloc (N,sizeof(double));
    double *z = (double*) calloc (N,sizeof(double));
    double *U = (double*) calloc (5*N3,sizeof(double));
    int ndims, nvars, dims[3];

    double *u,*v,*w;
    u = (double*) calloc (N3, sizeof(double));
    v = (double*) calloc (N3, sizeof(double));
    w = (double*) calloc (N3, sizeof(double));

    printf("Reading binary solution file op.bin\n");
    inp = fopen("op.bin","rb");

    if (!inp) {
      printf("Output file op.bin not found!\n");
      return(0);
    }

    fread(&ndims,sizeof(int),1,inp);
    fread(&nvars,sizeof(int),1,inp);
    fread(dims,sizeof(int),3,inp);
    fread(x,sizeof(double),N,inp);
    fread(y,sizeof(double),N,inp);
    fread(z,sizeof(double),N,inp);
    fread(U,sizeof(double),5*N3,inp);
    fclose(inp);

    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
          long p1 = 5*(i+N*(j+N*k));
          long p2 = k+N*(j+N*i);
          double rho = U[p1];
          u[p2] = U[p1+1]/rho;
          v[p2] = U[p1+2]/rho;
          w[p2] = U[p1+3]/rho;
        }
      }
    }

    free(x);
    free(y);
    free(z);
    free(U);

    double rms_velocity = 0;
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
          long p = k+N*(j+N*i);
          rms_velocity +=  (u[p]*u[p] + v[p]*v[p] + w[p]*w[p]);
        }
      }
    }
    rms_velocity = sqrt(rms_velocity / (3*((double)N3)));
    printf("RMS velocity (component-wise): %1.16E\n",rms_velocity);

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
    printf("Velocity divergence: %1.16E\n",DivergenceNorm);

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

    printf("Taylor microscales: %1.16E, %1.16E, %1.16E\n",
           TaylorMicroscale[0],TaylorMicroscale[1],TaylorMicroscale[2]);

    fftw_complex *uhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
    printf("Computing Fourier transform (u).\n");
    fourierTransform( N, u, uhat );
    free(u);

    fftw_complex *vhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
    printf("Computing Fourier transform (v).\n");
    fourierTransform( N, v, vhat );
    free(v);

    fftw_complex *what = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
    printf("Computing Fourier transform (w).\n");
    fourierTransform( N, w, what );
    free(w);

    energySpectrum(N, uhat, vhat, what);

    fftw_free(uhat);
    fftw_free(vhat);
    fftw_free(what);

  } else {

    printf("Error: Unsupported output file type. Use binary output only!\n");

  }

  return(0);
}
