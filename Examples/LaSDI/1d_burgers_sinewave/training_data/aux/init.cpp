/*
  Code to generate the initial and exact
  solutions for:
  Case: Sine Wave
  Model: Burger1D

  Needs: solver.inp
  Writes out: initial.inp

  If the final time is less than shock
  formation time, then the exact solution
  is also written out.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#define _MAX_STRING_SIZE_ 200

void GetStringFromInteger(int a, char *A, int width)
{
  int i;
  for (i=0; i<width; i++) {
    char digit = (char) (a%10 + '0');
    a /= 10;
    A[width-1-i] = digit;
  }
  A[width] = 0;
  return;
}

typedef struct _params_
{
  double ts;
} Params;


int main()
{
  /* shock formation time */
  std::vector<double> ts = {1.8, 2.2};
  int n_ts = ts.size();

  /* write out parameter set */
  {
    FILE* out;
    out = fopen("parameters.dat","w");
    fprintf(out, "begin\n");
    fprintf(out, "  t_shock ");
    for (int i = 0; i < n_ts; i++) {
      fprintf(out, "  %f", ts[i]);
    }
    fprintf(out, "\n");
    fprintf(out, "end");
    fclose(out);
  }

  int n_param_sets = n_ts;
  Params *params = (Params*) calloc( n_param_sets, sizeof(Params));
  {
    for (int i = 0; i < n_ts; i++) {
      int s = i;
      params[s].ts = ts[i];
    }
  }

  /* maximum number of iterations for computing exact solution */
  int MAX_ITER = 10000;

  /* Convergence tolerance for computing exact solution */
  double tolerance = 1e-15;

  /* value of pi */
  double pi = 4.0*atan(1.0);

  int nsims = 1;

  FILE* in;
  int ndims, niter;
  double dt;

  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  /* find out the number of simulation domains */
  printf("Reading file \"simulation.inp\" if available.\n");
  in = fopen("simulation.inp","r");
  if (in) {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "nsims")) fscanf(in,"%d",&nsims);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  if (nsims > n_param_sets) {
    printf("Error: parameter set is too small.");
    return 0;
  }
  printf("Number of simulation domains: %d\n", nsims);
  int* NI = (int*) calloc (nsims, sizeof(int));

  /* read solver.inp */
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) {
          fscanf(in,"%d",&ndims);
        } else if (!strcmp(word, "size")) {
          for (int ns = 0; ns < nsims; ns++) {
            fscanf(in, "%d", (NI+ns));
          }
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in,"%s",ip_file_type);
        } else if (!strcmp(word, "n_iter")) {
          fscanf(in,"%d", &niter);
        } else if (!strcmp(word, "dt")) {
          fscanf(in,"%lf", &dt);
        }
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n");
    return(0);
  }
  if (nsims == 1) {
    printf("Grid:\t\t\t%d\n",NI[0]);
  } else {
    printf("Grid sizes:\n");
    for (int ns = 0; ns < nsims; ns++) printf("\t%d\n", NI[ns]);
  }

  int i;
  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  for (int ns = 0; ns < nsims; ns++) {

    double ts = params[ns].ts;
    printf("(ts=%f)\n", ts);

    int N = NI[ns];
    double dx = 1.0 / ((double)N);

    /* initial solution */
    {
      char fname[_MAX_STRING_SIZE_] = "initial";
      if (nsims > 1) {
        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(ns, index, (int)log10(nsims)+1);
        strcat(fname, "_");
        strcat(fname, index);
      }
      strcat(fname, ".inp");

      double *x, *u;
      x = (double*) calloc (N, sizeof(double));
      u = (double*) calloc (N, sizeof(double));
      for (i = 0; i < N; i++){
        x[i] = i*dx;
        u[i] = sin(2.0*pi*x[i]) / (ts*2.0*pi);
      }
      FILE *out;
      if (!strcmp(ip_file_type,"ascii")) {
        printf("Writing ASCII initial solution file %s\n", fname);
        out = fopen(fname,"w");
        for (i = 0; i < N; i++)  fprintf(out,"%lf ",x[i]);
        fprintf(out,"\n");
        for (i = 0; i < N; i++)  fprintf(out,"%lf ",u[i]);
        fprintf(out,"\n");
        fclose(out);
      } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
        printf("Writing binary initial solution file %s\n", fname);
        out = fopen(fname,"wb");
        fwrite(x,sizeof(double),N,out);
        fwrite(u,sizeof(double),N,out);
        fclose(out);
      }
      free(x);
      free(u);
    }

    /* exact solution */
    if (tf < ts) {

      char fname[_MAX_STRING_SIZE_] = "exact";
      if (nsims > 1) {
        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(ns, index, (int)log10(nsims)+1);
        strcat(fname, "_");
        strcat(fname, index);
      }
      strcat(fname, ".inp");

      double *x, *u;
      x = (double*) calloc (N, sizeof(double));
      u = (double*) calloc (N, sizeof(double));
      for (i = 0; i < N; i++){
        x[i] = i*dx;
        u[i] = sin(2.0*pi*x[i]) / (ts*2.0*pi);
      }

      int k;
      printf("Computing exact solution iteratively...\n");
      for (k = 0; k < MAX_ITER; k++) {
        double maxres = 0;
        for (i = 0; i < N; i++){
          double new_u =  sin(2.0*pi*(x[i]-u[i]*tf)) / (ts*2.0*pi);
          double res = sqrt((new_u-u[i])*(new_u-u[i]));
          if (res > maxres) maxres = res;
          u[i] = new_u;
        }
        printf("  iter=%6d, max res=%1.6e\n", k, maxres);
        if (maxres < tolerance) break;
      }

      FILE *out;
      if (!strcmp(ip_file_type,"ascii")) {
        printf("Writing ASCII exact solution file %s\n", fname);
        out = fopen(fname,"w");
        for (i = 0; i < N; i++)  fprintf(out,"%lf ",x[i]);
        fprintf(out,"\n");
        for (i = 0; i < N; i++)  fprintf(out,"%lf ",u[i]);
        fprintf(out,"\n");
        fclose(out);
      } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
        printf("Writing binary exact solution file %s\n", fname);
        out = fopen(fname,"wb");
        fwrite(x,sizeof(double),N,out);
        fwrite(u,sizeof(double),N,out);
        fclose(out);
      }
      free(x);
      free(u);

    } else {

      printf("Final time (%f) greater than shock formation time (%f).\n", tf, ts);
      printf("Exact solution not available.\n");

    }
  }

  return(0);
}
