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


int main()
{
  /* shock formation time */
  double ts = 2.0;

  /* maximum number of iterations for computing exact solution */
  int MAX_ITER = 10000;

  /* Convergence tolerance for computing exact solution */
  double tolerance = 1e-15;

  /* value of pi */
  double pi = 4.0*atan(1.0);

  int NI, ndims, niter;
  double dt;

  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims"))             fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size"))         fscanf(in,"%d",&NI);
        else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
        else if (!strcmp(word, "n_iter"))   fscanf(in,"%d" , &niter);
        else if (!strcmp(word, "dt"))       fscanf(in,"%lf",    &dt);
      }
    } else {
      printf("Error: Illegal format in solver.inp. Crash and burn!\n");
    }
  }
  fclose(in);

  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n");
    return(0);
  }
  printf("Grid:\t\t\t%d\n",NI);

  int i;
  double dx = 1.0 / ((double)NI);
  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  /* initial solution */
  {
    double *x, *u;
    x = (double*) calloc (NI, sizeof(double));
    u = (double*) calloc (NI, sizeof(double));
    for (i = 0; i < NI; i++){
      x[i] = i*dx;
      u[i] = sin(2.0*pi*x[i]) / (ts*2.0*pi);
    }
    FILE *out;
    if (!strcmp(ip_file_type,"ascii")) {
      printf("Writing ASCII initial solution file initial.inp\n");
      out = fopen("initial.inp","w");
      for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
      fprintf(out,"\n");
      for (i = 0; i < NI; i++)  fprintf(out,"%lf ",u[i]);
      fprintf(out,"\n");
      fclose(out);
    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
      printf("Writing binary initial solution file initial.inp\n");
      out = fopen("initial.inp","wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(u,sizeof(double),NI,out);
      fclose(out);
    }
    free(x);
    free(u);
  }

  /* exact solution */
  if (tf < ts) {

    double *x, *u;
    x = (double*) calloc (NI, sizeof(double));
    u = (double*) calloc (NI, sizeof(double));
    for (i = 0; i < NI; i++){
      x[i] = i*dx;
      u[i] = sin(2.0*pi*x[i]) / (ts*2.0*pi);
    }

    int k;
    printf("Computing exact solution iteratively...\n");
    for (k = 0; k < MAX_ITER; k++) {
      double maxres = 0;
      for (i = 0; i < NI; i++){
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
      printf("Writing ASCII exact solution file exact.inp\n");
      out = fopen("exact.inp","w");
      for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
      fprintf(out,"\n");
      for (i = 0; i < NI; i++)  fprintf(out,"%lf ",u[i]);
      fprintf(out,"\n");
      fclose(out);
    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
      printf("Writing binary exact solution file exact.inp\n");
      out = fopen("exact.inp","wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(u,sizeof(double),NI,out);
      fclose(out);
    }
    free(x);
    free(u);

  } else {

    printf("Final time (%f) greater than shock formation time (%f).\n", tf, ts);
    printf("Exact solution not available.\n");

  }

  return(0);
}
