/*
  Code to generate the initial for:
  Case: 1D Sine Wave
  Model: Burger2D

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

  int NI, NJ, ndims, niter;
  double dt;

  FILE *in;

  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {

    fprintf(stderr,"Error: Input file \"solver.inp\" not found.\n");
    return(0);

  } else {

    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in,"%s",ip_file_type);
        } else if (!strcmp(word, "n_iter")) {
          fscanf(in,"%d" , &niter);
        } else if (!strcmp(word, "dt")) {
          fscanf(in,"%lf", &dt);
        }
      }

    } else {

      fprintf(stderr,"Error: Illegal format in solver.inp. Crash and burn!\n");
      return(0);

    }

  }

  fclose(in);

  if (ndims != 2) {
    fprintf(stderr,"ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n");
    return(0);
  }
  printf("Grid: %d, %d\n",NI,NJ);

  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  int i,j;
  double dx = 1.0 / ((double)NI);
  double dy = 1.0 / ((double)NJ);

  /* initial solution */
  {

    double *x, *y, *u;
    x = (double*) calloc (NI   , sizeof(double));
    y = (double*) calloc (NJ   , sizeof(double));
    u = (double*) calloc (NI*NJ, sizeof(double));

    for (i = 0; i < NI; i++){
      for (j = 0; j < NJ; j++){
        x[i] = i*dx;
        y[j] = j*dy;
        int p = NJ*i + j;
        u[p] = sin(2*pi*x[i]) / (2*ts*pi);
      }
    }

    FILE *out;
    if (!strcmp(ip_file_type,"ascii")) {

      out = fopen("initial.inp","w");
      for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
      fprintf(out,"\n");
      for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
      fprintf(out,"\n");
      for (j = 0; j < NJ; j++)  {
        for (i = 0; i < NI; i++)  {
          int p = NJ*i + j;
          fprintf(out,"%1.16e ",u[p]);
        }
      }
      fprintf(out,"\n");
      fclose(out);

    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

      printf("Writing binary exact solution file initial.inp\n");
      out = fopen("initial.inp","wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(y,sizeof(double),NJ,out);
      double *U = (double*) calloc (NI*NJ,sizeof(double));
      for (i=0; i < NI; i++) {
        for (j = 0; j < NJ; j++) {
          int p = NJ*i + j;
          int q = NI*j + i;
          U[q+0] = u[p];
        }
      }
      fwrite(U,sizeof(double),NI*NJ,out);
      free(U);
      fclose(out);

    }

    free(x);
    free(y);
    free(u);

  }

  /* exact solution */
  if (tf < ts) {

    double *x, *y, *u;
    x = (double*) calloc (NI   , sizeof(double));
    y = (double*) calloc (NJ   , sizeof(double));
    u = (double*) calloc (NI*NJ, sizeof(double));

    for (i = 0; i < NI; i++){
      for (j = 0; j < NJ; j++){
        x[i] = i*dx;
        y[j] = j*dy;
        int p = NJ*i + j;
        u[p] = sin(2*pi*x[i]) / (2*ts*pi);
      }
    }

    int k;
    printf("Computing exact solution iteratively...\n");
    for (k = 0; k < MAX_ITER; k++) {
      double maxres = 0;
      for (i = 0; i < NI; i++){
        for (j = 0; j < NJ; j++){
          int p = NJ*i + j;
          double new_u =  sin(2.0*pi*(x[i]-u[p]*tf)) / (ts*2.0*pi);
          double res = sqrt((new_u-u[p])*(new_u-u[p]));
          if (res > maxres) maxres = res;
          u[p] = new_u;
        }
      }
      printf("  iter=%6d, max res=%1.6e\n", k, maxres);
      if (maxres < tolerance) break;
    }

    FILE *out;
    if (!strcmp(ip_file_type,"ascii")) {

      out = fopen("exact.inp","w");
      for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
      fprintf(out,"\n");
      for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
      fprintf(out,"\n");
      for (j = 0; j < NJ; j++)  {
        for (i = 0; i < NI; i++)  {
          int p = NJ*i + j;
          fprintf(out,"%1.16e ",u[p]);
        }
      }
      fprintf(out,"\n");
      fclose(out);

    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

      printf("Writing binary exact solution file exact.inp\n");
      out = fopen("exact.inp","wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(y,sizeof(double),NJ,out);
      double *U = (double*) calloc (NI*NJ,sizeof(double));
      for (i=0; i < NI; i++) {
        for (j = 0; j < NJ; j++) {
          int p = NJ*i + j;
          int q = NI*j + i;
          U[q+0] = u[p];
        }
      }
      fwrite(U,sizeof(double),NI*NJ,out);
      free(U);
      fclose(out);

    }

    free(x);
    free(y);
    free(u);

  } else {

    printf("Final time (%f) greater than shock formation time (%f).\n", tf, ts);
    printf("Exact solution not available.\n");

  }

  return(0);
}
