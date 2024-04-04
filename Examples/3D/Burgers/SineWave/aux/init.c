/*
  Code to generate the initial for:
  Case: Sine Wave
  Model: Burger3D

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

  int NI, NJ, NK, ndims, niter;
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
          fscanf(in,"%d",&NK);
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

  if (ndims != 3) {
    fprintf(stderr,"ndims is not 3 in solver.inp. this code is to generate 3D initial conditions\n");
    return(0);
  }
  printf("Grid: %d, %d, %d\n", NI, NJ, NK);

  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  int i, j, k;
  double dx = 1.0 / ((double)NI);
  double dy = 1.0 / ((double)NJ);
  double dz = 1.0 / ((double)NK);

  /* initial solution */
  {

    double *x, *y, *z, *u;
    x = (double*) calloc (NI      , sizeof(double));
    y = (double*) calloc (NJ      , sizeof(double));
    z = (double*) calloc (NK      , sizeof(double));
    u = (double*) calloc (NI*NJ*NK, sizeof(double));

    for (i = 0; i < NI; i++) {
      for (j = 0; j < NJ; j++) {
        for (k = 0; k < NK; k++) {
          x[i] = (i+0.5)*dx;
          y[j] = (j+0.5)*dy;
          z[k] = (k+0.5)*dz;
          int p = i + NI*j + NI*NJ*k;
          u[p] = sin(2*pi*x[i]) * sin(2*pi*y[j]) * sin(2*pi*z[k]) / (2*ts*pi);
        }
      }
    }

    FILE *out;
    if (!strcmp(ip_file_type,"ascii")) {
      printf("Writing ASCII initial solution file initial.inp\n");
      out = fopen("initial.inp","w");
      for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
      fprintf(out,"\n");
      for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
      fprintf(out,"\n");
      for (k = 0; k < NK; k++)  fprintf(out,"%1.16e ",z[k]);
      fprintf(out,"\n");
      for (k = 0; k < NK; k++) {
        for (j = 0; j < NJ; j++)  {
          for (i = 0; i < NI; i++)  {
            int p = i + NK*j + NI*NJ*k;
            fprintf(out,"%1.16e ",u[p]);
          }
        }
      }
      fprintf(out,"\n");
      fclose(out);
    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
      printf("Writing binary initial solution file initial.inp\n");
      out = fopen("initial.inp","wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(y,sizeof(double),NJ,out);
      fwrite(z,sizeof(double),NK,out);
      fwrite(u,sizeof(double),NI*NJ*NK,out);
      fclose(out);
    }

    free(x);
    free(y);
    free(z);
    free(u);

  }

  /* exact solution */
  if (tf < ts) {

    double *x, *y, *z, *u;
    x = (double*) calloc (NI      , sizeof(double));
    y = (double*) calloc (NJ      , sizeof(double));
    z = (double*) calloc (NK      , sizeof(double));
    u = (double*) calloc (NI*NJ*NK, sizeof(double));

    for (i = 0; i < NI; i++) {
      for (j = 0; j < NJ; j++) {
        for (k = 0; k < NK; k++) {
          x[i] = (i+0.5)*dx;
          y[j] = (j+0.5)*dy;
          z[k] = (k+0.5)*dz;
          int p = i + NI*j + NI*NJ*k;
          u[p] = sin(2*pi*x[i]) * sin(2*pi*y[j]) * sin(2*pi*z[k]) / (2*ts*pi);
        }
      }
    }

    int iter;
    printf("Computing exact solution iteratively...\n");
    for (iter = 0; iter < MAX_ITER; iter++) {
      double maxres = 0;
      for (i = 0; i < NI; i++) {
        for (j = 0; j < NJ; j++) {
          for (k = 0; k < NK; k++) {
            int p = i + NI*j + NI*NJ*k;
            double new_u =  sin(2.0*pi*(z[k]-u[p]*tf)) * sin(2.0*pi*(y[j]-u[p]*tf)) * sin(2.0*pi*(x[i]-u[p]*tf)) / (ts*2.0*pi);
            double res = sqrt((new_u-u[p])*(new_u-u[p]));
            if (res > maxres) maxres = res;
            u[p] = new_u;
          }
        }
      }
      printf("  iter=%6d, max res=%1.6e\n", iter, maxres);
      if (maxres < tolerance) break;
    }

    FILE *out;
    if (!strcmp(ip_file_type,"ascii")) {
      printf("Writing ASCII exact solution file exact.inp\n");
      out = fopen("exact.inp","w");
      for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
      fprintf(out,"\n");
      for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
      fprintf(out,"\n");
      for (k = 0; k < NK; k++)  fprintf(out,"%1.16e ",z[k]);
      fprintf(out,"\n");
      for (k = 0; k < NK; k++) {
        for (j = 0; j < NJ; j++)  {
          for (i = 0; i < NI; i++)  {
            int p = i + NK*j + NI*NJ*k;
            fprintf(out,"%1.16e ",u[p]);
          }
        }
      }
      fprintf(out,"\n");
      fclose(out);
    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
      printf("Writing binary exact solution file exact.inp\n");
      out = fopen("exact.inp","wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(y,sizeof(double),NJ,out);
      fwrite(z,sizeof(double),NK,out);
      fwrite(u,sizeof(double),NI*NJ*NK,out);
      fclose(out);
    }

    free(x);
    free(y);
    free(z);
    free(u);

  } else {

    printf("Final time (%f) greater than shock formation time (%f).\n", tf, ts);
    printf("Exact solution not available.\n");

  }

  return(0);
}
