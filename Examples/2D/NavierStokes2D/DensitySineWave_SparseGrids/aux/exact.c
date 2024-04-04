#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double power(double x,double a)
{
  return(exp(a*log(x)));
}

int main(){
  const double pi     = 4.0*atan(1.0);
  const double GAMMA  = 1.4;

  int     NI,NJ,ndims,n_iter;
  double  tf, dt;
  char    ip_file_type[50]; strcpy(ip_file_type,"ascii");

  FILE *in, *out;

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
        } else if (!strcmp(word, "n_iter"))     fscanf(in,"%d",&n_iter);
        else if (!strcmp(word, "dt"))           fscanf(in,"%lf",&dt);
        else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D exact solution\n");
    return(0);
  }
  printf("Grid:\t\t\t%d X %d\n",NI,NJ);

  int i,j;
  double dx = 1.0 / ((double)NI);
  double dy = 1.0 / ((double)NJ);

  tf = (double)n_iter * dt;
  printf("Final time: %lf\n", tf);

  double *x, *y, *u0, *u1, *u2, *u3;
  x   = (double*) calloc (NI   , sizeof(double));
  y   = (double*) calloc (NJ   , sizeof(double));
  u0  = (double*) calloc (NI*NJ, sizeof(double));
  u1  = (double*) calloc (NI*NJ, sizeof(double));
  u2  = (double*) calloc (NI*NJ, sizeof(double));
  u3  = (double*) calloc (NI*NJ, sizeof(double));

  double rho_inf = 1.0;
  double drho = 0.1;
  double u_inf = 1.0;
  double v_inf = 1.0;
  double P_inf = 1.0/GAMMA;

  /* Initial solution */
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NJ; j++) {
      x[i] = i*dx;
      y[j] = j*dy;
      int p = NJ*i + j;
      double rho, u, v, P;
      rho = rho_inf + drho * sin(2*pi*x[i]) * cos(2*pi*y[j]);
      P   = P_inf;
      u   = u_inf;
      v   = v_inf;
      u0[p] = rho;
      u1[p] = rho*u;
      u2[p] = rho*v;
      u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
    }
  }
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
    out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%lf ",y[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u0[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u1[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u2[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u3[p]);
      }
    }
    fprintf(out,"\n");
    fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary initial solution file initial.inp\n");
    out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    double *U = (double*) calloc (4*NI*NJ,sizeof(double));
    for (i=0; i < NI; i++) {
      for (j = 0; j < NJ; j++) {
        int p = NJ*i + j;
        int q = NI*j + i;
        U[4*q+0] = u0[p];
        U[4*q+1] = u1[p];
        U[4*q+2] = u2[p];
        U[4*q+3] = u3[p];
      }
    }
    fwrite(U,sizeof(double),4*NI*NJ,out);
    free(U);
    fclose(out);
  }

  /* Exact solution */
  for (i = 0; i < NI; i++) {
    for (j = 0; j < NJ; j++) {
      x[i] = i*dx;
      y[j] = j*dy;
      int p = NJ*i + j;
      double rho, u, v, P;
      rho = rho_inf + drho * sin(2*pi*(x[i]-u_inf*tf)) * cos(2*pi*(y[j]-v_inf*tf));
      P   = P_inf;
      u   = u_inf;
      v   = v_inf;
      u0[p] = rho;
      u1[p] = rho*u;
      u2[p] = rho*v;
      u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
    }
  }
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII exact solution file exact.inp\n");
    out = fopen("exact.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%lf ",y[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u0[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u1[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u2[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
      for (i = 0; i < NI; i++)  {
        int p = NJ*i + j;
        fprintf(out,"%lf ",u3[p]);
      }
    }
    fprintf(out,"\n");
    fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary exact solution file exact.inp\n");
    out = fopen("exact.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    double *U = (double*) calloc (4*NI*NJ,sizeof(double));
    for (i=0; i < NI; i++) {
      for (j = 0; j < NJ; j++) {
        int p = NJ*i + j;
        int q = NI*j + i;
        U[4*q+0] = u0[p];
        U[4*q+1] = u1[p];
        U[4*q+2] = u2[p];
        U[4*q+3] = u3[p];
      }
    }
    fwrite(U,sizeof(double),4*NI*NJ,out);
    free(U);
    fclose(out);
  }

  free(x);
  free(y);
  free(u0);
  free(u1);
  free(u2);
  free(u3);

  return(0);
}
