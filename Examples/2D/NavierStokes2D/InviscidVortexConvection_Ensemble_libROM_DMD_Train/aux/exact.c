#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

double power(double x,double a)
{
  return(exp(a*log(x)));
}

int main(){
  const double pi     = 4.0*atan(1.0);
  const double GAMMA  = 1.4;

  int     *NI, *NJ, ndims, ns, nsims, n_iter;
  double  tf, dt;
  char    ip_file_type[50];
  FILE    *in, *out;

  /* default values */
  nsims = 1;
  strcpy(ip_file_type,"binary");

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
  printf("Number of simulation domains: %d\n", nsims);
  NI = (int*) calloc (nsims, sizeof(int));
  NJ = (int*) calloc (nsims, sizeof(int));

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) {
          fscanf(in,"%d",&ndims);
        } else if (!strcmp(word, "size")) {
          for (ns = 0; ns < nsims; ns++) {
            fscanf(in, "%d", (NI+ns));
            fscanf(in, "%d", (NJ+ns));
          }
        } else if (!strcmp(word, "n_iter")) {
          fscanf(in,"%d",&n_iter);
        } else if (!strcmp(word, "dt")) {
          fscanf(in,"%lf",&dt);
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in,"%s",ip_file_type);
        }
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D exact solution\n");
    return(0);
  }

  if (nsims == 1) {
    printf("Grid:\t\t\t%d\n",NI[0], NJ[0]);
  } else {
    printf("Grid sizes:\n");
    for (ns = 0; ns < nsims; ns++) printf("\t%d X %d\n", NI[ns], NJ[ns]);
  }

  tf = (double)n_iter * dt; double tff = tf;
  while (tf > 20) tf -= 20; // Time period

  for (ns = 0; ns < nsims; ns++) {

    int Ni = NI[ns];
    int Nj = NJ[ns];

    int i,j;
    double dx = 10.0 / ((double)Ni);
    double dy = 10.0 / ((double)Nj);

    double *x, *y, *u0, *u1, *u2, *u3;
    x   = (double*) calloc (Ni   , sizeof(double));
    y   = (double*) calloc (Nj   , sizeof(double));
    u0  = (double*) calloc (Ni*Nj, sizeof(double));
    u1  = (double*) calloc (Ni*Nj, sizeof(double));
    u2  = (double*) calloc (Ni*Nj, sizeof(double));
    u3  = (double*) calloc (Ni*Nj, sizeof(double));

    double u_inf = 0.5;
    double v_inf = 0.0;
    double b = u_inf;
    double x0, y0;

    {
      /* Initial solution */
      x0 = 5.0, y0 = 5.0;
      for (i = 0; i < Ni; i++){
        for (j = 0; j < Nj; j++){
          x[i] = i*dx;
          y[j] = j*dy;
          int p = Nj*i + j;

          double rx, ry;
          rx = (x[i] - x0);
          ry = (y[j] - y0);
          if (rx < -5)      { rx += 10; }
          else if (rx > 5)  { rx -= 10; }

          double rsq = rx*rx + ry*ry;
          double rho, u, v, P;
          double du, dv;
          rho = power(1.0 - ((GAMMA-1.0)*b*b)/(8.0*GAMMA*pi*pi) * exp(1.0-rsq), 1.0/(GAMMA-1.0));
          P   = power(rho,GAMMA);
          du  = - b/(2.0*pi) * exp(0.5*(1.0-rsq)) * ry;
          dv  =   b/(2.0*pi) * exp(0.5*(1.0-rsq)) * rx;
          u   = u_inf + du;
          v   = v_inf + dv;
          u0[p] = rho;
          u1[p] = rho*u;
          u2[p] = rho*v;
          u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
        }
      }

      char fname[_MAX_STRING_SIZE_] = "initial";
      if (nsims > 1) {
        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(ns, index, (int)log10(nsims)+1);
        strcat(fname, "_");
        strcat(fname, index);
      }
      strcat(fname, ".inp");
      if (!strcmp(ip_file_type,"ascii")) {
        printf("Writing ASCII initial solution file %s\n", fname);
        out = fopen(fname,"w");
        for (i = 0; i < Ni; i++)  fprintf(out,"%lf ",x[i]);
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  fprintf(out,"%lf ",y[j]);
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u0[p]);
          }
        }
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u1[p]);
          }
        }
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u2[p]);
          }
        }
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u3[p]);
          }
        }
        fprintf(out,"\n");
        fclose(out);
      } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
        printf("Writing binary initial solution file %s\n", fname);
        out = fopen(fname,"wb");
        fwrite(x,sizeof(double),Ni,out);
        fwrite(y,sizeof(double),Nj,out);
        double *U = (double*) calloc (4*Ni*Nj,sizeof(double));
        for (i=0; i < Ni; i++) {
          for (j = 0; j < Nj; j++) {
            int p = Nj*i + j;
            int q = Ni*j + i;
            U[4*q+0] = u0[p];
            U[4*q+1] = u1[p];
            U[4*q+2] = u2[p];
            U[4*q+3] = u3[p];
          }
        }
        fwrite(U,sizeof(double),4*Ni*Nj,out);
        free(U);
        fclose(out);
      }
    }


    {
      /* Exact solution */
      x0 = 5.0+tf*u_inf, y0 = 5.0;
      if (x0 > 10) x0 -= 10; //periodic domain
      printf("Final time: %lf, Vortex center: %lf, %lf\n",tff,x0,y0);
      for (i = 0; i < Ni; i++){
        for (j = 0; j < Nj; j++){
          x[i] = i*dx;
          y[j] = j*dy;
          int p = Nj*i + j;

          double rx, ry;
          rx = (x[i] - x0);
          ry = (y[j] - y0);
          if (rx < -5)      { rx += 10; }
          else if (rx > 5)  { rx -= 10; }

          double rsq = rx*rx + ry*ry;
          double rho, u, v, P;
          double du, dv;
          rho = power(1.0 - ((GAMMA-1.0)*b*b)/(8.0*GAMMA*pi*pi) * exp(1.0-rsq), 1.0/(GAMMA-1.0));
          P   = power(rho,GAMMA);
          du  = - b/(2.0*pi) * exp(0.5*(1.0-rsq)) * ry;
          dv  =   b/(2.0*pi) * exp(0.5*(1.0-rsq)) * rx;
          u   = u_inf + du;
          v   = v_inf + dv;
          u0[p] = rho;
          u1[p] = rho*u;
          u2[p] = rho*v;
          u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
        }
      }

      char fname[_MAX_STRING_SIZE_] = "exact";
      if (nsims > 1) {
        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(ns, index, (int)log10(nsims)+1);
        strcat(fname, "_");
        strcat(fname, index);
      }
      strcat(fname, ".inp");
      if (!strcmp(ip_file_type,"ascii")) {
        printf("Writing ASCII exact solution file %s\n", fname);
        out = fopen(fname,"w");
        for (i = 0; i < Ni; i++)  fprintf(out,"%lf ",x[i]);
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  fprintf(out,"%lf ",y[j]);
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u0[p]);
          }
        }
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u1[p]);
          }
        }
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u2[p]);
          }
        }
        fprintf(out,"\n");
        for (j = 0; j < Nj; j++)  {
          for (i = 0; i < Ni; i++)  {
            int p = Nj*i + j;
            fprintf(out,"%lf ",u3[p]);
          }
        }
        fprintf(out,"\n");
        fclose(out);
      } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
        printf("Writing binary exact solution file %s\n", fname);
        out = fopen(fname,"wb");
        fwrite(x,sizeof(double),Ni,out);
        fwrite(y,sizeof(double),Nj,out);
        double *U = (double*) calloc (4*Ni*Nj,sizeof(double));
        for (i=0; i < Ni; i++) {
          for (j = 0; j < Nj; j++) {
            int p = Nj*i + j;
            int q = Ni*j + i;
            U[4*q+0] = u0[p];
            U[4*q+1] = u1[p];
            U[4*q+2] = u2[p];
            U[4*q+3] = u3[p];
          }
        }
        fwrite(U,sizeof(double),4*Ni*Nj,out);
        free(U);
        fclose(out);
      }
    }

    free(x);
    free(y);
    free(u0);
    free(u1);
    free(u2);
    free(u3);

  }

  free(NI);
  free(NJ);

  return(0);
}
