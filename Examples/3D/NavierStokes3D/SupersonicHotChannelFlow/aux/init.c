#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{

  const double pi     = 4.0*atan(1.0);
  const double GAMMA  = 1.4;

  int NI, NJ, NK, ndims;
  char ip_file_type[50] = "ascii";

  printf("Reading file \"solver.inp\"...\n");
  FILE *in;
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
        if (!strcmp(word, "ndims")) {
          fscanf(in, "%d", &ndims);
        } else if (!strcmp(word, "size")) {
          fscanf(in, "%d", &NI);
          fscanf(in, "%d", &NJ);
          fscanf(in, "%d", &NK);
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in, "%s", ip_file_type);
        }
      }
    } else {
      printf("Error: Illegal format in solver.inp. Crash and burn!\n");
      return(0);
    }
  }
  fclose(in);

  if (ndims != 3) {
    printf("Error: ndims is not 3 in solver.inp. this code is to generate 3D exact solution\n");
    return(0);
  }
  printf("Grid:  %d, %d, %d.\n", NI, NJ, NK);

  double Lx = 1.0,
         Ly = 0.2;

  double  dx = Lx / ((double)NI-1),
          dy = Ly / ((double)NJ-1),
          dz = (dx < dy ? dx : dy);

  double  rho_inf = 1.0,
          p_inf = 1.0/GAMMA,
          u_inf = 2.0,
          v_inf = 0.0,
          w_inf = 0.0;

  double *x, *y, *z, *U;

   x   = (double*) calloc (NI        , sizeof(double));
  y   = (double*) calloc (NJ        , sizeof(double));
  z   = (double*) calloc (NK        , sizeof(double));
  U   = (double*) calloc (5*NI*NJ*NK, sizeof(double));

  int i,j,k;
  for (i = 0; i < NI; i++){
    for (j = 0; j < NJ; j++){
      for (k = 0; k < NK; k++){

        x[i] = i*dx;
        y[j] = j*dy;
        z[k] = k*dz;

        int p = i + NI*j + NI*NJ*k;

        double rho, u, v, w, P;
        rho = rho_inf;
        P   = p_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;

        U[5*p+0] = rho;
        U[5*p+1] = rho*u;
        U[5*p+2] = rho*v;
        U[5*p+3] = rho*w;
        U[5*p+4] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
    }
  }

  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
    out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%1.16E ",y[j]);
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)  fprintf(out,"%1.16E ",z[k]);
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)  {
      for (j = 0; j < NJ; j++)  {
        for (i = 0; i < NI; i++)  {
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+0]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)  {
      for (j = 0; j < NJ; j++)  {
        for (i = 0; i < NI; i++)  {
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+1]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)  {
      for (j = 0; j < NJ; j++)  {
        for (i = 0; i < NI; i++)  {
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+2]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)  {
      for (j = 0; j < NJ; j++)  {
        for (i = 0; i < NI; i++)  {
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+3]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)  {
      for (j = 0; j < NJ; j++)  {
        for (i = 0; i < NI; i++)  {
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+4]);
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
    fwrite(U,sizeof(double),5*NI*NJ*NK,out);
    fclose(out);

  }

  free(x);
  free(y);
  free(z);
  free(U);

  return(0);
}
