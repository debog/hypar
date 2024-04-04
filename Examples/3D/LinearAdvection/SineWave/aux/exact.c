#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{
  double pi = 4.0*atan(1.0);

  int NI, NJ, NK, ndims, n_iter;
  char ip_file_type[50];
  double tf, dt;
  FILE *in;

  strcpy(ip_file_type,"ascii");
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) {
          fscanf(in,"%d",&ndims);
        } else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
          fscanf(in,"%d",&NK);
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

  if (ndims != 3) {
    printf("ndims is not 3 in solver.inp. this code is to generate 3D initial conditions\n");
    return(0);
  }
  printf("Grid:\t\t\t%d X %d X %d\n", NI, NJ, NK);

  double ax, ay, az;
  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) {
    fprintf(stderr,"Error: Input file \"physics.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "advection")) {
          fscanf(in,"%lf",&ax);
          fscanf(in,"%lf",&ay);
          fscanf(in,"%lf",&az);
        }
      }
    } else {
      fprintf(stderr,"Error: Illegal format in physics.inp. Crash and burn!\n");
      return(0);
    }
  }
  fclose(in);
  printf("Advection: %3.1f, %3.1f, %3.1f\n", ax, ay, az);

  tf = (double)n_iter * dt;
  printf("dt: %lf, n_iter: %d, Final time: %lf\n", dt, n_iter, tf);

  double Lx, Ly, Lz;
  Lx = 1.0;
  Ly = 1.0;
  Lz = 1.0;

  int i,j,k;
  double dx = Lx / ((double)NI);
  double dy = Ly / ((double)NJ);
  double dz = Lz / ((double)NK);

  double *x, *y, *z, *u;
  x = (double*) calloc (NI      , sizeof(double));
  y = (double*) calloc (NJ      , sizeof(double));
  z = (double*) calloc (NK      , sizeof(double));
  u = (double*) calloc (NI*NJ*NK, sizeof(double));

  for (i = 0; i < NI; i++)  x[i] = (i+0.5)*dx;
   for (j = 0; j < NJ; j++)  y[j] = (j+0.5)*dy;
  for (k = 0; k < NK; k++)  z[k] = (k+0.5)*dz;

  {
    char filename[50] = "initial.inp";

    for (i = 0; i < NI; i++){
      for (j = 0; j < NJ; j++){
        for (k = 0; k < NK; k++) {
          int p = i + NI*j + NI*NJ*k;
          u[p] = sin(2*pi*x[i]) * sin(2*pi*y[j]) * sin(2*pi*z[k]);
        }
      }
    }

    FILE *out;
    if (!strcmp(ip_file_type,"ascii")) {
      fprintf(stderr, "Sorry, binary files only! Set \"ip_file_type\" to binary in solver.inp.\n");
    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
      printf("Writing binary initial solution file %s\n", filename);
      out = fopen(filename,"wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(y,sizeof(double),NJ,out);
      fwrite(z,sizeof(double),NK,out);
      fwrite(u,sizeof(double),NI*NJ*NK,out);
      fclose(out);
    }

  }

  {
    char filename[50] = "exact.inp";

    for (i = 0; i < NI; i++){
      for (j = 0; j < NJ; j++){
        for (k = 0; k < NK; k++) {
          int p = i + NI*j + NI*NJ*k;
          u[p] = sin(2*pi*(x[i]-ax*tf)) * sin(2*pi*(y[j]-ay*tf)) * sin(2*pi*(z[k]-az*tf));
        }
      }
    }

    FILE *out;
    if (!strcmp(ip_file_type,"ascii")) {
      fprintf(stderr, "Sorry, binary files only! Set \"ip_file_type\" to binary in solver.inp.\n");
    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
      printf("Writing binary initial solution file %s\n", filename);
      out = fopen(filename,"wb");
      fwrite(x,sizeof(double),NI,out);
      fwrite(y,sizeof(double),NJ,out);
      fwrite(z,sizeof(double),NK,out);
      fwrite(u,sizeof(double),NI*NJ*NK,out);
      fclose(out);
    }

  }

  free(x);
  free(y);
  free(z);
  free(u);

  return(0);
}
