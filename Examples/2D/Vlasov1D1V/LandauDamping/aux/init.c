/* Reference: Finn et al., A Numerical Study of Landau Damping with PETSc-PIC,
 *            https://doi.org/10.48550/arXiv.2303.12620 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main()
{
  double pi = 4.0 * atan(1.0);
  double two_pi = 8.0 * atan(1.0);

  double k = 1.0;
  double L = two_pi / k;
  double v_th = 1.0;
  double alpha = 0.01;
  double vmax = 10.0*v_th;

  int NI,NJ,ndims,n_iter;
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

  int i,j;
  double dx = L / ((double)NI);
  double dv = (2.0*vmax) / ((double)NJ);
  double start_x = 0.0;
  double start_v = -vmax;

  double *x, *v, *f;
  x = (double*) calloc (NI   , sizeof(double));
  v = (double*) calloc (NJ   , sizeof(double));
  f = (double*) calloc (NI*NJ, sizeof(double));

  for (i = 0; i < NI; i++){
    for (j = 0; j < NJ; j++){
      x[i] = start_x + i*dx;
      v[j] = start_v + j*dv;
      int p = NJ*i + j;
      double term1 = (1.0/sqrt(two_pi*v_th*v_th));
      double term2 = exp(-((v[j]*v[j])/(2*v_th*v_th)));
      double term3 = (1.0+alpha*cos(k*x[i]));
      f[p] = term1 * term2 * term3;
    }
  }

  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {

    out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",v[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  {
       for (i = 0; i < NI; i++)  {
          int p = NJ*i + j;
          fprintf(out,"%1.16e ", f[p]);
       }
    }
    fprintf(out,"\n");
    fclose(out);

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary exact solution file initial.inp\n");
    out = fopen("initial.inp","wb");
    fwrite(x, sizeof(double),NI,out);
    fwrite(v, sizeof(double),NJ,out);
    double *F = (double*) calloc (NI*NJ,sizeof(double));
    for (i=0; i < NI; i++) {
      for (j = 0; j < NJ; j++) {
        int p = NJ*i + j;
        int q = NI*j + i;
        F[q+0] = f[p];
      }
    }
    fwrite(F, sizeof(double),NI*NJ,out);
    free(F);
    fclose(out);

  }

  free(x);
  free(v);
  free(f);

  return(0);
}
