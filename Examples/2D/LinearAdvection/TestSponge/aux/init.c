#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(){

  int NI,NJ,ndims;
  double pi = 4.0*atan(1.0);

  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims"))     fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
        }
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n");
    return(0);
  }
  printf("Grid:\t\t\t%d x %d\n",NI,NJ);

  int i,j;
  double dx = 12.0 / ((double)NI);
  double dy = 6.0 / ((double)NJ);

  double *x, *y, *u;
  x = (double*) calloc (NI   , sizeof(double));
  y = (double*) calloc (NJ   , sizeof(double));
  u = (double*) calloc (NI*NJ, sizeof(double));

  for (i = 0; i < NI; i++){
    for (j = 0; j < NJ; j++){
      x[i] = -6 + i*dx;
      y[j] = -3 + j*dy;
      int p = NJ*i + j;
      u[p] = exp(-(((x[i]+3)*(x[i]+3)/2+y[j]*y[j]/2)));
    }
  }

  FILE *out;
  out = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%lf ",y[j]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  {
    for (i = 0; i < NI; i++)  {
      int p = NJ*i + j;
      fprintf(out,"%lf ",u[p]);
    }
  }
  fprintf(out,"\n");
  fclose(out);

  free(x);
  free(y);
  free(u);

  return(0);
}
