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

int main()
{

  const double pi = 4.0*atan(1.0);

	int *NI, ns, ndims, nsims, niter;
  char ip_file_type[50];
  double adv, dt;
  FILE *in;
  
  /* default values */
  nsims = 1;
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
  printf("Number of simulation domains: %d\n", nsims);
  NI = (int*) calloc (nsims, sizeof(int));

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
          for (ns = 0; ns < nsims; ns++) {
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

  /* read physics.inp */
  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) {
    printf("Error: Input file \"physics.inp\" not found. Default values will be used.\n");
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "advection")) {
          fscanf(in,"%lf",&adv);
        }
      }
    } else printf("Error: Illegal format in physics.inp. Crash and burn!\n");
  }
  fclose(in);

  /* some checks */
  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n");
    return(0);
  }
  if (nsims == 1) {
    printf("Grid:\t\t\t%d\n",NI[0]);
  } else {
    printf("Grid sizes:\n");
    for (ns = 0; ns < nsims; ns++) printf("\t%d\n", NI[ns]);
  }

  printf("Advection speed: %lf\n", adv);

  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  for (ns = 0; ns < nsims; ns++) {

    int N = NI[ns];

    char fname[_MAX_STRING_SIZE_] = "initial";
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(ns, index, (int)log10(nsims)+1);
      strcat(fname, "_");
      strcat(fname, index);
    }
    strcat(fname, ".inp");

  	int i;
  	double dx = 1.0 / ((double)N);
  
  	double *x, *u;
  	x = (double*) calloc (N, sizeof(double));
  	u = (double*) calloc (N, sizeof(double));
  
  	for (i = 0; i < N; i++){
  		x[i] = i*dx;
  		u[i] = sin(2*pi*x[i]);
  	}
  
    FILE *out;
  
    if (!strcmp(ip_file_type,"ascii")) {
      printf("Writing ASCII initial solution file %s\n", fname);
    	out = fopen(fname,"w");
      for (i = 0; i < N; i++)  fprintf(out,"%lf ",x[i]);
      fprintf(out,"\n");
  	  for (i = 0; i < N; i++)	fprintf(out,"%lf ",u[i]);						
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

  for (ns = 0; ns < nsims; ns++) {

    int N = NI[ns];

    char fname[_MAX_STRING_SIZE_] = "exact";
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(ns, index, (int)log10(nsims)+1);
      strcat(fname, "_");
      strcat(fname, index);
    }
    strcat(fname, ".inp");

  	int i;
  	double dx = 1.0 / ((double)N);
  
  	double *x, *u;
  	x = (double*) calloc (N, sizeof(double));
  	u = (double*) calloc (N, sizeof(double));
  
  	for (i = 0; i < N; i++){
  		x[i] = i*dx;
  		u[i] = sin(2*pi*(x[i]-adv*tf));
  	}
  
    FILE *out;
  
    if (!strcmp(ip_file_type,"ascii")) {
      printf("Writing ASCII exact solution file %s\n", fname);
    	out = fopen(fname,"w");
      for (i = 0; i < N; i++)  fprintf(out,"%lf ",x[i]);
      fprintf(out,"\n");
  	  for (i = 0; i < N; i++)	fprintf(out,"%lf ",u[i]);						
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
  }

  free(NI);

	return(0);
}
