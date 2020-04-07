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

	int *NI, ns, ndims, nsims;
  FILE *in;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

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
        }
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
    
  char fname[_MAX_STRING_SIZE_] = "initial";
  char adv_fname_root[_MAX_STRING_SIZE_] = "none";
  int flag = 0;

  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) {
    printf("Error: Input file \"physics.inp\" not found.\n");
    return 0;
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "advection_filename")) {
          fscanf(in,"%s", adv_fname_root);
          flag = 1;
        }
      }
    } else {
      printf("Error: Illegal format in solver.inp. Crash and burn!\n");
    }
  }
  fclose(in);

  if (flag == 0) {
    printf("Error: Input \"advection_filename\" not found in physics.inp file.\n");
    return 0;
  }

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

  for (ns = 0; ns < nsims; ns++) {

    int N = NI[ns];

    char fname[_MAX_STRING_SIZE_] = "initial";
    char adv_filename[_MAX_STRING_SIZE_];
    strcpy(adv_filename, adv_fname_root);
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(ns, index, (int)log10(nsims)+1);
      strcat(fname, "_");
      strcat(fname, index);
      strcat(adv_filename, "_");
      strcat(adv_filename, index);
    }
    strcat(fname, ".inp");
    strcat(adv_filename, ".inp");

  	int i;
  	double dx = 1.0 / ((double)N);
  
  	double *x, *u, *a;
  	x = (double*) calloc (N, sizeof(double));
  	u = (double*) calloc (N, sizeof(double));
  	a = (double*) calloc (N, sizeof(double));
  
  	for (i = 0; i < N; i++) {
  		x[i] = i*dx;
  		u[i] = sin(2*pi*x[i]);
  		a[i] = 1.0 + 0.5*sin(2*pi*x[i])*sin(2*pi*x[i]);
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
      printf("Writing ASCII advection field file %s\n", adv_filename);
  	  out = fopen(adv_filename,"w");
      for (i = 0; i < N; i++)  fprintf(out,"%1.16E ",x[i]);
      fprintf(out,"\n");
  	  for (i = 0; i < N; i++)	fprintf(out,"%1.16E ",a[i]);						
      fprintf(out,"\n");
  	  fclose(out);
    } else if (     (!strcmp(ip_file_type,"binary")) 
                ||  (!strcmp(ip_file_type,"bin"))     ) {
      printf("Writing binary initial solution file %s\n", fname);
      out = fopen(fname,"wb");
      fwrite(x,sizeof(double),N,out);
      fwrite(u,sizeof(double),N,out);
      fclose(out);
      printf("Writing binary advection field file %s\n", adv_filename);
  	  out = fopen(adv_filename,"wb");
      fwrite(x,sizeof(double),N,out);
      fwrite(a,sizeof(double),N,out);
  	  fclose(out);
    }
  
  	free(x);
  	free(u);
  	free(a);
  }

	return 0;
}
