#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main()
{  
  double pi = 4.0*atan(1.0);

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
	double dx = 1.0 / ((double)NI);
	double dy = 1.0 / ((double)NJ);

  double *x, *y, *u;
	x = (double*) calloc (NI   , sizeof(double));
	y = (double*) calloc (NJ   , sizeof(double));
	u = (double*) calloc (NI*NJ, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      int p = NJ*i + j;
		  u[p] = sin(2*pi*x[i]) * sin(2*pi*y[j]) / (4*pi);
	  }
	}

  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {

    out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
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

	return(0);
}
