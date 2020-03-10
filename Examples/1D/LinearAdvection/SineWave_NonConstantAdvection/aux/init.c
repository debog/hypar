#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main()
{
  const double pi = 4.0*atan(1.0);

	int NI, ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return 0;
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) {
          fscanf(in,"%d",&ndims);
        } else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in,"%s",ip_file_type);
        }
      }
    } else {
      printf("Error: Illegal format in solver.inp. Crash and burn!\n");
    }
  }
  fclose(in);

  char adv_filename[500] = "none";
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
          fscanf(in,"%s", adv_filename);
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
	printf("Grid:\t\t\t%d\n",NI);

	int i;
	double dx = 1.0 / ((double)NI);

	double *x, *u, *a;
	x = (double*) calloc (NI, sizeof(double));
	u = (double*) calloc (NI, sizeof(double));
	a = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
		u[i] = sin(2*pi*x[i]);
		a[i] = 1.0+0.5*sin(2*pi*x[i])*sin(2*pi*x[i]);
	}

  FILE *out;
  strcat(adv_filename, ".inp");

  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
  	out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%lf ",u[i]);						
    fprintf(out,"\n");
    fclose(out);
    printf("Writing ASCII advection field file %s\n", adv_filename);
	  out = fopen(adv_filename,"w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",a[i]);						
    fprintf(out,"\n");
	  fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary initial solution file initial.inp\n");
    out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(u,sizeof(double),NI,out);
    fclose(out);
    printf("Writing binary advection field file %s\n", adv_filename);
	  out = fopen(adv_filename,"wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(a,sizeof(double),NI,out);
	  fclose(out);
  }

	free(x);
	free(u);
	free(a);

	return(0);
}
