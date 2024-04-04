#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#define _MAX_STRING_SIZE_ 200

void GetStringFromInteger(int a, char *A, int width)
{
  for (int i=0; i<width; i++) {
    char digit = (char) (a%10 + '0');
    a /= 10;
    A[width-1-i] = digit;
  }
  A[width] = 0;
  return;
}

typedef struct _params_
{
  double rhoL;
  double rhoR;
  double pL;
  double pR;
} Params;

int main()
{
  std::vector<double> rhoL = {1.0};
  std::vector<double> rhoR = {0.125};
  std::vector<double> pL = {1.0};
  std::vector<double> pR = {0.1};

  int n_rhoL = rhoL.size();
  int n_rhoR = rhoR.size();
  int n_pL = pL.size();
  int n_pR = pR.size();

  /* write out parameter set */
  {
    FILE* out;
    out = fopen("parameters.dat","w");
    fprintf(out, "begin\n");
    fprintf(out, "  rhoL ");
    for (int i = 0; i < n_rhoL; i++) {
      fprintf(out, "  %f", rhoL[i]);
    }
    fprintf(out, "\n");
    fprintf(out, "  rhoR ");
    for (int i = 0; i < n_rhoR; i++) {
      fprintf(out, "  %f", rhoR[i]);
    }
    fprintf(out, "\n");
    fprintf(out, "  pL ");
    for (int i = 0; i < n_pL; i++) {
      fprintf(out, "  %f", pL[i]);
    }
    fprintf(out, "\n");
    fprintf(out, "  pR ");
    for (int i = 0; i < n_pR; i++) {
      fprintf(out, "  %f", pR[i]);
    }
    fprintf(out, "\n");
    fprintf(out, "end");
    fclose(out);
  }

  int n_param_sets = n_rhoL * n_rhoR * n_pL * n_pR;
  Params *params = (Params*) calloc( n_param_sets, sizeof(Params));
  {
    int i, j;
    for (int i = 0; i < n_rhoL; i++) {
      for (int j = 0; j < n_rhoR; j++) {
        for (int k = 0; k < n_pL; k++) {
          for (int l = 0; l < n_pR; l++) {
            int s = l + k*n_pR + j*n_pR*n_pL + i*n_pR*n_pL*n_rhoR;
            params[s].rhoL = rhoL[i];
            params[s].rhoR = rhoR[j];
            params[s].pL = pL[k];
            params[s].pR = pR[l];
          }
        }
      }
    }
  }


  int nsims = 1;
  int* NI;
  int ndims=1;

  FILE *in;
  char ip_file_type[50];
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
  if (nsims > n_param_sets) {
    printf("Error: parameter set is too small.");
    return 0;
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
          for (int ns = 0; ns < nsims; ns++) {
            fscanf(in, "%d", (NI+ns));
          }
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in,"%s",ip_file_type);
        }
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
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
    for (int ns = 0; ns < nsims; ns++) printf("\t%d\n", NI[ns]);
  }


  for (int ns = 0; ns < nsims; ns++) {

    int N = NI[ns];

    printf( "(rhoL=%f, rhoR=%f, pL=%f, pR=%f)\n",
            params[ns].rhoL,
            params[ns].rhoR,
            params[ns].pL,
            params[ns].pR );

    double dx = 1.0 / ((double)(N-1));

    double *x, *u;
    x = (double*) calloc (N, sizeof(double));
    u = (double*) calloc (3*N, sizeof(double));

    for (int i = 0; i < N; i++){
      x[i] = i*dx;
      double RHO,U,P;
      if (x[i] < 0.5) {
        RHO = params[ns].rhoL;
        U   = 0.0;
        P   = params[ns].pL;
      } else {
        RHO = params[ns].rhoR;
        U   = 0;
        P   = params[ns].pR;
      }
      u[3*i+0] = RHO;
      u[3*i+1] = RHO*U;
      u[3*i+2] = P/0.4 + 0.5*RHO*U*U;
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
      printf("Error: Writing ASCII initial solution file not implemented. ");
      printf("Please choose ip_file_type in solver.inp as \"binary\".\n");
    } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
      FILE* out;
      printf("Writing binary initial solution file %s\n", fname);
      out = fopen(fname,"wb");
      fwrite(x,sizeof(double),N,out);
      fwrite(u,sizeof(double),3*N,out);
      fclose(out);
    }

    free(x);
    free(u);
  }

  free(params);

  return(0);
}
