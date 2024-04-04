/* This code extracts the maximum E-field for the Landau
 * damping case from the electric field files efield_<nnnnn>.dat,
 * which are text files with three columns (grid index, x-coordinate,
 * and the e-field value). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void IncrementFilename(char *f)
{
  if (f[11] == '9') {
    f[11] = '0';
    if (f[10] == '9') {
      f[10] = '0';
      if (f[9] == '9') {
        f[9] = '0';
        if (f[8] == '9') {
          f[8] = '0';
          if (f[7] == '9') {
            f[7] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[7]++;
          }
        } else {
          f[8]++;
        }
      } else {
        f[9]++;
      }
    } else {
      f[10]++;
    }
  } else {
    f[11]++;
  }
}

double absVal(double x)
{
  return (x < 0 ? -x : x);
}

int readEField( const char* const a_fname,
                const int         a_size,
                double* const     a_E_max )
{
  FILE* in = fopen(a_fname,"r");

  if (!in) {

    printf("File %s not found.\n", a_fname);
    return 1;

  } else {

    printf("Reading file %s.\n", a_fname);

    int i, index;
    double xval, eval;

    *a_E_max = 0.0;
    for (i = 0; i < a_size; i++) {
      fscanf(in, "%d %lf %lf", &index, &xval, &eval);
      if (index != i) {
        printf("Error! Index mismatch in file %s.\n", a_fname);
        return 1;
      }
      if (absVal(eval) > *a_E_max) {
        *a_E_max = absVal(eval);
      }
    }

    fclose(in);
  }

  return 0;
}

int main()
{
  int NI,NJ,ndims;
  double dt, samay;
  int file_op_iter, restart_iter=0;
  char filename[50];

  FILE* in = fopen("solver.inp","r");
  if (!in) {

    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return 1;

  } else {

    char word[100];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
         } else if (!strcmp(word, "dt"))  fscanf(in,"%lf",&dt);
         else if (!strcmp(word, "file_op_iter")) fscanf(in,"%d",&file_op_iter);
         else if (!strcmp(word, "restart_iter")) fscanf(in,"%d",&restart_iter);
      }
    }
    fclose(in);

  }

  samay = (double) restart_iter*dt;
  FILE* out;
  if (!restart_iter) {
    printf("Writing to efield_history.dat (new file).\n");
    out = fopen("efield_history.dat","w");
  } else {
    printf("Writing to efield_history.dat (append).\n");
    out = fopen("efield_history.dat","a");
  }
  strcpy(filename,"efield_00000.dat");
  {
    int t;
    for (t=0; t<restart_iter; t++) {
      if ((t+1)%file_op_iter == 0) IncrementFilename(filename);
    }
  }

  while(1) {

    double maxE = 0.0;
    int status = readEField( filename, NI, &maxE );

    if (!status) {
      fprintf(out, "%1.16e  %1.16e\n", samay, maxE);
    } else {
      break;
    }

    samay += (dt * (double)file_op_iter);
    IncrementFilename(filename);
  }

  fclose(out);
  return 0;
}
