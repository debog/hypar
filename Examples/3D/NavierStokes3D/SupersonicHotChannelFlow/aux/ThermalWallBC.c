/* 
 * This code generates the temperature field and writes it out
 * in the format that HyPar wants.
 *
 * It needs the solver.inp file.
 *
 * For parallel runs, the data is split in this code.
 *
 * NOTE: To compile this file, you need to include hypar/Extras
 * in the include path, i.e., 
 *   gcc -I/path/to/hypar/Extras/ ThermalWallBC.c -lm -o <exec_name>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <thermalwallbc_utils.h>

/*
 * Main function
*/
int main() 
{
  const double GAMMA = 1.4;
  int ierr, i, j, k;

  /* Ask for the BC spatial dimension */
  int bc_dim = 1; /* along y */

  /* This file generates a steady temperature data */
  int n_time_levels = 1; /* number of time levels is 1 */
  double *time_levels =(double*) calloc (n_time_levels, sizeof(double));
  time_levels[0] = 0.0; /* specify the temperature field at t=0 */

  int dims_global[3], ndims, nvars, iproc[3];
  char init_file_format[1000];

  /* read simulation parameters */
  ierr = ReadSolverInp(dims_global, iproc, &ndims, &nvars, init_file_format);
  if (ierr) return(ierr);

  if (ndims != 3) {
    printf("Error: ndims needs to be 3!\n");
    return(0);
  }
  printf("Grid dimensions: %d x %d x %d\n",dims_global[0],dims_global[1],dims_global[2]);
  printf("Number of processors: %d x %d x %d\n",iproc[0],iproc[1],iproc[2]);

  double  rho_inf = 1.0,
          p_inf = 1.0/GAMMA,
          T_inf = p_inf/rho_inf;

  double  T_cold = T_inf,
          T_hot = 2*T_inf;

  double *x = (double*) calloc (dims_global[0],sizeof(double));
  double *y = (double*) calloc (dims_global[1],sizeof(double));
  double *z = (double*) calloc (dims_global[2],sizeof(double));
  ierr = ReadGrid(3, dims_global, x, y, z, init_file_format);
  if (ierr) return(ierr);

  /* size of temperature data */
  dims_global[bc_dim] = n_time_levels;
  int data_size = dims_global[0]*dims_global[1]*dims_global[2];

  /* generate the temperature field */
  double *T_global = (double*) calloc (data_size,sizeof(double));

  /* In this example, a hot temperature is specified for x > 0.5 */
  double  xc = 0.2;

  for (i = 0; i < dims_global[0]; i++) {
    for (j = 0; j < dims_global[1]; j++) { 
      for (k = 0; k < dims_global[2]; k++) {

        int p = i + dims_global[0]*j + dims_global[0]*dims_global[1]*k;
        if (x[i] < xc)  T_global[p] = T_cold;
        else            T_global[p] = T_hot;

      }
    }
  }

  /* write y-min BC data */
  ierr = WriteToFile(ndims, T_global, iproc, dims_global, bc_dim, 1, "temperature_ymin_bc.dat", time_levels);
  if (ierr) return(ierr);
  /* write y-max BC data */
  ierr = WriteToFile(ndims, T_global, iproc, dims_global, bc_dim, -1, "temperature_ymax_bc.dat", time_levels);
  if (ierr) return(ierr);

  free(T_global);  
  free(x);
  free(y);
  free(z);
  free(time_levels);
  return(0);
}

