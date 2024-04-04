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
  int ierr, i, j, k;

  /* Ask for the BC spatial dimension */
  int bc_dim, bc_dir;
  printf("Enter spatial dimension for this BC (0 - x, 1 - y, or 2 - z): ");
  scanf("%d", &bc_dim);
  printf("Enter face (1 - low side or -1 - high side: ");
  scanf("%d", &bc_dir);

  char filename[500] = "temperature_data.dat";

  /* This file generates a steady temperature data */
  int n_time_levels = 1; /* number of time levels is 1 */
  double *time_levels =(double*) calloc (n_time_levels, sizeof(double));
  time_levels[0] = 0.0; /* specify the temperature field at t=0 */

  int dims_global[3], ndims, nvars, iproc[3];
  char init_file_format[1000];

  /* read the simulation parameters */
  ierr = ReadSolverInp(dims_global, iproc, &ndims, &nvars, init_file_format);
  if (ierr) return(ierr);

  if (ndims != 3) {
    printf("Error: ndims needs to be 3!\n");
    return(0);
  }
  printf("Grid dimensions: %d x %d x %d\n",dims_global[0],dims_global[1],dims_global[2]);
  printf("Number of processors: %d x %d x %d\n",iproc[0],iproc[1],iproc[2]);

  double R, rho_ref, p_ref, T_ref;
  ierr = ReadPhysicsInp(&R, &rho_ref, &p_ref);
  if (ierr) return(ierr);
  T_ref = (p_ref/rho_ref)/R;

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

  /* In this example, a hot temperature is specified within a disk
   * of radius 100m, and a cold temperature is specified outside */

  double  rc = 100.0;
  double  xc = 500.0,
          yc = 500.0,
          zc = 500.0;

  double T_cold = T_ref;
  double T_hot = T_ref + 10.0; /* 10K hotter */

  for (i = 0; i < dims_global[0]; i++) {
    for (j = 0; j < dims_global[1]; j++) {
      for (k = 0; k < dims_global[2]; k++) {

        double r, rx, ry, rz;
        rx = x[i] - xc;
        ry = y[j] - yc;
        rz = z[k] - zc;
        if (bc_dim == 0)      r = sqrt(ry*ry+rz*rz);
        else if (bc_dim == 1) r = sqrt(rx*rx+rz*rz);
        else if (bc_dim == 2) r = sqrt(rx*rx+ry*ry);

        int p = i + dims_global[0]*j + dims_global[0]*dims_global[1]*k;
        if (r < rc) T_global[p] = R*T_hot;
        else        T_global[p] = R*T_cold;

      }
    }
  }

  ierr = WriteToFile(ndims, T_global, iproc, dims_global, bc_dim, bc_dir, filename, time_levels);
  if (ierr) return(ierr);

  free(T_global);
  free(x);
  free(y);
  free(z);
  free(time_levels);
  return(0);
}

