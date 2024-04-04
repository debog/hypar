#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * This function must be consistent with the same function
 * in hypar/src/MPIFunctions/MPIPartition1D.c
 */
int MPIPartition1D(int nglobal,int nproc,int rank)
{
  int nlocal;
  if (nglobal%nproc == 0) nlocal = nglobal/nproc;
  else {
    if (rank == nproc-1)  nlocal = nglobal/nproc + nglobal%nproc;
    else                  nlocal = nglobal/nproc;
  }
  return(nlocal);
}

/*
 * This function must be consistent with the same function
 * in hypar/src/MPIFunctions/MPILocalDomainLimits.c
 */
int MPILocalDomainLimits(int ndims,int* ip,int *iproc,int *dim_global,int *is, int *ie)
{
  int i;
  for (i=0; i<ndims; i++) {
    int imax_local, isize, root = 0;
    imax_local = MPIPartition1D(dim_global[i],iproc[i],root );
    isize      = MPIPartition1D(dim_global[i],iproc[i],ip[i]);
    if (is)  is[i] = ip[i]*imax_local;
    if (ie)  ie[i] = ip[i]*imax_local + isize;
  }
  return(0);
}

/* Read solver.inp for basic simulation inputs like grid size and number of MPI ranks */
int ReadSolverInp(int *a_dims_global, int *a_iproc, int *a_ndims, int *a_nvars, char *ip_format)
{
  FILE *in;
  in = fopen("solver.inp","r");

  if (!in) {
    printf("Error: solver.inp not found!\n");
    return(1);
  }

  char word[500];
  fscanf(in,"%s",word);
  if (!strcmp(word,"begin")) {

    while (strcmp(word,"end")) {

      fscanf(in,"%s",word);

      if (!strcmp(word,"ndims")) {
        fscanf(in,"%d", a_ndims);
      } else if (!strcmp(word,"nvars")) {
        fscanf(in,"%d", a_nvars);
      } else if (!strcmp(word,"size")) {
        fscanf(in,"%d", &a_dims_global[0]);
        fscanf(in,"%d", &a_dims_global[1]);
        fscanf(in,"%d", &a_dims_global[2]);
      } else if (!strcmp(word,"iproc")) {
        fscanf(in,"%d", &a_iproc[0]);
        fscanf(in,"%d", &a_iproc[1]);
        fscanf(in,"%d", &a_iproc[2]);
      } else if (!strcmp(word, "ip_file_type")) {
        fscanf(in,"%s",ip_format);
      }
    }

  } else {

    printf("Error: illegal format in solver.inp.\n");
    return(1);

  }

  fclose(in);
  return(0);
}

/* Read physics.inp for the universal gas constant, reference density, and reference pressure */
int ReadPhysicsInp(double *a_R, double *a_rho_ref, double *a_p_ref)
{
  printf("Reading file \"physics.inp\"...\n");
  FILE *in;
  in = fopen("physics.inp","r");

  if (!in) {

    printf("Error: Input file \"physics.inp\" not found.\n");
    return(1);

  } else {

    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if      (!strcmp(word, "rho_ref")) fscanf(in,"%lf",a_rho_ref);
        else if (!strcmp(word, "p_ref"  )) fscanf(in,"%lf",a_p_ref  );
        else if (!strcmp(word, "R"      )) fscanf(in,"%lf",a_R      );
      }
    } else printf("Error: Illegal format in physics.inp. Crash and burn!\n");
  }
  fclose(in);
  return(0);
}

int ReadGrid(
              int     ndims,          /*!< Number of spatial dimensions */
              int     *dim,           /*!< Integer array of size ndims with global grid size in each dimension */
              double  *x,             /*!< Grid associated with the array */
              double  *y,             /*!< Grid associated with the array */
              double  *z,             /*!< Grid associated with the array */
              char    *ip_file_type   /*!< Initial solution file type (ASCII or binary) */
            )
{
  int i, d, ferr, size,offset;

  /* allocate grid array */
  size = 0;
  for (d=0; d<ndims; d++) size += dim[d];
  double *xg = (double*) calloc(size,sizeof(double));

  if (!strcmp(ip_file_type,"ascii")) {

    FILE *in;
    in = fopen("initial.inp","r");

    if (!in) {

      printf("ERROR: initial.inp not found!\n");
      return(1);

    } else {

      printf("Reading from ASCII file initial.inp.\n");

      /* read grid */
      int offset = 0;
      for (d = 0; d < ndims; d++) {
        for (i = 0; i < dim[d]; i++) ferr = fscanf(in,"%lf",&xg[i+offset]);
        offset += dim[d];
      }

      fclose(in);

    }

  } else if ((!strcmp(ip_file_type,"bin")) || (!strcmp(ip_file_type,"binary"))) {

    FILE *in;
    in = fopen("initial.inp","rb");

    if (!in) {

      printf("ERROR: initial.inp not found!\n");
      return(1);

    } else {

      printf("Reading array from binary file initial.inp.\n");
      size_t bytes;
      int size;

      size = 0; for (d=0; d<ndims; d++) size += dim[d];
      xg = (double*) calloc(size,sizeof(double));

      /* read grid */
      size = 0;
      for (d = 0; d < ndims; d++) size += dim[d];
      bytes = fread(xg, sizeof(double), size, in);
      if ((int)bytes != size) {
        fprintf(stderr,"Error in ReadArray(): Unable to read grid. Expected %d, Read %d.\n",
                size, (int)bytes);
      }

      fclose(in);
    }

  }


  /* copy to x,y,z */
  for (i = 0; i < dim[0]; i++) x[i] = xg[i];
  for (i = 0; i < dim[1]; i++) y[i] = xg[i+dim[0]];
  for (i = 0; i < dim[2]; i++) z[i] = xg[i+dim[0]+dim[1]];

  free(xg);
  return(0);
}

/* This function takes in the temperature data as a global array, the number of MPI ranks along each spatial dimension,
 * the global dimension along each spatial dimension, the dimension along which the thermal wall BC is being applied,
 * and the face at which it is being applied (i.e. the min face (lo end) or the max face (hi end)).
 *
 * It then partitions the global data into subdomains for each participating MPI ranks and writes them out to
 * a binary file. The domain decomposition is consistent with that in HyPar, and the file is written in a format
 * that HyPar wants.
 *
 * Note that all MPI ranks are not participants; only those abutting the boundary face at which the thermal wall BC
 * is being applied need this data.
 *
 * Along all dimensions that are not a_bc_dim, a_dims_global[i] must be equal to the global grid size of the
 * HyPar simulation along that dimension.
 *
 * If the temperature field is steady, i.e., the data is available for only one time level (t=0),
 *  a_dims_global[a_bc_dim] must be 1.
 *
 * If the temperature field is unsteady, HyPar will want to read in the stacked data for all the time levels,
 *  a_dims_global[a_bc_dim] must be the number of time levels
 * During the simulation, HyPar will apply the temperature data from the index (along a_bc_dim) corresponding to
 * the highest time level that is <= to the current simulation/stage time.
 *
 * For example, if a_bc_dim = 1 (the second spatial dimension), HyPar implements:
 *  T[i][j_boundary][k](t) = T_input[i][t_index][k],
 * where t_index is the highest integer in [0,a_dims_global[a_bc_dim]-1] satisfying
 *  a_time_levels[t_index] <= t
*/
int  WriteToFile( int     ndims,          /* number of spatial dimensions */
                  double  *a_T_global,    /* global array of the temperature data */
                  int     *a_iproc,       /* int array with the number of MPI ranks along each spatial dimension */
                  int     *a_dims_global, /* int array with global dimension along each spatial dimension */
                  int     a_bc_dim,       /* dimension along which the thermal BC is being applied */
                  int     a_bc_dir,       /* face at which the thermal BC is being applied: 1 -> min face, -1 -> max face */
                  char    *a_filename,    /* filename to which to write the BC data */
                  double  *a_time_levels  /* the simulation time for each temperature data; this array must
                                             be of size a_dims_global[a_bc_dim] */
                )
{
  int ferr;
  FILE *out;
  out = fopen(a_filename,"wb");

  int ip[3];

  for (ip[0] = 0; ip[0] < a_iproc[0]; ip[0]++) {
    for (ip[1] = 0; ip[1] < a_iproc[1]; ip[1]++) {
      for (ip[2] = 0; ip[2] < a_iproc[2]; ip[2]++) {

        if (    ( (a_bc_dir ==  1) && (ip[a_bc_dim] == 0) )
            ||  ( (a_bc_dir == -1) && (ip[a_bc_dim] == a_iproc[a_bc_dim]-1) )  ) {

          int rank[ndims];
          rank[0] = ip[0];
          rank[1] = ip[1];
          rank[2] = ip[2];

          int rank_fake[ndims];
          rank_fake[0] = rank[0];
          rank_fake[1] = rank[1];
          rank_fake[2] = rank[2];
          rank_fake[a_bc_dim] = 0;

          int nproc[ndims];
          nproc[0] = a_iproc[0];
          nproc[1] = a_iproc[1];
          nproc[2] = a_iproc[2];
          nproc[a_bc_dim] = 1;

          ferr = fwrite(rank,sizeof(int),ndims,out);
          if (ferr != ndims) {
            printf("Error: (%d,%d,%d) Unable to write boundary data file (rank).\n",rank[0],rank[1],rank[2]);
            return(1);
          }
          int size[ndims], is[ndims], ie[ndims];
          MPILocalDomainLimits(ndims,rank_fake,nproc,a_dims_global,is,ie);
          size[0] = ie[0] - is[0];
          size[1] = ie[1] - is[1];
          size[2] = ie[2] - is[2];
          if (size[a_bc_dim] != a_dims_global[a_bc_dim]) {
            printf("Error: Something went wrong. size[a_bc_dim] is %d, but a_dims_global[a_bc_dim] is %d.\n",
                    size[a_bc_dim], a_dims_global[a_bc_dim]);
            return(1);
          }
          ferr = fwrite(size,sizeof(int),ndims,out);
          if (ferr != ndims) {
            printf("Error: (%d,%d,%d) Unable to write boundary data file (size).\n",rank[0],rank[1],rank[2]);
            return(1);
          }

          ferr = fwrite(a_time_levels,sizeof(double),size[a_bc_dim],out);
          if (ferr != size[a_bc_dim]) {
            printf("Error: (%d,%d,%d) Unable to write boundary data file (time levels).\n",rank[0],rank[1],rank[2]);
            return(1);
          }

          int local_size = size[0] * size[1] * size[2];
          double *T_local = (double*) calloc (local_size,sizeof(double));
          int i,j,k;
          for (i = 0; i < size[0]; i++) {
            for (j = 0; j < size[1]; j++) {
              for (k = 0; k < size[2]; k++) {
                int p1 = i + size[0]*j + size[0]*size[1]*k;
                int p2 =    (i+is[0])
                          + a_dims_global[0]*(j+is[1])
                          + a_dims_global[0]*a_dims_global[1]*(k+is[2]);
                T_local[p1] = a_T_global[p2];
              }
            }
          }
          ferr = fwrite(T_local,sizeof(double),local_size,out);
          if (ferr != local_size) {
            printf("Error: (%d,%d,%d) Unable to write boundary data file (T_local) (%d, expected %d).\n",
                   rank[0],rank[1],rank[2],ferr,local_size);
            return(1);
          }
          free(T_local);

        }
      }
    }
  }
  fclose(out);

  return(0);
}

