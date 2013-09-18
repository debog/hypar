#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int InitialSolution(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0,i,d;

  /* Only root process reads in initial solution file */
  double *ug,*xg,*dxinvg; /* global solution vector and grid arrays */
  if (!mpi->rank) {
    int size,offset;
    /* allocate global solution array */
    size = solver->npoints_global*solver->nvars;
    ug = (double*) calloc(size,sizeof(double));
    size = 0;
    for (d=0; d<solver->ndims; d++) size += solver->dim_global[d];
    xg      = (double*) calloc(size,sizeof(double));
    dxinvg  = (double*) calloc(size,sizeof(double));

    /* Reading grid and initial solution */
    printf("Reading grid and initial conditions from file \"initial.inp\".\n");
    FILE *in; in = fopen("initial.inp","r");
    if (!in) {
      fprintf(stderr,"Error: initial solution file \"initial.inp\" not found.\n");
      return(1);
    }

    /* read grid and calculate dxinv*/
    offset = 0;
    for (d = 0; d < solver->ndims; d++) {
      for (i = 0; i < solver->dim_global[d]; i++) ierr = fscanf(in,"%lf",&xg[i+offset]);
      for (i = 0; i < solver->dim_global[d]; i++) {
        if      (i == 0)                        dxinvg[i+offset] = 1.0/(xg[i+1+offset]-xg[i  +offset]);
        else if (i == solver->dim_global[d]-1)  dxinvg[i+offset] = 1.0/(xg[i  +offset]-xg[i-1+offset]);
        else                                    dxinvg[i+offset] = 2.0/(xg[i+1+offset]-xg[i-1+offset]);
      }
      offset += solver->dim_global[d];
    }

    /* read solution */
    for (i = 0; i < solver->nvars; i++) {
      int *index = solver->index;
      int done = 0; ierr = ArraySetValue_int(index,solver->ndims,0); CHECKERR(ierr);
      while (!done) {
        int p = ArrayIndex1D(solver->ndims,solver->dim_global,index,NULL,0);
        ierr = fscanf(in,"%lf",&ug[p*solver->nvars+i]);
        if (ierr != 1) return(1);
        done = ArrayIncrementIndex(solver->ndims,solver->dim_global,index);
      }
    }

    fclose(in);

  } else {
    ug      = NULL;
    xg      = NULL;
    dxinvg  = NULL;
  }

  ierr = MPIPartitionArraynD(solver->ndims,mpi,(mpi->rank?NULL:ug),solver->u,
                             solver->dim_global,solver->dim_local,
                             solver->ghosts,solver->nvars); CHECKERR(ierr);

  int offset_global, offset_local;
  offset_global = offset_local = 0;
  for (d=0; d<solver->ndims; d++) {
    ierr = MPIPartitionArray1D(mpi,(mpi->rank?NULL:&xg[offset_global]),&solver->x[offset_local],
                                    mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    ierr = MPIPartitionArray1D(mpi,(mpi->rank?NULL:&dxinvg[offset_global]),&solver->dxinv[offset_local],
                                    mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    offset_global += solver->dim_global[d];
    offset_local  += solver->dim_local [d];
  }

  if (!mpi->rank) {
    free(ug);
    free(xg);
    free(dxinvg);
  }
  return(0);
}