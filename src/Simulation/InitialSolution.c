/*! @file InitialSolution.c
    @author Debojyoti Ghosh
    @brief Read in initial solution from file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <io.h>
#include <simulation.h>

int VolumeIntegral(double*,double*,void*,void*);

/*! Read in initial solution from file, and compute grid spacing 
    and volume integral of the initial solution */
int InitialSolution ( void  *s,   /*!< Array of simulation objects of type #SimulationObject */
                      int   nsims /*!< Number of simulation objects */
                    )
{
  SimulationObject* simobj = (SimulationObject*) s;
  int n, flag, d, i, offset;
  _DECLARE_IERR_;

  for (n = 0; n < nsims; n++) {

    int ghosts = simobj[n].solver.ghosts;

    char fname_root[_MAX_STRING_SIZE_] = "initial";
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(n, index, (int)log10(nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }

    IERR ReadArray( simobj[n].solver.ndims,
                    simobj[n].solver.nvars,
                    simobj[n].solver.dim_global,
                    simobj[n].solver.dim_local,
                    simobj[n].solver.ghosts,
                    &(simobj[n].solver),
                    &(simobj[n].mpi),
                    simobj[n].solver.x,
                    simobj[n].solver.u,
                    fname_root,
                    &flag );
    if (!flag) {
      fprintf(stderr,"Error: initial solution file not found.\n");
      return(1);
    }
    CHECKERR(ierr);

    /* exchange MPI-boundary values of u between processors */
    MPIExchangeBoundariesnD(  simobj[n].solver.ndims,
                              simobj[n].solver.nvars,
                              simobj[n].solver.dim_local,
                              simobj[n].solver.ghosts,
                              &(simobj[n].mpi),
                              simobj[n].solver.u  );

    /* calculate dxinv */
    offset = 0;
    for (d = 0; d < simobj[n].solver.ndims; d++) {
      for (i = 0; i < simobj[n].solver.dim_local[d]; i++) {
        simobj[n].solver.dxinv[i+offset+ghosts] 
          = 2.0 / (simobj[n].solver.x[i+1+offset+ghosts]-simobj[n].solver.x[i-1+offset+ghosts]);
      }
      offset += (simobj[n].solver.dim_local[d] + 2*ghosts);
    }

    /* exchange MPI-boundary values of dxinv between processors */
    offset = 0;
    for (d = 0; d < simobj[n].solver.ndims; d++) {
      IERR MPIExchangeBoundaries1D( &(simobj[n].mpi),
                                    &(simobj[n].solver.dxinv[offset]),
                                    simobj[n].solver.dim_local[d],
                                    ghosts,
                                    d,
                                    simobj[n].solver.ndims ); CHECKERR(ierr);
      offset += (simobj[n].solver.dim_local[d] + 2*ghosts);
    }

    /* fill in ghost values of dxinv at physical boundaries by extrapolation */
    offset = 0;
    for (d = 0; d < simobj[n].solver.ndims; d++) {
      double *dxinv = &(simobj[n].solver.dxinv[offset]);
      int    *dim = simobj[n].solver.dim_local;
      if (simobj[n].mpi.ip[d] == 0) {
        /* fill left boundary along this dimension */
        for (i = 0; i < ghosts; i++) dxinv[i] = dxinv[ghosts];
      }
      if (simobj[n].mpi.ip[d] == simobj[n].mpi.iproc[d]-1) {
        /* fill right boundary along this dimension */
        for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) dxinv[i] = dxinv[dim[d]+ghosts-1];
      }
      offset  += (dim[d] + 2*ghosts);
    }

    /* calculate volume integral of the initial solution */
    IERR VolumeIntegral(  simobj[n].solver.VolumeIntegralInitial,
                          simobj[n].solver.u,
                          &(simobj[n].solver),
                          &(simobj[n].mpi) ); CHECKERR(ierr);
    if (!simobj[n].mpi.rank) {
      if (nsims > 1) printf("Volume integral of the initial solution on domain %d:\n", n);
      else           printf("Volume integral of the initial solution:\n");
      for (d=0; d<simobj[n].solver.nvars; d++) {
        printf("%2d:  %1.16E\n",d,simobj[n].solver.VolumeIntegralInitial[d]);
      }
    }
    /* Set initial total boundary flux integral to zero */
    _ArraySetValue_(simobj[n].solver.TotalBoundaryIntegral,simobj[n].solver.nvars,0);

  }

  return(0); 
}
