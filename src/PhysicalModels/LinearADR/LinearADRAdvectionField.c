/*! @file LinearADRAdvectionField.c
    @author Debojyoti Ghosh
    @brief Contains the function to read in a non-constant advection field
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <io.h>
#include <physicalmodels/linearadr.h>
#include <hypar.h>
#include <mpivars.h>

/*! Set the advection field over the domain - reads the spatially-varying
 * advection data from a file, if available. The array to store the field
 * \b must already be allocated.
 *
 * For a single simulation, it reads in the data from a file. The data
 * must have the same grid dimensions as the solution.
 *
 * For an ensemble simulation, it reads in the advection field from the
 * file corresponding to the \b idx (index of this simulation). The
 * data in this file must have the same grid dimensions as the domain
 * it is being read in for.
 *
 * For a sparse grids simulation, it reads in the advection field from
 * the file with data that has the grid dimensions as the full grid. The
 * field on the current grid is obtained by interpolation.
*/
int LinearADRAdvectionField(void *s,      /*!< Solver object of type #HyPar */
                            void *m,      /*!< MPI object of type #MPIVariables*/
                            int  idx,     /*!< Index of this simulation */
                            int  nsims,   /*!< Total number of simulations */
                            int *dim_adv  /*!< Grid dimensions of the advection field stored in file */
                           )
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  LinearADR     *param  = (LinearADR*) solver->physics;
  double        *adv    = param->a;

  /* sanity checks */
  if (param->constant_advection == 1) {
    fprintf(stderr,"Error in LinearADRAdvectionField(): param->constant_advection == 1!\n");
    return(1);
  }
  if (!strcmp(param->adv_filename,"none")) {
    fprintf(stderr,"Error in LinearADRAdvectionField(): param->adv_filename is \"none\"!\n");
    return(1);
  }

  int d, done, flag;
  int *dim = solver->dim_local;
  int ghosts = solver->ghosts;

  /* number of advection field components at a grid point is
   * number of spatial dimensions times number of solution
   * components */
  int adv_nvar = solver->ndims * solver->nvars;
  if (param->adv_arr_size != adv_nvar*solver->npoints_local_wghosts) {
    fprintf(stderr,"Error in LinearADRAdvectionField(): Incorrect adv_arr_size\n");
    return(1);
  }

  char fname_root[_MAX_STRING_SIZE_];
  strcpy(fname_root, param->adv_filename);
  if (idx >= 0) {
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(idx, index, (int)log10(nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }
  }

  /* read spatially-varying advection field from provided file */
  if (dim_adv == NULL) {
    int ierr = ReadArray( solver->ndims,
                          adv_nvar,
                          solver->dim_global,
                          solver->dim_local,
                          solver->ghosts,
                          solver,
                          mpi,
                          NULL,
                          adv,
                          fname_root,
                          &flag);
    if (ierr) {
      fprintf(stderr,"Error in LinearADRAdvectionField()\n");
      fprintf(stderr,"  ReadArray() returned error!\n");
      return ierr;
    }
  } else {
    int ierr = ReadArraywInterp(solver->ndims,
                                adv_nvar,
                                solver->dim_global,
                                solver->dim_local,
                                dim_adv,
                                solver->ghosts,
                                solver,
                                mpi,
                                NULL,
                                adv,
                                fname_root,
                                &flag);
    if (ierr) {
      fprintf(stderr,"Error in LinearADRAdvectionField()\n");
      fprintf(stderr,"  ReadArraywInterp() returned error!\n");
      return ierr;
    }
  }

  if (!flag) {
    /* if advection file not available, set it to zero */
    _ArraySetValue_(adv,solver->npoints_local_wghosts*adv_nvar,0.0);
  }

  /* if parallel, exchange MPI-boundary ghost point data */
  IERR MPIExchangeBoundariesnD( solver->ndims,
                                adv_nvar,
                                solver->dim_local,
                                solver->ghosts,
                                mpi,
                                adv );

  /* Along each dimension:
     - If boundaries are periodic, then set ghost points accordingly if
       number of MPI ranks along that dimension is 1 (if more than 1,
       MPIExchangeBoundariesnD() called above has taken care of it).
     - Else, extrapolate at the boundaries. */
  int indexb[solver->ndims],
      indexi[solver->ndims],
      bounds[solver->ndims],
      offset[solver->ndims];
  for (d = 0; d < solver->ndims; d++) {
    if (solver->isPeriodic[d] && (mpi->iproc[d] == 1)) {
      _ArrayCopy1D_(dim,bounds,solver->ndims); bounds[d] = ghosts;
      /* left boundary */
      done = 0; _ArraySetValue_(indexb,solver->ndims,0);
      _ArraySetValue_(offset,solver->ndims,0); offset[d] = -ghosts;
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,solver->ndims); indexi[d] = indexb[d] + dim[d] - ghosts;
        int p1; _ArrayIndex1DWO_(solver->ndims,dim,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (solver->ndims,dim,indexi,ghosts,p2);
        _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
        _ArrayIncrementIndex_(solver->ndims,bounds,indexb,done);
      }
      /* right boundary */
      done = 0; _ArraySetValue_(indexb,solver->ndims,0);
      _ArraySetValue_(offset,solver->ndims,0); offset[d] = dim[d];
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,solver->ndims);
        int p1; _ArrayIndex1DWO_(solver->ndims,dim,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (solver->ndims,dim,indexi,ghosts,p2);
        _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
        _ArrayIncrementIndex_(solver->ndims,bounds,indexb,done);
      }
    } else {
      /* left boundary */
      if (!mpi->ip[d]) {
        _ArrayCopy1D_(dim,bounds,solver->ndims); bounds[d] = ghosts;
        _ArraySetValue_(offset,solver->ndims,0); offset[d] = -ghosts;
        done = 0; _ArraySetValue_(indexb,solver->ndims,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,solver->ndims); indexi[d] = ghosts-1-indexb[d];
          int p1; _ArrayIndex1DWO_(solver->ndims,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (solver->ndims,dim,indexi,ghosts,p2);
          _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
          _ArrayIncrementIndex_(solver->ndims,bounds,indexb,done);
        }
      }
      /* right boundary */
      if (mpi->ip[d] == mpi->iproc[d]-1) {
        _ArrayCopy1D_(dim,bounds,solver->ndims); bounds[d] = ghosts;
        _ArraySetValue_(offset,solver->ndims,0); offset[d] = dim[d];
        done = 0; _ArraySetValue_(indexb,solver->ndims,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,solver->ndims); indexi[d] = dim[d]-1-indexb[d];
          int p1; _ArrayIndex1DWO_(solver->ndims,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (solver->ndims,dim,indexi,ghosts,p2);
          _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
          _ArrayIncrementIndex_(solver->ndims,bounds,indexb,done);
        }
      }
    }
  }

  return(0);
}
