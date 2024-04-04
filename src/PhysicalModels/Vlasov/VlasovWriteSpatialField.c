/*! @file VlasovWriteSpatialField.c
    @author Ping-Hsuan Tsai
    @brief Write out the field to file
*/

#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <io.h>
#include <mpivars.h>
#include <hypar.h>
#include <physicalmodels/vlasov.h>

/*! Write out a spatial field or variable to file */
int VlasovWriteSpatialField( void   *s,         /*!< Solver object of type #HyPar */
                             void   *m,         /*!< MPI object of type #MPIVariables */
                             double *a_field,   /*!< Vector field to write */
                             char   *fname_root /*!< Filename root (extension is added automatically).
                                                      For unsteady output, a numerical index is added
                                                      that is the same as for the solution
                                                      output files. */
                            )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  Vlasov        *param  = (Vlasov*)       solver->physics;

  if (!solver->WriteOutput) return 0;

  if (solver->nsims > 1) {
    char index[_MAX_STRING_SIZE_];
    GetStringFromInteger(solver->my_idx, index, (int)log10(solver->nsims)+1);
    strcat(fname_root, "_");
    strcat(fname_root, index);
    strcat(fname_root, "_");
  }

  char filename[_MAX_STRING_SIZE_] = "";
  strcat(filename,fname_root);
  if (!strcmp(solver->op_overwrite,"no")) {
    strcat(filename,"_");
    strcat(filename,solver->filename_index);
  }
  strcat(filename,".dat");

  int d,
      ghosts = solver->ghosts,
      ndims_x = param->ndims_x;

  int dim_global_x[ndims_x];
  _ArrayCopy1D_(solver->dim_global, dim_global_x, ndims_x);
  int dim_local_x[ndims_x];
  _ArrayCopy1D_(solver->dim_local, dim_local_x, ndims_x);

  /* gather the spatial coordinates into a global array */
  double *xg;
  {
    int size_g = param->npts_global_x;
    xg = (double*) calloc (size_g, sizeof(double));
    _ArraySetValue_(xg, size_g, 0.0);

    int offset_global, offset_local;
    offset_global = offset_local = 0;
    for (d=0; d<ndims_x; d++) {
      IERR MPIGatherArray1D(  mpi,
                              (mpi->rank?NULL:&xg[offset_global]),
                              &solver->x[offset_local+ghosts],
                              mpi->is[d],
                              mpi->ie[d],
                              solver->dim_local[d],
                              0);
      offset_global += dim_global_x[d];
      offset_local  += dim_local_x [d] + 2*ghosts;
    }
  }

  /* gather the field into a global array */
  double *field_g;
  {
    int size_g = param->npts_global_x * param->ndims_x;
    field_g = (double*) calloc (size_g, sizeof(double));
    _ArraySetValue_(field_g, size_g, 0.0);

    if (param->ndims_x > 1) {
      if (!mpi->rank) {
        fprintf(stderr,"Warning in VlasovWriteSpatialField():\n");
        fprintf(stderr,"  field writing not yet supported for >1 spatial dimensions.\n");
      }
    } else {
      IERR MPIGatherArray1D(  mpi,
                              (mpi->rank ? NULL : field_g),
                              (a_field+ghosts),
                              mpi->is[0],
                              mpi->ie[0],
                              solver->dim_local[0],
                              0);
    }
  }

  if (!mpi->rank) {
    int idx[ndims_x];
    WriteText(  ndims_x,
                1,
                dim_global_x,
                xg,
                field_g,
                filename,
                idx );
  }

  /* free up arrays */
  free(xg);
  free(field_g);

  return 0;
}
