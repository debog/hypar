/*! @file VlasovWriteEField.c
    @author Debojyoti Ghosh
    @brief Write out the self-consistent electric field to file
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

/*! Write out the self-consistent field to file */
int VlasovWriteEField( void* s, /*!< Solver object of type #HyPar */
                       void* m  /*!< MPI object of type #MPIVariables */
                     )
{
  HyPar         *solver = (HyPar*)          s;
  MPIVariables  *mpi    = (MPIVariables*)   m;
  Vlasov        *param  = (Vlasov*)         solver->physics;

  if (param->self_consistent_electric_field) {

    char fname_root[_MAX_STRING_SIZE_] = "efield";
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
        cfg_dim = param->cfg_dim;

    int cfg_dim_global[cfg_dim];
    _ArrayCopy1D_(solver->dim_global, cfg_dim_global, cfg_dim);
    int cfg_dim_local[cfg_dim];
    _ArrayCopy1D_(solver->dim_local, cfg_dim_local, cfg_dim);

    /* gather the spatial coordinates into a global array */
    double *xg;
    {
      int size_global = 0;
      for (d=0; d<cfg_dim; d++) size_global += cfg_dim_global[d];
      xg = (double*) calloc (size_global,sizeof(double));
      _ArraySetValue_(xg,size_global,0.0);
  
      int offset_global, offset_local;
      offset_global = offset_local = 0;
      for (d=0; d<cfg_dim; d++) {
        IERR MPIGatherArray1D(  mpi,
                                (mpi->rank?NULL:&xg[offset_global]),
                                &solver->x[offset_local+ghosts],
                                mpi->is[d],
                                mpi->ie[d],
                                solver->dim_local[d],
                                0);
        offset_global += cfg_dim_global[d];
        offset_local  += cfg_dim_local [d] + 2*ghosts;
      }
    }

    /* gather the field into a global array */
    double *fieldg;
    {
      /* hard-coded for 1D-1V for now */
      int maxdim = 1;

      int size_global = 0;
      for (d=0; d<maxdim; d++) size_global += cfg_dim_global[d];
      fieldg = (double*) calloc (size_global,sizeof(double));
      _ArraySetValue_(fieldg,size_global,0.0);
  
      int offset_global, offset_local;
      offset_global = offset_local = 0;
      for (d=0; d<maxdim; d++) {
#ifdef fftw
        IERR MPIGatherArray1D(  mpi,
                                (mpi->rank?NULL:&fieldg[offset_global]),
                                &param->field[offset_local+ghosts],
                                mpi->is[d],
                                mpi->ie[d],
                                solver->dim_local[d],
                                0);
#endif
        offset_global += cfg_dim_global[d];
        offset_local  += cfg_dim_local [d] + 2*ghosts;
      }
    }

    if (!mpi->rank && solver->WriteOutput) {
      int index_cfg[cfg_dim];
      WriteText(  cfg_dim,
                  1,
                  cfg_dim_global,
                  xg,
                  fieldg,
                  filename,
                  index_cfg );
    }

    /* free up arrays */
    free(xg);
    free(fieldg);

  }

  return(0);
}
