/*! @file VlasovWriteEFieldAndPotential.c
    @author Debojyoti Ghosh
    @brief Write out the electric field and potential to file
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
int VlasovWriteSpatialField( void*, void*, double*, char* );

/*! Write out the electric field and potential to file */
int VlasovWriteEFieldAndPotential(  void* s, /*!< Solver object of type #HyPar */
                                    void* m  /*!< MPI object of type #MPIVariables */ )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  Vlasov        *param  = (Vlasov*)       solver->physics;

  {
    char fname_root[_MAX_STRING_SIZE_] = "efield";
    VlasovWriteSpatialField( solver, mpi, param->e_field, fname_root );
  }

  if (param->self_consistent_electric_field) {
    char fname_root[_MAX_STRING_SIZE_] = "potential";
    VlasovWriteSpatialField( solver, mpi, param->potential, fname_root );
  }

  return 0;
}
