#ifdef with_librom

/*! @file CalculateROMDiff.c
    @author Debojyoti Ghosh
    @brief Computes the diff between PDE and ROM solutions.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <timeintegration.h>
#include <mpivars.h>
#include <hypar.h>

/*! Calculates the L1, L2, & Linf norms of the diff between the PDE
 * solution and the predicted solution by libROM.
 * error in the solution if the exact solution is
*/
int CalculateROMDiff( void *s, /*!< Solver object of type #HyPar */
                      void *m  /*!< MPI object of type #MPIVariables */ )
{
  HyPar* solver = (HyPar*) s;
  MPIVariables* mpi = (MPIVariables*) m;
  double sum = 0, global_sum = 0;

  static const double tolerance = 1e-15;

  int size = solver->m_npoints_local_wghosts * solver->m_nvars;
  double* u_diff = (double*) calloc (size, sizeof(double));

  /* calculate solution norms (for relative error) */
  double solution_norm[3] = {0.0,0.0,0.0};
  /* L1 */
  sum = ArraySumAbsnD ( solver->m_nvars,
                        solver->m_ndims,
                        solver->m_dim_local,
                        solver->m_ghosts,
                        solver->m_index,
                        solver->m_u );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
  solution_norm[0] = global_sum/((double)solver->m_npoints_global);
  /* L2 */
  sum = ArraySumSquarenD  ( solver->m_nvars,
                            solver->m_ndims,
                            solver->m_dim_local,
                            solver->m_ghosts,
                            solver->m_index,
                            solver->m_u );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
  solution_norm[1] = sqrt(global_sum/((double)solver->m_npoints_global));
  /* Linf */
  sum = ArrayMaxnD  ( solver->m_nvars,
                      solver->m_ndims,
                      solver->m_dim_local,
                      solver->m_ghosts,
                      solver->m_index
                      ,solver->m_u  );
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->m_world);
  solution_norm[2] = global_sum;

  /* set u_diff to PDE solution */
  _ArrayCopy1D_(solver->m_u, u_diff, size);
  /* subtract the ROM solutions */
  _ArrayAXPY_(solver->m_u_rom_predicted, -1.0, u_diff, size);

  /* calculate L1 norm of error */
  sum = ArraySumAbsnD ( solver->m_nvars,
                        solver->m_ndims,
                        solver->m_dim_local,
                        solver->m_ghosts,
                        solver->m_index,
                        u_diff  );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
  solver->m_rom_diff_norms[0] = global_sum/((double)solver->m_npoints_global);

  /* calculate L2 norm of error */
  sum = ArraySumSquarenD  ( solver->m_nvars,
                            solver->m_ndims,
                            solver->m_dim_local,
                            solver->m_ghosts,
                            solver->m_index,
                            u_diff  );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
  solver->m_rom_diff_norms[1] = sqrt(global_sum/((double)solver->m_npoints_global));

  /* calculate Linf norm of error */
  sum = ArrayMaxnD  ( solver->m_nvars,
                      solver->m_ndims,
                      solver->m_dim_local,
                      solver->m_ghosts,
                      solver->m_index,
                      u_diff  );
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->m_world);
  solver->m_rom_diff_norms[2] = global_sum;

  /* decide whether to normalize and report relative diff norms,
    or report absolute diff norms. */
  if (    (solution_norm[0] > tolerance)
      &&  (solution_norm[1] > tolerance)
      &&  (solution_norm[2] > tolerance) ) {
    solver->m_rom_diff_norms[0] /= solution_norm[0];
    solver->m_rom_diff_norms[1] /= solution_norm[1];
    solver->m_rom_diff_norms[2] /= solution_norm[2];
  }

  free(u_diff);
  return 0;
}

#endif
