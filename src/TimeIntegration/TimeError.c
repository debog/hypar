/*! @file TimeError.c
    @brief Compute time integration error
    @author Debojyoti Ghosh
*/

#include <math.h>
#include <string.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Computes the time integration error in the solution: If the time
  integration method chosen has a mechanism to compute the error
  in the numerical integration, it is computed in this function. In
  addition, if an exact solution is available (see function argument
  \a uex, the error of the computed solution (stored in #HyPar::m_u)
  with respect to this exact solution is also computed. These errors
  are written to a text file whose name depends on the time integration
  method being used.

  Time integration method with error estimation currently implemented:
  + #_GLM_GEE_ (General Linear Methods with Global Error Estimation)
    (see TimeGLMGEE(), TimeGLMGEEInitialize() )
*/
int TimeError(
                void    *a_s,   /*!< Solver object of type #HyPar */
                void    *a_m,   /*!< MPI object of type #MPIVariables */
                double  *a_uex  /*!< Exact solution (stored with the same
                                   array layout as #HyPar::m_u (may be
                                   NULL) */
             )
{
  HyPar               *solver = (HyPar*)           a_s;
  MPIVariables        *mpi    = (MPIVariables*)    a_m;
  TimeIntegration     *TS     = (TimeIntegration*) solver->m_time_integrator;
  int                 size    = solver->m_npoints_local_wghosts * solver->m_nvars;
  double              sum     = 0.0, global_sum = 0.0;
  static const double tolerance = 1e-15;
  if (!TS) return(0);

  if (!strcmp(solver->m_time_scheme,_GLM_GEE_)) {
    /* For GLM-GEE methods, calculate the norm of the estimated global error */
    GLMGEEParameters *params = (GLMGEEParameters*) solver->m_msti;
    double error[6] = {0,0,0,0,0,0}, *Uerr;
    if (!strcmp(params->ee_mode,_GLM_GEE_YEPS_)) Uerr = TS->m_U[params->r];
    else {
      Uerr = TS->m_U[0];
      _ArraySubtract1D_(Uerr,solver->m_u,TS->m_U[params->r],size);
      _ArrayScale1D_(Uerr,(1.0/(1.0-params->gamma)),size);
    }

    /* calculate solution norm for relative errors */
    double sol_norm[3] = {0.0,0.0,0.0};
    /* L1 */
    sum = ArraySumAbsnD   (solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                           solver->m_ghosts,solver->m_index,solver->m_u);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
    sol_norm[0] = global_sum/((double)solver->m_npoints_global);
    /* L2 */
    sum = ArraySumSquarenD(solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                           solver->m_ghosts,solver->m_index,solver->m_u);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
    sol_norm[1] = sqrt(global_sum/((double)solver->m_npoints_global));
    /* Linf */
    sum = ArrayMaxnD      (solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                           solver->m_ghosts,solver->m_index,solver->m_u);
    global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->m_world);
    sol_norm[2] = global_sum;

    /* calculate L1 norm of error */
    sum = ArraySumAbsnD   (solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                           solver->m_ghosts,solver->m_index,Uerr);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
    error[0] = global_sum/((double)solver->m_npoints_global);
    /* calculate L2 norm of error */
    sum = ArraySumSquarenD(solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                           solver->m_ghosts,solver->m_index,Uerr);
    global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
    error[1] = sqrt(global_sum/((double)solver->m_npoints_global));
    /* calculate Linf norm of error */
    sum = ArrayMaxnD      (solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                           solver->m_ghosts,solver->m_index,Uerr);
    global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->m_world);
    error[2] = global_sum;

    if (a_uex) {
      _ArrayAXBY_(TS->m_Udot[0],1.0,solver->m_u,-1.0,a_uex,size);
      _ArrayAXPY_(Uerr,-1.0,TS->m_Udot[0],size);
      /* calculate L1 norm of error */
      sum = ArraySumAbsnD   (solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                             solver->m_ghosts,solver->m_index,TS->m_Udot[0]);
      global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
      error[3] = global_sum/((double)solver->m_npoints_global);
      /* calculate L2 norm of error */
      sum = ArraySumSquarenD(solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                             solver->m_ghosts,solver->m_index,TS->m_Udot[0]);
      global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->m_world);
      error[4] = sqrt(global_sum/((double)solver->m_npoints_global));
      /* calculate Linf norm of error */
      sum = ArrayMaxnD      (solver->m_nvars,solver->m_ndims,solver->m_dim_local,
                             solver->m_ghosts,solver->m_index,TS->m_Udot[0]);
      global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->m_world);
      error[5] = global_sum;
    } else error[3] = error[4] = error[5] = -1;

    if (   (sol_norm[0] > tolerance)
        && (sol_norm[1] > tolerance)
        && (sol_norm[2] > tolerance) ) {
      error[0] /= sol_norm[0];
      error[1] /= sol_norm[1];
      error[2] /= sol_norm[2];
      if (a_uex) {
        error[3] /= sol_norm[0];
        error[4] /= sol_norm[1];
        error[5] /= sol_norm[2];
      }
    }

    /* write to file */
    if (!mpi->m_rank) {
      FILE *out;
      out = fopen("glm_err.dat","w");
      fprintf(out,"%1.16E  %1.16E  %1.16E  %1.16E  ",TS->m_dt,error[0],error[1],error[2]);
      fprintf(out,"%1.16E  %1.16E  %1.16E\n",error[3],error[4],error[5]);
      fclose(out);
      printf("Estimated time integration errors (GLM-GEE time-integration):-\n");
      printf("  L1         Error           : %1.16E\n",error[0]);
      printf("  L2         Error           : %1.16E\n",error[1]);
      printf("  Linfinity  Error           : %1.16E\n",error[2]);
    }
  }

  return(0);
}
