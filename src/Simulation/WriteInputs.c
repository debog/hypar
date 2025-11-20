/*! @file ReadInputs.c
    @author Debojyoti Ghosh
    @brief Read the input parameters from \b solver.inp
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <timeintegration.h>
#include <mpivars.h>
#include <simulation_object.h>

/*! Write the simulation inputs read from the file \b solver.inp. */
int WriteInputs ( void  *a_s,     /*!< Array of simulation objects of type #SimulationObject
                                     of size nsims */
                  int   a_nsims,  /*!< Number of simulation objects */
                  int   a_rank    /*!< MPI rank of this process */
                )
{
  SimulationObject *sim = (SimulationObject*) a_s;
  int n;

  if (sim == NULL)  return 0;

  if (!a_rank) {

    printf("  No. of dimensions                          : %d\n",sim[0].solver.m_ndims);
    printf("  No. of variables                           : %d\n",sim[0].solver.m_nvars);
    if (a_nsims > 1) {
      printf("  Domain sizes:\n");
      for (int n = 0; n < a_nsims; n++) {
        printf("    domain %3d - ", n);
        for (int i=0; i<sim[n].solver.m_ndims; i++) printf ("%d ",sim[n].solver.m_dim_global[i]);
        printf("\n");
      }
#ifndef serial
      printf("  Processes along each dimension:\n");
      for (int n = 0; n < a_nsims; n++) {
        printf("    domain %3d - ", n);
        for (int i=0; i<sim[n].solver.m_ndims; i++) printf ("%d ",sim[n].mpi.m_iproc[i]);
        printf("\n");
      }
#endif
      printf("  Exact solution domain sizes:\n");
      for (int n = 0; n < a_nsims; n++) {
        printf("    domain %3d - ", n);
        for (int i=0; i<sim[n].solver.m_ndims; i++) printf ("%d ",sim[n].solver.m_dim_global_ex[i]);
        printf("\n");
      }
    } else {
      printf("  Domain size                                : ");
      for (int i=0; i<sim[0].solver.m_ndims; i++) printf ("%d ",sim[0].solver.m_dim_global[i]);
      printf("\n");
#ifndef serial
      printf("  Processes along each dimension             : ");
      for (int i=0; i<sim[0].solver.m_ndims; i++) printf ("%d ",sim[0].mpi.m_iproc[i]);
      printf("\n");
#endif
      printf("  Exact solution domain size                 : ");
      for (int i=0; i<sim[0].solver.m_ndims; i++) printf ("%d ",sim[0].solver.m_dim_global_ex[i]);
      printf("\n");
    }
    printf("  No. of ghosts pts                          : %d\n"     ,sim[0].solver.m_ghosts              );
    printf("  No. of iter.                               : %d\n"     ,sim[0].solver.m_n_iter              );
    printf("  Restart iteration                          : %d\n"     ,sim[0].solver.m_restart_iter        );
#ifdef with_petsc
    if (sim[0].solver.m_use_petsc_ts)
      printf("  Time integration scheme                    : PETSc \n"                            );
    else {
      printf("  Time integration scheme                    : %a_s ",sim[0].solver.m_time_scheme             );
      if (strcmp(sim[0].solver.m_time_scheme,_FORWARD_EULER_)) {
        printf("(%a_s)",sim[0].solver.m_time_scheme_type                                                    );
      }
      printf("\n");
    }
#else
    printf("  Time integration scheme                    : %a_s ",sim[0].solver.m_time_scheme               );
    if (strcmp(sim[0].solver.m_time_scheme,_FORWARD_EULER_)) {
      printf("(%a_s)",sim[0].solver.m_time_scheme_type                                                      );
    }
    printf("\n");
#endif
    printf("  Spatial discretization scheme (hyperbolic) : %a_s\n"     ,sim[0].solver.m_spatial_scheme_hyp  );
    printf("  Split hyperbolic flux term?                : %a_s\n"     ,sim[0].solver.m_split_hyperbolic_flux );
    printf("  Interpolation type for hyperbolic term     : %a_s\n"     ,sim[0].solver.m_interp_type         );
    printf("  Spatial discretization type   (parabolic ) : %a_s\n"     ,sim[0].solver.m_spatial_type_par    );
    printf("  Spatial discretization scheme (parabolic ) : %a_s\n"     ,sim[0].solver.m_spatial_scheme_par  );
    printf("  Time Step                                  : %E\n"     ,sim[0].solver.m_dt                  );
    printf("  Check for conservation                     : %a_s\n"     ,sim[0].solver.m_conservation_check   );
    printf("  Screen output iterations                   : %d\n"     ,sim[0].solver.m_screen_op_iter      );
    printf("  File output iterations                     : %d\n"     ,sim[0].solver.m_file_op_iter        );
    printf("  Initial solution file type                 : %a_s\n"     ,sim[0].solver.m_ip_file_type        );
    printf("  Initial solution read mode                 : %a_s"       ,sim[0].solver.m_input_mode          );
    if (strcmp(sim[0].solver.m_input_mode,"serial"))    printf("  [%d file IO a_rank(a_s)]\n",sim[0].mpi.m_N_IORanks  );
    else                                        printf("\n");
    printf("  Solution file write mode                   : %a_s"       ,sim[0].solver.m_output_mode         );
    if (strcmp(sim[0].solver.m_output_mode,"serial"))   printf("  [%d file IO a_rank(a_s)]\n",sim[0].mpi.m_N_IORanks  );
    else                                        printf("\n");
    printf("  Solution file format                       : %a_s\n"     ,sim[0].solver.m_op_file_format      );
    printf("  Overwrite solution file                    : %a_s\n"     ,sim[0].solver.m_op_overwrite        );
#if defined(HAVE_CUDA)
    printf("  Use GPU                                    : %a_s\n"     ,(sim[0].solver.m_use_gpu == 1)? "yes" : "no");
    printf("  GPU device no                              : %d\n"     ,(sim[0].solver.m_gpu_device_no));
#endif
    printf("  Physical model                             : %a_s\n"     ,sim[0].solver.m_model               );
    if (sim[0].solver.m_flag_ib) {
      printf("  Immersed Body                              : %a_s\n"     ,sim[0].solver.m_ib_filename         );
    }
  }

  return 0;
}
