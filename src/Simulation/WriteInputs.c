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
int WriteInputs ( void  *s,     /*!< Array of simulation objects of type #SimulationObject
                                     of size nsims */
                  int   nsims,  /*!< Number of simulation objects */
                  int   rank    /*!< MPI rank of this process */
                )
{
  SimulationObject *sim = (SimulationObject*) s;
  int n;

  if (sim == NULL)  return 0;

  if (!rank) {

    printf("  No. of dimensions                          : %d\n",sim[0].solver.ndims);
    printf("  No. of variables                           : %d\n",sim[0].solver.nvars);
    if (nsims > 1) {
      printf("  Domain sizes:\n");
      for (int n = 0; n < nsims; n++) {
        printf("    domain %3d - ", n);
        for (int i=0; i<sim[n].solver.ndims; i++) printf ("%d ",sim[n].solver.dim_global[i]);
        printf("\n");
      }
#ifndef serial
      printf("  Processes along each dimension:\n");
      for (int n = 0; n < nsims; n++) {
        printf("    domain %3d - ", n);
        for (int i=0; i<sim[n].solver.ndims; i++) printf ("%d ",sim[n].mpi.iproc[i]);
        printf("\n");
      }
#endif
      printf("  Exact solution domain sizes:\n");
      for (int n = 0; n < nsims; n++) {
        printf("    domain %3d - ", n);
        for (int i=0; i<sim[n].solver.ndims; i++) printf ("%d ",sim[n].solver.dim_global_ex[i]);
        printf("\n");
      }
    } else {
      printf("  Domain size                                : ");
      for (int i=0; i<sim[0].solver.ndims; i++) printf ("%d ",sim[0].solver.dim_global[i]);
      printf("\n");
#ifndef serial
      printf("  Processes along each dimension             : ");
      for (int i=0; i<sim[0].solver.ndims; i++) printf ("%d ",sim[0].mpi.iproc[i]);
      printf("\n");
#endif
      printf("  Exact solution domain size                 : ");
      for (int i=0; i<sim[0].solver.ndims; i++) printf ("%d ",sim[0].solver.dim_global_ex[i]);
      printf("\n");
    }
    printf("  No. of ghosts pts                          : %d\n"     ,sim[0].solver.ghosts              );
    printf("  No. of iter.                               : %d\n"     ,sim[0].solver.n_iter              );
    printf("  Restart iteration                          : %d\n"     ,sim[0].solver.restart_iter        );
#ifdef with_petsc
    if (sim[0].solver.use_petscTS)
      printf("  Time integration scheme                    : PETSc \n"                            );
    else {
      printf("  Time integration scheme                    : %s ",sim[0].solver.time_scheme             );
      if (strcmp(sim[0].solver.time_scheme,_FORWARD_EULER_)) {
        printf("(%s)",sim[0].solver.time_scheme_type                                                    );
      }
      printf("\n");
    }
#else
    printf("  Time integration scheme                    : %s ",sim[0].solver.time_scheme               );
    if (strcmp(sim[0].solver.time_scheme,_FORWARD_EULER_)) {
      printf("(%s)",sim[0].solver.time_scheme_type                                                      );
    }
    printf("\n");
#endif
    printf("  Spatial discretization scheme (hyperbolic) : %s\n"     ,sim[0].solver.spatial_scheme_hyp  );
    printf("  Split hyperbolic flux term?                : %s\n"     ,sim[0].solver.SplitHyperbolicFlux );
    printf("  Interpolation type for hyperbolic term     : %s\n"     ,sim[0].solver.interp_type         );
    printf("  Spatial discretization type   (parabolic ) : %s\n"     ,sim[0].solver.spatial_type_par    );
    printf("  Spatial discretization scheme (parabolic ) : %s\n"     ,sim[0].solver.spatial_scheme_par  );
    printf("  Time Step                                  : %E\n"     ,sim[0].solver.dt                  );
    printf("  Check for conservation                     : %s\n"     ,sim[0].solver.ConservationCheck   );
    printf("  Screen output iterations                   : %d\n"     ,sim[0].solver.screen_op_iter      );
    printf("  File output iterations                     : %d\n"     ,sim[0].solver.file_op_iter        );
    printf("  Initial solution file type                 : %s\n"     ,sim[0].solver.ip_file_type        );
    printf("  Initial solution read mode                 : %s"       ,sim[0].solver.input_mode          );
    if (strcmp(sim[0].solver.input_mode,"serial"))    printf("  [%d file IO rank(s)]\n",sim[0].mpi.N_IORanks  );
    else                                        printf("\n");
    printf("  Solution file write mode                   : %s"       ,sim[0].solver.output_mode         );
    if (strcmp(sim[0].solver.output_mode,"serial"))   printf("  [%d file IO rank(s)]\n",sim[0].mpi.N_IORanks  );
    else                                        printf("\n");
    printf("  Solution file format                       : %s\n"     ,sim[0].solver.op_file_format      );
    printf("  Overwrite solution file                    : %s\n"     ,sim[0].solver.op_overwrite        );
#if defined(HAVE_CUDA)
    printf("  Use GPU                                    : %s\n"     ,(sim[0].solver.use_gpu == 1)? "yes" : "no");
    printf("  GPU device no                              : %d\n"     ,(sim[0].solver.gpu_device_no));
#endif
    printf("  Physical model                             : %s\n"     ,sim[0].solver.model               );
    if (sim[0].solver.flag_ib) {
      printf("  Immersed Body                              : %s\n"     ,sim[0].solver.ib_filename         );
    }
  }

  return 0;
}
