/*! @file TimePrintStep.c
    @brief Print to screen
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <simulation_object.h>
#include <timeintegration.h>

/*!
  Print information to screen (also calls any physics-specific
  printing function, if defined).
*/
int TimePrintStep(void *a_ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration* TS = (TimeIntegration*) a_ts;
  SimulationObject* sim = (SimulationObject*) TS->m_simulation;
  int ns, nsims = TS->m_nsims;

  if ((!TS->m_rank) && ((TS->m_iter+1)%sim[0].solver.m_screen_op_iter == 0)) {
    if (nsims > 1) {
      printf("--\n");
      printf("iter=%7d,  t=%1.3e\n", TS->m_iter+1, TS->m_waqt);
      if (TS->m_max_cfl >= 0) printf("  CFL=%1.3E\n", TS->m_max_cfl);
      if (TS->m_norm >= 0) printf("  norm=%1.4E\n", TS->m_norm);
      printf("  wctime=%1.1E (s)\n",TS->m_iter_wctime);
    } else {
      printf("iter=%7d  ",TS->m_iter+1 );
      printf("t=%1.3E  ",TS->m_waqt );
      if (TS->m_max_cfl >= 0) printf("CFL=%1.3E  ",TS->m_max_cfl );
      if (TS->m_norm >= 0) printf("norm=%1.4E  ",TS->m_norm );
      printf("wctime: %1.1E (s)  ",TS->m_iter_wctime);
    }

    /* calculate and print conservation error */
    if (!strcmp(sim[0].solver.m_conservation_check,"yes")) {

      double error = 0;
      for  (ns = 0; ns < nsims; ns++) {
        int v;
        for (v=0; v<sim[ns].solver.m_nvars; v++) {
          error += (sim[ns].solver.m_conservation_error[v] * sim[ns].solver.m_conservation_error[v]);
        }
      }
      error = sqrt(error);
      printf("  cons_err=%1.4E\n", error);

    } else {

      if (nsims == 1) printf("\n");

    }

    /* print physics-specific info, if available */
    for (ns = 0; ns < nsims; ns++) {
      if (sim[ns].solver.PrintStep) {
        if (nsims > 1) printf("Physics-specific output for domain %d:\n", ns);
        sim[ns].solver.PrintStep( &(sim[ns].solver),
                                  &(sim[ns].mpi),
                                  TS->m_waqt );
      }
    }
    if (nsims > 1) {
      printf("--\n");
      printf("\n");
    }
  }

  return(0);
}

