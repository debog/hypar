/*! @file TimePrintStep.c
    @brief Print to screen
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <simulation.h>
#include <timeintegration.h>

/*!
  Print information to screen (also calls any physics-specific 
  printing function, if defined).
*/
int TimePrintStep(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->simulation;
  int ns, nsims = TS->nsims;

  if ((!TS->rank) && ((TS->iter+1)%sim[0].solver.screen_op_iter == 0)) {
    printf("Iteration: %7d  "       ,TS->iter+1  );
    printf("Time: %1.3E  "          ,TS->waqt    );
    printf("Max CFL: %1.3E  "       ,TS->max_cfl );
    printf("Max Diff. No.: %1.3E  " ,TS->max_diff);
    printf("Norm: %1.4E  "          ,TS->norm    );

    /* calculate and print conservation error */
    if (!strcmp(sim[0].solver.ConservationCheck,"yes")) {

      double error = 0;
      for  (ns = 0; ns < nsims; ns++) {
        int v;
        for (v=0; v<sim[ns].solver.nvars; v++) {
          error += (sim[ns].solver.ConservationError[v] * sim[ns].solver.ConservationError[v]);
        }
      }
      error = sqrt(error);
      printf("Conservation loss: %1.4E",error);

    }
    printf("\n");

    /* print physics-specific info, if available */
    for (ns = 0; ns < nsims; ns++) {
      if (nsims > 1) printf("Domain %d:\n", ns);
      if (sim[ns].solver.PrintStep) {
        sim[ns].solver.PrintStep( &(sim[ns].solver),
                                  &(sim[ns].mpi),
                                  TS->waqt );
      }
    }
  }

  return(0);
}

