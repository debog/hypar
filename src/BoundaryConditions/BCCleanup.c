/*! @file BCCleanup.c
    @author Debojyoti Ghosh, Youngdae Kim
    @brief Function to clean up boundary-conditions related variables
*/
#include <stdlib.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#endif
#include <boundaryconditions.h>

/*! Cleans up a boundary object of type #DomainBoundary */
int BCCleanup(  void *b, /*!< Boundary object of type #DomainBoundary*/
                int flag_gpu /*!< Flag indicating if GPU is being used */ )
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  free(boundary->xmin);
  free(boundary->xmax);
  free(boundary->is);
  free(boundary->ie);
  if (boundary->DirichletValue) free(boundary->DirichletValue);
  if (boundary->SpongeValue   ) free(boundary->SpongeValue   );
  if (boundary->FlowVelocity  ) free(boundary->FlowVelocity  );

  if (boundary->UnsteadyDirichletSize) free(boundary->UnsteadyDirichletSize);
  if (boundary->UnsteadyDirichletData) free(boundary->UnsteadyDirichletData);

  if (boundary->UnsteadyTemperatureSize)  free(boundary->UnsteadyTemperatureSize);
  if (boundary->UnsteadyTimeLevels)       free(boundary->UnsteadyTimeLevels);
  if (boundary->UnsteadyTemperatureData)  free(boundary->UnsteadyTemperatureData);

#if defined(HAVE_CUDA)
  if (flag_gpu) {
    gpuFree(boundary->gpu_bounds);
    gpuFree(boundary->gpu_is);
    gpuFree(boundary->gpu_ie);
    if (boundary->FlowVelocity) gpuFree(boundary->gpu_FlowVelocity  );
  }
#endif

  return 0;
}
