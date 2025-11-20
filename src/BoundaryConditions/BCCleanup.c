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
  free(boundary->m_xmin);
  free(boundary->m_xmax);
  free(boundary->m_is);
  free(boundary->m_ie);
  if (boundary->m_DirichletValue) free(boundary->m_DirichletValue);
  if (boundary->m_SpongeValue   ) free(boundary->m_SpongeValue   );
  if (boundary->m_FlowVelocity  ) free(boundary->m_FlowVelocity  );

  if (boundary->m_UnsteadyDirichletSize) free(boundary->m_UnsteadyDirichletSize);
  if (boundary->m_UnsteadyDirichletData) free(boundary->m_UnsteadyDirichletData);

  if (boundary->m_UnsteadyTemperatureSize)  free(boundary->m_UnsteadyTemperatureSize);
  if (boundary->m_UnsteadyTimeLevels)       free(boundary->m_UnsteadyTimeLevels);
  if (boundary->m_UnsteadyTemperatureData)  free(boundary->m_UnsteadyTemperatureData);

#if defined(HAVE_CUDA)
  if (flag_gpu) {
    gpuFree(boundary->m_gpu_bounds);
    gpuFree(boundary->m_gpu_is);
    gpuFree(boundary->m_gpu_ie);
    if (boundary->m_FlowVelocity) gpuFree(boundary->m_gpu_FlowVelocity  );
  }
#endif

  return 0;
}
