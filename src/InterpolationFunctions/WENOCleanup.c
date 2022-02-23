/*! @file WENOCleanup.c
    @brief Cleans up allocations specific to WENO-type methods
    @author Debojyoti Ghosh, Youngdae Kim
*/

#include <stdlib.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#endif
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Cleans up all allocations related to the WENO-type methods.
*/
int WENOCleanup(void *s, /*!< WENO object of type #WENOParameters */
                int flag_gpu /*!< flag to indicate if on GPU */ )
{
  WENOParameters  *weno   = (WENOParameters*) s;

  if (weno->offset) free(weno->offset);
#if defined(HAVE_CUDA)
  if (flag_gpu) {
    if (weno->w1) gpuFree(weno->w1);
    if (weno->w2) gpuFree(weno->w2);
    if (weno->w3) gpuFree(weno->w3);
  } else {
#endif
    if (weno->w1) free(weno->w1);
    if (weno->w2) free(weno->w2);
    if (weno->w3) free(weno->w3);
#if defined(HAVE_CUDA)
  }
#endif

  return(0);
}
