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
int WENOCleanup(void *a_s, /*!< WENO object of type #WENOParameters */
                int flag_gpu /*!< flag to indicate if on GPU */ )
{
  WENOParameters  *weno   = (WENOParameters*) a_s;

  if (weno->m_offset) free(weno->m_offset);
#if defined(HAVE_CUDA)
  if (a_flag_gpu) {
    if (weno->m_w1) gpuFree(weno->m_w1);
    if (weno->m_w2) gpuFree(weno->m_w2);
    if (weno->m_w3) gpuFree(weno->m_w3);
  } else {
#endif
    if (weno->m_w1) free(weno->m_w1);
    if (weno->m_w2) free(weno->m_w2);
    if (weno->m_w3) free(weno->m_w3);
#if defined(HAVE_CUDA)
  }
#endif

  return(0);
}
