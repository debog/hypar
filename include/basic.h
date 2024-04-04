/*! @file basic.h
    @brief Some basic definitions and macros
    @author Debojyoti Ghosh
*/

#ifndef _BASIC_H_
#define _BASIC_H_

#if defined(HAVE_CONFIG_H)
# include <config.h>
#endif

/*! Maximum size of character strings */
#define _MAX_STRING_SIZE_ 500
#ifdef debug
  #define IERR ierr=
  #define _DECLARE_IERR_ int ierr = 0;
  #define CHECKERR(ierr) { if (ierr) return(ierr); }
#else
  #define IERR
  #define _DECLARE_IERR_
  #define CHECKERR(ierr)
#endif

/*! Numbers smaller than this will be treated as zero */
#define _MACHINE_ZERO_ 1.0e-14

/*! \def _GetCoordinate_
 * Macro to get the spatial coordinate \a coord along the \a dir dimension of a grid point whose index in that dimension is \a i, in a grid of size \a dim with \a ghosts ghost points in each dimension. The array containing the spatial coordinates of the entire grid is \a x.
*/
#define _GetCoordinate_(dir,i,dim,ghosts,x,coord) \
  { \
    int arraycounter,offset = 0; \
    for (arraycounter = 0; arraycounter < dir; arraycounter++) offset += (dim[arraycounter]+2*ghosts); \
    coord = (x[offset+ghosts+i]); \
  }

/*! Minimum number of grid points per dimension per MPI rank */
#define _MIN_GRID_PTS_PER_PROC_ 4

#endif
