/*! @file mpivars.h
    @brief MPI related function definitions.
    @author Debojyoti Ghosh
 */

#ifndef _MPIVARS_h_
#define _MPIVARS_h_

#include <mpivars_struct.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! Broadcast a double to all ranks */
int MPIBroadcast_double     (double*,int,int,void*);
/*! Broadcast an integer to all ranks */
int MPIBroadcast_integer    (int*,int,int,void*);
/*! Broadcast a character to all ranks */
int MPIBroadcast_character  (char*,int,int,void*);

/*! Create communicators required for the tridiagonal solver in compact schemes */
int MPICreateCommunicators  (int,void*);
/*! Free/destroy communicators created */
int MPIFreeCommunicators    (int,void*);

/*! Create I/O groups for file reading and writing -- Group leaders gather data from
 * all other ranks in the group and write to file, and read from file and sends the
 * data to the appropriate ranks in the group. Thus, number of ranks participating
 * in file I/O is equal to the number of groups (which is an input), and can be set
 * to the number of I/O nodes available. */
int MPICreateIOGroups       (void*);

/*! Exchange boundary (ghost point) values for an essentially 1D array (like grid
 * coordinates) */
int MPIExchangeBoundaries1D (void*,double*,int,int,int,int);

/*! Exchange boundary (ghost point) values for an n-dimensional array (like the
 * solution array) */
int MPIExchangeBoundariesnD (int,int,int*,int,void*,double*);

/*! Gather local arrays into a global array for an essentially 1D array */
int MPIGatherArray1D            (void*,double*,double*,int,int,int,int);
/*! Gather local arrays into a global array for an n-dimensional array */
int MPIGatherArraynD            (int,void*,double*,double*,int*,int*,int,int);
/*! Gather local arrays into a global array for an n-dimensional array (with ghosts) */
int MPIGatherArraynDwGhosts     (int,void*,double*,double*,int*,int*,int,int);
/*! Partition a global array into local arrays for an n-dimensional array */
int MPIPartitionArraynD         (int,void*,double*,double*,int*,int*,int,int);
/*! Partition a global array into local arrays for an n-dimensional array */
int MPIPartitionArraynDwGhosts  (int,void*,double*,double*,int*,int*,int,int);
/*! Partition a global array into local arrays for an essentially 1D array */
int MPIPartitionArray1D         (void*,double*,double*,int,int,int,int);

/*! fetch data from an n-dimensional local array on another rank */
int MPIGetArrayDatanD       (double*,double*,int*,int*,int*,int*,int,int,int,void*);

/*! Calculate the local domain limits/extend in terms of the global domain */
int MPILocalDomainLimits    (int,int,void*,int*,int*,int*);

/*! Find the maximum in an integer array over all ranks */
int MPIMax_integer          (int*,int*,int,void*);
/*! Find the maximum in a long integer array over all ranks */
int MPIMax_long             (long*,long*,int,void*);
/*! Find the maximum in a double array over all ranks */
int MPIMax_double           (double*,double*,int,void*);
/*! Find the minimum in an integer array over all ranks */
int MPIMin_integer          (int*,int*,int,void*);
/*! Find the minimum in a double array over all ranks */
int MPIMin_double           (double*,double*,int,void*);
/*! Calculate the sum of an array of doubles over all ranks */
int MPISum_double           (double*,double*,int,void*);
/*! Calculate the sum of an array of integers over all ranks */
int MPISum_integer          (int*,int*,int,void*);

/*! Partition (along a dimension) the domain given global size and number of ranks */
int MPIPartition1D          (int,int,int);

/*! Calculate 1D rank from the n-dimensional rank */
int MPIRank1D               (int,int*,int*);
/*! Calculate the n-dimensional rank from the 1D rank */
int MPIRanknD               (int,int,int*,int*);

/*! Generate a unique filename given the rank of the process to let that process
 * write to its own file */
void MPIGetFilename         (char*,void*,char*);

#if defined(HAVE_CUDA)
int gpuMPIExchangeBoundariesnD (int,int,const int*,int,void*,double*);
#endif

#ifdef __cplusplus
}
#endif

#endif
