/*! @file mpivars_cpp.h
    @brief C++ declarations for MPI-related functions
    @author Debojyoti Ghosh
 */

#ifndef _MPIVARS_CPP_H_
#define _MPIVARS_CPP_H_

#include <mpivars_struct.h>

/*! Broadcast a double to all ranks */
extern "C" int MPIBroadcast_double     (double*,int,int,void*);
/*! Broadcast an integer to all ranks */
extern "C" int MPIBroadcast_integer    (int*,int,int,void*);
/*! Broadcast a character to all ranks */
extern "C" int MPIBroadcast_character  (char*,int,int,void*);

/*! Create communicators required for the tridiagonal solver in compact schemes */
extern "C" int MPICreateCommunicators  (int,void*);
/*! Free/destroy communicators created */
extern "C" int MPIFreeCommunicators    (int,void*);

/*! Create I/O groups for file reading and writing -- Group leaders gather data from
 * all other ranks in the group and write to file, and read from file and sends the
 * data to the appropriate ranks in the group. Thus, number of ranks participating
 * in file I/O is equal to the number of groups (which is an input), and can be set
 * to the number of I/O nodes available. */
extern "C" int MPICreateIOGroups       (void*);

/*! Exchange boundary (ghost point) values for an essentially 1D array (like grid
 * coordinates) */
extern "C" int MPIExchangeBoundaries1D (void*,double*,int,int,int,int);

/*! Exchange boundary (ghost point) values for an n-dimensional array (like the
 * solution array) */
extern "C" int MPIExchangeBoundariesnD (int,int,int*,int,void*,double*);

/*! Gather local arrays into a global array for an essentially 1D array */
extern "C" int MPIGatherArray1D           (void*,double*,double*,int,int,int,int);
/*! Gather local arrays into a global array for an n-dimensional array */
extern "C" int MPIGatherArraynD           (int,void*,double*,double*,int*,int*,int,int);
/*! Gather local arrays into a global array for an n-dimensional array (with ghosts) */
extern "C" int MPIGatherArraynDwGhosts    (int,void*,double*,double*,int*,int*,int,int);
/*! Partition a global array into local arrays for an n-dimensional array */
extern "C" int MPIPartitionArraynD        (int,void*,double*,double*,int*,int*,int,int);
/*! Partition a global array into local arrays for an n-dimensional array */
extern "C" int MPIPartitionArraynDwGhosts (int,void*,double*,double*,int*,int*,int,int);
/*! Partition a global array into local arrays for an essentially 1D array */
extern "C" int MPIPartitionArray1D        (void*,double*,double*,int,int,int,int);

/*! fetch data from an n-dimensional local array on another rank */
extern "C" int MPIGetArrayDatanD       (double*,double*,int*,int*,int*,int*,int,int,int,void*);

/*! Calculate the local domain limits/extend in terms of the global domain */
extern "C" int MPILocalDomainLimits    (int,int,void*,int*,int*,int*);

/*! Find the maximum in an integer array over all ranks */
extern "C" int MPIMax_integer          (int*,int*,int,void*);
/*! Find the maximum in a long integer array over all ranks */
extern "C" int MPIMax_long             (long*,long*,int,void*);
/*! Find the maximum in a double array over all ranks */
extern "C" int MPIMax_double           (double*,double*,int,void*);
/*! Find the minimum in an integer array over all ranks */
extern "C" int MPIMin_integer          (int*,int*,int,void*);
/*! Find the minimum in a double array over all ranks */
extern "C" int MPIMin_double           (double*,double*,int,void*);
/*! Calculate the sum of an array of doubles over all ranks */
extern "C" int MPISum_double           (double*,double*,int,void*);
/*! Calculate the sum of an array of integers over all ranks */
extern "C" int MPISum_integer          (int*,int*,int,void*);

/*! Partition (along a dimension) the domain given global size and number of ranks */
extern "C" int MPIPartition1D          (int,int,int);

/*! Calculate 1D rank from the n-dimensional rank */
extern "C" int MPIRank1D               (int,int*,int*);
/*! Calculate the n-dimensional rank from the 1D rank */
extern "C" int MPIRanknD               (int,int,int*,int*);

/*! Generate a unique filename given the rank of the process to let that process
 * write to its own file */
void MPIGetFilename         (char*,void*,char*);

#endif
