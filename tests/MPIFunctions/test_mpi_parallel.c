/*! @file test_mpi_parallel.c
 *  @brief Unit tests for MPI utility functions with actual MPI (2 ranks)
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef serial
#include <mpi.h>
#endif

#include <mpivars.h>
#include <common.h>

#define TEST_PASS 0
#define TEST_FAIL 1
#define TOLERANCE 1e-10

/* Simple test harness */
typedef struct {
  int passed;
  int failed;
  char last_test[256];
} TestStats;

void print_test_result(TestStats *stats, const char *test_name, int result, int rank) {
  if (rank == 0) {
    if (result == TEST_PASS) {
      printf("[PASS] %s\n", test_name);
      stats->passed++;
    } else {
      printf("[FAIL] %s\n", test_name);
      stats->failed++;
      strcpy(stats->last_test, test_name);
    }
  }
}

/* Test MPIPartition1D with 2 ranks */
int test_partition1d_2ranks(int rank, void *comm) {
  int nproc = 2;
  int nglobal = 100;

  /* Each rank should get 50 */
  int nlocal = MPIPartition1D(nglobal, nproc, rank);
  int expected = 50;

  if (nlocal != expected) {
    printf("    Rank %d: got %d, expected %d\n", rank, nlocal, expected);
    return TEST_FAIL;
  }

  /* Verify total using MPI reduction */
  int total = 0;
  MPISum_integer(&total, &nlocal, 1, comm);

  if (rank == 0 && total != nglobal) {
    printf("    Total %d != global %d\n", total, nglobal);
    return TEST_FAIL;
  }

  return TEST_PASS;
}

/* Test MPIPartition1D with uneven division */
int test_partition1d_uneven_2ranks(int rank, void *comm) {
  int nproc = 2;
  int nglobal = 101;

  int nlocal = MPIPartition1D(nglobal, nproc, rank);

  /* Rank 0 gets 50, rank 1 gets 51 */
  int expected = (rank == 0) ? 50 : 51;

  if (nlocal != expected) {
    printf("    Rank %d: got %d, expected %d\n", rank, nlocal, expected);
    return TEST_FAIL;
  }

  /* Verify total */
  int total = 0;
  MPISum_integer(&total, &nlocal, 1, comm);

  if (rank == 0 && total != nglobal) {
    printf("    Total %d != global %d\n", total, nglobal);
    return TEST_FAIL;
  }

  return TEST_PASS;
}

/* Test MPIRank1D and MPIRanknD conversion with 2 ranks */
int test_rank_conversion_2ranks(int rank, void *comm) {
  int ndims = 2;
  int iproc[2] = {2, 1};  /* 2 ranks along x, 1 along y */
  int ip[2];

  /* Convert my 1D rank to nD */
  MPIRanknD(ndims, rank, iproc, ip);

  /* Verify bounds */
  if (ip[0] < 0 || ip[0] >= iproc[0] || ip[1] < 0 || ip[1] >= iproc[1]) {
    printf("    Rank %d: nD rank (%d, %d) out of bounds\n", rank, ip[0], ip[1]);
    return TEST_FAIL;
  }

  /* Convert back to 1D */
  int rank_back = MPIRank1D(ndims, iproc, ip);

  if (rank_back != rank) {
    printf("    Rank %d: converted to (%d, %d) and back to %d\n",
           rank, ip[0], ip[1], rank_back);
    return TEST_FAIL;
  }

  return TEST_PASS;
}

/* Test MPIGetFilename generates unique names per rank */
int test_get_filename_2ranks(int rank, void *comm) {
  char root[] = "output";
  char filename[256];

  MPIGetFilename(root, comm, filename);

  /* Expected filenames: "output.0000" and "output.0001" */
  char expected[256];
  sprintf(expected, "output.%04d", rank);

  if (strcmp(filename, expected) != 0) {
    printf("    Rank %d: got '%s', expected '%s'\n", rank, filename, expected);
    return TEST_FAIL;
  }

  return TEST_PASS;
}

/* Test MPIMax_integer across 2 ranks */
int test_mpimax_integer_2ranks(int rank, void *comm) {
  int size = 5;
  int var[5];
  int global[5];

  /* Rank 0: [1, 3, 5, 7, 9], Rank 1: [2, 6, 4, 8, 10] */
  if (rank == 0) {
    var[0] = 1; var[1] = 3; var[2] = 5; var[3] = 7; var[4] = 9;
  } else {
    var[0] = 2; var[1] = 6; var[2] = 4; var[3] = 8; var[4] = 10;
  }

  MPIMax_integer(global, var, size, comm);

  /* Expected max: [2, 6, 5, 8, 10] */
  int expected[5] = {2, 6, 5, 8, 10};

  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (global[i] != expected[i]) {
      printf("    Rank %d: global[%d] = %d, expected %d\n",
             rank, i, global[i], expected[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Test MPIMax_double across 2 ranks */
int test_mpimax_double_2ranks(int rank, void *comm) {
  int size = 5;
  double var[5];
  double global[5];

  /* Rank 0: [1.5, 3.5, 5.5, 7.5, 9.5], Rank 1: [2.5, 6.5, 4.5, 8.5, 10.5] */
  if (rank == 0) {
    var[0] = 1.5; var[1] = 3.5; var[2] = 5.5; var[3] = 7.5; var[4] = 9.5;
  } else {
    var[0] = 2.5; var[1] = 6.5; var[2] = 4.5; var[3] = 8.5; var[4] = 10.5;
  }

  MPIMax_double(global, var, size, comm);

  /* Expected max: [2.5, 6.5, 5.5, 8.5, 10.5] */
  double expected[5] = {2.5, 6.5, 5.5, 8.5, 10.5};

  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (fabs(global[i] - expected[i]) > TOLERANCE) {
      printf("    Rank %d: global[%d] = %f, expected %f\n",
             rank, i, global[i], expected[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Test MPIMin_integer across 2 ranks */
int test_mpimin_integer_2ranks(int rank, void *comm) {
  int size = 5;
  int var[5];
  int global[5];

  /* Rank 0: [1, 3, 5, 7, 9], Rank 1: [2, 6, 4, 8, 10] */
  if (rank == 0) {
    var[0] = 1; var[1] = 3; var[2] = 5; var[3] = 7; var[4] = 9;
  } else {
    var[0] = 2; var[1] = 6; var[2] = 4; var[3] = 8; var[4] = 10;
  }

  MPIMin_integer(global, var, size, comm);

  /* Expected min: [1, 3, 4, 7, 9] */
  int expected[5] = {1, 3, 4, 7, 9};

  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (global[i] != expected[i]) {
      printf("    Rank %d: global[%d] = %d, expected %d\n",
             rank, i, global[i], expected[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Test MPIMin_double across 2 ranks */
int test_mpimin_double_2ranks(int rank, void *comm) {
  int size = 5;
  double var[5];
  double global[5];

  /* Rank 0: [1.5, 3.5, 5.5, 7.5, 9.5], Rank 1: [2.5, 6.5, 4.5, 8.5, 10.5] */
  if (rank == 0) {
    var[0] = 1.5; var[1] = 3.5; var[2] = 5.5; var[3] = 7.5; var[4] = 9.5;
  } else {
    var[0] = 2.5; var[1] = 6.5; var[2] = 4.5; var[3] = 8.5; var[4] = 10.5;
  }

  MPIMin_double(global, var, size, comm);

  /* Expected min: [1.5, 3.5, 4.5, 7.5, 9.5] */
  double expected[5] = {1.5, 3.5, 4.5, 7.5, 9.5};

  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (fabs(global[i] - expected[i]) > TOLERANCE) {
      printf("    Rank %d: global[%d] = %f, expected %f\n",
             rank, i, global[i], expected[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Test MPISum_integer across 2 ranks */
int test_mpisum_integer_2ranks(int rank, void *comm) {
  int size = 5;
  int var[5];
  int global[5];

  /* Rank 0: [1, 2, 3, 4, 5], Rank 1: [10, 20, 30, 40, 50] */
  if (rank == 0) {
    var[0] = 1; var[1] = 2; var[2] = 3; var[3] = 4; var[4] = 5;
  } else {
    var[0] = 10; var[1] = 20; var[2] = 30; var[3] = 40; var[4] = 50;
  }

  MPISum_integer(global, var, size, comm);

  /* Expected sum: [11, 22, 33, 44, 55] */
  int expected[5] = {11, 22, 33, 44, 55};

  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (global[i] != expected[i]) {
      printf("    Rank %d: global[%d] = %d, expected %d\n",
             rank, i, global[i], expected[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Test MPISum_double across 2 ranks */
int test_mpisum_double_2ranks(int rank, void *comm) {
  int size = 5;
  double var[5];
  double global[5];

  /* Rank 0: [1.1, 2.2, 3.3, 4.4, 5.5], Rank 1: [10.1, 20.2, 30.3, 40.4, 50.5] */
  if (rank == 0) {
    var[0] = 1.1; var[1] = 2.2; var[2] = 3.3; var[3] = 4.4; var[4] = 5.5;
  } else {
    var[0] = 10.1; var[1] = 20.2; var[2] = 30.3; var[3] = 40.4; var[4] = 50.5;
  }

  MPISum_double(global, var, size, comm);

  /* Expected sum: [11.2, 22.4, 33.6, 44.8, 56.0] */
  double expected[5] = {11.2, 22.4, 33.6, 44.8, 56.0};

  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (fabs(global[i] - expected[i]) > TOLERANCE) {
      printf("    Rank %d: global[%d] = %f, expected %f\n",
             rank, i, global[i], expected[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Test MPIMax in-place operation */
int test_mpimax_inplace_2ranks(int rank, void *comm) {
  int size = 3;
  double var[3];

  /* Rank 0: [1.0, 5.0, 3.0], Rank 1: [2.0, 4.0, 6.0] */
  if (rank == 0) {
    var[0] = 1.0; var[1] = 5.0; var[2] = 3.0;
  } else {
    var[0] = 2.0; var[1] = 4.0; var[2] = 6.0;
  }

  /* In-place operation */
  MPIMax_double(var, var, size, comm);

  /* Expected max: [2.0, 5.0, 6.0] */
  double expected[3] = {2.0, 5.0, 6.0};

  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (fabs(var[i] - expected[i]) > TOLERANCE) {
      printf("    Rank %d: var[%d] = %f, expected %f\n",
             rank, i, var[i], expected[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Test that verifies results are identical on all ranks */
int test_allreduce_consistency(int rank, void *comm) {
  int size = 4;
  double var[4];
  double global[4];

  /* Each rank has different values */
  for (int i = 0; i < size; i++) {
    var[i] = (double)(rank * 10 + i);
  }

  /* All ranks should get the same result */
  MPISum_double(global, var, size, comm);

  /* Exchange result with other rank to verify consistency */
  double other_global[4];

#ifndef serial
  MPI_Comm mpi_comm = *((MPI_Comm*)comm);
  int other_rank = (rank == 0) ? 1 : 0;

  MPI_Sendrecv(global, size, MPI_DOUBLE, other_rank, 0,
               other_global, size, MPI_DOUBLE, other_rank, 0,
               mpi_comm, MPI_STATUS_IGNORE);
#else
  for (int i = 0; i < size; i++) other_global[i] = global[i];
#endif

  /* Verify both ranks got the same result */
  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (fabs(global[i] - other_global[i]) > TOLERANCE) {
      printf("    Rank %d: global[%d] = %f, other rank has %f\n",
             rank, i, global[i], other_global[i]);
      test_result = TEST_FAIL;
    }
  }

  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  int rank = 0, nproc = 1;

#ifndef serial
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nproc);

  if (nproc != 2) {
    if (rank == 0) {
      printf("Error: This test must be run with exactly 2 MPI ranks\n");
      printf("Usage: mpirun -n 2 ./test_mpi_parallel\n");
    }
    MPI_Finalize();
    return 1;
  }
#else
  void *comm = NULL;
  if (rank == 0) {
    printf("Warning: Running in serial mode - MPI tests will be limited\n\n");
  }
#endif

  TestStats stats = {0, 0, ""};

  if (rank == 0) {
    printf("========================================\n");
    printf("MPI Functions Unit Tests (2 MPI Ranks)\n");
    printf("========================================\n\n");
  }

  if (rank == 0) printf("Testing Domain Partitioning:\n");
  print_test_result(&stats, "MPIPartition1D - 2 Ranks Even",
                    test_partition1d_2ranks(rank, &comm), rank);
  print_test_result(&stats, "MPIPartition1D - 2 Ranks Uneven",
                    test_partition1d_uneven_2ranks(rank, &comm), rank);

  if (rank == 0) printf("\nTesting Rank Conversions:\n");
  print_test_result(&stats, "Rank Conversion - 2 MPI Ranks",
                    test_rank_conversion_2ranks(rank, &comm), rank);

  if (rank == 0) printf("\nTesting Utility Functions:\n");
  print_test_result(&stats, "MPIGetFilename - Unique per Rank",
                    test_get_filename_2ranks(rank, &comm), rank);

  if (rank == 0) printf("\nTesting Reduction Operations (2 Ranks):\n");
  print_test_result(&stats, "MPIMax - Integer Array Reduction",
                    test_mpimax_integer_2ranks(rank, &comm), rank);
  print_test_result(&stats, "MPIMax - Double Array Reduction",
                    test_mpimax_double_2ranks(rank, &comm), rank);
  print_test_result(&stats, "MPIMin - Integer Array Reduction",
                    test_mpimin_integer_2ranks(rank, &comm), rank);
  print_test_result(&stats, "MPIMin - Double Array Reduction",
                    test_mpimin_double_2ranks(rank, &comm), rank);
  print_test_result(&stats, "MPISum - Integer Array Reduction",
                    test_mpisum_integer_2ranks(rank, &comm), rank);
  print_test_result(&stats, "MPISum - Double Array Reduction",
                    test_mpisum_double_2ranks(rank, &comm), rank);
  print_test_result(&stats, "MPIMax - In-place Operation",
                    test_mpimax_inplace_2ranks(rank, &comm), rank);
  print_test_result(&stats, "Allreduce - Result Consistency",
                    test_allreduce_consistency(rank, &comm), rank);

  if (rank == 0) {
    printf("\n========================================\n");
    printf("Test Results:\n");
    printf("  Passed: %d\n", stats.passed);
    printf("  Failed: %d\n", stats.failed);
    printf("========================================\n");

    if (stats.failed > 0) {
      printf("\nLast failed test: %s\n", stats.last_test);
    }
  }

#ifndef serial
  MPI_Finalize();
#endif

  return (stats.failed > 0) ? 1 : 0;
}
