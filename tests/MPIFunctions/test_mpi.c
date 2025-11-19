/*! @file test_mpi.c
 *  @brief Unit tests for MPI utility functions (serial mode)
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

void print_test_result(TestStats *stats, const char *test_name, int result) {
  if (result == TEST_PASS) {
    printf("[PASS] %s\n", test_name);
    stats->passed++;
  } else {
    printf("[FAIL] %s\n", test_name);
    stats->failed++;
    strcpy(stats->last_test, test_name);
  }
}

/* Test MPIPartition1D with evenly divisible size */
int test_partition1d_even() {
  int nglobal = 100;
  int nproc = 4;
  int expected_size = 25;
  
  int test_result = TEST_PASS;
  
  /* All ranks should get equal size when evenly divisible */
  for (int rank = 0; rank < nproc; rank++) {
    int nlocal = MPIPartition1D(nglobal, nproc, rank);
    if (nlocal != expected_size) {
      printf("    Rank %d: got %d, expected %d\n", rank, nlocal, expected_size);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIPartition1D with non-evenly divisible size */
int test_partition1d_uneven() {
  int nglobal = 103;
  int nproc = 4;
  
  int test_result = TEST_PASS;
  int total = 0;
  
  /* First 3 ranks get floor(103/4)=25, last rank gets 25+3=28 */
  for (int rank = 0; rank < nproc; rank++) {
    int nlocal = MPIPartition1D(nglobal, nproc, rank);
    total += nlocal;
    
    if (rank < nproc - 1) {
      if (nlocal != 25) {
        printf("    Rank %d: got %d, expected 25\n", rank, nlocal);
        test_result = TEST_FAIL;
      }
    } else {
      if (nlocal != 28) {
        printf("    Rank %d: got %d, expected 28\n", rank, nlocal);
        test_result = TEST_FAIL;
      }
    }
  }
  
  /* Verify total equals global size */
  if (total != nglobal) {
    printf("    Total size %d != global size %d\n", total, nglobal);
    test_result = TEST_FAIL;
  }
  
  return test_result;
}

/* Test MPIPartition1D with single rank */
int test_partition1d_single() {
  int nglobal = 50;
  int nproc = 1;
  int rank = 0;
  
  int nlocal = MPIPartition1D(nglobal, nproc, rank);
  
  if (nlocal != nglobal) {
    printf("    Single rank: got %d, expected %d\n", nlocal, nglobal);
    return TEST_FAIL;
  }
  
  return TEST_PASS;
}

/* Test MPIRank1D and MPIRanknD bidirectional conversion (1D case) */
int test_rank_conversion_1d() {
  int ndims = 1;
  int iproc[1] = {8};  /* 8 ranks along 1 dimension */
  
  int test_result = TEST_PASS;
  
  for (int rank1d = 0; rank1d < iproc[0]; rank1d++) {
    int ip[1];
    
    /* Convert 1D to nD */
    MPIRanknD(ndims, rank1d, iproc, ip);
    
    /* Convert back to 1D */
    int rank1d_back = MPIRank1D(ndims, iproc, ip);
    
    if (rank1d_back != rank1d) {
      printf("    1D rank %d -> nD (%d) -> 1D %d (mismatch)\n", 
             rank1d, ip[0], rank1d_back);
      test_result = TEST_FAIL;
    }
    
    /* For 1D, ip[0] should equal rank1d */
    if (ip[0] != rank1d) {
      printf("    1D rank %d -> nD (%d) != %d\n", rank1d, ip[0], rank1d);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIRank1D and MPIRanknD bidirectional conversion (2D case) */
int test_rank_conversion_2d() {
  int ndims = 2;
  int iproc[2] = {4, 3};  /* 4 ranks along x, 3 along y = 12 total ranks */
  int total_ranks = iproc[0] * iproc[1];
  
  int test_result = TEST_PASS;
  
  for (int rank1d = 0; rank1d < total_ranks; rank1d++) {
    int ip[2];
    
    /* Convert 1D to nD */
    MPIRanknD(ndims, rank1d, iproc, ip);
    
    /* Verify nD rank is within bounds */
    if (ip[0] < 0 || ip[0] >= iproc[0] || ip[1] < 0 || ip[1] >= iproc[1]) {
      printf("    1D rank %d -> nD (%d, %d) out of bounds\n", 
             rank1d, ip[0], ip[1]);
      test_result = TEST_FAIL;
    }
    
    /* Convert back to 1D */
    int rank1d_back = MPIRank1D(ndims, iproc, ip);
    
    if (rank1d_back != rank1d) {
      printf("    1D rank %d -> nD (%d, %d) -> 1D %d (mismatch)\n", 
             rank1d, ip[0], ip[1], rank1d_back);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIRank1D and MPIRanknD bidirectional conversion (3D case) */
int test_rank_conversion_3d() {
  int ndims = 3;
  int iproc[3] = {2, 3, 4};  /* 2x3x4 = 24 total ranks */
  int total_ranks = iproc[0] * iproc[1] * iproc[2];
  
  int test_result = TEST_PASS;
  
  for (int rank1d = 0; rank1d < total_ranks; rank1d++) {
    int ip[3];
    
    /* Convert 1D to nD */
    MPIRanknD(ndims, rank1d, iproc, ip);
    
    /* Verify nD rank is within bounds */
    for (int d = 0; d < ndims; d++) {
      if (ip[d] < 0 || ip[d] >= iproc[d]) {
        printf("    1D rank %d -> nD (%d, %d, %d) dimension %d out of bounds\n", 
               rank1d, ip[0], ip[1], ip[2], d);
        test_result = TEST_FAIL;
      }
    }
    
    /* Convert back to 1D */
    int rank1d_back = MPIRank1D(ndims, iproc, ip);
    
    if (rank1d_back != rank1d) {
      printf("    1D rank %d -> nD (%d, %d, %d) -> 1D %d (mismatch)\n", 
             rank1d, ip[0], ip[1], ip[2], rank1d_back);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIGetFilename */
int test_get_filename() {
  char root[] = "output";
  char filename[256];
  
  /* In serial mode, rank is always 0 */
  MPIGetFilename(root, NULL, filename);
  
  /* Expected: "output.0000" */
  if (strcmp(filename, "output.0000") != 0) {
    printf("    Got '%s', expected 'output.0000'\n", filename);
    return TEST_FAIL;
  }
  
  return TEST_PASS;
}

/* Test MPIMax_integer in serial mode */
int test_mpimax_integer() {
  int size = 5;
  int var[5] = {3, 7, 2, 9, 5};
  int global[5];
  
  MPIMax_integer(global, var, size, NULL);
  
  int test_result = TEST_PASS;
  /* In serial mode, global should equal var */
  for (int i = 0; i < size; i++) {
    if (global[i] != var[i]) {
      printf("    global[%d] = %d, expected %d\n", i, global[i], var[i]);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIMax_double in serial mode */
int test_mpimax_double() {
  int size = 5;
  double var[5] = {3.14, 7.21, 2.71, 9.99, 5.55};
  double global[5];
  
  MPIMax_double(global, var, size, NULL);
  
  int test_result = TEST_PASS;
  /* In serial mode, global should equal var */
  for (int i = 0; i < size; i++) {
    if (fabs(global[i] - var[i]) > TOLERANCE) {
      printf("    global[%d] = %f, expected %f\n", i, global[i], var[i]);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIMin_integer in serial mode */
int test_mpimin_integer() {
  int size = 5;
  int var[5] = {3, 7, 2, 9, 5};
  int global[5];
  
  MPIMin_integer(global, var, size, NULL);
  
  int test_result = TEST_PASS;
  /* In serial mode, global should equal var */
  for (int i = 0; i < size; i++) {
    if (global[i] != var[i]) {
      printf("    global[%d] = %d, expected %d\n", i, global[i], var[i]);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIMin_double in serial mode */
int test_mpimin_double() {
  int size = 5;
  double var[5] = {3.14, 7.21, 2.71, 9.99, 5.55};
  double global[5];
  
  MPIMin_double(global, var, size, NULL);
  
  int test_result = TEST_PASS;
  /* In serial mode, global should equal var */
  for (int i = 0; i < size; i++) {
    if (fabs(global[i] - var[i]) > TOLERANCE) {
      printf("    global[%d] = %f, expected %f\n", i, global[i], var[i]);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPISum_integer in serial mode */
int test_mpisum_integer() {
  int size = 5;
  int var[5] = {3, 7, 2, 9, 5};
  int global[5];
  
  MPISum_integer(global, var, size, NULL);
  
  int test_result = TEST_PASS;
  /* In serial mode, global should equal var */
  for (int i = 0; i < size; i++) {
    if (global[i] != var[i]) {
      printf("    global[%d] = %d, expected %d\n", i, global[i], var[i]);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPISum_double in serial mode */
int test_mpisum_double() {
  int size = 5;
  double var[5] = {3.14, 7.21, 2.71, 9.99, 5.55};
  double global[5];
  
  MPISum_double(global, var, size, NULL);
  
  int test_result = TEST_PASS;
  /* In serial mode, global should equal var */
  for (int i = 0; i < size; i++) {
    if (fabs(global[i] - var[i]) > TOLERANCE) {
      printf("    global[%d] = %f, expected %f\n", i, global[i], var[i]);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test MPIMax in-place operation */
int test_mpimax_inplace() {
  int size = 5;
  double var[5] = {3.14, 7.21, 2.71, 9.99, 5.55};
  double original[5];
  
  /* Save original values */
  for (int i = 0; i < size; i++) {
    original[i] = var[i];
  }
  
  /* In-place operation: var is both input and output */
  MPIMax_double(var, var, size, NULL);
  
  int test_result = TEST_PASS;
  /* In serial mode, values should remain unchanged */
  for (int i = 0; i < size; i++) {
    if (fabs(var[i] - original[i]) > TOLERANCE) {
      printf("    var[%d] = %f, expected %f\n", i, var[i], original[i]);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test rank ordering property in 2D */
int test_rank_ordering_2d() {
  int ndims = 2;
  int iproc[2] = {3, 4};  /* 3x4 = 12 ranks */
  
  int test_result = TEST_PASS;
  
  /* Test specific known mappings */
  /* Rank 0 should be (0,0) */
  int ip[2];
  MPIRanknD(ndims, 0, iproc, ip);
  if (ip[0] != 0 || ip[1] != 0) {
    printf("    Rank 0 -> (%d, %d), expected (0, 0)\n", ip[0], ip[1]);
    test_result = TEST_FAIL;
  }
  
  /* Rank 1 should be (1,0) - x changes first */
  MPIRanknD(ndims, 1, iproc, ip);
  if (ip[0] != 1 || ip[1] != 0) {
    printf("    Rank 1 -> (%d, %d), expected (1, 0)\n", ip[0], ip[1]);
    test_result = TEST_FAIL;
  }
  
  /* Rank 3 should be (0,1) - wraps to next y */
  MPIRanknD(ndims, 3, iproc, ip);
  if (ip[0] != 0 || ip[1] != 1) {
    printf("    Rank 3 -> (%d, %d), expected (0, 1)\n", ip[0], ip[1]);
    test_result = TEST_FAIL;
  }
  
  /* Rank 11 should be (2,3) - last rank */
  MPIRanknD(ndims, 11, iproc, ip);
  if (ip[0] != 2 || ip[1] != 3) {
    printf("    Rank 11 -> (%d, %d), expected (2, 3)\n", ip[0], ip[1]);
    test_result = TEST_FAIL;
  }
  
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};
  
  printf("========================================\n");
  printf("MPI Functions Unit Tests (Serial Mode)\n");
  printf("========================================\n\n");
  
  printf("Testing Domain Partitioning:\n");
  print_test_result(&stats, "MPIPartition1D - Even Division", 
                    test_partition1d_even());
  print_test_result(&stats, "MPIPartition1D - Uneven Division", 
                    test_partition1d_uneven());
  print_test_result(&stats, "MPIPartition1D - Single Rank", 
                    test_partition1d_single());
  
  printf("\nTesting Rank Conversions:\n");
  print_test_result(&stats, "Rank Conversion - 1D Bidirectional", 
                    test_rank_conversion_1d());
  print_test_result(&stats, "Rank Conversion - 2D Bidirectional", 
                    test_rank_conversion_2d());
  print_test_result(&stats, "Rank Conversion - 3D Bidirectional", 
                    test_rank_conversion_3d());
  print_test_result(&stats, "Rank Ordering - 2D Grid Layout", 
                    test_rank_ordering_2d());
  
  printf("\nTesting Utility Functions:\n");
  print_test_result(&stats, "MPIGetFilename - Serial Mode", 
                    test_get_filename());
  
  printf("\nTesting Reduction Operations (Serial Mode):\n");
  print_test_result(&stats, "MPIMax - Integer Array", 
                    test_mpimax_integer());
  print_test_result(&stats, "MPIMax - Double Array", 
                    test_mpimax_double());
  print_test_result(&stats, "MPIMin - Integer Array", 
                    test_mpimin_integer());
  print_test_result(&stats, "MPIMin - Double Array", 
                    test_mpimin_double());
  print_test_result(&stats, "MPISum - Integer Array", 
                    test_mpisum_integer());
  print_test_result(&stats, "MPISum - Double Array", 
                    test_mpisum_double());
  print_test_result(&stats, "MPIMax - In-place Operation", 
                    test_mpimax_inplace());
  
  printf("\n========================================\n");
  printf("Test Results:\n");
  printf("  Passed: %d\n", stats.passed);
  printf("  Failed: %d\n", stats.failed);
  printf("========================================\n");
  
  if (stats.failed > 0) {
    printf("\nLast failed test: %s\n", stats.last_test);
    return 1;
  }
  
  return 0;
}
