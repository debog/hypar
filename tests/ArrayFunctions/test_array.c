/*! @file test_array.c
 *  @brief Unit tests for array functions
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>

#define TEST_PASS 0
#define TEST_FAIL 1
#define TOLERANCE 1e-12

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

/* Test _ArrayCopy1D_ macro */
int test_array_copy_1d() {
  int size = 10;
  double *x = (double*) malloc(size * sizeof(double));
  double *y = (double*) malloc(size * sizeof(double));
  
  for (int i = 0; i < size; i++) {
    x[i] = (double)(i + 1);
    y[i] = 0.0;
  }
  
  _ArrayCopy1D_(x, y, size);
  
  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (fabs(y[i] - x[i]) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d: expected %f, got %f\n", i, x[i], y[i]);
      break;
    }
  }
  
  free(x);
  free(y);
  return test_result;
}

/* Test _ArraySetValue_ macro */
int test_array_set_value() {
  int size = 20;
  double *x = (double*) malloc(size * sizeof(double));
  double value = 3.14159;
  
  _ArraySetValue_(x, size, value);
  
  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    if (fabs(x[i] - value) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d: expected %f, got %f\n", i, value, x[i]);
      break;
    }
  }
  
  free(x);
  return test_result;
}

/* Test _ArrayScale1D_ macro */
int test_array_scale_1d() {
  int size = 10;
  double *x = (double*) malloc(size * sizeof(double));
  double scale = 2.5;
  
  for (int i = 0; i < size; i++) {
    x[i] = (double)(i + 1);
  }
  
  double *original = (double*) malloc(size * sizeof(double));
  _ArrayCopy1D_(x, original, size);
  
  _ArrayScale1D_(x, scale, size);
  
  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    double expected = original[i] * scale;
    if (fabs(x[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d: expected %f, got %f\n", i, expected, x[i]);
      break;
    }
  }
  
  free(x);
  free(original);
  return test_result;
}

/* Test _ArrayAXPY_ macro (y = y + a*x) */
int test_array_axpy() {
  int size = 10;
  double *x = (double*) malloc(size * sizeof(double));
  double *y = (double*) malloc(size * sizeof(double));
  double a = 3.0;
  
  for (int i = 0; i < size; i++) {
    x[i] = (double)(i + 1);
    y[i] = (double)(2 * i);
  }
  
  double *y_orig = (double*) malloc(size * sizeof(double));
  _ArrayCopy1D_(y, y_orig, size);
  
  /* _ArrayAXPY_ is defined as: y[i] += a*x[i], so arguments are (x, a, y, size) */
  _ArrayAXPY_(x, a, y, size);
  
  int test_result = TEST_PASS;
  for (int i = 0; i < size; i++) {
    double expected = y_orig[i] + a * x[i];
    if (fabs(y[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d: expected %f, got %f\n", i, expected, y[i]);
      break;
    }
  }
  
  free(x);
  free(y);
  free(y_orig);
  return test_result;
}

/* Test _ArrayIndex1D_ and _ArrayIndexnD_ for 1D */
int test_array_index_1d() {
  int ndims = 1;
  int dim[1] = {10};
  int ghosts = 2;
  int index_nd[1];
  int index_1d;
  
  int test_result = TEST_PASS;
  
  /* Test converting from nD to 1D and back */
  for (int i = 0; i < dim[0]; i++) {
    index_nd[0] = i;
    _ArrayIndex1D_(ndims, dim, index_nd, ghosts, index_1d);
    
    int reconstructed[1];
    _ArrayIndexnD_(ndims, index_1d, dim, reconstructed, ghosts);
    
    if (reconstructed[0] != index_nd[0]) {
      test_result = TEST_FAIL;
      printf("    Failed for i=%d: got %d\n", i, reconstructed[0]);
      break;
    }
  }
  
  return test_result;
}

/* Test _ArrayIndex1D_ and _ArrayIndexnD_ for 2D */
int test_array_index_2d() {
  int ndims = 2;
  int dim[2] = {5, 8};
  int ghosts = 1;
  int index_nd[2];
  int index_1d;
  
  int test_result = TEST_PASS;
  
  /* Test converting from nD to 1D and back */
  for (int i = 0; i < dim[0]; i++) {
    for (int j = 0; j < dim[1]; j++) {
      index_nd[0] = i;
      index_nd[1] = j;
      _ArrayIndex1D_(ndims, dim, index_nd, ghosts, index_1d);
      
      int reconstructed[2];
      _ArrayIndexnD_(ndims, index_1d, dim, reconstructed, ghosts);
      
      if (reconstructed[0] != index_nd[0] || reconstructed[1] != index_nd[1]) {
        test_result = TEST_FAIL;
        printf("    Failed for (%d,%d): got (%d,%d)\n", 
               i, j, reconstructed[0], reconstructed[1]);
        break;
      }
    }
    if (test_result == TEST_FAIL) break;
  }
  
  return test_result;
}

/* Test _ArrayIndex1D_ and _ArrayIndexnD_ for 3D */
int test_array_index_3d() {
  int ndims = 3;
  int dim[3] = {4, 5, 6};
  int ghosts = 2;
  int index_nd[3];
  int index_1d;
  
  int test_result = TEST_PASS;
  
  /* Test a subset of points to keep test fast */
  for (int i = 0; i < dim[0]; i += 1) {
    for (int j = 0; j < dim[1]; j += 1) {
      for (int k = 0; k < dim[2]; k += 2) {
        index_nd[0] = i;
        index_nd[1] = j;
        index_nd[2] = k;
        _ArrayIndex1D_(ndims, dim, index_nd, ghosts, index_1d);
        
        int reconstructed[3];
        _ArrayIndexnD_(ndims, index_1d, dim, reconstructed, ghosts);
        
        if (reconstructed[0] != index_nd[0] || 
            reconstructed[1] != index_nd[1] || 
            reconstructed[2] != index_nd[2]) {
          test_result = TEST_FAIL;
          printf("    Failed for (%d,%d,%d): got (%d,%d,%d)\n", 
                 i, j, k, reconstructed[0], reconstructed[1], reconstructed[2]);
          break;
        }
      }
      if (test_result == TEST_FAIL) break;
    }
    if (test_result == TEST_FAIL) break;
  }
  
  return test_result;
}

/* Test _ArrayProduct1D_ macro */
int test_array_product_1d() {
  int size = 5;
  int x[5] = {2, 3, 4, 5, 6};
  int product;
  
  _ArrayProduct1D_(x, size, product);
  
  int expected = 2 * 3 * 4 * 5 * 6;
  if (product != expected) {
    printf("    Expected %d, got %d\n", expected, product);
    return TEST_FAIL;
  }
  
  return TEST_PASS;
}

/* Test ArrayCopynD function for 2D arrays */
int test_array_copy_nd_2d() {
  int ndims = 2;
  int nvars = 3;
  int dim[2] = {4, 5};
  int ghosts = 1;
  
  int size_with_ghosts = 1;
  for (int d = 0; d < ndims; d++) {
    size_with_ghosts *= (dim[d] + 2*ghosts);
  }
  
  double *x = (double*) calloc(size_with_ghosts * nvars, sizeof(double));
  double *y = (double*) calloc(size_with_ghosts * nvars, sizeof(double));
  int *index = (int*) calloc(ndims, sizeof(int));
  
  /* Initialize x with some values */
  for (int i = 0; i < size_with_ghosts * nvars; i++) {
    x[i] = (double)(i + 1);
  }
  
  int ierr = ArrayCopynD(ndims, x, y, dim, ghosts, ghosts, index, nvars);
  
  int test_result = TEST_PASS;
  if (ierr != 0) {
    printf("    ArrayCopynD returned error: %d\n", ierr);
    test_result = TEST_FAIL;
  } else {
    /* Check interior points (excluding ghosts) */
    int index[2];
    int done = 0;
    _ArraySetValue_(index, ndims, 0);
    while (!done) {
      int p;
      _ArrayIndex1D_(ndims, dim, index, ghosts, p);
      for (int v = 0; v < nvars; v++) {
        if (fabs(y[p*nvars+v] - x[p*nvars+v]) > TOLERANCE) {
          test_result = TEST_FAIL;
          printf("    Mismatch at index (%d,%d), var %d\n", index[0], index[1], v);
          break;
        }
      }
      if (test_result == TEST_FAIL) break;
      _ArrayIncrementIndex_(ndims, dim, index, done);
    }
  }
  
  free(x);
  free(y);
  free(index);
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};
  
  printf("========================================\n");
  printf("Array Functions Unit Tests\n");
  printf("========================================\n\n");
  
  printf("Testing Basic Array Operations:\n");
  print_test_result(&stats, "Array Copy 1D", test_array_copy_1d());
  print_test_result(&stats, "Array Set Value", test_array_set_value());
  print_test_result(&stats, "Array Scale 1D", test_array_scale_1d());
  print_test_result(&stats, "Array AXPY", test_array_axpy());
  
  printf("\nTesting Array Indexing:\n");
  print_test_result(&stats, "Array Index 1D", test_array_index_1d());
  print_test_result(&stats, "Array Index 2D", test_array_index_2d());
  print_test_result(&stats, "Array Index 3D", test_array_index_3d());
  print_test_result(&stats, "Array Product 1D", test_array_product_1d());
  
  printf("\nTesting Multi-dimensional Operations:\n");
  print_test_result(&stats, "Array Copy nD (2D)", test_array_copy_nd_2d());
  
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
