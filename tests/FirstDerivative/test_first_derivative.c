/*! @file test_first_derivative.c
 *  @brief Unit tests for first derivative functions
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <firstderivative.h>
#include <hypar.h>

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

/* Test FirstDerivativeFirstOrder on a linear function */
int test_first_derivative_first_order_linear() {
  int N = 20;
  int ghosts = 1;
  int nvars = 1;
  int dir = 0;
  int bias = 1;

  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;

  double *f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *Df = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));

  /* Initialize with linear function f(x) = 2*x */
  for (int i = 0; i < N+2*ghosts; i++) {
    double x = (double)(i - ghosts);
    f[i] = 2.0 * x;
  }

  int ierr = FirstDerivativeFirstOrder(Df, f, dir, bias, &solver, NULL);

  /* First derivative of 2*x is 2 (not divided by dx) */
  /* For first order: Df[i] = f[i+1] - f[i] = 2*(i+1) - 2*i = 2 */
  int test_result = TEST_PASS;
  for (int i = ghosts; i < N+ghosts-1; i++) {
    double expected = 2.0;
    if (fabs(Df[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d: expected %f, got %f\n", i, expected, Df[i]);
      break;
    }
  }

  free(f);
  free(Df);
  return test_result;
}

/* Test FirstDerivativeSecondOrderCentral on a quadratic function */
int test_first_derivative_second_order_quadratic() {
  int N = 20;
  int ghosts = 2;
  int nvars = 1;
  int dir = 0;
  int bias = 0;

  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;

  double *f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *Df = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));

  /* Initialize with quadratic function f(x) = x^2 */
  for (int i = 0; i < N+2*ghosts; i++) {
    double x = (double)(i - ghosts);
    f[i] = x * x;
  }

  int ierr = FirstDerivativeSecondOrderCentral(Df, f, dir, bias, &solver, NULL);

  /* Second order central: Df[i] = (f[i+1] - f[i-1])/2 (not divided by dx)
   * For f(x) = x^2: df/dx = 2*x
   * Df[i] = ((i+1)^2 - (i-1)^2)/2 = (i^2+2i+1 - i^2+2i-1)/2 = 2i */
  int test_result = TEST_PASS;
  for (int i = ghosts; i < N+ghosts; i++) {
    double x = (double)(i - ghosts);
    double expected = 2.0 * x;
    if (fabs(Df[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d (x=%f): expected %f, got %f\n", i, x, expected, Df[i]);
      break;
    }
  }

  free(f);
  free(Df);
  return test_result;
}

/* Test FirstDerivativeFourthOrderCentral on a cubic function */
int test_first_derivative_fourth_order_cubic() {
  int N = 20;
  int ghosts = 3;
  int nvars = 1;
  int dir = 0;
  int bias = 0;

  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;

  double *f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *Df = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));

  /* Initialize with cubic function f(x) = x^3 */
  for (int i = 0; i < N+2*ghosts; i++) {
    double x = (double)(i - ghosts);
    f[i] = x * x * x;
  }

  int ierr = FirstDerivativeFourthOrderCentral(Df, f, dir, bias, &solver, NULL);

  /* Fourth order central: Df[i] = (-f[i+2] + 8*f[i+1] - 8*f[i-1] + f[i-2])/12
   * For f(x) = x^3: df/dx = 3*x^2 */
  int test_result = TEST_PASS;
  for (int i = ghosts; i < N+ghosts; i++) {
    double x = (double)(i - ghosts);
    double xm2 = x - 2.0;
    double xm1 = x - 1.0;
    double xp1 = x + 1.0;
    double xp2 = x + 2.0;
    double expected = (-xp2*xp2*xp2 + 8.0*xp1*xp1*xp1 - 8.0*xm1*xm1*xm1 + xm2*xm2*xm2) / 12.0;
    if (fabs(Df[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d (x=%f): expected %f, got %f, error=%e\n",
             i, x, expected, Df[i], fabs(Df[i] - expected));
      break;
    }
  }

  free(f);
  free(Df);
  return test_result;
}

/* Test FirstDerivativeSecondOrderCentral on a 2D array */
int test_first_derivative_2d() {
  int ndims = 2;
  int nvars = 1;
  int dim_local[2] = {10, 12};
  int ghosts = 2;
  int bias = 0;

  HyPar solver;
  solver.ndims = ndims;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  solver.dim_local = dim_local;

  int size_with_ghosts = 1;
  for (int d = 0; d < ndims; d++) {
    size_with_ghosts *= (dim_local[d] + 2*ghosts);
  }

  double *f = (double*) calloc(size_with_ghosts * nvars, sizeof(double));
  double *Df_x = (double*) calloc(size_with_ghosts * nvars, sizeof(double));
  double *Df_y = (double*) calloc(size_with_ghosts * nvars, sizeof(double));

  /* Initialize with function f(x,y) = x^2 + y^2, including ghost points */
  for (int i = -ghosts; i < dim_local[0]+ghosts; i++) {
    for (int j = -ghosts; j < dim_local[1]+ghosts; j++) {
      int index[2] = {i, j};
      int p;
      _ArrayIndex1D_(ndims, dim_local, index, ghosts, p);
      double x = (double)i;
      double y = (double)j;
      f[p] = x*x + y*y;
    }
  }

  /* Compute derivatives in x direction */
  int ierr = FirstDerivativeSecondOrderCentral(Df_x, f, 0, bias, &solver, NULL);

  /* Compute derivatives in y direction */
  ierr = FirstDerivativeSecondOrderCentral(Df_y, f, 1, bias, &solver, NULL);

  /* Check derivatives: df/dx = 2x, df/dy = 2y */
  int test_result = TEST_PASS;
  for (int i = 0; i < dim_local[0]; i++) {
    for (int j = 0; j < dim_local[1]; j++) {
      int index[2] = {i, j};
      int p;
      _ArrayIndex1D_(ndims, dim_local, index, ghosts, p);
      double x = (double)i;
      double y = (double)j;
      double expected_x = 2.0 * x;
      double expected_y = 2.0 * y;

      if (fabs(Df_x[p] - expected_x) > TOLERANCE || fabs(Df_y[p] - expected_y) > TOLERANCE) {
        test_result = TEST_FAIL;
        printf("    At (%d,%d): expected Df_x=%f, got %f; expected Df_y=%f, got %f\n",
               i, j, expected_x, Df_x[p], expected_y, Df_y[p]);
        break;
      }
    }
    if (test_result == TEST_FAIL) break;
  }

  free(f);
  free(Df_x);
  free(Df_y);
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};

  printf("========================================\n");
  printf("First Derivative Functions Unit Tests\n");
  printf("========================================\n\n");

  printf("Testing First Order Scheme:\n");
  print_test_result(&stats, "First Order - Linear Function",
                    test_first_derivative_first_order_linear());

  printf("\nTesting Second Order Central Scheme:\n");
  print_test_result(&stats, "Second Order - Quadratic Function",
                    test_first_derivative_second_order_quadratic());

  printf("\nTesting Fourth Order Central Scheme:\n");
  print_test_result(&stats, "Fourth Order - Cubic Function",
                    test_first_derivative_fourth_order_cubic());

  printf("\nTesting Multi-dimensional:\n");
  print_test_result(&stats, "Second Order - 2D Array",
                    test_first_derivative_2d());

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
