/*! @file test_second_derivative.c
 *  @brief Unit tests for second derivative functions
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <secondderivative.h>
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

/* Test SecondDerivativeSecondOrderCentral on a quadratic function */
int test_second_derivative_second_order_quadratic() {
  int N = 20;
  int ghosts = 2;
  int nvars = 1;
  int dir = 0;

  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;

  double *f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *D2f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));

  /* Initialize with quadratic function f(x) = x^2 */
  for (int i = 0; i < N+2*ghosts; i++) {
    double x = (double)(i - ghosts);
    f[i] = x * x;
  }

  int ierr = SecondDerivativeSecondOrderCentral(D2f, f, dir, &solver, NULL);

  /* Second derivative of x^2 is 2 (not divided by dx^2)
   * Second order: D2f[i] = f[i+1] - 2*f[i] + f[i-1] = (i+1)^2 - 2*i^2 + (i-1)^2 = 2 */
  int test_result = TEST_PASS;
  for (int i = ghosts; i < N+ghosts; i++) {
    double expected = 2.0;
    if (fabs(D2f[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d: expected %f, got %f\n", i, expected, D2f[i]);
      break;
    }
  }

  free(f);
  free(D2f);
  return test_result;
}

/* Test SecondDerivativeSecondOrderCentral on a cubic function */
int test_second_derivative_second_order_cubic() {
  int N = 20;
  int ghosts = 2;
  int nvars = 1;
  int dir = 0;

  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;

  double *f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *D2f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));

  /* Initialize with cubic function f(x) = x^3 */
  for (int i = 0; i < N+2*ghosts; i++) {
    double x = (double)(i - ghosts);
    f[i] = x * x * x;
  }

  int ierr = SecondDerivativeSecondOrderCentral(D2f, f, dir, &solver, NULL);

  /* Second derivative of x^3 is 6*x
   * Second order: D2f[i] = f[i+1] - 2*f[i] + f[i-1] */
  int test_result = TEST_PASS;
  for (int i = ghosts; i < N+ghosts; i++) {
    double x = (double)(i - ghosts);
    double xp1 = x + 1.0;
    double xm1 = x - 1.0;
    double expected = xp1*xp1*xp1 - 2.0*x*x*x + xm1*xm1*xm1;
    if (fabs(D2f[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d (x=%f): expected %f, got %f\n", i, x, expected, D2f[i]);
      break;
    }
  }

  free(f);
  free(D2f);
  return test_result;
}

/* Test SecondDerivativeFourthOrderCentral on a quartic function */
int test_second_derivative_fourth_order_quartic() {
  int N = 20;
  int ghosts = 3;
  int nvars = 1;
  int dir = 0;

  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;

  double *f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *D2f = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));

  /* Initialize with quartic function f(x) = x^4 */
  for (int i = 0; i < N+2*ghosts; i++) {
    double x = (double)(i - ghosts);
    f[i] = x * x * x * x;
  }

  int ierr = SecondDerivativeFourthOrderCentral(D2f, f, dir, &solver, NULL);

  /* Fourth order: D2f[i] = (-f[i+2] + 16*f[i+1] - 30*f[i] + 16*f[i-1] - f[i-2])/12 */
  int test_result = TEST_PASS;
  for (int i = ghosts; i < N+ghosts; i++) {
    double x = (double)(i - ghosts);
    double xp2 = x + 2.0;
    double xp1 = x + 1.0;
    double xm1 = x - 1.0;
    double xm2 = x - 2.0;
    double expected = (-xp2*xp2*xp2*xp2 + 16.0*xp1*xp1*xp1*xp1 - 30.0*x*x*x*x
                      + 16.0*xm1*xm1*xm1*xm1 - xm2*xm2*xm2*xm2) / 12.0;
    if (fabs(D2f[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d (x=%f): expected %f, got %f, error=%e\n",
             i, x, expected, D2f[i], fabs(D2f[i] - expected));
      break;
    }
  }

  free(f);
  free(D2f);
  return test_result;
}

/* Test SecondDerivativeSecondOrderCentral on a 2D array */
int test_second_derivative_2d() {
  int ndims = 2;
  int nvars = 1;
  int dim_local[2] = {10, 12};
  int ghosts = 2;

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
  double *D2f_x = (double*) calloc(size_with_ghosts * nvars, sizeof(double));
  double *D2f_y = (double*) calloc(size_with_ghosts * nvars, sizeof(double));

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

  /* Compute second derivatives in x direction */
  int ierr = SecondDerivativeSecondOrderCentral(D2f_x, f, 0, &solver, NULL);

  /* Compute second derivatives in y direction */
  ierr = SecondDerivativeSecondOrderCentral(D2f_y, f, 1, &solver, NULL);

  /* Check derivatives: d2f/dx2 = 2, d2f/dy2 = 2 */
  int test_result = TEST_PASS;
  for (int i = 0; i < dim_local[0]; i++) {
    for (int j = 0; j < dim_local[1]; j++) {
      int index[2] = {i, j};
      int p;
      _ArrayIndex1D_(ndims, dim_local, index, ghosts, p);
      double expected_x = 2.0;
      double expected_y = 2.0;

      if (fabs(D2f_x[p] - expected_x) > TOLERANCE || fabs(D2f_y[p] - expected_y) > TOLERANCE) {
        test_result = TEST_FAIL;
        printf("    At (%d,%d): expected D2f_x=%f, got %f; expected D2f_y=%f, got %f\n",
               i, j, expected_x, D2f_x[p], expected_y, D2f_y[p]);
        break;
      }
    }
    if (test_result == TEST_FAIL) break;
  }

  free(f);
  free(D2f_x);
  free(D2f_y);
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};

  printf("========================================\n");
  printf("Second Derivative Functions Unit Tests\n");
  printf("========================================\n\n");

  printf("Testing Second Order Central Scheme:\n");
  print_test_result(&stats, "Second Order - Quadratic Function",
                    test_second_derivative_second_order_quadratic());
  print_test_result(&stats, "Second Order - Cubic Function",
                    test_second_derivative_second_order_cubic());

  printf("\nTesting Fourth Order Central Scheme:\n");
  print_test_result(&stats, "Fourth Order - Quartic Function",
                    test_second_derivative_fourth_order_quartic());

  printf("\nTesting Multi-dimensional:\n");
  print_test_result(&stats, "Second Order - 2D Array",
                    test_second_derivative_2d());

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
