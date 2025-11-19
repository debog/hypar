/*! @file test_tridiag.c
 *  @brief Unit tests for tridiagonal matrix solver functions
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tridiagLU.h>

#define TEST_PASS 0
#define TEST_FAIL 1
#define TOLERANCE 1e-8  /* Relaxed tolerance for tridiag solver */

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

/* Test tridiagonal LU solve for a simple system */
int test_tridiag_lu_simple() {
  int n = 5;  /* System size */

  /* Allocate arrays for tridiagonal matrix
   * a = subdiagonal, b = diagonal, c = superdiagonal */
  double *a = (double*) calloc(n, sizeof(double));
  double *b = (double*) calloc(n, sizeof(double));
  double *c = (double*) calloc(n, sizeof(double));
  double *x = (double*) calloc(n, sizeof(double));
  double *rhs = (double*) calloc(n, sizeof(double));

  /* Set up a simple tridiagonal system:
   * [2  1  0  0  0] [x0]   [1]
   * [1  2  1  0  0] [x1]   [2]
   * [0  1  2  1  0] [x2] = [3]
   * [0  0  1  2  1] [x3]   [4]
   * [0  0  0  1  2] [x4]   [5] */

  for (int i = 0; i < n; i++) {
    b[i] = 2.0;  /* diagonal */
    if (i < n-1) c[i] = 1.0;  /* superdiagonal */
    if (i > 0) a[i] = 1.0;  /* subdiagonal */
    rhs[i] = (double)(i + 1);
  }

  /* Save original matrix and rhs for verification (tridiagLU modifies inputs) */
  double *a_orig = (double*) calloc(n, sizeof(double));
  double *b_orig = (double*) calloc(n, sizeof(double));
  double *c_orig = (double*) calloc(n, sizeof(double));
  double *rhs_orig = (double*) calloc(n, sizeof(double));
  for (int i = 0; i < n; i++) {
    a_orig[i] = a[i];
    b_orig[i] = b[i];
    c_orig[i] = c[i];
    rhs_orig[i] = rhs[i];
  }

  /* Initialize TridiagLU structure */
  TridiagLU params;
  params.evaluate_norm = 0;
  params.maxiter = 10;
  params.atol = 1e-12;
  params.rtol = 1e-10;
  strcpy(params.reducedsolvetype, _TRIDIAG_GS_);

  /* Solve the system (note: tridiagLU modifies a, b, c, and rhs becomes the solution) */
  int ierr = tridiagLU(a, b, c, rhs, n, 1, &params, NULL);

  /* Copy solution from rhs to x */
  for (int i = 0; i < n; i++) {
    x[i] = rhs[i];
  }

  /* Verify solution by computing A*x and comparing with original rhs */
  int test_result = TEST_PASS;
  for (int i = 0; i < n; i++) {
    double Ax_i = b_orig[i] * x[i];
    if (i > 0) Ax_i += a_orig[i] * x[i-1];
    if (i < n-1) Ax_i += c_orig[i] * x[i+1];

    if (fabs(Ax_i - rhs_orig[i]) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At row %d: A*x = %f, expected rhs = %f\n", i, Ax_i, rhs_orig[i]);
      break;
    }
  }

  free(a_orig);
  free(b_orig);
  free(c_orig);
  free(rhs_orig);
  free(a);
  free(b);
  free(c);
  free(x);
  free(rhs);
  return test_result;
}

/* Test tridiagonal LU solve for a diagonally dominant system */
int test_tridiag_lu_diag_dominant() {
  int n = 10;

  double *a = (double*) calloc(n, sizeof(double));
  double *b = (double*) calloc(n, sizeof(double));
  double *c = (double*) calloc(n, sizeof(double));
  double *x = (double*) calloc(n, sizeof(double));
  double *rhs = (double*) calloc(n, sizeof(double));

  /* Set up a diagonally dominant system:
   * b[i] = 5, a[i] = c[i] = 1 */
  for (int i = 0; i < n; i++) {
    b[i] = 5.0;
    if (i < n-1) c[i] = 1.0;
    if (i > 0) a[i] = 1.0;
    rhs[i] = (double)(2 * i + 1);
  }

  /* Save original for verification */
  double *a_orig = (double*) calloc(n, sizeof(double));
  double *b_orig = (double*) calloc(n, sizeof(double));
  double *c_orig = (double*) calloc(n, sizeof(double));
  double *rhs_orig = (double*) calloc(n, sizeof(double));
  for (int i = 0; i < n; i++) {
    a_orig[i] = a[i];
    b_orig[i] = b[i];
    c_orig[i] = c[i];
    rhs_orig[i] = rhs[i];
  }

  TridiagLU params;
  params.evaluate_norm = 0;
  params.maxiter = 10;
  params.atol = 1e-12;
  params.rtol = 1e-10;
  strcpy(params.reducedsolvetype, _TRIDIAG_GS_);

  int ierr = tridiagLU(a, b, c, rhs, n, 1, &params, NULL);

  /* Copy solution from rhs to x */
  for (int i = 0; i < n; i++) {
    x[i] = rhs[i];
  }

  /* Verify solution */
  int test_result = TEST_PASS;
  for (int i = 0; i < n; i++) {
    double Ax_i = b_orig[i] * x[i];
    if (i > 0) Ax_i += a_orig[i] * x[i-1];
    if (i < n-1) Ax_i += c_orig[i] * x[i+1];

    if (fabs(Ax_i - rhs_orig[i]) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At row %d: A*x = %f, expected rhs = %f, error=%e\n",
             i, Ax_i, rhs_orig[i], fabs(Ax_i - rhs_orig[i]));
      break;
    }
  }

  free(a_orig);
  free(b_orig);
  free(c_orig);
  free(rhs_orig);
  free(a);
  free(b);
  free(c);
  free(x);
  free(rhs);
  return test_result;
}

/* Test tridiagonal system with known solution */
int test_tridiag_lu_known_solution() {
  int n = 8;

  double *a = (double*) calloc(n, sizeof(double));
  double *b = (double*) calloc(n, sizeof(double));
  double *c = (double*) calloc(n, sizeof(double));
  double *x = (double*) calloc(n, sizeof(double));
  double *rhs = (double*) calloc(n, sizeof(double));
  double *x_exact = (double*) calloc(n, sizeof(double));

  /* Set up a system with known solution x[i] = i + 1 */
  for (int i = 0; i < n; i++) {
    b[i] = 4.0;
    if (i < n-1) c[i] = 1.0;
    if (i > 0) a[i] = 1.0;
    x_exact[i] = (double)(i + 1);
  }

  /* Compute rhs = A * x_exact */
  for (int i = 0; i < n; i++) {
    rhs[i] = b[i] * x_exact[i];
    if (i > 0) rhs[i] += a[i] * x_exact[i-1];
    if (i < n-1) rhs[i] += c[i] * x_exact[i+1];
  }

  /* Solve the system (tridiagLU modifies inputs) */
  TridiagLU params;
  params.evaluate_norm = 0;
  params.maxiter = 10;
  params.atol = 1e-12;
  params.rtol = 1e-10;
  strcpy(params.reducedsolvetype, _TRIDIAG_GS_);

  int ierr = tridiagLU(a, b, c, rhs, n, 1, &params, NULL);

  /* Copy solution from rhs to x */
  for (int i = 0; i < n; i++) {
    x[i] = rhs[i];
  }

  /* Compare with exact solution */
  int test_result = TEST_PASS;
  for (int i = 0; i < n; i++) {
    if (fabs(x[i] - x_exact[i]) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At index %d: computed x = %f, expected = %f\n", i, x[i], x_exact[i]);
      break;
    }
  }

  free(a);
  free(b);
  free(c);
  free(x);
  free(rhs);
  free(x_exact);
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};

  printf("========================================\n");
  printf("Tridiagonal LU Solver Unit Tests\n");
  printf("========================================\n\n");

  printf("Testing Tridiagonal LU Solver:\n");
  print_test_result(&stats, "Simple 5x5 System",
                    test_tridiag_lu_simple());
  print_test_result(&stats, "Diagonally Dominant 10x10 System",
                    test_tridiag_lu_diag_dominant());
  print_test_result(&stats, "Known Solution Test",
                    test_tridiag_lu_known_solution());

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
