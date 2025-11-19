/*! @file test_interpolation.c
 *  @brief Unit tests for interpolation functions
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
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

/* Test FirstOrderUpwind interpolation with constant function */
int test_first_order_upwind_constant() {
  int N = 10;
  int ghosts = 1;
  int nvars = 1;
  int upw = 1;
  int dir = 0;
  
  /* Allocate arrays */
  double *fC = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *fI = (double*) calloc((N+1)*nvars, sizeof(double));
  double *u = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *x = (double*) calloc(N+2*ghosts, sizeof(double));
  
  /* Setup HyPar structure */
  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;
  
  /* Setup MPIVariables structure (for serial) */
  MPIVariables mpi;
  mpi.rank = 0;
  mpi.nproc = 1;
  
  /* Initialize with constant function f(x) = 2.0 */
  for (int i = 0; i < (N+2*ghosts)*nvars; i++) {
    fC[i] = 2.0;
    u[i] = 2.0;
  }
  
  /* Call interpolation function */
  int ierr = Interp1PrimFirstOrderUpwind(fI, fC, u, x, upw, dir, &solver, &mpi, 1);
  
  /* Check result - should be constant 2.0 at all interfaces */
  int test_result = TEST_PASS;
  for (int i = 0; i < (N+1)*nvars; i++) {
    if (fabs(fI[i] - 2.0) > TOLERANCE) {
      test_result = TEST_FAIL;
      break;
    }
  }
  
  free(fC);
  free(fI);
  free(u);
  free(x);
  
  return test_result;
}

/* Test FirstOrderUpwind with linear function */
int test_first_order_upwind_linear() {
  int N = 10;
  int ghosts = 1;
  int nvars = 1;
  int upw = 1;
  int dir = 0;
  
  double *fC = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *fI = (double*) calloc((N+1)*nvars, sizeof(double));
  double *u = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *x = (double*) calloc(N+2*ghosts, sizeof(double));
  
  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;
  
  MPIVariables mpi;
  mpi.rank = 0;
  mpi.nproc = 1;
  
  /* Initialize with linear function f(x) = x */
  for (int i = 0; i < N+2*ghosts; i++) {
    fC[i] = (double)(i - ghosts);
    u[i] = fC[i];
  }
  
  int ierr = Interp1PrimFirstOrderUpwind(fI, fC, u, x, upw, dir, &solver, &mpi, 1);
  
  /* For upwind=1 (left biased), interface values should be from left cell */
  int test_result = TEST_PASS;
  for (int i = 0; i < N+1; i++) {
    double expected = (double)(i - 1);
    if (fabs(fI[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      break;
    }
  }
  
  free(fC);
  free(fI);
  free(u);
  free(x);
  
  return test_result;
}

/* Test SecondOrderCentral interpolation with linear function */
int test_second_order_central_linear() {
  int N = 10;
  int ghosts = 2;
  int nvars = 1;
  int upw = 1;
  int dir = 0;
  
  double *fC = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *fI = (double*) calloc((N+1)*nvars, sizeof(double));
  double *u = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *x = (double*) calloc(N+2*ghosts, sizeof(double));
  
  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;
  
  MPIVariables mpi;
  mpi.rank = 0;
  mpi.nproc = 1;
  
  /* Initialize with linear function f(x) = x */
  for (int i = 0; i < N+2*ghosts; i++) {
    fC[i] = (double)(i - ghosts);
    u[i] = fC[i];
  }
  
  int ierr = Interp1PrimSecondOrderCentral(fI, fC, u, x, upw, dir, &solver, &mpi, 1);
  
  /* Second order central: fI[i] = (fC[i-1] + fC[i]) / 2 */
  /* For linear function f(x) = x, this gives exact interface values */
  int test_result = TEST_PASS;
  for (int i = 0; i < N+1; i++) {
    /* Interface i is between cells i-1 and i (in the interior numbering) */
    double fL = (double)(i - 1);
    double fR = (double)(i);
    double expected = (fL + fR) / 2.0;
    if (fabs(fI[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At interface %d: got %f, expected %f\n", i, fI[i], expected);
      break;
    }
  }
  
  free(fC);
  free(fI);
  free(u);
  free(x);
  
  return test_result;
}

/* Test SecondOrderCentral with quadratic function */
int test_second_order_central_quadratic() {
  int N = 10;
  int ghosts = 2;
  int nvars = 1;
  int upw = 1;
  int dir = 0;
  
  double *fC = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *fI = (double*) calloc((N+1)*nvars, sizeof(double));
  double *u = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *x = (double*) calloc(N+2*ghosts, sizeof(double));
  
  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;
  
  MPIVariables mpi;
  mpi.rank = 0;
  mpi.nproc = 1;
  
  /* Initialize with quadratic function f(x) = x^2 */
  for (int i = 0; i < N+2*ghosts; i++) {
    double xi = (double)(i - ghosts);
    fC[i] = xi * xi;
    u[i] = fC[i];
  }
  
  int ierr = Interp1PrimSecondOrderCentral(fI, fC, u, x, upw, dir, &solver, &mpi, 1);
  
  /* Second order central: fI[i] = (fC[i-1] + fC[i]) / 2 */
  /* For quadratic, this gives the average of two cell values */
  int test_result = TEST_PASS;
  for (int i = 0; i < N+1; i++) {
    double xL = (double)(i - 1);
    double xR = (double)(i);
    double expected = (xL*xL + xR*xR) / 2.0;
    if (fabs(fI[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At interface %d: got %f, expected %f\n", i, fI[i], expected);
      break;
    }
  }
  
  free(fC);
  free(fI);
  free(u);
  free(x);
  
  return test_result;
}

/* Test FourthOrderCentral with cubic function */
int test_fourth_order_central_cubic() {
  int N = 10;
  int ghosts = 2;
  int nvars = 1;
  int upw = 1;
  int dir = 0;
  
  double *fC = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *fI = (double*) calloc((N+1)*nvars, sizeof(double));
  double *u = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *x = (double*) calloc(N+2*ghosts, sizeof(double));
  
  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;
  
  MPIVariables mpi;
  mpi.rank = 0;
  mpi.nproc = 1;
  
  /* Initialize with cubic function f(x) = x^3 */
  for (int i = 0; i < N+2*ghosts; i++) {
    double xi = (double)(i - ghosts);
    fC[i] = xi * xi * xi;
    u[i] = fC[i];
  }
  
  int ierr = Interp1PrimFourthOrderCentral(fI, fC, u, x, upw, dir, &solver, &mpi, 1);
  
  /* Fourth order central: fI[i] = -1/12*fC[i-2] + 7/12*fC[i-1] + 7/12*fC[i] - 1/12*fC[i+1] */
  int test_result = TEST_PASS;
  for (int i = 1; i < N; i++) {  /* Skip boundaries where stencil extends outside */
    double xLL = (double)(i - 2);
    double xL = (double)(i - 1);
    double xR = (double)(i);
    double xRR = (double)(i + 1);
    double expected = (-1.0/12.0)*(xLL*xLL*xLL) + (7.0/12.0)*(xL*xL*xL) 
                    + (7.0/12.0)*(xR*xR*xR) + (-1.0/12.0)*(xRR*xRR*xRR);
    if (fabs(fI[i] - expected) > TOLERANCE) {
      test_result = TEST_FAIL;
      printf("    At interface %d: got %f, expected %f, error=%e\n", 
             i, fI[i], expected, fabs(fI[i] - expected));
      break;
    }
  }
  
  free(fC);
  free(fI);
  free(u);
  free(x);
  
  return test_result;
}

/* Test smooth sine wave - tests WENO limiting */
int test_weno_smooth_function() {
  int N = 20;
  int ghosts = 3;
  int nvars = 1;
  int upw = 1;
  int dir = 0;
  
  double *fC = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *fI = (double*) calloc((N+1)*nvars, sizeof(double));
  double *u = (double*) calloc((N+2*ghosts)*nvars, sizeof(double));
  double *x = (double*) calloc(N+2*ghosts, sizeof(double));
  
  HyPar solver;
  solver.ndims = 1;
  solver.nvars = nvars;
  solver.ghosts = ghosts;
  int dim_local[1] = {N};
  solver.dim_local = dim_local;
  int stride_with_ghosts[1] = {1};
  solver.stride_with_ghosts = stride_with_ghosts;
  
  /* Initialize WENO parameters */
  WENOParameters weno;
  weno.mapped = 0;
  weno.borges = 0;
  weno.yc = 0;
  weno.no_limiting = 0;
  weno.eps = 1e-6;
  weno.p = 2.0;
  weno.tol = 1e-10;
  weno.rc = 0.3;
  weno.xi = 0.001;
  
  /* Allocate and initialize WENO weights and offset arrays */
  weno.offset = (int*) calloc(solver.ndims, sizeof(int));
  weno.offset[0] = 0;
  
  /* Calculate total size for weight arrays */
  int total_size = nvars * (N + 1);  /* For 1D: nvars * (dim+1) */
  weno.size = total_size;
  
  /* Allocate weight arrays (4*size for upwind/conservative flags) */
  weno.w1 = (double*) calloc(4*total_size, sizeof(double));
  weno.w2 = (double*) calloc(4*total_size, sizeof(double));
  weno.w3 = (double*) calloc(4*total_size, sizeof(double));
  
  /* Initialize weights to optimal values */
  for (int i = 0; i < 4*total_size; i++) {
    weno.w1[i] = 0.1;  /* Optimal weight 1 */
    weno.w2[i] = 0.6;  /* Optimal weight 2 */
    weno.w3[i] = 0.3;  /* Optimal weight 3 */
  }
  
  solver.interp = &weno;
  solver.SetInterpLimiterVar = NULL;
  
  MPIVariables mpi;
  mpi.rank = 0;
  mpi.nproc = 1;
  
  /* Initialize with smooth sine wave f(x) = sin(2*pi*x/N) */
  double pi = 4.0 * atan(1.0);
  for (int i = 0; i < N+2*ghosts; i++) {
    double xi = 2.0 * pi * (double)(i - ghosts) / (double)N;
    fC[i] = sin(xi);
    u[i] = fC[i];
  }
  
  int ierr = Interp1PrimFifthOrderWENO(fI, fC, u, x, upw, dir, &solver, &mpi, 1);
  
  /* For smooth function, WENO should behave close to 5th order upwind */
  /* Just check that interpolation is reasonable and doesn't produce spurious oscillations */
  int test_result = TEST_PASS;
  for (int i = 1; i < N; i++) {
    /* Check that values are within reasonable bounds */
    if (fabs(fI[i]) > 1.5) {  /* sin(x) should be bounded by 1 */
      test_result = TEST_FAIL;
      printf("    Interface %d has unreasonable value: %f\n", i, fI[i]);
      break;
    }
  }
  
  /* Clean up allocated memory */
  free(weno.w1);
  free(weno.w2);
  free(weno.w3);
  free(weno.offset);
  free(fC);
  free(fI);
  free(u);
  free(x);
  
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};
  
  printf("========================================\n");
  printf("Interpolation Functions Unit Tests\n");
  printf("========================================\n\n");
  
  printf("Testing First Order Upwind Scheme:\n");
  print_test_result(&stats, "First Order Upwind - Constant Function", 
                    test_first_order_upwind_constant());
  print_test_result(&stats, "First Order Upwind - Linear Function", 
                    test_first_order_upwind_linear());
  
  printf("\nTesting Second Order Central Scheme:\n");
  print_test_result(&stats, "Second Order Central - Linear Function", 
                    test_second_order_central_linear());
  print_test_result(&stats, "Second Order Central - Quadratic Function", 
                    test_second_order_central_quadratic());
  
  printf("\nTesting Fourth Order Central Scheme:\n");
  print_test_result(&stats, "Fourth Order Central - Cubic Function", 
                    test_fourth_order_central_cubic());
  
  printf("\nTesting WENO Scheme:\n");
  print_test_result(&stats, "Fifth Order WENO - Smooth Sine Wave", 
                    test_weno_smooth_function());
  
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
