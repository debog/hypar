/*! @file test_limiters.c
 *  @brief Unit tests for slope limiter functions
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limiters.h>

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

/* Test MinMod limiter */
int test_minmod_limiter() {
  int test_result = TEST_PASS;
  
  /* MinMod properties:
   * phi(r) = max(0, min(1, r)) for r > 0
   * phi(r) = 0 for r <= 0 */
  
  double r, phi, expected;
  
  /* Test 1: r = 0.5 */
  r = 0.5;
  phi = LimiterMinMod(r);
  expected = 0.5;
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 2: r = 1.5 */
  r = 1.5;
  phi = LimiterMinMod(r);
  expected = 1.0;
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 3: r = -0.5 (opposite signs) */
  r = -0.5;
  phi = LimiterMinMod(r);
  expected = 0.0;
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 4: r = 1.0 */
  r = 1.0;
  phi = LimiterMinMod(r);
  expected = 1.0;
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  return test_result;
}

/* Test van Leer limiter */
int test_vanleer_limiter() {
  int test_result = TEST_PASS;
  
  /* van Leer: phi(r) = (r + |r|) / (1 + |r|) for r > 0
   *           phi(r) = 0 for r <= 0 */
  
  double r, phi, expected;
  
  /* Test 1: r = 1.0 */
  r = 1.0;
  phi = LimiterVanLeer(r);
  expected = (r + fabs(r)) / (1.0 + fabs(r));
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 2: r = 0.5 */
  r = 0.5;
  phi = LimiterVanLeer(r);
  expected = (r + fabs(r)) / (1.0 + fabs(r));
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 3: r = 2.0 */
  r = 2.0;
  phi = LimiterVanLeer(r);
  expected = (r + fabs(r)) / (1.0 + fabs(r));
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 4: r = -1.0 */
  r = -1.0;
  phi = LimiterVanLeer(r);
  expected = 0.0;
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  return test_result;
}

/* Test SuperBee limiter */
int test_superbee_limiter() {
  int test_result = TEST_PASS;
  
  /* SuperBee: phi(r) = max(0, min(2r, 1), min(r, 2)) */
  
  double r, phi, expected;
  
  /* Test 1: r = 0.5 */
  r = 0.5;
  phi = LimiterSuperBee(r);
  expected = 1.0;  /* max(0, min(1.0, 1), min(0.5, 2)) = max(0, 1, 0.5) = 1 */
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 2: r = 1.0 */
  r = 1.0;
  phi = LimiterSuperBee(r);
  expected = 1.0;  /* max(0, min(2.0, 1), min(1.0, 2)) = max(0, 1, 1) = 1 */
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 3: r = 2.0 */
  r = 2.0;
  phi = LimiterSuperBee(r);
  expected = 2.0;  /* max(0, min(4.0, 1), min(2.0, 2)) = max(0, 1, 2) = 2 */
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  /* Test 4: r = -0.5 */
  r = -0.5;
  phi = LimiterSuperBee(r);
  expected = 0.0;
  if (fabs(phi - expected) > TOLERANCE) {
    printf("    r=%.2f: expected %.4f, got %.4f\n", r, expected, phi);
    test_result = TEST_FAIL;
  }
  
  return test_result;
}

/* Test limiter properties: monotonicity */
int test_limiter_monotonicity() {
  int test_result = TEST_PASS;
  
  /* Test that limiters are monotone: phi(r) should increase with r for r > 0 */
  double r_values[] = {0.1, 0.5, 1.0, 1.5, 2.0, 3.0};
  int n = 6;
  
  /* MinMod */
  for (int i = 0; i < n-1; i++) {
    double phi1 = LimiterMinMod(r_values[i]);
    double phi2 = LimiterMinMod(r_values[i+1]);
    if (phi2 < phi1 - TOLERANCE) {
      printf("    MinMod: Not monotone at r=%.2f (phi1=%.4f) and r=%.2f (phi2=%.4f)\n",
             r_values[i], phi1, r_values[i+1], phi2);
      test_result = TEST_FAIL;
    }
  }
  
  /* van Leer */
  for (int i = 0; i < n-1; i++) {
    double phi1 = LimiterVanLeer(r_values[i]);
    double phi2 = LimiterVanLeer(r_values[i+1]);
    if (phi2 < phi1 - TOLERANCE) {
      printf("    VanLeer: Not monotone at r=%.2f (phi1=%.4f) and r=%.2f (phi2=%.4f)\n",
             r_values[i], phi1, r_values[i+1], phi2);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Test limiter bounds: 0 <= phi(r) <= 2r and 0 <= phi(r) <= 2 */
int test_limiter_bounds() {
  int test_result = TEST_PASS;
  
  double r_values[] = {0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
  int n = 7;
  
  for (int i = 0; i < n; i++) {
    double r = r_values[i];
    
    /* Test MinMod */
    double phi = LimiterMinMod(r);
    if (phi < -TOLERANCE || phi > 2.0 * r + TOLERANCE || phi > 2.0 + TOLERANCE) {
      printf("    MinMod violates bounds at r=%.2f: phi=%.4f\n", r, phi);
      test_result = TEST_FAIL;
    }
    
    /* Test van Leer */
    phi = LimiterVanLeer(r);
    if (phi < -TOLERANCE || phi > 2.0 * r + TOLERANCE || phi > 2.0 + TOLERANCE) {
      printf("    VanLeer violates bounds at r=%.2f: phi=%.4f\n", r, phi);
      test_result = TEST_FAIL;
    }
    
    /* Test SuperBee */
    phi = LimiterSuperBee(r);
    if (phi < -TOLERANCE || phi > 2.0 * r + TOLERANCE || phi > 2.0 + TOLERANCE) {
      printf("    SuperBee violates bounds at r=%.2f: phi=%.4f\n", r, phi);
      test_result = TEST_FAIL;
    }
  }
  
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};
  
  printf("========================================\n");
  printf("Limiter Functions Unit Tests\n");
  printf("========================================\n\n");
  
  printf("Testing Individual Limiters:\n");
  print_test_result(&stats, "MinMod Limiter", test_minmod_limiter());
  print_test_result(&stats, "van Leer Limiter", test_vanleer_limiter());
  print_test_result(&stats, "SuperBee Limiter", test_superbee_limiter());
  
  printf("\nTesting Limiter Properties:\n");
  print_test_result(&stats, "Monotonicity", test_limiter_monotonicity());
  print_test_result(&stats, "TVD Bounds", test_limiter_bounds());
  
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
