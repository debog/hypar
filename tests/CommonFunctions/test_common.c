/*! @file test_common.c
 *  @brief Unit tests for common utility functions
 *  @author Factory AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

/* Test GetStringFromInteger */
int test_get_string_from_integer() {
  char str[100];
  int test_result = TEST_PASS;

  /* Test 1: Convert 41 with width 5 */
  GetStringFromInteger(41, str, 5);
  if (strcmp(str, "00041") != 0) {
    printf("    Test 1 failed: expected '00041', got '%s'\n", str);
    test_result = TEST_FAIL;
  }

  /* Test 2: Convert 123 with width 3 */
  GetStringFromInteger(123, str, 3);
  if (strcmp(str, "123") != 0) {
    printf("    Test 2 failed: expected '123', got '%s'\n", str);
    test_result = TEST_FAIL;
  }

  /* Test 3: Convert 9 with width 1 */
  GetStringFromInteger(9, str, 1);
  if (strcmp(str, "9") != 0) {
    printf("    Test 3 failed: expected '9', got '%s'\n", str);
    test_result = TEST_FAIL;
  }

  /* Test 4: Convert 0 with width 4 */
  GetStringFromInteger(0, str, 4);
  if (strcmp(str, "0000") != 0) {
    printf("    Test 4 failed: expected '0000', got '%s'\n", str);
    test_result = TEST_FAIL;
  }

  return test_result;
}

/* Test takeLog function */
int test_take_log() {
  int n = 5;
  double *array = (double*) malloc(n * sizeof(double));

  /* Initialize with positive values */
  array[0] = 1.0;
  array[1] = 2.718281828;  /* e */
  array[2] = 10.0;
  array[3] = 100.0;
  array[4] = 0.5;

  /* Expected log values */
  double expected[5];
  expected[0] = log(1.0);
  expected[1] = log(2.718281828);
  expected[2] = log(10.0);
  expected[3] = log(100.0);
  expected[4] = log(0.5);

  takeLog(array, n);

  int test_result = TEST_PASS;
  for (int i = 0; i < n; i++) {
    if (fabs(array[i] - expected[i]) > TOLERANCE) {
      printf("    At index %d: expected %f, got %f\n", i, expected[i], array[i]);
      test_result = TEST_FAIL;
      break;
    }
  }

  free(array);
  return test_result;
}

/* Test takeExp function */
int test_take_exp() {
  int n = 5;
  double *array = (double*) malloc(n * sizeof(double));

  /* Initialize array */
  array[0] = 0.0;
  array[1] = 1.0;
  array[2] = -1.0;
  array[3] = 2.0;
  array[4] = -2.0;

  /* Expected exp values */
  double expected[5];
  expected[0] = exp(0.0);
  expected[1] = exp(1.0);
  expected[2] = exp(-1.0);
  expected[3] = exp(2.0);
  expected[4] = exp(-2.0);

  takeExp(array, n);

  int test_result = TEST_PASS;
  for (int i = 0; i < n; i++) {
    if (fabs(array[i] - expected[i]) > TOLERANCE) {
      printf("    At index %d: expected %f, got %f\n", i, expected[i], array[i]);
      test_result = TEST_FAIL;
      break;
    }
  }

  free(array);
  return test_result;
}

/* Test log and exp are inverse operations */
int test_log_exp_inverse() {
  int n = 10;
  double *array = (double*) malloc(n * sizeof(double));
  double *original = (double*) malloc(n * sizeof(double));

  /* Initialize with positive values */
  for (int i = 0; i < n; i++) {
    array[i] = (double)(i + 1) * 0.5;
    original[i] = array[i];
  }

  /* Take log then exp */
  takeLog(array, n);
  takeExp(array, n);

  /* Should get back original values */
  int test_result = TEST_PASS;
  for (int i = 0; i < n; i++) {
    if (fabs(array[i] - original[i]) > TOLERANCE) {
      printf("    At index %d: expected %f, got %f\n", i, original[i], array[i]);
      test_result = TEST_FAIL;
      break;
    }
  }

  free(array);
  free(original);
  return test_result;
}

/* Main test runner */
int main(int argc, char *argv[]) {
  TestStats stats = {0, 0, ""};

  printf("========================================\n");
  printf("Common Functions Unit Tests\n");
  printf("========================================\n\n");

  printf("Testing String Conversion:\n");
  print_test_result(&stats, "GetStringFromInteger",
                    test_get_string_from_integer());

  printf("\nTesting Mathematical Functions:\n");
  print_test_result(&stats, "Take Log", test_take_log());
  print_test_result(&stats, "Take Exp", test_take_exp());
  print_test_result(&stats, "Log-Exp Inverse", test_log_exp_inverse());

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
