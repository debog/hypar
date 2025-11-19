# HyPar Unit Tests

This directory contains unit tests for HyPar components.

## Directory Structure

```
tests/
├── InterpolationFunctions/   # Tests for interpolation schemes
│   ├── test_interpolation.c  # Main test driver
│   └── Makefile.am           # Build configuration
├── Makefile.am               # Parent Makefile
└── README.md                 # This file
```

## InterpolationFunctions Tests

The interpolation functions tests verify the correctness of various spatial interpolation schemes used in HyPar, including:

### Tested Schemes

1. **First Order Upwind** (`Interp1PrimFirstOrderUpwind`)
   - Constant function test
   - Linear function test

2. **Second Order Central** (`Interp1PrimSecondOrderCentral`)
   - Linear function test (exact for 2nd order)
   - Quadratic function test (exact for 2nd order)

3. **Fourth Order Central** (`Interp1PrimFourthOrderCentral`)
   - Cubic function test

4. **Fifth Order WENO** (`Interp1PrimFifthOrderWENO`)
   - Smooth sine wave test (checks no spurious oscillations)

### Test Coverage

The tests verify:
- **Accuracy**: Interpolation schemes reproduce polynomial functions of appropriate order
- **Stability**: No spurious oscillations for smooth functions
- **Interface Values**: Correct reconstruction at cell interfaces

### Running the Tests

#### Quick Method (after building HyPar)
```bash
make check
```

#### Building and Testing from Scratch
```bash
# From the HyPar root directory
autoreconf -i
./configure --enable-serial  # or without --enable-serial for MPI
make
make check
```

#### Running Individual Tests
```bash
cd tests/InterpolationFunctions
./test_interpolation
```

### Test Output

Tests produce output in the format:
```
========================================
Interpolation Functions Unit Tests
========================================

Testing First Order Upwind Scheme:
[PASS] First Order Upwind - Constant Function
[PASS] First Order Upwind - Linear Function

Testing Second Order Central Scheme:
[PASS] Second Order Central - Linear Function
[PASS] Second Order Central - Quadratic Function
...

========================================
Test Results:
  Passed: 6
  Failed: 0
========================================
```

### Adding New Tests

To add new interpolation scheme tests:

1. Add a new test function in `test_interpolation.c`:
```c
int test_my_new_scheme() {
  // Setup test data
  // Call interpolation function
  // Verify results
  return TEST_PASS or TEST_FAIL;
}
```

2. Register the test in `main()`:
```c
print_test_result(&stats, "My New Scheme Test", test_my_new_scheme());
```

3. Rebuild and run:
```bash
make clean
make check
```

## Continuous Integration

Unit tests are automatically run by GitHub Actions on every push and pull request. See `.github/workflows/unit-tests.yml` for the CI configuration.

The CI runs tests in two configurations:
- Serial mode (--enable-serial)
- MPI mode (default)

## Future Test Additions

Planned additions include:
- Tests for MUSCL schemes (Second/Third Order)
- Tests for remaining WENO variants (CRWENO, HCWENO)
- Tests for characteristic-based interpolation
- Tests for compact schemes
- Multi-variable (nvars > 1) tests
- Multi-dimensional tests

## Notes

- Tests use a simple custom test harness (no external testing framework required)
- Tolerance for floating-point comparisons: `1e-10`
- Tests are designed to be fast and self-contained
- No MPI communication required for current tests (can run in serial mode)
