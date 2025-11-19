#!/bin/bash
# Wrapper script to run MPI test with 2 ranks

# Check if mpirun/mpiexec is available
if command -v mpirun &> /dev/null; then
    MPI_CMD="mpirun"
elif command -v mpiexec &> /dev/null; then
    MPI_CMD="mpiexec"
else
    echo "Error: Neither mpirun nor mpiexec found"
    exit 1
fi

# Run the MPI tests with 2 ranks
$MPI_CMD -n 2 ./test_mpi_parallel
TEST_RESULT=$?

# Also run the serial test with 1 rank if it exists
if [ -f ./test_mpi ]; then
    echo ""
    echo "========================================="
    echo "Running serial unit tests (1 rank)..."
    echo "========================================="
    $MPI_CMD -n 1 ./test_mpi
    SERIAL_RESULT=$?
    
    # Return failure if either test failed
    if [ $TEST_RESULT -ne 0 ] || [ $SERIAL_RESULT -ne 0 ]; then
        exit 1
    fi
else
    # Only parallel test exists
    exit $TEST_RESULT
fi

exit 0
