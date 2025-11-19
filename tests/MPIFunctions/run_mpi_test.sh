#!/bin/bash
# Wrapper script to run MPI test with 2 ranks

# Check if mpirun/mpiexec is available and determine MPI implementation
MPI_FLAGS=""
if command -v mpirun &> /dev/null; then
    MPI_CMD="mpirun"
    # Check if this is OpenMPI (supports --oversubscribe)
    if mpirun --version 2>&1 | grep -iq "open mpi"; then
        MPI_FLAGS="--oversubscribe"
    fi
elif command -v mpiexec &> /dev/null; then
    MPI_CMD="mpiexec"
else
    echo "Error: Neither mpirun nor mpiexec found"
    exit 1
fi

# Check if test executable exists
if [ ! -f ./test_mpi_parallel ]; then
    echo "Error: test_mpi_parallel executable not found"
    exit 1
fi

# Run the MPI tests with 2 ranks
$MPI_CMD $MPI_FLAGS -n 2 ./test_mpi_parallel
exit $?
