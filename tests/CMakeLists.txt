# Unit tests for HyPar

enable_testing()

# Helper function to add a test executable
function(add_hypar_test test_name test_source)
    add_executable(${test_name} ${test_source})
    target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${test_name} ${ARGN})
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

# InterpolationFunctions tests
add_hypar_test(test_interpolation
    InterpolationFunctions/test_interpolation.c
    InterpolationFunctions
    ArrayFunctions
    CommonFunctions
)

# ArrayFunctions tests
add_hypar_test(test_array
    ArrayFunctions/test_array.c
    ArrayFunctions
    CommonFunctions
)

# FirstDerivative tests
add_hypar_test(test_first_derivative
    FirstDerivative/test_first_derivative.c
    FirstDerivative
    ArrayFunctions
    CommonFunctions
)

# SecondDerivative tests
add_hypar_test(test_second_derivative
    SecondDerivative/test_second_derivative.c
    SecondDerivative
    ArrayFunctions
    CommonFunctions
)

# TridiagLU tests
add_hypar_test(test_tridiag
    TridiagLU/test_tridiag.c
    TridiagLU
    CommonFunctions
)

# CommonFunctions tests
add_hypar_test(test_common
    CommonFunctions/test_common.c
    CommonFunctions
)

# LimiterFunctions tests
add_hypar_test(test_limiters
    LimiterFunctions/test_limiters.c
    LimiterFunctions
    CommonFunctions
)

# MathFunctions tests
add_hypar_test(test_math
    MathFunctions/test_math.c
    MathFunctions
    ArrayFunctions
    CommonFunctions
)

# MPIFunctions tests (parallel)
if(NOT ENABLE_SERIAL)
    add_executable(test_mpi_parallel
        MPIFunctions/test_mpi_parallel.c
    )
    target_include_directories(test_mpi_parallel PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(test_mpi_parallel
        MPIFunctions
        CommonFunctions
    )

    # Add MPI test with 2 processes
    if(MPIEXEC_EXECUTABLE)
        add_test(NAME test_mpi_parallel
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_mpi_parallel>
            ${MPIEXEC_POSTFLAGS}
        )
    else()
        # Fallback if MPIEXEC not found - try common names
        find_program(MPIRUN_EXECUTABLE NAMES mpirun mpiexec)
        if(MPIRUN_EXECUTABLE)
            add_test(NAME test_mpi_parallel
                COMMAND ${MPIRUN_EXECUTABLE} -n 2 $<TARGET_FILE:test_mpi_parallel>
            )
        else()
            message(WARNING "MPI executable (mpirun/mpiexec) not found. MPIFunctions tests will not be added.")
        endif()
    endif()
else()
    message(STATUS "MPIFunctions tests disabled in serial mode")
endif()
