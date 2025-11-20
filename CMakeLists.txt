cmake_minimum_required(VERSION 3.12)
project(HyPar VERSION 4.1 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

option(ENABLE_SERIAL "Enable serial mode" OFF)
option(ENABLE_CUDA "Enable CUDA support" OFF)
option(ENABLE_OMP "Enable OpenMP threads" OFF)
option(ENABLE_SCALAPACK "Enable ScaLAPACK" OFF)
option(ENABLE_FFTW "Enable FFTW" OFF)
option(ENABLE_PYTHON "Enable Python features" OFF)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${CMAKE_SOURCE_DIR}/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_SOURCE_DIR}/include")

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    link_libraries(m)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    link_libraries(m)
endif()

if(ENABLE_SERIAL)
    add_compile_definitions(serial)
endif()

if(ENABLE_OMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_compile_definitions(with_omp)
        link_libraries(OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
    endif()
endif()

if(NOT ENABLE_SERIAL)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        include_directories(${MPI_C_INCLUDE_DIRS})
        include_directories(${MPI_CXX_INCLUDE_DIRS})
        link_libraries(${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
    endif()
endif()

if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_compile_definitions(HAVE_CUDA)
    # Disable min/max macros that conflict with C++ std::min/std::max
    add_compile_definitions(NOMINMAX)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DHAVE_CUDA -DNOMINMAX -lineinfo -Xptxas='-dlcm=ca -v'")
    link_libraries(CUDA::cudart)
endif()

if(ENABLE_SCALAPACK)
    add_compile_definitions(with_scalapack)
    # Try to find ScaLAPACK library with different possible names
    find_library(SCALAPACK_LIB
        NAMES scalapack scalapack-openmpi scalapack-mpich libscalapack.so
        PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib
        PATH_SUFFIXES openmpi mpich
    )
    if(SCALAPACK_LIB)
        message(STATUS "Found ScaLAPACK: ${SCALAPACK_LIB}")
        link_libraries(${SCALAPACK_LIB} blas lapack)
    else()
        message(STATUS "ScaLAPACK library not found, trying common names")
        link_libraries(scalapack-openmpi blas lapack)
    endif()
endif()

if(ENABLE_FFTW)
    if(NOT ENABLE_SERIAL)
        find_package(PkgConfig)
        if(PKG_CONFIG_FOUND)
            pkg_check_modules(FFTW3 fftw3)
            if(FFTW3_FOUND)
                add_compile_definitions(fftw)
                include_directories(${FFTW3_INCLUDE_DIRS})
                link_directories(${FFTW3_LIBRARY_DIRS})
                link_libraries(fftw3_mpi fftw3)
            endif()
        endif()
    else()
        message(WARNING "Cannot enable FFTW when compiling in serial mode")
    endif()
endif()

if(DEFINED ENV{PETSC_DIR} AND DEFINED ENV{PETSC_ARCH})
    if(EXISTS "$ENV{PETSC_DIR}")
        message(STATUS "Compiling with PETSc interface")
        add_compile_definitions(with_petsc)

        set(PETSC_MAKEFILE_CONTENT
"PETSC_ARCH = $ENV{PETSC_ARCH}
PETSC_DIR  = $ENV{PETSC_DIR}
include \${PETSC_DIR}/lib/petsc/conf/variables
include \${PETSC_DIR}/lib/petsc/conf/rules_util.mk
")
        file(WRITE "${CMAKE_BINARY_DIR}/petscmake.tmp" "${PETSC_MAKEFILE_CONTENT}")

        execute_process(
            COMMAND make --no-print-directory -f ${CMAKE_BINARY_DIR}/petscmake.tmp getincludedirs
            OUTPUT_VARIABLE PETSC_INCLUDES
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        execute_process(
            COMMAND make --no-print-directory -f ${CMAKE_BINARY_DIR}/petscmake.tmp getlinklibs
            OUTPUT_VARIABLE PETSC_LIBS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        file(REMOVE "${CMAKE_BINARY_DIR}/petscmake.tmp")

        string(STRIP "${PETSC_INCLUDES}" PETSC_INCLUDES)
        string(STRIP "${PETSC_LIBS}" PETSC_LIBS)

        separate_arguments(PETSC_INCLUDES_LIST UNIX_COMMAND "${PETSC_INCLUDES}")
        separate_arguments(PETSC_LIBS_LIST UNIX_COMMAND "${PETSC_LIBS}")

        foreach(inc ${PETSC_INCLUDES_LIST})
            if(inc MATCHES "^-I")
                string(SUBSTRING ${inc} 2 -1 inc_path)
                include_directories(${inc_path})
            endif()
        endforeach()

        foreach(lib ${PETSC_LIBS_LIST})
            if(lib MATCHES "^-L")
                string(SUBSTRING ${lib} 2 -1 lib_path)
                link_directories(${lib_path})
            elseif(lib MATCHES "^-l")
                string(SUBSTRING ${lib} 2 -1 lib_name)
                link_libraries(${lib_name})
            else()
                link_libraries(${lib})
            endif()
        endforeach()
    endif()
endif()

if(DEFINED ENV{LIBROM_DIR})
    if(EXISTS "$ENV{LIBROM_DIR}")
        message(STATUS "Compiling with libROM interface")
        add_compile_definitions(with_librom)
        include_directories($ENV{LIBROM_DIR}/lib)
        link_directories($ENV{LIBROM_DIR}/build/lib)
        link_libraries(ROM)
    endif()
endif()

if(ENABLE_PYTHON)
    if(DEFINED ENV{PYTHON_INCLUDE_PATH} AND DEFINED ENV{PYTHON_LIB_PATH} AND
       DEFINED ENV{PYTHON_LIB_NAME} AND DEFINED ENV{PYTHON_BIN_PATH})
        message(STATUS "Compiling with Python features")
        add_compile_definitions(with_python)
        include_directories($ENV{PYTHON_INCLUDE_PATH})
        link_directories($ENV{PYTHON_LIB_PATH})
        link_libraries($ENV{PYTHON_LIB_NAME} util)

        if(DEFINED ENV{NUMPY_INCLUDE_PATH})
            message(STATUS "Compiling with Python-numpy features")
            add_compile_definitions(with_python_numpy)
            include_directories($ENV{NUMPY_INCLUDE_PATH})
        endif()
    else()
        message(WARNING "PYTHON_INCLUDE_PATH, PYTHON_LIB_PATH, PYTHON_LIB_NAME, or PYTHON_BIN_PATH not set. Skipping Python features.")
    endif()
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/physicalmodels)

add_subdirectory(src)

# Add tests if BUILD_TESTING is enabled (default ON)
option(BUILD_TESTING "Build the testing tree" ON)
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

install(DIRECTORY include/ DESTINATION include)
install(FILES README.md DESTINATION share/doc/hypar)
