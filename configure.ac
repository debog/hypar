AC_INIT([hypar],[4.1],[debojyoti.ghosh@gmail.com])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
echo -n AC_PACKAGE_VERSION>VERSION
AC_SUBST(PACKAGE_VERSION)
AC_MSG_NOTICE([hypar AC_PACKAGE_VERSION])

AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT('$(top_srcdir)')
includedir="${ac_pwd}/${srcdir}/include"
if ! test -d ${includedir} ; then
  includedir="${srcdir}/include"
fi

AC_ARG_ENABLE([serial],AS_HELP_STRING([--enable-serial],[Enable serial mode]))
AC_ARG_ENABLE([cuda],AS_HELP_STRING([--enable-cuda],[Enable CUDA]))
AC_ARG_WITH([cuda_dir],AS_HELP_STRING([--with-cuda-dir],[Specify path where CUDA is installed.]))
AC_ARG_WITH([mpi_dir],AS_HELP_STRING([--with-mpi-dir],[Specify path where MPI is installed.]))
AC_ARG_ENABLE([omp],AS_HELP_STRING([--enable-omp],[Enable OpenMP threads]))
AC_ARG_ENABLE([scalapack],AS_HELP_STRING([--enable-scalapack],[Enable ScaLAPACK]))
AC_ARG_WITH([blas_dir],AS_HELP_STRING([--with-blas-dir],[Specify path where BLAS libraries are installed.]))
AC_ARG_WITH([lapack_dir],AS_HELP_STRING([--with-lapack-dir],[Specify path where LAPACK libraries are installed.]))
AC_ARG_WITH([scalapack_dir],AS_HELP_STRING([--with-scalapack-dir],[Specify path where ScaLAPACK libraries are installed.]))
AC_ARG_ENABLE([fftw],AS_HELP_STRING([--enable-fftw],[Enable FFTW]))
AC_ARG_WITH([fftw_dir],AS_HELP_STRING([--with-fftw-dir],[Specify path where FFTW libraries are installed.]))
AC_ARG_ENABLE([python],AS_HELP_STRING([--enable-python],[Enable Python features (environment variables PYTHON_LIB_PATH, PYTHON_BIN_PATH, PYTHON_INCLUDE_PATH, NUMPY_INCLUDE_PATH, and PYTHON_LIB_NAME must be defined).]))

AC_PROG_CC
AC_PROG_CXX

CFLAGS="$CFLAGS -I$includedir"
if test "$GCC" = "yes" ; then
  CFLAGS="$CFLAGS -std=c99"
  LIBS="$LIBS -lm"
  MPI_CC=mpicc
else
  AC_MSG_WARN([Compilation with non-GNU compilers not supported.])
fi

CXXFLAGS="$CXXFLAGS -I$includedir"
if test "$GXX" = "yes" ; then
  CXXFLAGS="$CXXFLAGS -std=c++11"
  LIBS="$LIBS -lm"
  MPI_CXX=mpicxx
else
  AC_MSG_WARN([Compilation with non-GNU compilers not supported.])
fi

CUDA_CFLAGS="$CUDA_CFLAGS -I$includedir"

if test "x$enable_serial" = "xyes"; then
  CFLAGS="$CFLAGS -Dserial"
  CXXFLAGS="$CXXFLAGS -Dserial"
  CUDA_CFLAGS="$CUDA_CFLAGS -Dserial"
fi

if test "x$enable_omp" = "xyes" ; then
  AC_MSG_NOTICE([Compiling with OpenMP support.])
  if test "$CC" = "gcc"  ; then
    CFLAGS="$CFLAGS -fopenmp"
    LIBS="$LIBS -lgomp"
  fi
  CXXFLAGS="$CXXFLAGS -Dwith_omp"
  if test "$CXX" = "g++"  ; then
    CXXFLAGS="$CXXFLAGS -fopenmp"
    LIBS="$LIBS -lgomp"
  fi
  CXXFLAGS="$CXXFLAGS -Dwith_omp"
fi

if test "x$enable_serial" = "x"; then

  if test "x$with_mpi_dir" != "x" ; then
    AC_MSG_CHECKING([for ${MPI_CC} in $with_mpi_dir/bin/])
    if test -e $with_mpi_dir/bin/${MPI_CC} ; then
      AC_MSG_RESULT([yes])
      MPICC=$with_mpi_dir/bin/${MPI_CC}
      MPI_DIR=$with_mpi_dir
    else
      AC_MSG_RESULT([no])
    fi
    AC_MSG_CHECKING([for ${MPI_CXX} in $with_mpi_dir/bin/])
    if test -e $with_mpi_dir/bin/${MPI_CXX} ; then
      AC_MSG_RESULT([yes])
      MPICXX=$with_mpi_dir/bin/${MPI_CXX}
    else
      AC_MSG_RESULT([no])
    fi
  else
    if test "x$MPI_DIR" != "x" ; then
      AC_MSG_NOTICE([Found environment variable MPI_DIR=${MPI_DIR}.])
      AC_MSG_CHECKING([for ${MPI_CC} in ${MPI_DIR}/bin/])
      if test -e ${MPI_DIR}/bin/${MPI_CC} ; then
        AC_MSG_RESULT([yes])
        MPICC=${MPI_DIR}/bin/${MPI_CC}
      else
        AC_MSG_RESULT([no])
      fi
      AC_MSG_CHECKING([for ${MPI_CXX} in ${MPI_DIR}/bin/])
      if test -e ${MPI_DIR}/bin/${MPI_CXX} ; then
        AC_MSG_RESULT([yes])
        MPICXX=${MPI_DIR}/bin/${MPI_CXX}
      else
        AC_MSG_RESULT([no])
      fi
    fi
  fi

  if test "x$MPICC" = "x" ; then
    AC_CHECK_PROG(MPICC,${MPI_CC},${MPI_CC})
  fi
  if test "x$MPICXX" = "x" ; then
    AC_CHECK_PROG(MPICXX,${MPI_CXX},${MPI_CXX})
  fi

fi
  
AM_CONDITIONAL([ENABLE_CUDA], [test "x$enable_cuda" = "xyes"])
if test "x$enable_cuda" = "xyes"; then
  AC_MSG_NOTICE([Compiling with CUDA support.])
  NV_CC=nvcc

  #CFLAGS="$CFLAGS -DHAVE_CUDA -malign-double"
  #CXXFLAGS="$CXXFLAGS -DHAVE_CUDA -malign-double"
  CFLAGS="$CFLAGS -DHAVE_CUDA"
  CXXFLAGS="$CXXFLAGS -DHAVE_CUDA"

  CUDA_CFLAGS="$CUDA_CFLAGS -DHAVE_CUDA -lineinfo -Xptxas='-dlcm=ca -v'"
  if test "x$with_cuda_dir" != "x"; then
    AC_MSG_CHECKING([for ${NV_CC} in $with_cuda_dir/bin/])
    if test -e $with_cuda_dir/bin/${NV_CC}; then
      AC_MSG_RESULT([yes])
      NVCC=$with_cuda_dir/bin/${NV_CC}
    else
      AC_MSG_RESULT([no])
    fi
    LIBS="$LIBS -L$with_cuda_dir/lib64 -lcudart"
  else
    if test "x$CUDA_DIR" != "x"; then
      AC_MSG_NOTICE([Found environment variable CUDA_DIR=${CUDA_DIR}.])
      AC_MSG_CHECKING([for ${NV_CC} in ${CUDA_DIR}/bin/])
      if test -e ${CUDA_DIR}/bin/${NV_CC}; then
        AC_MSG_RESULT([yes])
        NVCC=${CUDA_DIR}/bin/${NV_CC}
      else
        AC_MSG_RESULT([no])
      fi
      LIBS="$LIBS -L$CUDA_DIR/lib64 -lcudart"
    else
      CUDA_DIR=$(which ${NVCC} | sed 's/\/bin\/nvcc//g')
      LIBS="$LIBS -L$CUDA_DIR/lib64 -lcudart"
    fi
  fi

  if test "x$NVCC" = "x"; then
    AC_CHECK_PROG(NVCC,${NV_CC},${NV_CC})
  fi

  if test "x$enable_serial" != "xyes"; then
    MPI_DIR=$(which ${MPICC} | sed 's/\/bin\/mpicc//g')
    if test -e ${MPI_DIR}/include/mpi.h ; then
      CUDA_CFLAGS="$CUDA_CFLAGS -I${MPI_DIR}/include -L${MPI_DIR}/lib -lmpi"
    else 
      CUDA_CFLAGS="$CUDA_CFLAGS -lmpi"
    fi
  fi

  AC_SUBST(CUDA_CFLAGS)
fi
  
if test "x$enable_serial" = "x"; then

  if test "x$MPICC" = "x" ; then
    AC_MSG_WARN([Cannot find ${MPI_CC}. Will compile with gcc. Use --with-mpi-dir to specify the location of ${MPI_CC}.])
    CFLAGS="$CFLAGS -Dserial"
    flag_serial="yes"
  else
    CC=$MPICC
  fi

  if test "x$MPICXX" = "x" ; then
    AC_MSG_WARN([Cannot find ${MPI_CXX}. Will compile with gcc. Use --with-mpi-dir to specify the location of ${MPI_CXX}.])
    CXXFLAGS="$CXXFLAGS -Dserial"
    flag_serial="yes"
  else
    CXX=$MPICXX
  fi

fi

if test "x$enable_scalapack" = "xyes" ; then
  AC_MSG_NOTICE([Compiling with ScaLAPACK support.])
  CFLAGS="$CFLAGS -Dwith_scalapack"
  CXXFLAGS="$CXXFLAGS -Dwith_scalapack"
  if test "x$with_scalapack_dir" != "x" ; then
    LDFLAGS="$LDFLAGS -Wl,-rpath,$with_scalapack_dir/lib -L$with_scalapack_dir/lib"
  fi
  if test "x$with_blas_dir" != "x" ; then
    LDFLAGS="$LDFLAGS -Wl,-rpath,$with_blas_dir/lib -L$with_blas_dir/lib"
  fi
  if test "x$with_lapack_dir" != "x" ; then
    LDFLAGS="$LDFLAGS -Wl,-rpath,$with_lapack_dir/lib -L$with_lapack_dir/lib"
  fi
  LIBS="$LIBS -lscalapack -lblas -llapack"
fi

if test "x$enable_fftw" = "xyes" ; then
  if test "x$enable_serial" = "x"; then
    AC_MSG_NOTICE([Compiling with FFTW support.])
    CFLAGS="$CFLAGS -Dfftw"
    CXXFLAGS="$CXXFLAGS -Dfftw"
    if test "x$with_fftw_dir" != "x" ; then
      CFLAGS="$CFLAGS -I$with_fftw_dir/include"
      CXXFLAGS="$CXXFLAGS -I$with_fftw_dir/include"
      LDFLAGS="$LDFLAGS -Wl,-rpath,$with_fftw_dir/lib -L$with_fftw_dir/lib"
    fi
    LIBS="$LIBS -lfftw3_mpi -lfftw3 -lm"
  else
    AC_MSG_WARN([Cannot enable FFTW when compiling in serial (--enable-serial)])
  fi
fi

if test "x$PETSC_DIR" = "x"; then
  AC_MSG_NOTICE([Environment variable PETSC_DIR not set. Compiling without PETSc interface.])
else
  if test [ ! -d "$PETSC_DIR" ]; then
    AC_MSG_NOTICE([Environment variable PETSC_DIR set, but PETSc not found there. Compiling without PETSc interface.])
  else
    if test "x$PETSC_ARCH" = "x"; then
      AC_MSG_NOTICE([Environment variable PETSC_ARCH not set. Compiling without PETSc interface.])
    else
      AC_MSG_NOTICE([Compiling with PETSc interface.])
cat <<MFILE > petscmake$$
PETSC_ARCH = ${PETSC_ARCH}
PETSC_DIR  = ${PETSC_DIR}
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules_util.mk
MFILE
      PETSC_INC="`make -f petscmake$$ getincludedirs`"
      PETSC_LIB="`make -f petscmake$$ getlinklibs`"
      rm petscmake$$
      CFLAGS="$CFLAGS -Dwith_petsc $PETSC_INC"
      CXXFLAGS="$CXXFLAGS -Dwith_petsc $PETSC_INC"
      LIBS="$LIBS $PETSC_LIB"
    fi
  fi
fi

if test "x$LIBROM_DIR" = "x"; then
  AC_MSG_NOTICE([Environment variable LIBROM_DIR not set. Compiling without libROM interface.])
else
  if test [ ! -d "$LIBROM_DIR" ]; then
    AC_MSG_NOTICE([Environment variable LIBROM_DIR set, but libROM not found there. Compiling without libROM interface.])
  else
    AC_MSG_NOTICE([Compiling with libROM interface.])
    CFLAGS="$CFLAGS -Dwith_librom -I$LIBROM_DIR/lib"
    CXXFLAGS="$CXXFLAGS -Dwith_librom -I$LIBROM_DIR/lib"
    LDFLAGS="$LDFLAGS -Wl,-rpath,$LIBROM_DIR/build/lib -L$LIBROM_DIR/build/lib"
    LIBS="$LIBS -lROM"
  fi
fi

AM_CONDITIONAL([ENABLE_PYTHON], [test "x$enable_python" = "xyes"])
if test "x$enable_python" = "xyes"; then
  if test "x$PYTHON_INCLUDE_PATH" != "x"; then
    AC_MSG_NOTICE([Environment variable PYTHON_INCLUDE_PATH found (${PYTHON_INCLUDE_PATH}).])
    if test "x$PYTHON_LIB_PATH" != "x"; then
      AC_MSG_NOTICE([Environment variable PYTHON_LIB_PATH found (${PYTHON_LIB_PATH}).])
      if test "x$PYTHON_LIB_NAME" != "x"; then
        AC_MSG_NOTICE([Environment variable PYTHON_LIB_NAME found (${PYTHON_LIB_NAME}).])
        if test "x$PYTHON_BIN_PATH" != "x"; then
          AC_MSG_NOTICE([Environment variable PYTHON_BIN_PATH found (${PYTHON_BIN_PATH}).])
          AC_MSG_NOTICE([Compiling with Python features.])
          CFLAGS="$CFLAGS -Dwith_python -I${PYTHON_INCLUDE_PATH}"
          CXXFLAGS="$CXXFLAGS -Dwith_python -I${PYTHON_INCLUDE_PATH}"
          LDFLAGS="$LDFLAGS -Wl,-rpath,$PYTHON_LIB_PATH -L$PYTHON_LIB_PATH"
          LIBS="$LIBS -${PYTHON_LIB_NAME} -lutil"
          if test "x$NUMPY_INCLUDE_PATH" != "x"; then
            AC_MSG_NOTICE([Environment variable NUMPY_INCLUDE_PATH found (${NUMPY_INCLUDE_PATH}).])
            AC_MSG_NOTICE([Compiling with Python-numpy features.])
            CFLAGS="$CFLAGS -Dwith_python_numpy -I${NUMPY_INCLUDE_PATH}"
            CXXFLAGS="$CXXFLAGS -Dwith_python_numpy -I${NUMPY_INCLUDE_PATH}"
          else
            AC_MSG_NOTICE([Environment variable NUMPY_INCLUDE_PATH not found. Skipping Python-numpy features.])
          fi
        else
          AC_MSG_NOTICE([Environment variable PYTHON_BIN_PATH not found. Skipping Python features.])
        fi
      else
        AC_MSG_NOTICE([Environment variable PYTHON_LIB_NAME not found. Skipping Python features.])
      fi
    else
      AC_MSG_NOTICE([Environment variable PYTHON_LIB_PATH not found. Skipping Python features.])
    fi
  else
    AC_MSG_NOTICE([Environment variable PYTHON_INCLUDE_PATH not found. Skipping Python features.])
  fi
fi

AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/ArrayFunctions/Makefile
 src/BandedMatrix/Makefile
 src/BoundaryConditions/Makefile
 src/CommonFunctions/Makefile
 src/FirstDerivative/Makefile
 src/HyParFunctions/Makefile
 src/ImmersedBoundaries/Makefile
 src/InterpolationFunctions/Makefile
 src/IOFunctions/Makefile
 src/LimiterFunctions/Makefile
 src/MathFunctions/Makefile
 src/MPIFunctions/Makefile
 src/PETScFunctions/Makefile
 src/PhysicalModels/Makefile
   src/PhysicalModels/LinearADR/Makefile
   src/PhysicalModels/Burgers/Makefile
   src/PhysicalModels/FPDoubleWell/Makefile
   src/PhysicalModels/FPPowerSystem/Makefile
   src/PhysicalModels/FPPowerSystem1Bus/Makefile
   src/PhysicalModels/FPPowerSystem3Bus/Makefile
   src/PhysicalModels/Euler1D/Makefile
   src/PhysicalModels/Euler2D/Makefile
   src/PhysicalModels/NavierStokes2D/Makefile
   src/PhysicalModels/NavierStokes3D/Makefile
   src/PhysicalModels/Numa2D/Makefile
   src/PhysicalModels/Numa3D/Makefile
   src/PhysicalModels/ShallowWater1D/Makefile
   src/PhysicalModels/ShallowWater2D/Makefile
   src/PhysicalModels/Vlasov/Makefile
 src/PlottingFunctions/Makefile
 src/ROM/Makefile
 src/SecondDerivative/Makefile
 src/Simulation/Makefile
 src/SparseGrids/Makefile
 src/TimeIntegration/Makefile
 src/TridiagLU/Makefile
])

AC_OUTPUT
