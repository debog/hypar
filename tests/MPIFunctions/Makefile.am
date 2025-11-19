check_PROGRAMS = test_mpi_parallel

test_mpi_parallel_SOURCES = test_mpi_parallel.c

test_mpi_parallel_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_mpi_parallel_LDADD = \
  $(top_builddir)/src/MPIFunctions/libMPIFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = run_mpi_test.sh

EXTRA_DIST = run_mpi_test.sh test_mpi.c
