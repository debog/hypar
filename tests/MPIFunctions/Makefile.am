check_PROGRAMS = test_mpi

test_mpi_SOURCES = test_mpi.c

test_mpi_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_mpi_LDADD = \
  $(top_builddir)/src/MPIFunctions/libMPIFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_mpi
