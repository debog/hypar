check_PROGRAMS = test_tridiag

test_tridiag_SOURCES = test_tridiag.c

test_tridiag_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_tridiag_LDADD = \
  $(top_builddir)/src/TridiagLU/libTridiagLU.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_tridiag
