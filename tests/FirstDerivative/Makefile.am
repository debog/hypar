check_PROGRAMS = test_first_derivative

test_first_derivative_SOURCES = test_first_derivative.c

test_first_derivative_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_first_derivative_LDADD = \
  $(top_builddir)/src/FirstDerivative/libFirstDerivative.a \
  $(top_builddir)/src/ArrayFunctions/libArrayFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_first_derivative
