check_PROGRAMS = test_second_derivative

test_second_derivative_SOURCES = test_second_derivative.c

test_second_derivative_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_second_derivative_LDADD = \
  $(top_builddir)/src/SecondDerivative/libSecondDerivative.a \
  $(top_builddir)/src/ArrayFunctions/libArrayFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_second_derivative
