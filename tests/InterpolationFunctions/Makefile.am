check_PROGRAMS = test_interpolation

test_interpolation_SOURCES = test_interpolation.c

test_interpolation_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_interpolation_LDADD = \
  $(top_builddir)/src/InterpolationFunctions/libInterpolationFunctions.a \
  $(top_builddir)/src/ArrayFunctions/libArrayFunctions.a \
  $(top_builddir)/src/MathFunctions/libMathFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_interpolation
