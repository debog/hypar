check_PROGRAMS = test_math

test_math_SOURCES = test_math.c

test_math_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_math_LDADD = \
  $(top_builddir)/src/MathFunctions/libMathFunctions.a \
  $(top_builddir)/src/ArrayFunctions/libArrayFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_math
