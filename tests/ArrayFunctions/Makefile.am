check_PROGRAMS = test_array

test_array_SOURCES = test_array.c

test_array_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_array_LDADD = \
  $(top_builddir)/src/ArrayFunctions/libArrayFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_array
