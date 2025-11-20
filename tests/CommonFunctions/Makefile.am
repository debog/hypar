check_PROGRAMS = test_common

test_common_SOURCES = test_common.c

test_common_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_common_LDADD = \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_common
