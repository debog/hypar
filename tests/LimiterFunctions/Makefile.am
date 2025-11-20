check_PROGRAMS = test_limiters

test_limiters_SOURCES = test_limiters.c

test_limiters_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

test_limiters_LDADD = \
  $(top_builddir)/src/LimiterFunctions/libLimiterFunctions.a \
  $(top_builddir)/src/CommonFunctions/libCommonFunctions.a

TESTS = test_limiters
