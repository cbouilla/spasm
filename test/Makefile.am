AM_CPPFLAGS = -I$(top_srcdir)/src
LDADD = $(top_builddir)/src/libspasm.la

TESTS = test_io.sh vec_perm test_perm.sh

AM_TESTS_ENVIRONMENT = \
	env AM_TAP_AWK='$(AWK)'
LOG_DRIVER = $(SHELL) $(top_srcdir)/build-aux/tap-driver.sh
EXTRA_DIST = $(TESTS) Expected

check_PROGRAMS = io vec_perm mat_perm

clean-local:
	rm -rf Output/*
