AM_CPPFLAGS = -I$(top_srcdir)/src
LDADD = $(top_builddir)/src/libspasm.la

TESTS = test_io.sh GFp vec_perm test_perm.sh test_gaxpy.sh test_reach.sh \
	test_triangular_solvers.sh test_lu.sh test_pluq.sh

check_PROGRAMS = io GFp gaxpy reach dense_solve sparse_solve lu vec_perm mat_perm pluq

AM_TESTS_ENVIRONMENT = \
	env AM_TAP_AWK='$(AWK)'
LOG_DRIVER = $(SHELL) $(top_srcdir)/build-aux/tap-driver.sh
EXTRA_DIST = $(TESTS) Expected


clean-local:
	rm -rf Output/*
