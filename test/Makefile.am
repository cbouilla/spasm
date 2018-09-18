AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = -lm
LDADD = $(top_builddir)/src/libspasm.la

TESTS = test_io.sh GFp vec_perm test_perm.sh test_gaxpy.sh test_reach.sh \
	test_triangular_solvers.sh test_lu.sh test_lu_solve.sh \
	test_transpose.sh test_matching.sh test_dm.sh \
	test_submatrix.sh test_cc.sh test_scc.sh test_dense_lu.sh

# known to fail
#test_schur.sh
#test_kernel.sh

check_PROGRAMS = io GFp gaxpy reach dense_usolve dense_lsolve \
	sparse_usolve lu permuted_lu lu_solve \
	vec_perm mat_perm transpose dm submatrix cc scc \
	kernel dense_lu matching schur


AM_TESTS_ENVIRONMENT = env AM_TAP_AWK='$(AWK)'
LOG_DRIVER = $(SHELL) $(top_srcdir)/build-aux/tap-driver.sh
EXTRA_DIST = $(TESTS) Expected Output

clean-local:
	rm -rf Output/*
