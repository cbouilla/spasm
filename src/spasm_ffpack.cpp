/*
 * This file encapsulates all the invocations of FFLAS-FFPACK
 * and exports simple C functions.
 */
#include <stdio.h>

#include <givaro/modular-balanced.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/ffpack/ffpack.h>

extern "C" {
#include "spasm.h"
}


/*
 * Q of size m (#cols). On output, Q[0:rank] is singificant.
 * It is implicit that U[i, 0:Q[i]] == 0 and U[i, 0:Q[i]] == 1.
 * M[i, Q[i]:m] contain the actual data.
 */
template<typename T>
static inline int spasm_ffpack_rref(i64 prime, int n, int m, T *A, int ldA, size_t *qinv)
{
	double start = spasm_wtime();
	fprintf(stderr, "[ffpack/rref] Matrix of dimension %d x %d mod %" PRId64"... ", n, m, prime);
	fflush(stderr);
	Givaro::ModularBalanced<T> GFp(prime);
	size_t *Qt = FFLAS::fflas_new<size_t>(m);
	size_t *P = FFLAS::fflas_new<size_t>(n);
	for (int i = 0; i < n; i++)
		P[i] = 0;
	for (int j = 0; j < m; j++)
		Qt[j] = 0;
	size_t rank = FFPACK::pReducedRowEchelonForm(GFp, n, m, A, ldA, P, Qt);
	fprintf(stderr, "done in %.1fs. Rank %zd\n", spasm_wtime() - start, rank);
	start = spasm_wtime();
	// FFPACK::getEchelonForm (GFp, FFLAS::FflasUpper, FFLAS::FflasUnit, n, m, rank, Qt, A, ldA, FFPACK::FfpackTileRecursive);
	// FFPACK::getReducedEchelonForm(GFp, FFLAS::FflasUpper,  n, m, rank, Qt, A, ldA, FFPACK::FfpackTileRecursive);
	/* Qt is in LAPACK representation; convert */
	FFPACK::LAPACKPerm2MathPerm (qinv, Qt, m);
	FFLAS::fflas_delete(P);
	FFLAS::fflas_delete(Qt);
	return rank;
}

/* force instatiations */

int spasm_ffpack_rref_double(i64 prime, int n, int m, double *A, int ldA, size_t *qinv)
{
	return spasm_ffpack_rref(prime, n, m, A, ldA, qinv);
}

int spasm_ffpack_rref_float(i64 prime, int n, int m, float *A, int ldA, size_t *qinv)
{
	return spasm_ffpack_rref(prime, n, m, A, ldA, qinv);
}

int spasm_ffpack_rref_i64(i64 prime, int n, int m, i64 *A, int ldA, size_t *qinv)
{
	return spasm_ffpack_rref(prime, n, m, A, ldA, qinv);
}
