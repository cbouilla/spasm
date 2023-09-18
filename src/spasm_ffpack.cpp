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

void spasm_ffpack_setzero(int prime, int n, int m, FFPACK_carrier *A, int ldA)
{
	Givaro::ModularBalanced<FFPACK_carrier> GFp(prime);
	FFLAS::fzero(GFp, n, m, A, ldA);
}

/*
 * Q of size m (#cols). On output, Q[0:rank] is singificant.
 * It is implicit that U[i, 0:Q[i]] == 0 and U[i, 0:Q[i]] == 1.
 * M[i, Q[i]:m] contain the actual data.
 */ 
int spasm_ffpack_echelonize(int prime, int n, int m, FFPACK_carrier *A, int ldA, size_t *qinv)
{
	double start = spasm_wtime();
	fprintf(stderr, "[ffpack/echelonize] Matrix of dimension %d x %d mod %d... ", n, m, prime);
	fflush(stderr);
	Givaro::ModularBalanced<FFPACK_carrier> GFp(prime);
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
