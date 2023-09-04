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

void spasm_ffpack_setzero(int prime, int n, int m, double *A, int ldA)
{
	Givaro::ModularBalanced<double> GFp(prime);
	FFLAS::fzero(GFp, n, m, A, ldA);
}

/*
 * Q of size m (#cols). On output, Q[0:rank] is singificant.
 * It is implicit that U[i, 0:Q[i]] == 0 and U[i, 0:Q[i]] == 1.
 * M[i, Q[i]:m] contain the actual data.
 */ 
int spasm_ffpack_echelonize(int prime, int n, int m, double *A, int ldA, size_t *Q)
{
	double start = spasm_wtime();
	fprintf(stderr, "[ffpack/echelonize] Matrix of dimension %d x %d mod %d... ", n, m, prime);
	fflush(stderr);
	Givaro::ModularBalanced<double> GFp(prime);
	size_t *P = FFLAS::fflas_new<size_t>(n);
	for (int j = 0; j < n; j++)
		P[j] = 0;
	for (int j = 0; j < m; j++)
		Q[j] = 0;
	size_t rank = FFPACK::RowEchelonForm(GFp, n, m, A, m, P, Q, 0, FFPACK::FfpackSlabRecursive);
	// FFPACK::getEchelonForm (GFp, FFLAS::FflasUpper, FFLAS::FflasUnit, n, m, rank, Q, A, m, FFPACK::FfpackSlabRecursive);
	FFLAS::fflas_delete(P);
	fprintf(stderr, "done in %.1fs. Rank %zd\n", spasm_wtime() - start, rank);
	return rank;
}