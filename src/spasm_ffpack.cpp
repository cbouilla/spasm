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
 * Q of size m (#cols). On output, Q[0:rank] is significant.
 * It is implicit that U[i, 0:Q[i]] == 0 and U[i, 0:Q[i]] == 1.
 * M[i, Q[i]:m] contain the actual data.
 */
template<typename T>
static int ffpack_rref(i64 prime, int n, int m, T *A, int ldA, size_t *qinv)
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

/*
 * Q of size m (#cols). On output, Q[0:rank] is significant.
 * It is implicit that U[i, 0:Q[i]] == 0 and U[i, 0:Q[i]] == 1.
 * M[i, Q[i]:m] contain the actual data.
 */
template<typename T>
static int ffpack_LU(i64 prime, int n, int m, T *A, int ldA, size_t *p, size_t *qinv)
{
	double start = spasm_wtime();
	fprintf(stderr, "[ffpack/LU] Matrix of dimension %d x %d mod %" PRId64"... ", n, m, prime);
	fflush(stderr);
	Givaro::ModularBalanced<T> GFp(prime);
	size_t *Qt = FFLAS::fflas_new<size_t>(m);
	size_t *P = FFLAS::fflas_new<size_t>(n);
	for (int i = 0; i < n; i++)
		P[i] = 0;
	for (int j = 0; j < m; j++)
		Qt[j] = 0;
	size_t rank = FFPACK::pPLUQ(GFp, FFLAS::FflasUnit, n, m, A, ldA, P, Qt);
	fprintf(stderr, "done in %.1fs. Rank %zd\n", spasm_wtime() - start, rank);
	start = spasm_wtime();
	// FFPACK::getEchelonForm (GFp, FFLAS::FflasUpper, FFLAS::FflasUnit, n, m, rank, Qt, A, ldA, FFPACK::FfpackTileRecursive);
	// FFPACK::getReducedEchelonForm(GFp, FFLAS::FflasUpper,  n, m, rank, Qt, A, ldA, FFPACK::FfpackTileRecursive);
	/* Qt is in LAPACK representation; convert */
	FFPACK::LAPACKPerm2MathPerm (p, P, n);
	FFPACK::LAPACKPerm2MathPerm (qinv, Qt, m);
	FFLAS::fflas_delete(P);
	FFLAS::fflas_delete(Qt);
	return rank;
}


int spasm_ffpack_rref(i64 prime, int n, int m, void *A, int ldA, spasm_datatype datatype, size_t *qinv)
{
	switch (datatype) {
	case SPASM_DOUBLE: return ffpack_rref<double>(prime, n, m, (double *) A, ldA, qinv);
	case SPASM_FLOAT: return ffpack_rref<float>(prime, n, m, (float *) A, ldA, qinv);
	case SPASM_I64: return ffpack_rref<i64>(prime, n, m, (i64 *) A, ldA, qinv);
	}
	assert(false);
}

int spasm_ffpack_LU(i64 prime, int n, int m, void *A, int ldA, spasm_datatype datatype, size_t *p, size_t *qinv)
{
	switch (datatype) {
	case SPASM_DOUBLE: return ffpack_LU<double>(prime, n, m, (double *) A, ldA, p, qinv);
	case SPASM_FLOAT: return ffpack_LU<float>(prime, n, m, (float *) A, ldA, p, qinv);
	case SPASM_I64: return ffpack_LU<i64>(prime, n, m, (i64 *) A, ldA, p, qinv);
	}
	assert(false);
}

/* code to make the "runtime type selection" mechanisme work */

spasm_ZZp spasm_datatype_read(const void *A, size_t i, spasm_datatype datatype)
{
	switch (datatype) {	
	case SPASM_DOUBLE: return ((double *) A)[i];
	case SPASM_FLOAT: return ((float *) A)[i];
	case SPASM_I64: return ((i64 *) A)[i];
	}	
	assert(false);
}

void spasm_datatype_write(void *A, size_t i, spasm_datatype datatype, spasm_ZZp value)
{
	switch (datatype) {	
	case SPASM_DOUBLE: ((double *) A)[i] = value; return;
	case SPASM_FLOAT: ((float *) A)[i] = value; return;
	case SPASM_I64: ((i64 *) A)[i] = value; return;
	}	
	assert(false);
}

size_t spasm_datatype_size(spasm_datatype datatype)
{
	switch (datatype) {	
	case SPASM_DOUBLE: return sizeof(double);
	case SPASM_FLOAT: return sizeof(float);
	case SPASM_I64: return sizeof(i64);
	}	
	assert(false);	
}

spasm_datatype spasm_datatype_choose(i64 prime)
{
	// return SPASM_DOUBLE;
	if (prime <= 8191)
		return SPASM_FLOAT;
	else if (prime <= 189812531)
		return SPASM_DOUBLE;
	else
		return SPASM_I64;
}

const char * spasm_datatype_name(spasm_datatype datatype)
{
	switch (datatype) {	
	case SPASM_DOUBLE: return "double";
	case SPASM_FLOAT: return "float";
	case SPASM_I64: return "i64";
	}	
	assert(false);	
}