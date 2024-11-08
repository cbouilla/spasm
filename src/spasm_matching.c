#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

void spasm_augment_matching(int head, const int *istack, const int *jstack, int *p, int *qinv)
{
	for (int px = head; px >= 0; px--) {
		int i = istack[px];
		int j = jstack[px];
		qinv[j] = i;
		p[i] = j;
	}
}

/* lookahead: search for unmatched column in row i. If found, it completes the
 * alternating path in istack/jstack, so we augment the matching. */
int spasm_lookahead(const struct spasm_csr *A, int i, int head, int *plookahead, int *istack, int *jstack, int *p, int *qinv)
{
	const i64 *Ap = A->p;
	const int *Aj = A->j;	
	for (i64 px = Ap[i] + plookahead[i]; px < Ap[i + 1]; px++) {
		int j = Aj[px];
		if (qinv[j] < 0) {
			plookahead[i] = px - Ap[i];
			jstack[head] = j;
			spasm_augment_matching(head, istack, jstack, p, qinv);
			return 1;
		}	
	}
	/* all column on A[i,:] are matched. start the DFS */
	plookahead[i] = Ap[i + 1] - Ap[i];
	return 0;
}

/**
 * Search an augmenting path w.r.t. the matching, starting from row istart 
 * (i.e. a path from an unmatched row to an unmatched column).
 *
 * This does a DFS starting from row k, and memorizes the path in (row_stack / 
 * col_stack). When looking for an unmatched column reachable from a row, the 
 * adjacent columns are examined first. This "lookeahead" amounts to do one step
 * of BFS inside the DFS.
 *
 * Because the matching increases monotonically (when row/column is matched, it
 * stays matched), is it useless to re-examine matched columns.
 */
int spasm_augmenting_path(const struct spasm_csr *A, int istart, int *istack, int *jstack, int *pstack, int *marks, int *plookahead, int *p, int *qinv)
{
	const i64 *Ap = A->p;
	const int *Aj = A->j;

	/* initialize the DFS */
	int head = 0;
	istack[head] = istart;

	/* stack empty ? */
	while (head >= 0) {
		/* search an unmatched column reachable from (unmatched) row i */
		int i = istack[head];

		if (marks[i] != istart) {
			marks[i] = istart;
			if (spasm_lookahead(A, i, head, plookahead, istack, jstack, p, qinv))
				return 1;
			/* nothing on row i: we have to start the DFS */
			pstack[head] = 0;
		}
		
		/* Depth-first-search of (matched) columns adjacent to row i */
		i64 px;
		for (px = Ap[i] + pstack[head]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			int inew = qinv[j];
			assert(inew != -1);   /* because we did the lookahead, we now that all columns adjacent to row i are matched */
			if (marks[inew] == istart)
				continue;
			/* pause DFS of row i, start DFS of row inew. */
			pstack[head] += 1;
			jstack[head] = j;
			istack[++head] = inew;
			break;
		}
		
		/* row i is done: pop it from stack */
		if (px == Ap[i + 1])
			head--;
	}
	return 0;
}

/** 
 * Computes a maximum matching using the Ford–Fulkerson algorithm.
 *
 * If the matrix is rectangular, it is a big advantage to transpose it so that n << m.
 * 
 * @param qinv[j] = row matched to column j (or -1) 
 *
 * @param p[i] = column matched to row i (or -1)
 * 
 * @return size of the matching
 */
int spasm_maximum_matching(const struct spasm_csr *A, int *p, int *qinv)
{
	int n = A->n;
	int m = A->m;
	int r = spasm_min(n, m); /* the matching cant' be bigger than this */

	/* get workspace */
	int *istack = spasm_malloc(n * sizeof(int));
	int *jstack = spasm_malloc(n * sizeof(int));
	int *pstack = spasm_malloc(n * sizeof(int));
	int *marks  = spasm_malloc(n * sizeof(int));
	int *plookahead = spasm_malloc(n * sizeof(int));
	for (int j = 0; j < m; j++)
		qinv[j] = -1;
	for (int i = 0; i < n; i++) {
		p[i] = -1;
		marks[i] = -1;
		plookahead[i] = 0;
	}
	int k = 0;
	double start = spasm_wtime();
	for (int i = 0; (i < n) && (k < r); i++) {
		if (p[i] < 0)
			k += spasm_augmenting_path(A, i, istack, jstack, pstack, marks, plookahead, p, qinv);
		fprintf(stderr, "\r[matching] %d / %d, size %d", i, n, k);
		fflush(stderr);
	}
	fprintf(stderr, " [%.1f s]\n", spasm_wtime() - start);
	free(istack);
	free(jstack);
	free(pstack);
	free(marks);
	free(plookahead);
	return k;
}

/*
 * given a row-matching of A, returns a row_matching of P*A*Q --- the result
 * of spasm_permute(A, p, q).
 */
int *spasm_permute_row_matching(int n, const int *jmatch, const int *p, const int *qinv) {
	int *jjmatch;
	int i;

	jjmatch = spasm_malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		if (jmatch[p[i]] == -1) {
			jjmatch[i] = -1;
		} else {
			jjmatch[i] = qinv[jmatch[p[i]]];
		}
	}
	return jjmatch;
}

int *spasm_permute_column_matching(int m, const int *imatch, const int *pinv, const int *q) {
	int *iimatch;
	int j;

	iimatch = spasm_malloc(m * sizeof(int));
	for (j = 0; j < m; j++) {
		if (imatch[q[j]] == -1) {
			iimatch[j] = -1;
		} else {
			iimatch[j] = pinv[imatch[q[j]]];
		}
	}
	return iimatch;
}


/*
 * returns (a copy of) the matching match restricted to the submatrix M[a:b,
 * c:d]
 */
int *spasm_submatching(const int *match, int a, int b, int c, int d) {
	int *pmatch;
	int i;

	pmatch = spasm_malloc((b - a) * sizeof(int));
	for (i = a; i < b; i++) {
		if (match[i] == -1) {
			pmatch[i - a] = -1;
		} else {
			pmatch[i - a] = match[i] - c;
			assert(pmatch[i - a] >= 0);
			assert(pmatch[i - a] < d);
		}
	}
	return pmatch;
}