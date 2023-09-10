#include <assert.h>
#include <stdlib.h>

#include "spasm.h"

/*
 * General convention: in U, the pivot is the first entry of the row.
 */

/* register a pivot in (i, j) ; return 1 iff it is new in both row i or col j */
static int register_pivot(int i, int j, int *p, int *qinv)
{
	int r = 1;
	int pi = p[i];
	int qinvj = qinv[j];
	if (pi != -1) {
		qinv[pi] = -1;
		r = 0;
	}
	if (qinvj != -1) {
		p[qinvj] = -1;
		r = 0;
	}
	p[i] = j;
	qinv[j] = i;
	return r;
}

/** Faugère-Lachartre pivot search.
 *
 * The leftmost entry of each row is a candidate pivot. Select the sparsest row
 * with a leftmost entry on the given column.
 *
 * update p/qinv and returns the number of pivots found. 
 */
static int spasm_find_FL_pivots(const spasm *A, int *p, int *qinv)
{
	int n = A->n;
	int m = A->m;
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	double start = spasm_wtime();
	int npiv = 0;

	for (int i = 0; i < n; i++) {
		int j = m + 1;         /* locate leftmost entry */
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++)
			if (Aj[px] < j)
				j = Aj[px];
		if (j == m + 1)            /* Skip empty rows */
			continue;
		/* check if it is a sparser pivot */
		if (qinv[j] == -1 || spasm_row_weight(A, i) < spasm_row_weight(A, qinv[j]))
			npiv += register_pivot(i, j, p, qinv);
	}
	fprintf(stderr, "[pivots] Faugère-Lachartre: %d pivots found [%.1fs]\n", npiv, spasm_wtime() - start);
	return npiv;
}


/*
 * Leftovers from FL. Column not occuring on previously selected pivot row
 * can be made pivotal, as this will not create alternating cycles.
 * 
 * w[j] = 1 <===> column j does not appear in a pivotal row
 * 
 */
static int spasm_find_FL_column_pivots(const spasm *A, int *pinv, int *qinv)
{
	int n = A->n;
	int m = A->m;
	i64 *Ap = A->p;
	int *Aj = A->j;
	int npiv = 0;
	int *w = spasm_malloc(m * sizeof(int));
	for (int j = 0; j < m; j++)
		w[j] = 1;
	double start = spasm_wtime();

	/* mark columns on pivotal rows as obstructed */
	for (int i = 0; i < n; i++) {
		if (pinv[i] < 0)
			continue;
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			w[j] = 0;
		}
	}

	/* find new pivots */
	for (int i = 0; i < n; i++) {
		if (pinv[i] >= 0)
			continue;

		/* does A[i,:] have an entry on an unobstructed column? */
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			if (w[j] == 0)
				continue;	/* this column is closed, skip this entry */
			if (qinv[j] >= 0)
				continue;       /* column j already pivotal */
			/* TODO: displace previous pivot on column j if this one is better */
			npiv += register_pivot(i, j, pinv, qinv);
			/* mark the columns occuring on this row as unavailable */
			for (i64 px = Ap[i]; px < Ap[i + 1]; px++) 
				w[Aj[px]] = 0;
			break; /* move on to the next row */
		}
	}
	free(w);
	fprintf(stderr, "[pivots] ``Faugère-Lachartre on columns'': %d pivots found [%.1fs]\n", 
		npiv, spasm_wtime() - start);
	return npiv;
}


/*
 * This implements the greedy parallel algorithm described in
 * https://doi.org/10.1145/3115936.3115944
 */
static inline void BFS_enqueue(char *w, int *queue, int *surviving, int *tail, int j)
{
	queue[(*tail)++] = j;
	*surviving -= w[j];
	w[j] = -1;
}

static inline void BFS_enqueue_row(char *w, int *queue, int *surviving, int *tail, const i64 *Ap, const int *Aj, int i) 
{
	for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
		/* this is the critical section */
		int j = Aj[px];
		if (w[j] >= 0)
			BFS_enqueue(w, queue, surviving, tail, j);
	}
}

static int spasm_find_cycle_free_pivots(const spasm *A, int *pinv, int *qinv)
{
	int n = A->n;
	int m = A->m;
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	int v = spasm_max(1, spasm_min(1000, n / 100));
	int processed = 0;
	int npiv = 0;
	double start = spasm_wtime();
	int *journal = spasm_malloc(n * sizeof(*journal));

	/*
	 * This uses "transactions". Rows with new pivots are appended to the journal.
	 * npiv is the number of such rows. If npiv has not changed since the beginning of
	 * a transaction, then it can be commited right away.
	 */
	#pragma omp parallel
	{
		char *w = spasm_malloc(m * sizeof(*w));
		int *queue = spasm_malloc(m * sizeof(*queue));

		/* workspace initialization */
		int tid = spasm_get_thread_num();
		for(int j = 0; j < m; j++)
			w[j] = 0;

		#pragma omp for schedule(dynamic, 1000)
		for (int i = 0; i < n; i++) {
			/*
			 * for each non-pivotal row, computes the columns reachable from its entries by alternating paths.
			 * Unreachable entries on the row can be chosen as pivots. 
			 * The w[] array is used for marking during the graph traversal. 
			 * Before the search: 
			 *   w[j] == 1 for each non-pivotal entry j on the row 
			 *   w[j] == 0 otherwise 
			 * After the search: 
			 *   w[j] ==  1  for each unreachable non-pivotal entry j on the row (candidate pivot) 
			 *   w[j] == -1  column j is reachable by an alternating path,
			 *                 or is pivotal (has entered the queue at some point) 
			 *   w[j] ==  0  column j was absent and is unreachable
			 */
			if ((tid == 0) && (i % v) == 0) {
				fprintf(stderr, "\r[pivots] %d / %d --- found %d new", processed, n, npiv);
				fflush(stderr);
			}
			if (pinv[i] >= 0)
				continue;   /* row is already pivotal */

			#pragma omp atomic update
			processed++;

			/* we will start reading qinv: begin the transaction by reading npiv */
			int npiv_local;
			#pragma omp atomic read
			npiv_local = npiv;

			/* scatters columns of A[i] into w, enqueue pivotal entries */
			int head = 0;
			int tail = 0;
			int surviving = 0;
			for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
				int j = Aj[px];
				if (qinv[j] < 0) {
					w[j] = 1;
					surviving += 1;
				} else {
					BFS_enqueue(w, queue, &surviving, &tail, j);
				}
			}

			/* BFS. This is where most of the time is spent */
	BFS:
			while (head < tail && surviving > 0) {
				int j = queue[head++];
				int I = qinv[j];
				if (I == -1)
					continue;	/* j is not pivotal: nothing to do */
				BFS_enqueue_row(w, queue, &surviving, &tail, Ap, Aj, I);
			}

			/* scan w for surviving entries */
			if (surviving == 0)
				goto cleanup;   /* no possible pivot */
			
			/* locate survivor in the row */
			int j = -1;
			for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
				j = Aj[px];
				if (w[j] == 1)  /* potential pivot */
					break;
			}
			assert(j != -1);

			/* try to commit the transaction */
			int npiv_target = -1;
			#pragma omp critical
			{
				if (npiv == npiv_local) {
					/* success */
					register_pivot(i, j, pinv, qinv);
					journal[npiv] = j;
					#pragma omp atomic update
					npiv += 1;
				} else {
					/* failure */
					#pragma omp atomic read
					npiv_target = npiv;
				}
			}

			if (npiv_target < 0)
				goto cleanup;  /* commit success */

			/* commit failure: new pivots have been found behind our back. Examine them */
			for (; npiv_local < npiv_target; npiv_local++) {
				int j = journal[npiv_local];
				if (w[j] == 0)	/* the new pivot plays no role here */
					continue;
				if (w[j] == 1) {
					/* a survivors becomes pivotal with this pivot */
					BFS_enqueue(w, queue, &surviving, &tail, j);
				} else {
					/* the new pivot has been hit */
					int i = qinv[j];
					BFS_enqueue_row(w, queue, &surviving, &tail, Ap, Aj, i);
				}
			}
			goto BFS;

	cleanup:
			/* reset w back to zero */
			for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
				int j = Aj[px];
				w[j] = 0;
			}
			for (int px = 0; px < tail; px++) {
				int j = queue[px];
				w[j] = 0;
			}
		}
		free(w);
		free(queue);
	}
	free(journal);
	fprintf(stderr, "\r[pivots] greedy alternating cycle-free search: %d pivots found [%.1fs]\n", 
		npiv, spasm_wtime() - start);
	return npiv;
}

/*
 * Find a permutation of rows/columns that selects pivots without arithmetic operations.
 * Return the number of pivots found. 
 * qinv[j] == i if (i, j) is a pivot or -1 if there is no pivot on column j.
 * pinv[i] == j if (i, j) is a pivot or -1 if there is no pivot on row i.
 *
 * p : row permutations. Pivotal rows are first, in topological order 
 * Both p, pinv and qinv must be preallocated
 */
static int spasm_pivots_find(const spasm *A, int *pinv, int *qinv, struct echelonize_opts *opts)
{
	int n = A->n;
	int m = A->m;
	for (int j = 0; j < m; j++)
		qinv[j] = -1;
	for (int i = 0; i < n; i++)
		pinv[i] = -1;
	int npiv = spasm_find_FL_pivots(A, pinv, qinv);
	npiv += spasm_find_FL_column_pivots(A, pinv, qinv);
	if (opts->enable_greedy_pivot_search)
		npiv += spasm_find_cycle_free_pivots(A, pinv, qinv);
	fprintf(stderr, "\r[pivots] %d pivots found\n", npiv);
	return npiv;
}

/*
 * build row permutation. Pivotal rows go first in topological order,
 * then non-pivotal rows
 */
static void spasm_pivots_reorder(const spasm *A, const int *pinv, const int *qinv, int *p)
{
	int n = A->n;
	int m = A->m;
	int k = 0;
	/* topological sort */
	int *xj = spasm_malloc(m * sizeof(*xj));
	int *marks = spasm_malloc(m * sizeof(*marks));
	for (int j = 0; j < m; j++)
		marks[j] = 0;
	int top = m;
	for (int j = 0; j < m; j++)
		if (qinv[j] != -1 && !marks[j])
			top = spasm_dfs(j, A, top, xj, p, marks, qinv);  /* use p as "pstack" */
	for (int px = top; px < m; px++) {
		int j = xj[px];
		int i = qinv[j];
		if (i != -1) {
			assert(pinv[i] == j);
			p[k] = i;
			k += 1;
		}
	}
	for (int i = 0; i < n; i++)
		if (pinv[i] == -1) {
			p[k] = i;
			k += 1;
		}
	assert(k == n);
	free(xj);
	free(marks);
}

/*
 * Identify stuctural pivots in A, and copy the relevant rows to U
 * Update Uqinv and p (pivotal rows of A first)
 * return the number of pivots found
 */
int spasm_pivots_extract_structural(const spasm *A, spasm *U, int *Uqinv, int *p, struct echelonize_opts *opts)
{
	int n = A->n;
	int m = A->m;
	int *qinv = spasm_malloc(m * sizeof(*qinv));  /* for pivot search */
	int *pinv = spasm_malloc(n * sizeof(*pinv));     /* for pivot search */

	/* find structural pivots in A */
	int npiv = spasm_pivots_find(A, pinv, qinv, opts);

	/* reorder pivots to make U upper-triangular (up to a column permutation) */
	spasm_pivots_reorder(A, pinv, qinv, p);

	/* compute total pivot nnz and reallocate U if necessary */
	i64 pivot_nnz = 0;
	for (int k = 0; k < npiv; k++) {
		int i = p[k];
		pivot_nnz += spasm_row_weight(A, i);
	}
	if (spasm_nnz(U) + pivot_nnz > U->nzmax)
		spasm_csr_realloc(U, spasm_nnz(U) + pivot_nnz);

	/* copy pivotal rows to U and make them unitary; update Uqinv */
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_GFp *Ax = A->x;
	i64 *Up = U->p;
	int *Uj = U->j;
	spasm_GFp *Ux = U->x;
	spasm_GFp prime = A->prime;
	i64 unz = spasm_nnz(U);      /* #entries in U */
	for (int k = 0; k < npiv; k++) {
		int i = p[k];
		int j = pinv[i];
		assert(j >= 0);
		assert(qinv[j] == i);
		Uqinv[j] = U->n;          /* register pivot in U */
		/* locate pivot in row */ 
		spasm_GFp pivot = 0;
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			if (Aj[px] == j) {
				pivot = Ax[px];
				break;
			}
		}
		/* make pivot unitary and add it first */
		spasm_GFp alpha = spasm_GFp_inverse(pivot, prime);
		Uj[unz] = j;
		Ux[unz] = 1;
		unz += 1;
		/* add the rest of the row */
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			if (j == Aj[px])
				continue;    /* skip pivot, already there */
			Uj[unz] = Aj[px];
			Ux[unz] = (alpha * Ax[px]) % prime;
			unz += 1;
		}
		U->n += 1;
		Up[U->n] = unz;
	}
	assert(unz <= U->nzmax);
	free(pinv);
	free(qinv);
	return npiv;
}

#if 0
/*
 * returns a permuted version of A where pivots are pushed to the top-left
 * and form an upper-triangular principal submatrix. qinv is modified.
 */
spasm *spasm_permute_pivots(const spasm *A, const int *p, int *qinv, int npiv)
{
	int m = A->m;
	const i64 *Ap = A->p;
	const int *Aj = A->j;

	/* pivotal columns first */
	int k = 0;
	for (int i = 0; i < npiv; i++) {
		/* the pivot is the first entry of each row */
		int inew = p[i];
		int j = Aj[Ap[inew]];
		qinv[j] = k;
		k += 1;
	}

	/* put remaining non-pivotal columns afterwards, in any order */
	for (int j = 0; j < m; j++)
		if (qinv[j] == -1) {
			qinv[j] = k;
			k += 1;
		}
	return spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
}
#endif