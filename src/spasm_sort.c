/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_sort.c */
#include <assert.h>
#include "spasm.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif

#define INSERT_SORT_THRESHOLD 42/* TODO : tune this value */


/* sort up to index right, excluded */
static void insertion_sort(const spasm * A, int *p, const int left, const int right) {
  int i, j, u, v;

  /* if (left <= 0) { */
  for (i = left + 1; i < right; i++) {
    u = p[i];
    v = spasm_row_weight(A, p[i]);
    j = i - 1;
    while (j >= 0 && spasm_row_weight(A, p[j]) > v) {
      p[j + 1] = p[j];
      j--;
    }
    p[j + 1] = u;
  }
  /*
   * } else { TODO, possible optimization : if ( left>0 ), then we know that
   * on the left of the current subfile, there is an element smaller than all
   * the elements of the subfile (because this element was a pivot).
   * Therefore, we don't have to check explicitly that we attained the left
   * boundary of the array...
   */
}


/* standard median-of-three pivoting strategy. Returns the pivot index */
static int choose_pivot(const spasm * A, int *p, const int left, const int right) {
  int mid = (left + right) / 2;

  if (spasm_row_weight(A, p[mid - 1]) > spasm_row_weight(A, p[mid])) {
    spasm_swap(p, mid - 1, mid);
  }
  if (spasm_row_weight(A, p[mid - 1]) > spasm_row_weight(A, p[mid + 1])) {
    spasm_swap(p, mid - 1, mid + 1);
  }
  if (spasm_row_weight(A, p[mid]) > spasm_row_weight(A, p[mid + 1])) {
    spasm_swap(p, mid, mid + 1);
  }
  return mid;
}


/* returns final position of pivot */
static int pivoting(const spasm * A, int *p, const int initial_left, const int initial_right, const int pivotIndex) {
  int pivotValue, left, right;

  spasm_swap(p, pivotIndex, initial_right - 1);
  pivotValue = spasm_row_weight(A, p[initial_right - 1]);

  right = initial_right - 2;
  left = initial_left;

  while (left < right) {
    while (spasm_row_weight(A, p[left]) < pivotValue) {
      left++;
    }
    while (spasm_row_weight(A, p[right]) > pivotValue) {
      right--;
    }

    if (left < right) {
      spasm_swap(p, left, right);
      left++;
    }
  }

  if (spasm_row_weight(A, p[right]) < pivotValue) {
    right++;
  }
  spasm_swap(p, right, initial_right - 1);
  return right;
}

static void spasm_quicksort(const spasm * A, int *p, const int left, const int right) {
  int pivotIndex, new_pivotIndex;

  if (right - left > INSERT_SORT_THRESHOLD) {

    pivotIndex = choose_pivot(A, p, left, right);
    new_pivotIndex = pivoting(A, p, left, right, pivotIndex);

    spasm_quicksort(A, p, left, new_pivotIndex);
    spasm_quicksort(A, p, new_pivotIndex + 1, right);
  } else {
    insertion_sort(A, p, left, right);
  }
}


int *spasm_row_sort(const spasm * A) {
  int *p;
  int i, n;

  n = A->n;
  p = spasm_malloc(n * sizeof(int));
  for (i = 0; i < n; i++) {
    p[i] = i;
  }
  spasm_quicksort(A, p, 0, n);
  return p;
}


int spasm_is_row_pivotal(const spasm *A, const int *qinv, const int i) {
  int *Ap, *Aj;

  Ap = A->p;
  Aj = A->j;
  return (qinv[Aj[Ap[i]]] == i);
}

/* make pivot the first entry of the row */
void spasm_prepare_pivot(const spasm *A, const int i, const int px) {
  int *Ap, *Aj;
  spasm_GFp *Ax;

  Ap = A->p;
  Aj = A->j;
  Ax = A->x;

  spasm_swap(Aj, Ap[i], px);
  if (Ax != NULL)
    spasm_swap(Ax, Ap[i], px);
}


/* Faugère-Lachartre pivot search.
 * The leftmost entry of each row is a candidate pivot. Select the sparsest row
 * with a leftmost entry on the given column. Selected pivots are moved to the
 * front of the row.
 * qinv must be initialized to -1
 * p can be arbitrary.
 * Return the number of pivots found. */
int spasm_find_FL_pivots(const spasm * A, int *p, int *qinv) {
  int n, m, idx_j, *Aj, *Ap, npiv;
  spasm_GFp *Ax;
  double start;

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;
  Ax = A->x;
  start = spasm_wtime();

  for (int i = 0; i < n; i++) {
    int j = -1;
    for (int px = Ap[i]; px < Ap[i + 1]; px++) {
      if (j == -1 || Aj[px] < j) {
	j = Aj[px];
	idx_j = px;
      }
    }
    
    if (j == -1) { /* Skip empty rows */
      continue;
    }

    /* check if it is a sparser pivot */
    if (qinv[j] == -1 || spasm_row_weight(A, i) < spasm_row_weight(A, qinv[j])) {
      qinv[j] = i;
      spasm_prepare_pivot(A, i, idx_j);
    }
  }

  /* build p */
  npiv = 0;
  for (int j = 0; j < m; j++)
    if (qinv[j] != -1)
      p[npiv++] = qinv[j];

  fprintf(stderr, "[pivots] Faugère-Lachartre: %d pivots found [%.1fs]\n", npiv, spasm_wtime() - start);
  return npiv;
}



/* Leftovers from FL.
 * Column not occuring on previously selected pivot row can be made pivotal, as
 * this will not create alternating cycles.
 *
 * w[j] = 1 <===> column j does not appear in a pivotal row
 * 
 */
int spasm_find_FL_column_pivots(const spasm * A, int *p, int *qinv, int npiv_fl) {
  int n, m, *Aj, *Ap, npiv, *w;
  spasm_GFp *Ax;
  double start;

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;
  Ax = A->x;
  npiv = npiv_fl;
  w = spasm_malloc(m * sizeof(int));
  spasm_vector_set(w, 0, m, 1);

  start = spasm_wtime();

  /* mark columns on previous pivot rows as obstructed */
  for (int i = 0; i < npiv; i++) {
    int inew = p[i];
    for (int px = Ap[inew]; px < Ap[inew + 1]; px++)
      w[Aj[px]] = 0;
  }

  /* find new pivots */
  for (int i = 0; i < n; i++) {
    if (spasm_is_row_pivotal(A, qinv, i))
      continue;
   
    /* does A[i,:] have an entry on an unobstructed column? */
    for (int px = Ap[i]; px < Ap[i + 1]; px++) {
      int j = Aj[px];
      if (w[j] == 0) {
	continue; /* this column is closed, skip this entry */
      }

      /* new pivot found! */
      if (qinv[j] == -1) {
	p[npiv++] = i;
	qinv[j] = i;
	spasm_prepare_pivot(A, i, px);
	/* mark the columns occuring on this row as unavailable */     
	for (int px = Ap[i]; px < Ap[i + 1]; px++)
	  w[Aj[px]] = 0;

	break;
      }
    }
  }
  free(w);

  fprintf(stderr, "[pivots] ``Faugère-Lachartre on columns'': %d pivots found [%.1fs]\n", npiv - npiv_fl, spasm_wtime() - start);
  return npiv;
}



int find_survivor(int *Ap, int *Aj, spasm_GFp * Ax, int i, int *w) {
  for (int px = Ap[i]; px < Ap[i + 1]; px++) {
    int j = Aj[px];
    if (w[j] == 1) {            /* potential pivot found */
      /* move potential pivot entry to the front */
      spasm_swap(Aj, Ap[i], px);
      if (Ax != NULL) {
	spasm_swap(Ax, Ap[i], px);
      }
      return j;
    }
  }
  return -1;
}

/*
 * provide already know pivots, and this looks for more. Updates qinv, but
 * DFS must be performed afterwards
 */
void BFS_enqueue(int *w, int *queue, int *surviving, int *tail, int j) {
  queue[(*tail)++] = j;
  *surviving -= w[j];
  w[j] = -1;
}

void BFS_enqueue_row(int *w, int *queue, int *surviving, int *tail, const int *Ap, const int *Aj, int i) {
  for (int px = Ap[i]; px < Ap[i + 1]; px++) {
    /* this is the critical section */
    int j = Aj[px];
    if (w[j] >= 0)
      BFS_enqueue(w, queue, surviving, tail, j);
  }
}

int spasm_find_cycle_free_pivots(const spasm * A, int *p, int *qinv, int npiv_start) {
  int n, m, *Aj, *Ap, processed, v, npiv, retries;
  double start;
  spasm_GFp *Ax;

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;
  Ax = A->x;
  v = spasm_max(1, spasm_min(1000, n / 100));
  processed = 0;
  retries = 0;
  npiv = npiv_start;
  start = spasm_wtime();

#pragma omp parallel
  {
    int *w = spasm_malloc(m * sizeof(int));
    int *queue = spasm_malloc(2 * m * sizeof(int));
    int head, tail, npiv_local, surviving, tid;

    /* workspace initialization */
    tid = 0;
    spasm_vector_set(w, 0, m, 0);

#ifdef USE_OPENMP
    tid = omp_get_thread_num();
    if (tid == 0)
	fprintf(stderr, "[pivots] Greedy pivot search starting on %d threads\n",  omp_get_num_threads());
#endif


#pragma omp for schedule(dynamic, 1000)
    for (int i = 0; i < n; i++) {
      /*
       * for each non-pivotal row, computes the columns reachable from its
       * entries by alternating paths. Unreachable entries on the row can be
       * chosen as pivots. The w[] array is used for marking during the graph
       * traversal. Before the search: w[j] ==  1  for each non-pivotal entry
       * j on the row w[j] ==  0  otherwise After the search: w[j] ==  1  for
       * each unreachable non-pivotal entry j on the row w[j] == -1  column j
       * is reachable by an alternating path, or is pivotal (has entered the
       * queue at some point) w[j] ==  0  column j was absent and is
       * unreachable
       */
      if (tid == 0 && (i % v) == 0) {
	fprintf(stderr, "\r[pivots] %d / %d --- found %d new --- %d retries", processed, n - npiv_start, npiv - npiv_start, retries);
	fflush(stderr);
      }
      if (qinv[Aj[Ap[i]]] == i) /* this row was already pivotal before the
				 * start: skip */
	continue;

#pragma omp atomic update
      processed++;

      /* we start reading qinv: begining of transaction */
#pragma omp atomic read
      npiv_local = npiv;

      /* scatters columns of A[i] into w, enqueue pivotal entries */
      head = 0;
      tail = 0;
      surviving = 0;
      for (int px = Ap[i]; px < Ap[i + 1]; px++) {
	int j = Aj[px];
	if (qinv[j] < 0) {
	  w[j] = 1;
	  surviving++;
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
	  continue;             /* j is not pivotal: nothing to do */
	BFS_enqueue_row(w, queue, &surviving, &tail, Ap, Aj, I);
      }

      /* scan w for surviving entries */
      if (surviving > 0) {
	int j = find_survivor(Ap, Aj, Ax, i, w);
	int npiv_target = -1;

	/* si aucun nouveau pivot n'est arrivé, ajouter */
#pragma omp critical
	{
	  if (npiv == npiv_local) {
	    qinv[j] = i;
	    p[npiv] = i;
#pragma omp atomic update
	    npiv++;
	  } else {
#pragma omp atomic read
	    npiv_target = npiv;
	    retries++;
	  }
	}

	if (npiv_target < 0)
	  goto cleanup;

	/* si on a découvert de nouveaux pivots à traiter... les traiter ! */
	for (; npiv_local < npiv_target; npiv_local++) {
	  int I = p[npiv_local];
	  int j = Aj[Ap[I]];
	  if (w[j] == 0)        /* the new pivot plays no role here */
	    continue;

	  if (w[j] == 1) {      /* a survivors becomes pivotal with this new
				 * pivot */
	    BFS_enqueue(w, queue, &surviving, &tail, j);
	  } else {
	    BFS_enqueue_row(w, queue, &surviving, &tail, Ap, Aj, I);
	  }
	}
	goto BFS;
      }
      /* reset w back to zero */
  cleanup:
      for (int px = Ap[i]; px < Ap[i + 1]; px++)
	w[Aj[px]] = 0;
      for (int px = 0; px < tail; px++)
	w[queue[px]] = 0;
    } /* end for */
    free(w);
    free(queue);
  } /* end of omp parallel */

  fprintf(stderr, "\r[pivots] greedy alternating cycle-free search: %d pivots found [%.1f]\n", npiv - npiv_start, spasm_wtime() - start);
  return npiv;
}

/*
 * return the number of pivots found. @param p : row permutations. Pivotal
 * rows are first. @param qinv : inverse column permutation. q[j] is the row
 * on which the pivot on column j is, or -1 if there is no pivot on column j.
 * both p and qinv must be preallocated
 */
int spasm_find_pivots(const spasm * A, int *p, int *qinv) {
  int n, m, i, j, k, npiv;
  int *Ap, *Aj;

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;

  spasm_vector_set(qinv, 0, m, -1);
  npiv = spasm_find_FL_pivots(A, p, qinv);
  /* npiv = spasm_find_FL_column_pivots(A, p, qinv, npiv); */
  npiv = spasm_find_cycle_free_pivots(A, p, qinv, npiv);

  /* --- build corresponding row permutation ---------------------- */
  int *xj = spasm_malloc(m * sizeof(int));
  int *marks = spasm_malloc(m * sizeof(int));
  int *pstack = spasm_malloc(n * sizeof(int));
  
  /* DFS */
  spasm_vector_set(marks, 0, m, 0);
  int top = m;
  for (j = 0; j < m; j++) {
    if (qinv[j] != -1 && !marks[j]) {
      top = spasm_dfs(j, A, top, xj, pstack, marks, qinv);
    }
  }

  /* reorders the first npiv rows of p */
  k = 0;
  for (j = top; j < m; j++) {
    i = qinv[xj[j]];
    if (i != -1)
      p[k++] = i;
  }
  assert(k == npiv);
  free(xj);
  free(pstack);
  free(marks);

  /* put other (non-empty) rows afterwards */
  for (i = 0; i < n; i++) {
    if (Ap[i] == Ap[i + 1]) {
      continue;
    }
    j = Aj[Ap[i]];
    if (qinv[j] != i) { /* row is non-pivotal */
      assert(k < n);
      p[k] = i;
      k++;
    }
  }

  /* put empty rows last */
  for (i = 0; i < n; i++) {
    if (Ap[i] == Ap[i + 1]) {
      p[k] = i;
      k++;
    }
  }

  fprintf(stderr, "\r[pivots] %d pivots found\n", npiv);
  return npiv;
}

/* returns a permuted version of A where pivots are pushed to the top-left
* and form an upper-triangular principal submatrix */
spasm * spasm_permute_pivots(const spasm *A, int *p, int *qinv, int npiv) {
  int i, j, k, m, *Ap, *Aj;

  m = A->m;
  Ap = A->p;
  Aj = A->j;

  /* pivotal column first, in row-order */
  k = 0;
  for (i = 0; i < npiv; i++) {
    j = Aj[Ap[p[i]]];         /* the pivot is the first entry of each row */
    qinv[j] = k++;
  }

  /* put remaining non-pivotal columns afterwards, in any order */
  for (j = 0; j < m; j++) {
    if (qinv[j] == -1) {
      qinv[j] = k++;
    }
  }

  return spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
}
