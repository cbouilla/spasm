/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_sort.c */
#include <assert.h>
#include "spasm.h"

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


int *spasm_cheap_pivots(const spasm * A, int *cheap_ptr) {
  int n, m, i, j, k, I, idx_j, px, n_cheap, pxI, head, tail;
  int *q, *p, *Ap, *Aj, *w, *queue;
  spasm_GFp *Ax;

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;
  Ax = A->x;

  q = spasm_malloc(m * sizeof(int));
  p = spasm_malloc(n * sizeof(int));
#if 0
  w = spasm_malloc(m * sizeof(int));
  queue = spasm_malloc(m * sizeof(int));
#endif

  /* --- Cheap pivot selection ----------------------------------- */
  for (j = 0; j < m; j++) {
    q[j] = -1;
  }
  for (i = 0; i < n; i++) {

    /* find leftmost entry */
    j = -1;
    for (px = Ap[i]; px < Ap[i + 1]; px++) {
      if (j == -1 || Aj[px] < j) {
        j = Aj[px];
        idx_j = px;
      }
    }
    /* Skip empty rows */
    if (j == -1) {
      continue;
    }
    /* make sure leftmost entry is the first of the row */
    spasm_swap(Aj, Ap[i], idx_j);
    spasm_swap(Ax, Ap[i], idx_j);

    /* check if it is a sparser pivot */
    if (q[j] == -1 || spasm_row_weight(A, i) < spasm_row_weight(A, q[j])) {
      q[j] = i;
    }
  }

  /* count the Faugère-Lachartre pivots, and store their rows in p */
  k = 0;
  for (j = 0; j < m; j++) {
    if (q[j] != -1) {
      p[k++] = q[j];
    }
  }
  fprintf(stderr, "[LU] found %d cheap pivots (stage1)\n", k);

#if 0
  /* --- transitive reduction ------------------------------------- */
  fprintf(stderr, "starting transitive reduction...\n");
  spasm *T = spasm_csr_alloc(k, m, A->nzmax, -1, SPASM_IGNORE_VALUES);
  int *Tp = T->p;
  int *Tj = T->j;
  int *qinv_red = spasm_malloc(m * sizeof(int));
  int tnz = 0, anz = 0;

  /*
   * workspace initialization. Marking scheme : -1 = unmarked 0 = marked by
   * original row (=kept in transitive reduction) 1 = marked by other,
   * reachable, rows (=removed in transitive reduction)
   */
  for (j = 0; j < m; j++) {
    w[j] = -1;
    qinv_red[j] = -1;
  }

  for (i = 0; i < k; i++) {
    /* perform reduction of row A[I], store result in T[i] */
    I = p[k - 1 - i];

    /* initialize BFS in T with A[I] */
    Tp[i] = tnz;

    head = 0;
    tail = 0;
    for (px = Ap[I]; px < Ap[I + 1]; px++) {
      j = Aj[px];
      queue[tail] = j;
      tail++;
      w[j] = 0;
    }
    anz += spasm_row_weight(A, I);

    /* start BFS */
    while (head < tail) {
      j = queue[head];
      head++;

      int Ired = qinv_red[j];
      if (Ired == -1) {
        continue;
      }
      /*
       * trick: the first entry is the pivot, we know it is already marked,
       * so we skip it
       */
      for (pxI = Tp[Ired] + 1; pxI < Tp[Ired + 1]; pxI++) {
        j = Tj[pxI];
        if (w[j] < 0) {         /* not marked : mark and add to queue */
          queue[tail] = j;
          tail++;
        }
        w[j] = 1;
      }
    }

    /* scan w for surviving entries, add them to T */
    for (px = Ap[I]; px < Ap[I + 1]; px++) {
      j = Aj[px];
      if (w[j] == 0) {          /* entry has not been "superseded", so we add
                                 * it to T */
        Tj[tnz] = j;
        tnz++;
      }
    }

    /* reset w */
    for (px = 0; px < tail; px++) {
      j = queue[px];
      w[j] = -1;
    }

    qinv_red[Aj[Ap[I]]] = i;
    fprintf(stderr, "\r%d / %d | %d vs %d", i, k, tnz, anz);
    fflush(stderr);
  }
  fprintf(stderr, "\r                                             \r");
  /* finalize the last row of T */
  Tp[k] = tnz;
  n_cheap = k;
  fprintf(stderr, "done. |A| = %d, |T| = %d\n", anz, tnz);
#endif

  /* --- find less-cheap pivots ----------------------------------- */
#if 0
  n_cheap = k;
  int n_cheap1 = k;
  int processed = 0;
  /* workspace initialization */
  for (j = 0; j < m; j++) {
    w[j] = 0;
  }

  for (i = 0; i < n; i++) {
    if (q[Aj[Ap[i]]] == i) {    /* this row is already pivotal: skip */
      continue;
    }
    if (i % (n / 100) == 0) {
      fprintf(stderr, "\rcheap : %d / %d --- found %d new", processed, n - n_cheap1, n_cheap - n_cheap1);
      fflush(stderr);
    }
    processed++;
    //printf("------------------------ %d\n", i);

    /* scatters non-pivotal columns of A[i] into w */
    for (px = Ap[i]; px < Ap[i + 1]; px++) {
      j = Aj[px];
      if (q[j] != -1) {         /* column is pivotal: skip */
        continue;
      }
      w[j] = 1;
    }

    head = 0;
    tail = 0;
    for (px = Ap[i]; px < Ap[i + 1]; px++) {
      j = Aj[px];
      if (q[j] == -1) {         /* column is not pivotal: skip */
        continue;
      }
      if (w[j] < 0) {           /* already marked: skip */
        continue;
      }
      queue[tail] = j;
      tail++;
      w[j] = -1;
      //printf("amorçage en %d\n", j);
      /* BFS */
      while (head < tail) {
        j = queue[head];
        //assert(w[j] < 0);
        head++;

        I = q[j];
        if (I == -1) {
          continue;             /* nothing to do */
        }
        for (pxI = Ap[I]; pxI < Ap[I + 1]; pxI++) {
          j = Aj[pxI];
          if (w[j] < 0) {
            continue;           /* already marked */
          }
          queue[tail] = j;
          tail++;
          w[j] = -1;
          /*
           * ne peut-on éviter d'empiler des trucs qui déclenchent le
           * "continue" ci-dessus ?
           */
        }
      }
    }

    /* scan w for surviving entries */
    k = -1;
    for (px = Ap[i]; px < Ap[i + 1]; px++) {
      j = Aj[px];
      if ((k == -1) && (w[j] == 1)) {
        k = j;
        idx_j = px;
      }
      w[j] = 0;                 /* utile ? */
    }

    /* reset w */
    for (px = 0; px < tail; px++) {
      j = queue[px];
      w[j] = 0;
    }

    if (k != -1) {
      q[k] = i;

      /* make sure leftmost entry is the first of the row */
      spasm_swap(Aj, Ap[i], idx_j);
      spasm_swap(Ax, Ap[i], idx_j);
      n_cheap++;
    }
  }
  fprintf(stderr, "[LU] found %d cheap pivots (stage2)\n", n_cheap);
#endif

  /* --- build corresponding row permutation ---------------------- */

  int *xj = spasm_malloc(m * sizeof(int));
  int *marks = spasm_malloc(m * sizeof(int));
  int *pstack = spasm_malloc(n * sizeof(int));

  for (j = 0; j < m; j++) {
    marks[j] = 0;
  }

  /* DFS */
  int top = m;
  for (j = 0; j < m; j++) {
    if (q[j] != -1 && !marks[j]) {
      top = spasm_dfs(j, A, top, xj, pstack, marks, q);
    }
  }

  /* reorders the first n_cheap rows of p */
  k = 0;
  for (j = top; j < m; j++) {
    if (q[xj[j]] != -1) {
      p[k++] = q[xj[j]];
    }
  }
  assert(k == n_cheap);

  free(xj);
  free(pstack);
  free(marks);

  n_cheap = k;
  *cheap_ptr = n_cheap;

  /* put other (non-empty) rows afterwards */
  for (i = 0; i < n; i++) {
    if (Ap[i] == Ap[i + 1]) {
      continue;
    }
    if (q[Aj[Ap[i]]] != i) {
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

  free(q);

  return p;
}
