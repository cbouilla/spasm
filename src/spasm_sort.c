#include <assert.h>
#include "spasm.h"

#define INSERT_SORT_THRESHOLD 42 // TODO : tune this value


// sort up to index right, excluded
static void insertion_sort(const spasm *A, int *p, const int left, const int right) {
  int i, j, u, v;

  //  if (left <= 0) {
    for(i = left + 1; i < right; i++) {
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
  } else {
    TODO, possible optimization : if ( left>0 ), then we know that
    on the left of the current subfile, there is an element smaller
    than all the elements of the subfile (because this element was a pivot).
    Therefore, we don't have to check explicitly that we attained the left
    boundary of the array...
    */
}


// standard median-of-three pivoting strategy. Returns the pivot index
static int choose_pivot(const spasm *A, int *p, const int left, const int right) {
  int mid = (left+right)/2;

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


// returns final position of pivot
static int pivoting(const spasm *A, int *p, const int initial_left, const int initial_right, const int pivotIndex) {
  int pivotValue, left, right;

  spasm_swap(p, pivotIndex, initial_right - 1);
  pivotValue = spasm_row_weight(A, p[initial_right - 1]);

  right = initial_right - 2;
  left = initial_left;

  while(left < right) {
    while(spasm_row_weight(A, p[left]) < pivotValue) {
      left++;
    }
    while(spasm_row_weight(A, p[right]) > pivotValue) {
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

static void spasm_quicksort(const spasm *A, int *p, const int left, const int right) {
  int pivotIndex, new_pivotIndex;

  if (right-left > INSERT_SORT_THRESHOLD) {

    pivotIndex = choose_pivot(A, p, left, right);
    new_pivotIndex = pivoting(A, p, left, right, pivotIndex);

    spasm_quicksort(A, p, left, new_pivotIndex);
    spasm_quicksort(A, p, new_pivotIndex + 1, right);
  } else {
    insertion_sort(A, p, left, right);
  }
}


int * spasm_row_sort (const spasm *A) {
  int *p;
  int i, n;

  n = A->n;
  p = spasm_malloc(n * sizeof(int));
  for(i = 0; i < n; i++) {
    p[i] = i;
  }
  spasm_quicksort(A, p, 0, n);
  return p;
}

int * spasm_cheap_pivots(const spasm *A) {
  int n, m, i, j, k, l, px;
  int *q, *p, *Ap, *Aj;

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;

  q = spasm_malloc(m * sizeof(int));
  p = spasm_malloc(n * sizeof(int));

  /* --- Cheap pivot selection ----------------------------------- */
  for(j = 0; j < m; j++) {
    q[j] = -1;
  }
  for(i = 0; i < n; i++) {
    j = -1;
    /* find leftmost entry */
    for(px = Ap[i]; px < Ap[i + 1]; px++) {
      if (j == -1 || Aj[px] > j) {
	j = Aj[px];
      }
    }

    if (j == -1) {
      continue; // empty row
    }
    if (q[j] == -1 || spasm_row_weight(A, i) < spasm_row_weight(A, q[j])) {
      q[j] = i;
    }
  }
  /* --- count cheap pivots ------------ */
  k = 0;
  for(j = 0; j < m; j++) {
    if (q[j] != -1) {
      k++;
    }
  }

  fprintf(stderr, "[LU] found %d cheap pivots\n", k);

  /* --- build corresponding row permutation ---------------------- */
  l = 0;
  for(i = 0; i < n; i++) {
    j = Aj[ Ap[i] ];
    if (q[j] == i) {
      p[l] = i;
      l++;
    } else {
      p[k] = i;
      k++;
    }
  }

  free(q);
  return p;
}
