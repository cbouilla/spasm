#include "spasm.h"

#define INSERT_SORT_THRESHOLD 42


// sort up to index right, excluded
static void insertion_sort(const spasm *A, int *p, const int left, const int right) {
  int i, j, u, v;

  //  printf("InsertionSort(%d, %d);\n", left, right);

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

    // if ( left>0 ), then we know that on the left of the current
    // subfile, there is an element smaller than all the elements of
    // the subfile (because this element was a pivot). Therefore, we
    // don't have to check explicitly that we attained the left
    // boundary of the array...

    for(i = left + 1; i <= right; i++) {
	v = a[i];
	j = i - 1;
	while (a[j] > v) {
	  a[j + 1] = a[j];
	  j--;
	}
	a[j + 1] = v;
      }
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

  //  printf("QuickSort(%d, %d);\n", left, right);

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
