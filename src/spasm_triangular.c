#include <assert.h>
#include "spasm.h"

#ifdef SPASM_TIMING
#include "cycleclock.h"

uint64_t reach = 0, scatter = 0;
#endif


int spasm_is_upper_triangular(const spasm *A) {
  int i, p, n, m, *Aj, *Ap;
  spasm_GFp *Ax;

  n = A->n;
  m = A->m;
  if (n > m) {
    return 0;
  }

  Ap = A->p;
  Aj = A->j;
  Ax = A->x;
  for(i = 0; i < n; i++) {
    /* check diagonal */
    if ( Aj[ Ap[i] ] != i ) {
      return 0;
    }
    if ( Ax[ Ap[i] ] == 0 ) {
      return 0;
    }

    /* check other entries */
    for (p = Ap[i] + 1; p < Ap[i + 1]; p++) {
      if ( Aj[p] < i ) {
	return 0;
      }
    }
  }
  return 1;
}


int spasm_is_lower_triangular(const spasm *A) {
  int i, n, m, p, *Aj, *Ap;
  spasm_GFp *Ax;

  n = A->n;
  m = A->m;
  if (n < m) {
    return 0;
  }

  Ap = A->p;
  Aj = A->j;
  Ax = A->x;
  for(i = 0; i < m; i++) {

    /* check diagonal */
    if (Aj[ Ap[i + 1] - 1 ] != i ) {
      return 0;
    }
    if ( Ax[ Ap[i + 1] - 1 ] == 0 ) {
      return 0;
    }

    /* check other entries */
    for (p = Ap[i]; p < Ap[i + 1] - 1; p++) {
      if ( Aj[p] > i ) {
	return 0;
      }
    }
  }
  return 1;
}


/*
 * Solving triangular systems, dense RHS
 */


/* dense backwards substitution solver. Solve x . L = b where x and b are dense.
 *
 * b is undefined on output
 *
 * L is assumed to be lower-triangular, with non-zero diagonal.
 *
 * The diagonal entry is the **last** of each row.
 * More precisely, L[j,j] is Lx[ Lp[j+1] - 1 ]
 *
 * p[j] == i indicates if the "diagonal" entry on column j is on row i
 *
 */
void spasm_dense_back_solve(const spasm * L, spasm_GFp *b, spasm_GFp * x, const int *p) {
  int i, j, n, m, *Lp, *Lj, prime;
    spasm_GFp *Lx;

    /* check inputs */
    assert(b != NULL);
    assert(x != NULL);
    assert(L != NULL);

    n = L->n;
    m = L->m;
    Lp = L->p;
    Lj = L->j;
    Lx = L->x;
    prime = L->prime;

    for(i = 0; i < n; i++) {
      x[i] = 0;
    }

    for (j = m - 1; j >= 0; j--) {
      i = (p != SPASM_IDENTITY_PERMUTATION) ? p[j] : j;

      /* pivot on the j-th column is on the i-th row */
      const spasm_GFp diagonal_entry = Lx[ Lp[i + 1] - 1 ];
      assert( diagonal_entry != 0 );

      // axpy - inplace
      x[i] = (b[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      spasm_scatter(Lj, Lx, Lp[i], Lp[i + 1] - 1, prime - x[i], b, prime);
    }
}

/* dense forwards substitution solver. Solve x . U = b where x and b are dense.
 *
 * b is undefined on output
 *
 * U is upper-triangular
 *
 * Assumption : the diagonal entry is always present, is always != 0.
 *
 * The diagonal entry is the first one of each row.
 * More precisely, U[i,i] is Ux[ Up[i] ]
 *
 * if q != SPASM_IDENTITY_PERMUTATION, then q[i] indicates the
 * column on which the i-th row pivot is.
 *
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_dense_forward_solve(const spasm * U, spasm_GFp *b, spasm_GFp * x, const int *q) {
  int i, j, n, m, *Up, *Uj, prime;
  spasm_GFp *Ux;

    /* check inputs */
    assert(b != NULL);
    assert(x != NULL);
    assert(U != NULL);

    n = U->n;
    m = U->m;
    assert(n <= m);

    Up = U->p;
    Uj = U->j;
    Ux = U->x;
    prime = U->prime;

    for(i = 0; i < n; i++) {
      x[i] = 0;
    }

    for (i = 0; i < n; i++) {
      j = (q != SPASM_IDENTITY_PERMUTATION) ? q[i] : i;
      if (b[j] != 0) {
      	/* check diagonal entry */
      	const  spasm_GFp diagonal_entry = Ux[ Up[i] ];
      	assert( diagonal_entry != 0 );

      	// axpy - inplace
      	x[i] = (b[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
      	spasm_scatter(Uj, Ux, Up[i] + 1, Up[i + 1], prime - x[i], b, prime);
      	b[j] = 0;
      }
    }

    for(i = 0; i < m; i++) {
      if (b[i] != 0) {
	      return SPASM_NO_SOLUTION;
      }
    }

    return SPASM_SUCCESS;
}



/*************** Triangular solving with sparse RHS
 *
 * solve x * U = B[k], where U is (permuted) upper triangular.
 *
 * x has size m (number of columns of U, paradoxically).
 *
 * when this function returns, the solution is scattered in x, and its pattern
 * is given in xi[top : m].
 *
 * top is the return value.
 *
 */
int spasm_sparse_forward_solve(const spasm *U, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv) {
  int i, I, p, px, top, m, prime, *Up, *Uj, *Bp, *Bj;
  spasm_GFp *Ux, *Bx;

#ifdef SPASM_TIMING
    uint64_t start;
#endif

    assert(U != NULL);
    assert(B != NULL);
    assert(xi != NULL);
    assert(x != NULL);

    m = U->m;
    Up = U->p;
    Uj = U->j;
    Ux = U->x;
    prime = U->prime;

    Bp = B->p;
    Bj = B->j;
    Bx = B->x;

#ifdef SPASM_TIMING
    start = spasm_ticks();
#endif

    /* xi[top : n] = Reach( U, B[k] ) */
    top = spasm_reach(U, B, k, m, xi, pinv);

#ifdef SPASM_TIMING
    reach += spasm_ticks() - start;
#endif

    /* clear x */
    for (p = top; p < m; p++) {
      x[ xi[p] ] = 0;
    }

    /* scatter B[k] into x */
    for (p = Bp[k]; p < Bp[k + 1]; p++) {
        x[ Bj[p] ] = Bx[p];
    }

    /* iterate over the (precomputed) pattern of x (= the solution) */
#ifdef SPASM_TIMING
    start = spasm_ticks();
#endif

    for (px = top; px < m; px++) {
      /* x[i] is nonzero */
      i = xi[px];

      /* i maps to row I of U */
      I = (pinv != NULL) ? (pinv[i]) : i;
      //printf("pinv[%d] = %d\n", i, I);
      if (I < 0) {
	/* row I is empty */
            continue;
      }

      /* get U[i,i] */
      const spasm_GFp diagonal_entry = Ux[ Up[I] ];
      assert( diagonal_entry != 0 );
      // axpy-in-place
      x[i] = (x[i] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      spasm_scatter(Uj, Ux, Up[I] + 1, Up[I + 1], prime - x[i], x, prime);
    }

#ifdef SPASM_TIMING
    scatter += spasm_ticks() - start;
#endif

    return top;
}

/*************** Triangular solving with sparse RHS
 *
 * solve x * U = y, where U is (permuted) upper triangular.
 *
 * x has size m (number of columns of U, paradoxically).
 *
 * when this function returns, the solution is scattered in x, and its pattern
 * is given in xi[top : m].
 *
 * top is the return value.
 *
 */
int spasm_sparse_forward_solve_scat(const spasm *U, int *y, int *yi, int ynz, int *xi, spasm_GFp *x, const int *pinv) {
  int i, I, p, px, top, m, prime, *Up, *Uj;
  spasm_GFp *Ux;

#ifdef SPASM_TIMING
    uint64_t start;
#endif

    assert(U != NULL);
    assert(y != NULL);
    assert(yi != NULL);
    assert(xi != NULL);
    assert(x != NULL);

    m = U->m;
    Up = U->p;
    Uj = U->j;
    Ux = U->x;
    prime = U->prime;

#ifdef SPASM_TIMING
    start = spasm_ticks();
#endif

    /* xi[top : n] = Reach( U, B[k] ) */
    top = spasm_scat_reach(U, yi, 0, ynz, m, xi, pinv);

#ifdef SPASM_TIMING
    reach += spasm_ticks() - start;
#endif

    /* clear x */
    for (p = top; p < m; p++) {
      x[ xi[p] ] = 0;
    }

    /* scatter y into x */
    for (p = 0; p < ynz; p++) {
        x[ yi[p] ] = y[yi[p]];
    }

    /* iterate over the (precomputed) pattern of x (= the solution) */
#ifdef SPASM_TIMING
    start = spasm_ticks();
#endif

    for (px = top; px < m; px++) {
      /* x[i] is nonzero */
      i = xi[px];

      /* i maps to row I of U */
      I = (pinv != NULL) ? (pinv[i]) : i;


      if (I < 0) {
	/* row I is empty */
            continue;
      }

      /* get U[i,i] */
    
      const spasm_GFp diagonal_entry = Ux[ Up[I] ];
      assert( diagonal_entry != 0 );
      // axpy-in-place
      x[i] = (x[i] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      spasm_scatter(Uj, Ux, Up[I] + 1, Up[I + 1], prime - x[i], x, prime);
    }

#ifdef SPASM_TIMING
    scatter += spasm_ticks() - start;
#endif

    return top;
}



/*************** Triangular solving with sparse RHS
 *
 * solve x * L = B[k], where L is (permuted) lower triangular.
 *
 * x has size m (number of columns of L).
 *
 * when this function returns, the solution scattered in x, and its pattern
 * is given in xi[top : n].
 *
 * top is the return value.
 *
 */
int spasm_sparse_backward_solve(const spasm *L, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv, int r_bound) {
  int i, I, p, px, top, n, m, prime, *Lp, *Lj, *Bp, *Bj, tmp;
  spasm_GFp *Lx, *Bx;

#ifdef SPASM_TIMING
    uint64_t start;
#endif

    assert(L != NULL);
    assert(B != NULL);
    assert(xi != NULL);
    assert(x != NULL);
    //assert(pinv == SPASM_IDENTITY_PERMUTATION);

    n = L->n;
    m = L->m;
    Lp = L->p;
    Lj = L->j;
    Lx = L->x;
    prime = L->prime;

    Bp = B->p;
    Bj = B->j;
    Bx = B->x;

#ifdef SPASM_TIMING
    start = spasm_ticks();
#endif
    /* xi[top : m] = Reach( L, B[k] ) */
    top = spasm_reach(L, B, k, n, xi, pinv);

#ifdef SPASM_TIMING
    reach += spasm_ticks() - start;
#endif

    /* clear x */
    for (p = top; p < n; p++) {
      x[ xi[p] ] = 0;
      //      printf("support : %d\n", xi[p]);
    }

    /* scatter B[k] into x */
   for (p = Bp[k]; p < Bp[k + 1]; p++) {
        x[ Bj[p] ] = Bx[p];
    }

    /* iterate over the (precomputed) pattern of x (= the solution) */
#ifdef SPASM_TIMING
    start = spasm_ticks();
#endif

    for (px = top; px < n; px++) {
      /* x[i] is nonzero */
      i = xi[px];

      /* i maps to row I of L */
      I = (pinv == SPASM_IDENTITY_PERMUTATION) ? i : pinv[i];

       if (i >= m) {
	/* column I is part of an implicit identity matrix */
	spasm_scatter(Lj, Lx, Lp[I], Lp[I + 1], prime - x[i], x, prime);
       } else if(i >= r_bound) {
	/* get L[i,i] */
	const spasm_GFp diagonal_entry = Lx[ Lp[I + 1] - 1];
	assert( diagonal_entry != 0 );
	// axpy-in-place
	x[i] = (x[I] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
	spasm_scatter(Lj, Lx, Lp[I], Lp[I + 1] - 1, prime - x[i], x, prime);
      }

       xi[px] = I;
       tmp = x[i];
       x[i] = 0;
       x[ xi[px] ] = tmp;

    }




#ifdef SPASM_TIMING
    scatter += spasm_ticks() - start;
#endif

    return top;
}


/*
 * Solve a system x*L = y with L "almost triangular" 
 * with (implicit) pivots on its diagonal.
 *
 * L : super_spasm "triangular" system.
 * y : right-hand side scattered in a dense vector of size L->n
 * yi : y pattern
 * start : begining of yi
 * x : solution (vector size L->n)
 * xi : x pattern
 *
 * retrun value top : begining of xi
 */
int super_spasm_sparse_solve(super_spasm *L, int *y, int *yi, int end, int *x, int *xi){
  int top, *pinv, n_big, *Lperm, i, n_small, p, i_new, *Lj, prime, *Lp;
  spasm_GFp *Lx;

  /* check inputs */  
  assert(L != NULL);
  assert(L->M != NULL);
  assert(y != NULL);
  assert(yi != NULL);
  assert(x != NULL);
  assert(xi != NULL);

  Lperm = L->p;
  n_big = L->n;
  n_small = L->M->n;
  Lj = L->M->j;
  Lp = L->M->p;
  Lx = L->M->x;
  prime = L->M->prime;

  assert(end > 0); //yi not empty.

  /* get workspace */
  // pinv[i] : corresponding row in L->M.
  pinv = spasm_malloc(n_big * sizeof(int)); 

  /* initialize pinv */
  for(i = 0; i < n_big; i++){
    pinv[i] = -1;
  }
  for(i = 0; i < n_small; i++){
    pinv[Lperm[i]] = i;
  }

  /* find x pattern xi */
  // xi[top : m] = Reach( L, y )
  top = spasm_scat_reach(L->M, yi, 0, end, n_big, xi, pinv);

  /* initialize x */
  for (p = top; p < n_big; p++) {
    // clear x
    x[ xi[p] ] = 0;
  }

  /* scatter y into x */
  for (p = 0; p < end; p++) {
    x[ yi[p] ] = y[yi[p]];
  }

   /* iterate over the (precomputed) pattern of x (= the solution) */
  for(p = top; p < n_big; p++){
    i = xi[p]; // x[i] non zero.
    i_new = pinv[i]; // corresponding row in L->M
    
    if(i_new == -1){
      // implicit identity row
      continue;
    }

    //update x
    spasm_scatter(Lj, Lx, Lp[i_new], Lp[i_new+1], prime - x[i], x, prime);

  }

  /* free workspace */
  free(pinv);

  return top;
}
