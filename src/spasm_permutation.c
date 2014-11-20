#include <assert.h>
#include "spasm.h"

/*
 * Permutations matrices are represented by vectors. p[k] = i means that
 * P[k,i] = 1
 */

/*
 * x <-- P.b, for dense vectors x and b; p=NULL denotes identity.
 *
 * This means that x[k] <--- b[ p[k] ]
 */
void cs_pvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n) {
    int k;
    assert(x != NULL);
    assert(b != NULL);

    for (k = 0; k < n; k++) {
        x[k] = b[(p != NULL) ? p[k] : k];
    }
}

/* x <--- P^{-1}.b, for dense vectors x and b; p=NULL denotes identity.
 *
 * This means that x[ p[k] ] <--- b[ k ]
 *
 * The function is given p, not p^{-1}.
 */
void cs_ipvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n) {
    int k;
    assert(x != NULL);
    assert(b != NULL);

    for (k = 0; k < n; k++) {
        x[(p != NULL) ? p[k] : k] = b[k];
    }
}

/* compute the inverse permutation */
int *spasm_pinv(int const *p, int n) {
    int k, *pinv;
    /* p = NULL denotes identity */
    if (p == NULL) {
        return (NULL);
    }
    /* allocate result */
    pinv = spasm_malloc(n * sizeof(int));
    /* invert the permutation */
    for (k = 0; k < n; k++) {
        pinv[ p[k] ] = k;
    }
    return pinv;
}


/*
 * C = P.A.Q where p and q are permutations of 0..m-1 and 0..n-1
 * respectively.
 *
 * beware that p is described by its inverse permutation.
 */
spasm *spasm_permute(const spasm * A, const int *pinv, const int *q, int values) {
    int t, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci;
    spasm_GFp *Cx, *Ax;
    spasm *C;

    /* check inputs */
    assert(spasm_is_csc(A));

    m = A->m;
    n = A->n;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;

    /* alloc result */
    C = spasm_spalloc(m, n, Ap[n], A->prime, values && (Ax != NULL), 0);
    Cp = C->p;
    Ci = C->i;
    Cx = C->x;

    for (k = 0; k < n; k++) {
        /* column k of C is column q[k] of A (denoted by j) */
        Cp[k] = nz;
        j = (q != NULL) ? q[k] : k;
        for (t = Ap[j]; t < Ap[j + 1]; t++) {
            /* row i of A is row pinv[i] of C */
            Ci[nz] = (pinv != NULL) ? pinv[Ai[t]] : Ai[t];
            if (Cx != NULL) {
                Cx[nz] = Ax[t];
            }
            nz++;
        }
    }
    /* finalize the last column of C */
    Cp[n] = nz;
    return C;
}
