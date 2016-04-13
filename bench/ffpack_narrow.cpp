/* Copyright (c) FFLAS-FFPACK
* Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
* ========LICENCE========
* This file is part of the library FFLAS-FFPACK.
*
* FFLAS-FFPACK is free software: you can redistribute it and/or modify
* it under the terms of the  GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
* ========LICENCE========
*/
#include <fflas-ffpack/fflas-ffpack-config.h>
#include <givaro/modular-balanced.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/utils/timer.h>
#include <fflas-ffpack/utils/Matio.h>
#include <fflas-ffpack/utils/args-parser.h>
#include <fflas-ffpack/ffpack/ffpack.h>


#include <iostream>

extern "C" {
#include "spasm.h"
}

using namespace FFLAS;

int main(int argc, char** argv) {

  int prime, n, m, n_cheap, schur_m, schur_n, rank;
  int *Aj, *Ap, *p, *q;
  spasm_GFp *Ax;
  spasm_triplet * T;
  spasm *A;
  typedef Givaro::Modular<int> Ring;
  Ring::Element * S;

  /* charge la matrice depuis l'entrée standard */
  prime = 42013;
  T = spasm_load_sms(stdin, prime);
  A = spasm_compress(T);
  spasm_triplet_free(T);
  n = A->n;
  m = A->m;
  Aj = A->j;
  Ap = A->p;
  Ax = A->x;

  /* find free pivots. */
  p = spasm_cheap_pivots(A, &n_cheap);
  printf("Found %d free pivots\n", n_cheap);

  q = (int *) spasm_malloc(m * sizeof(int));
  schur_m = m - n_cheap;
  schur_n = schur_m + 10;

  /* mark non-pivotal columns */
  for(int j=0; j<m; j++) {
    q[j] = 0;
  }
  for(int i = 0; i < n_cheap; i++) {
    int inew = p[i];
    int j = Aj[ Ap[inew] ]; /* the pivot is the first entry of each pivotal row */
    assert(j < m);
    q[j] = -1;
  }
  int k = 0;
  for(int j=0; j<m; j++) {
    if (q[j] == 0) {
      q[j] = k++;
    }
  }
  assert(k == schur_m);

  /* allocate dense matrix to hold the schur complement */
  Ring F(prime);
  S = fflas_new(F, schur_n, schur_m);

  int progress = 0;

  spasm_GFp *y = (spasm_GFp *) spasm_malloc(m * sizeof(spasm_GFp));

  for(int k=0; k<schur_n; k++) {

    /* compute a random linear combination of the non-pivotal rows */
    for(int j = 0; j < m; j++) {
      y[j] = 0;
    }
    for(int i = n_cheap; i < n; i++) {
      int inew = p[i];
      spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], rand() % prime, y, prime);
    }

    /* eliminate everything in y */
    for(int i = 0; i < n_cheap; i++) {
      int inew = p[i];
      int j = Aj[ Ap[inew] ];
      const spasm_GFp diagonal_entry = Ax[ Ap[inew] ];
      assert (diagonal_entry != 0);
      if (y[j] == 0) {
        continue;
      }
      const spasm_GFp d = (y[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
      spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], prime - d, y, prime);
    }

    /* copy y into A[k,:] */
    for(int j=0; j<m; j++) {    /* ceci est améliorable (pas besoin de scanner tout le vecteur, juste les cols non-pivot) */
      if (q[j] >= 0) {
        F.init(*(S + k*schur_m + q[j]), y[j]);
      }
    }
   
    progress++;
    fprintf(stderr, "\rBuilding S: %d / %d", progress, schur_n);
    fflush(stderr);
  }

  free(y);

  free(p);
  free(q);
  spasm_csr_free(A);
  
  printf("\nDense Rank...\n");
  rank = FFPACK::Rank(F, schur_n, schur_m, S, schur_m);
  printf("Rank of input matrix: %d (%d + %d)\n", rank + n_cheap, n_cheap, rank);
  
  fflas_delete(S);
  return 0;
}

