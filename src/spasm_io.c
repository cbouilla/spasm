#include <assert.h>
#include "spasm.h"

// set prime = -1 to avoid loading values
spasm_triplet * spasm_load_sms(FILE *f, int prime) {
    int i, j ;   /* use double for integers to avoid int conflicts */
    spasm_GFp x ;
    spasm_triplet *T;
    char type;
    assert(f != NULL);

    if (fscanf(f, "%d %d %c\n", &i, &j, &type) != 3) {
      fprintf(stderr, "[spasm_load_sms] bad SMS file (header)\n");
      exit(0);
    }

    if (prime != -1 && type != 'M') {
      fprintf(stderr, "[spasm_load_sms] only ``Modular'' type supported\n");
      exit(0);
    }

    /* allocate result */
    T = spasm_triplet_alloc (i, j, 1, prime, prime != -1);

    while (fscanf (f, "%d %d %d\n", &i, &j, &x) == 3) {
      if (i == 0 && j == 0 && x == 0) {
	break;
      }
      assert(i != 0);
      assert(j != 0);
      spasm_add_entry(T, i - 1, j - 1, x);
    }

    return T;
}


void spasm_save_csr(FILE *f, const spasm *A) {
  int i, n, m, p;
  int *Aj, *Ap;
  spasm_GFp *Ax;

    assert(f != NULL);
    assert(A != NULL);

    Aj = A->j;
    Ap = A->p;
    Ax = A->x;
    n  = A->n;
    m  = A->m;

    fprintf(f, "%d %d M\n", n, m);
    for(i = 0; i < n; i++) {
      for(p = Ap[i]; p < Ap[i + 1]; p++) {
	fprintf(f, "%d %d %d\n", i + 1, Aj[p] + 1, (Ax != NULL) ? Ax[p] : 1);
      }
    }

    fprintf(f, "0 0 0\n");
}

void spasm_save_triplet(FILE *f, const spasm_triplet *A) {
  int i, nz, n, m;
  int *Ai, *Aj;
  spasm_GFp *Ax;

    assert(f != NULL);
    assert(A != NULL);
    Ai = A->i;
    Aj = A->j;
    Ax = A->x;
    nz = A->nz;
    n  = A->n;
    m  = A->m;

    fprintf(f, "%d %d M\n", n, m);

    for(i = 0; i < nz; i++) {
      fprintf(f, "%d %d %d\n", Ai[i] + 1, Aj[i] + 1, (Ax != NULL) ? Ax[i] : 1);
    }

    fprintf(f, "0 0 0\n");
}

/* Saves a PBM (bitmap) file with one pixel per entry of A */
void spasm_save_pbm(FILE *f, const spasm *A) {
  int i, j, n, m, p;
  int *Aj, *Ap, *x;

  assert(f != NULL);
  assert(A != NULL);

  Aj = A->j;
  Ap = A->p;
  n  = A->n;
  m  = A->m;
  x = spasm_malloc(m * sizeof(int));
  for(j = 0; j < m; j++) {
    x[j] = 0;
  }

  fprintf(f, "P1\n");
  fprintf(f, "%d %d\n", n, m);
  for(i = 0; i < n; i++) {

    // scatters row i to x
    for(p = Ap[i]; p < Ap[i + 1]; p++) {
      x[ Aj[p] ] = 1;
    }

    // print row i
    for(j = 0; j < m; j++) {
      fprintf(f, "%d ", x[j]);
    }

    // reset x
    for(p = Ap[i]; p < Ap[i + 1]; p++) {
      x[ Aj[p] ] = 0;
    }
  }

  free(x);
}

/* Saves a PBM (graymap) of specified dimensions of A */
void spasm_save_pgm(FILE *f, int x, int y, const spasm *A) {
  int i, j, k, n, m, t, p;
  int *Aj, *Ap, *w;
  double max;

  assert(f != NULL);
  assert(A != NULL);

  Aj = A->j;
  Ap = A->p;
  n  = A->n;
  m  = A->m;
  x = spasm_min(x, m);
  y = spasm_min(y, n);

  w = spasm_malloc(x * sizeof(int));
  for(j = 0; j < x; j++) {
    w[j] = 0;
  }

  fprintf(f, "P2\n");
  fprintf(f, "%d %d\n", x, y);
  fprintf(f, "255\n");

  max = (1.0 * m / x) * (1.0 * n / y);
  t = 0;
  i = 0;
  while(i < n) {
    for(k = 0; k < spasm_max(1, n / y) && i < n; k++) {

      // scatters row i to x
      for(p = Ap[i]; p < Ap[i + 1]; p++) {
	w[ (Aj[p] * x) / m ]++;
      }
      i++;
    }

    // print row
    for(j = 0; j < x; j++) {
      double intensity = 1.0 - w[j] / max;
      assert( 0 <= intensity && intensity <= 1.0 );
      fprintf(f, "%.0f ", 255.0 * intensity);
      t++;
      if ((t & 31) == 0) {
	fprintf(f, "\n");
      }
    }

    // reset x
    for(j = 0; j < x; j++) {
      w[j] = 0;
    }
  }

  free(w);
}

/* Saves a PPM (color pixmap) of specified dimensions of A, with an optional DM decomposition */
void spasm_save_ppm(FILE *f, int x, int y, const spasm *A, const spasm_partition *P) {
  int i, j, jj, n, m, t, p, u, v;
  int *Aj, *Ap, *w, *rr, *cc;
  double max, r, g, b, k;

  int colors[16] = { 0xFF0000, 0xFF6633, 0xCC0000, 0x990000,
		     0xFFFFFF, 0xFFFFFF, 0xFFCC00, 0xCC9900,
		     0xFFFFFF, 0xFFFFFF, 0xFFFFFF, 0x99FF99,
		     0xFFFFFF, 0xFFFFFF, 0xFFFFFF, 0x33CC33 };

  assert(f != NULL);
  assert(A != NULL);

  Aj = A->j;
  Ap = A->p;
  n  = A->n;
  m  = A->m;

  x = spasm_min(x, m);
  y = spasm_min(y, n);

  w = spasm_malloc(x * sizeof(int));
  for(j = 0; j < x; j++) {
    w[j] = 0;
  }

  if (P != NULL) {
    rr = (int *) P->rr;
    cc = (int *) P->cc;
  }

  fprintf(f, "P3\n");
  fprintf(f, "%d %d\n", x, y);
  fprintf(f, "255\n");

  max = (1.0 * m / x) * (1.0 * n / y);
  t = 0;
  i = 0;
  while(i < n) {
    for(k = 0; k < spasm_max(1, n / y) && i < n; k++) {

      // scatters row i to x
      for(p = Ap[i]; p < Ap[i + 1]; p++) {
	w[ (Aj[p] * x) / m ]++;
      }
      i++;
    }

    // print row
    for(j = 0; j < x; j++) {

      jj = j * m / x;
      for(t = 0; t < 4; t++) {
	if (P != NULL && (i > rr[t])) {
	  u = t;
	}
	if (P != NULL && (jj >= cc[t])) {
	  v = t;
	}
      }

      r = (colors[4 * u + v] >> 16) & 0xff;
      g = (colors[4 * u + v] >> 8)  & 0xff;
      b =  colors[4 * u + v]        & 0xff;

      k = 1.0 - w[j] / max;
      if (k < 0.0 | k > 1.0 ) {
	printf("ça craint : %f (%d vs %f)\n", k, w[j], max);
	exit(1);
      }

      fprintf(f, "%.0f %.0f %.0f ", r * k, g * k, b * k);
      t++;
      if ((t & 7) == 0) {
	fprintf(f, "\n");
      }
    }

    // reset x
    for(j = 0; j < x; j++) {
      w[j] = 0;
    }
  }

  fprintf(f, "\n");
  free(w);
}
