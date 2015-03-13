#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "spasm.h"
#include "metis.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *B;
  int n, m, i, j, result, *Ap, *Aj, k, npart, edgecut, *row_part, *col_part, *p, *q, *qinv, px, ch, output;

  /* options descriptor */
  struct option longopts[2] = {
    { "matrix",        no_argument,       NULL,    'm' },
    { NULL,            0,                 NULL,     0  }
  };

  output = 0;

  while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
    switch (ch) {
    case 'm':
      output = 1;
      break;
    default:
      printf("Unknown option\n");
      exit(1);
    }
  }
  argc -= optind;
  argv += optind;

  T = spasm_load_sms(stdin, (output == 0) ? -1 : 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  // element = row
  // node = column

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;
  npart = 2;
  row_part = spasm_malloc(n * sizeof(int));
  col_part = spasm_malloc(m * sizeof(int));

  p = spasm_malloc(n * sizeof(int));
  q = spasm_malloc(m * sizeof(int));

  result = METIS_PartMeshNodal( &n, &m, Ap, Aj, NULL, NULL, &npart, NULL, NULL, &edgecut, row_part, col_part);

  if (result != METIS_OK) {
    fprintf(stderr, "something went wrong\n");
    exit(1);
  }

  /* find (row) separator */
  k = 0;
  for(i = 0; i < n; i++) {
    for(px = Ap[i]; px < Ap[i + 1]; px++) {
      // row i belongs to separator ?
      if (row_part[i] != col_part[ Aj[px] ] ) {
	row_part[i] = npart;
	k++;
	break;
      }
    }
  }

  /* for(i = 0; i < n; i++) {
    printf("%d : %d\n", i, row_part[i]);
  }
  printf("------------------------\n");
  for(i = 0; i < m; i++) {
    printf("%d : %d\n", i, col_part[i]);
    }*/
  

  if (output == 0) {
    printf("%d\t %.1f %%\n", k, 100.0 * k / n);
  } else {
    /* we output the permuted matrix */

    /* build p */
    px = 0;
    for(k = 0; k < npart + 1; k++) {
      for(i = 0; i < n; i++) {
	if (row_part[i] == k) {
	  p[px] = i;
	  px++;
	}
      }
    }

    // deal with empty rows
    for(i = 0; i < n; i++) {
      if (row_part[i] < 0) {
	p[px] = i;
	px++;
      }
    }
    assert(px == n);

    /* build q */
    px = 0;
    for(k = 0; k < npart; k++) {
      for(j = 0; j < m; j++) {
	if (col_part[j] == k) {
	  q[px] = j;
	  px++;
	}
      }
    }
    assert(px == m);

    free(row_part);
    free(col_part);

    qinv = spasm_pinv(q, m);
    B = spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);

    free(q);
    free(qinv);

    spasm_save_csr(stdout, B);
    spasm_csr_free(B);
  }

  spasm_csr_free(A);
  return 0;
}
