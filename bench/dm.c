#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#include "spasm.h"

int nontrivial_diag_size = 0;
int nontrivial_diag_rank = 0;
int trivial_diag_rank = 0;


int subrank(const spasm *M, int a, int b, int c, int d) {
  spasm *C;
  int *p;
  spasm_lu *LU;
  int r;

  C = spasm_submatrix(M, a, c, b, d, SPASM_WITH_NUMERICAL_VALUES);
  //  printf("M : %d, C : %d nnz, dim %d x %d\n", spasm_nnz(M), spasm_nnz(C), C->n, C->m);
  p = spasm_cheap_pivots(C);
  LU = spasm_LU(C, p, SPASM_DISCARD_L);
  free(p);
  r = LU->U->n;
  spasm_free_LU(LU);
  spasm_csr_free(C);
  return r;
}

void show(const spasm *M, spasm_cc *Y) {
  int i, j, a,b,c,d,e,f,g,h, r;

  for(i = 0; i < Y->CC->nr; i++) {
    a = Y->CC->rr[i];
    b = Y->CC->cc[i];
    c = Y->CC->rr[i + 1];
    d = Y->CC->cc[i + 1];

    //    printf("   *) Connected component (%d x %d) --- (%d, %d) to (%d, %d)\n", c-a, d-b, a,b,c,d);
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	e = Y->SCC[i]->rr[j];
	f = Y->SCC[i]->cc[j];
	g = Y->SCC[i]->rr[j + 1];
	h = Y->SCC[i]->cc[j + 1];
	r = subrank(M, e, f, g, h);
	if (g-e > 1) {
	  nontrivial_diag_size += g-e;
	  nontrivial_diag_rank += r;
	}
	if  (g-e > 1) {
	  trivial_diag_rank += 1;
	}
	//	printf("       *) SCC (%d x %d, deffect %d) --- (%d, %d) to (%d, %d)\n", g-e, h-f, spasm_min(g-e, h-f)-r, e, f, g, h);
      }
    }
  }
}

int largest_diagonal_block(spasm_cc *Y) {
  int i, j, r = -1;

  for(i = 0; i < Y->CC->nr; i++) {
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	int SCC_n = Y->SCC[i]->rr[j + 1] - Y->SCC[i]->rr[j];
	r = spasm_max(r, SCC_n);
      }
    }
  }
  return r;
}


int main(int argc, char **argv) {
    spasm_triplet *T;
    spasm *A, *B;
    spasm_dm *x;
    int n, m, i, *qinv, verbose, ch, *rr, *cc;
    char *pm_file, *img_file;
    FILE *f;

  /* options descriptor */
  struct option longopts[6] = {
    { "permuted",      required_argument, NULL,    'p' },
    { "verbose",       no_argument,       NULL,    'v' },
    { "tabulated",     no_argument,       NULL,    't' },
    { "image",         required_argument, NULL,    'i' },
    { NULL,            0,                 NULL,     0  }
  };

  pm_file = NULL;
  img_file = NULL;
  verbose = 0;
  B = NULL;

  while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
    switch (ch) {
    case 'p':
      pm_file = optarg;
      break;
    case 'v':
      verbose = 2;
      break;
    case 't':
      verbose = 1;
      break;
    case 'i':
      img_file = optarg;
      break;
    default:
      printf("Unknown option\n");
      exit(1);
    }
  }
  argc -= optind;
  argv += optind;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

    n = A->n;
    m = A->m;

    x = spasm_dulmage_mendelsohn(A);
    rr = x->DM->rr;
    cc = x->DM->cc;

    qinv = spasm_pinv(x->DM->q, m);
    B = spasm_permute(A, x->DM->p, qinv, SPASM_WITH_NUMERICAL_VALUES);
    free(qinv);

    /* terse output with just the size of the largest diagonal block */
    if (verbose == 1) {
      i = -1;
      if (x->H != NULL) {
	i = spasm_max(i, largest_diagonal_block(x->H));
      }
      if (x->S != NULL) {
	i = spasm_max(i, largest_diagonal_block(x->S));
      }
      if (x->V != NULL) {
	i = spasm_max(i, largest_diagonal_block(x->V));
      }
      printf("%5d \t %5d \t %6d \t %6d \t %.1f\n", n, m, spasm_nnz(A), i, 100.0 * i / spasm_min(n, m));
    }

    if (verbose == 2) {
      printf("structural rank = %d\n", rr[2] + cc[4] - cc[3]);
      if (x->H != NULL) {
	int h_n = rr[1] - rr[0];
	int h_m = cc[2] - cc[0];
	printf("*) H (%d x %d) : \n", h_n, h_m);
	show(B, x->H);
      }
      if (x->S != NULL) {
	int s_n = rr[2] - rr[1];
	int s_m = cc[3] - cc[2];
	printf("*) S (%d x %d) : \n", s_n, s_m);
	show(B, x->S);
      }
      if (x->V != NULL) {
	int v_n = rr[4] - rr[2];
	int v_m = cc[4] - cc[3];
	printf("*) V (%d x %d) : \n", v_n, v_m);
	show(B, x->V);
      }
    }


    if (pm_file != NULL) {
      f = fopen(pm_file, "w");
      spasm_save_csr(f, B);
      fclose(f);
    }

    if (img_file != NULL) {
      FILE *f = fopen(img_file, "w");
      spasm_save_ppm(f, B, x);
      fclose(f);
    }

    spasm_csr_free(B);
    spasm_csr_free(A);

    printf("non trivial diagonal blocks size : %d\n", nontrivial_diag_size);
    printf("non trivial diagonal blocks rank : %d\n", nontrivial_diag_rank);
    printf("trivial diagonal blocks : %d\n", trivial_diag_rank);
    return 0;
}
