/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_dm.c */
#include <assert.h>
#include "spasm.h"


/*
 * R0 : unmatched rows C0 : unmatched columns
 * 
 * R1 : rows reachable from C0 (by alternating paths) C1 : column matched to R1
 * 
 * C3 : columns reachable from R0 (by alternating paths) R3 : rows matched to C3
 * 
 * R2 : other rows C2 : other columns
 */


/*
 * breadth-first search for coarse decomposition. This determines R0,R3,C3
 * (or C0,C1,R1 when given the transpose of A).
 */
static void spasm_bfs(const spasm * A, int *wi, int *wj, int *queue, const int *imatch, const int *jmatch, int mark) {
  int *Ap, *Aj;
  int n, head, tail, j, i, p, i2;

  head = 0;
  tail = 0;
  Ap = A->p;
  Aj = A->j;
  n = A->n;

  /* place all unmatched nodes in queue */
  for (i = 0; i < n; i++) {
    /* skip row i if matched */
    if (jmatch[i] >= 0) {
      continue;
    }
    /* i in set R0 (C0 if transpose) */
    wi[i] = 0;

    /* place unmatched row i in queue */
    queue[tail] = i;
    tail++;
  }

  /* quick return if no unmatched nodes */
  if (tail == 0) {
    return;
  }
  /* while queue is not empty */
  while (head < tail) {

    /* get the head of the queue */
    i = queue[head];
    head++;

    for (p = Ap[i]; p < Ap[i + 1]; p++) {
      j = Aj[p];

      /* skip if j is marked */
      if (wj[j] >= 0) {
        continue;
      }
      /* j in set C3 (R1 if transpose) */
      wj[j] = mark;

      /* traverse alternating path to i2 */
      i2 = imatch[j];

      /* skip i2 if it is marked */
      if (wi[i2] >= 0) {
        continue;
      }
      /* i2 in set R3 (C1 if transpose) */
      wi[i2] = mark;

      /* add i2 to queue */
      queue[tail] = i2;
      tail++;
    }
  }
}


/* collect unmatched rows into the permutation vector p */
static void spasm_unmatched(int n, const int *wi, int *p, int *rr, int set) {
  int i, kr;

  kr = rr[set];
  for (i = 0; i < n; i++) {
    if (wi[i] == 0) {
      p[kr] = i;
      kr++;
    }
  }
  rr[set + 1] = kr;
}



/* collect matched rows and columns into p and q */
static void spasm_matched(int n, const int *wj, const int *imatch, int *p, int *q, int *cc, int *rr, int set, int mark) {
  int kc, kr, j;

  kc = cc[set];
  kr = rr[set - 1];
  for (j = 0; j < n; j++) {
    if (wj[j] != mark) {
      continue;                 /* skip if j is not in C set */
    }
    p[kr++] = imatch[j];
    q[kc++] = j;
  }
  cc[set + 1] = kc;
  rr[set] = kr;
}


/* Given A and a maximum matching, compute coarse dm decomposition */
static spasm_partition *spasm_dulmage_mendelson_coarse(const spasm * A, const spasm * A_t, const int *jmatch, const int *imatch) {
  int m, n, i, j;
  int *wi, *wj, *p, *q, *cc, *rr, *queue;
  spasm_partition *P;

  /* check inputs */
  assert(A != NULL);
  assert(A_t != NULL);
  n = A->n;
  m = A->m;

  /* allocate result */
  wi = spasm_malloc(n * sizeof(int));
  wj = spasm_malloc(m * sizeof(int));
  queue = spasm_malloc(spasm_max(n, m) * sizeof(int));

  P = spasm_partition_alloc(n, m, 4, 4);
  p = P->p;
  q = P->q;
  cc = P->cc;
  rr = P->rr;
  P->nr = 4;
  P->nc = 4;
  for (i = 0; i <= 4; i++) {
    rr[i] = 0;
    cc[i] = 0;
  }

  /* --- Reachability --------------------------------------------- */

  /* unmark everything for bfs */
  for (j = 0; j < m; j++) {
    wj[j] = -1;
  }
  for (i = 0; i < n; i++) {
    wi[i] = -1;
  }

  /* find R3, C3 from R0 */
  spasm_bfs(A, wi, wj, queue, imatch, jmatch, 3);

  /* find R1, C1 from C0 */
  spasm_bfs(A_t, wj, wi, queue, jmatch, imatch, 1);

  spasm_unmatched(m, wj, q, cc, 0);     /* unmatched set C0 */
  spasm_matched(m, wj, imatch, p, q, cc, rr, 1, 1);     /* set R1 and C1 */
  spasm_matched(m, wj, imatch, p, q, cc, rr, 2, -1);    /* set R2 and C2 */
  spasm_matched(m, wj, imatch, p, q, cc, rr, 3, 3);     /* set R3 and C3 */
  spasm_unmatched(n, wi, p, rr, 3);     /* unmatched set R0 */

  /* cleanup */
  free(queue);
  free(wi);
  free(wj);

  return P;
}

static spasm_partition *process_square_part(const spasm * B, int rx, int ry, int cx, int cy, int *p, int *q, int ra, int ca) {
  spasm *C;
  spasm_partition *SCC;
  int i, k;

  assert(ry - rx == cy - cx);
  C = spasm_submatrix(B, rx, ry, cx, cy, SPASM_IGNORE_VALUES);
  SCC = spasm_strongly_connected_components(C);
  spasm_csr_free(C);
  k = SCC->nr;

  /* update permutations of B */
  spasm_range_pvec(p, rx, ry, SCC->p);
  spasm_range_pvec(q, cx, cy, SCC->q);

  /* shift SCC */
  for (i = 0; i <= k; i++) {
    SCC->rr[i] += rx + ra;
    SCC->cc[i] += cx + ca;
  }

  return SCC;
}

static spasm_cc *process_rectangular_part(const spasm * B, int ra, int rb, int ca, int cb, int *p, int *q, const int *jmatch) {
  spasm *M, *MM;
  int n, m, CC_k, i, k, rx, ry, cx, cy, C_n, C_m;
  int *M_jmatch;
  int *CC_qinv;
  spasm_partition *CC;
  spasm_cc *result;


  M = spasm_submatrix(B, ra, rb, ca, cb, SPASM_IGNORE_VALUES);
  M_jmatch = spasm_submatching(jmatch, ra, rb, ca, cb);

  n = M->n;
  m = M->m;

  if (n == 0 || m == 0) {
    return NULL;
  }
  /* --- connected components of M */
  CC = spasm_connected_components(M, NULL, M_jmatch);
  CC_k = CC->nr;

  /* permute M to expose the connected components */
  CC_qinv = spasm_pinv(CC->q, m);
  MM = spasm_permute(M, CC->p, CC_qinv, SPASM_IGNORE_VALUES);
  free(CC_qinv);

  result = spasm_malloc(sizeof(spasm_cc));
  result->CC = CC;
  result->SCC = spasm_malloc(CC_k * sizeof(spasm_partition *));

  for (i = 0; i < CC_k; i++) {

    /* process i-th connected component of M */
    result->SCC[i] = NULL;
    C_n = CC->rr[i + 1] - CC->rr[i];
    C_m = CC->cc[i + 1] - CC->cc[i];

    if (C_n == 0 || C_m == 0) {
      continue;
    }
    /* extract the (square) perfectly-matched part */
    k = spasm_min(C_n, C_m);
    cx = CC->cc[i];
    ry = CC->rr[i + 1];
    if (C_n == C_m) {
      /* square case: the matching is perfect */
      rx = CC->rr[i];
      cy = CC->cc[i + 1];
    } else if (C_n < C_m) {
      /* horizontal case: matched columns are on the left */
      rx = CC->rr[i];
      cy = CC->cc[i] + k;
    } else {
      /* vertical case: matched rows are on the bottom */
      rx = CC->rr[i + 1] - k;
      cy = CC->cc[i + 1];
    }

    result->SCC[i] = process_square_part(MM, rx, ry, cx, cy, CC->p, CC->q, ra, ca);
  }

  /* update permutations of B */
  spasm_range_pvec(p, ra, rb, CC->p);
  spasm_range_pvec(q, ca, cb, CC->q);

  /* translate CC */
  for (i = 0; i <= CC_k; i++) {
    CC->rr[i] += ra;
    CC->cc[i] += ca;
  }

  /* cleanup */
  free(M_jmatch);
  spasm_csr_free(MM);
  spasm_csr_free(M);

  return result;
}

spasm_dm *spasm_dulmage_mendelsohn(const spasm * A) {
  spasm *B, *A_t;
  spasm_partition *DM;
  int n, m;
  int *p, *q, *rr, *cc;
  int *qinv;
  int *imatch, *jmatch, *Bjmatch;
  spasm_dm *result;

  assert(A != NULL);
  n = A->n;
  m = A->m;

  A_t = spasm_transpose(A, SPASM_IGNORE_VALUES);

  /* --- Maximum matching ------------------------------------------------- */
  jmatch = spasm_malloc(n * sizeof(int));
  imatch = spasm_malloc(m * sizeof(int));

  if (n < m) {
    spasm_maximum_matching(A, jmatch, imatch);
  } else {
    spasm_maximum_matching(A_t, imatch, jmatch);
  }

  /* --- coarse DM decomposition ----------------------------------------- */
  DM = spasm_dulmage_mendelson_coarse(A, A_t, jmatch, imatch);
  p = DM->p;
  q = DM->q;
  rr = DM->rr;
  cc = DM->cc;

  spasm_csr_free(A_t);
  free(imatch);

  result = spasm_malloc(sizeof(spasm_dm));
  result->DM = DM;
  result->H = NULL;
  result->S = NULL;
  result->V = NULL;

  qinv = spasm_pinv(q, m);
  B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);

  /* optimization : this could be done in-place */
  Bjmatch = spasm_permute_row_matching(n, jmatch, p, qinv);
  free(qinv);
  free(jmatch);

  /* ------------------- H --------------------- */
  if (cc[2] != cc[0]) {
    result->H = process_rectangular_part(B, rr[0], rr[1], cc[0], cc[2], p, q, Bjmatch);
  }
  /* --------------- S ----------------------- */
  if (rr[2] != rr[1]) {
    result->S = process_rectangular_part(B, rr[1], rr[2], cc[2], cc[3], p, q, Bjmatch);
  }
  /* ------------------- V --------------------- */
  if (rr[4] != rr[2]) {
    result->V = process_rectangular_part(B, rr[2], rr[4], cc[3], cc[4], p, q, Bjmatch);
  }
  /* cleanup */
  spasm_csr_free(B);
  free(Bjmatch);

  return result;
}
