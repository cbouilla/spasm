#include <assert.h>
#include "spasm.h"

/* allocate a cs_dmperm or cs_scc result */
spasm_dm *spasm_dm_alloc(int n, int m) {
  int i;
    spasm_dm *D;

    D = spasm_malloc(sizeof(spasm_dm));
    D->p = spasm_malloc (n     * sizeof(int)) ;
    D->r = spasm_malloc ((n+6) * sizeof(int)) ;
    D->q = spasm_malloc (m     * sizeof(int)) ;
    D->s = spasm_malloc ((m+6) * sizeof(int)) ;
    D->nb = 0;
    for(i = 0; i < 5; i++) {
      D->rr[i] = 0;
      D->cc[i] = 0;
    }
    return D;
}

void spasm_dm_free(spasm_dm *D) {
  if (D == NULL) {
    return;
  }
  free(D->p);
  free(D->q);
  free(D->r);
  free(D->s);
  free(D);
}


/* R0 : unmatched rows
   C0 : unmatched columns

   R1 : rows reachable from C0 (by alternating paths)
   C1 : column matched to R1

   C3 : columns reachable from R0 (by alternating paths)
   R3 : rows matched to C3

   R2 : other rows
   C2 : other columns
*/


/* breadth-first search for coarse decomposition. This determines R0,R3,C3
 * (or C0,C1,R1 when given the transpose of A).
 */
static void spasm_bfs (const spasm *A, int *wi, int *wj, int *queue, const int *imatch, const int *jmatch, int mark) {
  int *Ap, *Aj;
  int n, head, tail, j, i, p, i2 ;

    head = 0;
    tail = 0;
    Ap = A->p;
    Aj = A->j;
    n = A->n;

    /* place all unmatched nodes in queue */
    for (i = 0 ; i < n ; i++) {
      /* skip row i if matched */
      if (jmatch[i] >= 0) {
	continue;
      }
      /* i in set R0 (C0 if transpose) */
      wi[i] = 0 ;

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
  for (i = 0 ; i < n ; i++) {
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
  for (j = 0 ; j < n ; j++) {
    if (wj[j] != mark) {
      continue ;      /* skip if j is not in C set */
    }
    p [kr++] = imatch[j];
    q [kc++] = j;
  }
  cc[set + 1] = kc;
  rr[set] = kr;
}


/* Given A, compute coarse and then fine dmperm */
spasm_dm * spasm_dulmage_mendelson(const spasm *A) {
  int m, n, i, j, k;
  int *jmatch, *imatch, *wi, *wj, *p, *q, *cc, *rr, *queue;

    spasm_dm *D;
    spasm *A_t;

    /* check inputs */
    assert(A != NULL);
    n = A->n;
    m = A->m;

    /* allocate result */
    jmatch = spasm_malloc(n * sizeof(int));
    imatch = spasm_malloc(m * sizeof(int));
    wi = spasm_malloc(n * sizeof(int));
    wj = spasm_malloc(m * sizeof(int));
    queue = spasm_malloc(spasm_max(n, m) * sizeof(int));

    D = spasm_dm_alloc(n, m);
    p = D->p;
    q = D->q;
    cc = D->cc;
    rr = D->rr;

    /* --- Maximum matching ------------------------------------------------- */
    spasm_maximum_matching(A, imatch, jmatch);

    /* --- Coarse decomposition --------------------------------------------- */

    /* unmark everything for bfs */
    for (j = 0 ; j < m ; j++) {
      wj[j] = -1;
    }
    for (i = 0 ; i < n ; i++) {
      wi[i] = -1;
    }

    A_t = spasm_transpose(A, 0);

    /* find R3, C3 from R0 */
    spasm_bfs(A,    wi, wj, queue, imatch, jmatch, 3);

    /* find R1, C1 from C0 */
    spasm_bfs(A_t,  wj, wi, queue, jmatch, imatch, 1);

    spasm_unmatched(m, wj, q, cc, 0) ;                      /* unmatched set C0 */
    spasm_matched  (m, wj, imatch, p, q, cc, rr, 1, 1) ;    /* set R1 and C1 */
    spasm_matched  (m, wj, imatch, p, q, cc, rr, 2, -1) ;   /* set R2 and C2 */
    spasm_matched  (m, wj, imatch, p, q, cc, rr, 3, 3) ;    /* set R3 and C3 */
    spasm_unmatched(n, wi, p, rr, 3) ;                      /* unmatched set R0 */

    /* cleanup */
    free(imatch);
    free(jmatch);
    free(queue);
    free(wi);
    free(wj);
    spasm_csr_free(A_t);

    return D;
}
