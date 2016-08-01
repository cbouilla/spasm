/*
 * indent -nfbs -i2 -nip -npsl -di0 -nut spasm_lu.c uetree.c
 * 
 * Created on: July, 2011 updated, August 2012 Author: Kamer Kaya and Bora Ucar
 * 
 * implements the algorithm given in On constructing elimination trees for
 * sparse unsymmetric matrices by Kamer Kaya and Bora Ucar
 * 
 */

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

#define max(x, y)  (((x) > (y)) ? (x) : (y))
#define min(x, y)  (((x) > (y)) ? (y) : (x))

/* returns the number of SCC.
   perm comes straight from the input
   cptrs is uninitialized on firs call;
   last 4 arguments are temporary storage */
int scc(int *Aptrs, int *Aids, int n, int *perm, int *cptrs, int *lowl, int *spos, int *prev, int *tedges) {
  int i, iv, iw, j, i2, ptr, lcnt, ist, stp, icnt, num, control, tnm;

  control = icnt = num = 0;
  tnm = 2 * n - 1;

  /* initialization */
  for (i = 0; i < n; i++) {
    spos[i] = -1;
  }
  memcpy(tedges, Aptrs, sizeof(int) * n);

  /* main loop */
  for (i = 0; i < n; i++) {
    /* if node i is already marked, skip */
    if (spos[i] != -1) {
      continue;
    }

    /* else, start a DFS from node i */

    iv = i;
    ist = 1;

    /* put the node to the stack */
    lowl[iv] = spos[iv] = 0;
    cptrs[n - 1] = iv;

    /* the body of this loop either puts a new node to the stack or backtrack */
    for (j = 0; j < tnm; j++) {
      ptr = tedges[iv];

      /* if there exists an edge to visit */
      if (ptr < Aptrs[iv + 1]) {
        i2 = Aptrs[iv + 1];
        control = 0;

        /*
         * search the edges leaving iv until one enters a new node or all
         * edges are exhausted
         */
        for (; ptr < i2; ptr++) {
          iw = Aids[ptr];
          if (iw < n) {
            /* check if node iw has not been on stack already */
            if (spos[iw] == -1) {
              /* put a new node to the stack */
              tedges[iv] = ptr + 1;
              prev[iw] = iv;
              iv = iw;

              lowl[iv] = spos[iv] = ist = ++ist;
              cptrs[n - ist] = iv;
              control = 1;
              break;
            }
            /* update lowl[iw] if necessary */
            if (lowl[iw] < lowl[iv]) {
              lowl[iv] = lowl[iw];
            }
          }
        }

        if (control == 1) {
          control = 0;
          continue;
        }
      }
      /* is node iv the root of a block */
      if (lowl[iv] >= spos[iv]) {
        /* order the nodes in the block */
        lcnt = icnt;

        /*
         * peel block off the top of the stack starting at the top and
         * working down to the root of the block
         */
        for (stp = n - ist; stp < n; stp++) {
          iw = cptrs[stp];
          lowl[iw] = n;
          spos[iw] = icnt++;
          if (iw == iv) {
            break;
          }
        }

        ist = (n - 1) - stp;
        cptrs[num++] = lcnt;

        /* are there any nodes left on the stack */
        if (ist == 0) {
          /* if all the nodes have been ordered */
          if (icnt == n) {
            control = 1;
          }
          break;
        }
      }
      /* backtrack to previous node on a path */
      iw = iv;
      iv = prev[iv];

      /* update the value of lowl(iv) if necessary */
      if (lowl[iw] < lowl[iv]) {
        lowl[iv] = lowl[iw];
      }
    }
    if (control == 1) {
      break;
    }
  }

  /* put permutation in the required form */
  for (i = 0; i < n; i++) {
    perm[spos[i]] = i;
  }
  cptrs[num] = n;

  return num;
}

/* #define DEBUG */

#ifdef DEBUG
int nnz(int *a, int size) {
  int i;
  int nz = 0;
  for (i = 0; i < size; i++) {
    if (a[i] != 0) {
      nz++;
    }
  }
  return nz;
}

int neq(int *a, int val, int size) {
  int i;
  int ne = 0;
  for (i = 0; i < size; i++) {
    if (a[i] == val) {
      ne++;
    }
  }
  return ne;
}

int npos(int *a, int size) {
  int i;
  int np = 0;
  for (i = 0; i < size; i++) {
    if (a[i] > 0) {
      np++;
    }
  }
  return np;
}
#endif


/* Aptrs --> Ap
   Aids  --> Ai CSC 
   parent --> the output 
   n --> matrix size
   s --> FL pivots to skip 
   vids --> ???? [permutation des lignes ? Vertex ID ?]
   perm --> ???
   blocks --> output ? 
   work, work2, work3, work4 --> temporary storage (size ?)
   xAids --> ??
   xShrPtrs --> ??
   */
void setparentv2(int *Aptrs, int *Aids, int *parent, int n, int s, int *vids,
      int *perm, int *blocks, int *work, int *work2, int *work3, int *work4,
                      int *xAids, int *xShrPtrs) {
  int i, root, k, num, ik, j, nnz, bi, bj, iv, jv, temp, temp2, ptr, ns,
      nk, sval, sctemp, totnnz, shrtemp;
  int *xAptrs, *cptrs, *dupflag;

  /*
   * printf("calling with n = %d nnz = %d s = %d rootsvid = %d\n", n,
   * Aptrs[n],s, root);
   */
  if (n == 2) {
    parent[vids[0]] = vids[1];
    return;
  } else if (s == n - 1) {
    root = vids[n - 1];
    for (i = 0; i < n - 1; i++) {
      parent[vids[i]] = root;
    }
    return;
  } else {
    nnz = Aptrs[n];
    cptrs = (int *)malloc((n + 1) * sizeof(int)); /* column pointers ? */

    k = s;
    while (k == s) {
      k = (int)(ceil((s + n) / 2.0));
      /*
       * mexPrintf("calling SCC for Gk (containing the first %d vertices) ->
       * ", k);
       */
      num = scc(Aptrs, Aids, k, perm, cptrs, work, work2, work3, work4);
      /* mexPrintf("number of SCs: %d\n", num); */

      if (num == k) {
        s = k;
      }
      if (s == n - 1) {
        root = vids[n - 1];
        for (i = 0; i < n - 1; i++) {
          parent[vids[i]] = root;
        }
        free(cptrs);
        return;
      }
    }

    nk = num - k;
    ns = n + nk;

    /* partition starts here ******************************************** */
    /* find block of each vertex and max edge count of each block        */
    xAptrs = (int *)malloc((n + num) * sizeof(int));

    /*
     * used to mark existence of an edge from some vertex to another one in
     * the shrinked graph
     */
    for (i = 0; i < num; i++) {
      for (j = cptrs[i]; j < cptrs[i + 1]; j++) {
        blocks[perm[j]] = num - i - 1;
      }
    }
    for (i = k; i < n; i++) {
      blocks[i] = i + nk;
    }                           /* this is the vertex id in the shrinked
                                 * graph for each vertex >=k */

    /* find new vertex ids and partition perm */
    memset(work, 0, sizeof(int) * num);
    for (i = 0; i < k; i++) {
      bi = (num - 1) - blocks[i];
      perm[work[bi] + cptrs[bi]] = i;
      work3[i] = work[bi]++;
    }

    /* set new edge list */
    dupflag = work2;
    for (i = 0; i < ns; i++)
      dupflag[i] = -1;

    xShrPtrs[ns] = shrtemp = nnz;
    for (i = n - 1; i >= k; i--) {
      for (ptr = Aptrs[i]; ptr < Aptrs[i + 1]; ptr++) {
        bj = blocks[Aids[ptr]];
        if (bj < num) {         /* >=k to SC */
          if (dupflag[bj] != i) {
            xAids[--shrtemp] = bj;
            dupflag[bj] = i;
          }
        } else {                /* >=k to >=k */
          xAids[--shrtemp] = bj;
        }
      }
      xShrPtrs[i + nk] = shrtemp;
    }

    sctemp = 0;
    totnnz = 0;
    for (i = 0; i < num; i++) { /* for each component */
      temp = (num - 1) - i;
      dupflag[temp] = i;
      for (j = cptrs[i]; j < cptrs[i + 1]; j++) {       /* visit the vertices in
                                                         * the component */
        iv = perm[j];           /* vertex id */
        xAptrs[j + i] = sctemp; /* ptrs for the first vertex of the ith block */
        for (ptr = Aptrs[iv]; ptr < Aptrs[iv + 1]; ptr++) {
          bj = blocks[Aids[ptr]];
          if (bj < num) {       /* SC to SC */
            if (temp == bj) {   /* inside a SC */
              xAids[totnnz + sctemp++] = work3[Aids[ptr]];
            } else if (dupflag[bj] != i) {      /* to other SCs */
              xAids[--shrtemp] = bj;
              dupflag[bj] = i;
            }
          } else {              /* SC to >= k */
            if (dupflag[bj] != i) {
              xAids[--shrtemp] = bj;
              dupflag[bj] = i;
            }
          }
        }
      }
      xShrPtrs[(num - 1) - i] = shrtemp;        /* ptrs for shrinked graph */
      totnnz += sctemp;
      xAptrs[j + i] = sctemp;
      sctemp = 0;
    }
    /* update shrink ptrs */
    for (i = 0; i <= ns; i++) {
      xShrPtrs[i] -= shrtemp;
    }
    /* partition ends here ******************************************** */

    /* recursive call for each strong component************************ */
    /* partition the original ids now */
    memcpy(work, vids, k * sizeof(int));
    for (i = 0; i < k; i++) {
      vids[i] = work[perm[i]];
    }

    totnnz = 0;
    for (i = 0; i < num; i++) {
      if (cptrs[i + 1] - cptrs[i] > 1) {
        sval = 0;
        for (j = cptrs[i]; j < cptrs[i + 1]; j++) {
          if (perm[j] < s)
            sval++;
          else
            break;
        }
        /*
         * printf("%dth recursive call with %d vertices s = %d nnz = %d\n",
         * i, cptrs[i+1] - cptrs[i], sval, xAptrs[i + cptrs[i+1]]);
         */
        setparentv2(&(xAptrs[cptrs[i] + i]), &(xAids[totnnz]), parent,
                    cptrs[i + 1] - cptrs[i], sval, &(vids[cptrs[i]]),
                    &(perm[cptrs[i]]), blocks, work, work2, work3, work4,
                    &(Aids[totnnz]), Aptrs);
      }
      totnnz += xAptrs[i + cptrs[i + 1]];
    }
    /******************************************************************/

    for (i = 0; i < num; i++)
      vids[i] = vids[cptrs[i + 1] - 1];

    free(cptrs);
    free(xAptrs);

    /* recursive call for condensed graph **************************** */
    /* mexPrintf("Shrinked recursive call with %d vertices\n", n - k + num); */
    memcpy(work, vids, num * sizeof(int));
    for (i = num - 1; i >= 0; i--)
      vids[i] = work[(num - 1) - i];
    for (i = k; i < n; i++) {
      vids[i + nk] = vids[i];
    }
    setparentv2(xShrPtrs, &(xAids[shrtemp]), parent, ns, num, vids, perm, blocks, work, work2, work3, work4, Aids, Aptrs);
    /******************************************************************/
  }
}
