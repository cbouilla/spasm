#include "spasm.h"

int spasm_is_upper_triangular(const spasm *A)
{
        int n = A->n;
        int m = A->m;
        if (n > m)
          return 0;
        const i64 *Ap = A->p;
        const int *Aj = A->j;
        const spasm_GFp*Ax = A->x;
        for (int i = 0; i < n; i++) {
                if (Ap[i] == Ap[i + 1])
                        return 0;
                /* check diagonal */
                if (Aj[Ap[i]] != i)
                        return 0;
                if (Ax[Ap[i]] != 1)
                        return 0;
                /* check other entries */
                for (i64 p = Ap[i] + 1; p < Ap[i + 1]; p++)
                        if (Aj[p] < i)
                    return 0;
        }
        return 1;
}


int spasm_is_lower_triangular(const spasm *A)
{  
        int n = A->n;
        int m = A->m;
        if (n < m)
                return 0;
        const i64 *Ap = A->p;
        const int *Aj = A->j;
        const spasm_GFp *Ax = A->x;
        for (int i = 0; i < m; i++) {
                if (Ap[i] == Ap[i + 1])
                        return 0;
                /* check diagonal */
                if (Aj[Ap[i + 1] - 1] != i)
                        return 0;
                if (Ax[Ap[i + 1] - 1] == 0)
                        return 0;
                /* check other entries */
                for (i64 p = Ap[i]; p < Ap[i + 1] - 1; p++)
                  if (Aj[p] > i)
                    return 0;
        }
        return 1;
}
