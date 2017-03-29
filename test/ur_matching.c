#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(){
 spasm_triplet *T;
 spasm *A;
 spasm_rc *match;
 int i, *col_TO, countcol;

 T = spasm_load_sms(stdin, 42013);
 A = spasm_compress(T);
 spasm_triplet_free(T);

 match = spasm_rc_alloc(A->n, A->m);
 col_TO = spasm_malloc(A->m *sizeof(int));
 
 spasm_ur_matching(A, match, col_TO, &countcol);


 for(i = 0; i < A->n; i++){
   if((match->r[i] != -1) && (match->c[match->r[i]] != i)){
     printf("not ok row %d: col %d, error %d\n", i, match->r[i], match->c[match->r[i]]);
   }
 }

 for(i = 0; i < A->m; i++){
   if((match->c[i] != -1) && (match->r[match->c[i]] != i)){
     printf("not ok col %d: row %d, error %d\n", i, match->c[i], match->r[match->c[i]]);
   }
 }
 
 printf("Ok\n");
 spasm_csr_free(A);
 spasm_rc_free(match);
 free(col_TO);
}
