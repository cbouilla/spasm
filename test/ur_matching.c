#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(){
 spasm_triplet *T;
 spasm *A;
 spasm_rc *match;
 int i; 

 T = spasm_load_sms(stdin, 42013);
 A = spasm_compress(T);
 spasm_triplet_free(T);

 match = spasm_ur_matching(A);

 for(i = 0; i < A->n; i++){
   if((match->r[i] != -1) && (match->c[match->r[i]] != i)){
     printf("not ok row %d: col %d, error %d\n", i, match->r[i], match->c[match->r[i]]);
   }
 }
 
 printf("Ok\n");
 spasm_csr_free(A);
 spasm_rc_free(match);
}
