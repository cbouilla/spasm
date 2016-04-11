#include <assert.h>
#include <stdio.h>
#include "spasm.h"

int main(){
  
  // charge la matrice depuis l'entrée standard
  int prime = 42013;
  spasm_triplet * T = spasm_load_sms(stdin, prime);
  spasm * A = spasm_compress(T);
  spasm_triplet_free(T);

  // permutation p : cheap pivots.



}
