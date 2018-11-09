// -----------------------------------------------------------------
// Clear a vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void vec_conjug(vector *v,vector *u) {

  register int i;
  for (i = 0; i < DIMF; i++){
    CONJG( v->c[i], u->c[i]);
  }

}
// ----------------------------------------------------------------
