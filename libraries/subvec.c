// -----------------------------------------------------------------
// Subtract two vectors
// c <-- a - b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void sub_vec(vector *a, vector *b, vector *c) {

   CSUB(a->c[0],b->c[0],c->c[0]);
   CSUB(a->c[1],b->c[1],c->c[1]);
  
}
// -----------------------------------------------------------------
