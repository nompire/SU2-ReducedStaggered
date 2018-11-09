// -----------------------------------------------------------------
// Add two vectors
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void add_vec(vector *a, vector *b, vector *c) {

  CADD(a->c[0],b->c[0],c->c[0]);
  CADD(a->c[1],b->c[1],c->c[1]);
  
}
// -----------------------------------------------------------------
