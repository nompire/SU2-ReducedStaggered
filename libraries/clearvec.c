// -----------------------------------------------------------------
// Clear a vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void clearvec(vector *v) {

  register int i;
  for (i = 0; i < DIMF; i++){
    v->c[i].real = 0.0;
    v->c[i].imag = 0.0;}

}
// ----------------------------------------------------------------
