// -----------------------------------------------------------------
// Print the given vector
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su2.h"

void dumpvec(vector *v) {

  register int i;
  
  for (i = 1; i < DIMF; i++)
    printf("  %.4g", v->c[i].real);
    printf("  %.4g", v->c[i].imag);
    printf("\n");

}
// -----------------------------------------------------------------
