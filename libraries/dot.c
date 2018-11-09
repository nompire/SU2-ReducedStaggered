// -----------------------------------------------------------------
// Return the dot product of two vectors: adag b
#include "../include/config.h"
#include "../include/su2.h"
#include "../include/complex.h"
#include <stdio.h>
#include <stdlib.h>

double dot(vector *a, vector *b) {
register double temp1,temp2;
    temp2 = a->c[0].real * b->c[0].real;
    temp1 = a->c[0].imag * b->c[0].imag; temp2 += temp1;
    temp1 = a->c[1].real * b->c[1].real; temp2 += temp1;
    temp1 = a->c[1].imag * b->c[1].imag; temp2 += temp1;
    return(temp2);
}
// -----------------------------------------------------------------
