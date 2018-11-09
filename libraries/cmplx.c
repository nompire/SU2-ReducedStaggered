/********************** cmplx.c (in complex.a) **********************/
/* Subroutines for operations on complex numbers */
/* make a complex number from two real numbers */
#include "../include/config.h"
#include "../include/complex.h"

complex cmplx( Real x, Real y )  {
    complex c;
    c.real = x; c.imag = y;
    return(c);
}
