/****************  s_m_a_vec.c  (in su2.a) ******************************
*									*
* void scalar_mult_add_su2_vector( su2_vector *a, su2_vector *b,	*
*	Real s, su2_vector *c)						*
* C <- A + s*B,   A,B and C vectors 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

/* c <- a + s*b, vectors */

void scalar_mult_add_vec(vector *a, vector *b, double s,vector *c){

register int i;
  for(i=0;i<2;i++){
    c->c[i].real = a->c[i].real + s*b->c[i].real;
    c->c[i].imag = a->c[i].imag + s*b->c[i].imag;
  }
  
}
