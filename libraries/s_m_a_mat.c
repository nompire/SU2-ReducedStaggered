/****************  s_m_a_mat.c  (in su2.a) ******************************
*									*
* void scalar_mult_add_su2_matrix( su2_matrix *a, su2_matrix *b,	*
*	Real s, su2_matrix *c)						*
* C <- A + s*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

/* c <- a + s*b, matrices */
void scalar_mult_add_su2_matrix(su2_matrix *a,su2_matrix *b,double s,
	su2_matrix *c){


register int i,j;
    for(i=0;i<2;i++){
    for(j=0;j<2;j++){
	c->e[i][j].real = a->e[i][j].real + s*b->e[i][j].real;
	c->e[i][j].imag = a->e[i][j].imag + s*b->e[i][j].imag;
       }
    }

}
