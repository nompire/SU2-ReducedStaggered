/****************  cs_m_s_mat.c  (in su2.a) *****************************
*									*
* void c_scalar_mult_sub_su2mat( su2_matrix *a, su2_matrix *b,		*
*	complex *s, su2_matrix *c)					*
* C <- A - s*B,   A,B and C matrices 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

/* c <- a - s*b, matrices */
void c_scalar_mult_sub_su2mat( su2_matrix *a, su2_matrix *b, complex *s,su2_matrix *c){


register int i,j;
register double sr,si;
sr = (*s).real; si = (*s).imag;
 
    for(i=0;i<2;i++){
    for(j=0;j<2;j++){
	
         c->e[i][j].real = a->e[i][j].real -  (sr*b->e[i][j].real -  si*b->e[i][j].imag);
	 c->e[i][j].imag = a->e[i][j].imag -  (sr*b->e[i][j].imag +  si*b->e[i][j].real); 

    }

  }

}
