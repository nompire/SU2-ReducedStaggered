/******************  s_m_mat.c  (in su2.a) ******************************
*									*
* void scalar_mult_su2_matrix( su2_matrix *a, Real s, su2_matrix *b)	*
* B <- s*A								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

/* b <- s*a, matrices */
void scalar_mult_su2_matrix( su2_matrix *a, double s, su2_matrix *b ){

register int i,j;
    for(i=0;i<2;i++)for(j=0;j<2;j++){
	b->e[i][j].real = s*a->e[i][j].real;
	b->e[i][j].imag = s*a->e[i][j].imag;
    }
}


