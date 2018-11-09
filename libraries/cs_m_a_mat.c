/******************  cs_m_a_mat.c  (in su3.a) ***************************
*									*
*  c_scalar_mult_add_su3mat( su3_matrix *ma, su3_matrix *m2,		*
*	complex *phase, su3_matrix *m3)					*
*  multiply an su3 matrix by a complex scalar and add it to another	*
*  matrix:   m3 <- m1 + number*m2 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void c_scalar_mult_add_su2mat( su2_matrix *m1, su2_matrix *m2,
	complex *phase, su2_matrix *m3){


register int i,j;
register double sr,si;
sr = (*phase).real; si = (*phase).imag;

    for(i=0;i<2;i++){
    for(j=0;j<2;j++){
	 m3->e[i][j].real = m1->e[i][j].real  +   (sr*m2->e[i][j].real - si*m2->e[i][j].imag);
	 m3->e[i][j].imag = m1->e[i][j].imag +   (sr*m2->e[i][j].imag +  si*m2->e[i][j].real);
        }


   }
}
