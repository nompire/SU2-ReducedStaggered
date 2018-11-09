/************  uncmp_ahmat.c  (in su2.a) ********************************
*									*
* void uncompress_anti_hermitian( anti_hermitmat *mat_antihermit,	*
*	su2_matrix *mat_su2 )						*
* uncompresses an anti_hermitian matrix to make a 2x2 complex matrix	*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void uncompress_anti_hermitian( anti_hermitmat *mat_antihermit,
	su2_matrix *mat_su2 ) {
/* uncompresses an anti_hermitian su2 matrix */
        Real temp1;
	mat_su2->e[0][0].imag=mat_antihermit->m00im;
	mat_su2->e[0][0].real=0.0;
	mat_su2->e[1][1].imag=mat_antihermit->m11im;
	mat_su2->e[1][1].real=0.0;
        mat_su2->e[0][1].imag=mat_antihermit->m01.imag;
	temp1=mat_antihermit->m01.real;
	mat_su2->e[0][1].real=temp1;
	mat_su2->e[1][0].real= -temp1;
	mat_su2->e[1][0].imag=mat_antihermit->m01.imag;
	
}



/*uncompress_anti_hermitian*/
