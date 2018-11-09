/*****************  make_ahmat.c  (in su2.a) ****************************
*									*
* void make_anti_hermitian( su2_matrix *m3, anti_hermitmat *ah3)	*
* take the traceless and anti_hermitian part of an su2 matrix 		*
* and compress it 							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

//#ifndef FAST
void make_anti_hermitian( su2_matrix *m3, anti_hermitmat *ah3 ) {
Real temp;
	
	temp = (m3->e[0][0].imag + m3->e[1][1].imag)*0.5;
	ah3->m00im = m3->e[0][0].imag - temp;
	ah3->m11im = m3->e[1][1].imag - temp;
	ah3->m01.real = (m3->e[0][1].real - m3->e[1][0].real)*0.5;
	ah3->m01.imag = (m3->e[0][1].imag + m3->e[1][0].imag)*0.5;
	

}/* make_anti_hermitian */


