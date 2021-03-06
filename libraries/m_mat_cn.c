/******************  m_mat_an.c  (in su3.a) *****************************
*									*
* void mult_su2_cn( su2_matrix *a,*b,*c )				*
* matrix multiply, first matrix is adjoint 				*
* C <-  A_conjug*B							*
*/


#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void mult_su2_cn( su2_matrix *a, su2_matrix *b, su2_matrix *c ){
register int i,j,k;
register complex x,y;
    for(i=0;i<2;i++)for(j=0;j<2;j++){
	x.real=x.imag=0.0;
	for(k=0;k<2;k++){
	    CMULJ_( a->e[i][k] , b->e[k][j], y );
	    CSUM( x , y );
	}
	c->e[i][j] = x;
    }
}


