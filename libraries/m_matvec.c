/****************  m_matvec.c  (in su2.a) *******************************
*									*
* void mult_su2_mat_vec( su2_matrix *a, su2_vector *b,*c )		*
* matrix times vector multiply, no adjoints 				*
*  C  <-  A*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"


void mult_su2_mat_vec( su2_matrix *a,vector *b,vector *c  ){
register int i,j;
register complex x,y;
    for(i=0;i<2;i++){
	x.real=x.imag=0.0;
	for(j=0;j<2;j++){
	    CMUL( a->e[i][j] , b->c[j] , y )
	    CSUM( x , y );
	}
	c->c[i] = x;
    }
}


