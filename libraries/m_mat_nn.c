/*******************  m_mat_nn.c  (in su2.a) ****************************
*									*
* void mult_su2_nn( su2_matrix *a,*b,*c )				*
* matrix multiply, no adjoints 						*
* C  <-  A*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

//#ifndef FAST
void mult_su2_nn( su2_matrix *a, su2_matrix *b, su2_matrix *c ){
register int i,j,k;
register complex x,y;
    for(i=0;i<2;i++)for(j=0;j<2;j++){
	x.real=x.imag=0.0;
	for(k=0;k<2;k++){
	    CMUL( a->e[i][k] , b->e[k][j] , y );
	    CSUM( x , y );
	}
	c->e[i][j] = x;
    }
}


