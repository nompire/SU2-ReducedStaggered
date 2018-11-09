/******************  realtr.c  (in su2.a) *******************************
*									*
* Real realtrace_su2( su2_matrix *a,*b)				*
* return Re( Tr( A_adjoint*B )  					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

double realtrace_su2(  su2_matrix *a, su2_matrix *b ){
register int i,j;
register double sum =0.0;
    for(i=0;i<2;i++){
    for(j=0;j<2;j++){
	sum+= a->e[i][j].real*b->e[i][j].real + a->e[i][j].imag*b->e[i][j].imag;}


        }
return(sum);
}
