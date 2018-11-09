/******************  su2_trans.c  (in su2.a) **************************
*									*
* void su2_trans( su2_matrix *a, su2_matrix *b )			*
* B  <- A_trans,  trans of an su2 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

/* adjoint of an su2 matrix */
void su2_trans( su2_matrix *a, su2_matrix *b ){
register int i,j;
    for(i=0;i<2;i++){
    for(j=0;j<2;j++){
        b->e[i][j].real = a->e[j][i].real;
        b->e[i][j].imag = a->e[j][i].imag;
    }
  }
}
