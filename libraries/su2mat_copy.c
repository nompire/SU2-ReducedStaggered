/*****************  su2mat_copy.c  (in su2.a) ***************************
*									*
* void su2mat_copy( su2_matrix *a, su2_matrix *b )			*
* Copy an su2 matrix:  B <- A   						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

/* Copy a su2 matrix:  b <- a   */
void su2mat_copy( su2_matrix *a, su2_matrix *b ){
register int i,j;
    for(i=0;i<2;i++){
     for(j=0;j<2;j++){
	b->e[i][j].real = a->e[i][j].real;
	b->e[i][j].imag = a->e[i][j].imag;
      }
   }
}
