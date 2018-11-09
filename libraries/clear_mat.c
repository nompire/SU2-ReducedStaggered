/********************  clear_mat.c  (in su2.a) ********************
*
*void clear_su2mat( su2_matrix *dest )
*  clear an su2 matrix
* dest  <-  zero_matrix
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void clear_su2mat( su2_matrix *dest ){
register int i,j;
    for(i=0;i<2;i++)for(j=0;j<2;j++){
	dest->e[i][j].real = dest->e[i][j].imag = 0.0;
    }
}
