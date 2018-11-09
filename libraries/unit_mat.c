/********************  clear_mat.c  (in su2.a) ********************
*
*void unit_su2mat( su2_matrix *dest )
*  clear an su2 matrix
* dest  <-  unit_matrix
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void unit_su2mat( su2_matrix *dest ){
register int i,j;
    for(i=0;i<2;i++){
    for(j=0;j<2;j++){
    if( i != j){dest->e[i][j].real = 0.0;  dest->e[i][j].imag = 0.0;}
    else{dest->e[i][j].real = 1.0; dest->e[i][j].imag = 0.0;}
    }
 }

}
