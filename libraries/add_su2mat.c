/********************  addmat.c (in su3.a)  *****************************
*									*
*  Add two SU2 matrices 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void add_matrix( su2_matrix *a, su2_matrix *b) {
register int i,j;
    for(i=0;i<2;i++)for(j=0;j<2;j++){
	CSUM( a->e[i][j], b->e[i][j]);
    }
}
