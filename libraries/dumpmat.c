/******************  dumpmat.c  (in su2.a) ******************************
*									*
*  void dumpmat( su2_matrix *mat )					*
*  print out a 2x2 complex matrix					*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su2.h"

void dumpmat( su2_matrix *m ){
int i,j;
    for(i=0;i<2;i++){
	for(j=0;j<2;j++)printf("(%.1g,%.1g)\t",
	    m->e[i][j].real,m->e[i][j].imag);
	
    }
    printf("\n");
}
