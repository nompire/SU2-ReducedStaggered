/******************  su3_dot.c  (in su3.a) ******************************
*									*
* complex cmplx_dot(vector *a, vector *b )			*
* return dot product of two su2_vectors: a^dagger b			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

complex cmplx_dot(vector *a,vector *b ){


complex temp1,temp2;
    CMULJ_(a->c[0],b->c[0],temp1) 
    CMULJ_(a->c[1],b->c[1],temp2)
    CSUM(temp1,temp2);
    return(temp1);


}
