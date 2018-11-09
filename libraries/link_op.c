
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"


complex link_op(su2_matrix *a , su2_matrix *b){

register int i,j;
complex x,y;
x.real=0.0; x.imag=0.0;
for(i=0 ; i < DIMF ; i++){
for(j=0 ; j < DIMF ; j++){

   CMUL( a->e[i][j] , b->e[i][j] , y);
   
   CSUM(x,y);
      }

  }   
return (x);
}
