#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"
#include "../include/macros.h"
#include "../RHMC/su2_includes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void exp_su2_matrix(su2_matrix *u ,su2_matrix *v){

        su2_matrix *del = malloc(sizeof(*del));
        su2_matrix *prod = malloc(sizeof(*prod));
	double fac=1.0;
        int i=1;
        int j,k,l;
        for( j=0; j < DIMF ; j++){
        for( k=0; k < DIMF ; k++){
        prod->e[j][k] = cmplx(0.0,0.0);  
        v->e[j][k] =cmplx(0.0,0.0); 
              }
        }
 
        for(l=0 ; l < DIMF ; l++){
        prod->e[l][l] = cmplx(1.0,0.0);
        v->e[l][l] = cmplx(1.0,0.0);
        }
        
register int sum=0,counter=0;
	
        do{
        fac=fac*(double)i;
        mult_su2_nn(prod,u,prod);
        //prod=prod*u;
        scalar_mult_su2_matrix(prod, (1.0/fac),del);
        //del=prod*(1.0/fac);
        add_su2_matrix(v,del,v);
        //v=v+del;
        i++;}
        while(sqrt(realtrace_su2(del,del))>GAUGETOL);
        
        sum+=i;
        counter++;
     if(counter==100000){ 
     node0_printf("mean no. of terms in exp()=%.1f\n" ,(double)(sum/counter) ); counter=0;sum=0;}
     free(del);
     free(prod);
}
