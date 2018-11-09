/******************  rand_ahmat.c  (in su3.a) ***************************
*									*
* void random_anti_hermitian( anti_hermitmat *mat_antihermit, passthru *prn_pt)*
* Creates gaussian random anti-hermitian matrices			*
* Normalization is < |m01|^2 > = 1, or < m01.real*m01.real > = 1/2	*
* The argument "prn_pt" is a pointer to be passed to gaussian_rand_no() *
* RS6000 may choke on void *						*
*/
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"
#include "../include/su2.h"

void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt) {
Real r3;

        r3=gaussian_rand_no(prn_pt);
	//r8=gaussian_rand_no(prn_pt);
	mat_antihermit->m00im= (1.0 /sqrt(2.0)) * r3;
	mat_antihermit->m11im= -(1.0 /sqrt(2.0)) * r3;
	
	mat_antihermit->m01.real=(1.0 /sqrt(2.0)) * gaussian_rand_no(prn_pt);
	
	mat_antihermit->m01.imag=(1.0 /sqrt(2.0)) * gaussian_rand_no(prn_pt);
	
}


/*random_anti_hermitian_*/
