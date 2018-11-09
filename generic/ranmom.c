/*************************** ranmom.c *******************************/
/* Produce Gaussian random momenta for the gauge fields. */

#include "generic_includes.h"
//#include <defines.h>                 /* For SITERAND */


void ranmom() {
  register int i, j, mu;
  register site *s;
  double grn;

  FORALLSITES(i, s) {
    for(mu = XUP ; mu <=TUP ; mu++) {
      clear_su2mat(&(s->mom[mu]));
      for (j = 0; j < NUMGEN; j++) {
#ifdef SITERAND
        grn = (double)gaussian_rand_no(&(s->site_prn));
       
#else
        
        grn = (double)gaussian_rand_no(&node_prn);
        
#endif
        scalar_mult_add_su2_matrix(&(s->mom[mu]) ,&(Lambda[j]), grn, &(s->mom[mu]));
      }
    }
  }
}


