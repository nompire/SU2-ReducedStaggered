/*********************** rand_gauge.c ***************************/
/* original code by UMH */
/* 2/19/98 Version 5 port CD */

/* Makes a random gauge transformation on the gauge fields
   and then reunitarizes them */
/* Warning KS fermion applications: Must be done with KS phases OUT! */
/* Requires site workspace G */

#include "../generic/generic_includes.h"

void randomize();
void gauge_trans();

void rand_gauge()
	/* G holds the gauge transformation matrices */
{
    randomize();
    gauge_trans();
    node0_printf("Applying a random gauge transformation\n");     
}

void randomize()
{
  

  register int i;
  register site *s;
  register int m;
  FORALLSITES(i, s) {
#ifdef SITERAND
   
  
  for (m=0 ; m < NUMGEN ; m++){
  
  
    clear_su2mat(&(s->tmp));
    scalar_mult_add_su2_matrix(&(s->tmp),&(Lambda[m]), gaussian_rand_no(&(s->site_prn)) , &(s->tmp)); 
   }
   exp_su2_matrix(&(s->tmp),&(s->G));
   

#else
 
  for(m=0 ; m < NUMGEN ; m++){
  
    clear_su2mat(&(s->tmp));
    scalar_mult_add_su2_matrix(&(s->tmp), &(Lambda[m]) , gaussian_rand_no(&node_prn) , &(s->tmp));
  }
   exp_su2_matrix(&(s->tmp),&(s->G));
   

#endif
}
}

void gauge_trans()
{
  register int i,mu;
  site *s;
  su2_matrix temp;
  msg_tag *tag[4];

  FORALLUPDIR(mu) 
    tag[mu] = start_gather_site(F_OFFSET(G),sizeof(su2_matrix),mu,EVENANDODD,
		       gen_pt[mu]);

  FORALLUPDIR(mu) {
    wait_gather(tag[mu]);
    FORALLSITES(i,s) {

       mult_su2_an(&(s->G), &(s->link[mu]), &temp);
       mult_su2_nn(&temp, (su2_matrix *)gen_pt[mu][i],
		       &(s->link[mu]));

    }
    cleanup_gather(tag[mu]);
  }
}
