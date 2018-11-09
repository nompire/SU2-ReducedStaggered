// -----------------------------------------------------------------
// Update the link matrices
// Go to eighth order in the exponential of the momentum matrices,
// since higher order integrators use large steps
// Evaluation is done as:
//   exp(H) * U = ( 1 + H + H^2/2 + H^3/6 ...) * U
//              = U + H*(U + (H/2)*(U + (H/3)*( ... )))

#include "su2_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_u(Real eps) {
  
  register int i, dir;
  register site *s;
  
  su2_matrix *link,temp1, temp2, htemp;
  register Real t2, t3, t4, t5, t6, t7, t8;
  
  //Take divisions out of site loop (can't be done by compiler)
  t2 = eps / 2.0;
  t3 = eps / 3.0;
  t4 = eps / 4.0;
  t5 = eps / 5.0;
  t6 = eps / 6.0;
  t7 = eps / 7.0;
  t8 = eps / 8.0;

  FORALLSITES(i, s) {
    for (dir = XUP; dir <= TUP; dir++) {
      


      su2mat_copy(&(s->mom[dir]), &htemp);
      
      link = &(s->link[dir]);

      /*scalar_mult_su2_matrix(&htemp,eps,&temp);
      exp_su2_matrix(&temp,&temp1);
      
      mult_su2_nn(&temp1,link,&temp2);
      su2mat_copy(&temp2,link); */

      mult_su2_nn(&htemp, link, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, t8, &temp2);

      mult_su2_nn(&htemp, &temp2, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, t7, &temp2);

      mult_su2_nn(&htemp, &temp2, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, t6, &temp2);

      mult_su2_nn(&htemp, &temp2, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, t5, &temp2);

      mult_su2_nn(&htemp, &temp2, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, t4, &temp2);

      mult_su2_nn(&htemp, &temp2, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, t3, &temp2);

      mult_su2_nn(&htemp, &temp2, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, t2, &temp2);

      mult_su2_nn(&htemp, &temp2, &temp1);
      scalar_mult_add_su2_matrix(link, &temp1, eps, &temp2);

      su2mat_copy(&temp2, link);
      
      //det_su2(link);
      
    }
  }
   
}
// -----------------------------------------------------------------
