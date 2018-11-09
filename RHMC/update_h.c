// -----------------------------------------------------------------
// Update the momentum matrices
#include "su2_includes.h"



// -----------------------------------------------------------------



// -----------------------------------------------------------------

// -----------------------------------------------------------------
double gauge_force(Real eps) {
  register int i, dir1, dir2;
  register site *st;
  register Real eb3 = eps;
  msg_tag *tag0, *tag1, *tag2;
  int start;
  su2_matrix tmat1, tmat2,tmat3,tmat4,tmat5;
  complex ctmp ;
  double norm = 0.0;
  


   



  // Loop over directions, update mom[dir1]
  for (dir1 = XUP; dir1 <= TUP; dir1++) {
    start = 1; // Indicates staple sum not initialized
    FORALLSITES(i, st)
      clear_su2mat(&(st->staple));    // Initialize staple

    // Loop over other directions
    // Compute force from plaquettes in the dir1, dir2 plane
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir2 != dir1) {
        // Get link[dir2] from direction dir1
        tag0 = start_gather_site(F_OFFSET(link[dir2]),
                                 sizeof(su2_matrix),
                                 dir1, EVENANDODD, gen_pt[0]);

        // Start gather for the "upper staple"
       tag2 = start_gather_site(F_OFFSET(link[dir1]),
                                sizeof(su2_matrix),
                                dir2, EVENANDODD, gen_pt[2]);

        // Begin the computation "at the dir2DOWN point"
        // We will later gather the intermediate result "to the home point"
        wait_gather(tag0);
        FORALLSITES(i, st) {
          mult_su2_an(&(st->link[dir2]), &(st->link[dir1]), &tmat1);
          mult_su2_nn(&tmat1, (su2_matrix *)gen_pt[0][i],
                      (su2_matrix *)&(st->tempmat1));
        }

        // Gather lower staple "up to home site"
        tag1 = start_gather_site(F_OFFSET(tempmat1), sizeof(su2_matrix),
                                 OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

        // The "upper" staple
        // One of the links has already been gathered,
        // since it was used in computing
        // the "lower" staple of the site above (in dir2)
        wait_gather(tag2);
        if (start) {  // This is the first contribution to staple
          FORALLSITES(i, st) {
            mult_su2_nn(&(st->link[dir2]), (su2_matrix *)gen_pt[2][i], &tmat1);
            mult_su2_na(&tmat1, (su2_matrix *)gen_pt[0][i], &(st->staple));
            scalar_mult_su2_matrix(&(st->staple), 1.0 ,&(st->staple));
          }
          start = 0;
        }
        else {
          FORALLSITES(i, st) {
            mult_su2_nn(&(st->link[dir2]), (su2_matrix *)gen_pt[2][i], &tmat1);
            mult_su2_na(&tmat1, (su2_matrix *)gen_pt[0][i], &tmat2);
            scalar_mult_add_su2_matrix(&(st->staple),&tmat2, 1.0 ,&(st->staple));
          }
        }

        wait_gather(tag1);
        FORALLSITES(i, st) {
          
          scalar_mult_add_su2_matrix(&(st->staple),(su2_matrix *)gen_pt[1][i], 1.0 ,&(st->staple));
        }
        cleanup_gather(tag0);
        cleanup_gather(tag1);
        cleanup_gather(tag2);
      }
    } // End of loop over other directions

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, st) {
      mult_su2_na(&(st->staple), &(st->link[dir1]) ,&tmat3);
      su2_adjoint(&(st->staple) ,&tmat5);
      mult_su2_nn(&(st->link[dir1]), &tmat5, &tmat1);
      
      scalar_mult_add_su2_matrix(&tmat1 ,&tmat3 , -1.0 , &(st->staple));
      unit_su2mat(&tmat4);
      scalar_mult_su2_matrix(&tmat4 , -0.5 ,&tmat4);
      ctmp =trace_su2(&(st->staple)) ;
      c_scalar_mult_add_su2mat(&(st->staple),&tmat4,&ctmp,&(st->staple));

      scalar_mult_su2_matrix(&(st->staple) , -1.0*(BETA/4.0) ,&tmat1);
      
      //add_su2_matrix(&(sigma[dir1][i]) , &tmat1 , &(sigma[dir1][i]));
      //Finally updating the momenta
     
      scalar_mult_add_su2_matrix(&(st->mom[dir1]), &tmat1, eb3, &(st->mom[dir1]));
      norm += (double)realtrace_su2(&tmat1, &tmat1);
    }
  } // End of loop over dir1

  g_doublesum(&norm);
  
  return (eb3 * sqrt(norm) / volume);
}
// -----------------------------------------------------------------


//Assume CG has been run and the solution is in psim[n]
double fermion_force(Real eps, vector *src, vector **sol){

Real ferm_epsilon = eps;
register site *st;
register int i,dir;
int n;
double norm =0.0;

//Zero the fermion force collectors

for (dir = XUP; dir <= TUP; dir++) {
     FORALLSITES(i, st) {
      
      clear_su2mat(&(sigma[dir][i]));
      
        }

}


//Hit the solution vectors with the fermion operator


for (n = 0; n < Norder; n++) {
//Hit the solution vector with M
fermion_op(sol[n], tempvec , PLUS); 

//Compute the fermion force

compute_fermion_force(tempvec,sol[n]);

for(dir = XUP ; dir <=TUP ; dir++){
FORALLSITES(i,st) { 

scalar_mult_add_su2_matrix(&(sigma[dir][i]),&(st->f_U[dir]),-1.0 * amp4[n], &(sigma[dir][i]));
//clear_su2mat(&(sigma[dir][i])); //Uncommenting this line sets the fermion force to zero

         }

    }

}

for(dir = XUP ; dir <=TUP ; dir++){

FORALLSITES(i,st){     
//Update the momenta with the fermion force 
scalar_mult_add_su2_matrix(&(st->mom[dir]), &sigma[dir][i],ferm_epsilon, &(st->mom[dir]));
      
norm += (double)realtrace_su2(&sigma[dir][i], &sigma[dir][i]);
     
    }

}


  
  g_doublesum(&norm);
  
  return (ferm_epsilon * sqrt(norm) / volume );
  
}
// -----------------------------------------------------------------

void compute_fermion_force(vector *psol,vector *sol){


 register int i, dir;
 register site *s;
 int dumcoords[NDIMS];
 double tr;
 complex dum;
 int c,d;
 vector tvec_psol, tvec_sol;
 su2_matrix tmat, tmat1,tmat2,tmat3,temp,temp1,tmp,tmp1,tmp2,tmp3;
 msg_tag *tag0[NDIMS] ,*tag1[NDIMS];
 Real link_mass =  site_mass;



for (dir = XUP; dir <= TUP; dir++) {
FORALLSITES(i,s) { 
         
       clear_su2mat(&(s->f_U[dir])); 
      
      } 

}

FORALLSITES(i,s){

vec_copy(&(psol[i]),&(src[i]));
vec_copy(&(sol[i]),&(dest[i]));

}



for (dir =XUP ;dir <=TUP ; dir++){
   FORALLSITES(i,s){
    clear_su2mat(&(s->temp_link[dir]));

    if(lattice[i].parity == EVEN){su2mat_copy(&(s->link[dir]),&(s->temp_link[dir])) ;}

    else{su2_conjug(&(s->link[dir]), &(s->temp_link[dir]));}


   }

}


for (dir = XUP; dir <= TUP; dir++) {

  //Start gathers for psol[n] = M sol[n]  , and sol[n]
   

   
   tag0[dir]= start_gather_field(dest, sizeof(vector), dir,
                                  EVENANDODD, gen_pt[dir]); 

   tag1[dir] = start_gather_field(src, sizeof(vector), dir,
                                  EVENANDODD, gen_pt[4+dir]); 

}


 

for (dir = XUP; dir <= TUP; dir++) {
wait_gather(tag0[dir]);
wait_gather(tag1[dir]);

FORALLSITES(i,s){

clear_su2mat(&tmp);clear_su2mat(&tmp2);
dumcoords[0]=s->x;
dumcoords[1]=s->y;
dumcoords[2]=s->z;
dumcoords[3]=s->t;

 vec_copy((vector *)gen_pt[dir][i], &tvec_sol);
      
 if (dir == TUP && PBC < 0 && s->t == nt - 1){  scalar_mult_vec(&tvec_sol, -1.0, &tvec_sol);}
 su2_projector(&tvec_sol, &(psol[i]),&temp);
 
  
  /*if(s->parity == EVEN){
  su2_trans(&temp,&tmat); 
  mult_su2_nn(&(s->temp_link[dir]),&tmat,&tmat1);
  add_matrix(&tmp,&tmat1);
  }
  else{
  su2_trans(&(s->temp_link[dir]),&tmat);
  mult_su2_nn(&temp,&tmat,&tmat2);
  sub_matrix(&tmp,&tmat2);
  }

  if(dumcoords[dir]%2==0){
  scalar_mult_su2_matrix(&tmp,0.5*(s->phase[dir] + link_mass*s->phase[dir]),&tmp1);
  }
  else{
  scalar_mult_su2_matrix(&tmp,0.5*(s->phase[dir] - link_mass*s->phase[dir]),&tmp1);
  }*/
  
 //Generator version for forces
 //-----------------------------------------------------------------
 for( c=0; c < NUMGEN; c++){
  
 if(s->parity == EVEN){
  mult_su2_nn(&(s->link[dir]),&temp,&tmat);
  mult_su2_nn(&(Lambda[c]),&tmat,&tmat1);
  dum = trace_su2(&tmat1);
  }
  else{ 
  mult_su2_cn(&(s->link[dir]),&temp,&tmat);
  mult_su2_cn(&(Lambda[c]), &tmat,&tmat1);
  dum =trace_su2(&tmat1);}
  
  if(dumcoords[dir]%2==0){
  tr = -0.5 * dum.real * (s->phase[dir] + link_mass * s->phase[dir]); 
  }
  else{
  tr = -0.5 * dum.real * (s->phase[dir] - link_mass * s->phase[dir]);
  }
  
  scalar_mult_add_su2_matrix(&tmp,&(Lambda[c]),tr,&tmp);
 }
 //-----------------------------------------------------------------
  vec_copy((vector *)gen_pt[4+dir][i], &tvec_psol);

  if (dir == TUP && PBC < 0 && s->t == nt - 1){scalar_mult_vec(&tvec_psol, -1.0, &tvec_psol); }
  su2_projector(&(sol[i]),&tvec_psol,&temp);

  /*if(s->parity == EVEN){
  mult_su2_nn(&(s->temp_link[dir]),&temp,&tmat1);
  add_matrix(&tmp2,&tmat1);
  }
  else{
  su2_trans(&(s->temp_link[dir]),&tmat);
  su2_trans(&temp,&tmat3);
  mult_su2_nn(&tmat3,&tmat,&tmat2);
  sub_matrix(&tmp2,&tmat2);
  }
 
  if(dumcoords[dir]%2==0){
  scalar_mult_su2_matrix(&tmp2,0.5*(s->phase[dir] + link_mass*s->phase[dir]),&tmp3);
  }
  else{
  scalar_mult_su2_matrix(&tmp2,0.5*(s->phase[dir] - link_mass*s->phase[dir]),&tmp3);
  }*/


  su2_trans(&temp,&temp1);
  for( d=0; d < NUMGEN; d++){
  if(s->parity == EVEN){
  mult_su2_nn(&(s->link[dir]),&temp1,&tmat);
  mult_su2_nn(&(Lambda[d]),&tmat,&tmat1);
  dum = trace_su2(&tmat1);}

  else{ 
  mult_su2_cn(&(s->link[dir]),&temp1,&tmat);
  mult_su2_cn(&(Lambda[d]),&tmat,&tmat1);
  dum = trace_su2(&tmat1);}
  

  if(dumcoords[dir]%2==0){
  tr = 0.5 * dum.real * (s->phase[dir] + link_mass * s->phase[dir]);
  }
 
  else{
  tr =  0.5 * dum.real * (s->phase[dir] - link_mass * s->phase[dir]);
  }
  



   scalar_mult_add_su2_matrix(&tmp,&(Lambda[d]),tr,&tmp);
   }
   


  /*sub_su2_matrix(&tmp1,&tmp3,&tmp);
  su2_adjoint(&tmp,&tmat);
  sub_su2_matrix(&tmp,&tmat,&temp); 
  dum=trace_su2(&temp);
  CMULREAL(dum,0.5,dum);
  unit_su2mat(&tmat);
  c_scalar_mult_sub_su2mat(&temp,&tmat,&dum,&tmp);
  */
  scalar_mult_su2_matrix(&tmp,2.0,&(s->f_U[dir]));
  
 }
  cleanup_gather(tag0[dir]);
  cleanup_gather(tag1[dir]);
} 

   


}
  

