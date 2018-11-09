// -----------------------------------------------------------------
// Evaluate Polyakov loops in arbitrary direction
#include "../generic/generic_includes.h"
// -----------------------------------------------------------------
// -----------------------------------------------------------------
complex ploop(int dir) {
  register int i;
  register site *s;
  msg_tag *tag;
  complex sum, plp;
  su2_matrix *staple,*tempmat1;
  int t, N = nt;

  sum = cmplx(0, 0);
  
  staple = (su2_matrix *)malloc(sites_on_node*sizeof(su2_matrix));
  if(staple == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  tempmat1 = (su2_matrix *)malloc(sites_on_node*sizeof(su2_matrix));
  if(tempmat1 == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  switch(dir) {
    case XUP: N = nx; break;
    case YUP: N = ny; break;
    case ZUP: N = nz; break;
    case TUP: N = nt; break;
    default:
      node0_printf("ERROR: UNRECOGNIZED DIRECTION IN PLOOP()\n");
      terminate(1);
  }

  
  

  FORALLSITES(i,s){
  tempmat[i]=lattice[i].link[dir];
  }


  for (t = 1; t < N; t++) {
   
    tag = start_gather_field(tempmat, sizeof(su2_matrix),
                                    dir,EVENANDODD, gen_pt[0]);
    wait_gather(tag);
  
  FORALLSITES(i,s){

  if(s->t != 0) continue;
  
  if(t==1){
  
    mult_su2_nn(&(s->link[dir]),(su2_matrix *)gen_pt[0][i],&staple[i]);
    

  }
  else{
  
    mult_su2_nn(&staple[i],(su2_matrix *)gen_pt[0][i],&(tempmat[i]));
    staple[i] = tempmat[i];
  }
  }
  
  FORALLSITES(i,s){
   
   tempmat1[i] = *(su2_matrix *)(gen_pt[0][i]);
   tempmat[i] = tempmat1[i];
  
  }

    cleanup_gather(tag);
  }

  FORALLSITES(i,s){
    if(s->t != 0) continue;
    plp = trace_su2(&staple[i]);
    CSUM(sum, plp);   // Running complex sum
  }

  // Average all the loops we just calculated
  g_complexsum(&sum);
  plp.real = sum.real  /((Real) (nx*ny*nz));
  plp.imag = sum.imag  /((Real) (nx*ny*nz));
  free(tempmat1);
  free(staple);
  return(plp);
}
// -----------------------------------------------------------------
