// -----------------------------------------------------------------
// Update lattice
// Leapfrog integrator
// Begin at "integral" time, with H and U evaluated at the same time

// Uncomment to print out debugging messages
#define UPDATE_DEBUG
#include "su2_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>         // For "finite"
#endif
// -----------------------------------------------------------------

// -----------------------------------------------------------------
int update_step(Real *fnorm, Real *gnorm, vector **src, vector ***psim) {
  int step, iters = 0, n;
  int STEP;
  Real final_rsq, eps = traj_length / (Real)nsteps[0];
  double tr;
   
  
  
  node0_printf("eps %.4g\n", eps);
 for(step = 0; step < (int) (0.5 * nsteps[0] ) ; step++){
  for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
      tr = fermion_force(0.5*eps, src[n], psim[n]); //Updates the gauge_momenta at the end of function
      
    fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
  }

 //Pure bosonic update where for each fermion update we have 
 for(STEP=0; STEP < nsteps[1] ; STEP ++) { 
 
   // Inner steps p(t) u(t)
    tr = gauge_force(0.5*(eps/nsteps[1])); // Updates the gauge_momenta at the end of function gauge_force(Real eps);
    *gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

    
  
   
   update_u(eps/nsteps[1]);

   
   tr = gauge_force(0.5*(eps/nsteps[1])); // Updates the gauge_momenta at the end of function gauge_force(Real eps);
    *gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;
   
  }
  
 //pure fermionic update 
  for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
      tr = fermion_force(0.5*eps, src[n], psim[n]); //Updates the gauge_momenta at the end of function
      
    fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
   }

}
     
  
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int j, n, iters = 0;
  Real final_rsq;
  
  static int first_time=1,accept=0,no_calls=0;
  double startaction, endaction, change;
  vector **src = malloc(Nroot * sizeof(**src));
  vector ***psim = malloc(Nroot * sizeof(***psim));

  for (n = 0; n < Nroot; n++) {
    src[n] = malloc(sites_on_node * sizeof(vector));
    psim[n] = malloc(Norder * sizeof(vector*));
    for (j = 0; j < Norder; j++)
      psim[n][j] = malloc(sites_on_node * sizeof(vector));
  }

  
  // Refresh the momenta
  ranmom();

  // Set up the fermion variables
  // Compute g and src = (Mdag M)^(1 / 8) g
  for (n = 0; n < Nroot; n++)
    iters += grsource(src[n]);

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (j = 0; j < Norder; j++)
    shift[j] = shift4[j];


#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_leapfrog -- original action\n");
#endif
  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
    
  
  node0_printf("No of iteration in CG = %.1d\n",iters);

  // Find initial action
  startaction = action(src, psim);
  gnorm = 0.0;
  max_gf = 0.0;
  for (n = 0; n < Nroot; n++) {
    fnorm[n] = 0.0;
    max_ff[n] = 0.0;
  }

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy gauge field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));
 
#endif
  // Do microcanonical updating
  

  
 iters += update_step(fnorm, &gnorm, src, psim);


  // Find ending action
  // Since update_step ended on a gauge update,
  // need to do conjugate gradient to get (Mdag M)^(-1 / 4) chi
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_leapfrog -- new action\n");
#endif
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
    
  endaction = action(src, psim);
  change = endaction - startaction;
#ifdef HMC_ALGORITHM
  // Reject configurations giving overflow
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change) > 1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting Apparent Overflow: Delta S = %.4g\n",
                 change);
    change = 1.0e20;
  }

  // Decide whether to accept, if not, copy old link field back
  // Careful -- must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  
  node0_printf("Exp(-change) = %.8g\n",exp(-change));
  
  node0_printf("deltaHovH = %.8g\n",fabs((double)(change/startaction)));
  
  if (exp(-change) < (double)xrandom) {
    if (traj_length > 0.0) {
      gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));
    }
    node0_printf("REJECT: delta S = %.4g start S = %.12g :end S = %.12g\n",
                 change, startaction, endaction);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
    accept++;
  }
#else
  // Only print check if not doing HMC
  node0_printf("CHECK: delta S = %.4g\n", (double)(change));
#endif // ifdef HMC

  for (n = 0; n < Nroot; n++) {
    free(src[n]);
    for (j = 0; j < Norder; j++)
      free(psim[n][j]);
    free(psim[n]);
  }
  free(src);
  free(psim);

  if (traj_length > 0) {
    node0_printf("MONITOR_FORCE_GAUGE   %.4g %.4g\n",
                 gnorm / (double)(2 * nsteps[0]), max_gf);
    for (n = 0; n < Nroot; n++) {
      node0_printf("MONITOR_FORCE_FERMION%d %.4g %.4g\n",
                   n, fnorm[n] / (double)(2 * nsteps[0]), max_ff[n]);
    }

    if((no_calls%10==0)&&(!first_time)){
    node0_printf("ACCEPTANCE RATE  %.4g\n",(double)accept/(double)no_calls );
    no_calls=0;
    accept=0;
   }
    first_time=0;
    no_calls++;
    return iters;
  }
  else
    return -99;
}
// -----------------------------------------------------------------
