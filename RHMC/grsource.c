// -----------------------------------------------------------------
// Construct a gaussian random vector R, return src = (Ddag.D)^{1 / 8} R
// Need to invert despite the positive power, since it is fractional
#include "su2_includes.h"

// Return the number of iterations from the inversion
int grsource(vector *src) {

  register int i, j;
  register site *s;
  int avs_iters;
  Real size_r;
  vector **psim = malloc(Norder * sizeof(**psim));

  // Allocate psim (will be zeroed in congrad_multi)
  for (i = 0; i < Norder; i++)
    psim[i] = malloc(sites_on_node * sizeof(vector));

  // Begin with pure gaussian random numbers
 
  
  FORALLSITES(i, s) {
#ifdef SITERAND

  
    src[i].c[0].real = 1.0/sqrt(2.0)*gaussian_rand_no(&(s->site_prn));
    src[i].c[0].imag = 1.0/sqrt(2.0)*gaussian_rand_no(&(s->site_prn));
    src[i].c[1].real = 1.0/sqrt(2.0)*gaussian_rand_no(&(s->site_prn));
    src[i].c[1].imag = 1.0/sqrt(2.0)*gaussian_rand_no(&(s->site_prn));
    
#else
     
    src[i].c[0].real = 1.0/sqrt(2.0)*gaussian_rand_no(&node_prn);
    src[i].c[0].imag = 1.0/sqrt(2.0)*gaussian_rand_no(&node_prn);
    src[i].c[1].real = 1.0/sqrt(2.0)*gaussian_rand_no(&node_prn);
    src[i].c[1].imag = 1.0/sqrt(2.0)*gaussian_rand_no(&node_prn);

#endif
  }


  double source_norm = 0.0;
  FORALLSITES(i, s) {

    source_norm += (double)magsq_vec(&(src[i]));
    
  }
  g_doublesum(&source_norm);
  
  node0_printf("source_norm in grsource %.4g\n", source_norm);


  // We now compute (Mdag M)^{1 / 8}.src
  for (i = 0; i < Norder; i++)
    shift[i] = shift8[i];
  
  avs_iters = congrad_multi(src, psim, niter, rsqmin, &size_r);

  

  node0_printf("Iters for source %d\n", avs_iters);

  // Reconstruct (Mdag M)^{1 / 8}.src from multi-mass CG solution psim
  FORALLSITES(i, s) {
    scalar_mult_vec(&(src[i]), ampdeg8, &(src[i]));
    for (j = 0; j < Norder; j++)
      scalar_mult_add_vec(&(src[i]), &(psim[j][i]), amp8[j], &(src[i]));
     
  }



  for (i = 0; i < Norder; i++)
    free(psim[i]);
  free(psim);
  return avs_iters;
}
// -----------------------------------------------------------------
