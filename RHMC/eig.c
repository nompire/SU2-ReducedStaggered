// -----------------------------------------------------------------
// Eigenvalue computation and helper functions
// !!!Path to primme.h may need local customization
#include "su2_includes.h"
#include "../PRIMME/PRIMMESRC/COMMONSRC/primme.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gaussian random vector
void rand_src(vector *src) {
  register int i;
  register site *s;

  // Begin with pure gaussian random numbers
  FORALLSITES(i, s) {
#ifdef SITERAND
    src[i].c[0].real = gaussian_rand_no(&(s->site_prn));
    src[i].c[0].imag = gaussian_rand_no(&(s->site_prn));
    src[i].c[1].real = gaussian_rand_no(&(s->site_prn));
    src[i].c[1].imag = gaussian_rand_no(&(s->site_prn));
#else
    src[i].c[0].real = gaussian_rand_no(&node_prn);
    src[i].c[0].imag = gaussian_rand_no(&node_prn);
    src[i].c[1].real = gaussian_rand_no(&node_prn);
    src[i].c[1].imag = gaussian_rand_no(&node_prn);
#endif
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get Nvec vectors (stored consecutively) and hit them by the matrix
void av_ov (void *x, void *y, int *Nvec, primme_params *primme) {
  int i, j, iter, ivec;
  Complex_Z *xx;

  for (ivec = 0; ivec < *Nvec; ivec++) {
    // Copy double-precision vector x into Real precision vector src
    // Each vector has DIMF components
    xx = ((Complex_Z*) x) + DIMF * ivec * sites_on_node;  // This vector in x
    iter = 0;
    for (i = 0; i < sites_on_node; i++) {
      for (j = 0; j < DIMF; j++) {
        src[i].c[j].real = xx[iter].r;
        src[i].c[j].imag = xx[iter].i;
        iter++;
      }
    }

#ifdef DEBUG_CHECK
    if (iter != DIMF * sites_on_node)
      printf("av_ov: iter = %d after source\n", iter);

    // Check that src is the same magnitude as x[ivec]
    register site *s;
    double xmag = 0.0, src_mag = 0.0;
    xx = ((Complex_Z*) x) + DIMF * ivec * sites_on_node;   // This vector in x
    for (i = 0; i < sites_on_node * DIMF; i++)
      xmag += xx[i].r * xx[i].r + xx[i].i * xx[i].i;
    FORALLSITES(i, s)
      src_mag += magsq_vec(&(src[i]));
    if (fabs(xmag - src_mag) > eig_tol * eig_tol) {
      node0_printf("av_ov: |x[%d]|^2 = %.4g but |src|^2 = %.4g (%.4g)\n",
                   ivec, xmag, src_mag, fabs(xmag - src_mag));
    }
#endif

#ifdef DEBUG_CHECK
    // Check that src is being copied appropriately
    node0_printf("eigVec[0] copy check:\n");
    dump_vec(&(src[0]));
#endif

    DSq(src, dest);    // D^2

    // Copy the resulting vector dest back to double-precision vector y
    // Each vector has DIMF components
    xx = ((Complex_Z*) y) + DIMF * ivec * sites_on_node;  // This vector in y
    iter = 0;
    for (i = 0; i < sites_on_node; i++) {
      for (j = 0; j < DIMF; j++) {
        xx[iter].r = (double)dest[i].c[j].real;
        xx[iter].i = (double)dest[i].c[j].imag;
        iter++;
      }
    }

#ifdef DEBUG_CHECK
    if (iter != DIMF * sites_on_node)
      printf("av_ov: iter = %d after source\n", iter);

    // Check that dest is the same magnitude as y[ivec]
    double ymag = 0.0, dest_mag = 0.0;
    xx = ((Complex_Z*) y) + DIMF * ivec * sites_on_node;   // This vector in x
    for (i = 0; i < sites_on_node * DIMF; i++)
      ymag += xx[i].r * xx[i].r + xx[i].i * xx[i].i;
    FORALLSITES(i, s)
      dest_mag += magsq_vec(&(dest[i]));
    if (fabs(ymag - dest_mag) > eig_tol * eig_tol) {
      node0_printf("av_ov: |y[%d]|^2 = %.4g but |dest|^2 = %.4g (%.4g)\n",
                   ivec, ymag, dest_mag, fabs(ymag - dest_mag));
    }
#endif
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function par_GlobalSumDouble is set as primme.globalSumDouble
void par_GlobalSumDouble(void *sendBuf, void *recvBuf,
                         int *count, primme_params *primme) {

  int i;
  for (i = 0; i < *count; i++)
    *((double*)recvBuf + i) = *((double*)sendBuf + i);

  g_vecdoublesum((double*)recvBuf, *count);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function make_evs computes eigenvalues through PRIMME
// Prints them with a quick check of |D^dag D phi - lambda phi|^2
// If flag==1 we calculate the smallest eigenvalues
// If flag==-1 we calculate the largest eigenvalues
int make_evs(int Nvec, vector **eigVec, double *eigVal, int flag) {
  register site* s;
  int i, j, ivec, iter = 0, ret, maxn = sites_on_node * DIMF;
  int total_iters = 0;
  double check, *rnorms = malloc(Nvec * sizeof(*rnorms)), dtime;
  Complex_Z *workVecs = malloc(Nvec * maxn * sizeof(*workVecs));
  static primme_params primme;
  vector tvec, *tmpVec = malloc(sites_on_node * sizeof(*tmpVec));

  // Check memory allocations
  if (workVecs == NULL) {
    node0_printf("ERROR in make_evs: couldn't allocate workVecs\n");
    exit(1);
  }
  if (rnorms == NULL) {
    node0_printf("ERROR in make_evs: couldn't allocate rnorms\n");
    exit(1);
  }

  // Initialize all the eigenvectors to random vectors
  for (ivec = 0; ivec < Nvec; ivec++) {
    eigVal[ivec] = 1e16;
    rand_src(eigVec[ivec]);
  }

  // Copy initial guesses into double-precision temporary fields
  // Each vector has DIMF components
  for (ivec = 0; ivec < Nvec; ivec++) {
    iter = DIMF * ivec * sites_on_node;   // This vector in workvecs
    for (i = 0; i < sites_on_node; i++) {
      for (j = 0; j < DIMF; j++) {
        workVecs[iter].r = eigVec[ivec][i].c[j].real;
        workVecs[iter].i = eigVec[ivec][i].c[j].imag;
        iter++;
      }
    }
  }
#ifdef DEBUG_CHECK
  if (iter != Nvec * maxn)
    printf("make_evs: iter = %d after input\n", iter);
#endif

  // Set the parameters of the EV finder
  primme_initialize(&primme);
  primme.n = maxn * number_of_nodes;              // Global size of matrix
  primme.nLocal = maxn;                           // Local volume
  primme.maxOuterIterations = maxIter;
  primme.maxMatvecs = maxIter + 5;
  primme.numProcs = number_of_nodes;
  primme.procID = this_node;
  primme.globalSumDouble = par_GlobalSumDouble;   // Wrapped function
  primme.matrixMatvec = av_ov;                    // Mat-vec wrapper

  primme_set_method(DEFAULT_MIN_MATVECS, &primme);
  primme.printLevel = 1;
  primme.eps = eig_tol;                   // Maximum residual
  primme.numEvals = Nvec;
  primme.initSize = 0;                    // Number of initial guesses
  if (flag == 1)
    primme.target = primme_smallest;
  else if (flag == -1)
    primme.target = primme_largest;
  else {
    node0_printf("make_evs: Unrecognized flag %d\n", flag);
    terminate(1);
  }
  primme.stats.numOuterIterations = 0;    // Just to make sure
//  primme_display_params(primme);

  // Call the actual EV finder and check return value
  dtime = -dclock();
  ret = zprimme(eigVal, workVecs, rnorms, &primme);
  iter = primme.stats.numOuterIterations;
  total_iters += iter;
  primme.stats.numOuterIterations = 0;
  while (ret != 0) {
    // Try again with looser residual
    primme.eps *= 10;
    dtime += dclock();
    node0_printf("%d iterations saturated in %.4g seconds, ", iter, dtime);
    node0_printf("loosening stopping condition to %.4g\n", primme.eps);

    dtime = -dclock();
    ret = zprimme(eigVal, workVecs, rnorms, &primme);
    iter = primme.stats.numOuterIterations;
    total_iters += iter;
    primme.stats.numOuterIterations = 0;
  }
  dtime += dclock();
  node0_printf("Converged after %d iterations in %.4g seconds\n", iter, dtime);

  // Copy double-precision temporary fields back into output
  // Each vector DIMF components
  for (ivec = 0; ivec < Nvec; ivec++) {
    iter = DIMF * ivec * sites_on_node;   // This vector in workvecs
    for (i = 0; i < sites_on_node; i++) {
      for (j = 0; j < DIMF; j++) {
        eigVec[ivec][i].c[j].real = workVecs[iter].r;
        eigVec[ivec][i].c[j].imag = workVecs[iter].i;
        iter++;
      }
    }
  }
#ifdef DEBUG_CHECK
  if (iter != Nvec * maxn)
    printf("make_evs: iter = %d after output\n", iter);
#endif

  // Print results and check |D^dag D phi - lambda phi|^2
  // !!! Assume Nvec has been increased by 2 to help quartet formation
  for (ivec = 0; ivec < Nvec - 2; ivec++) {
    check = 0.0;
    DSq(eigVec[ivec], tmpVec);
    FORALLSITES(i, s) {
      // tvec = tmpVec - eigVal[ivec] * eigVec[ivec]
      scalar_mult_add_vec(&(tmpVec[i]), &(eigVec[ivec][i]),
                                       -1.0 * eigVal[ivec], &tvec);
      check += magsq_vec(&tvec);
    }
    g_doublesum(&check);    // Accumulate across all nodes
    if (flag == 1)  {       // Braces suppress compiler warning
      node0_printf("EIGENVALUE %d %.8g %.8g\n", ivec, eigVal[ivec], check);
    }
    else if (flag == -1)
      node0_printf("BIGEIGVAL  %d %.8g %.8g\n", ivec, eigVal[ivec], check);
  }
  fflush(stdout);

  // Clean up
  free(workVecs);
  free(rnorms);
  primme_Free(&primme);
  return total_iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Check matrix elements and eigenvalues of <psi_j | D | psi_i>
// where the psi are eigenvectors of DDag.D
// Have checked that <psi_j | Ddag | psi_i> produces conjugate eigenvalues
void check_Dmat(int Nvec, vector **eigVec) {
  register int i;
  register site *s;
  char N = 'N';
  int ivec, jvec, stat = 0, unit = 1, doub = 2 * Nvec;
  double *store, *work, *eigs, *dum;
  complex tc, check;
  double real, imag;
  vector *tmpTF = malloc(sites_on_node * sizeof(*tmpTF));

  // Allocate double arrays expected by LAPACK
  store = malloc(2 * Nvec * Nvec * sizeof(*store));
  work = malloc(4 * Nvec * sizeof(*work));
  eigs = malloc(2 * Nvec * sizeof(*eigs));
  dum = malloc(2 * sizeof(*dum));

  // Hit each eigVec with D, then contract with every other eigVec
  for (ivec = 0; ivec < Nvec; ivec++) {
    fermion_op(eigVec[ivec], tmpTF, PLUS);

    for (jvec = 0; jvec < Nvec; jvec++) {
      check = cmplx(0.0, 0.0);
      FORALLSITES(i, s) {
        tc = cmplx_dot(&(eigVec[jvec][i]), &(tmpTF[i]));
        CSUM(check, tc);
        
      }
      real = check.real;
      imag = check.imag;
      g_doublesum(&real);
      g_doublesum(&imag);    // Accumulate across all nodes
//      node0_printf("D[%d, %d] (%.8g, %.4g)\n",
//                   ivec, jvec, check.real, check.imag);

      // Save in column-major double array expected by LAPACK
      store[2 * (jvec + Nvec * ivec)] = check.real;
      store[2 * (jvec + Nvec * ivec) + 1] = check.imag;
    }
  eigs[ivec] =0.0;
  }
  free(tmpTF);

  // Diagonalize <psi|D|psi> using LAPACK
  // Arguments summarized in susy_includes.h
  node0_printf("Using LAPACK to diagonalize <psi_j | D | psi_i>\n");
  zgeev_(&N, &N, &Nvec, store, &Nvec, eigs,
         dum, &unit, dum, &unit, work, &doub, work, &stat);

  // Print resulting eigenvalues
  for (ivec = 0; ivec < Nvec; ivec++)
    node0_printf("D_eig %d (%.6g, %.6g)\n",
                 ivec, eigs[2 * ivec], eigs[2 * ivec + 1]);

  // Free double arrays expected by LAPACK
  free(store);
  free(work);
  free(eigs);
  free(dum);
}
// -----------------------------------------------------------------
