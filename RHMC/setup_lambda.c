#include "../include/complex.h"
#include "lattice.h"
#include "math.h"
#include "su2_includes.h"


#define IMAG_TOL 1.0e-8

#define DEBUC_CHECK

void setup_lambda() {
  int i, j, k, l, count;
  complex inv_sqrt = cmplx(1.0 / sqrt(2.0), 0.0);
  complex i_inv_sqrt = cmplx(0.0, 1.0 / sqrt(2.0));

  // Make sure Lambda matrices are initialized
  for (i = 0; i < NUMGEN; i++)
    clear_su2mat(&(Lambda[i]));

  // N * (N - 1) off-diagonal SU(N) generators
  // (T^{ij, +})_{kl} = i * (de_{ki} de_{lj} + de_{kj} de_{li}) / sqrt(2)
  // (T^{ij, -})_{kl} = (de_{ki} de_{lj} - de_{kj} de_{ki}) / sqrt(2)
  // Sign in second chosen to match previous values
  count = 0;
  for (i = 0; i < DIMF; i++) {
    for (j = i + 1; j < DIMF; j++) {
      for (k = 0; k < DIMF; k++) {
        for (l = 0; l < DIMF; l++) {
          if (k == i && l == j) {
            CSUM(Lambda[count].e[k][l], i_inv_sqrt);
            CSUM(Lambda[count + 1].e[k][l], inv_sqrt);
          }
          else if (k == j && l == i) {
            CSUM(Lambda[count].e[k][l], i_inv_sqrt);
            CDIF(Lambda[count + 1].e[k][l], inv_sqrt);
          }
        }
      }
      count += 2;
    }
  }
  if (count != DIMF * (DIMF - 1)) {
    node0_printf("ERROR: Wrong number of off-diagonal generators, ");
    node0_printf("%d vs. %d\n", count, DIMF * (DIMF - 1));
    terminate(1);
  }

  // N - 1 diagonal SU(N) generators
  // T^k = i * diag(1, 1, ..., -k, 0, ..., 0) / sqrt(k * (k + 1))
  for (i = 0; i < DIMF - 1; i++) {
    j = DIMF * (DIMF - 1) + i;    // Index after +/- above
    k = i + 1;
    i_inv_sqrt = cmplx(0.0, 1.0 / sqrt(k * (k + 1.0)));
    for (l = 0; l <= k; l++)
      Lambda[j].e[l][l] = i_inv_sqrt;
    CMULREAL(Lambda[j].e[k][k], -1.0 * k, Lambda[j].e[k][k]);
  }





  
  // Print Lambdas
  node0_printf("Computing generators for SU(N)\n");
  for (i = 0; i < NUMGEN; i++){
    node0_printf("Lambda[%d]\n",i);
    if (this_node == 0)
      dumpmat(&(Lambda[i]));
  }


  int a;
  complex tc;
  su2_matrix tmat;
  // Test group theory (useful reference: arXiv:1310.5353)
  node0_printf("Check group theory ");
  node0_printf("Sum_a Lambda^a_{kl} Lambda^a_{ij} = -delta_kj delta_il\n");
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      for (k = 0; k < DIMF; k++) {
        for (l = 0; l < DIMF; l++) {
          tc.real = Lambda[0].e[k][l].real * Lambda[0].e[i][j].real
                  - Lambda[0].e[k][l].imag * Lambda[0].e[i][j].imag;
          tc.imag = Lambda[0].e[k][l].imag * Lambda[0].e[i][j].real
                  + Lambda[0].e[k][l].real * Lambda[0].e[i][j].imag;
          for (a = 1; a < DIMF; a++) {
            tc.real += Lambda[a].e[k][l].real * Lambda[a].e[i][j].real
                     - Lambda[a].e[k][l].imag * Lambda[a].e[i][j].imag;
            tc.imag += Lambda[a].e[k][l].imag * Lambda[a].e[i][j].real
                     + Lambda[a].e[k][l].real * Lambda[a].e[i][j].imag;
          }
          if (cabs_sq(&tc) > IMAG_TOL)
            node0_printf("Sum_a La^a_{%d%d} La^a_{%d%d} = (%.4g, %.4g)\n",
                         k, l, i, j, tc.real, tc.imag);
        }
      }
    }
  }

  // Test orthogonality of products of Lambdas
  for (i = 0; i < NUMGEN; i++) {
    for (j = 0; j < NUMGEN; j++) {
      mult_su2_nn(&(Lambda[i]), &(Lambda[j]), &tmat);
      tc = trace_su2(&tmat);
      if (tc.real * tc.real > IMAG_TOL)
        node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n",
                     i, j, tc.real, tc.imag);
    }
  }

}
