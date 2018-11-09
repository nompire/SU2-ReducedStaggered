// -----------------------------------------------------------------
// Main procedure for SU(2) gauge model with red. stagg'd fermions measurements
#define CONTROL
#include "su2_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, avm_iters=0;
  double dtime;
  complex plp;
  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
 
  setup_lambda();

  setup_rhmc();



  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  
  // This is the only local measurement for now
  

  // Main measurements
#ifdef CORR
  // Correlator measurements
  avm_iters=0;
  //avm_iters += condensates();
  //int j;
  //for (j = 0; j < Nsrc; j++) {
  //  node0_printf("Source point %d %d %d %d\n",
    //             pnts[j][0], pnts[j][1], pnts[j][2], pnts[j][3]);
  //  avm_iters += correlators(pnts[j]);
  //}


  plp = ploop(TUP);
  node0_printf("POLYAKOV LINE TUP %.8g %.8g\n",plp.real, plp.imag);

  plp = ploop(XUP);
  node0_printf("POLYAKOV LINE XUP %.8g %.8g\n",plp.real, plp.imag);

  plp = ploop(YUP);
  node0_printf("POLYAKOV LINE YUP %.8g %.8g\n",plp.real, plp.imag);

  plp = ploop(ZUP);
  node0_printf("POLYAKOV LINE ZUP %.8g %.8g\n",plp.real, plp.imag);
  
#endif

  node0_printf("RUNNING COMPLETED\n");
#ifdef CORR
  node0_printf("CG iters for measurements: %d\n", avm_iters);
  // total_iters is accumulated in the multiCG itself
  // Should equal total for correlator measurements
  node0_printf("total_iters = %d\n", total_iters);
#endif
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
