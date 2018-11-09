// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "../include/defines.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/su2.h"
#include "../include/random.h"    // For double_prn
#include "../include/dirs.h"      // For NDIMS
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

#ifdef SITERAND
  // The state information for a random number generator
  double_prn site_prn;
#endif

  // Staggered phases not including temporal boundary conditions
  Real phase[NDIMS];



 //Gauge field
 su2_matrix link[4];

#ifdef HMC_ALGORITHM
 su2_matrix old_link[4]; // For accept/reject
#endif

//Antihermitian momentum matrices

//anti_hermitmat mom[4];

su2_matrix mom[4];

//#ifdef CORR
//su2_matrix prop; //Propagators for fermion correlators
//#endif

//Link matrices for defining curly U_{mu}(x)
su2_matrix temp_link[4];

//Temporary matrices for use in computing fermion forces
su2_matrix f_U[4];


//Wilson staple
su2_matrix staple,tempmat1,tempmat2;

//Random gauge transformation
su2_matrix G, tmp;

} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

// Global variables
EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int PBC;             // Temporal fermion boundary condition
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int warms, trajecs, niter, propinterval;
EXTERN Real traj_length;



// Epsilon tensor and matrices to translate (mu, nu) to linear index
// (either anti-symmetric or self-dual)
//EXTERN Real perm[DIMF][DIMF][DIMF][DIMF];
//EXTERN int as_index[DIMF][DIMF];
//EXTERN int sd_index[DIMF][DIMF];

// More global parameters
EXTERN Real rsqmin, rsqprop;
EXTERN Real BETA;
EXTERN Real site_mass;
EXTERN double sigmasum;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag;     // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int saveflag;      // 1 if we will save the lattice;
EXTERN int total_iters;   // To be incremented by the multi-mass CG

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Stuff for multi-mass CG and RHMC
EXTERN int nsteps[2];           // Fermion and scalar steps
EXTERN Real ampdeg, *amp, *shift;
EXTERN Real ampdeg4, *amp4, *shift4;
EXTERN Real ampdeg8, *amp8, *shift8;
EXTERN int Nroot, Norder;
EXTERN Real gnorm, *fnorm, max_gf, *max_ff;

// Momenta and forces for the scalars

//EXTERN su2_matrix *temp_link1[4];
EXTERN su2_matrix *sigma[4];

EXTERN su2_matrix Lambda[NUMGEN];


// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

// Persistent fermions for matrix--vector operation
// Used in fermion_op and assemble_fermion_force
// Also used by correlator calculation
EXTERN vector *src;
EXTERN vector *dest;
//EXTERN vector **psol;

// Temporary vectors, matrices

EXTERN vector *tempvec;
EXTERN vector *tempvec1;
EXTERN vector *tempvec2;







EXTERN su2_matrix *tempmat;
//EXTERN su2_matrix *tmpmat2; // used in Polyakov loop calculation

//Structures defined for computing the generator matrices




EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// Probably need at least 8=2NDIMS
#define N_POINTERS 12
EXTERN char **gen_pt[N_POINTERS];


#ifdef CORR
EXTERN int Nstoch;                  // Number of stochastic sources
EXTERN int Nsrc;                    // Number of point sources
EXTERN int pnts[MAX_SRC][NDIMS];    // Point sources
EXTERN su2_matrix *prop;
EXTERN su2_matrix *prop2;
EXTERN su2_matrix *prop3;
#endif

#ifdef EIG
// Eigenvalue stuff
EXTERN int Nvec;
EXTERN double *eigVal;
EXTERN vector **eigVec;
EXTERN Real eig_tol;          // Tolerance for the eigenvalue computation
EXTERN int maxIter;           // Maximum iterations
#endif

#ifdef PFAFF
// Pfaffian stuff
EXTERN int Nmatvecs;                // For timing/counting
EXTERN int ckpt_load, ckpt_save;    // For checkpointing
#endif

#endif // _LATTICE_H
// -----------------------------------------------------------------
