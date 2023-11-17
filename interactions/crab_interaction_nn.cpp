/* These first two numbers describe the number of ident numbers used to
   screen phase space pt.s from files that are to be used in the generation
   of correlation functions. These codes may be different than those of the
    particles you are investigating. For instance, you may wish to do pi-pi-
    correlations and use the phase space pt.s for all 3 species of pions */
#define N1TYPES 1
#define N2TYPES 1
/* These are the identifications. See pdg.lbl.gov/rpp/mcdata/all.mc for a
   list of identification numbers */
const int IDENT1[N1TYPES]={2112};
const int IDENT2[N2TYPES]={2112};

#define MASS1 939.6
#define MASS2 939.6

/* Define if particles are identical */
#define IDENTICAL
/* Turn off and on the Coulomb Interaction */
//#define COULOMB
/* Turn off and on the Strong Interaction */
#define STRONG_INTERACTION

#define Q1Q2 1

#define INTERACTION_WSYM 0.25
#define INTERACTION_WANTI 0.75
#define INTERACTION_WNOSYM 0.0
/* fractions of symmetric and antisym weights of the various spin channels */

#define INTERACTION_DELK  4.0
/* spacing of mom. mesh for calc. of strong/coul. int. given in MeV/c
   and mom. is reduced momentum (1/2 Qinv for m1=m2) */
/* these are used to define both the strong and coulomb meshes */
#define INTERACTION_NKMAX 25
/* number of momentum points in mesh for strong/coul. interaction */

/* DEFINE THE STRONG INTERACTION POTENTIAL */
#ifdef STRONG_INTERACTION

/* Potentials can be found the file crab_potentials.cpp */
#define POTENTIAL vreid

/* These are definitions used for calculating wave func.s */
#define STRONG_NSPINS 2
#define STRONG_NETPARTIALS 2

/* these are used to define both the strong and coulomb meshes */
/* number of spatial points is mesh */
#define STRONG_NRMAX 150
/* The strong interaction is calculated for r<rmax.  delr=RMAX/NRMAX */
#define STRONG_RMAX 6.0
/* For calculating wave func.s when  r > RMAX, asymptotic forms are used */

/* ************************************************************ */

/* The Following are global variables used for the calculation of wave
   functions */

const int STRONG_NPARTIALS[STRONG_NSPINS]={1,1};
/* for various spin channels, this is the number of partial waves in
   each spin channel. These numbers should sum up to NETPARTIALS*/

const int STRONG_SYMM[STRONG_NSPINS]={1,-1};
/* parity required by spin choice, choose -1 for antisym, 1 for symm.,
   and 0 for none. for example, in pp, the s=1 waves would have negative,
   while the s=0 wave would have positive. for non-identical particles,
   one might choose zero. */

const double STRONG_WEIGHT[STRONG_NSPINS]={0.25,0.75};
/* weight of spin channels */

const int STRONG_L[STRONG_NETPARTIALS]={0,1};

/* ************************************************************ */

double vreid(double r,int ipart);

double vreid(double r,int ipart){
  double pmux,f1,f2,f4,f6,f7,vr;
  /* See the appendix of B.D. Day, PRC24, p. 1203 (1981).
     with Tensor forces neglected */
  if(ipart==0){
    /* S=0 */
    pmux=r*0.7;
    f1=exp(-pmux);
    f4=(f1*f1*f1*f1);
    f7=f4*(f1*f1*f1);
    vr=-10.463*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
  }
  if(ipart>0){
    /* S=1 */
    pmux=r*0.7;
    f1=exp(-pmux);
    f2=f1*f1;
    f4=f2*f2;
    f6=f4*f2;
    vr=((-10.463/3.0)*f1-933.48*f4+4152.1*f6)/pmux;
  }
  return vr;
}

#endif


