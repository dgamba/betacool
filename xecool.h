//---------------------------------------------------------------------------
#ifndef xEcoolH
#define xEcoolH
//---------------------------------------------------------------------------
#include "bolideu.h"
#include "xbeam.h"
#include "xring.h"
#include "xeffect.h"
#include "xoptics.h"
#include "xebeam.h"
#include "xforce.h"
#include "xrunge.h"

//---------------------------------------------------------------------------

class xEcool : public xEffect, public xOpticsLibrary
{
 public:

   int E_model;                                  // model of friction force calculation
   int integr_step;                              // number of integration steps for Euler or RK

   int I_model;                                  // model of cooler (thin lens, Euler, R-K)
   int step_tr;                                  // number of integration steps in transverse direction
   int step_long;                                // number of integration steps in long direction (for bunch)
   int num_MC;                                   // number of particles for Monte Carlo
   //bool capture;                                 // switch of electron capture
   int SectionNumber;                            // ECOOL Section numbers
   bool FringeField;                             // frienge filed of solenoid
   int Lload;

   xEcool();
   int OnGet();
   int OnSet();
   vectorU GetForces (xTime&t, U_Energy&e, vectorU ion);              // calculates Forces (raght parts) for ECOOL effect
   vectorU Rates(xTime&t, xBeam&, xRing&);                            // calculates rates for ECOOL effect
   xLattice FinalLattice(xLattice& lat);                              // lattice function in the end of ECOOL system
   vectorU SingleParticle(xTime&t, xBeam& beam, xRing& ring);         // calculates ECOOL as single particle model
   doubleU New_Coord(xTime&tE, xBeam& beam, xRing& ring, vectorU&X);  // new coordinates for Monte Carlo & Kick
   vectorU MonteCarlo(xTime&t, xBeam& beam, xRing& ring);             // calculates ECOOL as Monte Carlo model
   doubleU LifeTime(xTime&, xBeam&, xRing&);                          // calculates beam lifetime accordingly to electron cooling influence

// these functions are now in xdistributor:
   static vectorU Courant_Snyder(vectorU, xLattice&, doubleU, bool);  // calculates Courant Snyder Invariants for current Lattice
//   vectorU CS_Coord2Emit(xBeam&, xLattice&, doubleU, bool);           // calculates beam emittances from the Invariants

   void Kick (xTime&, xBeam&, xRing&);
};
//---------------------------------------------------------------------------

extern xEcool iEcool; 

#endif
