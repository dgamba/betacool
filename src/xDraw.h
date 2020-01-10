//---------------------------------------------------------------------------
#ifndef xDrawH
#define xDrawH
#include "bolideU.h"
#include "xBeam.h"
#include "xRing.h"
#include "xEffect.h"
#include "xForce.h"
#include "xDynamic.h"
#include "xEcool.h"    
#include "pellets.h"
#include "xIBS.h"
#include "bpData.h"

#define ABCDD_
class xForce;
extern char* CalcTime;

class xDrawGraf
{public:
   xGraf Graf;
   int x_axis;
   int y_axis;
};

//---------------------------------------------------------------------------

class xDraw : public xData
{
 public:
	xDraw();
   int OnRun();
   int OnGet();

   xGraf Space;
   xGraf Barrier;
   xGraf LongProfile;
   xGraf GSystems;
   //xGraf ElectronShiftH[6];
   //xGraf ElectronShiftV[6];
   xTimer SpaceTimer;
   bool DrawRealSpace(xBeam&, xRing&, xLattice&, int);  // draws all beam real space plots

   double Sigmas;
   int Divisions;
   int EvolutionFor;
   int NormalisedOn;
   int EvolutionSlices;
   double EvolutionStep;
   //xGraf Distribution[9];
   xGraf Coordinate;
   xGraf Profile;
   xGraf Invariant;
   xSurf Evolution;
   bool  evolution;
   bool EvolutionAver;
   bool BiGaussian;
   bool SumBiGaussian;
   bool GFitD;
   bool GFitE;
   int cut;                              // 0-center, 1 - top (Effective lumminosity calc)  Krestnikov 14. 01. 09

   int MaximumPoints;                        // maximum number of points in curve (according to skipped)
   int CurrentNumber;                        // current number of point in curve
   int CrossNumber;
   xGraf TrajX;
   xGraf TrajY;
   xGraf Cross;
   xGraf Target;
   void SaveTrajectory();
   void ResetTrajectory(int, int);
   void Trajectory(xTime&, xBeam&, bool);

   xGraf Gamma2;
   xGraf Gamma3;
   void DrawGammaPi(xBeam&);
   void DrawGamma3(xBeam&);

   xGraf InvV;                                // plot for horizontal Invariant
   xGraf InvH;                                // plot for vertical Invariant
   xGraf InvdP;                               // plot for longitudinal Invariant

   xGraf BeamTxy;                            // plot for beam rms evolution: 3D plot emittance(momentum spread)
   xGraf BeamPrn;                            // footpront of invariants
   bool Emittances;
   bool Spread;
   xGraf Evolut;
   xGraf Emittance;
   xTimer EvolutionTimer;
   xGraf LaserLight;
   bool DrawBeamEvolution(xBeam&, xTime&, xColl&, int);  // draws all beam rms evolution plots

   xGraf Rates1;
   xTimer RatesTimer;
   bool BeamRates(vectorU, xTime&, int);     // draws all rate evolution plots

   xGraf Lattice;
   void LatticeFunctions(xRing&, xTime&);    // draws all lattice functions

   xSurf RateEh;                             // 3D plot for rates of horizontal emittance on momentum spread
   xSurf RateEv;                             // 3D plot for rates of vertical emittance on momentum spread
   xSurf RatedP;                             // 3D plot for rates of momentum spread on transverse emittance
   void Rates(xTaskRates&, xTime&, xBeam&, xRing&);   // draws 3D plots for rates

   xSurf FFtr;                               // 3D plot for transverse friction force of electron cooling effect
   xSurf FFlong;                             // 3D plot for longitudinal friction force of electron cooling effect
   xGraf FCtr;                              // plot transverse friction force alon line
   xGraf FClong;                             // plot for longitudinal friction force along line
   xGraf FAtr;                               //Force plots along angle
   xGraf FAlong;
   void FF(xForce&, xTime&, xEbeam&, xRing&);  // draws 3D plots for friction forces

   xGraf SpaceCharge;
   xGraf DriftVelocity;
   xGraf InEcool;
   void SpaceChargeDraw(xBeam& beam, xEbeam& ebeam);

   xGraf LuminosityTest;
   xSurf Luminosity3D;

   xGraf LaserForce;
   void DrawLaserForce(xTime&time, xRing&ring, xLaserCooling&lcool);

   xGraf EBunch;
   void DrawEBunch(xEbeam& ebeam);
   xSurf EDensity;                               // 3D plot for electron density calculated from array
   void DrawEDensity(xEbeam& ebeam);

//13.12.2007
   xGraf StochForce;
   xGraf StochDiff;
   void DrawStochForce(xTime&time, xRing&ring, xBeam& beam, xStochastic&scool);

   xGraf Pellet;
   xGraf Average;
   xGraf Pellets;
   xGraf BeamSize;
   xGraf BeamCros;
   double DrawPellet(double);
   void Oscilograph();

   xGraf Gain;
   void DrawGain();
   xGraf LStoch;
   void DrawLStoch();

   xSurf Potentials;
   xSurf ElectricField;
   void TuneShift();
};

extern xDraw iDraw;
#endif
