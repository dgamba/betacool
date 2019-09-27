//---------------------------------------------------------------------------
#ifndef xDynamicH
#define xDynamicH
//#include <condefs.h>
#include "doubleu.h"
#include "xibs.h"
#include "xbeam.h"
#include "xring.h"
#include "xeffect.h"
#include "xrunge.h"
#include "xecool.h"
#include "xtarget.h"
#include "xrestgas.h"
#include "xstoch.h"
#include "xdraw.h"
#include "xbucket.h"
#include "pellets.h"

//---------------------------------------------------------------------------
class xDynamics : public xData
{public:

   xLoader Loader;

//Algorithm
   int Algoritm;          // type of the calculation model 0..3
   int cycles;
//Output
   int SavingInterval;      // interval of redrawing, sec
//Evolution
   int SkipPoints;         // number of every skipped points
   int CurveSize;          // size of the curve array
   bool AutoSkip;          // smitch for AutoSkip mode
   int SkipCount;          // counter for skipped points
// Dynamics
   doubleU dt;             // integration step for dynamics
   doubleU Maxdt;          // maximal allowed integration step
   double ifmore;          // step multiplier for variable step integration
   double ratio;           // growth tempo for beam rms parameters, %
// Model Beam
   doubleU IntegrationStep;// time step
   doubleU TurnNumber;      // turn number step
   double TimeTurn;
   int Synchrotron;        // Model of synchrotron motion
   bool LongSC;            // Calculate longitudinal space charge
   Tray<double>LongHyst;   // Longitudinal space charge distribution
   bool LongOK;
   int betatron;           // common betatron tune
   bool rotate;             // phase rotation in transverse plan
   doubleU angle;           // rotation angle per second
   doubleU Time_Stop;       // stop time
//Tracking
   doubleU ds;             // step over  optics element
   bool CheckCrossing;     // check naumber of crossing in longotudinal direction
   bool Crystall;          // use for crystall distribution
   bool FORCES;
//Monte-Carlo
/*
   int TurnStep25;         // integration step for dynamics
   bool TimeBox25;         // enabling variable step integration
   doubleU TimeStep25;     // maximal allowed integration step
   double MinInv25;        // growth tempo for invariants, %
*/
   xDynamics();             // xDynamics constructor
   int OnGet();
   int OnSet();
   int OnRun();
   void Distribution(vectorU&,xLattice&,int,bool);  // generates rms accordingly to Lattice (mean or current), matching of the beam
   bool Drawing(xLattice&);            // drawing for rms beam parameters (Evolution) and rate (Rates)

   void Dynamics();                    // rms beam dynamics
   void ModelBeam();                    // testbeam algorithm
	void Tracking();                    // beam tracking
   //void MonteCarlo();                  // Monte Carlo algorithm

   void Algorithms();                  // calls algorithm function in accordance to the chosen model
   void Lattice();                     // draws Lattices to Ring|Optic structure
   void Rates();                       // calls SummaryRates function and draws rates to Task|Rates
   void BeamTest();                    // test function for beam rms calculation
   void FFTest();                      // test function for Friction Force (ECOOL) calculation, draws 3D plots of Fr.force
   void LuminosityCalculation();       // test of Luminosity calculation in realtime
};
 
extern xDynamics iDynamics;

#endif
