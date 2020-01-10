//---------------------------------------------------------------------------
#ifndef xEffectH
#define xEffectH
#include "bolideU.h"
#include "xBeam.h"
#include "xRing.h"
//---------------------------------------------------------------------------

class xEffect : public ATemplate<xEffect>, public xData
{
 public:
   chars EffectName;
   int  Multi;                 // Multi > 0 - increase step, Multi < 0 - decrease
   bool Loss;
   int ID;                                               // effect ID (label in List of effects)
   //static int Count;                                   // Counter for effects
   xEffect();
   virtual vectorU Rates  (xTime&, xBeam&, xRing&) = 0;  // calculates current effect's rates
   static  vectorU Summary(xTime&, xBeam&, xRing&);      // calculates summary rates from all effectS
   static xEffect* GetEffect(int);
   static xEffect* GetEffect(const char*); 

   xLattice Lattice;
   static void SetLattice(xLattice&);
   virtual void Kick (xTime&, xBeam&, xRing&) { ; }
};
//---------------------------------------------------------------------------

class xTaskRates : public xData
{
 public:
   bool show3D, E_log, dP_log;                               // markers for rates plots
   doubleU Eh_min, Eh_max, Ev_min, Ev_max, dP_min, dP_max;   // boundaries for rates plots
   int E_div, dP_div;                                        // divisions for rates plots
	xTaskRates();
	int OnGet();
	int OnRun();
};
//---------------------------------------------------------------------------

class xLosses : public xEffect
{
 public:
   bool Decay;
   bool Acceptance;
   bool Separatrix;
   bool Injection;


   xLosses();
   int OnGet();
   int OnRun();
   vectorU Rates(xTime&, xBeam&, xRing&);     // calculates effect's rates
   void Kick (xTime&, xBeam&, xRing&);
};
//---------------------------------------------------------------------------

class xHeat : public xEffect
{
 public:
   bool rate;
   doubleU Rate[3];
   double tapered;                          // tapered force
   bool rnum;

   bool linear;                             // switcher for linear type of heating
   doubleU Linear[3];                       // measure of linear heating of emittances and dP/P
   bool lnum;

   bool power;
   doubleU Power[3];
   bool pnum;

   bool diffusion;                           // switcher for diffusion type of heating
   doubleU Diffusion[3];                     // measure of diffusion heating of emittances and dP/P
   doubleU DiffusP;                          // measure of diffusion heating of dP/P (coasting beam)
   bool difP;                                // switcher for this tipe
   bool dnum;

   xHeat();
   int OnGet();
   int OnRun();
   vectorU Rates(xTime&, xBeam&, xRing&);    // calculates effect's rates
   void Kick (xTime&, xBeam&, xRing&);

   bool coherent;
   double k1;
   double k2;
   doubleU tau;
   doubleU position;
   doubleU length;
};
//---------------------------------------------------------------------------

class xColl : public xEffect
{
 public:
   doubleU Points;               //number of crossing points
   doubleU Cross;                //cross-section of the paticle loss
   doubleU Luminosity;
   doubleU Events;
   //bool EnableLoss;
   int LuminosityModel;
   bool Realtime;
   bool CollBeam;
   bool Local_BB;

   int Divisions;
   double From;
   double Upto;
   double Steps;
   int BB_EmitDef;
   double percent;
   doubleU Kappa_x;
   doubleU Kappa_y;

// 11.12.2007
   bool Include;           //Include diffusion
   doubleU Diffusion;      //Noice power

   xColl();
   int OnGet();
   int OnRun();
   int OnSet();
   double hourglass;
   double HourGlass();
   vectorU Rates(xTime&t, xBeam&, xRing&);               // calculates effect's rates
   void Kick (xTime&, xBeam&, xRing&);
   void CalcLuminosity(xTime&, xBeam&, xRing&);          // test of Luminosity calculation in realtime
   void BeamBeam(xTime& time, xBeam& beam, xRing& ring); // calculates beam-beam parameter
   // --- 05.06
   doubleU CollBeamDensity(double*);                     //local density of colliding beam
   doubleU MalakhovDensity(double*);                     //local density by Malakhov
};
//---------------------------------------------------------------------------

class xOptStoch : public xEffect
{
 public:
   doubleU R51;
   doubleU R52;
   doubleU R56;
   doubleU Go;
   doubleU Zo;
   doubleU SigmaZ;
   doubleU lambda;

   xOptStoch();
   int OnGet();
   int OnRun() { return OnGet(); }
   int OnSet();
   vectorU Rates(xTime&t, xBeam&b, xRing&r) { return vectorU(s_^-1); }
   void Kick (xTime&, xBeam&, xRing&);
};
//---------------------------------------------------------------------------

class xLaserCooling : public xEffect, public xOpticsLibrary
{
 public:
   doubleU ResonantFreq;
   doubleU LaserFreq;
   doubleU CoolLength;
   doubleU NuturalWidth;
   doubleU Saturation;
   doubleU SpotSize;
   doubleU Delta;
   doubleU S;
   double  Sweep0;
   double  Sweep;

   bool    SweepOn1;
   doubleU DeltaStart1;
   doubleU DeltaStop1;
   doubleU To;
   doubleU Shift1;
   doubleU Shift1fin;       

   bool    SweepOn2;
   doubleU DeltaStart2;
   doubleU DeltaStop2;
   int     SweepStep2;
   doubleU Shift2;
   doubleU Shift2fin;

   Tine<double> hyst;
   double TotalLight;

   xLaserCooling();
   int OnGet();
   int OnRun() { return OnGet(); }
   int OnSet();
   vectorU Rates(xTime&t, xBeam&b, xRing&r) { return vectorU(s_^-1); }
   void Kick (xTime&, xBeam&, xRing&);
   vectorU GetF      (xTime&t, U_Energy&e, vectorU ion);
   vectorU GetForces (xTime&t, U_Energy&e, vectorU ion);
};
//---------------------------------------------------------------------------

class xLastEffect : public xEffect
{
 public:
   xLastEffect(){ Multi = 1; }
   vectorU Rates(xTime&t, xBeam&b, xRing&r) { return vectorU(s_^-1); }
   void Kick (xTime&t, xBeam&b, xRing&r) { ; }
};
//---------------------------------------------------------------------------

extern xTaskRates iTaskRates;
extern xLosses iLosses;
extern xHeat iHeat;
extern xColl iColl;
extern xOptStoch iOptStoch;
extern xLaserCooling iLaserCooling;

#endif
