//---------------------------------------------------------------------------
#ifndef xBeamH       
#define xBeamH
#include "bolideu.h"
#include "xring.h"
#include "xoptics.h"
#include "xdistributor.h"
//---------------------------------------------------------------------------
class xDraw;
class xTarget;

class xBeam : public xDistributor //xGener
{friend class xTarget;
 protected:
   U_Energy *pEnergy;
   xRing* pRing;
   doubleU *pCirc;   // ring circumference
   doubleU *pBetaH;  // ring horizontal beta-function
   doubleU *pBetaV;  // ring vertical beta-function
   doubleU *pTrev;   // ring revolution period
   doubleU *pEta;    // ring Off momentum compaction
   doubleU *pR;      // mean radius
   doubleU *pRF;     // RF system parameters
   BData Cross;

 public:
   double percent;
   double Nrms;
   int iMomentum;
   Units_ uMomentum;
   
   doubleU CellSize;      // number of particles in cell (for simulation)
   doubleU CellNum;       // number of cells
   int IniNumber;         // initial number of model particles
   int MaxNumber;         // maximum number of model particles
   doubleU Impact;        // minimum interparticle distance of collision
   int NImpact;           // number of impact events per turn;
   double rr, ss, Ir, Is, L;
   //doubleU Huge;
   //doubleU GammaPi;
   //doubleU Over;
   int Index;
   int GenerateOn;
   doubleU Ttrans;        //initial transverse temperature [K] - for MD
   doubleU Tlong;         //initial longitudinal temperature [K] - for MD
   doubleU Gamma1;        // gamma1 parameter - for MD
   doubleU Gamma2;        // gamma2 parameter - for MD
   doubleU Lambda;        // linear density - for MD
   //doubleU Tmin;
   //bool CanCooling;

   bool inject;           // inceject ions
   doubleU interval;      // injection interval
   int cycles;            // cycle number
   int ModelParticles;    // number of initial modeling particles
   vectorU InitialEmit;   // initial emittances

	xBeam(xRing&, int n = 1);

   void SpaceChargeND(xTime& time);          // space charge with Molecular Dynamics
   void SpaceChargeMD(xTime& time);          // space charge with MD (Hiroshi Tsutsui)
   void SpaceCharge2C(xTime& time);          // space charge for 2-cell approximation (A.Smirnov)
   void SpaceChargeBunch(xTime& time);       // space charge for bunch beam

   void Temperature();
   doubleU GetGammaPi(doubleU tlong);
   doubleU GetGamma2 (doubleU momentum);
   doubleU GetGamma3 (doubleU emit);
   doubleU GetEquillibrium(doubleU emit);
   void KeepCross();
   bool Crossed();
   bool Loss(xLattice&, int, double,bool);

// from Grisha
//---------------------------------------------------------------------------
protected:

   doubleU* pQ_x;       // pointer on ring horizontal tune
   doubleU* pKsy_x;     // pointer on ring chromaticity
   doubleU tz_0;        // koefficient 377 Ohm  - necessary to add to doubleU class **********


 public:

   //bool bunched;
   doubleU* pz_0; // pointer on tz_0 (= 377)
   doubleU* pRF_h; // pointer on harmonic number tRF_h
   doubleU* pRF_V; // pointer on RF Voltage tRF_V
   doubleU* pRF_Q_s; // pointer on Synchrotron tune tRF_Q_s
   doubleU* pRF_L_s; // pointer on Separatrix length tRF_L_s
   doubleU* pVac_a; // mean radius of vacuum chamber
   //20.11.06 For barrier bucket
   doubleU* pBB_V_B;                             //Voltage amplitude
   doubleU* pBB_T1;                               //RF duration in T0
   doubleU* pBB_T2;                               //Gap duration in T0
   doubleU* pBB_BarHeight;                       //Barrier height
   //
   //xLumBeam Lum;
   doubleU a; 				// beam radius
   doubleU F_sc;        // Image force correction factor
   doubleU F_l;         //Longitudinal distribution factor
   doubleU Z_l_i;
   doubleU Z_l;         //Longitudinal coupling impedance
   doubleU F_t;         //Transverse distribution factor
   doubleU Z_t_i;
   doubleU Z_t;         //Transverse coupling impedance
   doubleU G;           //Longitudinal form-factor
   doubleU Z_l_sc;      //Longitudinal space charge impedance
   doubleU Z_t_sc;      //Transverse space charge impedance
   doubleU I;           //Beam current, A
   doubleU D_Q_x;       //Tune shift horizontal
   doubleU D_Q_z;       //Tune shift vertical
   doubleU KS;          //Keil-Schnell parameter
   doubleU S_Q;         //Tune spread
   doubleU DM;          //Dipole mode parameter

   doubleU s_s;         //Rms bunch length, sm
   doubleU N_max;       //Maximum number of particles
   doubleU Q_s;         //Synchrotron  tune
   doubleU B_f;         //Bunching factor
   doubleU N_b;         //Number of bunches
   doubleU Ep;          // Longitudinal emittance
   //bool collider;       // switch for collider regime

//******************* functions *******************

   int CalcStability();    //calculates criterium for beam stability
   void Keil_Shneile();    //calculates Kei-Shneil criterium (beam instability)
   void Tune_shift();      // calculates betatron (X and Y) tune shifts
   void Shneil_Zotter();   //calculates Kei-Shneil criterium (beam instability)
   void CalcBunch();       // calculates bunch parameters: Synchrotron tune, Bunching factor, Number of particles in bunch
   //----------------20.11.06--
   double CalcBucket(doubleU dp2); //calculates the bucket length
   void SetBucket(vectorU& X);
   void GetBucket(vectorU& X, doubleU& T_int);
   doubleU Size_Bac;                        //Bucket length
   doubleU T_s;                             //Synchrotron period
   doubleU Amp2;                            //Square of momentum amplitude
   doubleU t_1;
   doubleU t_2;
   doubleU ksi2;

   doubleU Bdp2(doubleU,doubleU);
/*
   void Luminosity();        // calculates luminosity
   void SimCoolLum();        // calculate luminosity for arbitrary distribution
   doubleU Lcross(doubleU);  // calculates luminosity for crossed collision (with angle)
   doubleU Lhead(doubleU);   // calculates luminosity for head-on collision
   doubleU integral();       // calculates integral luminosity
*/
};

class xBeam1 : public xBeam, public xData
{
 public:
   xBeam1(xRing& ring, int n = 1) : xBeam(ring, n) { ; }
   int OnGet();
   int OnSet();
   int OnRun();
};

extern xBeam1 iBeam;


//---------------------------------------------------------------------------
// Colliding Beam

class CollidingBeam : public xData
{
   Units_ UKinetic;                         // keep Energy definition

 public:
   U_Energy Energy;                         // structure of U_Energy
   xLattice Lattice;                        // Lattice functions in collision points

   bool bunched;
   doubleU Sigma_x;
   doubleU Sigma_y;
   doubleU Sigma_s;
   vectorU Emit;     // vector of Emittances


   CollidingBeam();
   int OnGet();
   int OnSet();
   //int OnRun();
};

extern CollidingBeam cBeam;

#endif
