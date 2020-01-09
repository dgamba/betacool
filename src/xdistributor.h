//---------------------------------------------------------------------------
#include "bolideu.h"
#include "xbeam.h"
#include <time.h>

#ifndef xDistributorH
#define xDistributorH
//---------------------------------------------------------------------------
class xBunch
{public:

   // Global bunch in LRF
   Tray<double>Bunch;    // array of particles

   // Local bunch in PRF
   Tray<double>Local;    // local array
   Tine<int>Minst;       // minimim distance index array
   Tine<double>Dist2;    // distances of particles to ion
   double Mean[6];       // local vector of RMS: x, x', y, y', s, dp/p
   double Mean2[6];      // local squre of the RMS bunch dimensions

   //xBunch();
   void GetMean(vectorU ion, Tray<double>&bunch);
   void LRF_PRF(U_Energy& e, Tray<double>&local);
   void PRF_LRF(U_Energy& e, Tray<double>&local);

};

enum BunchEnum{COASTING, BUNCHED, BUCKET};

class xDistributor : public xBunch
{public:
   //Tray<double>Ions;
   Tray<double>Invs;

   xDistributor(int n = 0);
   void Destructor();
   void Number(int n);
   int  Number() { return Bunch.Number; }
   bool index(int);

   double* operator[] (int);
   vectorU operator() (int);
   void    operator() (int, vectorU);
   void    operator() (int, matrixU&);
   void    Matrix_2x2(int i, matrixU&);
   doubleU operator() (int, int);
   void    operator() (int, int, doubleU);

   bool PRF;
/*
   void PRF2LRF(U_Energy&);
   void LRF2PRF(U_Energy&);
*/
   void Copy(int, int);
   void Add (int, int);
   void Add (int, vectorU);
   void Add (int, int, doubleU);

   double* InvI(int);
   vectorU Inv(int);
   void    Inv(int, vectorU);
   doubleU Inv(int, int);
   void    Inv(int, int, doubleU);
//=============================================================================
   BunchEnum benum;
   virtual double CalcBucket(doubleU dp2){return 0;}
   virtual void SetBucket(vectorU& X){;}

   vectorU Sigma;        // vector of RMS: x, x', y, y', s, dp/p
   vectorU Sigma_2;      // Squre of the RMS bunch dimensions
   doubleU Sigma_s;      // Bunch length
   vectorU Emit;         // vector of Emittances
                         // [0] - horizontal emittance   [pi m rad]
                         // [1] - vertical emittance     [pi m rad]
                         // [2] - momentum spread
                         // [3] - particle number
                         // [4] - longitudinal emittance [pi mm rad]
                         // [5] - momentum deviation
   int EmitDef;          // emittance type presentation (RMS, Inv, FWHM, %)
   bool UseForIBS;
   bool MeanPercents;    // mean coordinates for Enclosed Pesents
   double Sigmas;        // number of sigma for distribution processing
   int Division;         // number of division for distribution processing
   double percentage;    // percentage of particles taken into account for invsriants or emittance
   double percentlong;   // percentage of particles for longitudinal degree of freedom
   Tray<double> Hyst3;
   //-------------------------------------------21.11.05------------------------
   Tine<double>CurrentAv; // Averaged profile for evolution
   bool HourGlass;
   int Array_FWHM;
   int X_FWHM;            // number of particles in FWHM X emittance
   int Y_FWHM;            // number of particles in FWHM Y emittance

   static double Flatten();      // generates Flattened distribution [-1..1]
   static double ran1(long *);// generates Flattened distribution [0..1]
   static double Gaussian();     // generates Gaussian distribution  [-3..3]
   static int Poisson(double);   // generates integer number at Poisson distribution
   int initial;                  // initial distribution
   double dpshift;
   double horshift;
   double vershift;
   double longshift;

   void RMS(vectorU emit, xLattice&);   // function calculates rms beam parameters
   void SetIon(int);                    // distribution for one ion
   void Matching(xLattice&,int);        // matching with lattice function
   void PlusDisp(xLattice&,int);        // matching with dispesrion
   void MinusDisp(xLattice&,int);       // matching without dispesrion
   void Distribution(vectorU emit, xLattice&, int n0=0);// distribution for beam

   void QuickSort (Tray<double>&array, int size, int col); // quick sort
   void SimpleSort(Tray<double>&array, int size, int col); // simple sort
   void Hystogra (double *hyst,  int div,  int col); // Invs hystorgram
   void Profile   (double *hyst,  int col, double emit, int div); // profile
   void Amplitudes(xLattice&);
   double Width(double* , int);            // calculates distribution width on half maximum

   void Courant_Snyder(xLattice&);              // calculates C-S Invariants for current Lattice
   void Courant_Snyder(xLattice&, int index);   // calculates C-S Invariants for current Lattice
   vectorU Emit_RMS(xLattice&);                 // calculates RMS emittances
   vectorU Emit_Tail(xLattice&,Tine<bool>&,int);// calculates RMS emittances in tail
   vectorU CS_Coord2Emit(xLattice&);            // calculates beam emittances from the Invariants
   vectorU EmitSort(double, double, xLattice&); // calculates emittance occupied by N(%) of particles
   vectorU Emittance_FWHM(xLattice&, int);      // calculates FWHM emittance
   vectorU Emit_FWHM(xLattice&);                // calculates FWHM emittance
   vectorU Emit_Def(xLattice&);                 // get emittance definition

   doubleU LocalDensity(double* X, int Count, xLattice&); // local density for coord.X inside Count of particles
   void RadialDensity(double*, int, double, double);      // radial density
   void ProfileDensity(double*, int, double**, double, double);    // profile to radial density
   //---05.06
   doubleU LocalDensity_G(double* X, xLattice&lat);      //local density for Gaussian beam
   doubleU LocalDensity_R(double* X, double part, xLattice&lat);

   double Fitting[3][5];
   Tray<double> xi;
   void Powel1(double *p, double *hys, double sig, double div);
   void Powel2(double *p, double *hys, double sig, double div);

   Tine<double> LProfile;                       // Longitidinal profile
   Tine<double> VProfile;                       // Potential profile
   Tine<double> Integral;                       // Integral of profile
   Tine<double> Momentum;                       // Momentum distribution
   Tine<double> Average;                        // Average momentum
   Tine<double> RatesDp;                        // Rates in slices
   int  index(int j, int num);
   void LongProfile();

};

#endif
