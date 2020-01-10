//---------------------------------------------------------------------------
#ifndef xbucketH
#define xbucketH
//#include "xBeam.h"
//---------------------------------------------------------------------------
class xBarrier
{
 public:

// input
   doubleU s1;          // initial coordinate
   doubleU s2;          // final coordinate
   doubleU Ub;          // barrier voltage
   doubleU Us;			// space charge potential
   xBarrier();

// calculated
   doubleU A;           // acceleration coefficient
   doubleU V;           // velocity coefficient
   doubleU H;           // barrier height [dP/P]
   doubleU UP;           // barrier potential
   void SetBarrier(double x1, double x2, double u);

   void SetS(vectorU& X, doubleU t);         // calculate coordinate
   doubleU GetT(vectorU& X, doubleU s);      // calculate time flight

   double D;          // particle density between barriers
   double N;          // particle number between barriers
   double P;          // momentum spread between barriers
   double M;          // average momentum
   double I;          // integral of particle number
};

class xBeam;

class xBucket
{
 public:

   doubleU SyncPeriod;                         // synchrotron period
   Tine<xBarrier> barrier;                     // series of barriers
   double Ucoeff;                              // potential scale coefficient
   double UPcoeff;                             // integral cossficient
   bool   bCoeff;
   bool IBSnorma;                              // IBS normalisation
   bool Show3D;
   double HourGlass;
   bool Periodical;

   xBucket();
   vectorU GetPeriod(vectorU X, doubleU dt);

   BData BarData;
   bool Analytic;
   bool Stationary;
   bool Moving;
   int SetBucket(BData& bdata);
   int SetMoving(doubleU t);
   void CalcParticle(xBeam&, xRing&);

   bool Harmonic1;
   int Number1;
   bool Harmonic2;
   int Number2;
   bool Induction;
   doubleU TimeStart;
   doubleU TimeFinish;
   complex<double>Impedance;
   bool LongSpaceCharge;

   //double xBucket::Density(xBeam& beam, xRing& ring, int j);
   int Pos(xBeam& beam, int j);
   double DroDs(xBeam& beam, int j);

   bool   cut;
   double vertex;
   Tray<double> LumiTray;
   double Luminosity(int i1, int i2, double beta, double ds);
   void LuminosityTest();
   void Generate();
};


#endif
