//---------------------------------------------------------------------------
#ifndef xrestgasH
#define xrestgasH
#include "xeffect.h"
#include "xtarget.h"
//---------------------------------------------------------------------------

class xRestGas : public xEffect
{public:
   xRestGas();                                           // constructor
   xMaterial Comp[4];                                    // array of vacuum components
   bool ElectronCapture;
   bool SingleScattering;
   bool NuclearReaction;
   doubleU Pressure;                                     // pressure, Torr
   vectorU Rate0;
   int OnGet();
   int OnSet();
   vectorU Rates(xTime&, xBeam&, xRing&);               // calculates rates for residual gas effect
   void GetRate0(xBeam&,xRing&);
   void Kick (xTime&, xBeam&, xRing&);
};
extern xRestGas iRestGas;

class xPoint
{public:
   doubleU p;                                           // position [m]
   doubleU l;                                           // length [m]
   doubleU d;                                           // diameter [m]
   doubleU S;                                           // pumping speed [l/s]
   doubleU Sl;                                          // effective pumping left [l/s]
   doubleU Sr;                                          // effective pumping right [l/s]
   doubleU Sn;                                          // relative pumping speed [l/s]
   doubleU U;                                           // conductivity [l/s]
   doubleU Q;                                           // leak flow [l Torr/s]
   doubleU F;                                           // total flow [l Torr/s]
   doubleU G;                                           // outgasing [m Torr/s]
   doubleU P;                                           // vacuum pressure [Torr
   xPoint(double v = 0);
};

class xComponent : public LTemplate<xPoint>
{public:
   void Number(int n) { LSetSize(n); }
   int  Number(){return LGetSize( ); }
   void Plot(); 
   void StaticVacuum();

};
extern xComponent iComponent;


#endif
