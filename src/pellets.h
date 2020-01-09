#ifndef PelletsH
#define PelletsH
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "bolideU.h"


class pTime
{
 public:
	double t, s, dt, ds, dh;                   // parameters of the pTime structure
   double Velocity;                       // pointer on Ring.velocity
   int Steps;
   pTime();
   void TimeStep    (double, double);     // makes step over time
   void DistStep    (double, double);     // makes step over ring
   void operator += (double);
   void operator  = (pTime&);
};


class Coordinates
{
   double X[6];
 public:
   Coordinates();
   double&     operator [](int i);
   void        operator +=(Coordinates x0);
   Coordinates operator + (Coordinates f2);
   Coordinates operator * (double dt);
   Coordinates operator / (double h);
};


class Pellet
{ public:
   Coordinates position;
   double   ds, PxPy;
   void     setPosition(pTime);
   double   getP(double px [], int Np, double, double R);
   void     set(Coordinates, Coordinates);
};


class Flux
{
   Pellet *pellets;
 public:
   int Nplts;
   double L, Detector;
   Coordinates C0, C1;
   Pellet& operator [](int i);
   void set();
   void init(int);
   double getP(double[], int, double, double);
   double getPxPy(Pellet&, double px[], double ps[], int, double, double);
  ~Flux();
};


class xPellet : public xData
{  public:
   double PxPy, Average, step, tt, Rx, Rs;
   int    Nplts, num, N;
   pTime t;
   Coordinates X0, X1;
   Flux pellets;
   double getP( double[], double, double);
   int OnGet();
   void Run(double*, double*);
};
#endif
