//---------------------------------------------------------------------------
#ifndef xRungeH
#define xRungeH
#include "bolideU.h"
#include "doubleU.h"

//---------------------------------------------------------------------------
class xRing;

class xTime : public xData
{
public:
   doubleU t, dt, so, sr, ds; // parameters of the xTime structure
   doubleU *pVelocity;        // pointer on Ring.velocity
   doubleU *pCirc;            // pointer on Ring.circumeference
   double Turns;              // number of turns
   int Steps;
   xTime(xRing &);
   xTime(xTime &);
   int OnGet();
   void TimeStep(doubleU); // makes step over time
   void DistStep(doubleU); // makes step over ring
   void operator*=(double);
   void operator+=(double);
   void operator++();
   void operator=(xTime &);
};

typedef vectorU f_Methods(xTime &t, vectorU &X);

vectorU Kutta4(xTime &t, doubleU &h, vectorU X, f_Methods f1); // calculates integral with method of Runge-Kutta of 4t order
vectorU Euler(xTime &t, doubleU &h, vectorU X, f_Methods f1);  // calculates integral with method of Euler with variable step
vectorU Symple(xTime &t, doubleU &h, vectorU X, f_Methods f1); // calculates integral with method of simplectic integral
/*
class xRunge : public xData
{
 public:

   matrixU EigenV;                          // matrix of eigen vectors
   matrixU NormV;                           // matrix of normalized eigen vectors
   matrixU Canonic;                         // matrix of canonical eigen vectors
   vectorU Eigen;                           // eigen vector
   vectorU NormE;                           // normalized vector

// ***************  for coupled motion ***********************
   vectorU Squareq(vectorU);
   vectorU qubstep(complex<double> a);
   vectorU Qubeq  (vectorU);
   vectorU Quateq (vectorU);
   vectorU Koeff  (matrixU m4x4);
   bool QuickEigen (matrixU m4x4, double accuracy = 1e-8);
   void EigenValue (matrixU);
   void EigenVector(matrixU, bool);
   void NormVector();
   void SetCanonic(matrixU space, double ro);

//------------------------------------------------------

   typedef void SNU (matrixU& result, matrixU& variable, PVOID pvoid);
   typedef SNU* SNUP;
   bool SNEquation (SNUP function, matrixU& unknown, PVOID pvoid, double accuracy = 1e-8);
};
*/
extern xTime iTime;

#endif
