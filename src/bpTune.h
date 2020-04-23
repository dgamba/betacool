#ifndef bpTuneH
#define bpTuneH
#include "bpData.h"

class bpTune
{
public:
   bool grid;
   Mass<double> *Potent;
   Mass<double> *FieldX;
   Mass<double> *FieldY;
   Mass<double> *FieldZ;
   int xyz[3];
   double abc[3];
   bpTune() { grid = false; }
   bpTune(int x, int y, int z, double a = 1, double b = 1, double c = 1);
   virtual ~bpTune();
   void SetGrid(int x, int y, int z, double a, double b, double c);
   void ResetGrid(double a, double b, double c);
   void SetValue(double v = 0);
   void AddCharge(double *particle, double gammaRel);
   //double** operator[] (int i) { return Grid[i]; }
   void SpaceCharge(double *particle, double K, double cell, bool benum);
   void Field(double *potential, int i, int j, int m, double cell, bool benum);
};

#endif
