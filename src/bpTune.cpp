#include "stdafx.h"
#include <iostream>
#include "bptune.h"
#include "math.h"
#include <omp.h>
using namespace std;

bpTune::bpTune (int x, int y, int z, double a, double b, double c) {
   grid = false;
   SetGrid(x, y, z, a, b, c);
}

bpTune::~bpTune() {
   grid =false;
   //delete Potent;
}

void bpTune::SetGrid(int x, int y, int z, double a, double b, double c) {
   if (!grid) {
      Potent = new Mass<double>[x];
      FieldX = new Mass<double>[x];
      FieldY = new Mass<double>[x];
      FieldZ = new Mass<double>[x];
      for (int i = 0; i < x; i++) {
         Potent[i].set(y,z); 
         FieldX[i].set(y,z); 
         FieldY[i].set(y,z); 
         FieldZ[i].set(y,z); 
      }
      xyz[0] = x;
      xyz[1] = y;
      xyz[2] = z;
      abc[0] = a;
      abc[1] = b;
      abc[2] = c;
   }
   grid = true;
}

void bpTune::ResetGrid(double a, double b, double c) {
   abc[0] = a;
   abc[1] = b;
   abc[2] = c;
}

void bpTune::SetValue(double v) {
   if (grid) {
      for (int i = 0; i < xyz[0]; i++) 
      for (int j = 0; j < xyz[1]; j++)
      for (int k = 0; k < xyz[2]; k++) {
         Potent[i][j][k] = v;
         FieldX[i][j][k] = v;
         FieldY[i][j][k] = v;
         FieldZ[i][j][k] = v;
      }
   }
}

void bpTune::AddCharge(double*particle, double gammaRel) {
   for (int p = 0; p < 3; p++)
   if (fabs(particle[p*2]) > (xyz[p]-1) * abc[p] / 2.) {
      cout << "Particle outside of Potent: coordinate[" << p*2 << "] = " << particle[p*2] << "\n";
      return;
   }
   double diagonal2 = abc[0]*abc[0] + abc[1]*abc[1] + abc[2]*abc[2]*gammaRel*gammaRel;
   double x = abc[0] * (xyz[0]-1)/2 + particle[0];
   double y = abc[1] * (xyz[1]-1)/2 + particle[2];
   double z = abc[2] * (xyz[2]-1)/2 + particle[4];

#pragma omp parallel for
   for (int i = 0; i < xyz[0]; i++) 
   for (int j = 0; j < xyz[1]; j++)
   for (int k = 0; k < xyz[2]; k++) {
      double dist2 = pow(abc[0]*i-x,2) + pow(abc[1]*j-y,2) + pow(abc[2]*k-z,2);
      if (dist2 < diagonal2) dist2 = diagonal2;
      Potent[i][j][k] += 1. / sqrt(dist2);
      
      dist2 = pow(abc[0]*i-x,2);
      if (dist2 < abc[0]*abc[0]) dist2 = abc[0]*abc[0];
      if (abc[0]*i-x < 0) dist2 *= -1;
      FieldX[i][j][k] += 1. / dist2;
      
      dist2 = pow(abc[1]*i-y,2);
      if (dist2 < abc[1]*abc[1]) dist2 = abc[1]*abc[1];
      if (abc[1]*i-y < 0) dist2 *= -1;
      FieldY[i][j][k] += 1. / dist2;
      
      dist2 = pow(abc[2]*i-z,2);
      if (dist2 < abc[2]*abc[2]*gammaRel*gammaRel) dist2 = abc[2]*abc[2]*gammaRel*gammaRel;
      if (abc[2]*i-z < 0) dist2 *= -1;
      FieldZ[i][j][k] += 1. / dist2;
   }
}

void bpTune::SpaceCharge(double*particle, double K, double cell, bool benum) {

   int i  = floor(particle[0] / abc[0] + 0.5) + (xyz[0]-1) / 2;
   int j  = floor(particle[2] / abc[1] + 0.5) + (xyz[1]-1) / 2;
   int k  = floor(particle[4] / abc[2] + 0.5) + (xyz[2]-1) / 2;
   if (i < 0 || i >= xyz[0] || j < 0 || j >= xyz[1] || k < 0 || k >= xyz[2]) {
      cout << "Particle outside of Potent: i= " <<i<<", j= " <<j<< ", k= " <<k<< "\n";
      return;
   }

   double field[3];
   Field(field, i, j, k, cell, benum);

   particle[1] += (field[0] * K);
   particle[3] += (field[1] * K);
   particle[5] += (field[2] * K);
}

void bpTune::Field(double*field, int i, int j, int m, double cell, bool benum) {
   int dk;
   if (benum) dk = 100000;
       else   dk = floor(cell / abc[2] + 0.5); 
   for (int p = 0; p < 3; p++)
      field[p] = 0;

   for (int p = -10; p <= 10; p++) 
   {
      int k = p * dk + m;
      //int k = m;
      if ((k >= 0) && (k < xyz[2]-1)) {

/*         field[0] += ((Potent[i]  [j]  [k]   - Potent[i+1][j]  [k]) +
                      (Potent[i]  [j+1][k]   - Potent[i+1][j+1][k]) +
                      (Potent[i]  [j]  [k+1] - Potent[i+1][j]  [k+1]) +
                      (Potent[i]  [j+1][k+1] - Potent[i+1][j+1][k+1])) / (4. * abc[0]);

         field[1] += ((Potent[i]  [j]  [k]   - Potent[i]  [j+1][k]) +
                      (Potent[i+1][j]  [k]   - Potent[i+1][j+1][k]) +
                      (Potent[i]  [j]  [k+1] - Potent[i]  [j+1][k+1]) +
                      (Potent[i+1][j]  [k+1] - Potent[i+1][j+1][k+1])) / (4. * abc[1]);
    
         field[2] += ((Potent[i]  [j]  [k]   - Potent[i]  [j]  [k+1]) +
                      (Potent[i+1][j]  [k]   - Potent[i+1][j]  [k+1]) +
                      (Potent[i]  [j+1][k]   - Potent[i]  [j+1][k+1]) +
                      (Potent[i+1][j+1][k]   - Potent[i+1][j+1][k+1])) / (4. * abc[2]);*/
         field[0] += (Potent[i]  [j]  [k] - Potent[i+1][j]  [k])   / abc[0];
         field[1] += (Potent[i]  [j]  [k] - Potent[i]  [j+1][k])   / abc[1];
         field[2] += (Potent[i]  [j]  [k] - Potent[i]  [j]  [k+1]) / abc[2];
      }
   }
}
