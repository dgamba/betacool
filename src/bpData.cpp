#include "stdafx.h"
#include <iostream>
#include "bpData.h"
#include <cstring>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
const double NULLdata = -1e-300;

double Flatten()
{
   return 2. * rand() / RAND_MAX - 1.;
}

double Gaussian()
{
   double rand1 = double(rand()) / RAND_MAX;
   double rand2 = double(rand()) / RAND_MAX;
   if (rand1 < 1e-9)
      rand1 = 1e-9;
   double g = sqrt(-2.0 * log(rand1)) * sin(2.0 * 3.1415926 * rand2);
   return g;
}

int bpData::LoadFileSep(char *fn, char sep)
{
   FILE *fp;
   fp = fopen(fn, "rt");
   if (!fp)
   {
      cout << "Cannot read file : " << fn << "\n";
      return -1;
   }
   char lin[5000];
   char buf[500];
   int ro = 0;
   while ((fgets(lin, 5000, fp)) && (ro < row))
   {
      int i = 0;
      int co = 0;
      while ((lin[i] != '\n') && (lin[i] != '\0') && (co < col))
      {
         int j = 0;
         while ((lin[i] != sep) && (lin[i] != '\n') && (lin[i] != '\0'))
            buf[j++] = lin[i++];
         buf[j] = '\0';
         double d;
         char *endptr;
         d = strtod(buf, &endptr);
         if ((!strlen(buf)) || (*endptr != NULL))
            mass[ro][co] = NULLdata;
         else
            mass[ro][co] = d;
         if (lin[i] == sep)
            i++;
         co++;
      }
      ro++;
   }
   fclose(fp);
   return ro;
}

void bpData::SaveFileSep(char *fn, char *sep, int RN)
{
   char buf[5000];
   char val[50];
   FILE *fp;
   fp = fopen(fn, "wt");
   if (!fp)
   {
      cout << "Cannot save to file : " << fn << "\n";
      return;
   }
   if (!RN)
      RN = row;
   for (int r = 0; r < RN; r++)
   {
      strcpy(buf, "");
      for (int c = 0; c < col; c++)
      {
         if (mass[r][c] != NULLdata)
         {
            gcvt(mass[r][c], 9, val);
            strcat(buf, val);
            if (c < col - 1)
               strcat(buf, sep);
         }
         else
         {
            strcat(buf, "\t");
         }
      }
      strcat(buf, "\n");
      fputs(buf, fp);
   }
   fclose(fp);
}
//------------------------------------bpRing-----------------------------------
void bpRing::Init(bpBeam &beam, bool Vrf_OR_sigmab)
{
   Trev = Circ / beam.velocity;
   eta = 1. / beam.gamma / beam.gamma - Imagenery / GammaTr / GammaTr;
   lattice.betax = Circ / (M_2PI * Q_h);
   lattice.betay = Circ / (M_2PI * Q_h);
   lattice.Dx = lattice.betax / Q_h;
   lattice.Dy = lattice.betay / Q_v;
   L_s = Circ / Hrf;

   if (Vrf_OR_sigmab)
   { // define sigmab via Vrf
      Q_s = sqrt(Hrf * abs(eta) * beam.q0 * Vrf / (M_2PI * beam.m0 * 931500 * beam.gamma)) / beam.beta;
      B_s = Circ * abs(eta) / (M_2PI * Q_s);
      beam.emit.sigmab = B_s * sqrt(beam.emit.momentum2);
   }
   else
   { //define Vrf via sigmab
      B_s = beam.emit.sigmab / sqrt(beam.emit.momentum2);
      Q_s = Circ * abs(eta) / (M_2PI * B_s);
      Vrf = pow(Q_s * beam.beta, 2) * M_2PI * beam.m0 * 931500 * beam.gamma / (Hrf * abs(eta) * beam.q0);
   }
   if (beam.bunched != 1)
      beam.emit.sigmab = Circ;
}

//------------------------------------pbBeam-----------------------------------

void bpBeam::Init()
{
   gamma = 1. + (energy / 931.5);
   beta = sqrt(((pow(gamma, 2)) - 1) / pow(gamma, 2));
   velocity = beta * 299792458;
}

void bpBeam::MinusD(bpLattice &lat)
{
   for (int i = 0; i < col; i++)
   {
      mass[0][i] -= lat.Dx * mass[5][i];
      mass[1][i] -= lat.Dpx * mass[5][i];
      mass[2][i] -= lat.Dy * mass[5][i];
      mass[3][i] -= lat.Dpy * mass[5][i];
   }
}

void bpBeam::PlusD(bpLattice &lat)
{
   for (int i = 0; i < col; i++)
   {
      mass[0][i] += lat.Dx * mass[5][i];
      mass[1][i] += lat.Dpx * mass[5][i];
      mass[2][i] += lat.Dy * mass[5][i];
      mass[3][i] += lat.Dpy * mass[5][i];
   }
}

void bpBeam::RmsEmittance(bpLattice &lat)
{
   double MO[6];
   double SIG[6];
   double Corr[3];
   MinusD(lat);
   for (int j = 0; j < 6; j++)
   {
      MO[j] = 0;
      SIG[j] = 0;
   }
   for (int j = 0; j < 3; j++)
      Corr[j] = 0;
   for (int j = 0; j < 6; j++)
   {

      for (int i = 0; i < col; i++)
         MO[j] += mass[j][i];
      MO[j] /= col;
   }
   for (int j = 0; j < 6; j++)
   {
      for (int i = 0; i < col; i++)
         SIG[j] += pow(mass[j][i] - MO[j], 2);
      SIG[j] /= col;
   }
   for (int j = 0; j < 3; j++)
   {
      for (int i = 0; i < col; i++)
         Corr[j] += (mass[j * 2][i] - MO[j * 2]) * (mass[j * 2 + 1][i] - MO[j * 2 + 1]);
      Corr[j] /= col;
   }
   emit.emittancex = sqrt((SIG[0] * SIG[1]) - (Corr[0] * Corr[0]));
   emit.emittancey = sqrt((SIG[2] * SIG[3]) - (Corr[1] * Corr[1]));
   emit.emittancez = sqrt((SIG[4] * SIG[5]) - (Corr[2] * Corr[2]));
   emit.sigma_s = sqrt(SIG[4]);
   emit.momentum2 = SIG[5];
   PlusD(lat);
}

void bpBeam::Distribution(bpLattice &lat)
{
   lat.gammax = (((lat.alphax) * (lat.alphax)) + 1) / lat.betax;
   lat.gammay = (((lat.alphay) * (lat.alphay)) + 1) / lat.betay;

   Sigma[0] = sqrt(emit.emittancex * lat.betax);
   Sigma[1] = sqrt(emit.emittancex * lat.gammax);
   Sigma[2] = sqrt(emit.emittancey * lat.betay);
   Sigma[3] = sqrt(emit.emittancey * lat.gammay);
   Sigma[4] = emit.sigmab;
   Sigma[5] = sqrt(emit.momentum2);

   for (int i = 0; i < col; i++)
   {
      for (int j = 0; j < 6; j++)
      {
         mass[j][i] = Sigma[j];
         if ((j == 4) && (bunched == 0))
            mass[j][i] *= Flatten() / 2.;
         else
            mass[j][i] *= Gaussian();
      }
      mass[0][i] *= sqrt(1 / (lat.betax * lat.gammax));
      mass[2][i] *= sqrt(1 / (lat.betay * lat.gammay));
      mass[0][i] -= lat.alphax * lat.betax * mass[1][i] / (1 + lat.alphax * lat.alphax); // rotate on x coord by thin lens
      mass[2][i] -= lat.alphay * lat.betay * mass[3][i] / (1 + lat.alphay * lat.alphay); // rotate on y coord by thin lens
   }
   PlusD(lat);
}

void bpBeam::Invariants(bpLattice &lat, double B_s, int i, double Invs[3])
{
   lat.gammax = (1. + (lat.alphax * lat.alphax)) / lat.betax;
   lat.gammay = (1. + (lat.alphay * lat.alphay)) / lat.betay;

   double Xb = mass[0][i] - lat.Dx * mass[5][i];
   double Xs = mass[1][i] - lat.Dpx * mass[5][i];
   double Yb = mass[2][i] - lat.Dy * mass[5][i];
   double Ys = mass[3][i] - lat.Dpy * mass[5][i];

   Invs[0] = ((lat.betax * Xs * Xs) + (2. * lat.alphax * Xb * Xs) + (lat.gammax * Xb * Xb));
   Invs[1] = ((lat.betay * Ys * Ys) + (2. * lat.alphay * Yb * Ys) + (lat.gammay * Yb * Yb));
   Invs[2] = mass[5][i] * mass[5][i];
   if (bunched == 1)
      Invs[2] += mass[4][i] * mass[4][i] / (B_s * B_s);
}

void bpBeam::TransRotation(bpLattice &L1, bpLattice &L2, double mu_x, double mu_y, int i)
{
   double TM[4][4];

   TM[0][0] = sqrt(L2.betax / L1.betax) * (cos(mu_x) + L1.alphax * sin(mu_x));
   TM[1][1] = sqrt(L1.betax / L2.betax) * (cos(mu_x) - L2.alphax * sin(mu_x));
   TM[0][1] = sqrt(L1.betax * L2.betax) * sin(mu_x);
   TM[1][0] = -1. * ((L2.alphax - L1.alphax) * cos(mu_x) + (1. + (L2.alphax * L1.alphax)) * sin(mu_x)) / sqrt(L1.betax * L2.betax);

   TM[2][2] = sqrt(L2.betay / L1.betay) * (cos(mu_y) + L1.alphay * sin(mu_y));
   TM[3][3] = sqrt(L1.betay / L2.betay) * (cos(mu_y) - L2.alphay * sin(mu_y));
   TM[2][3] = sqrt(L1.betay * L2.betay) * sin(mu_y);
   TM[3][2] = -1. * ((L2.alphay - L1.alphay) * cos(mu_y) + (1. + (L2.alphay * L1.alphay)) * sin(mu_y)) / sqrt(L1.betay * L2.betay);

   //for (int i = 0; i < col; i++)
   {
      double ion[4];
      ion[0] = TM[0][0] * mass[0][i] + TM[0][1] * mass[1][i];
      ion[1] = TM[1][0] * mass[0][i] + TM[1][1] * mass[1][i];
      ion[2] = TM[2][2] * mass[2][i] + TM[2][3] * mass[3][i];
      ion[3] = TM[3][2] * mass[2][i] + TM[3][3] * mass[3][i];
      for (int j = 0; j < 4; j++)
         mass[j][i] = ion[j];
   }
}

void bpBeam::LongRotation(bpRing &ring, double mu, double dt)
{
   if (bunched == 0)
   {
      for (int j = 0; j < col; j++)
      {
         mass[4][j] += mass[5][j] * velocity * abs(ring.eta) * dt;
         mass[4][j] -= ring.Circ * floor(mass[4][j] / ring.Circ);
         mass[4][j] -= ring.Circ / 2;
      }
   }
   else if (bunched == 1)
   {
      for (int j = 0; j < col; j++)
      {
         mass[4][j] = cos(mu) * mass[4][j] + ring.B_s * sin(mu) * mass[5][j];
         mass[5][j] = cos(mu) * mass[5][j] - (((mass[4][j] - ring.B_s * sin(mu) * mass[5][j]) / cos(mu)) / ring.B_s) * sin(mu);
      }
   }
}

void bpBeam::operator+=(bpBeam &beam1)
{
   bpBeam beam0(col);
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < col; j++)
         beam0.mass[i][j] = mass[i][j];
   setsize(col + beam1.col);
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < beam0.col; j++)
         mass[i][j] = beam0.mass[i][j];
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < beam1.col; j++)
         mass[i][j + beam0.col] = beam1.mass[i][j];
}

void bpBeam::operator-=(int i)
{
   if (col > 2)
   {
      col--;
      if (i < col)
         for (int j = 0; j < 6; j++)
            mass[j][i] = mass[j][col];
   }
}
