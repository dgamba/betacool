//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xRestgas.h"
//---------------------------------------------------------------------------

xRestGas::xRestGas()
{
   ID = 2;
   EffectName = "XRESTGAS";

   Pressure._(Torr_);
   Rate0[0]._(m_);
   Rate0[1]._(m_);
   Rate0[3]._(s_^-1);
}

int xRestGas::OnGet()
{
   Pressure = Data[43][2];

   Comp[0].Percent = Data.OnGet(43,3,100);
   Comp[1].Percent = Data.OnGet(43,4,0);
   Comp[2].Percent = Data.OnGet(43,5,0);
   Comp[3].Percent = Data.OnGet(43,6,0);

   Comp[0].A = Data.OnGet(43,7,1);
   Comp[1].A = Data.OnGet(43,8,1);
   Comp[2].A = Data.OnGet(43,9,1);
   Comp[3].A = Data.OnGet(43,10,1);

   Comp[0].Z = Data.OnGet(43,11,1);
   Comp[1].Z = Data.OnGet(43,12,1);
   Comp[2].Z = Data.OnGet(43,13,1);
   Comp[3].Z = Data.OnGet(43,14,1);

   ElectronCapture   = Data.OnGetB(43, 15, true);
   SingleScattering  = Data.OnGetB(43, 16, true);
   NuclearReaction   = Data.OnGetB(43, 17, true);

   for (int i = 0; i < 4; i++)
   {
      Comp[i].ElectronCapture = ElectronCapture;
      Comp[i].SingleScattering = SingleScattering;
      Comp[i].NuclearReaction = NuclearReaction;
   }

   return 0;
}

int xRestGas::OnSet()
{
   doubleU V_st(22400, (cm_^3));
   doubleU den(cm_^-3);
   den = Pressure()*U_NA/(V_st*760);
   for (int j = 0; j < 4; j++)                   
      Data.OnSet(43,j+18,den(cm_^-3)*Comp[j].Percent/100);
   return 0;
}

void xRestGas::GetRate0(xBeam&beam, xRing& ring)
{
   xLattice Lat;
   doubleU V_st(22400, (cm_^3));

   for (int i = 0; i < ring.Number(); i++)
   {  Lat = ring.GetLattice(i);
      for (int j = 0; j < 4; j++)
      if(Comp[j].Percent > 0 && ring[i].Length() > 0) // -- Kyoto --
      {  Comp[j].l = ring[i].Length;
         Comp[j].N = (Pressure()*U_NA) /(V_st*760);
         Comp[j].Nx = Comp[j].N * Comp[j].l;
         Comp[j].Ro = Comp[j].N*Comp[j].A/(U_NA*Comp[j].gr);
         Comp[j].Rox = Comp[j].Ro*Comp[j].l;
         Comp[j].RMSANGLE(ring.Energy);
         Comp[j].Bethe(ring.Energy);
         Comp[j].Deviations(ring.Energy, Lat);

         Rate0[0] += Comp[j].dEhor * Comp[j].Percent;
         Rate0[1] += Comp[j].dEver * Comp[j].Percent;
         Rate0[2] +=(Comp[j].dPP2  - Comp[j].dPB2) * Comp[j].Percent;
         if (Loss)
         Rate0[3] -= Comp[j].Lifetime(ring, beam, Lat)  * Comp[j].Percent;
      }
   }
   Rate0[0] /= 100 * ring.Trev;
   Rate0[1] /= 100 * ring.Trev;
   Rate0[2] /= 100 * ring.Trev;
   Rate0[3] /= 100;
}

vectorU xRestGas::Rates(xTime&time, xBeam&beam, xRing&ring)
{
   vectorU rates(s_^-1);
/*
   vectorU rates_tmp(s_^-1);
   U_Energy& Energy = ring.Energy;
   xLattice Lat;
   doubleU V_st(22400, (cm_^3));
   int kol;
   kol = 0;

    for (int j = 0; j < 4; j++)
		  if (Comp[j].Percent != 0)  kol++;

   for (int i = 0; i < ring.Number(); i++)
   {
      Lat = ring.GetLattice(i);
    for (int j = 0; j < kol; j++)
      {
       Comp[j].l = ring[i].Length;
       Comp[j].N = (Pressure()*U_NA) /(V_st*760);
       Comp[j].Nx = Comp[j].N * Comp[j].l;
       Comp[j].Ro = Comp[j].N*Comp[j].A/(U_NA*Comp[j].gr);
       Comp[j].Rox = Comp[j].Ro*Comp[j].l;
       Comp[j].RMSANGLE(Energy);
       Comp[j].Bethe(Energy);
       Comp[j].Deviations(Energy, Lat);

      rates[0] += Comp[j].dEhor * Comp[j].Percent / (100*beam.Emit[0]*ring.Trev);
      rates[1] += Comp[j].dEver * Comp[j].Percent / (100*beam.Emit[1]*ring.Trev);
      rates[2] += Comp[j].dPP2 * Comp[j].Percent  / (100*beam.Emit[2]*ring.Trev);
      rates[2] -= Comp[j].dPB2 * Comp[j].Percent /  (100*beam.Emit[2]*ring.Trev);
      rates[3] -= Comp[j].Lifetime(ring, beam)  * Comp[j].Percent / 100;
     }
   }
*/
   static bool once = true;
   if(once)
   {  once = false;
      GetRate0(beam, ring);
   }

   rates[0] = Rate0[0] / beam.Emit[0];
   rates[1] = Rate0[1] / beam.Emit[1];
   rates[2] = Rate0[2] / beam.Emit[2];
   rates[3] = Rate0[3];

   return rates;
}

void xRestGas::Kick (xTime&time, xBeam&beam, xRing&ring)
{
   beam.RMS(beam.Emit,Lattice);
   vectorU rates;
   rates = Rates(time, beam, ring);
   double loss = 0;
   if (Loss)
      loss = (-time.dt * rates[3])(U1_);

   for (int j = 0; j < beam.Number(); j++)
   {  if (!beam.Loss(Lattice, j, loss,false))
      {
         beam.Add(j,1,(U_Abs((beam.Emit[0] / Lattice.betax) *
         rates[0] * time.dt * 2.0)^0.5)*xDistributor::Gaussian());
         beam.Add(j,3,(U_Abs((beam.Emit[1] / Lattice.betay) *
         rates[1] * time.dt * 2.0)^0.5)*xDistributor::Gaussian());
         beam.Add(j,5,(U_Abs( beam.Emit[2]                  *
         rates[2] * time.dt * 2.0)^0.5)*xDistributor::Gaussian());
      }
   }
}

//---------------------------------------------------------------------------
xPoint::xPoint(double v)
{
   p = 0;
   l = 0;
   d = 0;
   S = 0;
   Sl= 0;
   Sr= 0;
   Sn= 0;
   U = 0;
   Q = 0;
   F = 0;
   G = 0;
   P = 0;
}

void xComponent::Plot()
{
   Number(10);
   xComponent& c = *this;

   for (int i = 0; i < Number(); i++)
   {
      c[i].p = i;
      c[i].l = 1;
      c[i].d = 1;
      //c[i].Q = 1;
   }
      c[2].S = 1;
      //c[9].S = 1;
      c[5].Q = 1;

   StaticVacuum();

   xGraf graf;
   graf.SetName("vac.cur");
   graf.Size(102,3);
   double X[3];                      

   c[0].F = c[Number()-1].F;
   for (int i = 0; i < Number(); i++)
   {
      X[0] = c[i].p();
      X[1] = c[i].F();
      X[2] = c[i].P();
      graf.Points(X, 3);
   }
   graf.Save();
}



void xComponent::StaticVacuum()
{
   xComponent& c = *this;

   for (int i = 1; i < Number(); i++)
   {
      c[i].l = c[i].p - c[i-1].p;
      c[i].Q += c[i].d * c[i].l * M_PI * c[i].G;
      c[i].U = (c[i].d^3) / c[i].l;
   }

   doubleU Rl, Rr;
   for (int i = 1; i < Number(); i++)
   if (c[i].Q() != 0)
   {
      c[i].Sl = 0;
      c[i].Sr = 0;
      for (int j = 0; j < Number(); j++)
         c[j].Sn = 0;

      for (int j = 1; j < Number(); j++)
      if (c[j].S() != 0)
      {
         Rl = 0; //1. / c[j].S;
         Rr = 0; //1. / c[j].S;
         //if (i == j) c[i].Sn = c[j].S;
         if (i > j)
         {  for (int k = j; k < i; k++)
               Rl += 1. / c[k].U;
            for (int k = 1; k < j; k++)
               Rr += 1. / c[k].U;
            for (int k = i; k < Number(); k++)
               Rr += 1. / c[k].U;

            for (int k = j; k < i; k++)
               c[k].Sn += 1. / Rl;
            for (int k = 1; k < j; k++)
               c[k].Sn -= 1. / Rr;
            for (int k = i; k < Number(); k++)
               c[k].Sn -= 1. / Rr;
         }
         if (i < j)
         {  for (int k = i; k < j; k++)
               Rr += 1. / c[k].U;
            for (int k = 1; k < i; k++)
               Rl += 1. / c[k].U;
            for (int k = j; k < Number(); k++)
               Rl += 1. / c[k].U;

            for (int k = i; k < j; k++)
               c[k].Sn -= 1. / Rr;
            for (int k = 1; k < i; k++)
               c[k].Sn += 1. / Rl;
            for (int k = j; k < Number(); k++)
               c[k].Sn += 1. / Rl;
         }
         c[i].Sl += (1. / Rl);//  + (1. / c[j].S);
         c[i].Sr += (1. / Rr);//  - (1. / c[j].S);
      }
      //c[0].Sn = c[Number()-1].Sn;
      for (int j = 1; j < Number(); j++)
         c[j].F += c[i].Q * c[j].Sn / (c[i].Sl + c[i].Sr);
   }

   for (int i = 1; i < Number(); i++)
      c[i].P = c[i-1].P + c[i].F;
}




