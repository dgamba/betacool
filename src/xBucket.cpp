//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xDynamic.h"
#include "xBucket.h"
//---------------------------------------------------------------------------

xBarrier::xBarrier()
{
   s1._(m_);
   s2._(m_);
   Ub._(V_);
   Us._(0, V_);
   V._(m_ / s_);
   A._(m_/(s_^2));
   UP._(V_*m_);
}
//shen modified:  delete absolute of Eta
void xBarrier::SetBarrier(double x1, double x2, double u)
{
   s1 = x1;
   s2 = x2;
   Ub = u;
   //Us = 0;
   V = iRing.Energy.Velocity * iRing.Eta;

   doubleU U1(0, V_);
   if ((Ub)(V_) >  1e-9) U1 =-1000; else
   if ((Ub)(V_) < -1e-9) U1 = 1000;

   A = iRing.Energy.Z * U_e * (Ub+Us) * iRing.Energy.Velocity * iRing.Eta /
      (iRing.Circ * iRing.Energy.Momentum);
   H = (2 * (s2-s1) * iRing.Energy.Z * U_Abs(Ub+Us) * U_e /
      (iRing.Circ * V * iRing.Energy.Momentum))^0.5;
   if ((Ub+Us)() < 0) H = -H;
}

void xBarrier::SetS(vectorU& X, doubleU t)
{
   X[4] += X[5]*V*t + A*t*t/2;
   X[5] += A * t / (iRing.Energy.Velocity * U_Abs(iRing.Eta));
}

doubleU xBarrier::GetT(vectorU& X, doubleU s)
{
   doubleU t(s_);
   double sign = 1;
   if (X[5] < 0) sign = -1;
   t  = U_Abs(( -X[5]*V + sign * (( ((X[5]*V)^2) + 2*A*(s-X[4]) )^0.5) ) / A);
   SetS(X, t);
   return t;
}

//---------------------------------------------------------------------------

xBucket::xBucket()
{
   SyncPeriod._(1, s_);
   Analytic = false;
   Stationary = false;
   Moving = false;
   bCoeff = false;
   Ucoeff = 1;
   UPcoeff = -1.;
   vertex = 0;
   cut = false;
   Periodical = false;

   TimeStart ._(s_);
   TimeFinish._(s_);
   LongSpaceCharge = false;
}

vectorU xBucket::GetPeriod(vectorU X0, doubleU dt)
{
   if (X0[5] == 0) X0[5] = 1e-199;
   vectorU X1, X = X0;
   doubleU s(m_), ds(m_), Time(s_);
   int i0 = -1;
   int i1 = 0, inc;
   int Nreflect = 0;
   bool stop = false;
   bool once = true;
   bool Reflected;
   SyncPeriod = 0;

   for(int i = 0; i < barrier.Number; i++)
   {  if(X[4] <= barrier[i].s2)
      {  if (X[5]() > 0)
            inc = 1;
         else
            inc = -1;
         i1 = i;
         break;
      }
   }
//                   ^ dP/P
//               8   |   1
//         7 /=======|=======\ 2
//          /        |        \      S-So
//      ---(---------|---------)----->
//          \        |        /
//         6 \=======|=======/ 3
//               5   |   4
   do
   {  X1 = X;        
      Reflected = false;
      if (i1==i0 && ((X[5] / X0[5])() >= 0) && !Moving)
      {  s = X0[4];
         stop = true;
      }else
         if (X[5]() > 0) s = barrier[i1].s2;
         else            s = barrier[i1].s1;

      if (barrier[i1].A() == 0)                                // 1,4,5,8
      {
         if (X[5] == 0) X[5] = 1e-199;
         Time = U_Abs((s - X[4]) / (barrier[i1].V * X[5]));
         X[4] = s;
      }else
      if ( (X[5] / barrier[i1].A)() >= 0)                      // 3,7
      {  Time = barrier[i1].GetT(X, s);
      }else                                                    // 2,6
      {  ds = ((X[5]*barrier[i1].V)^2)/(barrier[i1].A * 2);
         if (U_Abs(ds) / U_Abs(s - X[4]) > 0.999999999)        // not reflected
         {  Time = barrier[i1].GetT(X, s);
         }else                                                 // reflected
         {  Time = -X[5] * barrier[i1].V / barrier[i1].A;
            X[4] -= ds;
            if (X[5]() > 0) X[5] =-1e-199;
            else            X[5] = 1e-199;
            Reflected = true;
         }
      }

      SyncPeriod += Time;
      if (dt() > 0 && dt < Time)
      {  barrier[i1].SetS(X1, dt);
         return X1;
      }
      dt -= Time;

      if(once)
      {  once = false;
         i0 = i1;
      }


      if (Reflected)
      {  inc = -inc;
         Nreflect++;
      }
      else
         i1 += inc;
      if (i1 < 0)
      {  i1 = barrier.Number-1;
         X[4] = barrier[i1].s2;
         Nreflect++;
      }
      if (i1 > barrier.Number-1)
      {  i1 = 0;
         X[4] = barrier[i1].s1;
         Nreflect++;
      }
   }while(!stop && Nreflect < 4);
   return X;
}

int xBucket::SetBucket(BData& bdata)
{
   if (bdata[bdata.Row()-1][0] != 0.5)
   {   bdata.Size(bdata.Row()+1, 2);
       bdata[bdata.Row()-1][0] = 0.5;
       bdata[bdata.Row()-1][1] = -bdata[0][1];
   }
   int n = bdata.Row();
   barrier.SetNumber(n);

   for (int j = 0; j < n; j++)
   if (bdata[j][0] <= -0.5 || bdata[j][0] > 0.5)
   {
      Warning("Position of barrier [ ", j, "]= ", bdata[j][0], "should be in range -0.5 < x <= 0.5", PressEnter);
      return 1;
   }

   for (int k = 0; k < n-1; k++)
   if (bdata[k][0] >= bdata[k+1][0])
   {
      Warning("Position of barrier [ ", k, "] =", bdata[k][0], "should be less than next barrier position=", bdata[k+1][0], PressEnter);
      return 1;
   }

   Ucoeff = 1;
   UPcoeff = -1.;
   barrier[0].SetBarrier(-iRing.Circ()/2, bdata[0][0] * iRing.Circ(), bdata[0][1]);
   barrier[0].UP = (barrier[0].s1 - barrier[0].s2) * barrier[0].Ub;

   for (int i = 1; i < n; i++)
   {  barrier[i].SetBarrier(bdata[i-1][0] * iRing.Circ(), bdata[i][0] * iRing.Circ(), bdata[i][1]);
      barrier[i].UP = barrier[i-1].UP + (barrier[i].s1 - barrier[i].s2) * barrier[i].Ub;
      if(Ucoeff < fabs(barrier[i].Ub(k_3*V_)))
         Ucoeff = fabs(barrier[i].Ub(k_3*V_));
      if(UPcoeff < fabs(barrier[i].UP(k_3*V_*m_)))
         UPcoeff = fabs(barrier[i].UP(k_3*V_*m_));
   }
   if (!Harmonic1)
      if (fabs(barrier[n-1].UP(V_*m_)) > 1e-8)
         Warning("Integral of barrier potentials [V m] = ", barrier[n-1].UP(V_*m_));
   return 0;
}

int xBucket::SetMoving(doubleU t)
{  int i = 0;

   if (Periodical)
   {
      doubleU period(BarData[BarData.Row()-1][0] - BarData[0][0], s_);
      while (t >= period)
         t -= period;
   }
   for (i = 0; i < BarData.Row()-1; i++)
      if (BarData[i][0] <= t(s_) && t(s_) < BarData[i+1][0])
         break;

   BData bdata;
   bdata.Size( (BarData.Col(i)-1)/2, 2);

   if (i == BarData.Row()-1)
   {
      for (int j = 0; j < bdata.Row(); j++)
      {
         bdata[j][0] = BarData[i][1+j*2];
         bdata[j][1] = BarData[i][2+j*2];
      }
   }else
   {
      double dt = (t(s_)-BarData[i][0])/(BarData[i+1][0]-BarData[i][0]);
      for (int j = 0; j < bdata.Row(); j++)
      {
         bdata[j][0] = BarData[i][1+j*2] + dt * (BarData[i+1][1+j*2]-BarData[i][1+j*2]);
         bdata[j][1] = BarData[i][2+j*2] + dt * (BarData[i+1][2+j*2]-BarData[i][2+j*2]);
      }
   }
   SetBucket(bdata);
   return 0;
}

void xBucket::CalcParticle(xBeam& beam, xRing& ring)
{
   for (int i = 0; i < barrier.Number; i++)
   {  barrier[i].N = 0;
      barrier[i].M = 0;
      barrier[i].P = 0;
   }

   for (int j = 0; j < beam.Number(); j++)
   for (int i = 0; i < barrier.Number; i++)
   if  (beam(j,4) <= barrier[i].s2)
   {  barrier[i].N += 1;
      barrier[i].M += beam[j][5];
      break;
   }
   double Dmax = 0;
   for (int i = 0; i < barrier.Number; i++)
   {  barrier[i].D = barrier[i].N / beam.Number() *
                     (ring.Circ/(barrier[i].s2-barrier[i].s1))(U1_);
      if(Dmax < barrier[i].D)
         Dmax = barrier[i].D;
      if(barrier[i].N)
         barrier[i].M /= barrier[i].N;
   }
   for (int j = 0; j < beam.Number(); j++)
   for (int i = 0; i < barrier.Number; i++)
   if  (beam(j,4) <= barrier[i].s2)
   {  barrier[i].P += pow(beam[j][5] - barrier[i].M, 2);
      break;
   }
   for (int i = 0; i < barrier.Number; i++)
      if(barrier[i].N)
         barrier[i].P /= barrier[i].N;

   barrier[0].I = barrier[0].N;
   for (int i = 1; i < barrier.Number; i++)
      barrier[i].I = barrier[i-1].I + barrier[i].N;
   Dmax /= barrier[barrier.Number-1].I;
   for (int i = 0; i < barrier.Number; i++)
      barrier[i].I *= Dmax;

   if (iRing.Bucket.LongSpaceCharge)
   {
	   doubleU a(((iBeam.Emit[0] * iRing.BetaH) + (iBeam.Emit[1] * iRing.BetaV)) ^ 0.5);
	   double g = 1. + 2. * U_Ln(iRing.a_m / a);
	   for (int i = 0; i < barrier.Number; i++)
	   {
		   doubleU Usc(V_);
		   Usc = (iRing.Energy.Z * U_e * iBeam.Emit[3] * g) / (iRing.Energy.Gamma ^ 2) / (barrier[i].s2 - barrier[i].s1);
		   double p = 0;
		   if ((U_Abs(barrier[i].Ub))(V_) > 1e-9)
		   {
			   if (i == 0)
				   barrier[i].Us = Usc * (barrier[barrier.Number - 2].D - barrier[i + 1].D);
			   else if (i == barrier.Number - 1)
				   barrier[i].Us = Usc * (barrier[i - 1].D - barrier[1].D);
			   else
				   barrier[i].Us = Usc * (barrier[i - 1].D - barrier[i + 1].D);
			   p = barrier[i].Us(V_);
		   }
		   else
			   barrier[i].Us = 0;
	   }
   }
/*
   int index;
   double Pmax = 0;
   for (int i = 0; i < Profile.Number; i++)
      Profile[i] = 0;
   for (int j = 0; j < beam.Number(); j++)
   {  index = Round( (beam[j][4] + iRing.Circ(m_)/2) / iRing.Circ(m_) * Profile.Number );
      if (index >= Profile.Number) index = Profile.Number - 1;
      Profile[index] += double(Profile.Number) / beam.Number();
      if(Pmax < Profile[index])
         Pmax = Profile[index];
   }
   Integral[0] = Profile[0];
   for (int i = 1; i < Profile.Number; i++)
      Integral[i] = Integral[i-1] + Profile[i];
   Pmax /= Integral[Profile.Number-1];
   for (int i = 0; i < Profile.Number; i++)
      Integral[i] *= Pmax;
*/
}
/*
double xBucket::Density(xBeam& beam, xRing& ring, int j)
{
   if (beam.benum != BUCKET) return 1.;
   if (!IBSnorma)
   {
      int index = Round( (beam[j][4] + iRing.Circ(m_)/2) / iRing.Circ(m_) * Profile.Number );
      if (index >= Profile.Number) index = Profile.Number - 1;
      return Profile[index];
   }else
      return barrier[Pos(beam, ring, j)].D;
}
*/
int xBucket::Pos(xBeam& beam, int j)
{
   for (int i = 0; i < barrier.Number; i++)
   if  (beam(j,4) <= barrier[i].s2)
   {  return i;
   }
   return 0;
}

double xBucket::DroDs(xBeam& beam, int j)
{
	int i = Pos(beam, j);
	if ((U_Abs(barrier[i].Ub))(V_) > 1e-9)
	{	if (i == 0)     return barrier[barrier.Number - 2].D - barrier[i + 1].D; else
		if (i == barrier.Number - 1) return barrier[i - 1].D - barrier[1].D; else
			return barrier[i - 1].D - barrier[i + 1].D;
	}
	else
		return 0;
}


double xBucket::Luminosity(int i1, int i2, double beta, double ds)
{
   LumiTray.SetNumber(i2-i1, i2-i1);
   double b, s;
   double L0 = 0, L1 = 0;
   for (int i = 0; i < i2-i1; i++)
   {
      for (int j = 0; j < i2-i1; j++)
      {
         s = ds * (j - i);
         b = beta + s*s/beta;
         if (cut)
         {
            if (vertex/100. < fabs(s))
               LumiTray[i][j] = 0;
            else
               LumiTray[i][j] = barrier[i+i1].D * barrier[j+i1].D / b;
         }else
            LumiTray[i][j] = barrier[i+i1].D * barrier[j+i1].D / b;

         L0 += barrier[i+i1].D * barrier[j+i1].D / beta;
         L1 += LumiTray[i][j];
      }
   }
   return L1 / L0;
}

void xBucket::LuminosityTest()
{

   HourGlass = Luminosity(1, barrier.Number-1,
               sqrt(iColl.Lattice.betax(m_)*iColl.Lattice.betay(m_)),
               (barrier[2].s2-barrier[2].s1)(m_)/2.);
   iDraw.Luminosity3D.SetXaxis(barrier.Number-3, barrier[1].s1()/2., barrier[barrier.Number-2].s2()/2.);
   iDraw.Luminosity3D.SetYaxis(barrier.Number-3, barrier[1].s1()/2., barrier[barrier.Number-2].s2()/2.);
   for (int i = 0; i < barrier.Number-2; i++)
   for (int j = 0; j < barrier.Number-2; j++)
      iDraw.Luminosity3D.SetPoint(i+1,j+1,LumiTray[i][j]);
   iDraw.Luminosity3D.Save();   

}

void xBucket::Generate()
{
	/*
   if (GenerateStep > TimeFinish1)
   {  Warning("Finish time smaller than generation step");
      return;
   }   
   double MaxTime = TimeFinish1();
   if (TimeFinish1 < TimeFinish2) MaxTime = TimeFinish2();
   int steps = MaxTime / GenerateStep();

   BData HarmonicRF(steps+1, 5);
   doubleU k1(s_^-1), k3(s_^-1), t(s_), k2(V_ / s_);
   k1 = (1 - ((AmplStart1 / AmplFinish1) ^ 0.5)) / (TimeFinish1 - TimeStart1);
   k3 = (1 - ((AmplStart2 / AmplFinish2) ^ 0.5)) / (TimeFinish2 - TimeStart2);
   k2 = (AmplFinish2 - AmplStart2) / (TimeFinish2 - TimeStart2);

   for (int i = 0; i < steps+1; i++)
   {  t = GenerateStep * i;
      HarmonicRF[i][0] = t(s_);
      HarmonicRF[i][1] = 0;
      HarmonicRF[i][3] = 0.2;

      if (t <= TimeStart1)
          HarmonicRF[i][2] = AmplStart1(V_);
      else
      if (t >= TimeFinish1)
          HarmonicRF[i][2] = AmplFinish1(V_);
      else
          HarmonicRF[i][2] = (AmplStart1 / ( (1 - ((t-TimeStart1)*k1) )^2))(V_);
      if (Harmonic2)
      {
		  if (t <= TimeStart2)
			  HarmonicRF[i][4] = 0; // AmplStart2(V_);
         else
         if (t >= TimeFinish2)
             HarmonicRF[i][4] = AmplFinish2(V_);
         else
			 HarmonicRF[i][4] = (AmplStart2 / ((1 - ((t - TimeStart2)*k3)) ^ 2))(V_);
		 //HarmonicRF[i][4] = (AmplStart2 + (k2*(t-TimeStart2)))(V_);
       }else
         HarmonicRF[i][4] = 0;
   }
   HarmonicRF.Save("harmonicrf.mov");
   */
}
