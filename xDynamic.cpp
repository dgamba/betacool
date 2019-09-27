//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xdynamic.h"
#include "bptune.h"
#include "bpdata.h"
#include <iostream>
#include <time.h>
#include <omp.h>
//---------------------------------------------------------------------------

xDynamics::xDynamics()
{
   dt._(s_);
   Maxdt._(s_);
   ds._(cm_);
   IntegrationStep._(s_);
   Algoritm = -1;
   Synchrotron = 0;
   LongSC = false;
   LongOK = false;
   betatron = true;
   rotate = true;
   angle._(1, s_^-1);
   Time_Stop._(s_);
}

int xDynamics::OnGet()
{
   //Algoritm       = Data.OnGet(21, 1, 0);
   SavingInterval = Data.OnGetI(21, 2, 5);
   SkipPoints     = Data.OnGetI(21, 4, 0);
   SkipCount      = 0;
   CurveSize      = Data.OnGetI(21, 5, 2000);
   AutoSkip       = Data.OnGetB(21, 6, true);
   xGraf::ResetAll(CurveSize, AutoSkip);

   dt     = Data.OnGet(22, 1, 1.);
   Maxdt  = Data.OnGet(22, 2, 1.);
   ifmore = Data.OnGet(22, 3, 2.);
   ratio  = Data.OnGet(22, 4, 10.);
   ratio = ratio / 100. + 1.;

   IntegrationStep = Data.OnGet(23, 1, 1.);
   TurnNumber = Data.OnGet(23, 3, 1.);
   TimeTurn = Data.OnGetB(23, 4, true);
   Synchrotron = Data.OnGetI(23, 5, 0);
   LongSC = Data.OnGetB(23, 12, false);
   Time_Stop = Data.OnGet(23, 7, 0);
   betatron = Data.OnGetI(23, 8, 0);
   rotate = Data.OnGetB(23, 10, true);
   angle = Data.OnGet(23, 11, 1.);

   ds            = Data.OnGet(25, 1, 1.);
   CheckCrossing = Data.OnGetB(25, 2, false);
   Crystall      = Data.OnGetB(25, 3, false);
   FORCES        = Data.OnGetB(25, 4, false);

   return 0;
}

int xDynamics::OnSet()
{
   if (TimeTurn == 0)
   {  if (IntegrationStep() < 0)
          IntegrationStep = iRing.Trev;
      TurnNumber = IntegrationStep / iRing.Trev;
   }else
   {  if (TurnNumber() < 0)
          TurnNumber = 1;
      IntegrationStep = TurnNumber * iRing.Trev;
   }

   Data.OnSet(23,1,IntegrationStep);
   Data.OnSet(23,3,TurnNumber);

   return 0;
}

int xDynamics::OnRun()
{
   //Algoritm       = Data.OnGet(21, 1, 0);
   SavingInterval = Data.OnGetI(21, 2, 5);
   SkipPoints     = Data.OnGetI(21, 4, 0);
   SkipCount      = 0;
   CurveSize      = Data.OnGetI(21, 5, 2000);
   AutoSkip       = Data.OnGetB(21, 6, true);
   xGraf::ResetAll(CurveSize, AutoSkip);

   //dt     = Data.OnGet(22, 1, 1.);
   Maxdt  = Data.OnGet(22, 2, 1.);
   ifmore = Data.OnGet(22, 3, 2.);
   ratio  = Data.OnGet(22, 4, 10.);
   ratio = ratio / 100. + 1.;

   IntegrationStep = Data.OnGet(23, 1, 1.);
   TurnNumber = Data.OnGet(23, 3, 1.);
   TimeTurn = Data.OnGetB(23, 4, true);
   Synchrotron = Data.OnGetI(23, 5, 0);
   LongSC = Data.OnGetB(23, 12, false);
   Time_Stop = Data.OnGet(23, 7, 0);
   betatron = Data.OnGetI(23, 8, 0);
   rotate = Data.OnGetB(23, 10, true);
   angle = Data.OnGet(23, 11, 1.);
   /*
   ds            = Data.OnGet(25, 1, 1);
   CheckCrossing = Data.OnGet(25, 2, false);
   Crystall      = Data.OnGet(25, 3, false);
   FORCES        = Data.OnGet(25, 4, false);
   */
   OnSet();

   return 0;
}
//---------------------------------------------------------------------------
void xDynamics::Distribution(vectorU& emit, xLattice& L, int num, bool crystall)
{
   iBeam.Distribution(emit, L, num);
   if (iBeam.initial == 2)
   {
      BData InjFile;
      InjFile.Load(Data.OnGetC(14,8,"NoFile"));
      iBeam.Number(num + InjFile.Row());
      int count = num;
      for (int i = 0; i < InjFile.Row(); i++)
      {  if (InjFile[i][0] !=EmptyData)
         {  for(int j = 0; j < 6; j++)
               iBeam[num+i][j] = InjFile[i][j];
            count++;
         }
      }
      iBeam.Number(count);
   }else
   {
      if (crystall)
         for (int i = 0; i < iBeam.Number(); i++)
            iBeam[i][4]=iBeam.CellSize() * i / iBeam.Number();
   }
   for (int i = num; i < iBeam.Number(); i++)
   {  iBeam[i][0] += iBeam.horshift;
      iBeam[i][2] += iBeam.vershift;
      iBeam[i][4] += iBeam.longshift;
      iBeam[i][5] += iBeam.dpshift;
   }                                       
}

bool xDynamics::Drawing(xLattice& L1)
{
   static vectorU emit1;
   vectorU rates(s_^-1);

   static doubleU time1(s_);
   bool rtrn = false;   
   if (!SkipCount || (SkipCount > SkipPoints))
   {  	   
	   SkipCount = 0;
      iBeam.Emit = iBeam.Emit_Def(L1);
      if(iBeam.benum == BUNCHED)
         iBeam.CalcBunch();
      else
      if(iBeam.benum == BUCKET)
         iBeam.CalcBucket(iBeam.Emit[2]);	  
      iDraw.DrawRealSpace (iBeam, iRing, L1, SavingInterval);	  
      rtrn = iDraw.DrawBeamEvolution (iBeam, iTime, iColl, SavingInterval);   	  
      if (iTime.t())
      {
         for (int i = 0; i < 4; i++)
            rates[i] = (iBeam.Emit[i] - emit1[i])/(emit1[i] * (iTime.t-time1));
         rates[2] /= 2.;   
         iDraw.BeamRates(rates, iTime, SavingInterval);		 
      }	  
      emit1 = iBeam.Emit;
      time1 = iTime.t;
      if (rtrn && Loader.Check("bolide.run"))
      {  Warning("BETACOOL parameters have been changed !");
         xData::Run(xData::File);
         if(iBeam.benum == BUNCHED)
            iBeam.CalcBunch();
         else
         if(iBeam.benum == BUCKET)
            iBeam.CalcBucket(iBeam.Emit[2]);
      }
   }      
   SkipCount++;
   if(((Time_Stop() > 0) && (iTime.t >= Time_Stop)) || xTimer::finish)
   {  ShowTime("STOP -");
      exit(11);
   }
   return rtrn;
}

//================================ RMS Dynamics ================================

void xDynamics::Dynamics()
{  Algoritm = 0;
   iDraw.DrawGammaPi(iBeam);
   iDraw.DrawGamma3(iBeam);
   Warning("Dynamics: Initial step [sec] =", dt(), ", Step multiplier =",ifmore,
           ", Max growth =", ratio);
   if ((iRing.LatticeFile == 2) && iIBS.Model)
   {  Warning("<Ring | Optics structure | No optics file> - only for <IBS | Piwiski>", PressEnter);
      return;
   }

   vectorU rate0;
   vectorU emit0;
   bool negative;
   iTime.TimeStep(dt);
   while (true)
   {
      rate0 = xEffect::Summary(iTime, iBeam, iRing);
      emit0 = iBeam.Emit;
      do
      {  negative = false;
         iBeam.Emit *=  rate0 * iTime.dt + 1;            // Euler step
         for (int i = 0; i < 3; i++)
         {  if (iBeam.Emit[i]() <= 0)                    // Check on negative values
               negative = true;
            else
            if (((iBeam.Emit[i]/emit0[i])() > ratio) ||  // Check on too large growth
                ((emit0[i]/iBeam.Emit[i])() > ratio) )
               negative = true;
         }
         if (negative)
         {  iTime *= 1./ifmore;                         // If any then decrease step
            iBeam.Emit = emit0;
         }
      }while(negative);

      iTime.TimeStep(iTime.dt);
      if (!SkipCount || (SkipCount > SkipPoints))
      {  SkipCount = 0;
         //iBeam.Luminosity();
         iDraw.BeamRates(rate0, iTime, SavingInterval);
         if (iDraw.DrawBeamEvolution(iBeam, iTime, iColl, SavingInterval))
         {  if (Loader.Check("bolide.run"))
            {  Warning("BETACOOL parameters have been changed !");
               xData::Run(xData::File);
            }
                   }
      }
      SkipCount++;
      ++iTime;
      if(iTime.dt > Maxdt/2)
         iTime.dt = Maxdt;
      else
         iTime *= ifmore;            // Increase step
      if(xTimer::finish)
      {  ShowTime("STOP TIME -");
         exit(11);
      }
   }
}
//================================ Model Beam ==================================

void xDynamics::ModelBeam()
{  Algoritm = 1;
   iDraw.DrawGammaPi(iBeam);
   iDraw.DrawGamma3(iBeam);
   Warning("Model Beam: Integration Step [sec] =", IntegrationStep());

   if ((iRing.LatticeFile == 2) && iIBS.Model)
   {  Warning("<Ring | Optics structure | No optics file> - only for <IBS | Piwinski>", PressEnter);
      return;
   }
   if ((iIBS.IBSkickModel == 5) || (iIBS.IBSkickModel == 7))
   {
	  if (iRing.LatticeFile != 1)
      {  Warning("Choose <Ring | Optics structure | Input MAD file>", PressEnter);
         return;
      }
      iIBS.Lattice = iRing.GetLattice(0);
   }
   iDraw.SpaceChargeDraw(iBeam,iEbeam);
   LongHyst.SetNumber(6,iDraw.Divisions+1);
   matrixU Matrix(4,4);
   matrixU Vector(4,1);
   xLattice Lattice2;

   Distribution(iBeam.InitialEmit,iRing.LATTICE, 0, false);
   if (iBeam.initial == 2)
   {  iBeam.Emit = iBeam.Emit_Def(iRing.LATTICE);
      iBeam.InitialEmit = iBeam.Emit;
   }
   iTime.TimeStep(IntegrationStep);
   xTime jtime = iTime;          // injection timer
   xTime etime = iTime;          // 3D evolution timer
   xTime dtime = iTime;          // time divider
   etime.t = iDraw.EvolutionStep;
   cycles = 1;      
   while(true)
   {	 	   
      if(iBeam.inject && (cycles < iBeam.cycles) && (jtime.t >= iBeam.interval))
      {  jtime.t = 0;
         cycles++;

         if (iLosses.Multi && iLosses.Injection)
         {
			iBeam.Courant_Snyder(iRing.LATTICE);
			
			for (int j = 0; j < iBeam.Number(); j++)
			   {
				if (((sqrt(iBeam.Invs[j][2])*fabs(iRing.LATTICE.Dx())+sqrt(iBeam.Invs[j][0]*iRing.LATTICE.betax()))> sqrt(iRing.AcceptHinj()*iRing.LATTICE.betax())) ||
					((sqrt(iBeam.Invs[j][2])*fabs(iRing.LATTICE.Dy())+sqrt(iBeam.Invs[j][1]*iRing.LATTICE.betay()))> sqrt(iRing.AcceptVinj()*iRing.LATTICE.betay())) ||
					((iBeam.benum == BUCKET) && (iBeam.longshift-iRing.S1(m_) < iBeam(j,4)(m_)) && (iBeam(j,4)(m_) < iBeam.longshift+iRing.S1(m_))) )
					iBeam.Loss(iRing.LATTICE, j, 1.);
			   }
			   
		   }
         int num = iBeam.Number();
		   //iBeam.Number(num + iBeam.IniNumber);
		   iBeam.Number(num * ((iBeam.InitialEmit[3]/iBeam.Emit[3])(U1_) + 1.));
		   Distribution(iBeam.InitialEmit, iRing.LATTICE, num, false);
         iBeam.Emit[3] += iBeam.InitialEmit[3];
      }
      if (iBeam.MaxNumber)
	   while (iBeam.MaxNumber < iBeam.Number())
      {  int k = 0;
         for (int i = 0; i < iBeam.Number(); i++)
         if (i % 5)
         {  for (int j = 0; j < 6; j++)
               iBeam.Bunch[k][j] = iBeam.Bunch[i][j];
            k++;
         }
         iBeam.Number(k);
      }	  
      if(iDraw.Evolution.Enabled &&
         (etime.t()-iDraw.EvolutionStep >= iDraw.EvolutionStep * 1e-6))
      {  etime.t = 0;
         iDraw.evolution = true;
      }	  
      if(iBeam.benum == BUNCHED)
         iBeam.CalcBunch();
      if (iBeam.benum == BUCKET)
      {  iRing.Bucket.CalcParticle(iBeam, iRing);
         if (iRing.Bucket.Analytic && (iRing.Bucket.Stationary || iRing.Bucket.Moving))
            for (int j4 = 0; j4 < iBeam.Number(); j4++)
            {  while (iBeam(j4,4) > iRing.Circ/(2*iRing.h()))
                  iBeam[j4][4] -= iRing.Circ(m_)/iRing.h();
			   while (iBeam(j4,4) <-iRing.Circ/(2*iRing.h()))
                  iBeam[j4][4] += iRing.Circ(m_)/iRing.h();
            }
      }   	  
      iBeam.LongProfile();	  	  
      Drawing(iRing.LATTICE);	  
      Lattice2 = iRing.LATTICE;
      iTime.TimeStep(IntegrationStep);
      jtime.TimeStep(IntegrationStep);
      etime.TimeStep(IntegrationStep);	  
      double mu;
      bool once = true;
      ExpKick = 0;
      for (int j2 = 0; j2 < iBeam.Number(); j2++)
         iBeam.MinusDisp(Lattice2, j2);

      if(rotate)
      {  double a = M_PI*(iTime.dt*angle)(U1_);
         double ca = cos(a);
         double sa = sin(a);
         vectorU v;
         for (int i = 0; i < iBeam.Number(); i++)
         {  v = iBeam(i);
            iBeam(i,0, v[0]*ca - v[2]*sa);
            iBeam(i,1, v[1]*ca - v[3]*sa);
            iBeam(i,2, v[2]*ca + v[0]*sa);
            iBeam(i,3, v[3]*ca + v[1]*sa);
         }
      }	  
/*
      double dp;
      if (iRing.Tsweep()== 0)
         dp = iRing.dPini();
      else if (iRing.Tsweep < iTime.t)
         dp = iRing.dPfin();
      else
		 dp = ((iRing.dPfin - iRing.dPini) * iTime.t / iRing.Tsweep + iRing.dPini)(U1_);
*/

// Longitudinal rotation

      if (iBeam.benum == COASTING)
      {  for (int j4 = 0; j4 < iBeam.Number(); j4++)
         {  iBeam.Add(j4,4,iBeam(j4)[5]*iRing.Energy.Velocity*U_Abs(iRing.Eta)*dtime.dt);
            iBeam[j4][4] -= iBeam.CellSize() * floor(iBeam[j4][4] / iBeam.CellSize());
            iBeam[j4][4] -= iBeam.CellSize() / 2;
         }
      }else

      if (iBeam.benum == BUNCHED)
      {  
         if (Synchrotron == 2)
         {  xBeam& b = iBeam;
            doubleU s_s((xLattice::B_s) * (b.InitialEmit[2]^0.5));
            doubleU d_t(dtime.dt);
            if (LongSC) 
            {  LongOK = true;
               doubleU a( ((b.Emit[0]*iRing.BetaH) + (b.Emit[1]*iRing.BetaV))^0.5);
               double g = 1. + 2. * U_Ln(iRing.a_m/a); 
               for (int i = 0; i < b.Number(); i++)
                  b.Invs[i][5] = sqrt(b[i][5]*b[i][5] + b[i][4]*b[i][4]/(xLattice::B_s()*xLattice::B_s()));
               b.Profile(LongHyst[5], 5, (b.InitialEmit[2]^0.5)(), iDraw.Divisions);
            
               double Xo, h = 2 * iDraw.Sigmas / iDraw.Divisions;
               for (int j = 0; j <= iDraw.Divisions; j++)
               {  LongHyst[0][j] = -iDraw.Sigmas + h*j;
                  LongHyst[1][j] =  iDraw.Divisions*LongHyst[5][j]/iDraw.Sigmas/b.Number();
                  if(j)LongHyst[2][j] = g * (Xo - LongHyst[1][j]);
                  else LongHyst[2][j] = 0;
                  Xo = LongHyst[1][j];
                  LongHyst[3][j] = (iRing.Energy.Z * U_e * b.Emit[3] * LongHyst[2][j] / s_s)(V_);
               }
            }
#pragma omp parallel for
            for (int j4 = 0; j4 < b.Number(); j4++)
            {  doubleU Vsc(0, V_);
               if (LongSC) {  
                  double s2 = 2. * iDraw.Sigmas * s_s(m_);
                  int ind = Round((b[j4][4] + s2/2.) / s2 * iDraw.Divisions);
                  Vsc = iRing.Energy.Z * U_e * b.Emit[3] * LongHyst[2][ind] / s_s;
               }
               b[j4][5] -= ( iRing.Energy.Z * U_e * ( iRing.V * U_Sin ( b(j4,4)*iRing.h*2*M_PI/iRing.Circ) - Vsc) /
                           (iRing.Circ * iRing.Energy.Momentum) * d_t)(U1_);
               b[j4][4] += (U_Abs(iRing.Eta) * iRing.Energy.Beta * U_c * d_t)(m_) * b[j4][5];
            }
         }else 
         if (iRing.NonLinear)
         {
            doubleU step(dtime.dt);
            doubleU period(iRing.Trev / iRing.Q_s);
            if (step > period) step = period;
            int divi = (int)(step * 36 / period)(U1_);
            if (divi == 0) divi = 1;
            doubleU d_t(step / divi);

            static bool oncediv = true;
            if(oncediv && dtime.dt > period)
               divi *= 36;
            oncediv = false;

            for (int d = 0; d < divi; d ++)
            for (int j4 = 0; j4 < iBeam.Number(); j4++)
            {
               iBeam[j4][5] -= (iRing.Energy.Z * U_e *
               (iRing.V    *U_Sin(iBeam(j4,4)*iRing.h    *2*M_PI/iRing.Circ) +
                iRing.V_low*U_Sin(iBeam(j4,4)*iRing.h_low*2*M_PI/iRing.Circ))/
               (iRing.Circ * iRing.Energy.Momentum) * d_t)(U1_);
               iBeam[j4][4] += (U_Abs(iRing.Eta) * iRing.Energy.Beta * U_c *
                  d_t)(m_) * iBeam[j4][5];
            }
         }else
         {
            if (Synchrotron == 1)
               mu = 2 * M_PI * fabs(xDistributor::Flatten());
            else
               mu = (2*M_PI*iRing.Q_s*dtime.dt/iRing.Trev)();
            for (int j4 = 0; j4 < iBeam.Number(); j4++)
            {  //if (dp) iBeam[j4][5] -= dp;
               iBeam[j4][4] = cos(mu)*iBeam[j4][4]+ xLattice::B_s() * sin(mu)*iBeam[j4][5];
               iBeam[j4][5] = cos(mu)*iBeam[j4][5]- (((iBeam[j4][4]-xLattice::B_s() * sin(mu)*iBeam[j4][5])/cos(mu))/xLattice::B_s())*sin(mu);
              //if (dp) iBeam[j4][5] += dp;
            }
         }
      }else
      if (iBeam.benum == BUCKET) 
      {
		  if (!iRing.Bucket.Analytic)
		  {  if (iRing.Bucket.Moving)
				iRing.Bucket.SetMoving(dtime.t);
			else
				iRing.Bucket.SetBucket(iRing.Bucket.BarData);
//#pragma omp parallel for - not ok !
			  for (int j = 0; j < iBeam.Number(); j++)
			  {  
				  doubleU dt(0,s_);
				  vectorU X;
				  X = iBeam(j);
				  if (iRing.Bucket.Moving)
					  dt = dtime.dt; 
				  else
				  {  iRing.Bucket.GetPeriod(X, dt);
				     if (Synchrotron)
					     dt = iRing.Bucket.SyncPeriod * fabs(xDistributor::Flatten());
				     else
					     dt = dtime.dt - iRing.Bucket.SyncPeriod * floor((dtime.dt/iRing.Bucket.SyncPeriod)());
				  }	
				  X = iRing.Bucket.GetPeriod(X, dt);
				  iBeam[j][4] = X[4](m_);
				  iBeam[j][5] = X[5](U1_);
			  }	
       }else
       if (iRing.Bucket.Moving || iRing.Bucket.Stationary)  // 2 sin harmonics (Ring | Barrier Bucket | Harmonic RF)
       {  if (iRing.Bucket.Moving)
         {  iRing.Bucket.SetMoving(dtime.t);
            iRing.V = iRing.Bucket.barrier[0].Ub;
            if (iRing.Bucket.Harmonic2)
                iRing.V_low = iRing.Bucket.barrier[1].Ub;
            else
                iRing.V_low = 0;
         }  

		complex<double> Ah(0, 0), Z(260, 2.64e4), Uind(0, 0);
		if (iRing.Bucket.Induction && (iTime.t >= iRing.Bucket.TimeStart) && (iTime.t <= iRing.Bucket.TimeFinish))
	    {   
		   for (int j = 0; j < iBeam.Number(); j++)
		   {
			   complex<double>ih(0, (iRing.h_low * iBeam(j, 4) * 2 * M_PI / iRing.Circ)(U1_));
			   Ah += exp(-ih);
		   }
		   Ah /= iBeam.Number();
		   Uind = -Ah;
		   Uind = -Ah * iRing.Bucket.Impedance * (iRing.Energy.Z * U_e * iBeam.Emit[3] / iRing.Trev)(A_);
	    }
#pragma omp parallel for
         for (int j4 = 0; j4 < iBeam.Number(); j4++)
		 {
			 doubleU V(0, V_), V1(0, V_) ;
			 double p1;
			 doubleU  xx1(0, V_), xx2(0, V_);
			 if (iRing.Bucket.Induction && (iTime.t >= iRing.Bucket.TimeStart) && (iTime.t <= iRing.Bucket.TimeFinish))
			 {
				 complex<double>ih(0, (iRing.h_low * iBeam(j4, 4) * 2 * M_PI / iRing.Circ)(U1_));
				 ih = Uind * (exp(ih) - 1.);
				 V = 2. * ih.real();
			 }
			 xx1 = iRing.V;
			 xx2 = U_Sin(1.57);
			 V1 = iRing.V    *U_Sin(iBeam(j4, 4)*iRing.h * 2 * M_PI / iRing.Circ) +
				 iRing.V_low*U_Sin(iBeam(j4, 4)*iRing.h_low * 2 * M_PI / iRing.Circ) + V;

			 p1 = iBeam[j4][5];
			 iBeam[j4][5] -= (iRing.Energy.Z * U_e *
            (iRing.V    *U_Sin(iBeam(j4,4)*iRing.h    *2*M_PI/iRing.Circ) +
             iRing.V_low*U_Sin(iBeam(j4,4)*iRing.h_low*2*M_PI/iRing.Circ) + V)/
            (iRing.Circ * iRing.Energy.Momentum) * dtime.dt)(U1_);

  //         iBeam[j4][4] += (U_Abs(iRing.Eta) * iRing.Energy.Beta * U_c * dtime.dt)(m_) *( iBeam[j4][5]+p1)*0.5;
			iBeam[j4][4] += (U_Abs(iRing.Eta) * iRing.Energy.Beta * U_c * dtime.dt)(m_)*iBeam[j4][5];

            while (iBeam(j4,4) > iRing.Circ/(2*iRing.h()))
               iBeam[j4][4] -= iRing.Circ(m_)/iRing.h();
            while (iBeam(j4,4) <-iRing.Circ/(2*iRing.h()))
               iBeam[j4][4] += iRing.Circ(m_)/iRing.h();
         }
       }else                                                // 2 analytic barriers
       { doubleU RFlength(m_), T_int(s_), T_curr(s_);
         vectorU X(m_, U1_);
         double rand1;
         RFlength = 0.5* iRing.T2 * iRing.Circ/iRing.Trev;
         for (int j5 = 0; j5 < iBeam.Number(); j5++)
         {  iBeam.Amp2 = iBeam[j5][5]*iBeam[j5][5];
            if(fabs(iBeam[j5][4])>(RFlength(m_) + iRing.T1(s_) * iRing.Circ(m_)/iRing.Trev(s_)))
               iBeam.Amp2 += iRing.BarHeight*iRing.BarHeight;
            else
            {  if(iBeam[j5][4] > RFlength(m_))
               iBeam.Amp2 += (2.0* U_e *iRing.Energy.Z * iRing.V_B *(iBeam[j5][4]-RFlength(m_)))/
                       (U_Abs(iRing.Eta)*iRing.Energy.Velocity*iRing.Circ(m_)*iRing.Energy.Momentum);
               if(iBeam[j5][4] < (-1.0*RFlength(m_)))
               iBeam.Amp2 += (-2.0* U_e *iRing.Energy.Z * iRing.V_B *(iBeam[j5][4]+RFlength(m_)))/
                       (U_Abs(iRing.Eta)*iRing.Energy.Velocity*iRing.Circ(m_)*iRing.Energy.Momentum);
            }
            if((iBeam.Amp2 < (iRing.BarHeight*iRing.BarHeight))&&(iBeam.Amp2()>0.0))
            {  iBeam.t_1 =  iRing.T2/(2.0*U_Abs(iRing.Eta)*(iBeam.Amp2^0.5));
               iBeam.t_2 =  iRing.Circ*iRing.Energy.Momentum*(iBeam.Amp2^0.5)/( U_e *iRing.Energy.Z * iRing.V_B);
               iBeam.T_s = (2.0*iRing.T2/((U_Abs(iRing.Eta)*(iBeam.Amp2^0.5))))+
                     (4.0*iRing.Circ*iRing.Energy.Momentum*(iBeam.Amp2^0.5)/( U_e *iRing.Energy.Z * iRing.V_B));
               if(Synchrotron)
               {  rand1 = double(rand()) / RAND_MAX;
                  T_int = iBeam.T_s*rand1;
               }else
               {   if(fabs(iBeam[j5][4]) < RFlength(m_))
                   {  if(iBeam[j5][4] > 0)
                      {  T_curr = iBeam.t_1 * (iBeam[j5][4]/RFlength(m_));
                         if(iBeam[j5][5] < 0) T_curr = (2.0*iBeam.t_1) + (2.0*iBeam.t_2) - T_curr;
                      }else
                      {  T_curr = iBeam.T_s + ((iBeam[j5][4]/RFlength(m_))*iBeam.t_1);
                         if(iBeam[j5][5] < 0) T_curr = (2.0*iBeam.t_1) + (2.0*iBeam.t_2) + (iBeam.T_s - T_curr);
                      }
                   }else
                   {  if(iBeam[j5][4] > 0)
                      {  if(iBeam[j5][5] > 0)
                            T_curr = iBeam.t_1 + (((iBeam.Amp2^0.5) - iBeam[j5][5])*iRing.Circ*iRing.Energy.Momentum/
                            (U_e *iRing.Energy.Z * iRing.V_B));
                         else
                            T_curr = iBeam.t_1 + (2.0*iBeam.t_2) - (((iBeam.Amp2^0.5) + iBeam[j5][5])*iRing.Circ*iRing.Energy.Momentum/
                            (U_e *iRing.Energy.Z * iRing.V_B));
                      }else
                      {  if(iBeam[j5][5] > 0)
                            T_curr = (3.0*iBeam.t_1) + (4.0*iBeam.t_2) - (((iBeam.Amp2^0.5) - iBeam[j5][5])*iRing.Circ*iRing.Energy.Momentum/
                            (U_e *iRing.Energy.Z * iRing.V_B));
                         else
                            T_curr = (3.0*iBeam.t_1) + (2.0*iBeam.t_2) + (((iBeam.Amp2^0.5) + iBeam[j5][5])*iRing.Circ*iRing.Energy.Momentum/
                            (U_e *iRing.Energy.Z * iRing.V_B));
                      }
                   }
                  T_int = T_curr + dtime.dt;
                  while(T_int > iBeam.T_s) T_int -= iBeam.T_s;
               }
               iBeam.GetBucket(X, T_int);
               iBeam[j5][4] = X[4](m_);
               iBeam[j5][5] = X[5](U1_);
            }
            //if(Amp2() > (iRing.BarHeight()*iRing.BarHeight()))
            //   iBeam[j5][4] = 0.5*iRing.Circ(m_)*xDistributor::Flatten();
         }
      }}

// Kicks from Effects
    
      for (int i = 0; i < xEffect::ACount; i++)
      if (xEffect::AItems[i]->Multi != 0)
      {
         dtime = iTime;
         int divider = abs(xEffect::AItems[i]->Multi);

         if (xEffect::AItems[i]->Multi > 0)
         {  dtime *= divider;
            divider = 1;
         }else
            dtime *= 1. / divider;

         if ((iTime.Steps % xEffect::AItems[i]->Multi == 0) ||
             (xEffect::AItems[i]->Multi < 0))
         for (int d = 0; d < divider; d++)
         {
            if (betatron == 0)
            {  Lattice2.mux = iRing.TunesH * TurnNumber;
               Lattice2.muy = iRing.TunesV * TurnNumber;
               Matrix = iRing.CalcTransformMatrix(Lattice2, xEffect::AItems[i]->Lattice);
            }else
            if (betatron == 1)
            {  Lattice2.mux = M_PI * (1. + xDistributor::Flatten());
               Lattice2.muy = M_PI * (1. + xDistributor::Flatten());
               Matrix = iRing.CalcTransformMatrix(Lattice2, xEffect::AItems[i]->Lattice);
            }
            for (int j0 = 0; j0 < iBeam.Number(); j0++)
            {  if (betatron == 2 && once)
               {
                  Lattice2.mux = M_PI * (1. + xDistributor::Flatten());
                  Lattice2.muy = M_PI * (1. + xDistributor::Flatten());
                  Matrix = iRing.CalcTransformMatrix(Lattice2, xEffect::AItems[i]->Lattice);
               }
               iBeam.Matrix_2x2(j0, Matrix);
            }  once = false;

            if (iBeam.benum == BUNCHED && !iRing.NonLinear && Synchrotron == 1)
            {  mu = M_PI * (1. + xDistributor::Flatten());
               for (int j4 = 0; j4 < iBeam.Number(); j4++)
               {  //if (dp) iBeam[j4][5] -= dp;
                  iBeam[j4][4] = cos(mu)*iBeam[j4][4]+ xLattice::B_s() * sin(mu)*iBeam[j4][5];
                  iBeam[j4][5] = cos(mu)*iBeam[j4][5]- (((iBeam[j4][4]-xLattice::B_s() * sin(mu)*iBeam[j4][5])/cos(mu))/xLattice::B_s())*sin(mu);
                  //if (dp) iBeam[j4][5] += dp;
               }
            }

            for (int j1 = 0; j1 < iBeam.Number(); j1++)
               iBeam.PlusDisp (xEffect::AItems[i]->Lattice, j1);
            xEffect::AItems[i]->Kick(dtime, iBeam, iRing);
            for (int j2 = 0; j2 < iBeam.Number(); j2++)
               iBeam.MinusDisp(xEffect::AItems[i]->Lattice, j2);
            Lattice2 = xEffect::AItems[i]->Lattice;
         }
      }

      for (int j1 = 0; j1 < iBeam.Number(); j1++)
      {
         iBeam.PlusDisp (Lattice2, j1);
         if (iRing.Inductor)
            iBeam.Add(j1, 5, U_e * iRing.Vind * iTime.dt / iRing.Trev /
            (iRing.Energy.Momentum * iRing.Energy.Velocity * iRing.Energy.Gamma));
      }
      ++iTime;
      ++jtime;
      ++etime;
   }
}

//================================ Tracking ====================================

vectorU f1(xTime&t, vectorU&X)
{
   return iRing.GetForces(t, X);
}

void xDynamics::Tracking()
{  Algoritm = 2;
   iDraw.DrawGammaPi(iBeam);
   iDraw.DrawGamma3(iBeam);
   iBeam.NImpact = 0;
   if (FORCES)
      iTime.DistStep(ds);
   else
      iTime.TimeStep(iRing.Trev);
   Warning("Tracking: Integration step [sec] =", iTime.dt(s_));

   if (iRing.LatticeFile != 1)
   {  Warning("Choose <Ring | Optics structure | MAD input file>",
      PressEnter);
      return;
   }
   if (iRing.ReduceExtend == 0)
   {  Warning("Tracking is not allowed for Reduce Lattice", PressEnter);
      return;
   }

   //xLattice::B_s = iBeam.CellSize;
   iBeam.CellSize = iRing.Circ / iBeam.CellNum;
   Distribution(iBeam.InitialEmit,iRing.GetLattice(0), 0, Crystall);
   matrixU matrix;
   xRFcavity cavity;
   //bpTune Potentials(81, 81, 81);
   double BumpX = 0;
   double BumpY = 0;
   double ShiftX = 0.007;
   double ShiftY = 0.006;
   double septumx = - 0.0059;
   double septumy = - 0.0042; 
   double angle = -1;
   cycles = 0;

   while (!iBeam.inject || (cycles < iBeam.cycles))
   {
      if (iBeam.inject)
      {  
         // bump and septum position 
         for (int i = 0; i < iIBS.Bumpfile.Row(); i++) 
         if (iIBS.Bumpfile[i][0] >= cycles)
         {
            ShiftX = iIBS.Bumpfile[i][1] / 1000.;
            ShiftY = iIBS.Bumpfile[i][2] / 1000.; 
            septumx  = iIBS.Bumpfile[i][3] / 1000.;
            septumy  = iIBS.Bumpfile[i][4] / 1000.; 
            angle = -tan(M_PI*iIBS.Bumpfile[i][5]/180);
            break;
         }
         
         BumpX += ShiftX;
         BumpY += ShiftY;
         double SeptX = BumpX + septumx;
         double SeptY = BumpY + septumy;
   
         for (int i = 0; i < iBeam.Number(); i++)           // position of stored beam
         {
            iBeam[i][0] += ShiftX;
            iBeam[i][2] += ShiftY;
         }

         int num = 0;
         if (cycles) {                                      // no new injection for 1st cycle
            if (iLosses.Multi && iLosses.Injection) 
               for (int i = 0; i < iBeam.Number(); i++)     // particle losses of stored beam above septum
                  if (iBeam[i][2]  > angle * (iBeam[i][0]  - SeptX) + SeptY)  // (-1.) means -45 degree
					       iBeam.Loss(iRing.LATTICE, i--, 1.);
            num = iBeam.Number();                       
		      iBeam.Number(num + iBeam.IniNumber);            // injection
		      //iBeam.Number(num * ((iBeam.InitialEmit[3]/iBeam.Emit[3])(U1_) + 1.));
		      Distribution(iBeam.InitialEmit, iRing.LATTICE, num, Crystall);
            iBeam.Emit[3] += iBeam.InitialEmit[3];
   
            for (int i = num; i < iBeam.Number(); i++)         // position of injected beam
            {
               iBeam[i][0] += BumpX;
               iBeam[i][2] += BumpY;
            }
         }  cycles++;

         if (iLosses.Multi) {
            if (iLosses.Injection) 
               for (int i = num; i < iBeam.Number(); i++)     // particle losses of injected beam below septum
                  if (iBeam[i][2]  < angle * (iBeam[i][0]  - SeptX) + SeptY)   // (-1.) means -45 degree
					       iBeam.Loss(iRing.LATTICE, i--, 1.);
            if (iLosses.Acceptance) {
			      iBeam.Courant_Snyder(iRing.LATTICE);
			      for (int j = 0; j < iBeam.Number(); j++)
				   if (((sqrt(iBeam.Invs[j][2])*fabs(iRing.LATTICE.Dx())+sqrt(iBeam.Invs[j][0]*iRing.LATTICE.betax()))> sqrt(iRing.AcceptHinj()*iRing.LATTICE.betax())) ||
					    ((sqrt(iBeam.Invs[j][2])*fabs(iRing.LATTICE.Dy())+sqrt(iBeam.Invs[j][1]*iRing.LATTICE.betay()))> sqrt(iRing.AcceptVinj()*iRing.LATTICE.betay()))) 
					   // || ((iBeam.benum == BUCKET) && (iBeam.longshift-iRing.S1(m_) < iBeam(j,4)(m_)) && (iBeam(j,4)(m_) < iBeam.longshift+iRing.S1(m_))) )
				      iBeam.Loss(iRing.LATTICE, j--, 1.);
            }
         }
         bpData Septum(1001,2); 
         for (int i = 0; i < 1000; i++)
         {
            Septum[i][0] = -0.2 + 0.001*(i);
            Septum[i][1] = (-1.) * (Septum[i][0] - SeptX) + SeptY;
         }
         Septum.SaveFileSep("septum.cur", "\t");
      }

      if (iBeam.benum == COASTING)
      {  for (int j4 = 0; j4 < iBeam.Number(); j4++)
         {  iBeam[j4][4] -= iBeam.CellSize() * floor(iBeam[j4][4] / iBeam.CellSize());
            iBeam[j4][4] -= iBeam.CellSize() / 2;
         }
      }
      bool traj = false;
      if (Drawing(iRing.LATTICE))
      {  traj = CheckCrossing;
         if (traj)
            iDraw.ResetTrajectory(CurveSize, SkipPoints);
      }
      bool bunch_once = false;
      iBeam.NImpact = 0;
     
      bpData beam4(iBeam.Number(),4);
      if (iIBS.Plottuneshift)
      for (int i = 0; i < iBeam.Number(); i++)
      {  beam4[i][0] = iBeam[i][0] - iRing.LATTICE.Dx()  * iBeam[i][5]; 
         beam4[i][1] = iBeam[i][1] - iRing.LATTICE.Dpx() * iBeam[i][5]; 
         beam4[i][2] = iBeam[i][2] - iRing.LATTICE.Dy()  * iBeam[i][5]; 
         beam4[i][3] = iBeam[i][3] - iRing.LATTICE.Dpy() * iBeam[i][5]; 
      }

      for (int i = 0; i < iRing.Number(); i++)
      {  iTime.so = 0;
         iRing.Index = i;
         if (iBeam.benum == COASTING)
         {  for (int j4 = 0; j4 < iBeam.Number(); j4++)
            {  iBeam[j4][4] -= iBeam.CellSize() * floor(iBeam[j4][4] / iBeam.CellSize());
               iBeam[j4][4] -= iBeam.CellSize() / 2;
            }
         }

         if(bunch_once && iBeam.benum == BUNCHED && iTime.sr >= iRing.Circ / 2.)
         {  bunch_once = false;
            //for (int j2 = 0; j2 < iBeam.Number(); j2++) iBeam.MinusDisp(iRing.GetLattice(i), j2);
            matrix = cavity.GetMatrix(iTime, iRing.Energy);
            for (int j = 0; j < iBeam.Number(); j++)
               iBeam(j, matrix);
            //for (int j2 = 0; j2 < iBeam.Number(); j2++) iBeam.PlusDisp(iRing.GetLattice(i), j2);
         }

         if(FORCES)
            for (int j = 0; j < iBeam.Number(); j++)
               iBeam(j, iRing.OnEnter(iTime, iBeam(j)));

         xRingOptics& opt = iRing[i];

         int Steps = (int)((opt.Length / ds)(U1_) + 0.5);
         if ((Steps < 1) || !FORCES) Steps = 1;
         iTime.DistStep(opt.Length / Steps);

         for (int steps = 0; steps < Steps; steps++)
         {  if (traj) iBeam.KeepCross();

            if(iIBS.Multi == 1) {
               /*
               if (iIBS.TuneShift) {
                  Potentials.ResetGrid(sqrt(iBeam.Emit[0](m_)*opt.Lattice.betax(m_)) / 5, 
                                       sqrt(iBeam.Emit[1](m_)*opt.Lattice.betay(m_)) / 5, 
                                       iBeam.CellSize(m_) / 20);
                  Potentials.SetValue();
                  for (int n = 0; n < iBeam.Number(); n++)
                     Potentials.AddCharge(iBeam[n],iRing.Energy.Gamma(U1_));

                  double K = (((iRing.Energy.Z*U_e)^2) * iTime.ds /
                                 (iRing.Energy.Momentum*iRing.Energy.Velocity*iRing.Energy.Gamma))(m_^2);
                  #pragma omp parallel for
                  for (int n = 0; n < iBeam.Number(); n++) 
                     Potentials.SpaceCharge(iBeam.Bunch[n], K, iBeam.CellSize(m_),iBeam.benum);
               } else */
                  if (iBeam.benum == BUNCHED)
                     iBeam.SpaceChargeBunch(iTime);
                  else
                     iBeam.SpaceChargeMD(iTime);
            }
            if(iHeat.Multi == 1) {
               iHeat.Kick (iTime, iBeam, iRing);
            }
            if (FORCES)
            {  for (int j = 0; j < iBeam.Number(); j++)
                  iBeam(j, Symple(iTime, iTime.ds, iBeam(j), f1));
            }else
            {  if (opt.LGetSize())
               {
                  if ((iEcool.Multi != 0) && (opt[0] == &iEcool))
                  {  for (int j = 0; j < iBeam.Number(); j++)
                        iBeam.Add(j, iEcool.GetForces(iTime,iRing.Energy,iBeam(j))* opt.Length);
                  }else
                  if ((iLaserCooling.Multi != 0) && (opt[0] == &iLaserCooling))
                  {  for (int j = 0; j < iBeam.Number(); j++)
                        iBeam.Add(j, iLaserCooling.GetForces(iTime,iRing.Energy,iBeam(j))*opt.Length);
                  }
               }
               #pragma omp parallel for
               for (int j = 0; j < iBeam.Number(); j++)
                  iBeam(j, opt.Matrix);
            }

            if (traj)
            {  if (iBeam.Crossed()) iDraw.CrossNumber++;
               iDraw.Trajectory(iTime, iBeam, AutoSkip);
            }
            ++iTime;
         }
         if (FORCES)
            for (int k = 0; k < iBeam.Number(); k++)
                iBeam(k, iRing.OnExit(iTime, iBeam(k)));
      }
      if (iRing.Inductor)
      for (int j1 = 0; j1 < iBeam.Number(); j1++)
            iBeam.Add(j1, 5, U_e * iRing.Vind /
            (iRing.Energy.Momentum * iRing.Energy.Velocity * iRing.Energy.Gamma));
      if (traj)
         iDraw.SaveTrajectory();
      if (iIBS.Multi == -1)
      {
         iTime.TimeStep(iRing.Trev);
         iIBS.Kick(iTime, iBeam, iRing);
      }  

// Tune Shift Diagram
      if (iIBS.Plottuneshift) {
         double Qx_av=0;
         double Qy_av=0;
         xLattice&lat = iRing.LATTICE;
         bpData PhaseAdvance2D(iBeam.Number()+1, 2);
         bpData PhaseShift(10,2);
         PhaseShift[0][0] = lat.Qx();
         PhaseShift[0][1] = lat.Qy();
         for (int i = 0; i < iBeam.Number(); i++)
         {  double ion[4];
            ion[0] = iBeam[i][0] - lat.Dx()  * iBeam[i][5]; 
            ion[1] = iBeam[i][1] - lat.Dpx() * iBeam[i][5]; 
            ion[2] = iBeam[i][2] - lat.Dy()  * iBeam[i][5]; 
            ion[3] = iBeam[i][3] - lat.Dpy() * iBeam[i][5]; 

            double sum1 = beam4[i][0]*lat.gammax() + beam4[i][1]*lat.alfax();
            double sum2 = beam4[i][0]*lat.alfax()  + beam4[i][1]*lat.betax();
            double cos_x = (ion[0]*sum1 + ion[1]*sum2) / (beam4[i][0]*sum1 + beam4[i][1]*sum2);

            double sin_x = ((ion[0] * beam4[i][1] - ion[1] * beam4[i][0]) /
            (beam4[i][0]*beam4[i][0]*lat.gammax() + 2*beam4[i][0]*beam4[i][1]*lat.alfax() + beam4[i][1]*beam4[i][1]*lat.betax()));
            double Q_x = atan(sin_x/cos_x) / (2.*M_PI);
            if (Q_x < iIBS.TSminX) Q_x += 0.5;
            //while (Q_x > lat.Qx()) Q_x -= 0.5;
            //while (Q_x < lat.Qx() - 0.4999) Q_x += 0.5;

            sum1 = beam4[i][2]*lat.gammay() + beam4[i][3]*lat.alfay();
            sum2 = beam4[i][2]*lat.alfay()  + beam4[i][3]*lat.betay();
            double cos_y = (ion[2]*sum1 + ion[3]*sum2) / (beam4[i][2]*sum1 + beam4[i][3]*sum2);

            double sin_y = ((ion[2] * beam4[i][3] - ion[3] * beam4[i][2]) /
            (beam4[i][2]*beam4[i][2]*lat.gammay() + 2*beam4[i][2]*beam4[i][3]*lat.alfay() + beam4[i][3]*beam4[i][3]*lat.betay()));
            double Q_y = atan(sin_y/cos_y) / (2.*M_PI);
            if (Q_y < iIBS.TSminY) Q_y += 0.5;
            //while (Q_y > lat.Qy()) Q_y -= 0.5;
            //while (Q_y < lat.Qy() - 0.4999) Q_y += 0.5;

            PhaseAdvance2D[i][0] = Q_x;
            PhaseAdvance2D[i][1] = Q_y;
            Qx_av += Q_x;
            Qy_av += Q_y;
         }
         PhaseAdvance2D.SaveFileSep("tuneshift.cur", "\t");
         PhaseShift[1][0] = Qx_av / iBeam.Number();
         PhaseShift[1][1] = Qy_av / iBeam.Number();
         PhaseShift[3][0] = iIBS.TSminX;
         PhaseShift[3][1] = iIBS.TSminY;
         PhaseShift[4][0] = iIBS.TSminX+0.5;
         PhaseShift[4][1] = iIBS.TSminY;
         PhaseShift[5][0] = iIBS.TSminX+0.5;
         PhaseShift[5][1] = iIBS.TSminY+0.5;
         PhaseShift[6][0] = iIBS.TSminX;
         PhaseShift[6][1] = iIBS.TSminY+0.5;
         PhaseShift[7][0] = iIBS.TSminX;
         PhaseShift[7][1] = iIBS.TSminY;
         PhaseShift.SaveFileSep("phaseshift.cur", "\t");
      }
   }
}

//---------------------------------------------------------------------------

void xDynamics::Algorithms()
{
   iDraw.DrawGammaPi(iBeam);
   iDraw.DrawGamma3(iBeam);

   switch(Algoritm)
   {
    case 0:
      Warning("Dynamics: Initial step [sec] =", dt(), ", Step multiplier =",ifmore,
              ", Max growth =", ratio);
      Dynamics();
      break;
    case 1:
      Warning("Model Beam: Integration Step [sec] =", IntegrationStep());
      ModelBeam();
      break;
    case 2:
      Warning("Tracking: Distance step [cm] =", ds());
      Tracking();
      break;
   }
}

//---------------------------------------------------------------------------

void xDynamics::Lattice()
{
   if (iRing.LatticeFile == 2)
      Warning("Choose optics file", PressEnter);
   else
      iDraw.LatticeFunctions(iRing, iTime);
}

void xDynamics::Rates()
{
   if ((iRing.LatticeFile == 2) && iIBS.Model)
   {  Warning("<Ring|Optics structure|No optics file> - only for <IBS|Piwiski>", PressEnter);
      return;
   }

   vectorU rates(s_^-1);
   rates = xEffect::Summary(iTime, iBeam, iRing);
   iTime.DistStep(ds);

   for (int i = 0; i < 4; i++)
      Data[32][i+1] = rates[i]();
}
//---------------------------------------------------------------------------
void xDynamics::BeamTest()
{
   xLattice latt;
   latt.betax = 1;
   latt.alfax = 3;
   latt.Dx = 10;
   latt.Dpx = 10;
   latt.betay = 1;
   latt.alfay = 3;

   iBeam.Number(5000);

   doubleU B_s(100,m_);
   //bool bunched = true;
/*
        iBeam.RMS(latt, B_s, bunched);
   iBeam.Round(bunched);
   iBeam.beam_matching_wD(latt);
   iBeam.beam_MO();
   iBeam.beam_sigma();
   iBeam.beam_correlation();
   iBeam.Emit = iBeam.Emitt();
   iBeam.beam_D(latt);
   Warning("Ex=",iBeam.Emit[0](),"Ey=",iBeam.Emit[1](),
           "dP/P=",sqrt(iBeam.Emit[2]()));
*/
   xDistributor beam(5000);
   beam.Emit[0] = iBeam.Emit[0];
   beam.Emit[1] = iBeam.Emit[1];
   beam.Emit[2] = iBeam.Emit[2];
   beam.Distribution(beam.Emit,latt);

// from 22.04.04
//   beam.Emit = beam.Emittance(latt);
   beam.Emit = beam.Emit_Def(latt);

   Warning("Ex=",beam.Emit[0](),"Ey=",beam.Emit[1](),
           "dP/P=",sqrt(beam.Emit[2]()));

}
//-----------------------------------------------------
void xDynamics::FFTest()
{
  iDraw.FF(iForce, iTime, iEbeam, iRing);
}

//-----------------------------------------------------
#include "xpowell.h"

void xDynamics::LuminosityCalculation()
{
   //xEffect* effect = xEffect::GetEffect("IBS");
   //if (effect)
   //   Warning("IBS betax =", effect->Lattice.betax());
   iColl.CalcLuminosity(iTime, iBeam, iRing);
}

