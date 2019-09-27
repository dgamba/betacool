//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xdraw.h"
#include "xlibrary.h"
#include <omp.h>
#include <iostream>

xDraw::xDraw()
{
   Space.SetName("space.cur");
   Space.Size(1000,6);
   Emittance.SetName("emittance.cur");
   Emittance.Size(1000,10);
   Evolut.SetName("evolution.cur");
   Evolut.Size(1000,7);
   Rates1.SetName("rates.cur");
   Rates1.Size(1000,5);
   Lattice.SetName("lattice.cur");
   Lattice.Size(1000,11);
   Barrier.SetName("barrier.cur");
   Barrier.Size(1000,7);
   LongProfile.SetName("longprof.cur");
   LongProfile.Size(1000,4);
   GSystems.SetName("gated.cur");
   GSystems.Size(1000,2);
/*
   ElectronShiftH[0].SetName("shift_xh.cur");
   ElectronShiftH[1].SetName("shiftpxh.cur");
   ElectronShiftH[2].SetName("shift_yh.cur");
   ElectronShiftH[3].SetName("shiftpxh.cur");
   ElectronShiftH[4].SetName("shift_sh.cur");
   ElectronShiftH[5].SetName("shiftdph.cur");
   ElectronShiftV[0].SetName("shift_xv.cur");
   ElectronShiftV[1].SetName("shiftpxv.cur");
   ElectronShiftV[2].SetName("shift_yv.cur");
   ElectronShiftV[3].SetName("shiftpxv.cur");
   ElectronShiftV[4].SetName("shift_sv.cur");
   ElectronShiftV[5].SetName("shiftdpv.cur");
*/
   Coordinate.SetName("coordinate.cur");
   Coordinate.Size(1000,4);
   Profile.SetName("profile.cur");
   Profile.Size(1000,5);
   Invariant.SetName("invariant.cur");
   Invariant.Size(1000,4);
   Evolution. SetName("evolution.sur");
   evolution = false;
   EvolutionAver = false;
   BiGaussian = false;
   SumBiGaussian = false;
   GFitD = false;
   GFitE = false;

   SpaceCharge.SetName("charge.cur");
   DriftVelocity.SetName("vdrift.cur");
   InEcool.SetName("inecool.cur");

   BeamTxy.SetName("txy2t.cur");
   BeamPrn.SetName("footprn.cur");
   Gamma2.SetName("gamma2.cur");
   Gamma3.SetName("gamma3.cur");
   Emittances = false;
   Spread = false;

   TrajX.SetName("trajX.cur");
   TrajY.SetName("trajY.cur");
   Cross.SetName("cross.cur");
   Target.SetName("target.cur");

   RateEh.SetName("rateeh.sur");
   RateEv.SetName("rateev.sur");
   RatedP.SetName("ratedp.sur");

   FFtr.SetName("fftr.sur");
   FFlong.SetName("fflong.sur");
   FCtr.SetName("ftr.cur");
   FClong.SetName("flong.cur");
   FAtr.SetName("fatr.cur");
   FAlong.SetName("falong.cur");

   LuminosityTest.SetName("lumitest.cur");
   Luminosity3D.SetName("lumi3d.sur");
   LaserForce.SetName("laserf.cur");
   LaserLight.SetName("laserl.cur");

   CrossNumber = 0;

   EBunch.SetName("ebunchl.cur");
   EDensity.SetName("density.sur");

//17.12.2007
   StochForce.SetName("stochrhic.cur");
   StochDiff.SetName("diffrhic.cur");

   Pellet.SetName("pellet.cur");
   Average.SetName("average.cur");
   Pellets.SetName("pellets.cur");
   BeamSize.SetName("beamsize.cur");
   BeamCros.SetName("beamcros.cur");

   Potentials.SetName("potentials.sur");
   ElectricField.SetName("electricfield.sur");

}

int xDraw::OnGet()
{
   Evolution.Enabled = Data.OnGetB(30,1,false);
   Sigmas         = Data.OnGet(30,2,5);
   Divisions      = Data.OnGetI(30,3,50);
   EvolutionFor   = Data.OnGetI(30,4,0);
   NormalisedOn   = Data.OnGetI(30,5,1);
   EvolutionSlices= Data.OnGetI(30,6,10);
   EvolutionStep  = Data.OnGet(30,7,1);
   Coordinate.Enabled = Data.OnGetB(30,8,false);
   Profile.Enabled = Data.OnGetB(30,9,false);
   Invariant.Enabled = Data.OnGetB(30,10,false);
   EvolutionAver = Data.OnGetB(30,11,false);
   GFitD = Data.OnGetB(30,12,false);
   //SumBiGaussian = Data.OnGetB(30,19,false);
   if (Evolution.Enabled)
   {
      Evolution.SetXaxis(Divisions, -Sigmas, Sigmas, 0);
      Evolution.SetYaxis(EvolutionSlices, 0, EvolutionStep*EvolutionSlices, 0);
   }

   Emittances = Data.OnGetB(29,1,false);
   Spread = Data.OnGetB(29,2,false);
   Gamma2.Enabled = Data.OnGetB(29,3,false);
   Gamma3.Enabled = Data.OnGetB(29,4,false);
   BeamPrn.Enabled = Data.OnGetB(29,5,false);
   BeamTxy.Enabled = Data.OnGetB(29,6,false);
   GFitE = Data.OnGetB(29,7,false);

   Rates1.Enabled = Data.OnGetB(32,8,false);

   InEcool.Enabled = Data.OnGetB(62,6,false);

   xGraf::SaveAllGraf();

   return 0;
}

int xDraw::OnRun()
{

   //Evolution.Enabled = (bool)Data[30][1];
   Sigmas         = Data.OnGet(30,2,5);
   Divisions      = Data.OnGetI(30,3,50);
   //EvolutionFor   = Data.OnGet(30,4,0);
   //NormalisedOn   = Data.OnGet(30,5,1);
   //EvolutionSlices= Data.OnGet(30,6,10);
   //EvolutionStep  = Data.OnGet(30,7,1);
   Coordinate.Enabled = Data.OnGetB(30,8,false);
   Profile.Enabled = Data.OnGetB(30,9,false);
   Invariant.Enabled = Data.OnGetB(30,10,false);
   GFitD = Data.OnGetB(30,12,false);
/*
   if (Evolution.Enabled)
   {
      Evolution.SetXaxis(Divisions, -Sigmas, Sigmas, 0);
      Evolution.SetYaxis(EvolutionSlices, 0, EvolutionStep*EvolutionSlices, 0);
   }
*/
   Gamma2.Enabled = Data.OnGetB(29,3,false);
   Gamma3.Enabled = Data.OnGetB(29,4,false);
   BeamPrn.Enabled = Data.OnGetB(29,5,false);
   BeamTxy.Enabled = Data.OnGetB(29,6,false);
   GFitE = Data.OnGetB(29,7,false);

   Rates1.Enabled = Data.OnGetB(32,8,false);

   InEcool.Enabled = Data.OnGetB(62,6,false);

   return 0;
}

void xDraw::SaveTrajectory()
{
   TrajX.Save();
   TrajY.Save();
   Cross.Save();
}

void xDraw::ResetTrajectory(int CurveSize, int SkipPoint)
{
   TrajX.Reset(SkipPoint);
   TrajY.Reset(SkipPoint);
   MaximumPoints = CurveSize * SkipPoint;
   CurrentNumber = 0;
   CrossNumber = 0;
}

void xDraw::Trajectory(xTime& time, xBeam& beam, bool askip)
{
   if ((CurrentNumber < MaximumPoints) || askip)
   {
      TrajX.Point(time.sr(), beam[0][0]);
      TrajY.Point(time.sr(), beam[0][2]);
      //Cross.Point(time.sr(), CrossNumber);
      CurrentNumber++;
   }
}

void xDraw::DrawGammaPi(xBeam& beam)
{
   for (double i= -200; i < 1; i++)
   {  doubleU momentum(pow(10, i/10.), U1_);
      Gamma2.Point(momentum(),beam.GetGamma2(momentum)(u_6*m_));
   }
   Gamma2.Save();
}

void xDraw::DrawGamma3(xBeam& beam)
{  /*
   for (double i= -200; i <= 100; i++)
   {  doubleU emit(pow(10, i/10.), m_);
      BeamTxy.Point(beam.GetGamma3(emit)(U1_),emit(m_));
   }
   BeamTxy.Save();
   BeamTxy.Reset(0);
   */
   for (double i= -200; i <= 100; i++)
   {  doubleU emit(pow(10, i/10.), m_);
      Gamma3.Point(beam.GetEquillibrium(emit)(U1_),emit(u_6*m_));
   }
   Gamma3.Save();
}
//---------------------------------------------------------------------------

bool xDraw::DrawRealSpace(xBeam& beam, xRing& ring, xLattice& Lattice, int period)
{
   bool cansave = SpaceTimer.Timer(period);
   
   if(evolution || cansave)
   {
      double emit[4];
      if(NormalisedOn)
      {  emit[0] =((beam.InitialEmit[0]*Lattice.betax)^0.5)();
         emit[1] =((beam.InitialEmit[1]*Lattice.betay)^0.5)();
         emit[2] = (beam.InitialEmit[2]               ^0.5)();
         emit[3] = (beam.Emit[3] / beam.InitialEmit[3])();
      }else
      {  emit[0] =((beam.Emit[0]*Lattice.betax)^0.5)();
         emit[1] =((beam.Emit[1]*Lattice.betay)^0.5)();
         emit[2] = (beam.Emit[2]               ^0.5)();
         emit[3] = 1;
      }

      if (iBeam.benum == BUNCHED) 
      {  Barrier.Reset(0);
         double X[7];
         for (int i = 0; i < 6; i++)
            X[i] = 0;
         
         if (iDynamics.Synchrotron == 2 && iDynamics.LongSC) 
         {  if (iDynamics.LongOK)
            {  doubleU s_s((xLattice::B_s) * (iBeam.InitialEmit[2]^0.5));
               for (int j = 0; j < Divisions; j++)
               {  X[0] = iDynamics.LongHyst[0][j] * s_s(m_);
                  X[3] = iDynamics.LongHyst[3][j] / 1000.;
                  X[5]-= X[3] * M_PI * Sigmas / Divisions / 2.;
                  Barrier.Points(X, 7);
               }
            }
         } else
         {
            int r = 200;
            int a = iRing.h() * r / 2;  // tottal ring
            //int a = r / 2;                // one harmonic
            for (int j = -a; j <= a; j++)
            {  X[0] = iRing.Circ() * j / (iRing.h() * r);
               X[3] = - iRing.V() * sin(X[0] * iRing.h() * 2 * M_PI / iRing.Circ());
               if (iRing.NonLinear) 
                   X[3] -= iRing.V_low()* sin(X[0] * iRing.h_low()* 2 * M_PI / iRing.Circ());
               X[5]-= X[3] * M_PI / r;
               Barrier.Points(X, 7);
            }
         }
      }
// Barrier bucket
	  if (beam.benum == BUCKET)
	  {
		  //ring.Bucket.CalcParticle(beam, ring);
		  if (iRing.Bucket.Analytic && (iRing.Bucket.Stationary || iRing.Bucket.Moving))
		  {
			  complex<double> Ah(0, 0), Z(260, 2.64e4), Uind(0, 0);
			  if (iRing.Bucket.Induction)
			  {
				  for (int j = 0; j < iBeam.Number(); j++)
				  {
					  complex<double>ih(0, (iRing.h_low * iBeam(j, 4) * 2 * M_PI / iRing.Circ)(U1_));
					  Ah += exp(-ih);
				  }
				  Ah /= iBeam.Number();
				  doubleU Io(iRing.Energy.Z * U_e * iBeam.Emit[3] / iRing.Trev);
				  double I = Io(A_);
				  Uind = -Ah * iRing.Bucket.Impedance * I;
			  }
			  Barrier.Reset(0);
			  double X[7];
			  for (int i = 0; i < 6; i++)
				  X[i] = 0;
			  int r = 400;
			  int a = r / 2;                
//***********************************************************************************************************
			  double WWmax = 0;
			  for (int j = -a; j <= a; j++)
			  {
				  X[0] = iRing.Circ(m_) * j / (iRing.h() * r);
				  doubleU V(0, V_);
				  if (iRing.Bucket.Induction)
				  {
					  complex<double>ih(0, (iRing.h_low * X[0] * 2 * M_PI / iRing.Circ)(m_^-1));
					  ih = Uind * (exp(ih));
					  V = 2. * ih.real();
				  }
				  X[3] = -iRing.V(k_3*V_)    * sin(X[0] * iRing.h() * 2 * M_PI / iRing.Circ())
					  - iRing.V_low(k_3*V_)* sin(X[0] * iRing.h_low() * 2 * M_PI / iRing.Circ()) - V(k_3*V_);
				  X[5] -= X[3] * iRing.Circ(m_) / (iRing.h() * r);
				  if ((j > -a) && (X[5] > WWmax)) WWmax = X[5];
			  }
//---
				for (int i = 0; i < 6; i++)
				  X[i] = 0;
			  for (int j = -a; j <= a; j++)
			  {
				  X[0] = iRing.Circ(m_) * j / (iRing.h() * r);
				  doubleU V(0, V_);
				  if (iRing.Bucket.Induction)
				  {
					  complex<double>ih(0, (iRing.h_low * X[0] * 2 * M_PI / iRing.Circ)(m_^-1));
					  ih = Uind * (exp(ih));
					  V = 2. * ih.real();
				  }
				  X[3] = -iRing.V(k_3*V_)    * sin(X[0] * iRing.h()     * 2 * M_PI / iRing.Circ())
					     -iRing.V_low(k_3*V_)* sin(X[0] * iRing.h_low() * 2 * M_PI / iRing.Circ()) - V(k_3*V_);
				  X[5] -= X[3] * iRing.Circ(m_) / (iRing.h() * r);
				  X[1] = -sqrt(fabs((iRing.Energy.Z * U_e *(WWmax - X[5]) * 2 / 
					  (iRing.Circ * iRing.Energy.Energy * iRing.Eta * iRing.Energy.Beta2))((k_3*V_*m_)^-1)));
				  X[4] = -X[1];
				  Barrier.Points(X, 7);
			  }
//**********************************************************************************************************
		  }
		  else
		  {
			  double X[7];
			  Barrier.Reset(0);
			  for (int i = 0; i < ring.Bucket.barrier.Number; i++)
			  {
				  X[0] = ring.Bucket.barrier[i].s1();
				  X[1] = ring.Bucket.barrier[i].H();
				  X[2] = ring.Bucket.barrier[i].D;
				  X[3] = (ring.Bucket.barrier[i].Ub)(k_3*V_);// +ring.Bucket.barrier[i].Us)(k_3*V_);
				  X[4] = sqrt(ring.Bucket.barrier[i].P);
				  if (i)
				  {
					  X[5] = ring.Bucket.barrier[i - 1].UP(k_3*V_*m_);
					  if (ring.Bucket.bCoeff && ring.Bucket.UPcoeff)
						  X[5] *= ring.Bucket.Ucoeff / ring.Bucket.UPcoeff;
				  }
				  else
					  X[5] = 0;
				  //X[6] = ring.Bucket.barrier[i].I;
				  X[6] = ring.Bucket.barrier[i].Us(k_3*V_);


				  Barrier.Points(X, 7);

				  X[0] = ring.Bucket.barrier[i].s2();
				  X[5] = ring.Bucket.barrier[i].UP(k_3*V_*m_);
				  if (ring.Bucket.bCoeff && ring.Bucket.UPcoeff)
					  X[5] *= ring.Bucket.Ucoeff / ring.Bucket.UPcoeff;
				  Barrier.Points(X, 7);
			  }
			  /* Remove Shen modofication
				double X[7];
				Barrier.Reset(0);
				for (int i = 0; i < ring.Bucket.barrier.Number; i++)
				{
				doubleU sstep(0,m_);
				sstep = (ring.Bucket.barrier[i].s2-ring.Bucket.barrier[i].s1)/20;
				X[1] = ring.Bucket.barrier[i].H();
				X[2] = ring.Bucket.barrier[i].D;
				X[4] = sqrt(ring.Bucket.barrier[i].P);
				X[6] = ring.Bucket.barrier[i].I;
				for (int j=0;j<21;j++)
				{
				X[0] = (ring.Bucket.barrier[i].s1 + sstep*j)();
				X[3] = ring.Bucket.barrier[i].U(k_3*V_)*(ring.Bucket.Issin?sin(j/20.0*M_PI):1);
				X[5] = (i?ring.Bucket.barrier[i-1].UP(k_3*V_*m_):0.0)*(cos(j/20.0*M_PI)+1)/2.0 + ring.Bucket.barrier[i].UP(k_3*V_*m_)*(1-cos(j/20.0*M_PI))/2.0;
				if (ring.Bucket.bCoeff && ring.Bucket.UPcoeff)
				X[5] *= ring.Bucket.Ucoeff / ring.Bucket.UPcoeff;
				Barrier.Points(X, 7);
				if((ring.Bucket.barrier[i].U() == 0 || !ring.Bucket.Issin) && !j) j=19;
				}
				}*/
		  }
	  }
// Longitudinal profile      
      LongProfile.Reset(0);
      double X[4];
      double circ;
      if (iBeam.benum == BUNCHED)
         circ = iRing.L_s(m_);
      else
         circ = iRing.Circ(m_);
      vectorU rates(s_^-1);
      vectorU emittance;
      emittance = iBeam.Emit;
      for (int i = 0; i < iBeam.LProfile.Number; i++)
      {
         X[0] = circ*i/beam.LProfile.Number-circ/2;
         X[1] = beam.LProfile[i];
         X[2] = beam.VProfile[i];
         //X[2] = beam.RatesDp[i];
         X[3] = sqrt(beam.Momentum[i]);
         LongProfile.Points(X, 4);
         /*
         if (i == beam.LProfile.Number-1)
         {
            X[1] = beam.LProfile[0];
            X[3] = sqrt(beam.Momentum[0]);
         }else
         {
            X[1] = beam.LProfile[i+1];
            X[2] = beam.Integral[i+1];
            X[3] = sqrt(beam.Momentum[i+1]);
         }
         LongProfile.Points(X, 4);
         */
      } 
      iBeam.Emit = emittance;

/*
         for (int i = 0; i < ring.Bucket.Profile.Number-1; i++)
         {
            LongProfile.Point((ring.Circ*i/ring.Bucket.Profile.Number-ring.Circ/2)(m_),
               ring.Bucket.Profile[i]);
            LongProfile.Point((ring.Circ*i/ring.Bucket.Profile.Number-ring.Circ/2)(m_),
               ring.Bucket.Profile[i+1]);
         }
         LongProfile.Point((ring.Circ/2)(m_), ring.Bucket.Profile[ring.Bucket.Profile.Number-1]);
         LongProfile.Point((ring.Circ/2)(m_), ring.Bucket.Profile[0]);
*/
	  
      GSystems.Reset(0);
      if (iGated.Multi)
      {
         for (int i = 0; i < iGated.Number; i++)
         if (iGated.System[i].use)
         {
            GSystems.Point((iGated.System[i].initial * ring.Circ)(m_), 0.25*(i+1.)*ring.Bucket.Ucoeff);
            GSystems.Point((iGated.System[i].final * ring.Circ)(m_), 0.25*(i+1.)*ring.Bucket.Ucoeff);
            GSystems.Point();
         }
      }
      
//Electron beam shifts
/*
      if (iEbeam.EnableShifts)
      {
         for (int f = 0; f < 6; f++)
         {
            ElectronShiftH[f].Size(10);
            ElectronShiftH[f].Enabled = false;
            ElectronShiftV[f].Size(10);
            ElectronShiftV[f].Enabled = false;
         }

         for (int n = 0; n < 4; n++)
         if(RealSpace[n].Graf.Enabled)
         {
            int xa = RealSpace[n].x_axis;
            int ya = RealSpace[n].y_axis;

            ElectronShiftH[xa].Enabled = true;
            ElectronShiftV[ya].Enabled = true;

            ElectronShiftH[xa].Point(iEbeam.Shift[xa](), -1e6);
            ElectronShiftH[xa].Point(iEbeam.Shift[xa](),  1e6);

            ElectronShiftV[ya].Point(-1e6, iEbeam.Shift[ya]());
            ElectronShiftV[ya].Point( 1e6, iEbeam.Shift[ya]());
         }
      }
*/
	  
// Coordinate distributions
      Tray<double>Hyst(6,Divisions+1);
      Space.Size(beam.Number()+1,6);
      if (Coordinate.Enabled)
      {
         Coordinate.Size(Divisions*2+5,4);
         for (int n = 0; n < 3; n++)
         for (int k = 0; k <= Divisions; k++)
            Hyst[n][k] = 0;
      }
      double h = 2 * Sigmas / Divisions;

      int index, coor;
      for (int i = 0; i < beam.Number(); i++)
      {
         double X[6];
         for (int j = 0; j < 6; j++)
            X[j] = beam[i][j];
         Space.Points(X, 6);

         if(Coordinate.Enabled)
         {
            for (int n4 = 0; n4 < 3; n4++)
            {
               switch(n4)
               {case  0: coor = 0; break;
                case  1: coor = 2; break;
                default: coor = 5;
               }
               index=Round( (beam[i][coor] / emit[n4]+Sigmas) / h );
               if (index < 0) index = 0;
               if (index > Divisions) index = Divisions;
               Hyst[n4][index] += 1;
            }
         }
      }
	  
      if(Coordinate.Enabled)
      {
         double Y[4];
         for (int j = 0; j <= Divisions; j++)
         {  Y[0] = -Sigmas + h*j;
            for (int n = 0; n < 3; n++)
               Y[n+1] = Divisions*Hyst[n][j]/Sigmas/beam.Number() * emit[3];
            Coordinate.Points(Y, 4);
            if (j < Divisions-1)
            {
               Y[0] = -Sigmas + h*(j+1);
               Coordinate.Points(Y, 4);
            }   
         }
      }
	  
// Profiles
      if (Profile.Enabled)
      {
         Profile.Size(Divisions*2+10,5);
         beam.Amplitudes(Lattice);
         for (int n5 = 3; n5 < 6; n5++)
            beam.Profile(Hyst[n5], n5, emit[n5-3], Divisions);

         double X[5], Xo;
         for (int j = 0; j <= Divisions; j++)
         {  X[0] = -Sigmas + h*j;
            for (int n6 = 3; n6 < 6; n6++)
               X[n6-2] = Divisions*Hyst[n6][j]/Sigmas/beam.Number() * emit[3];
            if (j)
               X[4] = (Xo - X[3]) * 2. * M_PI;
            else
               X[4] = 0;
            Xo = X[3];
            Profile.Points(X, 5);
         }
         Profile.Point();
		 
         if (GFitE || GFitD)
         {
            for (int n6 = 3; n6 < 6; n6++)
            {  beam.Powel1(beam.Fitting[n6-3], Hyst[n6], beam.Sigmas, Divisions);
               beam.Fitting[n6-3][0] = beam.Fitting[n6-3][1] * beam.Fitting[n6-3][1];
               if(NormalisedOn)
                  beam.Fitting[n6-3][0] /= (beam.Emit[n6-3] / beam.InitialEmit[n6-3])(U1_);
            }
            if (GFitD)
            for (int j = 0; j <= Divisions; j++)
            {  X[0] = -Sigmas + h*j;
               for (int n6 = 3; n6 < 6; n6++)
                  X[n6-2] = Divisions * emit[3] / Sigmas / beam.Number() *
                    beam.Fitting[n6-3][2] * exp(-0.5 * pow(X[0] / beam.Fitting[n6-3][1], 2.));
               Profile.Points(X, 4);
            }
         }
/*
            Profile[i].Point();
            for (int j = 0; j <= Divisions; j++)
            {  x = -Sigmas + h*j;
               Distribution[n6].Point(x, 100. / beam.Number() *
                 beam.Fitting[n6-3][3] * exp(-0.5 * pow(x / beam.Fitting[n6-3][1], 2.)) );
            }
            Distribution[n6].Point();
            for (int j = 0; j <= Divisions; j++)
            {  x = -Sigmas + h*j;
               Distribution[n6].Point(x, 100. / beam.Number() *
                 beam.Fitting[n6-3][4] * exp(-0.5 * pow(x / beam.Fitting[n6-3][2], 2.)) );
            }
         }
*/

/*
            if (BiGaussian && iIBS.IBSkickModel == 1)
            {
               Distribution[n6].Point();
               if (SumBiGaussian)
               for (int j = 0; j <= Divisions; j++)
               {  x = -Sigmas + h*j;
                  Distribution[n6].Point(x, 100. / beam.Number() *
                  ( beam.Fitting[n6-3][3] * exp(-0.5 * pow(x / beam.Fitting[n6-3][1], 2.)) +
                    beam.Fitting[n6-3][4] * exp(-0.5 * pow(x / beam.Fitting[n6-3][2], 2.)) ));
               }else
               {
                  for (int j = 0; j <= Divisions; j++)
                  {  x = -Sigmas + h*j;
                     Distribution[n6].Point(x, 100. / beam.Number() *
                       beam.Fitting[n6-3][3] * exp(-0.5 * pow(x / beam.Fitting[n6-3][1], 2.)) );
                  }
                  Distribution[n6].Point();
                  for (int j = 0; j <= Divisions; j++)
                  {  x = -Sigmas + h*j;
                     Distribution[n6].Point(x, 100. / beam.Number() *
                       beam.Fitting[n6-3][4] * exp(-0.5 * pow(x / beam.Fitting[n6-3][2], 2.)) );
                  }
               }
            }
         }
*/
      }
	  
// Invariant
      if(Invariant.Enabled)
      {
         Invariant.Reset(0);
         beam.EmitSort(100, 100, Lattice);
         double X[4];
         for (int i = 0; i < beam.Number(); i++)
         {  X[3] = 100.*i/beam.Number();
            for (int n = 0; n < 3; n ++)
               X[n]   = beam.InvI(i)[n];
            Invariant.Points(X, 4);
         }		 
      }
	  
// Footprint
      if(BeamPrn.Enabled)
      {  BeamPrn.Reset(0);
         beam.Courant_Snyder(Lattice);
         for (int i = 0; i < beam.Number(); i++)
            BeamPrn.Point(sqrt(beam.InvI(i)[2]), beam.InvI(i)[0]*1e6);
         BeamPrn.Save();
      }
	  
// 3D Evolution
      static int Slices = 0;
      static int AvN = 0;
      static bool canevolution = false;
      if ((Evolution.Enabled) && (Slices < EvolutionSlices))
      {
        if(( iTime.t() - Slices * EvolutionStep < EvolutionStep * -1e-9) && EvolutionAver)
           {
              AvN++;
              for (int j = 1; j <= Divisions; j++)
                   beam.CurrentAv[j] += Divisions*Hyst[EvolutionFor][j-1]/Sigmas/beam.Number()*emit[3];
           }
        if (iTime.t() - Slices * EvolutionStep >= EvolutionStep * -1e-9)
           {
              canevolution = true;
              Slices++;
              Evolution.SetPoint(0, Slices, iTime.t());
              if(EvolutionAver)
                {
                if(AvN == 0) AvN = 1;
                for (int j = 1; j <= Divisions; j++)
                   {
                   Evolution.SetPoint(j, Slices, beam.CurrentAv[j]/AvN);
                   beam.CurrentAv[j] = 0;
                   }
                AvN = 0;
                }
              else
                {
                for (int j = 1; j <= Divisions; j++)
                   Evolution.SetPoint(j, Slices, Divisions*Hyst[EvolutionFor][j-1]/Sigmas/beam.Number()*emit[3]);
                }
           }
      }	  
      if (cansave)
      {  Space.Save();
         InEcool.Save();
         LongProfile.Save();
         if (beam.benum == BUCKET)
            Barrier.Save();
         if (beam.benum == BUNCHED && (iRing.Bucket.Stationary || iRing.Bucket.Moving))
            Barrier.Save();
         GSystems.Save();
         Coordinate.Save();
         Profile.Save();
         Invariant.Save();
         if(canevolution)
         {  canevolution = false;
            Evolution.Save();
         }
      }
      evolution = false;
      return true;
	}
   return false;
}

bool xDraw::DrawBeamEvolution(xBeam& beam,xTime& itime,xColl&coll,int period)
{
   double X[10];
   X[0] = itime.t(s_);

   if (Emittances)
   {  X[1] = (beam.Emit[0]* iRing.Energy.Beta*iRing.Energy.Gamma*beam.Nrms)(u_6*m_);
      X[2] = (beam.Emit[1]* iRing.Energy.Beta*iRing.Energy.Gamma*beam.Nrms)(u_6*m_);
   }else
   {  X[1] = beam.Emit[0](u_6*m_);
      X[2] = beam.Emit[1](u_6*m_);
   }
   if (Spread)
   {  if (beam.iMomentum < 4)
      {  X[3] = ((beam.Emit[2]^0.5)*iRing.Energy.Momentum)(beam.uMomentum);
         X[4] = ( beam.Emit[5]     *iRing.Energy.Momentum)(beam.uMomentum);
      }else
      {  X[3] = ((beam.Emit[2]^0.5)*iRing.Energy.Kinetic *
                   (iRing.Energy.Gamma + 1) / iRing.Energy.Gamma)(beam.uMomentum);
         X[4] = ( beam.Emit[5]     *iRing.Energy.Kinetic *
                   (iRing.Energy.Gamma + 1) / iRing.Energy.Gamma)(beam.uMomentum);
      }
   }else
   {  X[3] = (beam.Emit[2]^0.5)();
      X[4] = beam.Emit[5]();
   }
   X[5] = beam.Emit[3]();       
   if (GFitE)
   {  X[6] = X[1] * beam.Fitting[0][0];
      X[7] = X[2] * beam.Fitting[1][0];
   }else
   {  X[6] = EmptyData;
      X[7] = EmptyData;
   }                               
   if (iBeam.inject)
      X[8] = iBeam.Emit[3](U1_) * iBeam.cycles / iDynamics.cycles;
   else
      X[8] = coll.Events(U1_);
   X[9] = beam.Emit[4](u_6*m_);
   Emittance.Points(X, 10);

   if (itime.t() > 0)
   {  if (iLaserCooling.Multi)
         X[1] = iLaserCooling.TotalLight / iBeam.Number();
      else
         X[1] = coll.Luminosity((cm_^-2)/s_);
      X[2] = coll.Kappa_x();
      X[3] = coll.Kappa_y();
   }else
   {  X[1] = EmptyData;
      X[2] = EmptyData;
      X[3] = EmptyData;
   }
   X[4] = beam.s_s(m_);
   X[5] = beam.Sigma_s(m_);
   X[6] = EmptyData;
   BeamTxy.Point((beam.Emit[2]^0.5)(),beam.Emit[0](u_6*m_));

	if (EvolutionTimer.Timer(period))
   {
      Warning("Tcal[m.s]=", (double)EvolutionTimer.minutes +
               (double)(EvolutionTimer.seconds%60)/100.,
               " Tref[s]=", itime.t(s_)," Ecool[%]=", 100.*ExpKick/beam.Number(), false);
      Warning(" Core[%]=", 100.*iIBS.CoreNumber/beam.Number()," Impact=", iBeam.NImpact);

      if (iTarget.Multi && iTarget.itype == 2 && itime.t() > 0)
         X[6] = DrawPellet(X[1]);
      Evolut.Points(X, 7);

      if (iLaserCooling.Multi)
      {  LaserLight.Reset(0);
         for (int i = 0; i < iLaserCooling.hyst.Number; i++)
            LaserLight.Point(double(i)/iLaserCooling.hyst.Number, iLaserCooling.hyst[i]);
         LaserLight.Save();
      }

      Emittance.Save();
      Evolut.Save();
      BeamTxy.Save();
      return true;
   }else
      Evolut.Points(X, 7);

	return false;
}

bool xDraw::BeamRates(vectorU rates, xTime& itime, int period)
{
   if (Rates1.Enabled)
   {
      double X[5];
      X[0] = itime.t(s_);
      for (int i = 0; i < 4; i++)
         X[i+1] = rates[i](s_^-1);
      Rates1.Points(X, 5);
      if (RatesTimer.Timer(period))
      {  Rates1.Save();
         return true;
      }
   }
	return false;
}

void xDraw::LatticeFunctions(xRing& ring, xTime& time)
{
   if (ring.Autoskip)
      Lattice.AutoSkip = true;
   else
      Lattice.Size(ring.Number() + 2, 11);
   double dist = 0;
   double X[11];
   for (int i = 0; i < ring.Number(); i++)
   {
      xLattice& lat = ring.GetLattice(i);
      X[0] = dist;
      X[1] = lat.betax();
      X[2] = lat.betay();
      X[3] = lat.alfax();
      X[4] = lat.alfay();
      X[5] = lat.Dx();
      X[6] = lat.Dpx();
      X[7] = lat.Dy();
      X[8] = lat.Dpy();
      X[9] = lat.mux();
      X[10]= lat.muy();
      Lattice.Points(X, 11);
      dist += (ring[i].Length)(m_);
   }
   Lattice.Skip = 0;

   xLattice& lat = ring.GetLattice(0);
   X[0] = dist;
   X[1] = lat.betax();
   X[2] = lat.betay();
   X[3] = lat.alfax();
   X[4] = lat.alfay();
   X[5] = lat.Dx();
   X[6] = lat.Dpx();
   X[7] = lat.Dy();
   X[8] = lat.Dpy();
   X[9] = lat.mux();
   X[10]= lat.muy();
   Lattice.Points(X, 11);
   Lattice.Save();
}

//****************************************************************

void xDraw::Rates(xTaskRates& rate, xTime& time, xBeam& beam, xRing& ring)
{
   RateEh.SetXaxis(rate.dP_div,rate.dP_min(),rate.dP_max(),rate.dP_log);
   RateEh.SetYaxis(rate.E_div, rate.Eh_min(),rate.Eh_max(),rate.E_log);

   RateEv.SetXaxis(rate.dP_div,rate.dP_min(),rate.dP_max(),rate.dP_log);
   RateEv.SetYaxis(rate.E_div, rate.Ev_min(),rate.Ev_max(),rate.E_log);

   RatedP.SetXaxis(rate.dP_div,rate.dP_min(),rate.dP_max(),rate.dP_log);
   RatedP.SetYaxis(rate.E_div, rate.Eh_min(),rate.Eh_max(),rate.E_log);

   vectorU rates(s_^-1), emit;
   for (int i = 1; i < rate.dP_div + 2; i++)
   for (int j = 1; j < rate.E_div  + 2; j++)
   { 	if (RatesTimer.Timer(60) && RatesTimer.seconds)
  		   Warning("Calc.time [min] =", RatesTimer.minutes, ", Left :",
         floor(100.*(double(j)/(rate.E_div+2)+i) /(rate.dP_div+2)+0.5),"%");

      beam.Emit[0] = RateEh.GetPoint(0, j) * 1e-6;
      beam.Emit[1] = RateEv.GetPoint(0, j) * 1e-6;
      beam.Emit[2] = RateEh.GetPoint(i, 0);
      beam.Emit[2] = beam.Emit[2]^2;

      rates = xEffect::Summary(time, beam, ring);
      RateEh.SetPoint (i, j, rates[0]());
      RateEv.SetPoint (i, j, rates[1]());
      RatedP.SetPoint (i, j, rates[2]());

   }
   RateEv.Save();
   RateEh.Save();
   RatedP.Save();
}
//***************************************************************
void xDraw::FF(xForce& Fr, xTime& time, xEbeam& Ebeam, xRing& ring)
{
  //---- 05.07
  doubleU V (m_/s_);
  doubleU E_e (M_6*eV_ * 0.5110034);
  doubleU A(U1_);
  //int div = 1;
  vectorU ion_new(m_, U1_);

   ion_new = iEbeam.Ibeam2Ebeam(time, ion_new);
if(iEbeam.inside)
{
   //electron transverse spread
   //iEbeam.F.V_tr_e = ((iEbeam.F.Ttemp/E_e)^0.5)*U_c;
   // electron longitudinal spread
   //iEbeam.F.V_long_e = ((iEbeam.F.Ltemp/E_e)^0.5)*U_c;
   //effective electron spread
   //iEbeam.F.V_eff_e = ((iEbeam.F.TempEff/E_e)^0.5)*U_c;


   if(Fr.Surf)
   {
      FFtr.SetXaxis(Fr.div(),Fr.Vtr_min(),Fr.Vtr_max(),0);
      FFtr.SetYaxis(Fr.div(),Fr.Vlong_min(),Fr.Vlong_max(),0);

      FFlong.SetXaxis(Fr.div(),Fr.Vtr_min(),Fr.Vtr_max(),0);
      FFlong.SetYaxis(Fr.div(),Fr.Vlong_min(),Fr.Vlong_max(),0);

     for (int i = 1; i < Fr.div() + 2; i++)    // over Vtr
      for (int j = 1; j < Fr.div() + 2; j++)    // over Vlong
      {
//     iForce.Vtr = FFtr.GetPoint(0, j);
//     iForce.Vlong = FFlong.GetPoint(i, 0);
        iForce.Vtr = FFtr.GetPoint(i, 0);
        iForce.v[2] = FFlong.GetPoint(0, j);
        iForce.v[1] = 0.0;
        iForce.v[0] = iForce.Vtr;

     switch(iForce.type)
      {
       case 0: iForce.Budker(iEbeam.F); break;
       case 1: iForce.NonMag(iEbeam.F); break;
       case 2: iForce.DerSkr(iEbeam.F); break;
       case 3: iForce.Parhom(iEbeam.F); break;
       case 4: iForce.Toepffer(iEbeam.F); break;
       case 5: iForce.Table (iEbeam.F); break;
       case 6: iForce.D3 (iEbeam.F); break;
      }

//     FFtr.SetPoint(j, i, iForce.Ftr());
//     FFlong.SetPoint(j, i, iForce.f[2]());
        if(iForce.type == 6)iForce.Ftr = iForce.f[0];
        FFtr.SetPoint(i, j, iForce.Ftr());
        FFlong.SetPoint(i, j, iForce.f[2]());
      }
      FFtr.Save();
      FFlong.Save();
   }
   if(Fr.Line)
   {
      FCtr.Size(Fr.div2() + 2);
      FClong.Size(Fr.div2() + 2);

      for (int j = 0; j < Fr.div2()+1; j++)    // over Vlong
      {
        V = Fr.V_min + ((Fr.V_max-Fr.V_min)/Fr.div2*j);
        iForce.Vtr = V*U_Sin(Fr.Angle);
        iForce.v[2] = V*U_Cos(Fr.Angle);
        iForce.v[0] = iForce.Vtr;
        iForce.v[1] = 0.0;

        switch(iForce.type)
        {
           case 0: iForce.Budker(iEbeam.F); break;
           case 1: iForce.NonMag(iEbeam.F); break;
           case 2: iForce.DerSkr(iEbeam.F); break;
           case 3: iForce.Parhom(iEbeam.F); break;
           case 4: iForce.Toepffer(iEbeam.F); break;
           case 5: iForce.Table (iEbeam.F); break;
           case 6: iForce.D3 (iEbeam.F);  break;
        }
        if(iForce.type == 6)iForce.Ftr = iForce.f[0];
        FCtr.Point(V(), iForce.Ftr());
        FClong.Point(V(), iForce.f[2]());
      }

      FCtr.Save();
      FClong.Save();
   }
   if(Fr.Circle)
   {
      FAtr.Size(Fr.div3() + 2);
      FAlong.Size(Fr.div3() + 2);

      for (int j = 0; j < Fr.div3()+1; j++)    // over Vlong
      {
        A = (U_pi/(2.0*Fr.div3))*j;
        iForce.Vtr = Fr.Velocity*U_Sin(A);
        iForce.v[2] = Fr.Velocity*U_Cos(A);
        iForce.v[0] = iForce.Vtr;
        iForce.v[1] = 0.0;

        switch(iForce.type)
        {
           case 0: iForce.Budker(iEbeam.F); break;
           case 1: iForce.NonMag(iEbeam.F); break;
           case 2: iForce.DerSkr(iEbeam.F); break;
           case 3: iForce.Parhom(iEbeam.F); break;
           case 4: iForce.Toepffer(iEbeam.F); break;
           case 5: iForce.Table (iEbeam.F); break;
           case 6: iForce.D3 (iEbeam.F); break;
        }
        if(iForce.type == 6)iForce.Ftr = iForce.f[0];
        FAtr.Point(A(), iForce.Ftr());
        FAlong.Point(A(), iForce.f[2]());
      }

      FAtr.Save();
      FAlong.Save();
   }
}else
Warning("Electron beam is shifted outside the axis", PressEnter);
}

void xDraw::SpaceChargeDraw(xBeam& beam, xEbeam& ebeam)
{
   SpaceCharge.Reset(0);
   DriftVelocity.Reset(0);

   doubleU r(m_);
   doubleU Vdr(U1_);
   int steps = 200;
   doubleU dP;

   for (int i1 = 0; i1 <= steps*2; i1++)
   {
      if (ebeam.model == 0)
      {  r = ebeam.bradius/steps*(i1-steps);
         SpaceCharge.Point(r(m_), ebeam.UC_dp_P(r, ebeam.bcurrent, ebeam.bradius)(U1_));
         DriftVelocity.Point(r(m_), ebeam.UC_Vdrift(r, ebeam.bcurrent, ebeam.bradius)(U1_));
      }
      if (ebeam.model == 4)
      {  r = ebeam.r_circle/steps*(i1-steps);

         dP = ebeam.UC_dp_P(r, ebeam.I_circle, ebeam.r_circle)
            - ebeam.UC_dp_P(r, ebeam.I_hocle,  ebeam.r_hole)
            + ebeam.UC_dp_P(r, ebeam.I_hole,   ebeam.r_hole);
         SpaceCharge.Point(r(m_), dP());

         Vdr = ebeam.UC_Vdrift(r, ebeam.I_circle, ebeam.r_circle)
             - ebeam.UC_Vdrift(r, ebeam.I_hocle,  ebeam.r_hole)
             + ebeam.UC_Vdrift(r, ebeam.I_hole,   ebeam.r_hole);
         DriftVelocity.Point(r(m_), Vdr());
      }

   }

   SpaceCharge.Point();
   double dp, dp_max = (beam.Emit[2]^0.5)()*30;
   for (int i2 = 0; i2 < steps*2+1; i2++)
   {
      dp = dp_max/steps*(i2-steps);
      SpaceCharge.Point((iEcool.Lattice.Dx*dp-((ebeam.Initial[0]+ebeam.Final[0])/2))(m_),
                          dp-ebeam.Initial[5](U1_));
   }

   SpaceCharge.Point();
   SpaceCharge.Point(-1, (-ebeam.Initial[5])(U1_));
   SpaceCharge.Point( 1, (-ebeam.Initial[5])(U1_));
   SpaceCharge.Save();
   DriftVelocity.Save();
}

void xDraw::DrawLaserForce(xTime&time, xRing&ring, xLaserCooling&lcool)
{
      vectorU force, ion(m_, U1_);

      for (int i = -lcool.SweepStep2; i <= lcool.SweepStep2 ; i++)
      {
         ion[5] = (iBeam.Emit[2]^0.5)*10 * i / lcool.SweepStep2;

         force = lcool.GetForces(time, ring.Energy, ion);
         
         //LaserForce.Point((ion[5] * ring.Energy.Velocity)(m_/s_),
         LaserForce.Point(ion[5](U1_),
         (ring.Energy.Momentum*ring.Energy.Velocity*ring.Energy.Gamma*force[5])(eV_/m_));
      }
      LaserForce.Save();
}

void xDraw::DrawEBunch(xEbeam&ebeam)
{
   Space.Size(ebeam.Bunch.Number+1,6);
   double X[6];
   for (int i = 0; i < ebeam.Bunch.Number; i++)
   {
      if (!ebeam.Slice ||
          ((ebeam.From < ebeam.Bunch[i][4]) && (ebeam.Upto > ebeam.Bunch[i][4])))
      {
         for (int j = 0; j < 6; j++)
            X[j] = ebeam.Bunch[i][j];
         Space.Points(X, 6);
      }
   }
   Space.Save();
/*

   Tray<double>Hyst(6,Divisions+1);
   double h = 2 * Sigmas / Divisions;
   int coor, index;

   for (int n = 0; n < 6; n++)
   if(Distribution[n].Enabled)
   {
      Distribution[n].Size(Divisions+5);
      for (int k = 0; k <= Divisions; k++)
          Hyst[n][k] = 0;

      switch(n)
      {
       case  0: coor = 0; break;
       case  1: coor = 2; break;
       case  2: coor = 4; break;
       case  3: coor = 1; break;
       case  4: coor = 3; break;
       case  5: coor = 5; break;
      }

      for (int i = 0; i < ebeam.Bunch.Number; i++)
      {
         if (!ebeam.Slice ||
             ((ebeam.From < ebeam.Bunch[i][4]) && (ebeam.Upto > ebeam.Bunch[i][4])))
         {
            index=Round( (ebeam.Bunch[i][coor] / ebeam.Asigma[coor]()+Sigmas) / h );
            if (index < 0) index = 0;
            if (index > Divisions) index = Divisions;
            Hyst[n][index] += 1;
         }   
      }

      for (int j = 0; j <= Divisions; j++)
         Distribution[n].Point(-Sigmas + h*j, 100.*Hyst[n][j]/ebeam.Bunch.Number);

      Distribution[n].Save();
   }
*/
}


void xDraw::DrawEDensity(xEbeam&ebeam)
{
   if(ebeam.model != 5)
   {
      Warning("This procedure works for electron array only", PressEnter);
      return;
   }
   if(iForce.type != 6)
   {
      Warning("This procedure works with 3D force only", PressEnter);
      return;
   }

   vectorU ion_new(m_, U1_);
   ebeam.EnableShifts = false;

   ion_new = ebeam.Ibeam2Ebeam(iTime, ion_new);

   EDensity.SetXaxis(ebeam.N_hor,-3.0,3.0,0);
   EDensity.SetYaxis(ebeam.N_long,-3.0,3.0,0);

   for (int i = 1; i < ebeam.N_hor + 2; i++)    // over horizontal coord
      for (int j = 1; j < ebeam.N_long + 2; j++)    // over longitudinal coord
      {
        ion_new[0] = EDensity.GetPoint(i, 0)*ebeam.Asigma[0];
        ion_new[4] = EDensity.GetPoint(0, j)*ebeam.Asigma[4];

        ion_new = ebeam.Ibeam2Ebeam(iTime, ion_new);
        EDensity.SetPoint(i, j, ebeam.F.n_e());

      }
   EDensity.Save();

}

void xDraw::DrawStochForce(xTime&time, xRing&ring, xBeam&beam, xStochastic&scool)
{
      vectorU force, ion(m_, U1_);
      doubleU diffusion(m_^-1);
      scool.eta = ring.Eta;

      for (int i = 0; i <= 100 ; i++)
      {
         ion[5] = (beam.Emit[2]^0.5) * (i - 50) * 4.0 / 50.0;
         force = scool.GetForces(time, ring.Energy, ion);
         diffusion = scool.GetDiffusion(time, ring.Energy, beam, ion);
         StochForce.Point(ion[5](U1_), force[5](m_^-1));
         StochDiff.Point(ion[5](U1_), diffusion(m_^-1));
      }
      StochForce.Save();
      StochDiff.Save();
}

double xDraw::DrawPellet(double L)
{
   if (L == EmptyData) return L;

   double a[2];
   a[0] =((iBeam.Emit[0]*iTarget.Lattice.betax)^0.5)();
   a[1] =((iBeam.Emit[1]*iTarget.Lattice.betay)^0.5)();

   Tray<double>Hyst(6, Divisions+1);
   iBeam.Amplitudes(iTarget.Lattice);
   for (int n5 = 3; n5 < 5; n5++)
      iBeam.Profile(Hyst[n5], n5, a[n5-3], Divisions);

   xPellet iPellet;
   iPellet.OnGet();
   double X[3];

   BeamSize.Size(5,3);
   X[0] = iPellet.pellets[iPellet.Nplts-1].position[4]; X[1] = -a[0]*1e3;  X[2] = -a[1]*1e3;
   BeamSize.Points(X, 3);
   X[0] = 0; X[1] = -a[0]*1e3;  X[2] = -a[1]*1e3;
   BeamSize.Points(X, 3);
   X[0] = 0; X[1] = a[0]*1e3;  X[2] = a[1]*1e3;
   BeamSize.Points(X, 3);
   X[0] = iPellet.pellets[iPellet.Nplts-1].position[4]; X[1] = a[0]*1e3;  X[2] = a[1]*1e3;
   BeamSize.Points(X, 3);
   BeamSize.Save();

   BeamCros.Size(362, 2);
   for (int c = 0; c < 361; c++)
      BeamCros.Point(a[0]*1e3*cos(c*M_PI/180), a[1]*1e3*sin(c*M_PI/180));
   BeamCros.Save();     

   double Z[6];
   Pellets.Size(iPellet.Nplts + 3, 6);
   for (int p = 0; p < iPellet.Nplts; p++)
   {  for (int z = 0; z < 6; z++)
         Z[z] = iPellet.pellets[p].position[z];
      Pellets.Points(Z, 6);
   }
   Pellets.Save();

   iPellet.Average = 0;
   iPellet.tt = fabs(iPellet.pellets.L * iPellet.pellets.Nplts / iPellet.X0[5]);
   iPellet.N = (int)(iPellet.tt / iPellet.t.dt);

   Pellet.Size(iPellet.N + 1, 2);

   for (int i = 0; i < iPellet.N; i++)
   {
      iPellet.PxPy = 0;

      for (int j = 0; j < iPellet.pellets.Nplts; j++)
      {
         iPellet.PxPy += iPellet.pellets.getPxPy(iPellet.pellets[j], Hyst[3], Hyst[4],
                   Divisions+1, a[0]*Sigmas*1e3, a[1]*Sigmas*1e3);
         iPellet.pellets[j].setPosition(iPellet.t);
      }

      double k = Divisions / Sigmas / iBeam.Number();
      iPellet.PxPy *= k*k;
      Pellet.Point(iPellet.t.t, iPellet.PxPy);
      iPellet.Average += iPellet.PxPy;
      iPellet.t+=1;
   }
   iPellet.Average /= iPellet.N;
//---------------------- Effective luminosity ----------------------------------
   double Ldetector = 0;
   double Kdetector = EmptyData;
   if (L)
   {  Kdetector = iPellet.Average * iPellet.pellets.Detector / L;
      for (int i = 0; i < iPellet.N; i++)
      {
         if (Pellet.Data[i][1] > Kdetector)
            {
                  if (cut == 1) // (cutTop)                  //Krestnikov 14.01.2009
                        Ldetector += Kdetector;
            }
         else
            Ldetector += Pellet.Data[i][1];
      }
      Ldetector /= iPellet.N * iPellet.Average;
   }
//------------------------------------------------------------------------------
   Average.Size(5,3);
   X[0] = 0; X[1] = iPellet.Average; X[2] = Kdetector;
   Average.Points(X, 3);
   X[0] = iPellet.tt; X[1] = iPellet.Average; X[2] = Kdetector;
   Average.Points(X, 3);
   Average.Save();
   Pellet.Save();
   return Ldetector * L;
}

void xDraw::Oscilograph()      
{
   BData data;
   data.Load(xData::Data.OnGetC(35, 18, "no"));
   xGraf osc;
   osc.SetName("osilograph.cur");
   osc.Size(data.Row());
   double step = data[1][0];
   double time = 0;
   for (int i = 4; i < data.Row(); i++)
   {  osc.Point(time, -data[i][0]);
      time += step;
   }
   osc.Save();   
}

void xDraw::DrawGain()
{
   Gain.Size(1005,4);
   Gain.SetName("gain.cur");
   vectorU rates;
   double X[4];

   for (double g = 0; g < 1.; g += 0.001)
   {
      for (int i = 0; i < iGated.Number; i++)
         iGated.System[i].gain = g;
      rates = iGated.Rates(iTime, iBeam, iRing);
      X[0] = g;
      X[1] = rates[0](s_^-1);
      X[2] = rates[1](s_^-1);
      X[3] = rates[2](s_^-1);
      Gain.Points(X, 4);
   }
   Gain.Save();

}

void xDraw::DrawLStoch()
{
   LStoch.Size(2002,4);
   LStoch.SetName("lstoch.cur");
   double X[3];

   for (double p = iGated.dpF() * (-2); p < iGated.dpF() * 2; p += iGated.dpF()/500)
   {
      X[0] = p;
      X[1] = iGated.F() * iGated.getForce   ( p / iGated.dpF() );
      X[2] = sqrt(iGated.D() * iGated.getDifusion( p / iGated.dpD(), ((iBeam.Emit[2]^0.5)/iGated.sigma)()));
      LStoch.Points(X, 3);
   }
   LStoch.Save();
    
}

void xDraw::TuneShift() {
      double CellSize;
      if (iBeam.benum == COASTING)
      {  
         CellSize = (iRing.Circ / iBeam.CellNum)(m_);
         for (int j4 = 0; j4 < iBeam.Number(); j4++)
         {  iBeam[j4][4] -= CellSize * floor(iBeam[j4][4] / CellSize);
            iBeam[j4][4] -= CellSize / 2;
         }
      }
      else CellSize = iBeam.Emit[4](m_);
      iIBS.iTune.SetGrid ( iIBS.CellNum[0]+1,iIBS.CellNum[1]+1,iIBS.CellNum[2]+1,
                           sqrt(iBeam.Emit[0](m_)*iRing.LATTICE.betax(m_)) / (iIBS.CellNum[0] / iIBS.GridSize[0]), 
                           sqrt(iBeam.Emit[1](m_)*iRing.LATTICE.betay(m_)) / (iIBS.CellNum[1] / iIBS.GridSize[1]), 
                           CellSize / (iIBS.CellNum[2] / iIBS.GridSize[2]));
      iIBS.iTune.SetValue();
      for (int n = 0; n < iBeam.Number(); n++)
         iIBS.iTune.AddCharge(iBeam[n], iRing.Energy.Gamma(U1_));

      bpData Potentials(iIBS.iTune.xyz[0]+1, iIBS.iTune.xyz[2]+1);
      bpData ElectricField(iIBS.iTune.xyz[0]+1, iIBS.iTune.xyz[2]+1);

      for (int i = 0; i < iIBS.iTune.xyz[0]; i++) {
         Potentials[i+1][0] = (i-(iIBS.iTune.xyz[0]-1)/2) * iIBS.iTune.abc[0]; 
         ElectricField[i+1][0] = Potentials[i+1][0]; }
      for (int k = 0; k < iIBS.iTune.xyz[2]; k++) {
         Potentials[0][k+1] = (k-(iIBS.iTune.xyz[2]-1)/2) * iIBS.iTune.abc[2]; 
         ElectricField[0][k+1] = Potentials[0][k+1]; }

      int j = iIBS.VerticalCell + (iIBS.iTune.xyz[1]-1)/2;
      if (j < 0) j = 0;
      if (j >= iIBS.iTune.xyz[1]) j = iIBS.iTune.xyz[1]-1;
      for (int i = 0; i < iIBS.iTune.xyz[0]-1; i++)
      for (int k = 0; k < iIBS.iTune.xyz[2]-1; k++) {
         Potentials[i+1][k+1] = iIBS.iTune.Potent[i][j][k];
         //ElectricField[i+1][k+1] = iIBS.iTune.FieldZ[i][j][k]; 
         double field[3];
         iIBS.iTune.Field(field,i,j,k,CellSize,bool(iBeam.benum)); 
         ElectricField[i+1][k+1] = field[iIBS.FieldComponent]; 
      }

      Potentials.SaveFileSep("potentials.sur", "\t");
      ElectricField.SaveFileSep("electricfield.sur", "\t");

}
