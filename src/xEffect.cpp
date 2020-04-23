//---------------------------------------------------------------------------
#include "stdafx.h"
#define _CRT_SECURE_NO_WARNINGS
#include "xEffect.h"
#include "xDraw.h"

//---------------------------------------------------------------------------
// Effects
//---------------------------------------------------------------------------

xEffect::xEffect()
{
   Multi = 0;
   Loss = false;
   ID = 0;
}

vectorU xEffect::Summary(xTime &t, xBeam &beam, xRing &ring)
{
   vectorU rates(s_ ^ -1);
   for (int i = 0; i < xEffect::ACount; i++)
      if (xEffect::AItems[i]->Multi != 0)
         rates += xEffect::AItems[i]->Rates(t, beam, ring);
   return rates;
}

xEffect *xEffect::GetEffect(int id)
{
   for (int i = 0; i < xEffect::ACount; i++)
      if (xEffect::AItems[i]->ID == id)
         return xEffect::AItems[i];
   return NULL;
}

xEffect *xEffect::GetEffect(const char *name)
{
   char str[20];
   strcpy(str, name);

   for (int i = 0; i < xEffect::ACount; i++)
      if (xEffect::AItems[i]->EffectName == str)
         return xEffect::AItems[i];

   Warning("No Effect found. BETACOOL code has a following Effects:");
   for (int j = 0; j < xEffect::ACount; j++)
      Warning(xEffect::AItems[j]->EffectName.get(), EmptyData,
              ": betax [m] =", xEffect::AItems[j]->Lattice.betax(),
              ", betay [m] =", xEffect::AItems[j]->Lattice.betay());

   return NULL;
}

void xEffect::SetLattice(xLattice &L)
{
   for (int i = 0; i < xEffect::ACount; i++)
      //if (!xEffect::AItems[xEffect::ACount-1]->ID)
      xEffect::AItems[i]->Lattice = L;
}

//---------------------------------------------------------------------------
// Task Rates
//---------------------------------------------------------------------------

xTaskRates::xTaskRates()
{
   Eh_min._(m_);
   Eh_max._(m_);
   Ev_min._(m_);
   Ev_max._(m_);
   dP_min._(U1_);
   dP_max._(U1_);
}

int xTaskRates::OnGet()
{
   xEffect *xeffect;
   for (int i = 1; i < Data.Col(31); i++)
   {
      xeffect = xEffect::GetEffect(i);
      if (xeffect)
         xeffect->Multi = Data.OnGetI(31, i, 0);
   }
   show3D = Data.OnGetB(32, 9, false);

   Eh_min = Data.OnGet(33, 1, 1);
   Eh_max = Data.OnGet(33, 2, 10);
   Ev_min = Data.OnGet(33, 3, 1);
   Ev_max = Data.OnGet(33, 4, 10);
   E_div = Data.OnGetI(33, 5, 10);
   E_log = Data.OnGetB(33, 6, true);
   dP_min = Data.OnGet(33, 7, 1);
   dP_max = Data.OnGet(33, 8, 10);
   dP_div = Data.OnGetI(33, 9, 10);
   dP_log = Data.OnGetB(33, 10, true);

   return 0;
}

int xTaskRates::OnRun()
{
   xEffect *xeffect;
   for (int i = 1; i < Data.Col(31); i++)
   {
      xeffect = xEffect::GetEffect(i);
      if (xeffect)
         xeffect->Multi = Data.OnGetI(31, i, 0);
   }
   return 0;
}

//---------------------------------------------------------------------------
// Particle losses
//---------------------------------------------------------------------------

xLosses::xLosses()
{
   ID = 5;
   Decay = false;
   Acceptance = false;
   Separatrix = false;
   Injection = false;
   EffectName = "XLOSSES";
}

int xLosses::OnGet()
{
   xEffect *xeffect;
   for (int i = 1; i < 5; i++)
   {
      xeffect = xEffect::GetEffect(i);
      if (Multi != 0)
         xeffect->Loss = Data.OnGetB(9, i, false);
      else
         xeffect->Loss = false;
   }

   Decay = Data.OnGetB(9, 5, false);
   Acceptance = Data.OnGetB(9, 6, false);
   Separatrix = Data.OnGetB(9, 7, false);
   Injection = Data.OnGetB(9, 8, false);

   return 0;
}

int xLosses::OnRun()
{
   return OnGet();
}

vectorU xLosses::Rates(xTime &time, xBeam &beam, xRing &ring)
{
   vectorU rates(s_ ^ -1);
   if (Decay)
      rates[3] -= 1 / ring.TLife;
   return rates;
}

void xLosses::Kick(xTime &time, xBeam &beam, xRing &ring)
{
   double loss;
   beam.RMS(beam.Emit, Lattice);

   if (Decay)
   {
      loss = (time.dt / ring.TLife)(U1_);
      for (int j = 0; j < beam.Number(); j++)
         beam.Loss(Lattice, j, loss, false);
   }
   if (Acceptance)
   {
      beam.Courant_Snyder(Lattice);
      for (int j = beam.Number() - 1; j >= 0; j--)
      {
         //if ((U_Abs(beam(j,0)) > ((ring.AcceptH*ring.BetaH)^0.5)) || (U_Abs(beam(j,2)) > ((ring.AcceptV*ring.BetaV)^0.5)) )
         if (((beam.Inv(j)[0] * ring.BetaH) ^ 0.5) + ((beam.Inv(j)[2] ^ 0.5) * ring.DispH) > ((ring.AcceptH * ring.BetaH) ^ 0.5))
            beam.Loss(Lattice, j, 1., false);
         else if (((beam.Inv(j)[1] * ring.BetaV) ^ 0.5) + ((beam.Inv(j)[2] ^ 0.5) * ring.DispV) > ((ring.AcceptV * ring.BetaV) ^ 0.5))
            beam.Loss(Lattice, j, 1., false);
         else
             //if (U_Abs(beam(j,5)) > ring.AcceptDp)
             if (beam.Inv(j)[2] > (ring.AcceptDp * ring.AcceptDp))
            beam.Loss(Lattice, j, 1., false);
      }
   }
   if (Separatrix)
   {
      for (int j = 0; j < beam.Number(); j++)
      {
         if (beam.benum == BUCKET)
            if (U_Abs(beam(j, 4)) > ring.AcceptS / 2.)
               beam.Loss(Lattice, j, 1., false);
         if (beam.benum == BUNCHED)
            if (U_Abs(beam(j, 4)) > ring.L_s * ring.Size_s / 200.)
               beam.Loss(Lattice, j, 1., false);
      }
   }
}

//---------------------------------------------------------------------------
// Additional heating
//---------------------------------------------------------------------------

xHeat::xHeat()
{
   ID = 7;
   EffectName = "XADDHEAT";

   rate = false;
   rnum = false;
   Rate[0]._(s_ ^ -1);
   Rate[1]._(s_ ^ -1);
   Rate[2]._(s_ ^ -1);

   linear = false;
   lnum = false;
   Linear[0]._(u_6 * m_ / s_);
   Linear[1]._(u_6 * m_ / s_);
   Linear[2]._(u_6 * (s_ ^ -1));

   power = false;
   pnum = false;
   Power[0]._(s_ ^ -1);
   Power[1]._(s_ ^ -1);
   Power[2]._(s_ ^ -1);

   diffusion = false;
   dnum = false;
   Diffusion[0]._(p_12 * (m_ ^ 2) / s_);
   Diffusion[1]._(p_12 * (m_ ^ 2) / s_);
   Diffusion[2]._(p_12 * (s_ ^ -1));
   //------------------------------23.11.05--------------------------------------
   DiffusP._(s_ ^ -1);
   //----------------------------
   tau._(s_);
   position._(m_);
   length._(m_);
}

int xHeat::OnGet()
{
   rate = Data.OnGetB(41, 1, false);
   linear = Data.OnGetB(41, 5, false);
   power = Data.OnGetB(41, 9, false);
   diffusion = Data.OnGetB(41, 13, false);

   for (int i = 0; i < 3; i++)
   {
      Rate[i] = Data.OnGet(41, i + 2, 1.);
      Linear[i] = Data.OnGet(41, i + 6, 1.);
      Power[i] = fabs(Data.OnGet(41, i + 10, 1.));
      Diffusion[i] = fabs(Data.OnGet(41, i + 14, 1.));
   }
   difP = Data.OnGetB(41, 17, false);
   DiffusP = 1e-6 * Data.OnGet(41, 18, 1.0);
   tapered = Data.OnGet(41, 19, 0);
   rnum = Data.OnGetB(41, 20, false);
   lnum = Data.OnGetB(41, 21, false);
   pnum = Data.OnGetB(41, 22, false);
   dnum = Data.OnGetB(41, 23, false);

   coherent = Data.OnGetB(41, 24, false);
   k1 = Data.OnGet(41, 25, 1);
   k2 = Data.OnGet(41, 26, 1);
   tau = Data.OnGet(41, 27, 1);
   position = Data.OnGet(41, 28, 1);
   length = Data.OnGet(41, 29, 1);
   return 0;
}

int xHeat::OnRun()
{
   return OnGet();
}

vectorU xHeat::Rates(xTime &t, xBeam &beam, xRing &ring)
{
   vectorU rates(s_ ^ -1);
   double num, No = (beam.Emit[3] / beam.InitialEmit[3])();

   for (int i = 0; i < 3; i++)
   {
      if (rate)
      {
         if (rnum)
            num = No;
         else
            num = 1;
         rates[i] += Rate[i] * num;
      }
      if (linear)
      {
         if (lnum)
            num = No;
         else
            num = 1;
         rates[i] += Linear[i] / beam.Emit[i] * num;
      }
      if (power)
      {
         if (pnum)
            num = No;
         else
            num = 1;
         rates[i] += Power[i] * 2 * num;
      }
      if (diffusion)
      {
         if (dnum)
            num = No;
         else
            num = 1;

         if (i == 2 && difP)
            rates[i] += (DiffusP / beam.Emit[i]) * num;
         else
            rates[i] += Diffusion[i] * num / (beam.Emit[i] ^ 2) * 2;
      }
   }
   return rates;
}

void xHeat::Kick(xTime &time, xBeam &beam, xRing &ring)
{
   beam.Emit = beam.Emit_Def(Lattice);
   double r[3], k = 1;
   if (beam.benum == BUNCHED)
      k = 2;
   double num, No = (beam.Emit[3] / beam.InitialEmit[3])();

   if (rate)
   {
      if (rnum)
         num = No;
      else
         num = 1;
      for (int i = 0; i < 3; i++)
      {
         r[i] = (Rate[i] * time.dt)(U1_)*num;
         if (i == 2)
            r[i] *= k / 2.;
         if (r[i] > -1.)
            r[i] += 1.;
         else
            r[i] = exp(r[i]);
      }

      for (int j = 0; j < beam.Number(); j++)
      {
         for (int i = 0; i < 2; i++)
            beam(j, i * 2 + 1, beam(j, i * 2 + 1) * r[i]);
         beam(j, 5, beam(j, 5) * r[2] + (tapered * beam(j, 0) * ring.Energy.Gamma2 / ring.DispH));
      }
   }
   if (linear)
   {
      if (lnum)
         num = No;
      else
         num = 1;
      for (int i = 0; i < 3; i++)
      {
         r[i] = (Linear[i] * time.dt / beam.Emit[i])(U1_)*num;
         if (i == 2)
            r[i] *= k / 2.;
         if (r[i] > -1.)
            r[i] += 1.;
         else
            r[i] = exp(r[i]);
      }

      for (int j = 0; j < beam.Number(); j++)
         for (int i = 0; i < 3; i++)
            beam(j, i * 2 + 1, beam(j, i * 2 + 1) * r[i]);
   }
   if (power)
   {
      if (pnum)
         num = No;
      else
         num = 1;
      r[0] = ((2 * Power[0] * time.dt * 2 * beam.Emit[0] / Lattice.betax) ^ 0.5)(U1_);
      r[1] = ((2 * Power[1] * time.dt * 2 * beam.Emit[1] / Lattice.betay) ^ 0.5)(U1_);
      r[2] = ((k * Power[2] * time.dt * 2 * beam.Emit[2]) ^ 0.5)(U1_);

      for (int j = 0; j < beam.Number(); j++)
         for (int i = 0; i < 3; i++)
            beam.Add(j, i * 2 + 1, r[i] * xDistributor::Gaussian() * sqrt(num));
   }
   if (diffusion)
   {
      if (dnum)
         num = No;
      else
         num = 1;
      r[0] = ((2 * Diffusion[0] * time.dt * 2 / beam.Emit[0] / Lattice.betax) ^ 0.5)(U1_);
      r[1] = ((2 * Diffusion[1] * time.dt * 2 / beam.Emit[1] / Lattice.betay) ^ 0.5)(U1_);
      if (difP)
         r[2] = ((k * DiffusP * time.dt) ^ 0.5)(U1_);
      else
         r[2] = ((k * Diffusion[2] * time.dt * 2 / beam.Emit[2]) ^ 0.5)(U1_);

      for (int j = 0; j < beam.Number(); j++)
         for (int i = 0; i < 3; i++)
            beam.Add(j, i * 2 + 1, r[i] * xDistributor::Gaussian() * sqrt(num));
   }
   if (coherent)
   {
      for (int j = 0; j < beam.Number(); j++)
         if ((beam(j, 4) >= position - (length / 2)) && beam(j, 4) <= position + (length / 2))
            beam.Add(j, 5, sin(k1 * beam(j, 5)(U1_)) * time.dt / k2 / tau);
   }
}

//---------------------------------------------------------------------------
//- Effect collisions
//---------------------------------------------------------------------------

xColl::xColl()
{
   ID = 4;
   EffectName = "XCOLLISION";
   Cross._(barn_);
   Luminosity._((m_ ^ -2) / s_);
   Events._(0, U1_);
   Realtime = false;
   CollBeam = false;
   // 11.12.2007
   Include = false;
   Diffusion._(U1_);
   hourglass = 1;
}

int xColl::OnGet()
{
   Points = Data.OnGetI(10, 1, 1);
   Cross = Data.OnGet(10, 2, 10);
   LuminosityModel = Data.OnGetI(10, 6, 2);
   CollBeam = Data.OnGetB(10, 7, false);
   Divisions = Data.OnGetI(10, 8, 100);
   From = Data.OnGet(10, 9, 0);
   Upto = Data.OnGet(10, 10, 100);
   Steps = Data.OnGet(10, 11, 10);
   BB_EmitDef = Data.OnGetI(10, 12, 0);
   percent = Data.OnGet(10, 13, 60);
   Local_BB = Data.OnGetB(10, 14, false);
   Include = Data.OnGetB(10, 15, false);
   Diffusion = Data.OnGet(10, 16, 0.0001);
   hourglass = Data.OnGet(10, 17, 1);

   return 0;
}

int xColl::OnRun()
{
   OnGet();
   OnSet();
   return 0;
}

int xColl::OnSet()
{
   Lattice.betax = Data.OnGet(10, 3, 1);
   Lattice.betay = Data.OnGet(10, 4, 1);
   if (iBeam.benum == BUNCHED)
      hourglass = HourGlass();
   Data.OnSet(10, 17, hourglass);
   return 0;
}

double xColl::HourGlass()
{
   double s_b = (Lattice.B_s * ((iBeam.Emit[2] / Lattice.betax / Lattice.betay) ^ 0.5))(U1_);
   double integral = 0;
   double du = 0.1;
   for (int i = -100; i <= 100; i++)
      integral += exp(-(du * i) * (du * i)) / (1. + (du * i * s_b) * (du * i * s_b));
   return integral * du / sqrt(M_PI);
}

vectorU xColl::Rates(xTime &t, xBeam &beam, xRing &ring)
{
   vectorU rates(s_ ^ -1);
   beam.RMS(beam.Emit, Lattice);
   if (CollBeam)
      Luminosity = beam.N_b * beam.Emit[3] * cBeam.Emit[3] /
                   ((((cBeam.Sigma_x * cBeam.Sigma_x + beam.Sigma_2[0]) *
                      (cBeam.Sigma_y * cBeam.Sigma_y + beam.Sigma_2[2])) ^
                     0.5) *
                    (ring.Trev) * U_pi * 2.0);
   else
      Luminosity = beam.N_b * (beam.Emit[3] ^ 2) /
                   (((beam.Sigma_2[0] * beam.Sigma_2[2]) ^ 0.5) * (ring.Trev) * U_pi * 4.0);

   if (beam.benum == BUNCHED && beam.HourGlass)
      Luminosity *= HourGlass();
   if (Loss)
      rates[3] = -Points * Cross * Luminosity / (beam.Emit[3] * beam.N_b);

   BeamBeam(t, beam, ring);
   //11.12.2007
   if (Include)
   {
      rates[0] += (0.25 * Diffusion * Diffusion * Kappa_x * Kappa_x * U_pi * U_pi) / (2.0 * ring.Trev);
      rates[1] += (0.25 * Diffusion * Diffusion * Kappa_y * Kappa_y * U_pi * U_pi) / (2.0 * ring.Trev);
      if (CollBeam)
      {
         rates[0] *= beam.Sigma_2[0] / (cBeam.Sigma_x * cBeam.Sigma_x);
         rates[1] *= beam.Sigma_2[2] / (cBeam.Sigma_y * cBeam.Sigma_y);
      }
   }

   return rates;
}
//---- 05.06 ----
void xColl::Kick(xTime &time, xBeam &beam, xRing &ring)
{
   if (Realtime)
      CalcLuminosity(time, beam, ring);

   BeamBeam(time, beam, ring);

   beam.RMS(beam.Emit, Lattice);
   doubleU lumi((m_ ^ -2) / s_);
   Luminosity = 0;

   if (!CollBeam)
   {
      cBeam.Energy = ring.Energy;
      cBeam.Lattice = Lattice;
      cBeam.Sigma_x = (beam.Emit[0] * Lattice.betax) ^ 0.5;
      cBeam.Sigma_y = (beam.Emit[1] * Lattice.betay) ^ 0.5;
      cBeam.Sigma_s = Lattice.B_s * (beam.Emit[2] ^ 0.5);
      cBeam.Emit[3] = beam.Emit[3];
   }
   //double part_sig = 0.05;

   //   double density2;
   double sum = 0;
   double hx, hy;
   Tine<double> density(Divisions);

   doubleU rc(m_ ^ -2);
   doubleU rt(m_ ^ -2);
   doubleU sigma_x_C(m_);
   doubleU sigma_y_C(m_);
   doubleU sigma_s_C(m_);
   doubleU sigma_x_T(m_);
   doubleU sigma_y_T(m_);
   doubleU sigma_s_T(m_);
   vectorU Ec, Et;
   double Nc[3];
   //   double Nc_mean;
   double CoreNumber;
   doubleU x(m_);
   doubleU y(m_);
   doubleU s(m_);

   hx = beam.Sigma[0]() * beam.Sigmas / (Divisions - 0.0);
   hy = beam.Sigma[2]() * beam.Sigmas / (Divisions - 0.0);

   if ((LuminosityModel == 4) && (!CollBeam))
   {
      ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~14/01/2008~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      beam.Emit = beam.Emit_RMS(Lattice);
      Et = beam.Emit;
      Ec = beam.Emit_FWHM(Lattice);
      beam.Courant_Snyder(Lattice);
      for (int i = 0; i < 3; i++)
      {
         beam.Powel2(beam.Fitting[i], beam.Hyst3[i], beam.Sigmas, beam.Division);
         if (fabs(beam.Fitting[i][1]) > fabs(beam.Fitting[i][2]))
         {
            Ec[i] = beam.Emit[i] * beam.Fitting[i][2] * beam.Fitting[i][2];
            Et[i] = beam.Emit[i] * beam.Fitting[i][1] * beam.Fitting[i][1];
            Nc[i] = beam.Number() / (1. + fabs(beam.Fitting[i][1] / beam.Fitting[i][2] *
                                               beam.Fitting[i][3] / beam.Fitting[i][4]));
         }
         else
         {
            Ec[i] = beam.Emit[i] * beam.Fitting[i][1] * beam.Fitting[i][1];
            Et[i] = beam.Emit[i] * beam.Fitting[i][2] * beam.Fitting[i][2];
            Nc[i] = beam.Number() / (1. + fabs(beam.Fitting[i][2] / beam.Fitting[i][1] *
                                               beam.Fitting[i][4] / beam.Fitting[i][3]));
         }
      }

      sigma_x_C = ((Ec[0] * Lattice.betax) ^ 0.5);
      sigma_y_C = ((Ec[1] * Lattice.betay) ^ 0.5);
      sigma_x_T = ((Et[0] * Lattice.betax) ^ 0.5);
      sigma_y_T = ((Et[1] * Lattice.betay) ^ 0.5);

      CoreNumber = 0;

      for (int i = 0; i < beam.Number(); i++)
      {
         bool InsideCore = true;
         for (int j = 0; j < 3; j++)
            if (Ec[j] * 2 < beam.Inv(i, j))
               InsideCore = false;

         if (InsideCore)
         {
            CoreNumber += 1;
         }
      }

      //Nc_mean = 0.5 * (Nc[1] + Nc[2]);
      if (beam.benum == BUNCHED)
      {
         sigma_s_C = Lattice.B_s * (Ec[2] ^ 0.5);
         sigma_s_T = Lattice.B_s * (Et[2] ^ 0.5);
      }
      ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~14/01/2008~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      beam.Emit = beam.Emit_Def(Lattice);
   }

   if ((LuminosityModel == 1) && (!CollBeam))
   {
      beam.RadialDensity(density(), Divisions, hx, hy);
      for (int k = 0; k < Divisions; k++)
         sum += density[k] * hx * hy * M_PI * (2. * k + 1.);
   }
   if ((LuminosityModel == 2) && (!CollBeam))
   {
      Tray<double> Profile(2, Divisions * 2 + 1);
      beam.Amplitudes(Lattice);
      beam.Profile(Profile[0], 3, beam.Sigma[0](), Divisions * 2);
      beam.Profile(Profile[1], 4, beam.Sigma[2](), Divisions * 2);
      beam.ProfileDensity(density(), Divisions, Profile(), hx, hy);

      for (int k = 0; k < Divisions; k++)
         sum += density[k] * hx * hy * M_PI * (2. * k + 1.);
   }
   /*
   double dnst[250];
   for (int d = 0; d < 250; d++)
      dnst[d] = density[d] / 1e9;

   xGraf Density;
   Density.SetName("density.cur");
   Density.Reset(0);
   for (int i = 0; i < 250; i++)
      Density.Point(beam.Sigmas * (1. / Divisions * i - 0.), dnst[i]);
   Density.Save();
*/
   for (int j = 0; j < beam.Number(); j++)
   {
      lumi = 0.0;
      if (!CollBeam)
      {
         cBeam.Energy = ring.Energy;
         cBeam.Lattice = Lattice;

         if (LuminosityModel == 0)
         {
            lumi = beam.Emit[3] / beam.Number() / ring.Trev *
                   beam.LocalDensity(beam[j], Divisions, Lattice);
         }
         else if (LuminosityModel == 1)
         {
            for (int k = 1; k <= Divisions; k++)
               if ((((beam[j][0] / (hx * k)) * (beam[j][0] / (hx * k))) + ((beam[j][2] / (hy * k)) * (beam[j][2] / (hy * k)))) <= 1.)
               {
                  //lumi = beam.Emit[3] / sum / ring.Trev *
                  //doubleU(density[k-1], m_^-2);
                  lumi = beam.Emit[3] / beam.Number() / ring.Trev *
                         doubleU(density[k - 1], m_ ^ -2);
                  break;
               }
         }
         else if (LuminosityModel == 2)
         {
            for (int k = 1; k <= Divisions; k++)
               if (pow((beam[j][0] / (hx * k)), 2.0) + pow((beam[j][2] / (hy * k)), 2.0) < 1.)
               {
                  lumi = beam.Emit[3] / sum / ring.Trev *
                         doubleU(density[k - 1], m_ ^ -2);
                  //lumi = beam.Emit[3] / beam.Number() / ring.Trev *
                  //doubleU(density[k-1], m_^-2);
                  break;
               }
         }
         else if (LuminosityModel == 3)
         {
            /*
            lumi = beam.Emit[3] / beam.Number()  / ring.Trev *
                   beam.LocalDensity_G(beam[j], Lattice);
            */
            //lumi = beam.Emit[3] / beam.Number()  / ring.Trev *
            //beam.LocalDensity_R(beam[j], part_sig, Lattice);
            lumi = MalakhovDensity(beam[j]) / ring.Trev;
         }
         else if (LuminosityModel == 4)
         {
            /*
             x = doubleU(beam[j][0],m_);
             y = doubleU(beam[j][2],m_);
             s = doubleU(beam[j][4],m_);
             rc = 0;
             rt = 0;

             if((x < 6.0*sigma_x_C)&&(y < 6.0*sigma_y_C))
                 {
                  rc = Nc_mean/(2.*U_pi*sigma_x_C*sigma_y_C);
                  rc *= 1./U_Exp(((x/sigma_x_C)^2.)/2.);
                  rc *= 1./U_Exp(((y/sigma_y_C)^2.)/2.);
                 }
               if((x < 6.0*sigma_x_T)&&(y < 6.0*sigma_y_T))
                 {
                  rt = (beam.Number() - Nc_mean)/(2.*U_pi*sigma_x_T*sigma_y_T);
                  rt *= 1./U_Exp(((x/sigma_x_T)^2.)/2.);
                  rt *= 1./U_Exp(((y/sigma_y_T)^2.)/2.);
                 }
             */
            /*
             if(!beam.bunched)
             {
               if((x < 6.0*sigma_x_C)&&(y < 6.0*sigma_y_C))
                 {
                  rc = Nc_mean/(2.*U_pi*sigma_x_C*sigma_y_C);
                  rc *= 1./U_Exp(((x/sigma_x_C)^2.)/2.);
                  rc *= 1./U_Exp(((y/sigma_y_C)^2.)/2.);
                 }
               if((x < 6.0*sigma_x_T)&&(y < 6.0*sigma_y_T))
                 {
                  rt = (beam.Number() - Nc_mean)/(2.*U_pi*sigma_x_T*sigma_y_T);
                  rt *= 1./U_Exp(((x/sigma_x_T)^2.)/2.);
                  rt *= 1./U_Exp(((y/sigma_y_T)^2.)/2.);
                 }
             }else
             {
               if((x < 6.0*sigma_x_C)&&(y < 6.0*sigma_y_C)&&(s < 6.0*sigma_s_C))
                 {
                  rc = Nc_mean/(2.*U_pi*sigma_x_C*sigma_y_C*
                      sqrt(beam[j][1]*beam[j][1]+beam[j][3]*beam[j][3]+1.0));
                  rc *= 1./U_Exp(((x/sigma_x_C)^2.)/2.);
                  rc *= 1./U_Exp(((y/sigma_y_C)^2.)/2.);
                  rc *= 1./U_Exp(((s/sigma_s_C)^2.)/2.);
                  rc *= U_Exp((((x*beam[j][1]+y*beam[j][3]-s)/sigma_s_C)^2.)/
                       (2.0 *(beam[j][1]*beam[j][1]+beam[j][3]*beam[j][3]+1.0)));
                 }
                 if((x < 6.0*sigma_x_T)&&(y < 6.0*sigma_y_T)&&(s < 6.0*sigma_s_T))
                 {
                  rt = (beam.Number() - Nc_mean)/(2.*U_pi*sigma_x_T*sigma_y_T*
                      sqrt(beam[j][1]*beam[j][1]+beam[j][3]*beam[j][3]+1.0));
                  rc *= 1./U_Exp(((x/sigma_x_T)^2.)/2.);
                  rc *= 1./U_Exp(((y/sigma_y_T)^2.)/2.);
                  rc *= 1./U_Exp(((s/sigma_s_T)^2.)/2.);
                  rc *= U_Exp((((x*beam[j][1]+y*beam[j][3]-s)/sigma_s_T)^2.)/
                       (2.0 *(beam[j][1]*beam[j][1]+beam[j][3]*beam[j][3]+1.0)));
                 }
             }
            */
            //lumi = beam.Emit[3] / beam.Number()  / ring.Trev * (rc + rt);
            cBeam.Sigma_x = sigma_x_C;
            cBeam.Sigma_y = sigma_y_C;
            cBeam.Sigma_s = sigma_s_C;
            cBeam.Emit[3] = CoreNumber * beam.Emit[3] / beam.Number();
            lumi = MalakhovDensity(beam[j]) / ring.Trev;
            cBeam.Sigma_x = sigma_x_T;
            cBeam.Sigma_y = sigma_y_T;
            cBeam.Sigma_s = sigma_s_T;
            cBeam.Emit[3] = (beam.Number() - CoreNumber) * beam.Emit[3] / beam.Number();
            lumi += MalakhovDensity(beam[j]) / ring.Trev;
         }
      }
      else
      {
         //lumi =  CollBeamDensity(beam[j])/ ring.Trev;
         lumi = MalakhovDensity(beam[j]) / ring.Trev;
      }
      if (beam.HourGlass && LuminosityModel != 0)
      {
         if (beam.benum == BUNCHED)
            lumi *= HourGlass();
         else if (beam.benum == BUCKET)
            lumi *= ring.Bucket.Luminosity(1, ring.Bucket.barrier.Number - 1,
                                           sqrt(iColl.Lattice.betax(m_) * iColl.Lattice.betay(m_)),
                                           (ring.Bucket.barrier[2].s2 - ring.Bucket.barrier[2].s1)(m_) / 2.);
      }
      if (Loss)
         beam.Loss(Lattice, j, (lumi * Points * time.dt * Cross)(U1_), false);
      Luminosity += lumi * beam.N_b * beam.Emit[3] / beam.Number();
   }
   /*
   if(LuminosityModel == 4)
   {
      Luminosity = 2.0*beam.N_b * beam.Emit[3] * beam.Emit[3] *
      (beam.Number() - Nc_mean)*Nc_mean/(beam.Number()*beam.Number())/
      ((((sigma_x_T*sigma_x_T + sigma_x_C*sigma_x_C) *
         (sigma_y_T*sigma_y_T + sigma_y_C*sigma_y_C))^0.5) *
         ring.Trev * U_pi* 2.0);
      Luminosity += beam.N_b * (beam.Emit[3]^2) * Nc_mean*Nc_mean/
      (beam.Number()*beam.Number())/
      (sigma_x_C*sigma_y_C * ring.Trev * U_pi* 4.0);
      Luminosity += beam.N_b * (beam.Emit[3]^2) * (beam.Number()-Nc_mean) *
      (beam.Number()-Nc_mean)/(beam.Number()*beam.Number())/
      (sigma_x_T*sigma_y_T * ring.Trev * U_pi* 4.0);
   }
   */
   //11.12.2007
   if (Include)
   {
      double r1, r2;
      r1 = (((0.25 * Diffusion * Diffusion * Kappa_x * Kappa_x * U_pi * U_pi) * time.dt * 2 *
             beam.Emit[0] / (Lattice.betax * ring.Trev)) ^
            0.5)(U1_);
      r2 = (((0.25 * Diffusion * Diffusion * Kappa_y * Kappa_y * U_pi * U_pi) * time.dt * 2 *
             beam.Emit[1] / (Lattice.betay * ring.Trev)) ^
            0.5)(U1_);
      if (CollBeam)
      {
         r1 *= (beam.Sigma[0] / cBeam.Sigma_x)(U1_);
         r2 *= (beam.Sigma[2] / cBeam.Sigma_y)(U1_);
      }
      for (int j = 0; j < beam.Number(); j++)
      {
         beam.Add(j, 1, r1 * xDistributor::Gaussian());
         beam.Add(j, 3, r2 * xDistributor::Gaussian());
      }
   }
}
//------------------------------------------------------------
//Local density for Gaussian beam 05.06
doubleU xColl::CollBeamDensity(double *X)
{
   doubleU r(m_ ^ -2);

   doubleU x(X[0], m_);
   doubleU y(X[2], m_);

   r = 0;
   if ((x < 6.0 * cBeam.Sigma_x) && (y < 6.0 * cBeam.Sigma_y))
   {
      r = cBeam.Emit[3] / (2. * U_pi * cBeam.Sigma_x * cBeam.Sigma_y);
      r *= 1. / U_Exp(((x / cBeam.Sigma_x) ^ 2.) / 2.);
      r *= 1. / U_Exp(((y / cBeam.Sigma_y) ^ 2.) / 2.);
   }

   return r;
}
//---------------------------------------------
//------------------------------------------------------------
//Local density for Gaussian beam by V.Malakhov 02.07
doubleU xColl::MalakhovDensity(double *X)
{
   double R;
   double v;
   double u;
   double beta_x;
   double beta_y;
   double sigma_x;
   double sigma_y;
   double sigma_z;
   double N;

   double C;

   v = iRing.Energy.Beta();
   u = cBeam.Energy.Beta();

   N = cBeam.Emit[3]();

   beta_x = cBeam.Lattice.betax(m_);
   beta_y = cBeam.Lattice.betay(m_);

   sigma_x = cBeam.Sigma_x(m_);
   sigma_y = cBeam.Sigma_y(m_);
   sigma_z = cBeam.Sigma_s(m_);

   C = N / (2.0 * M_PI * sqrt(M_PI) * sigma_x * sigma_y);

   double c[8] = {1.996040609e-4, 1.707798308e-2, 0.207802326, 0.661147013, 0.661147013, 0.207802326, 1.707798308e-2, 1.996040609e-4};
   double w[8] = {-2.9306374203, -1.98165675667, -1.1571937124, -0.3811869902, 0.3811869902, 1.1571937124, 1.98165675667, 2.9306374203};

   double Ax;
   double Ay;
   double Bx;
   double By;
   double powex;

   R = 0;

   for (int i = 0; i < 8; i++)
   {
      Ax = X[0] + X[1] * X[4] + X[1] * (v / (u + v)) * (sqrt(2.0) * sigma_z * w[i] - X[4]);
      Ay = X[2] + X[3] * X[4] + X[3] * (v / (u + v)) * (sqrt(2.0) * sigma_z * w[i] - X[4]);
      Bx = 1.0 + (u * u / ((u + v) * (u + v) * beta_x * beta_x)) * (sqrt(2.0) * sigma_z * w[i] - X[4]) * (sqrt(2.0) * sigma_z * w[i] - X[4]);
      By = 1.0 + (u * u / ((u + v) * (u + v) * beta_y * beta_y)) * (sqrt(2.0) * sigma_z * w[i] - X[4]) * (sqrt(2.0) * sigma_z * w[i] - X[4]);
      powex = Ax * Ax / (2.0 * sigma_x * sigma_x * Bx);
      powex += Ay * Ay / (2.0 * sigma_y * sigma_y * By);
      if (powex < 20.0)
         R += c[i] / (exp(powex) * sqrt(Bx * By));
   }

   R *= C;
   doubleU r(R, m_ ^ -2);
   return r;
}
//---------------------------------------------
void xColl::CalcLuminosity(xTime &time, xBeam &beam, xRing &ring)
{
   if (!Realtime)
      beam.Distribution(beam.Emit, Lattice);
   Realtime = false;

   iDraw.LuminosityTest.Reset(0);
   double keep_divisions = Divisions;
   bool keep_enebleloss = Loss;
   Loss = false;

   int percents = 0;
   for (int n = (int)From; n <= (int)Upto; n += (int)Steps)
   {
      Divisions = n;
      if (Divisions < 1)
         Divisions = 1;
      Kick(time, beam, ring);
      iDraw.LuminosityTest.Point(n, Luminosity((cm_ ^ -2) / s_));
      if (percents <= (100. * n / (Upto - From)))
      {
         Warning("Real-time Luminosity calculation. Left ", percents, "%");
         percents += 10;
      }
   }
   Divisions = (int)keep_divisions;
   Loss = keep_enebleloss;
   iDraw.LuminosityTest.Save();
}
//------------------------------------------------------------

void xColl::BeamBeam(xTime &time, xBeam &beam, xRing &ring)
{
   vectorU emit1;
   emit1[0]._(m_);
   emit1[1]._(m_);

   double A[6];
   doubleU l_density(m_ ^ -2);
   double ratio;

   doubleU sigma_x(m_);
   doubleU sigma_y(m_);

   int N = 0;
   double Nx = 0;
   double Ny = 0;
   //int Em;
   int ind = 0;

   if (Local_BB && iDynamics.Algoritm == 1) // if B-B is calculated via local density
   {
      for (int j = 0; j < 6; j++)
         A[j] = 0;
      l_density = beam.LocalDensity(A, Divisions, Lattice);
      //l_density = beam.LocalDensity_G(A, Lattice);
      ratio = beam.Emit[3]() / beam.Number();

      Kappa_x = (Lattice.betax * l_density * ratio / (4.)) *
                (ring.Energy.Z * ring.Energy.Z * U_e * U_e / (ring.Energy.A * U_amu * U_c * U_c)) *
                ((1. + ring.Energy.Beta2) / (ring.Energy.Gamma * ring.Energy.Beta));

      Kappa_y = (Lattice.betay * l_density * ratio / (4.)) *
                (ring.Energy.Z * ring.Energy.Z * U_e * U_e / (ring.Energy.A * U_amu * U_c * U_c)) *
                ((1. + ring.Energy.Beta2) / (ring.Energy.Gamma * ring.Energy.Beta));
   }
   else
   {
      if (CollBeam)
      {
         Kappa_x = (cBeam.Emit[3] * Lattice.betax / (4. * U_pi)) *
                   (ring.Energy.Z * cBeam.Energy.Z * U_e * U_e / (ring.Energy.A * U_amu * U_c * U_c)) *
                   ((ring.Energy.Beta * cBeam.Energy.Beta + 1) / (ring.Energy.Gamma * ring.Energy.Beta)) *
                   (1. / (cBeam.Sigma_x * (cBeam.Sigma_x + cBeam.Sigma_y)));

         Kappa_y = (cBeam.Emit[3] * Lattice.betay / (4. * U_pi)) *
                   (ring.Energy.Z * cBeam.Energy.Z * U_e * U_e / (ring.Energy.A * U_amu * U_c * U_c)) *
                   ((ring.Energy.Beta * cBeam.Energy.Beta + 1) / (ring.Energy.Gamma * ring.Energy.Beta)) *
                   (1. / (cBeam.Sigma_y * (cBeam.Sigma_x + cBeam.Sigma_y)));
      }
      else
      {
         if (iDynamics.Algoritm == 1)
         {
            int j;
            switch (BB_EmitDef)
            {
            case 0:
               emit1 = iBeam.Emit_RMS(Lattice);
               Nx = iBeam.Emit[3]();
               Ny = iBeam.Emit[3]();
               break;
            case 1:
               //          emit1 = iBeam.Emit_FWHM(Lattice);
               emit1 = iBeam.Emittance_FWHM(Lattice, 100);
               N = iBeam.Number();
               //Em = emit1[3]();
               Nx = iBeam.X_FWHM * emit1[3]() / N;
               Ny = iBeam.Y_FWHM * emit1[3]() / N;
               break;
            case 2:
               ind = (int)(floor(percent * iBeam.Number() / 100.) - 1);
               for (j = 0; j < 3; j++)
                  emit1[j] = iBeam.Inv(ind, j);
               //          emit1 = iBeam.EmitSort(percent, Lattice);
               emit1[3] = iBeam.Emit[3] * 2 * U_Ln(1. / (1 - (percent / 100.)));
               Nx = emit1[3]();
               Ny = emit1[3]();
               break;
            default:
               Warning("Wrong Emittance Definition");
            }
         }
         else //________________________________27.09.2005_____________________________//
         {
            emit1[0] = beam.Emit[0];
            emit1[1] = beam.Emit[1];
            Nx = beam.Emit[3]();
            Ny = beam.Emit[3]();
         }
         sigma_x = U_Pow((emit1[0] * Lattice.betax)(), 0.5);
         sigma_y = U_Pow((emit1[1] * Lattice.betay)(), 0.5);

         Kappa_x = (Nx * Lattice.betax / (4. * U_pi)) *
                   (ring.Energy.Z * ring.Energy.Z * U_e * U_e / (ring.Energy.A * U_amu * U_c * U_c)) *
                   ((1. + ring.Energy.Beta2) / (ring.Energy.Gamma * ring.Energy.Beta)) *
                   (1. / (sigma_x * (sigma_x + sigma_y)));

         Kappa_y = (Ny * Lattice.betay / (4. * U_pi)) *
                   (ring.Energy.Z * ring.Energy.Z * U_e * U_e / (ring.Energy.A * U_amu * U_c * U_c)) *
                   ((1. + ring.Energy.Beta2) / (ring.Energy.Gamma * ring.Energy.Beta)) *
                   (1. / (sigma_y * (sigma_x + sigma_y)));
      }
   }
}

//---------------------------------------------------------------------------
//- Effect Optical Stochastic Cooling
//---------------------------------------------------------------------------

xOptStoch::xOptStoch()
{
   //   ID = 9;
   R52._(m_);
   R56._(m_);
   Go._(s_ ^ -1);
   Zo._(m_);
   SigmaZ._(m_);
   lambda._(m_);
}

int xOptStoch::OnGet()
{
   R51 = Data.OnGet(80, 1, 1);
   R52 = Data.OnGet(80, 2, 1);
   R56 = Data.OnGet(80, 3, 1);
   Go = Data.OnGet(80, 4, 0);
   Zo = Data.OnGet(80, 5, 1);
   SigmaZ = Data.OnGet(80, 6, 1);
   lambda = Data.OnGet(80, 7, 1);

   return 0;
}

int xOptStoch::OnSet()
{
   Lattice.betax = Data.OnGet(81, 1, 1);
   Lattice.betay = Data.OnGet(81, 2, 1);
   Lattice.alfax = Data.OnGet(81, 3, 0);
   Lattice.alfay = Data.OnGet(81, 4, 0);
   Lattice.Dx = Data.OnGet(81, 5, 0);
   Lattice.Dy = Data.OnGet(81, 6, 0);
   Lattice.Dpx = Data.OnGet(81, 7, 0);
   Lattice.Dpy = Data.OnGet(81, 8, 0);

   return 0;
}

void xOptStoch::Kick(xTime &time, xBeam &beam, xRing &ring)
{
   doubleU phi;
   doubleU G;

   for (int j = 0; j < beam.Number(); j++)
   {
      phi = -(R51 * beam(j, 0) + R52 * beam(j, 1) + R56 * beam(j, 5)) * 2 * M_PI / lambda;
      G = Go * time.dt * U_Exp(-(((beam(j, 4) - Zo) / (SigmaZ * 2)) ^ 2));
      beam.Add(j, 5, G * U_Sin(phi));
   }
}

//---------------------------------------------------------------------------
//- Effect Laser Cooling
//---------------------------------------------------------------------------

xLaserCooling::xLaserCooling()
{
   ID = 10;
   OpticName = "XLCOOL";
   EffectName = "XLCOOL";

   ResonantFreq._(G_9 * Hz_);
   LaserFreq._(G_9 * Hz_);
   CoolLength._(1, m_);
   NuturalWidth._(285.7, M_6 * Hz_);
   Saturation._(1);
   SpotSize._(1, m_3 * m_);
   Delta._(G_9 * Hz_);
   To._(0, s_);
   DeltaStart1._(G_9 * Hz_);
   DeltaStop1._(G_9 * Hz_);
   DeltaStart2._(G_9 * Hz_);
   DeltaStop2._(G_9 * Hz_);
   Shift1._(0, m_3 * m_);
   Shift2._(0, m_3 * m_);
   Shift1fin._(0, m_3 * m_);
   Shift2fin._(0, m_3 * m_);
}

int xLaserCooling::OnGet()
{
   ResonantFreq = Data.OnGet(85, 1, 1);
   CoolLength = Data.OnGet(85, 3, 1);
   NuturalWidth = Data.OnGet(85, 4, 285.7);
   Saturation = Data.OnGet(85, 5, 1);
   SpotSize = Data.OnGet(85, 6, 1);
   Delta = Data.OnGet(85, 7, 1);

   LaserFreq = ResonantFreq + Delta;
   Data.OnSet(85, 2, LaserFreq(G_9 * Hz_));

   SweepOn1 = Data.OnGetB(87, 1, false);
   DeltaStart1 = Data.OnGet(87, 2, 0);
   ;
   DeltaStop1 = Data.OnGet(87, 3, 0);
   Shift1 = Data.OnGet(87, 4, 0);
   Shift1fin = Data.OnGet(87, 5, 0);

   SweepOn2 = Data.OnGetB(87, 6, false);
   DeltaStart2 = Data.OnGet(87, 7, 0);
   ;
   DeltaStop2 = Data.OnGet(87, 8, 0);
   Shift2 = Data.OnGet(87, 9, 0);
   Shift2fin = Data.OnGet(87, 10, 0);

   To = Data.OnGet(87, 11, 0);
   SweepStep2 = Data.OnGetI(87, 12, 100);

   return 0;
}

int xLaserCooling::OnSet()
{

   Lattice.betax = Data.OnGet(86, 1, 1);
   Lattice.betay = Data.OnGet(86, 2, 1);
   Lattice.alfax = Data.OnGet(86, 3, 0);
   Lattice.alfay = Data.OnGet(86, 4, 0);
   Lattice.Dx = Data.OnGet(86, 5, 0);
   Lattice.Dy = Data.OnGet(86, 6, 0);
   Lattice.Dpx = Data.OnGet(86, 7, 0);
   Lattice.Dpy = Data.OnGet(86, 8, 0);

   return 0;
}

vectorU xLaserCooling::GetF(xTime &t, U_Energy &e, vectorU ion)
{
   vectorU forces(U1_, m_ ^ -1);

   Delta = LaserFreq * iRing.Energy.Gamma *
               (1 - iRing.Energy.Beta * (1 + ion[5] / iRing.Energy.Gamma / iRing.Energy.Gamma)) -
           (ResonantFreq * iRing.Energy.Gamma * (1 - iRing.Energy.Beta));
   /*
   Delta = 2 * M_PI * U_c / WaveLength *
          ((1 + e.Beta * (1 + ion[5] - Sweep)) / (1 + e.Beta) - 1);
          ((1 + e.Beta * (1 + ion[5] - Sweep)) / (1 + e.Beta * (1 - Sweep)) - 1);
*/
   S = Saturation * U_Exp(-0.5 * (ion[0] * ion[0] + ion[2] * ion[2]) /
                          (SpotSize * SpotSize));

   forces[5] = 0.5 * U_hbar * ResonantFreq / U_c * NuturalWidth * S /
               (1 + S + ((2 * Delta / NuturalWidth) ^ 2)) /
               (e.Momentum * e.Velocity * e.Gamma);

   return forces;
}

vectorU xLaserCooling::GetForces(xTime &time, U_Energy &e, vectorU ion)
{
   vectorU ions, force(U1_, m_ ^ -1);

   if (SweepOn2)
   {
      ions = ion;
      if (To() == 0)
      {
         LaserFreq = DeltaStart2 + ResonantFreq;
         ions[0] -= Shift2;
      }
      else if (time.t < To)
      {
         LaserFreq = DeltaStart2 + ResonantFreq + ((DeltaStop2 - DeltaStart2) * (time.t / To)());
         ions[0] -= (Shift2fin - Shift2) * (time.t / To)() + Shift2;
      }
      else
      {
         LaserFreq = DeltaStop2 + ResonantFreq;
         ions[0] -= Shift2fin;
      }
      force -= GetF(time, iRing.Energy, ions);
   }

   if (SweepOn1)
   {
      ions = ion;
      if (To() == 0)
      {
         LaserFreq = DeltaStart1 + ResonantFreq;
         ions[0] -= Shift1;
      }
      else if (time.t < To)
      {
         LaserFreq = DeltaStart1 + ResonantFreq + ((DeltaStop1 - DeltaStart1) * (time.t / To)());
         ions[0] -= (Shift1fin - Shift1) * (time.t / To)() + Shift1;
      }
      else
      {
         LaserFreq = DeltaStop1 + ResonantFreq;
         ions[0] -= Shift1fin;
      }
      force += GetF(time, iRing.Energy, ions);
   }
   else
   {
      force = GetF(time, iRing.Energy, ion);
   }
   return force;
}

void xLaserCooling::Kick(xTime &time, xBeam &beam, xRing &ring)
{
   vectorU force;
   int div = beam.Number() / 10;
   hyst.SetNumber(div);
   for (int k = 0; k < div; k++)
      hyst[k] = 0;
   TotalLight = 0;

   for (int j = 0; j < beam.Number(); j++)
   {
      force = GetForces(time, ring.Energy, beam(j));
      beam.Add(j, 5, force[5] * CoolLength * time.dt / ring.Trev);
      int index = Round((beam(j, 0) * div / (SpotSize * 10))(U1_)) + div / 2;
      if (index < 0)
         index = 0;
      if (index > div - 1)
         index = div - 1;
      hyst[index] += force[5](m_ ^ -1);
      TotalLight += force[5](m_ ^ -1);
   }
}
