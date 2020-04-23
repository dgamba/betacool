//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xEbeam.h"
#include "xDraw.h"
//---------------------------------------------------------------------------

xFrParam::xFrParam()
{
   V_tr_e._(m_ / s_);
   // ----- 05.06  for 3D force
   V_tr_x._(m_ / s_);
   V_tr_y._(m_ / s_);
   //-----
   V_long_e._(m_ / s_);
   Ttemp._(eV_);
   Ttemp_centre._(eV_);
   Ltemp._(eV_);
   TempEff._(eV_);
   Theta_Eff._(U1_);
   V_eff_e._(m_ / s_);
   Smoos._(U1_);
   mfield._(k_3 * G_);
   n_e._(1 / (m_ ^ 3));
   tau._(s_);
   T_plasma._(s_);
   //Undulator parameters
   lambda._(cm_);
   B_field._(G_);
   r_0._(u_6 * m_);
   Theta_U._(U1_);
   V_und._(m_ / s_);
}

xEbeam::xEbeam(xRing &ring)
{
   pBeam = ring.pBeam;
   LinearRates[0]._(Hz_);
   LinearRates[1]._(Hz_);
   LinearRates[2]._(Hz_);
   e_Energy.Massa._(electron_);
   pZ = &ring.Energy.Z;
   pA = &ring.Energy.A;
   pBeta = &ring.Energy.Beta;
   e_emit_tr._(m_);
   e_dpp._(U1_);
   CoLength._(m_);
   bradius._(cm_);
   bcurrent._(A_);
   dVdr._(s_ ^ -1);
   for (int i = 0; i < 6; i += 2)
   {
      Initial[i]._(0, m_);
      Initial[i + 1]._(0, U1_);
      Final[i]._(0, m_);
      Final[i + 1]._(0, U1_);
      Shift[i]._(0, m_);
      Shift[i + 1]._(0, U1_);
   }

   // 04.12.06
   Bdistance._(m_);
   //
   size_x._(cm_);
   size_y._(cm_);
   size_s._(cm_);
   //dist_uni._(cm_);
   Ne_uni._(U1_);
   Ie_uni._(A_);

   //        06.07
   Gsigma[0]._(m_);  //RMS horizontal dimension
   Gsigma[1]._(U1_); //RMS horizontal angular spread
   Gsigma[2]._(m_);  //RMS vertical dimension
   Gsigma[3]._(U1_); //RMS vertical angular spread
   Gsigma[4]._(m_);  //RMS bunch length
   Gsigma[5]._(U1_); //RMS momentum spread
   Ne._(U1_);
   Ieb._(A_);
   //
   sigma_x_cil._(cm_);
   sigma_y_cil._(cm_);
   Ne_cil._(1 / cm_);
   Ie_cil._(A_);
   r_hole._(cm_);
   n_hole._(1 / (cm_ ^ 3));
   r_circle._(cm_);
   n_circle._(1 / (cm_ ^ 3));
   I_circle._(A_);
   I_hocle._(A_);
   I_hole._(A_);
   Asigma[0]._(m_);  //RMS horizontal dimension
   Asigma[1]._(U1_); //RMS horizontal angular spread
   Asigma[2]._(m_);  //RMS vertical dimension
   Asigma[3]._(U1_); //RMS vertical angular spread
   Asigma[4]._(m_);  //RMS bunch length
   Asigma[5]._(U1_); //RMS momentum spread
   //Adist._(cm_);                                  //Distance between ion and electron bunch centers
   ANe._(U1_); //Electron number in the bunch
   //m_x._(U1_);
   //m_y._(U1_);
   //m_p._(U1_);                                     //components of mean electron velocity

   //    06.07
   // parameters of parabolic cylinder
   pradius._(cm_);
   pcurrent._(A_);
   dVdrp._(s_ ^ -1);
   n_centre._(1 / (cm_ ^ 3));
   ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07
   dVdrd._(s_ ^ -1);
   dradius._(cm_);
}
//---------------------------------------------------------------------------

int xEbeam::OnGet()
{
   F.Z = (*pZ);
   F.A = (*pA);

   e_Energy.A = 1;
   e_Energy.Z = 1;
   e_Energy.Beta = (*pBeta);
   e_Energy.Set(U_BETA);

   CoLength = Data.OnGet(50, 3, 1);
   F.mfield = Data.OnGet(50, 4, 1);

   bradius = Data.OnGet(55, 2, 1);
   bcurrent = Data.OnGet(55, 3, 1);
   Neutralization = Data.OnGet(55, 4, 100) / 100.;
   //dpshift        = Data.OnGet(55,5,0);
   //xshift         = Data.OnGet(55,6,0);
   dVdr = Data.OnGet(55, 7, 0);
   LongPulse = Data.OnGet(55, 8, 0);
   ;
   PulseFrom = Data.OnGet(55, 9, -0.1);
   ;
   PulseUpTo = Data.OnGet(55, 10, 0.1);
   ;

   emit = Data.OnGetI(60, 2, 1);
   temp = Data.OnGetI(60, 3, 0);
   velocity = Data.OnGetI(60, 8, 0);
   e_emit_tr = Data.OnGet(60, 4, 1);
   F.Ttemp = Data.OnGet(60, 5, 0.1);
   e_dpp = Data.OnGet(60, 6, 1);
   F.Ltemp = Data.OnGet(60, 7, 0.001);
   F.V_tr_e = Data.OnGet(60, 9, 100000.);
   F.V_long_e = Data.OnGet(60, 10, 10000.);

   EnableShifts = Data.OnGetB(53, 1, false);
   EnableFinal = Data.OnGetB(53, 2, false);
   if (EnableShifts)
   {
      Initial[0] = Data.OnGet(53, 3, 0);
      Final[0] = Data.OnGet(53, 4, 0);
      Initial[2] = Data.OnGet(53, 5, 0);
      Final[2] = Data.OnGet(53, 6, 0);
      Initial[4] = Data.OnGet(53, 7, 0);
      Final[4] = Data.OnGet(53, 8, 0);
      Initial[5] = Data.OnGet(53, 9, 0);
      Final[5] = Data.OnGet(53, 10, 0);
      for (int e = 0; e < 6; e++)
      {
         if (!EnableFinal)
            Final[e] = Initial[e];
      }
   }
   // 04.12.06
   Bdistance = Data.OnGet(50, 8, 0.1);
   BunchNumber = Data.OnGetI(50, 7, 1);

   UsePeriod = Data.OnGetB(53, 11, false);
   PaintingPeriod = Data.OnGetI(53, 12, 1);
   FromFile = Data.OnGetB(53, 13, false);
   if (FromFile)
   {
      Solenoid.AddSeparators(",");
      Solenoid.Load(Data.OnGetC(53, 14, "none"));
   }
   PaintingTable = Data.OnGetB(53, 15, false);
   if (PaintingTable)
   {
      PaintingData.AddSeparators(",");
      PaintingData.Load(Data.OnGetC(53, 16, "none"));
   }
   ScallingFactor = Data.OnGetI(53, 17, 1);
   if (ScallingFactor == 0)
      ScallingFactor = 1;

   model = Data.OnGetI(55, 1, 0);

   size_x = Data.OnGet(56, 1, 1);
   size_y = Data.OnGet(56, 2, 1);
   size_s = Data.OnGet(56, 3, 1);
   //dist_uni= Data.OnGet(56,4,1);
   Ne_uni = Data.OnGet(56, 5, 1);

   //       06.07
   Gsigma[0] = Data.OnGet(57, 1, 0.01);
   Gsigma[2] = Data.OnGet(57, 2, 0.01);
   Gsigma[4] = Data.OnGet(57, 3, 0.01);
   Gsigma[1] = Data.OnGet(57, 8, 0.00001);
   Gsigma[3] = Data.OnGet(57, 9, 0.00001);
   Gsigma[5] = Data.OnGet(57, 10, 0.0001);
   From_model = Data.OnGetB(57, 11, true);
   Ne = Data.OnGet(57, 5, 1);

   sigma_x_cil = Data.OnGet(58, 1, 1);
   sigma_y_cil = Data.OnGet(58, 2, 1);
   Ne_cil = Data.OnGet(58, 3, 1);

   r_hole = Data.OnGet(59, 1, 1);
   n_hole = Data.OnGet(59, 2, 1);
   r_circle = Data.OnGet(59, 3, 1);
   n_circle = Data.OnGet(59, 4, 1);
   SpaceCharge = Data.OnGetB(59, 6, false);

   F.Smoos = Data.OnGet(61, 1, 1);
   spread_temp = Data.OnGetI(61, 2, 1);
   F.Theta_Eff = Data.OnGet(61, 4, 1);
   F.TempEff = Data.OnGet(61, 5, 1);

   F.undulator = Data.OnGetI(74, 1, 0);
   F.lambda = Data.OnGet(74, 2, 10);
   F.B_field = Data.OnGet(74, 3, 50);

   // 06.07
   BiGauss = Data.OnGetB(54, 27, false);
   Uniform = Data.OnGetB(54, 30, false);

   sigma_core = Data.OnGet(54, 29, 0.0);
   GaussN = Data.OnGetB(75, 5, false);
   AFromFile = Data.OnGetB(54, 1, false);
   ANe = Data.OnGet(54, 11, 1);
   //Adist      = Data.OnGet (54,10,1);
   Minst.Number = Data.OnGetI(54, 13, 100);
   //      06.07
   N_hor = Data.OnGetI(54, 25, 10);
   N_long = Data.OnGetI(54, 26, 10);
   Box = Data.OnGet(54, 14, 0.5);
   Slice = Data.OnGetB(54, 15, 1);
   From = Data.OnGet(54, 16, 0);
   Upto = Data.OnGet(54, 17, 0.01);

   if (model == 5)
   {
      double sum[6];
      doubleU lengthU(1, m_);

      if (AFromFile)
      {
         EArray.AddSeparators(",");
         EArray.Load(Data.OnGetC(54, 3, "NoFile"));
         Bunch.Number = (int)EArray[0][0];
         AEnergy = EArray[0][3];
         AGamma = (AEnergy + 0.5110034) / (AEnergy + 2.0 * 0.5110034);
      }
      else
      {
         for (int x = 0; x < 6; x++)
            Asigma[x] = Data.OnGet(54, x + 4, 1);
         Bunch.Number = Data.OnGetI(54, 12, 1000);
         N_core = Data.OnGetI(54, 28, 0) * Bunch.Number;
      }

      if (Minst.Number < 1)
         Minst.Number = 1;
      if (Minst.Number > Bunch.Number)
         Minst.Number = Bunch.Number;

      Minst.SetNumber(Minst.Number);
      Bunch.SetNumber(Bunch.Number, 6);
      Dist2.SetNumber(Bunch.Number);

      if (AFromFile)
      {
         for (int j = 0; j < 6; j++)
            sum[j] = 0;

         for (int k = 1; k <= Bunch.Number; k++)
            for (int j = 0; j < 6; j++)
               sum[j] += EArray[k][j];
         for (int j = 0; j < 6; j++)
            sum[j] /= Bunch.Number;

         for (int k = 1; k <= Bunch.Number; k++)
            for (int j = 0; j < 6; j++)
               Bunch[k - 1][j] = EArray[k][j] - sum[j];

         for (int k = 0; k < Bunch.Number; k++)
         {
            Bunch[k][4] *= 0.4263 / 360.0;
            Bunch[k][5] *= AGamma / AEnergy;
         }
         for (int j = 0; j < 6; j++)
            sum[j] = 0;

         for (int k = 0; k < Bunch.Number; k++)
            for (int j = 0; j < 6; j++)
               sum[j] += Bunch[k][j] * Bunch[k][j];
         for (int j = 0; j < 6; j++)
            sum[j] /= Bunch.Number;

         Asigma[0] = sqrt(sum[0]) * lengthU;
         Asigma[1] = sqrt(sum[1]);
         Asigma[2] = sqrt(sum[2]) * lengthU;
         Asigma[3] = sqrt(sum[3]);
         Asigma[4] = sqrt(sum[4]) * lengthU;
         Asigma[5] = sqrt(sum[5]);
      }
      else
      {
         if (BiGauss)
         {
            for (int k = 0; k < N_core; k++)
            {
               Bunch[k][0] = sigma_core * Asigma[0](m_) * xDistributor::Gaussian();
               Bunch[k][1] = sigma_core * Asigma[1](U1_) * xDistributor::Gaussian();
               Bunch[k][2] = sigma_core * Asigma[2](m_) * xDistributor::Gaussian();
               Bunch[k][3] = sigma_core * Asigma[3](U1_) * xDistributor::Gaussian();
               Bunch[k][4] = sigma_core * Asigma[4](m_) * xDistributor::Gaussian();
               Bunch[k][5] = sigma_core * Asigma[5](U1_) * xDistributor::Gaussian();
            }
            for (int k = N_core; k < Bunch.Number; k++)
            {
               Bunch[k][0] = Asigma[0](m_) * xDistributor::Gaussian();
               Bunch[k][1] = Asigma[1](U1_) * xDistributor::Gaussian();
               Bunch[k][2] = Asigma[2](m_) * xDistributor::Gaussian();
               Bunch[k][3] = Asigma[3](U1_) * xDistributor::Gaussian();
               Bunch[k][4] = Asigma[4](m_) * xDistributor::Gaussian();
               Bunch[k][5] = Asigma[5](U1_) * xDistributor::Gaussian();
            }
         }
         else
         {
            if (Uniform)
            {
               for (int k = 0; k < Bunch.Number; k++)
               {
                  Bunch[k][0] = 2.0 * Asigma[0](m_) * xDistributor::Flatten();
                  Bunch[k][1] = Asigma[1](U1_) * xDistributor::Gaussian();
                  Bunch[k][2] = 2.0 * Asigma[2](m_) * xDistributor::Flatten();
                  Bunch[k][3] = Asigma[3](U1_) * xDistributor::Gaussian();
                  Bunch[k][4] = Asigma[4](m_) * xDistributor::Gaussian();
                  Bunch[k][5] = Asigma[5](U1_) * xDistributor::Gaussian();
               }
            }
            else
            {
               for (int k = 0; k < Bunch.Number; k++)
               {
                  Bunch[k][0] = Asigma[0](m_) * xDistributor::Gaussian();
                  Bunch[k][1] = Asigma[1](U1_) * xDistributor::Gaussian();
                  Bunch[k][2] = Asigma[2](m_) * xDistributor::Gaussian();
                  Bunch[k][3] = Asigma[3](U1_) * xDistributor::Gaussian();
                  Bunch[k][4] = Asigma[4](m_) * xDistributor::Gaussian();
                  Bunch[k][5] = Asigma[5](U1_) * xDistributor::Gaussian();
               }
            }
         }
      }
      Local = Bunch;
      LRF_PRF(e_Energy, Local);
   }

   // parabolic cylinder
   pradius = Data.OnGet(69, 1, 1);
   pcurrent = Data.OnGet(69, 2, 1);
   dVdrp = Data.OnGet(69, 3, 0);

   ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07
   if (model == 7)
   {
      Density.SetSeparators("\t,");
      if (!Density.Load(Data.OnGetC(69, 5, "NoFile")))
         return 1;
      dradius = Density[Density.Row() - 1][0];
   }
   dVdrd = Data.OnGet(69, 7, 0);

   return 0;
}
//---------------------------------------------------------------------------
int xEbeam::OnSet()
{
   doubleU M_e(U_me * U_c * U_c);

   if (emit == 1)
   {
      switch (model)
      {
      case 0:
         F.Ttemp = M_e * 4 * e_emit_tr * e_emit_tr / (bradius * bradius);
         break;
      case 1:
         F.Ttemp = M_e * 4 * e_emit_tr * e_emit_tr / (size_x * size_y);
         break;
      case 2:
         F.Ttemp = M_e * e_emit_tr * e_emit_tr / (sigma_x_cil * sigma_y_cil);
         break;
      case 3:
         //         06.07
         if (From_model)
         {
            F.Ttemp = M_e * e_emit_tr * e_emit_tr / (Gsigma[0] * Gsigma[2]);
            Gsigma[1] = ((F.Ttemp / M_e) ^ 0.5) / (e_Energy.Beta * e_Energy.Gamma);
            Gsigma[3] = Gsigma[1];
         }
         else
         {
            F.Ttemp = M_e * e_Energy.Beta * e_Energy.Beta * e_Energy.Gamma *
                      e_Energy.Gamma * 0.5 * ((Gsigma[1] * Gsigma[1]) + (Gsigma[3] * Gsigma[3]));
            e_emit_tr = (F.Ttemp * Gsigma[0] * Gsigma[2] / M_e) ^ 0.5;
            e_dpp = Gsigma[5];
         }

         break;
      case 4:
         F.Ttemp = M_e * 4 * e_emit_tr * e_emit_tr / (r_circle * r_circle);
         break;
      //   05.06
      case 5:
         F.Ttemp = M_e * e_Energy.Beta * e_Energy.Beta * e_Energy.Gamma *
                   e_Energy.Gamma * 0.5 * ((Asigma[1] * Asigma[1]) + (Asigma[3] * Asigma[3]));
         e_emit_tr = (F.Ttemp * Asigma[0] * Asigma[2] / M_e) ^ 0.5;
         e_dpp = Asigma[5];
         break;
      //  06.07
      case 6:
         F.Ttemp = M_e * 4 * e_emit_tr * e_emit_tr / (pradius * pradius);
         break;
      }
      F.Ltemp = M_e * e_dpp * e_dpp * e_Energy.Beta * e_Energy.Beta;
      F.V_tr_e = ((F.Ttemp / U_me) ^ 0.5);
      F.V_long_e = ((F.Ltemp / U_me) ^ 0.5);

      Data.OnSet(60, 5, F.Ttemp);
      Data.OnSet(60, 7, F.Ltemp);
      Data.OnSet(60, 9, F.V_tr_e);
      Data.OnSet(60, 10, F.V_long_e);
   }
   else if (temp == 1)
   {
      switch (model)
      {
      case 0:
      case 1:
         e_emit_tr = (F.Ttemp * bradius * bradius / (M_e * 4)) ^ 0.5;
         break;
      case 2:
         e_emit_tr = (F.Ttemp * sigma_x_cil * sigma_x_cil / M_e) ^ 0.5;
         break;
      case 3:
         if (From_model)
         {
            e_emit_tr = (F.Ttemp * Gsigma[0] * Gsigma[0] / M_e) ^ 0.5;
            Gsigma[1] = ((F.Ttemp / M_e) ^ 0.5) / (e_Energy.Beta * e_Energy.Gamma);
            Gsigma[3] = Gsigma[1];
         }
         else
         {
            F.Ttemp = M_e * e_Energy.Beta * e_Energy.Beta * e_Energy.Gamma *
                      e_Energy.Gamma * 0.5 * ((Gsigma[1] * Gsigma[1]) + (Gsigma[3] * Gsigma[3]));
            F.Ltemp = M_e * Gsigma[5] * Gsigma[5] * e_Energy.Beta * e_Energy.Beta;
            e_emit_tr = (F.Ttemp * Gsigma[0] * Gsigma[2] / M_e) ^ 0.5;
         }
         break;
      case 4:
         e_emit_tr = (F.Ttemp * r_circle * r_circle / (M_e * 4)) ^ 0.5;
         break;
      //  05.06
      case 5:
         F.Ttemp = M_e * e_Energy.Beta * e_Energy.Beta * e_Energy.Gamma *
                   e_Energy.Gamma * 0.5 * ((Asigma[1] * Asigma[1]) + (Asigma[3] * Asigma[3]));
         F.Ltemp = M_e * Asigma[5] * Asigma[5] * e_Energy.Beta * e_Energy.Beta;
         e_emit_tr = (F.Ttemp * Asigma[0] * Asigma[2] / M_e) ^ 0.5;
         break;
      //  06.07
      case 6:
         e_emit_tr = (F.Ttemp * pradius * pradius / (M_e * 4)) ^ 0.5;
         break;
      }
      e_dpp = (F.Ltemp / (M_e * e_Energy.Beta * e_Energy.Beta)) ^ 0.5;
      F.V_tr_e = ((F.Ttemp / U_me) ^ 0.5);
      F.V_long_e = ((F.Ltemp / U_me) ^ 0.5);

      Data.OnSet(60, 4, e_emit_tr);
      Data.OnSet(60, 6, e_dpp);
      Data.OnSet(60, 9, F.V_tr_e);
      Data.OnSet(60, 10, F.V_long_e);
   }
   else if (velocity == 1)
   {

      F.Ttemp = U_me * F.V_tr_e * F.V_tr_e;
      F.Ltemp = U_me * F.V_long_e * F.V_long_e;
      switch (model)
      {
      case 0:
      case 1:
         e_emit_tr = (F.Ttemp * bradius * bradius / (M_e * 4)) ^ 0.5;
         break;
      case 2:
         e_emit_tr = (F.Ttemp * sigma_x_cil * sigma_x_cil / M_e) ^ 0.5;
         break;
      case 3:
         if (From_model)
         {
            e_emit_tr = (F.Ttemp * Gsigma[0] * Gsigma[0] / M_e) ^ 0.5;
            Gsigma[1] = ((F.Ttemp / M_e) ^ 0.5) / (e_Energy.Beta * e_Energy.Gamma);
            Gsigma[3] = Gsigma[1];
         }
         else
         {
            F.Ttemp = M_e * e_Energy.Beta * e_Energy.Beta * e_Energy.Gamma *
                      e_Energy.Gamma * 0.5 * ((Gsigma[1] * Gsigma[1]) + (Gsigma[3] * Gsigma[3]));
            F.Ltemp = M_e * Gsigma[5] * Gsigma[5] * e_Energy.Beta * e_Energy.Beta;
            e_emit_tr = (F.Ttemp * Gsigma[0] * Gsigma[2] / M_e) ^ 0.5;
            F.V_tr_e = ((F.Ttemp / U_me) ^ 0.5);
            F.V_long_e = ((F.Ltemp / U_me) ^ 0.5);
         }
         break;
      case 4:
         e_emit_tr = (F.Ttemp * r_circle * r_circle / (M_e * 4)) ^ 0.5;
         break;
      case 5:
         F.Ttemp = M_e * e_Energy.Beta * e_Energy.Beta * e_Energy.Gamma *
                   e_Energy.Gamma * 0.5 * ((Asigma[1] * Asigma[1]) + (Asigma[3] * Asigma[3]));
         F.Ltemp = M_e * Asigma[5] * Asigma[5] * e_Energy.Beta * e_Energy.Beta;
         e_emit_tr = (F.Ttemp * Asigma[0] * Asigma[2] / M_e) ^ 0.5;
         F.V_tr_e = ((F.Ttemp / U_me) ^ 0.5);
         F.V_long_e = ((F.Ltemp / U_me) ^ 0.5);
         break;
      //  06.07
      case 6:
         e_emit_tr = (F.Ttemp * pradius * pradius / (M_e * 4)) ^ 0.5;
         break;
      }
      e_dpp = (F.Ltemp / (M_e * e_Energy.Beta * e_Energy.Beta)) ^ 0.5;

      Data.OnSet(60, 4, e_emit_tr);
      Data.OnSet(60, 6, e_dpp);
      Data.OnSet(60, 5, F.Ttemp);
      Data.OnSet(60, 7, F.Ltemp);
   }

   F.Ttemp_centre = F.Ttemp;
   F.V_tr_e = ((F.Ttemp / U_me) ^ 0.5);
   F.V_long_e = ((F.Ltemp / U_me) ^ 0.5);
   F.tau = CoLength / (e_Energy.Gamma * e_Energy.Velocity);
   //Data[Tag][8] = F.tau();
   //---  05.06
   F.V_tr_x = F.V_tr_e;
   F.V_tr_y = F.V_tr_e;
   //----
   //     06.07
   F.n_e = bcurrent / (U_pi * bradius * bradius * e_Energy.Gamma * e_ * e_Energy.Velocity); //in PRF
   F.T_plasma = 2 * U_pi * ((U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);
   //Data[Tag][9] = F.T_plasma();

   Ie_uni = e_Energy.Velocity * U_e * Ne_uni / (size_s * ((2. * U_pi) ^ 0.5));
   Data.OnSet(56, 6, Ie_uni);

   Ieb = e_Energy.Velocity * U_e * Ne / (Gsigma[4] * ((2. * U_pi) ^ 0.5));
   Data.OnSet(57, 6, Ieb);
   if (model == 3)
   {
      F.n_e = Ne / (((2 * U_pi) ^ 1.5) * Gsigma[0] * Gsigma[2] * Gsigma[4] * e_Energy.Gamma);
      Data.OnSet(57, 7, F.n_e);
   }
   Ie_cil = e_Energy.Velocity * U_e * Ne_cil;
   Data.OnSet(58, 4, Ie_cil);

   I_circle = e_Energy.Velocity * U_e * U_pi * r_circle * r_circle * n_circle;
   I_hocle = e_Energy.Velocity * U_e * U_pi * r_hole * r_hole * n_circle;
   I_hole = e_Energy.Velocity * U_e * U_pi * r_hole * r_hole * n_hole;
   Data.OnSet(59, 5, I_circle - I_hocle + I_hole);

   if (spread_temp == 1)
   {
      F.TempEff = F.Theta_Eff * F.Theta_Eff * M_e * e_Energy.Gamma * e_Energy.Gamma;
      Data.OnSet(61, 5, F.TempEff);
   }
   else
   {
      F.Theta_Eff = (F.TempEff / (M_e * e_Energy.Gamma * e_Energy.Gamma)) ^ 0.5;
      Data.OnSet(61, 4, F.Theta_Eff);
   }
   //------------Undulator
   F.Theta_U = U_e * F.B_field * F.lambda / (2. * U_pi * M_e * e_Energy.Gamma * e_Energy.Beta);
   F.r_0 = F.Theta_U * F.lambda / (2. * U_pi);
   F.V_und = e_Energy.Gamma * e_Energy.Beta * U_c * F.Theta_U;
   Data.OnSet(74, 4, F.r_0);

   //   05.06
   //m_x = 0;
   //m_y = 0;
   //m_p = 0;

   if (model == 5)
   {
      for (int x = 0; x < 6; x++)
         Data.OnSet(54, x + 4, Asigma[x]);
      Data.OnSet(54, 12, Bunch.Number);

      Data.OnSet(60, 4, e_emit_tr);
      Data.OnSet(60, 6, e_dpp);
      Data.OnSet(60, 5, F.Ttemp);
      Data.OnSet(60, 7, F.Ltemp);
      Data.OnSet(60, 9, F.V_tr_e);
      Data.OnSet(60, 10, F.V_long_e);
   }
   //     06.07
   if ((model == 3) && (!From_model))
   {
      Data.OnSet(60, 4, e_emit_tr);
      Data.OnSet(60, 6, e_dpp);
      Data.OnSet(60, 5, F.Ttemp);
      Data.OnSet(60, 7, F.Ltemp);
      Data.OnSet(60, 9, F.V_tr_e);
      Data.OnSet(60, 10, F.V_long_e);
   }

   if (model == 6)
   {
      ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07
      //       n_centre = 3.0* pcurrent /(U_pi* pradius*pradius* e_Energy.Gamma* e_ *e_Energy.Velocity);  //in PRF
      n_centre = 2.0 * pcurrent / (U_pi * pradius * pradius * e_Energy.Gamma * e_ * e_Energy.Velocity); //in PRF
      F.n_e = n_centre;
      Data.OnSet(69, 4, F.n_e);
   }
   if (model == 7)
   {
      doubleU I(0, A_);
      doubleU J(A_ / (cm_ ^ 2));
      doubleU r0(m_3 * m_), r1(m_3 * m_);
      r0 = Density[0][0];
      J = Density[0][1];
      I = J * U_pi * r0 * r0;
      Data.OnSet(69, 8, (J / (e_Energy.Gamma * e_ * e_Energy.Velocity))(m_ ^ -3));

      for (int i = 1; i < Density.Row(); i++)
      {
         r0 = Density[i - 1][0];
         r1 = Density[i][0];
         J = Density[i][1];
         I += J * U_pi * (r1 * r1 - r0 * r0);
      }
      Data.OnSet(69, 6, I);
   }

   if (model == 5)
   {
      vectorU ion_new(m_, U1_);
      ion_new = Ibeam2Ebeam(iTime, ion_new);
      if (inside)
      {
         Data.OnSet(54, 18, F.n_e);
         Data.OnSet(54, 19, Mean[1]);
         Data.OnSet(54, 20, Mean[3]);
         Data.OnSet(54, 21, Mean[5]);
         Data.OnSet(54, 22, sqrt(Mean2[1]));
         Data.OnSet(54, 23, sqrt(Mean2[3]));
         Data.OnSet(54, 24, sqrt(Mean2[5]));
      }
      else
      {
         Data.OnSet(54, 18, 0.0);
         Data.OnSet(54, 19, 0.0);
         Data.OnSet(54, 20, 0.0);
         Data.OnSet(54, 21, 0.0);
         Data.OnSet(54, 22, 0.0);
         Data.OnSet(54, 23, 0.0);
         Data.OnSet(54, 24, 0.0);
      }
   }
   return 0;
}

int xEbeam::OnRun()
{
   CoLength = Data.OnGet(50, 3, 1);
   F.mfield = Data.OnGet(50, 4, 1);

   bradius = Data.OnGet(55, 2, 1);
   bcurrent = Data.OnGet(55, 3, 1);
   Neutralization = Data.OnGet(55, 4, 100) / 100.;
   //dpshift        = Data.OnGet(55,5,0);
   //xshift         = Data.OnGet(55,6,0);
   LongPulse = Data.OnGet(55, 8, 0);
   ;
   PulseFrom = Data.OnGet(55, 9, -0.1);
   ;
   PulseUpTo = Data.OnGet(55, 10, 0.1);
   ;

   //   emit_temp = Data.OnGetI(60,2,1);
   e_emit_tr = Data.OnGet(60, 4, 1);
   F.Ttemp = Data.OnGet(60, 5, 0.1);
   e_dpp = Data.OnGet(60, 6, 1);
   F.Ltemp = Data.OnGet(60, 7, 0.001);

   EnableShifts = Data.OnGetB(53, 1, false);
   EnableFinal = Data.OnGetB(53, 2, false);
   if (EnableShifts)
   {
      Initial[0] = Data.OnGet(53, 3, 0);
      Final[0] = Data.OnGet(53, 4, 0);
      Initial[2] = Data.OnGet(53, 5, 0);
      Final[2] = Data.OnGet(53, 6, 0);
      Initial[4] = Data.OnGet(53, 7, 0);
      Final[4] = Data.OnGet(53, 8, 0);
      Initial[5] = Data.OnGet(53, 9, 0);
      Final[5] = Data.OnGet(53, 10, 0);
      for (int e = 0; e < 6; e++)
      {
         if (!EnableFinal)
            Final[e] = Initial[e];
      }
   }
   UsePeriod = Data.OnGetB(53, 11, false);
   PaintingPeriod = Data.OnGetI(53, 12, 1);
   FromFile = Data.OnGetB(53, 13, false);
   if (FromFile)
   {
      Solenoid.AddSeparators(",");
      Solenoid.Load(Data.OnGetC(53, 14, "none"));
   }
   PaintingTable = Data.OnGetB(53, 15, false);
   if (PaintingTable)
   {
      PaintingData.AddSeparators(",");
      PaintingData.Load(Data.OnGetC(53, 16, "none"));
   }
   ScallingFactor = Data.OnGetI(53, 17, 1);
   if (ScallingFactor == 0)
      ScallingFactor = 1;

   model = Data.OnGetI(55, 1, 0);
   size_x = Data.OnGet(56, 1, 1);
   size_y = Data.OnGet(56, 2, 1);
   size_s = Data.OnGet(56, 3, 1);
   //dist_uni= Data.OnGet(56,4,1);
   Ne_uni = Data.OnGet(56, 5, 1);

   //         06.07
   Gsigma[0] = Data.OnGet(57, 1, 0.01);
   Gsigma[2] = Data.OnGet(57, 2, 0.01);
   Gsigma[4] = Data.OnGet(57, 3, 0.01);
   Gsigma[1] = Data.OnGet(57, 8, 0.00001);
   Gsigma[3] = Data.OnGet(57, 9, 0.00001);
   Gsigma[5] = Data.OnGet(57, 10, 0.0001);
   From_model = Data.OnGetB(57, 11, true);
   Ne = Data.OnGet(57, 5, 1);
   //
   sigma_x_cil = Data.OnGet(58, 1, 1);
   sigma_y_cil = Data.OnGet(58, 2, 1);
   Ne_cil = Data.OnGet(58, 3, 1);

   r_hole = Data.OnGet(59, 1, 1);
   n_hole = Data.OnGet(59, 2, 1);
   r_circle = Data.OnGet(59, 3, 1);
   n_circle = Data.OnGet(59, 4, 1);
   SpaceCharge = Data.OnGetB(59, 6, false);

   F.Smoos = Data.OnGet(61, 1, 1);
   spread_temp = Data.OnGetI(61, 2, 1);
   F.Theta_Eff = Data.OnGet(61, 4, 1);
   F.TempEff = Data.OnGet(61, 5, 1);
   /*
   F.undulator = Data.OnGetI(74,1,0);
   F.lambda    = Data.OnGet(74,2,10);
   F.B_field   = Data.OnGet(74,3,50);
   */
   // 05.06
   //Adist    = Data.OnGet(54,10,1);

   OnSet();

   return 0;
}

//---------------------------------------------------------------------------

vectorU xEbeam::Ibeam2Ebeam(xTime &t, vectorU ion)
{
   for (int i = 0; i < 6; i++)
      Shift[i] = 0;
   if (EnableShifts)
      ion = CoordinateShifts(t, ion);
   if (iDraw.InEcool.Enabled)
      iDraw.InEcool.Point(ion[0](), ion[5]());

   switch (model)
   {
   case 0:
      ion = UniCilinder(t, ion);
      break;
   case 1:
      ion = UniBunch(t, ion);
      break;
   case 2:
      ion = GaussCilinder(t, ion);
      break;
   case 3:
      ion = GaussBunch(t, ion);
      break;
   case 4:
      ion = HollowBeam(t, ion);
      break;
   case 5:
      ion = Array(t, ion);
      break;
   case 6:
      ion = Parabolic(t, ion);
      break;
   case 7:
      ion = DensityFile(t, ion);
      break;
   }
   return ion;
}
//---------------------------------------------------------------------------

doubleU xEbeam::UC_dp_P(doubleU r, doubleU Ie, doubleU Re) // calculates momentum shift for Uniformely distr. cilinder
{
   doubleU I0(U_me * (U_c ^ 3) / U_e);
   if (U_Abs(r) <= Re)
      return (1. - Neutralization) * Ie * r * r / (I0 * Re * Re * e_Energy.Gamma * e_Energy.Beta3);
   else
      return (1. - Neutralization) * Ie * (1 + 2 * U_Ln(U_Abs(r / Re))) / (I0 * e_Energy.Gamma * e_Energy.Beta3);
}
//---------------------------------------------------------------------------

doubleU xEbeam::UC_Vdrift(doubleU r, doubleU Ie, doubleU Re) // calculates drift angle (drift velocity over long. velocity) for Uniformely distr. cilinder
{
   if (U_Abs(r) <= Re)
      return (1. - (Neutralization * e_Energy.Gamma2)) *
             2. * U_Abs(r) * Ie / (U_c * F.mfield * Re * Re * e_Energy.Gamma2 * e_Energy.Beta2);
   else
      return (1. - (Neutralization * e_Energy.Gamma2)) *
             2. * Ie / (U_c * F.mfield * U_Abs(r) * e_Energy.Gamma2 * e_Energy.Beta2);
}
//---------------------------------------------------------------------------
/*
doubleU xEbeam::GC_dp_P(doubleU r) // calculates momentum shift for Gaussian distr. cilinder
{
  doubleU dP_P(U1_);
  doubleU I0(A_);
  doubleU lambda(m_^-1);

  I0 = U_me*(U_c^3)/U_e;
//   lambda = 0.01/ring.Circ;
//   shift = 1000.0*con2*lambda*(0.0-Form(r/(1.41*beam.a)));
//  dP_P += (1.-Neutralization) * bcurrent* r * r /(I0*bradius*bradius*e_Energy.Gamma*e_Energy.Beta3);

//shift *=con1*ring.ion.beta;
//}
//return shift;


  return dP_P;
}
//---------------------------------------------------------------------------

doubleU xEbeam::GC_Vdrift(doubleU r) // calculates drift angle (drift velocity over long. velocity) for Gaussian distr. cilinder
{
  doubleU Vd(U1_);
   Vd =(1.0 - Neutralization)*2.*U_Abs(r)*bcurrent/(U_c*F.mfield*bradius*bradius*e_Energy.Gamma2*e_Energy.Beta2);
  return Vd;
}
*/
//---------------------------------------------------------------------------

vectorU xEbeam::CoordinateShifts(xTime &t, vectorU &ion)
{
   doubleU s(m_), l(m_);

   if (FromFile)
   {
      for (int i = 1; i < Solenoid.Row(); i++)
      {
         if (Solenoid[i][0] >= t.so())
         {
            s = t.so() - Solenoid[i - 1][0];
            l = Solenoid[i][0] - Solenoid[i - 1][0];
            Initial[0] = Solenoid[i - 1][1];
            Initial[2] = Solenoid[i - 1][2];
            Final[0] = Solenoid[i][1];
            Final[2] = Solenoid[i][2];
            break;
         }
      }
      if (l() <= 0)
      {
         s = 1;
         l = 1;
         Initial[0] = 0;
         Initial[2] = 0;
         Final[0] = 0;
         Final[2] = 0;
      }
   }
   else
   {
      s = t.so;
      l = CoLength;
   }

   Shift[1] = (Final[0] - Initial[0]) / l;
   Shift[3] = (Final[2] - Initial[2]) / l;
   Shift[0] = Initial[0] + Shift[1] * s;
   Shift[2] = Initial[2] + Shift[3] * s;
   Shift[4] = Initial[4];
   Shift[5] = Initial[5];

   if (PaintingTable)
   {
      PaintingPeriod = PaintingData.Row();
      int period;
      if (ScallingFactor < 0)
         period = (t.Steps / abs(ScallingFactor)) % PaintingPeriod;
      else
         period = (t.Steps * ScallingFactor) % PaintingPeriod;
      for (int i = 0; i < 6; i++)
         Shift[i] = PaintingData[period][i];
   }
   else if (UsePeriod && PaintingPeriod)
   {
      double period = t.Steps % PaintingPeriod;
      period /= PaintingPeriod;

      for (int i = 0; i < 4; i++)
         Shift[i] *= period;

      if (EnableFinal)
      {
         Shift[4] += (Final[4] - Initial[4]) * period;
         Shift[5] += (Final[5] - Initial[5]) * period;
      }
      else
         for (int i = 4; i < 6; i++)
            Shift[i] *= period;
   }

   return ion - Shift;
}

//---------------------------------------------------------------------------

vectorU xEbeam::UniCilinder(xTime &t, vectorU &ion) // models electron beam like a uniformly charged cilinder
{
   vectorU ion_new;
   ion_new = ion;
   doubleU rc = ((ion_new[0] * ion_new[0]) + (ion_new[2] * ion_new[2])) ^ 0.5;

   if (rc > bradius)
   {
      inside = false;
      return ion_new;
   }

   F.n_e = bcurrent / (U_pi * e_Energy.Gamma * bradius * bradius * e_ * e_Energy.Velocity);
   F.Ttemp = F.Ttemp_centre;

   F.V_tr_e = (F.Ttemp / U_me) ^ 0.5;
   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;

   if (Neutralization < 2.)
   {
      doubleU Vdr(U1_);
      doubleU E_e(M_6 * eV_ * 0.5110034);
      ion_new[5] -= UC_dp_P(rc, bcurrent, bradius);
      Vdr = UC_Vdrift(rc, bcurrent, bradius);

      F.Ttemp += Vdr * Vdr * E_e / 2.;
      F.V_tr_e = (F.Ttemp / U_me + (dVdr * rc * dVdr * rc)) ^ 0.5;

      if (rc() > 0.)
      {
         ion_new[1] += Vdr * ion_new[2] / rc;
         ion_new[3] += Vdr * ion_new[0] / rc;
      }
   }

   inside = true;
   return ion_new;
}
//---------------------------------------------------------------------------

vectorU xEbeam::Parabolic(xTime &t, vectorU &ion) // models electron beam like a uniformly charged cilinder
{
   vectorU ion_new;
   ion_new = ion;
   doubleU rc = ((ion_new[0] * ion_new[0]) + (ion_new[2] * ion_new[2])) ^ 0.5;

   if (rc > pradius)
   {
      inside = false;
      return ion_new;
   }

   F.n_e = n_centre * (1.0 - (rc * rc / (pradius * pradius)));
   F.T_plasma = ((U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);
   F.Ttemp = F.Ttemp_centre;

   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;
   F.V_tr_e = (F.Ttemp / U_me + (dVdrp * rc * dVdrp * rc)) ^ 0.5;

   inside = true;
   return ion_new;
}
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07

vectorU xEbeam::DensityFile(xTime &t, vectorU &ion) //density distribution from file
{
   vectorU ion_new;
   ion_new = ion;
   doubleU rc = ((ion_new[0] * ion_new[0]) + (ion_new[2] * ion_new[2])) ^ 0.5;

   if (rc > dradius)
   {
      inside = false;
      return ion_new;
   }

   int i = 0;
   doubleU r(m_3 * m_);
   for (i = 0; i < Density.Row() - 1; i++)
   {
      r = Density[i][0];
      if (r >= rc)
         break;
   }
   doubleU J(Density[i][1], A_ / (cm_ ^ 2));
   F.n_e = J / (e_Energy.Gamma * e_ * e_Energy.Velocity);
   F.T_plasma = ((U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);
   F.Ttemp = F.Ttemp_centre;

   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;

   doubleU dV = dVdrd;
   if (Density.Col(i) == 3)
      dV = Density[i][2];
   F.V_tr_e = (F.Ttemp / U_me + (dV * rc * dV * rc)) ^ 0.5;

   inside = true;
   return ion_new;
}
//---------------------06.07---------------------------------------------------

vectorU xEbeam::HollowBeam(xTime &t, vectorU &ion)
{
   vectorU ion_new;
   ion_new = ion;
   doubleU rc = ((ion_new[0] * ion_new[0]) + (ion_new[2] * ion_new[2])) ^ 0.5;

   if (rc > r_circle)
   {
      inside = false;
      return ion_new;
   }
   if (rc > r_hole)
      F.n_e = n_circle / e_Energy.Gamma;
   else
      F.n_e = n_hole / e_Energy.Gamma;
   F.Ttemp = F.Ttemp_centre;

   F.V_tr_e = (F.Ttemp / U_me) ^ 0.5;
   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;

   if (SpaceCharge)
   {
      doubleU Vdr(U1_);
      doubleU E_e(M_6 * eV_ * 0.5110034);

      ion_new[5] -= UC_dp_P(rc, I_circle, r_circle) - UC_dp_P(rc, I_hocle, r_hole) + UC_dp_P(rc, I_hole, r_hole);

      Vdr = UC_Vdrift(rc, I_circle, r_circle) - UC_Vdrift(rc, I_hocle, r_hole) + UC_Vdrift(rc, I_hole, r_hole);

      F.Ttemp += Vdr * Vdr * E_e / 2.;
      F.V_tr_e = (F.Ttemp / U_me) ^ 0.5;

      if (rc() > 0.)
      {
         ion_new[1] += Vdr * ion_new[2] / rc;
         ion_new[3] += Vdr * ion_new[0] / rc;
      }
   }
   inside = true;
   return ion_new;
}
//---------------------------------------------------------------------------

vectorU xEbeam::GaussCilinder(xTime &t, vectorU &ion) // models electron beam like a Gaussian charged cilinder
{
   vectorU ion_new;
   ion_new = ion;
   if ((ion_new[0] > 3. * sigma_x_cil) || (ion_new[2] > 3. * sigma_y_cil))
   {
      inside = false;
      return ion_new;
   }

   F.V_tr_e = (F.Ttemp / U_me) ^ 0.5;
   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;

   doubleU E_e(M_6 * eV_ * 0.5110034);
   F.n_e = Ne_cil / (((2 * U_pi) ^ 1.5) * sigma_x_cil * sigma_y_cil * e_Energy.Gamma);
   F.n_e *= 1. / U_Exp(((ion_new[0] / sigma_x_cil) ^ 2.) / 2. + ((ion_new[2] / sigma_y_cil) ^ 2.) / 2.);
   F.T_plasma = ((U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);

   inside = true;
   return ion_new;
}
//---------------------------------------------------------------------------
// ----- 05.06 -- dist between bunches
vectorU xEbeam::UniBunch(xTime &t, vectorU &ion)
{
   vectorU ion_new;
   ion_new = ion;
   //ion_new[4] += dist_uni;
   if ((U_Abs(ion_new[0]) > size_x) || (U_Abs(ion_new[2]) > size_y) ||
       (U_Abs(ion_new[4]) > 3. * size_s))
   {
      inside = false;
      return ion_new;
   }

   F.V_tr_e = (F.Ttemp / U_me) ^ 0.5;
   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;

   doubleU E_e(M_6 * eV_ * 0.5110034);
   F.n_e = Ne_uni / (((2 * U_pi) ^ 0.5) * U_pi * size_x * size_y * size_s * e_Energy.Gamma);
   F.n_e /= U_Exp(((ion_new[4] / size_s) ^ 2.) / 2.);
   F.T_plasma = ((U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);

   inside = true;
   return ion_new;
}
//---------------------------------------------------------------------------
// ----- 06.07 -- sigmas from form
vectorU xEbeam::GaussBunch(xTime &t, vectorU &ion) // models electron beam like a bunch with Gaussian distribution
{
   vectorU ion_new;
   ion_new = ion;
   //ion_new[4] += dist;
   if ((U_Abs(ion_new[0]) > 3. * Gsigma[0]) || (U_Abs(ion_new[2]) > 3. * Gsigma[2]) ||
       (U_Abs(ion_new[4]) > 3. * Gsigma[4]))
   {
      inside = false;
      return ion_new;
   }

   F.V_tr_e = (F.Ttemp / U_me) ^ 0.5;
   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;
   F.V_tr_x = e_Energy.Gamma * e_Energy.Beta * U_c * Gsigma[1];
   F.V_tr_y = e_Energy.Gamma * e_Energy.Beta * U_c * Gsigma[3];

   doubleU E_e(M_6 * eV_ * 0.5110034);
   F.n_e = Ne / (((2 * U_pi) ^ 1.5) * Gsigma[0] * Gsigma[2] * Gsigma[4] * e_Energy.Gamma);
   F.n_e /= U_Exp((((ion_new[0] / Gsigma[0]) ^ 2.) / 2.) +
                  (((ion_new[2] / Gsigma[2]) ^ 2.) / 2.) + (((ion_new[4] / Gsigma[4]) ^ 2.) / 2.));
   F.T_plasma = ((U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);

   inside = true;
   return ion_new;
}
//---------------------------------------------------------------------------

vectorU xEbeam::Array(xTime &t, vectorU &ion) // models electron beam like a bunch with Gaussian distribution
{
   doubleU M_e(U_me * U_c * U_c);
   doubleU lengthU(1, m_);

   vectorU ion_new, ion_local;
   ion_new = ion;
   doubleU n_loc_an(m_ ^ -3);
   doubleU diference(U1_);

   int Loc_n;

   if ((U_Abs(ion_new[0]) > 3. * Asigma[0]) || (U_Abs(ion_new[2]) > 3. * Asigma[2]) ||
       (U_Abs(ion_new[4]) > 3. * Asigma[4]))
   {
      inside = false;
      return ion_new; // ion out of ebunch in LRF
   }

   ion_local = ion_new; // ion in PRF
   ion_local[1] *= U_c * e_Energy.Beta * e_Energy.Gamma;
   ion_local[3] *= U_c * e_Energy.Beta * e_Energy.Gamma;
   ion_local[4] *= e_Energy.Gamma;
   ion_local[5] *= U_c * e_Energy.Beta;
   GetMean(ion_local, Local); // get Mean in PRF

   F.Ttemp = M_e * 0.5 * (Mean2[1] + Mean2[3]) / (U_c() * U_c()); // F in PRF
   F.Ltemp = M_e * Mean2[5] / (U_c() * U_c());
   F.V_tr_e = (F.Ttemp / U_me) ^ 0.5;
   F.V_tr_x = U_c * sqrt(Mean2[1]) / (U_c());
   F.V_tr_y = U_c * sqrt(Mean2[3]) / (U_c());
   F.V_long_e = (F.Ltemp / U_me) ^ 0.5;
   F.V_eff_e = (F.TempEff / U_me) ^ 0.5;

   Loc_n = 0;
   for (int i = 0; i < Minst.Number; i++)
      if (
          (
              ((Local[Minst[i]][4] - ion_local[4]()) * (Local[Minst[i]][4] - ion_local[4]())) < (1.0 * Mean2[4])) &&
          ((
               ((Local[Minst[i]][0] - ion_local[0]()) * (Local[Minst[i]][0] - ion_local[0]()) / (Mean2[0] * Box * Box)) +
               ((Local[Minst[i]][2] - ion_local[2]()) * (Local[Minst[i]][2] - ion_local[2]()) / (Mean2[2] * Box * Box))) < 1.0))
         Loc_n++;

   F.n_e = Loc_n * ANe / (2.0 * M_PI * Box * Box * sqrt(1.0 * Mean2[0] * Mean2[2] * Mean2[4]) * Bunch.Number * lengthU * lengthU * lengthU);

   //2D + Gauss
   /*
   F.n_e = Minst.Number*ANe /(2*sqrt(2.)*2.*M_PI*sqrt(2.*M_PI* Mean2[0]* Mean2[2]* Mean2[4])*
           Bunch.Number*lengthU*lengthU*lengthU);
   */
   /*
   //3D + Gauss corrected

   F.n_e = Minst.Number*ANe /(4.*M_PI*sqrt(2.*M_PI*Mean2[0]*Mean2[2]*Mean2[4])*
           Bunch.Number*lengthU*lengthU*lengthU);
   */

   //ellipsoid
   /*
   F.n_e = 2.0*Minst.Number*ANe /(5.*sqrt(5.)*(4./3.)*M_PI*sqrt(Mean2[0]*Mean2[2]*Mean2[4])*
           Bunch.Number*lengthU*lengthU*lengthU);

   */
   /*
   // box from max
   if ((fabs(ion_local[0]()-Mean[0]) > sqrt(Mean2[0]))||(fabs(ion_local[2]()-Mean[2]) > sqrt(Mean2[2]))||
       (fabs(ion_local[4]()-Mean[4]) > sqrt(Mean2[4])))
   {  inside = false;
      return ion_new;                                      // ion out of ebunch in LRF
   }

   F.n_e = Minst.Number*ANe /((4./3.)*M_PI*sqrt(Mean2[0]*Mean2[2]*Mean2[4])*
           Bunch.Number*lengthU*lengthU*lengthU);
   */
   // analytical from the ion

   n_loc_an = ANe * U_Exp(((-1.0 * ion_local[0] * ion_local[0]) / (2.0 * Asigma[0] * Asigma[0])) + ((-1.0 * ion_local[2] * ion_local[2]) / (2.0 * Asigma[2] * Asigma[2])) + ((-1.0 * ion_local[4] * ion_local[4]) / (2.0 * Asigma[4] * Asigma[4] * e_Energy.Gamma * e_Energy.Gamma))) /
              (2.0 * M_PI * sqrt(2.0 * M_PI) * Asigma[0] * Asigma[2] * Asigma[4] * e_Energy.Gamma);
   if (BiGauss)
   {
      n_loc_an *= (Bunch.Number - N_core) / double(Bunch.Number);
      n_loc_an += (N_core * ANe / Bunch.Number) * U_Exp(((-1.0 * ion_local[0] * ion_local[0]) / (2.0 * sigma_core * sigma_core * Asigma[0] * Asigma[0])) + ((-1.0 * ion_local[2] * ion_local[2]) / (2.0 * sigma_core * sigma_core * Asigma[2] * Asigma[2])) + ((-1.0 * ion_local[4] * ion_local[4]) / (2.0 * sigma_core * sigma_core * Asigma[4] * Asigma[4] * e_Energy.Gamma * e_Energy.Gamma))) /
                  (2.0 * M_PI * sqrt(2.0 * M_PI) * sigma_core * sigma_core * sigma_core * Asigma[0] * Asigma[2] * Asigma[4] * e_Energy.Gamma);
   }
   if (Uniform)
   {
      n_loc_an = 0.0;
      if ((U_Abs(ion_new[0]) < 2. * Asigma[0]) && (U_Abs(ion_new[2]) < 2. * Asigma[2]))
         n_loc_an = ANe * U_Exp((-1.0 * ion_local[4] * ion_local[4]) / (2.0 * Asigma[4] * Asigma[4] * e_Energy.Gamma * e_Energy.Gamma)) /
                    (4.0 * M_PI * sqrt(2.0 * M_PI) * Asigma[0] * Asigma[2] * Asigma[4] * e_Energy.Gamma);
   }
   //
   // analytical from the array
   /*
   F.n_e = ANe* U_Exp(((-1.0* Mean[0] *Mean[0])/(2.0*Asigma[0]*Asigma[0])) +
           ((-1.0* Mean[2] *Mean[2])/(2.0*Asigma[2]*Asigma[2])) +
           ((-1.0* Mean[4] *Mean[4])/(2.0*Asigma[4]*Asigma[4]*e_Energy.Gamma*e_Energy.Gamma)) )/
           (2.0*M_PI*sqrt(2.0*M_PI)*Asigma[0]*Asigma[2]*Asigma[4]*e_Energy.Gamma);
   */
   //
   //diference = F.n_e/n_loc_an;

   if (GaussN)
      F.n_e = n_loc_an;
   if (F.n_e() == 0.0)
   {
      inside = false;
      return ion_new; // ion out of ebunch in LRF
   }
   F.T_plasma = ((U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);

   Shift[1] += Mean[1] / (e_Energy.Beta * e_Energy.Gamma * U_c)(m_ / s_); // mean velocity in LRF
   Shift[3] += Mean[3] / (e_Energy.Beta * e_Energy.Gamma * U_c)(m_ / s_);
   Shift[5] += Mean[5] / (e_Energy.Beta * U_c)(m_ / s_);

   inside = true;
   return ion_new;
}
//---------------------------------------------------------------------------
/*
vectorU xEbeam::Array2(xTime&t, vectorU&ion) // models electron beam like a bunch with Gaussian distribution
{
   doubleU M_e(U_me*U_c*U_c);
   vectorU ion_new;
   ion_new = ion;
   ion_new[4] += Adist;

   if ((U_Abs(ion_new[0]) > 3.*Asigma[0])||(U_Abs(ion_new[2]) > 3.*Asigma[2])||
       (U_Abs(ion_new[4]) > 3.*Asigma[4]))
   {  inside = false;
      return ion_new;
   }

   doubleU lengthU(1,m_);
   double x2, y2, s2;
   double max, min;
   int test12;
   double test22;
   doubleU test33(m_);

   double mean_r2;
   mean_r2 = 0;
   double mean_x, mean_y, mean_s;
   mean_x = 0;
   mean_y = 0;
   mean_s = 0;
   double mean_x_s, mean_y_s, mean_p;
   mean_x_s = 0;
   mean_y_s = 0;
   mean_p = 0;
   double mean_x2, mean_y2, mean_s2;
   mean_x2 = 0;
   mean_y2 = 0;
   mean_s2 = 0;
   double mean_x_s2, mean_y_s2, mean_p2;
   mean_x_s2 = 0;
   mean_y_s2 = 0;
   mean_p2 = 0;

   for (int i0 = 0; i0 < Bunch.Number; i0++)
   {
      x2  = (Bunch[i0][0] - ion_new[0]())*(Bunch[i0][0] - ion_new[0]());
      y2  = (Bunch[i0][2] - ion_new[2]())*(Bunch[i0][2] - ion_new[2]());
      s2  = (Bunch[i0][4] - ion_new[4]())*(Bunch[i0][4] - ion_new[4]());

      Dist2[i0] = x2 + y2 + s2;
   }

   for(int j = 0; j < Minst.Number; j++)
   {  max = 1e100;                                      // start from bigest value
      for(int i = 0; i < Dist2.Number; i++)
      {  if((max > Dist2[i]) && (Dist2[i] > min))       // search for min exept last
         {  max = Dist2[i];                             // keep minimum distance
            Minst[j] = i;                               // keep index
         }
      }
      min = max;                                        // keep last value
   }

   for (int i0 = 0; i0 < Minst.Number; i0++)
     {
       mean_r2  += Dist2[Minst[i0]];
       mean_x   += Bunch[Minst[i0]][0];
       mean_y   += Bunch[Minst[i0]][2];
       mean_s   += Bunch[Minst[i0]][4];
       mean_x_s += Bunch[Minst[i0]][1];
       mean_y_s += Bunch[Minst[i0]][3];
       mean_p   += Bunch[Minst[i0]][5];
     }
   mean_r2 /= Minst.Number;
   mean_x /= Minst.Number;
   mean_y /= Minst.Number;
   mean_s /= Minst.Number;
   mean_x_s /= Minst.Number;
   mean_y_s /= Minst.Number;
   mean_p /= Minst.Number;
   m_x = mean_x_s;
   m_y = mean_y_s;
   m_p = mean_p;
   //if(Mean)
   {
      ion_new[1] += mean_x_s;
      ion_new[3] += mean_y_s;
      ion_new[5] += mean_p;
   }

   for (int i0 = 0; i0 < Minst.Number; i0++)
     {
       mean_x_s2 += (Bunch[Minst[i0]][1]- mean_x_s)*(Bunch[Minst[i0]][1]- mean_x_s);
       mean_y_s2 += (Bunch[Minst[i0]][3]- mean_y_s)*(Bunch[Minst[i0]][3]- mean_y_s);
       mean_p2 += (Bunch[Minst[i0]][5]- mean_p)*(Bunch[Minst[i0]][5]- mean_p);
       mean_x2 += (Bunch[Minst[i0]][0] - mean_x)*(Bunch[Minst[i0]][0] - mean_x);
       mean_y2 += (Bunch[Minst[i0]][2] - mean_y)*(Bunch[Minst[i0]][2] - mean_y);
       mean_s2 += (Bunch[Minst[i0]][4] - mean_s)*(Bunch[Minst[i0]][4] - mean_s);
     }

   mean_x_s2 /= Minst.Number;
   mean_y_s2 /= Minst.Number;
   mean_p2 /= Minst.Number;
   mean_x2 /= Minst.Number;
   mean_y2 /= Minst.Number;
   mean_s2 /= Minst.Number;

   //test22 = pow(mean_x2+mean_y2+mean_s2,1.5);
   //test22 = pow(mean_r2,1.5);
   F.Ttemp = M_e*e_Energy.Beta* e_Energy.Beta*e_Energy.Gamma*
                e_Energy.Gamma*0.5*(mean_x_s2 + mean_y_s2);
   F.Ltemp = M_e * mean_p2 *e_Energy.Beta* e_Energy.Beta;


   F.V_tr_e   = (F.Ttemp  /U_me)^0.5;
   F.V_tr_x   = U_c*e_Energy.Beta*e_Energy.Gamma*sqrt(mean_x_s2);
   F.V_tr_y   = U_c*e_Energy.Beta*e_Energy.Gamma*sqrt(mean_y_s2);
   F.V_long_e = (F.Ltemp  /U_me)^0.5;
   F.V_eff_e  = (F.TempEff/U_me)^0.5;

   F.n_e = Minst.Number*ANe /(2.3*sqrt(2.0)*2.0*M_PI*sqrt(2.0*M_PI*mean_x2*mean_y2*mean_s2)*
           Bunch.Number *e_Energy.Gamma *lengthU*lengthU*lengthU);

   F.T_plasma = ((U_me/(4*U_pi*F.n_e*U_e*U_e))^0.5);

   inside = true;
   return ion_new;
}
*/
//---------------------------------------------------------------------------

double xEbeam::FormFactor(double t) // calculates field distribution in radial direction, returns dimensionless result
                                    // this function calculates integral by rectangulars method and is used for
                                    // dP/P shift calculation for gaussian distribution in cylinder or bunch
{
   double F = 1;
   int i = 100;
   double step = fabs(t / (double)i);
   double temp;
   int j;

   if (t == 0.)
      return F;
   else
   {
      for (j = 1; j < i; j++)
      {
         temp = step * j;
         F += step * (1.0 - exp(-temp * temp)) / temp;
      }
      return F;
   }
}
//---------------------------------------------------------------------------

vectorU xEbeam::Ebeam2Ibeam(vectorU &ion)
{
   vectorU ion_new(m_, U1_);
   ion_new = ion;

   return ion_new;
}
//---------------------------------------------------------------------------
