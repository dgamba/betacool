//---------------------------------------------------------------------------
#include "stdafx.h"               
#include "xbeam.h"
#include "xdynamic.h"
#include <iostream>
//using namespace std;
//---------------------------------------------------------------------------

xBeam::xBeam(xRing& ring, int n)
{
   pRing = &ring;
   pEnergy = &ring.Energy;
   pCirc   = &ring.Circ;
   pTrev   = &ring.Trev;
   pBetaH  = &ring.BetaH;
   pBetaV  = &ring.BetaV;
   pEta    = &ring.Eta;
   //20.11.06
   pBB_V_B = &ring.V_B;                             //Voltage amplitude
   pBB_T1 = &ring.T1;                               //RF duration in T0
   pBB_T2 = &ring.T2;                               //Gap duration in T0
   pBB_BarHeight = &ring.BarHeight;                       //Barrier height
   //
   Emit[0]._(m_);
   Emit[1]._(m_);
   Emit[2]._(U1_);
   Emit[3]._(U1_);
   Ep._(eV_ * s_);
   Emit[5]._(U1_);
   CellSize._(m_);
   Impact._(m_);

   percent = 38;
   Nrms = 1;
   iMomentum = 0;
   uMomentum(eV_c_);

   Number(n);
   BParam[0] = &rr;
   BParam[1] = &ss;
   BParam[2] = &Ir;
   BParam[3] = &Is;
   BParam[4] = &L;
   Index = -1;
   Ttrans._(K_);
   Tlong._(K_);

//*********************   from Grisha ****************************
   pQ_x = &ring.TunesH;
   pKsy_x = &ring.HromatH;
   pR = &ring.R_m;


   benum = COASTING;  // coasting regime

   F_sc._(U1_);       // Image force correction factor
   F_l._(U1_);        // Longitudinal distribution factor
   Z_l_i._(Ohm_);
   Z_l._(Ohm_);       // Longitudinal coupling impedance
   F_t._(U1_);        // Transverse distribution factor
   Z_t_i._(Ohm_/m_);
   Z_t._(Ohm_/m_);    // Transverse coupling impedance
   G._(U1_);          // Longitudinal form-factor
   B_f._(U1_);        // Bunching factor
   a._(m_);           // Mean Beam radius
   Z_l_sc._(Ohm_);    // Longitudinal space charge impedance
   Z_t_sc._(Ohm_/m_); // Transverse space charge impedance
   I._(A_);           // Beam current, A
   D_Q_x._(U1_);      // Tune shift
   D_Q_z._(U1_);      // Tune shift
   KS._(U1_);         // Keil-Schnell parameter
   S_Q._(U1_);        // Tune spread
   DM._(U1_);         // Dipole mode parameter
   s_s._(cm_);        // Rms bunch length, cm
   N_max._(U1_);      // Maximum number of particles
   Q_s._(U1_);        // Synchrotron  tune
   N_b._(U1_);        // Number of bunches

// declaration for re-named pointers (from Ring and Energy):
//***************************************
   tz_0._(377, Ohm_);

//***************************************
   pRF_h = &ring.h;
   pRF_V = &ring.V;
   pRF_Q_s = &ring.Q_s;
   pRF_L_s = &ring.L_s;

   pz_0 = &tz_0;
   pVac_a = &ring.a_m;
//****************
//   N_b = 20;
//   Huge = 1;
//****************
   //collider = false;

   inject = false;
   interval._(s_);
   cycles = 1;

   T_s._(s_);                             //Synchrotron period
   Size_Bac._(m_);                        //Bucket length
   t_1._(s_);
   t_2._(s_);
}

int xBeam1::OnGet()
{
   Emit[0] = Data.OnGet(1,1,1)*1e-6;
   Emit[1] = Data.OnGet(1,2,1)*1e-6;
   Emit[2] = pow(Data.OnGet(1,3,1e-3), 2);
   Emit[3] = Data.OnGet(1,4,1e3);
   InitialEmit = Emit;

   IniNumber = Data.OnGetI(2,1,1000);
   Number(IniNumber);
   CellNum = Emit[3] / Number();
   ModelParticles = Number();
   MaxNumber = Data.OnGetI(2,2,0);
   if (MaxNumber < 10) MaxNumber = 0;
   Impact = Data.OnGet(2,3,1e-6);

   benum = (BunchEnum)((int)Data[1][5]);
   EmitDef =   Data.OnGetI(1,12,0);
   UseForIBS = Data.OnGetB(1,13,true);
   percent = Data.OnGet(1,17,40);
   iMomentum = Data.OnGetI(1,22,0);
   switch (iMomentum)
   {case 0: uMomentum(eV_c_);      break;
    case 1: uMomentum(k_3*eV_c_);  break;
    case 2: uMomentum(M_6*eV_c_);  break;
    case 3: uMomentum(G_9*eV_c_);  break;
    case 4: uMomentum(eV_);        break;
    case 5: uMomentum(k_3*eV_);    break;
    case 6: uMomentum(M_6*eV_);    break;
    case 7: uMomentum(G_9*eV_);    break;
   }

   Sigmas    = Data.OnGet(30,2,5);
   Division  = Data.OnGetI(30,3,50);
   percentage= Data.OnGet(1,14,50);
   percentlong = Data.OnGet(1,15,50);
   MeanPercents = Data.OnGetB(1,16,false);

   //collider = Data.OnGetI(1,6,0);
   N_b = Data.OnGet(4,1,1);

   F_sc  = Data.OnGet(3,1,1);
   F_l   = Data.OnGet(3,4,1);
   Z_l_i = Data.OnGet(3,5,1);
   F_t   = Data.OnGet(3,8,1);
   Z_t_i = Data.OnGet(3,9,1);

   GenerateOn = Data.OnGetI(9,9,0);
   HourGlass  = Data.OnGetB(10,5,false);

   Hyst3.SetNumber(4, Division+1);
   CurrentAv.SetNumber(Division+1);

   initial  = Data.OnGetI(14,1,0);
   inject   = Data.OnGetB(14,2,false);
   interval = Data.OnGet(14,3,1);
   cycles   = Data.OnGetI(14,4,1);
   dpshift  = Data.OnGet(14,5,0);
   horshift  = Data.OnGet(14,6,0);
   vershift  = Data.OnGet(14,7,0);
   longshift  = Data.OnGet(14,11,0);

   LProfile.SetNumber(Data.OnGetI(26,14,100));
   VProfile.SetNumber(Data.OnGetI(26,14,100));
   Integral.SetNumber(Data.OnGetI(26,14,100));
   Average. SetNumber(Data.OnGetI(26,14,100));
   Momentum.SetNumber(Data.OnGetI(26,14,100));
   RatesDp .SetNumber(Data.OnGetI(26,14,100));

   return 0;
}

int xBeam1::OnRun()
{
   Sigmas     = Data.OnGet(30,2,5);
   Division   = Data.OnGetI(30,3,50);

   GenerateOn = Data.OnGetI(9,9,0);
   HourGlass  = Data.OnGetB(10,5,false);

   Hyst3.SetNumber(4, Division+1);
   CurrentAv.SetNumber(Division+1);

   //initial  = Data.OnGetI(14,1,0);
   inject   = Data.OnGetB(14,2,false);
   interval = Data.OnGet(14,3,1);
   cycles   = Data.OnGetI(14,4,1);
   dpshift  = Data.OnGet(14,5,0);
   horshift  = Data.OnGet(14,6,0);
   vershift  = Data.OnGet(14,7,0);
   longshift  = Data.OnGet(14,11,0);

   LProfile.SetNumber(Data.OnGetI(26,14,100));
   VProfile.SetNumber(Data.OnGetI(26,14,100));
   Integral.SetNumber(Data.OnGetI(26,14,100));
   Average. SetNumber(Data.OnGetI(26,14,100));
   Momentum.SetNumber(Data.OnGetI(26,14,100));
   RatesDp .SetNumber(Data.OnGetI(26,14,100));

   return 0;
}

int xBeam1::OnSet()
{
   for (int i = 0; i < 3; i++)
   {
      Fitting[i][0] = 0;
      Fitting[i][1] = 1;
      Fitting[i][2] = sqrt(2*M_PI) * Number() * Sigmas / Division;
      if (iIBS.IBSkickModel)
         Fitting[i][2] = 1;
      Fitting[i][3] = sqrt(2*M_PI) * Number() * Sigmas / Division;
      Fitting[i][4] = sqrt(2*M_PI) * Number() * Sigmas / Division;
   }

   if (initial == 2)
   {
      Data.OnSet(1,1,Emit[0](u_6*m_));
      Data.OnSet(1,2,Emit[1](u_6*m_));
      Data.OnSet(1,3,sqrt(Emit[2]()));
      Data.OnSet(2,1,Number());
   }
   if (CalcStability())
   {
//   RMS(x_Lattice);
      Keil_Shneile();
      Tune_shift();
      Shneil_Zotter();
      //Luminosity();
   }
   //20.11.06
   if (benum == BUNCHED) CalcBunch();
   if (benum == BUCKET)  CalcBucket(Emit[2](U1_));
   //
// for General Beam Parameters:
   Data.OnSet(1,7,a);
   Data.OnSet(1,8,G);
   Data.OnSet(1,9,Z_l_sc);
   Data.OnSet(1,10,Z_t_sc);
   Data.OnSet(1,11,I);

   if (percent > 0 && percent < 100)
   {
      Nrms = - 2. * log(1. - percent/100.);
      Data.OnSet(1,18,Nrms);
      Data.OnSet(1,19,(Emit[0]*pEnergy->Beta*pEnergy->Gamma*Nrms)(u_6*m_));
      Data.OnSet(1,20,(Emit[1]*pEnergy->Beta*pEnergy->Gamma*Nrms)(u_6*m_));
   }
   if (iMomentum < 4)
      Data.OnSet(1,21,((Emit[2]^0.5) * pEnergy->Momentum)(uMomentum));
   else
      Data.OnSet(1,21,((Emit[2]^0.5) * pEnergy->Kinetic *
                   (pEnergy->Gamma + 1) / pEnergy->Gamma)(uMomentum));



// for Beam Stability characteristics:
   Data.OnSet(3,2,D_Q_x);
   Data.OnSet(3,3,D_Q_z);
   Data.OnSet(3,6,KS);
   Data.OnSet(3,7,S_Q);
   Data.OnSet(3,10,DM);
// for Bunched Beam parameters:
   Data.OnSet(4,2,s_s);
   Data.OnSet(4,3,N_max);
   Data.OnSet(4,4,Q_s);
   Data.OnSet(4,5,B_f);
   Data.OnSet(4,6,Ep);
// for Luminosity:
   //Data[1+6][10] = Lum.L()/10000.;

   CellSize = *pCirc;
   s_s = CellSize;
   Data.OnSet(2,4,*pCirc / CellNum);
   Lambda = ( Emit[3] / *pCirc *
   (( (pEnergy->Z^2) * 1.5 * U_rp * 0.5 *
     ((*pBetaH^2)+(*pBetaV^2)) /
     (pEnergy->A * (pEnergy->Gamma^5) * (pEnergy->Beta^2))
    )^(1./3.))    )(U1_);
   Data.OnSet(2,5,Lambda);

   Temperature();
   Data.OnSet(2,6,Tlong(K_));
   Data.OnSet(2,7,Ttrans(K_));
   Data.OnSet(2,8,Gamma1);
   Data.OnSet(2,9,Gamma2);
   //20.11.06
   Data.OnSet(26,6,Size_Bac);
   Data.OnSet(26,5,Amp2^0.5);
   Data.OnSet(26,7,T_s);

   return 0;
}
//---------------------------------------------------------------------------
void Bessel(double**);
void Bessel(double r, double s, double &dudr, double &duds);
/*
void xBeam::SpaceChargeND(xTime& time)
{
   static double dx, dy, ds, dr, a3, F[3], Kr, impact;
   L = CellSize()*pEnergy->Gamma();
   if (Impact() > 0)
      impact = pow(Impact(m_), -3.);
   else
      impact = 1e300;

   static doubleU m_2(1, m_^-2);
   doubleU K(((pEnergy->Z*U_e)^2) * time.ds * m_2 /
          (pEnergy->Momentum*pEnergy->Velocity*pEnergy->Gamma));
   double K1 = K(U1_);

   for (int i = 0; i < Bunch.Number-1; i++)
   {	for (int j = i+1; j < Bunch.Number; j++)
      {  dx = Bunch[i][0] - Bunch[j][0];
         dy = Bunch[i][2] - Bunch[j][2];
         ds = (Bunch[i][4] - Bunch[j][4]) * pEnergy->Gamma();
         ds-= L * floor(ds / L + 0.5);
         dr = dx*dx + dy*dy;
         a3 = pow(ds*ds + dr, -1.5);
         if (a3 < impact)
         {  dr = sqrt(dr);
            rr = dr / L;
            ss = ds / L;
            Bessel(BParam);
            if (dr)
               Kr = a3 - Ir * 2 / (dr * L * L);
            else
               Kr = 0;

            F[0] = dx * Kr;
            F[1] = dy * Kr;
            F[2] = ds * a3 - Is * 2 / L;

            for (int k = 0; k < 3; k++)
            {  Bunch[i][k*2+1] += K1 * F[k];
               Bunch[j][k*2+1] -= K1 * F[k];
            }
         }else
            Warning("Particle distance less than Impact Parameter : ", pow(a3, -1./3.),"[ m ]");
         }
      }
}
*/
void xBeam::SpaceChargeMD(xTime& time)
{
   double L, impact;
   L = CellSize()*pEnergy->Gamma();
   if (Impact() > 0)
      impact = pow(Impact(m_), -3.);
   else
      impact = 1e300;

   static doubleU m_2(1, m_^-2);
   doubleU K(((pEnergy->Z*U_e)^2) * time.ds * m_2 /
          (pEnergy->Momentum*pEnergy->Velocity*pEnergy->Gamma));
   double K1 = K(U1_);

   int nimpact = 0;
   #pragma omp parallel for reduction(+: nimpact)
   for (int i = 0; i < Bunch.Number-1; i++)
   {	for (int j = i+1; j < Bunch.Number; j++)
      {  
         double dx, dy, ds, dr, a3, F[3], Ir, Is, Kr;         
         dx = Bunch[i][0] - Bunch[j][0];
         dy = Bunch[i][2] - Bunch[j][2];
         ds = (Bunch[i][4] - Bunch[j][4]) * pEnergy->Gamma();
         ds-= L * floor(ds / L + 0.5);
         dr = dx*dx + dy*dy;
         a3 = pow(ds*ds + dr, -1.5);
         
         if (a3 > impact) { // 2016/04/18 by Alex
            a3 = impact;
            nimpact ++; 
         } 
         //if (a3 < impact)
         {  dr = sqrt(dr);
            Bessel(dr / L, ds / L, Ir, Is);
            if (dr)
               Kr = a3 - Ir * 2 / (dr * L * L);
            else
               Kr = 0;

            F[0] = dx * Kr;
            F[1] = dy * Kr;
            F[2] = ds * a3 - Is * 2 / (L * L);

            for (int k = 0; k < 3; k++)
            {  Bunch[i][k*2+1] += K1 * F[k];
               Bunch[j][k*2+1] -= K1 * F[k];
            }
         }
         //else Warning("Particle distance less than Impact Parameter : ", pow(a3, -1./3.),"[ m ]");
      }
   }
   NImpact += nimpact;
}

// 2-cell approximation of space charge (by Alexander Smirov)

void xBeam::SpaceCharge2C(xTime& time)
{
   static double a1, a2, a0, L;
   static double s1, s2, F[3];
   static doubleU m_2(1, m_^-2);

   L = CellSize()*pEnergy->Gamma();
   doubleU K(((pEnergy->Z*U_e)^2) * time.ds * m_2 /
        (pEnergy->Momentum*pEnergy->Velocity*pEnergy->Gamma));
   double K1 = K(U1_);

   for (int i = 0; i < Bunch.Number - 1; i++)
   {	for (int j = i+1; j < Bunch.Number; j++)
      {
         s1 = (Bunch[i][4] - Bunch[j][4]) * pEnergy->Gamma();
         s1-= L * floor(s1 / L + 0.5);
         s2 = s1 - L;
         a0 = pow(Bunch[i][0] - Bunch[j][0],2) + pow(Bunch[i][2] - Bunch[j][2],2);
         a1 = pow(pow(s1 * pEnergy->Gamma(), 2) + a0, -1.5);
         a2 = pow(pow(s2 * pEnergy->Gamma(), 2) + a0, -1.5);

         F[0] = (Bunch[i][0] - Bunch[j][0]) * (a1+a2);
         F[1] = (Bunch[i][2] - Bunch[j][2]) * (a1+a2);
         F[2] = ((s1*a1) + (s2*a2));

         for (int k = 0; k < 3; k++)
         {  Bunch[i][k*2+1] += K1 * F[k];
            Bunch[j][k*2+1] -= K1 * F[k];
         }
      }
   }
}

void xBeam::SpaceChargeBunch(xTime& time)
{
   static double dx, dy, ds, a3, F[3], impact;
   static doubleU m_2(1, m_^-2);
   doubleU K(((pEnergy->Z*U_e)^2) * time.ds * m_2 /
          (pEnergy->Momentum*pEnergy->Velocity*pEnergy->Gamma));
   double K1 = K(U1_);
   impact = pow(Impact(m_), -3.);

   for (int i = 0; i < Bunch.Number-1; i++)
   {	for (int j = i+1; j < Bunch.Number; j++)
      {  dx = Bunch[i][0] - Bunch[j][0];
         dy = Bunch[i][2] - Bunch[j][2];
         ds = (Bunch[i][4] - Bunch[j][4]) * pEnergy->Gamma();
         a3 = pow(ds*ds + dx*dx + dy*dy, -1.5);
         if (a3 < impact)
         {
            F[0] = dx * a3;
            F[1] = dy * a3;
            F[2] = ds * a3;

            for (int k = 0; k < 3; k++)
            {  Bunch[i][k*2+1] += K1 * F[k];
               Bunch[j][k*2+1] -= K1 * F[k];
            }
         }
         else Warning("Particle distance less than Impact Parameter : ", pow(a3, -1./3.),"[ m ]");
      }
   }
}

//---------------------------------------------------------------------------

void xBeam::Temperature()
{
   Tlong = pEnergy->A * pEnergy->Massa / U_k *
           (U_c^2) * (pEnergy->Beta^2);
   Ttrans= Tlong * (pEnergy->Gamma^2) * (Emit[0] / *pBetaH);
   Tlong *= Emit[2];

   Gamma1 = ((pEnergy->Z * pEnergy->Charge)^2) * Emit[3] /
      ( pEnergy->Gamma * (*pCirc) * ((Ttrans * Tlong)^0.5) );

   Gamma2 = ((pEnergy->Z * pEnergy->Charge)^2) / ( Tlong *
   ((Emit[0]* *pBetaH)^0.5) );
}

doubleU xBeam::GetGammaPi(doubleU tlong)
{
   doubleU teta;
   teta = (pEnergy->Z^4) * (pEnergy->Charge^4) /
          (U_pi^2) / (tlong^2) /
          ((*pBetaH^2)+(*pBetaV^2));

   return pEnergy->A * pEnergy->Massa / U_k * teta *
          (U_c^2) * (pEnergy->Beta^2) * (pEnergy->Gamma^2);
}

doubleU xBeam::GetGamma2(doubleU momentum)
{
   doubleU a(m_);
   a = (pEnergy->Z^2) * (pEnergy->Charge^2) /
   (pEnergy->A * pEnergy->Massa * (U_c^2) * (pEnergy->Beta^2)) /
   ((momentum^2) * U_pi);
   return (a^2) * 2. / (*pBetaH + *pBetaV);
}

doubleU xBeam::GetGamma3(doubleU emit)
{  /*
   return
   ( U_rp * 2 * (pEnergy->Z^2) / (pEnergy->A) / (pEnergy->Beta^2) *
     ( ( ( (*pBetaH * emit)+(*pBetaV * emit))^-0.5) -
       ( ( (*pBetaH * emit)+(*pBetaV * emit) +
           ((pEnergy->Gamma * *pCirc / Emit[3])^2))^-0.5)
     )
   )^0.5;
   */
   return
   ( U_rp * 2 * (Emit[3]^3) * (pEnergy->Z^2) * ((*pBetaH * emit)^0.5) /
     ( pEnergy->A * (pRing->TunesH^2) * (pEnergy->Gamma^3) *
      (pEnergy->Beta^2) * ((*pCirc)^2) * U_pi
     )
   );
}
doubleU xBeam::GetEquillibrium(doubleU emit)
{
   return ((pEnergy->Gamma^2) * emit * 2. / (*pBetaH + *pBetaV))^0.5;
}

void xBeam::KeepCross()
{
        int N = Number();
   if (N > 1000) return;
        Cross.Size(N, 4);
        xBeam& beam = *this;

   for (int i = 0; i < N; i++)
        Cross[i][0] = beam[i][4];
}

bool xBeam::Crossed()
{
        int N = Number();
   if (N > 1000) return true;
        Cross.Size(N, 4);
        xBeam& beam = *this;

   for (int i = 0; i < N; i++)
   {
      Cross[i][0]-= CellSize() * floor(Cross[i][0] / CellSize());
      Cross[i][1] = i;
        Cross[i][2] = beam[i][4];
      Cross[i][2]-= CellSize() * floor(Cross[i][2] / CellSize());
      Cross[i][3] = i;
   }

   bool trans;
   double d;
   for (int j = 0; j < 4; j += 2) do
   {
        trans = false;
      for (int i = 0; i < N-1; i++)
      {
        if (Cross[i][j] > Cross[i+1][j])
         {
                for (int k = 0; k < 2; k++)
            {
                d = Cross[i][j+k];
               Cross[i][j+k] = Cross[i+1][j+k];
               Cross[i+1][j+k] = d;
                                }
                trans = true;
         }
      }
   }while (trans);

   bool crossed = true;
   for (int k = -1; k < 2; k++)
   {
        bool equal = true;
      for (int i = 0; i < N; i++)
        if (i+k >= 0 && i+k < N)
                if (Cross[i][1] != Cross[i+k][3])
                equal = false;
      if (equal) crossed = false;
   }
   return crossed;
}

bool xBeam::Loss(xLattice& lat, int j, double loss, bool linj_loss)
{
	
   if (loss > 1) loss = 1;
    doubleU dn1;
      if (linj_loss)
   {
	   dn1 = Emit[3] / Bunch.Number;
	   Emit[3] -= dn1 * loss;
   }
   else
   {
	   Emit[3] *= exp(-loss / Bunch.Number);
   }
     
    
   if (Emit[3] < 10) Emit[3] = 10;

   if (fabs(xDistributor::Flatten()) < loss)
   {  int r;
      double mu;

      switch (GenerateOn)
      {case 0:
         SetIon(j);
         Matching(lat, j);
         PlusDisp(lat, j);
         Bunch[j][5] += Emit[5]();
         break;
       case 1:
         r = (int)(Number() * fabs(xDistributor::Flatten()));
         if (r == Number()) r = 0;
//         for (int n = 0; n < 6; n++) Bunch[j][n] = Bunch[r][n];
         Courant_Snyder(lat, r);

         mu = fabs(xDistributor::Flatten()) * 2. * M_PI;
         Bunch[j][0] = sqrt(Invs[r][0] * lat.betax()) * sin(mu)  + lat.Dx() * Bunch[j][5];
         Bunch[j][1] = sqrt(Invs[r][0] / lat.betax()) * (cos(mu) - lat.alfax()*sin(mu)) +
                                                                 lat.Dpx() * Bunch[j][5];
         mu = fabs(xDistributor::Flatten()) * 2. * M_PI;
         Bunch[j][2] = sqrt(Invs[r][1] * lat.betay()) * sin(mu)  + lat.Dy() * Bunch[j][5];
         Bunch[j][3] = sqrt(Invs[r][1] / lat.betay()) * (cos(mu) - lat.alfay()*sin(mu)) +
                                                                 lat.Dpy() * Bunch[j][5];
         if (benum == BUNCHED)
         {
            mu = fabs(xDistributor::Flatten()) * 2. * M_PI;
            Bunch[j][4] = sqrt(Invs[r][2]) * lat.B_s() * sin(mu);
            Bunch[j][5] = sqrt(Invs[r][2]) * cos(mu);
         } else
         {
            Bunch[j][4] = Bunch[r][4];
            Bunch[j][5] = Bunch[r][5];
         }
         break;
       case 2:
         if (Bunch.Number > 10)
         {  Bunch.Number--;
            if (j != Bunch.Number)
            for (r = 0; r < 6; r++)
               Bunch[j][r] = Bunch[Bunch.Number][r];
         }
         break;
      }
      return true;
   }
   return false;
}

//---------------------------------------------------------------------------
//************** from GRISHA***********************************************
//---------------------------------------------------------------------------
int xBeam::CalcStability()
{
   doubleU Z_t_sc_p(0, Ohm_*m_);
   G = 1.;
   B_f = 1.;
   doubleU a_x(m_);
   doubleU a_z(m_);
   a_x = ((Emit[0]*(*pBetaH))^0.5);
        a_z = ((Emit[1]*(*pBetaV))^0.5);
   a = 2.0*((a_x*a_z)^0.5);

   Z_t_sc_p = (*pz_0)*(*pR) /(pEnergy->Beta2 * pEnergy->Gamma2);

 if(a < (*pVac_a))
   {
      G = 1 + (2*U_Ln((*pVac_a)/a));
      Z_t_sc = Z_t_sc_p * ( (1/(a*a)) - (1/((*pVac_a)*(*pVac_a))) );
   }
 else Warning("Mean beam radius= ",a(),"  is more than chamber ", (*pVac_a)());

   Z_l_sc = ((*pz_0) * G * 0.5) / (pEnergy->Beta * pEnergy->Gamma2);
   Z_l = Z_l_i + Z_l_sc;
   Z_t = Z_t_i + Z_t_sc;

   return 1;
}
//------------------------------------------------------------
//---------    Beam stability parameters   -------------------
//------------------------------------------------------------

void xBeam::Keil_Shneile()
{
   doubleU e_p(938.2796, M_6*eV_);

   if(benum == BUNCHED)
      CalcBunch();

   I = (pEnergy->Z * Emit[3] * U_e) / (B_f * (*pTrev));

   if(benum == BUNCHED) I = I * (*pRF_h);
   KS =  (I * Z_l * pEnergy->Z * e_) / (F_l* e_p * 4.0 * pEnergy->A * pEnergy->Gamma*
         pEnergy->Beta2 * U_Abs((*pEta))* Emit[2]);
}
//------------------------------------------------------------

void xBeam::Tune_shift()
{
   B_f = 1.;

   if(benum == BUNCHED)
      CalcBunch();

   doubleU np = Emit[3];
   D_Q_x = (F_sc * pEnergy->Z * pEnergy->Z * U_rp * Emit[3])/(pEnergy->A * pEnergy->Beta2 * 2. * U_pi * pEnergy->Gamma3 * Emit[0]*(((Emit[1]/Emit[0])^0.5) + 1.0) * B_f);

   if(benum == BUNCHED)
      D_Q_x *= (*pRF_h);

   if(benum == BUNCHED)
      CalcBunch();

   D_Q_z = (F_sc * pEnergy->Z * pEnergy->Z * U_rp * Emit[3])/(pEnergy->A * pEnergy->Beta2 * 2. * U_pi * pEnergy->Gamma3 * Emit[1]*(((Emit[0]/Emit[1])^0.5) + 1.0) * B_f);

   if(benum == BUNCHED)
      D_Q_z *= (*pRF_h);

}
//------------------------------------------------------------

void xBeam::Shneil_Zotter()
{
   doubleU e_p(938.2796, M_6*eV_);
   doubleU temp = ( ((*pQ_x)-1) * (*pEta) + (*pKsy_x) ) * (Emit[2]^0.5);

   S_Q = (((temp*temp) + (D_Q_z*D_Q_z))^0.5);

   DM = (I*Z_t*(*pR)*pEnergy->Z*U_e)/(F_t * e_p * 8 * pEnergy->A * pEnergy->Gamma * pEnergy->Beta * (*pQ_x)* S_Q);
}
//------------------------------------------------------------

void xBeam::CalcBunch()
{
  doubleU e_p(938.2796, M_6*eV_);
  doubleU L_b (0, m_);
   s_s = (xLattice::B_s) * (Emit[2]^0.5);
   L_b = s_s * ((2*U_pi)^0.5);
   B_f = L_b /(*pRF_L_s);

   N_max = (L_b* L_b* L_b* pEnergy->Gamma2* (*pRF_h)* U_e *(*pRF_V)) /
          (pEnergy->Z* G* (*pR)* (*pR)* e_p *U_rp* 3.* U_pi);

   if(Emit[3] < N_max) Q_s = (*pRF_Q_s) * (((-Emit[3]/N_max) +1)^0.5);
   else Q_s = 0.0;
   Ep = U_pi * pRing->Energy.A * U_amu * pRing->Energy.Beta *
        pRing->Energy.Gamma * (Emit[2]^0.5) * s_s * U_c;
   //Ep = (U_pi * pRing->Energy.Beta *
   //     pRing->Energy.Gamma * (Emit[2]^0.5) * s_s)(m_);
}


double xBeam::CalcBucket(doubleU dp2)
{

   ksi2 = (U_e *pEnergy->Z * (*pBB_V_B)*(*pBB_T2))/
          (2.0*pEnergy->Momentum*U_Abs((*pEta))*(*pCirc));
   Amp2 = (((9.0*(ksi2-dp2)*(ksi2-dp2)+(12.0*ksi2*dp2))^0.5)- (3.0*(ksi2-dp2)))/2.0;

   Size_Bac = ((*pBB_T2)*pEnergy->Velocity) +
              ((U_Abs((*pEta))*pEnergy->Velocity*(*pCirc)*pEnergy->Momentum*Amp2)/
              ( U_e *pEnergy->Z * (*pBB_V_B)));
   s_s = Size_Bac;
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07
   if (Amp2())
   {
      T_s = (2.0*(*pBB_T2)/((U_Abs((*pEta))*(Amp2^0.5))))+
         (4.0*(*pCirc)*pEnergy->Momentum*(Amp2^0.5)/( U_e *pEnergy->Z * (*pBB_V_B)));

      t_1 =  (*pBB_T2)/(2.0*U_Abs((*pEta))*(Amp2^0.5));
      t_2 =  (*pCirc)*pEnergy->Momentum*(Amp2^0.5)/( U_e *pEnergy->Z * (*pBB_V_B));
   }
   return Size_Bac(m_);
}

void xBeam::SetBucket(vectorU& X)
{
   doubleU dp2;
   do
   {
      X[5] = Sigma[5]*Gaussian();
      dp2 = X[5] * X[5];
      CalcBucket(dp2);
   }while(Amp2 > (pRing->BarHeight*pRing->BarHeight));

   if(Amp2() == 0.0)
   {
      X[4] = iRing.T2(s_)*Flatten()/(2.0*iRing.Trev(s_));
      X[5] = 0.0;
   }else
   {
      double rand1 = double(rand()) / RAND_MAX;
      doubleU T_int(s_);
      T_int = T_s*rand1;
      GetBucket(X, T_int);
   }
}

void xBeam::GetBucket(vectorU& X, doubleU& T_int)
{
   if (T_int < t_1)
   {
      X[4] = U_Abs(pRing->Eta)*pRing->Energy.Velocity*(Amp2^0.5)*T_int;
      X[5] = Amp2^0.5;
   }else
   {if (T_int < (t_1 + (2.0*t_2)))
   {
      X[4] = (U_Abs(pRing->Eta)*pRing->Energy.Velocity*(Amp2^0.5)*T_int)
         -(U_Abs(pRing->Eta)*pRing->Energy.Velocity*U_e*pRing->Energy.Z*pRing->V_B
         *(T_int - t_1)*(T_int - t_1)/
         (2.0*pRing->Circ*pRing->Energy.Momentum));
      X[5] = (Amp2^0.5) - (U_e*pRing->Energy.Z*pRing->V_B*(T_int - t_1)/
         (pRing->Circ*pRing->Energy.Momentum));
   }else
   {if(T_int < ((3.0*t_1) + (2.0*t_2)))
   {
      X[4] = U_Abs(pRing->Eta)*pRing->Energy.Velocity*(Amp2^0.5)*
         ((2.0*t_1) + (2.0*t_2) - T_int);
      X[5] = -1.0*(Amp2^0.5);
   }else
   {if(T_int < ((3.0*t_1) + (4.0*t_2)))
   {
      X[4] = (-0.5*pRing->Energy.Velocity*pRing->T2)-
           (U_Abs(pRing->Eta)*pRing->Energy.Velocity*(Amp2^0.5)*(T_int - (3.0*t_1) - (2.0*t_2)))
           +(U_Abs(pRing->Eta)*pRing->Energy.Velocity*U_e*pRing->Energy.Z*pRing->V_B
           *(T_int - (3.0*t_1)-(2.0*t_2))*(T_int - (3.0*t_1)-(2.0*t_2))/
           (2.0*pRing->Circ*pRing->Energy.Momentum));
      X[5] = (U_e*pRing->Energy.Z*pRing->V_B*(T_int - (3.0*t_1) - (2.0*t_2))/
           (pRing->Circ*pRing->Energy.Momentum)) - (Amp2^0.5);
   }else
   {
      X[4] = U_Abs(pRing->Eta)*pRing->Energy.Velocity*(Amp2^0.5)*(T_int - T_s);
      X[5] = Amp2^0.5;
   }}}}
}

doubleU xBeam::Bdp2(doubleU s,doubleU dp)
{
doubleU Bp2 (0, U1_);
doubleU Centr(m_);
doubleU Am2;
//double As;
//As = fabs(s());
Centr = (*pBB_T2)*pEnergy->Velocity / 2.0;
if(fabs(s()) > Centr())
{
   Am2 = dp*dp + (2.0* U_e *pEnergy->Z * (*pBB_V_B) * (U_Abs(s) - Centr))/
         (U_Abs((*pEta))*pEnergy->Velocity*(*pCirc)*pEnergy->Momentum);

}else
{
   Am2 = dp*dp;
}

Bp2 = ((Am2*Am2) +(3.0*ksi2*Am2))/
       (3.0*(Am2+ksi2));

return Bp2;
}
//---------------------------------------------------------------------------
/*
void xBeam::SimCoolLum()
{
//Lum.L = N_b * (Emit[3]^2)/(((sigma_x2*sigma_z2)^0.5)* (*pTrev) * U_pi* 4.0);
Lum.L = 0.0;
int i, j;
xLattice lattice;
lattice = iRing.CalcLattice(iRing.Matrix);


doubleU r_max(m_), r1(m_), r2(m_), r_step(m_), r_i(m_);

r_max = 4.0 * Sigma[0];

r_step = r_max/Lum.step;

int N_slice;           // Number of particles in slice

for (j = 0; j < Lum.step; j++)
  {
   N_slice = 0;
   r1 = r_step*j;
   r2 = r_step*(j+1);

   xBeam& beam = *this;
   for (i = 0; i < Bunch.Number; i++)
     {
        r_i = ((beam(i,0)*beam(i,0)* Lum.beta_x_1/lattice.betax +
              beam(i,2)*beam(i,2)*Lum.beta_z_1/lattice.betay)^0.5);
      if ((r_i>r1)&&(r_i<=r2))
      N_slice++;
          }
   Lum.L += ((Emit[3]*N_slice/Bunch.Number)^2)/(U_pi * r_step*r_step*(2.0*j +1.0) * (*pTrev));

   }
Lum.L *= N_b;
}
//-----------------------------------------------
void xBeam::Luminosity()
{
Lum.L = 1.0;
Lum.head = true;
Lum.head_bunched = true;
Lum.thesame = true;
Lum.cross = false;
Lum.merge = false;
Lum.horizontal = true;
Lum.vertical = false;

if(Lum.horizontal)
{
Sigma_2[0] = Emit[0]*Lum.beta_x_1;
Sigma_2[2] = Emit[1]*Lum.beta_z_1;

Lum.sigma_x2 = 10000.* Lum.Eh*Lum.beta_x_2/2.0;
Lum.sigma_z2 = 10000.* Lum.Ev*Lum.beta_z_2/2.0;
}
if(Lum.vertical)
{
Sigma_2[2] = 10000.* Lum.Eh*Lum.beta_x_1/2.0;
Sigma_2[0] = 10000.* Lum.Ev*Lum.beta_z_1/2.0;

Lum.sigma_z2 = 10000.* Lum.Eh*Lum.beta_x_2/2.0;
Lum.sigma_x2 = 10000.* Lum.Ev*Lum.beta_z_2/2.0;
}


doubleU pi2 =2 *U_pi;
doubleU LL(s_^-1);

//if(Lum.head) Lum.phy = 0.0;

if(Lum.head)
{
    if(bunched&&Lum.head_bunched)
    {
     if(Lum.thesame)  //if(Lum.thesame&& !Lum.s)
        {
        Lum.L = N_b * (Emit[3]^2)/(((Sigma_2[0]*Sigma_2[2])^0.5)* (*pTrev) * U_pi* 4.0);
        if(iDynamics.Algoritm == 1)
        SimCoolLum();
        }
    else
        {
        Sigma_2[4] = s_s * s_s;
        Lum.sigma_s2 = Lum.s * Lum.s;
        Lum.L = integral() * N_b * Emit[3] * Lum.N * ((pEnergy->Beta2 + Lum.beta2 +
                (pEnergy->Beta * Lum.beta * 2.0) )^0.5)/((*pTrev) * pi2 * pi2 * pi2);        }
    }
}
if(Lum.cross)
{
    if(bunched&&Lum.head_bunched)
    {
     if(Lum.thesame)  //if(Lum.thesame&&!Lum.s)
        {
        Lum.L = N_b * (Emit[3]^2) / (((Sigma_2[0]*Sigma_2[2])^0.5) * (*pTrev) * U_pi * 4.0);

        Lum.L *= ( ((-pEnergy->Beta2 * U_Sin(Lum.phy)*U_Sin(Lum.phy) + 1) /
                 ((s_s * U_Sin(Lum.phy) / ( ( (Sigma_2[0])^0.5) * U_Cos(Lum.phy) ) ) * (s_s*U_Sin(Lum.phy) /
                 ( ( (Sigma_2[0])^0.5) * U_Cos(Lum.phy) ) ) + 1.) )^0.5)/U_Cos(Lum.phy);
        }
    else
        {
        Lum.L = N_b * Emit[3] * Lum.N * ((pEnergy->Beta2 + Lum.beta2 +
                pEnergy->Beta * Lum.beta * U_Cos(Lum.phy * 2.0) * 2.0 -
                pEnergy->Beta2 * Lum.beta2 * U_Sin(Lum.phy * 2.0) * U_Sin(Lum.phy * 2.0))^0.5) /
                ((*pTrev) * pi2 * pi2 * pi2);
        Sigma_2[4] = s_s * s_s;
        Lum.sigma_s2 = Lum.s*Lum.s;
        Lum.L *= integral();
        }
    }
}
if(Lum.merge)
    {
    Sigma_2[4] = s_s * s_s;
    Lum.sigma_s2 = Lum.s*Lum.s;
    Lum.L = Lcross(1.0);
    }
}
//----------------------integral=--------------------------------

doubleU xBeam::integral()
{
doubleU I (m_^-2);
doubleU interval(m_);
doubleU step;
interval = Lum.interval*(Lum.s + s_s);
step = interval/Lum.step;
doubleU z = interval*-1.;
do
{
    if(Lum.cross)I += Lcross(z)*step;
    if(Lum.head)I += Lhead(z)*step;
    z += step;
}while(z < interval);
return I;
}
//---------------------------------------------------------------

doubleU xBeam::Lcross(doubleU z)
{
doubleU Lz = 1.0;
doubleU sig1_x2,sig1_y2,sig2_x2,sig2_y2;

sig1_x2 = Sigma_2[0]*(z*z/(Lum.beta_x_1*Lum.beta_x_1) + 1);
sig1_y2 = Sigma_2[2]*(z*z/(Lum.beta_z_1*Lum.beta_z_1) + 1);
sig2_x2 = Lum.sigma_x2*(z*z/(Lum.beta_x_2*Lum.beta_x_2) + 1);
sig2_y2 = Lum.sigma_z2*(z*z/(Lum.beta_z_2*Lum.beta_z_2) + 1);
//---------Calculation of multiplicator-----------------
doubleU MUL, COS, SIN;
MUL = (U_pi*2)^1.5;
MUL /= ((sig1_y2 + sig2_y2)^0.5);
COS = (sig1_x2 + sig2_x2) * U_Cos(Lum.phy) * U_Cos(Lum.phy)*
      (pEnergy->Beta2 * Lum.sigma_s2 + Lum.beta2 * Sigma_2[4]);
SIN = sig1_x2*sig2_x2*(pEnergy->Beta-Lum.beta)*
      (pEnergy->Beta-Lum.beta)*U_Sin(Lum.phy)*U_Sin(Lum.phy);
MUL /= ((COS+SIN)^0.5);
//----------------exp--------------------
doubleU A1, B1, B2, C1, D1, E1, ExP;
A1 = pEnergy->Beta2 / Sigma_2[4] + Lum.beta2 / Lum.sigma_s2;
B1 = (pEnergy->Beta / Sigma_2[4] + Lum.beta / Lum.sigma_s2) * U_Sin(Lum.phy);
B2 = (pEnergy->Beta * -1. / Sigma_2[4] + Lum.beta / Lum.sigma_s2) * U_Cos(Lum.phy);
C1 = ((Sigma_2[0]^-1) + (Lum.sigma_x2^-1))*U_Cos(Lum.phy) * U_Cos(Lum.phy) +
     ((Sigma_2[4]^-1) + (Lum.sigma_s2^-1))*U_Sin(Lum.phy) * U_Sin(Lum.phy);
D1 = ((Sigma_2[0]^-1) + (Lum.sigma_x2^-1) + (Sigma_2[4]^-1) - (Lum.sigma_s2^-1))*
        U_Sin(Lum.phy) * U_Cos(Lum.phy);
E1 = ((Sigma_2[4]^-1) + (Lum.sigma_s2^-1)) * U_Cos(Lum.phy) * U_Cos(Lum.phy) +
     ((Sigma_2[0]^-1) + (Lum.sigma_x2^-1)) * U_Sin(Lum.phy) * U_Sin(Lum.phy);

ExP =  z * z * (E1 - B2 * B2 / A1 - (B1 * B2 + A1 * D1) * (B1 * B2 + A1 * D1)/
       (A1 * (A1 * C1 - B1 * B1))) * -0.5;

Lz = MUL * U_Exp(ExP);
return Lz;
}
//
//---------------------------------------------------------------
doubleU xBeam::Lhead(doubleU z)
{
doubleU Lz = 1.0;
doubleU sig1_x2,sig1_y2,sig2_x2,sig2_y2;
sig1_x2 = Sigma_2[0] * (z * z / (Lum.beta_x_1 * Lum.beta_x_1) + 1);
sig1_y2 = Sigma_2[4] * (z * z / (Lum.beta_z_1 * Lum.beta_z_1) + 1);
sig2_x2 = Lum.sigma_x2 * (z * z / (Lum.beta_x_2 * Lum.beta_x_2) + 1);
sig2_y2 = Lum.sigma_z2 * (z * z / (Lum.beta_z_2 * Lum.beta_z_2) + 1);


doubleU MUL,ExP;
MUL = ((U_pi*2)^3./2.);
MUL /= (((sig1_y2 + sig2_y2) * (sig1_x2 + sig2_x2) *
      (pEnergy->Beta2 * Lum.sigma_s2 + Lum.beta2 * Sigma_2[4]))^0.5);

ExP = z * z * (pEnergy->Beta + Lum.beta) * (pEnergy->Beta + Lum.beta) * -0.5/
      (Lum.beta2 * Sigma_2[4] + pEnergy->Beta2 * Lum.sigma_s2);

Lz = MUL * U_Exp(ExP);
return Lz;
}
*/

//---------------------------------------------------------------

CollidingBeam::CollidingBeam()
{
   bunched = false;
   Sigma_x._(m_);
   Sigma_y._(m_);
   Sigma_s._(m_);
   Emit[0]._(u_6*m_);
   Emit[1]._(u_6*m_);
}

int CollidingBeam::OnGet()
{
  int uk;
  int Reference;
  double Kinetic;

  Energy.A = Data.OnGet(65, 6, 1);
  Energy.Z = Data.OnGet(65, 7, 1);
  uk       = Data.OnGetI(65, 9, 0);

  if (uk > 3)
     Energy.PerNucleon = false;

   switch(uk)
   {case 0: case 4: UKinetic =     eV_; break;
    case 1: case 5: UKinetic = k_3*eV_; break;
    case 2: case 6: UKinetic = M_6*eV_; break;
    case 3: case 7: UKinetic = G_9*eV_; break;
   }

   Reference = Data.OnGetI(65, 1, 0);
   switch(Reference)
   {case 0:
      Energy.Gamma = Data.OnGet(65, 2, 1.1);
   	if (Energy.Set(U_GAMMA)) return 2;
   	break;
    case 1:
      Energy.Beta = Data.OnGet(65, 3, 0.9);
   	if (Energy.Set(U_BETA)) return 2;
   	break;
    case 2:
      Kinetic = Data.OnGet(65, 4, 1.);
   	if (Energy.Set(U_KINETIC, Kinetic, UKinetic)) return 2;
   	break;
    case 3:
      Energy.Momentum = Data.OnGet(65, 5, 1.);
   	if (Energy.Set(U_MOMENTUM)) return 2;
   	break;
	}

   Emit[0] = Data.OnGet(66, 1, 1.);
   Emit[1] = Data.OnGet(66, 2, 1.);
   Sigma_s = Data.OnGet(66, 3, 1.);
   Emit[3] = Data.OnGet(66, 4, 1.);
   bunched = Data.OnGetB(66, 5, false);

   Lattice.betax = Data.OnGet(67, 1, 1.);
   Lattice.betay = Data.OnGet(67, 2, 1.);

   return 0;
}

int CollidingBeam::OnSet()
{
   Data.OnSet(65, 2, Energy.Gamma);
   Data.OnSet(65, 3, Energy.Beta);
   Data.OnSet(65, 4, Energy.Kinetic(UKinetic));
   Data.OnSet(65, 5, Energy.Momentum);

   Sigma_x = (Emit[0] * Lattice.betax )^0.5;
   Sigma_y = (Emit[1] * Lattice.betay )^0.5;

   return 0;
}
