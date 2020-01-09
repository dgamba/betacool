//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xTarget.h"
//---------------------------------------------------------------------------

xMaterial::xMaterial()
{
   gr._(g_^-1);
   Z._(U1_);
   A._(U1_);
   Ro._(g_/(cm_^3));
   l._(m_3*m_);
   N._(1/(cm_^3));
   Rox._(g_/(cm_^2));
   Nx._(1/(cm_^2));
   Omega0._(U1_);
   ChiA._(U1_);

   Emax._(eV_);
   I._(eV_);

   //Urban model
   R._(U1_);
   Nex1._(U1_);
   Eex1._(eV_);
   Nex2._(U1_);
   Eex2._(eV_);
   Nion._(U1_);
   GGG._(U1_);
   //
   Teta2._(U1_);
   dE._(eV_);
   Estr2._(eV_^2);
   Use_dE = true;

   Hhor._(m_);
   Hver._(m_);
   dEhor._(m_);
   dEver._(m_);

   ElectronCapture   = true;
   SingleScattering  = true;
   NuclearReaction   = true;
   InteractionEvents = false;
   CrossSection._(barn_);
}

void xMaterial::Density(int ud)
{
   if (ud)
      Ro = N*A/(U_NA*gr);
   else
      N = Ro*(U_NA*gr)/A;
   Rox = Ro*l;
   Nx = N * l;
}

void xMaterial::RMSANGLE(U_Energy&e)
{
   doubleU db(U1_);
   doubleU Uin(U1_);
   doubleU a1(cm_);
   doubleU a2(U1_);
   doubleU Ksi(U1_);
   doubleU r0(1.3E-15,m_);
   doubleU a0(0.529E-8,cm_);

   //Moliere parameters
   doubleU ZSP(U1_);
   doubleU ZEP(U1_);
   doubleU ZXP(U1_);
   doubleU BC(cm_^-1);
   doubleU BCC(6702.23,cm_^-1);
   doubleU ChiCC(U1_);
   doubleU ChiC(U1_);

   ZSP = Z*(Z+1.0)/A;
   ZEP = ZSP*U_Ln(Z^(-2.0/3.0));
   ZXP = ZSP*U_Ln(1.0+(3.34*e.Z*Z/(137.0*e.Beta)));

   BC = BCC*Ro(g_/(cm_^3))*ZSP*U_Exp((ZEP-ZXP)/ZSP);

   Omega0 = BC*l*e.Z*e.Z/(e.Beta*e.Beta);

   ChiCC = 0.39612*U_Pow((ZSP* Ro(g_/(cm_^3))),0.5);
   ChiC = ChiCC*e.Z*U_Pow(l(cm_),0.5)/(e.A*931.49432*e.Gamma*e.Beta);
   ChiA = ChiC/U_Pow((1.167*Omega0),0.5);

   if (A <= 7)  Uin = -3.6 - ((A-1)/6);
   else
	if (A <= 16) Uin = -4.6 - (0.4*(A-7)/9);
   else
                Uin = -5.8;

   db =( U_Ln((1130*e.Beta2)/((Z^(4./3.))*((1-e.Beta)^2))) +
         Uin - (e.Beta2/2) ) / Z;
   a1 = 0.885*a0*(((Z^(2./3))+(e.Z^(2./3)))^0.5);
   a2 = a1 / (r0* ((A^(1./3))+(e.A^(1./3))));
   Ksi = (a2^2) / (1.13*(1+(3.33*((Z*e.Z/(137*e.Beta))^2))));
   Teta2 = 2*U_pi*Nx*((Z*e.Z*U_rp/e.A/(e.Beta2*e.Gamma))^2)*
           (U_Ln(Ksi)-1+db);
}

void xMaterial::Bethe(U_Energy&e)
{
   doubleU Ra;
   Ra._(M_6 * eV_*(cm_^2)/g_);
   doubleU Ra1;
   Ra1._(eV_);
   doubleU Psi;
   Psi._(eV_);
   doubleU F1(1.0,U1_);
   doubleU F2(0.0,U1_);
   doubleU temp(U1_);
//   double initrand;
   long randinit = -1;

   Eex2 = 10*Z*Z*Ra1;
   /*
   if (Z == 1)
      I = 13.6*Ra1;
   else
      I = 10*Z*Ra1;
   */
   I = 16.*U_Pow(Z,0.9)*Ra1;
   if (Z > 2)
   {
      F2 = 2.0/Z;
      F1 = 1.0 - F2;
      temp = (U_Ln(I/Ra1) - (F2*U_Ln(Eex2/Ra1)))/F1;
      Eex1 = Ra1 * U_Exp(temp);
   }
   else Eex1 = I;

   Emax = 2* U_me*(U_c^2)*e.Beta*e.Beta*e.Gamma*e.Gamma/
          (1+(2*e.Gamma*U_me/(e.A*e.Massa)+((U_me/(e.A*e.Massa))^2)));
   GGG = Emax/(Emax + I);
   Psi = 0.1535*Ra*e.Z*e.Z*Z*Rox/(e.Beta*e.Beta*A);
   dE = 2*Psi*(U_Ln(Emax/I)-e.Beta*e.Beta);
   Estr2 = Psi*Emax*(1-e.Beta*e.Beta/2);
   Nex1 = dE * (1.0 - R) * (F1/Eex1) *
          (U_Ln(Emax/Eex1) - e.Beta*e.Beta)/(U_Ln(Emax/I)-e.Beta*e.Beta);
   Nex2 = dE * (1.0 - R) * (F2/Eex2) *
          (U_Ln(Emax/Eex2) - e.Beta*e.Beta)/(U_Ln(Emax/I)-e.Beta*e.Beta);
   Nion = dE * R * Emax/
          (I*(Emax + I)* U_Ln((Emax + I)/I));
   xDistributor::ran1(&randinit);
}

void xMaterial::Deviations(U_Energy&Energy, xLattice&L)
{
   Hhor = (((1+L.alfax)^2)*(L.Dx^2)/L.betax)+(2*L.alfax*L.Dx*L.Dpx) + (L.betax*(L.Dpx^2));
   Hver = (((1+L.alfay)^2)*(L.Dy^2)/L.betay)+(2*L.alfay*L.Dy*L.Dpy) + (L.betay*(L.Dpy^2));
   dPP2 = ((Energy.Gamma/(1+Energy.Gamma))^2)*(Estr2 / ((Energy.A * Energy.Kinetic)^2));
   dPB2 = ((Energy.Gamma/(1+Energy.Gamma))^2)*(dE*dE / ((Energy.A *Energy.Kinetic)^2));
   dEhor = (L.betax*Teta2/2) + (Hhor * dPB2);
   dEver = (L.betay*Teta2/2) + (Hver * dPB2);
}

doubleU xMaterial::Lifetime(xRing&ring, xBeam&beam, xLattice&L)
{

   doubleU sec(0, cm_^2); //Electron capture life-time in accordance with Schlachter
   doubleU sss(0, cm_^2); //Single scattering life time
   doubleU snr(0, cm_^2); //Nuclear reaction life time
   doubleU sie(0, cm_^2); //Total interaction cross section

   if(ElectronCapture)
   {
      doubleU E_eff(eV_);
      E_eff = ring.Energy.Kinetic/((ring.Energy.Z^0.7)*(Z^1.25));
      sec = 1.1e-8*(1.0-exp(-0.037*pow(E_eff(k_3*eV_),2.2))) *
            (1.0 - exp(-2.44e-5*pow(E_eff(k_3*eV_),2.6))) /
            pow(E_eff(k_3*eV_),4.8);
   }
   if (SingleScattering)
   {
      doubleU Tetacc;
      Tetacc = ring.AcceptV/(U_pi * L.betay)*(1-((beam.Emit[1]/ring.AcceptV)^2));
      sss = 4*U_pi * ((U_rp * Z * ring.Energy.Z / ring.Energy.A)^2) /
            (ring.Energy.Beta3 * ring.Energy.Gamma2 * Tetacc);
   }
   if (NuclearReaction)
      snr = U_pi() * 1e-26 * pow(pow(ring.Energy.A(),1./3.)+ pow(A(),1./3.),2) /
         ring.Energy.Beta2();

   if (InteractionEvents)
      sie = CrossSection;

   return (sss + sec + snr + sie) * Nx / ring.Trev;
}
//---------------------------------------------------------------------------

xTarget::xTarget(xBeam& beam)
{
   ID = 3;
//   OpticName = "XTARGET";
   EffectName = "XTARGET";
   CrossNumber = 0;

//   Length = 0;
   pBeam = &beam;
   Cross._(cm_^2);
   probability._(1, U1_);

   PShiftX._(m_3*m_);
   XWidth._(m_3*m_);
   PelletX._(m_3*m_);
   PelletY._(m_3*m_);
   Velosity._(m_3*m_/s_);
   Interval._(m_3*m_);
   Ntar = 0;
}

int xTarget::OnGet()
{
   dens = Data.OnGetI(36,1,0);
   A  = Data.OnGet(36,3,1);
   Z  = Data.OnGet(36,4,1);
   l  = Data.OnGet(36,5,1);
   Ro = Data.OnGet(36,6,1e10);
   N  = Data.OnGet(36,7,1e10);
   mat =       Data.OnGetI(36,8,0);
   Nx    =     Data.OnGet(36,9,1e10);
   Teta2 = pow(Data.OnGet(36,10,1), 2);
   dE    =     Data.OnGet(36,11,1);
   Estr2 = pow(Data.OnGet(36,12,1), 2);
   Use_dE=     Data.OnGetB(36,13,1);

   itype   = Data.OnGetI(38,1,0);
   //imethod = Data.OnGet(38,2,0);
   //Cross   = Data.OnGet(38,3,0);
   //num_MC  = Data.OnGet(38,4,1000);
   //Gaussian = Data.OnGet(38,5,1);
   Urban  = Data.OnGetI(38,5,0);
   Plural = Data.OnGetI(38,6,0);
   R = Data.OnGet(38,7,0);
   probability = 1.0;
/*
   PShiftX = Data.OnGet(40,1,0);
   XWidth =  Data.OnGet(40,7,2);
   PelletX = Data.OnGet(40,2,0.01);
   if(XWidth() == 0) XWidth = 0.5*PelletX;
   PelletY = Data.OnGet(40,3,0.01);
   Velosity= Data.OnGet(40,4,10);
   Interval= Data.OnGet(40,5,1);
*/
   if (itype == 2)
   {
      l *= (U_pi / 6.)^(1./3.);
      PelletX = l;
      PelletY = l;

      PShiftX = Data.OnGet(35,1,0);
      XWidth =  Data.OnGet(35,7,2);
      if(XWidth() == 0) XWidth = 0.5*PelletX;
      Velosity= fabs(Data.OnGet(35,6,10));

      Interval= Data.OnGet(35,13,1);
      if(Interval < PelletY)
         Interval = PelletY;
   }

   ElectronCapture  = Data.OnGetB(39,1,1);
   SingleScattering = Data.OnGetB(39,2,1);
   NuclearReaction  = Data.OnGetB(39,3,1);
   InteractionEvents= Data.OnGetB(39,4,1);
   CrossSection     = Data.OnGet(39,5,1);

   return 0;
}

int xTarget::OnSet()
{
   Lattice.betax = Data.OnGet(37,1,1);
   Lattice.alfax = Data.OnGet(37,2,1);
   Lattice.Dx    = Data.OnGet(37,3,0);
   Lattice.Dpx   = Data.OnGet(37,4,0);
   Lattice.betay = Data.OnGet(37,5,0);
   Lattice.alfay = Data.OnGet(37,6,0);
   Lattice.Dy    = Data.OnGet(37,7,0);
   Lattice.Dpy   = Data.OnGet(37,8,0);

   U_Energy& Energy = *pBeam->pEnergy;
   if (mat)
   {  Density(dens);
      RMSANGLE(Energy);
      Bethe(Energy);
   }else
      Rox = Nx*A/(U_NA*gr);
   Deviations(Energy, Lattice);

   vectorU X(m_,U1_ );
   xTime t(iTime);
   if (itype > 0)
   {
      X[4] = 0;
      X[5] = ((2.*iBeam.Emit[2])^0.5);
      X[0] = X[5]*Lattice.Dx;
      X[1] = ((iBeam.Emit[0]*2./Lattice.betax)^0.5)+ (X[5]*Lattice.Dpx);
      X[2] = X[5]*Lattice.Dy;
      X[3] = ((iBeam.Emit[1] *2./Lattice.betay)^0.5)+ (X[5]*Lattice.Dpy);
      if (itype == 1)
         probability = CylinderProbability(t,X,iRing);
      else
         probability = FiberProbability(t,X,iRing);
   }else
      probability = 1;

   Data.OnSet(36,6,Ro);
   Data.OnSet(36,7,N);
   Data.OnSet(36,9,Nx);
   Data.OnSet(36,10,sqrt(Teta2()));
   Data.OnSet(36,11,dE);
   Data.OnSet(36,12,sqrt(Estr2()));
   Data.OnSet(36,14,(Nx * probability)(cm_^-2));
   Data.OnSet(39,6,(iBeam.N_b * iBeam.Emit[3] * Nx* probability / iRing.Trev)((cm_^-2)*(s_^-1)));

   return 0;
}

vectorU xTarget::Rates(xTime&time, xBeam&beam, xRing&ring)
{
   probability = U_1;
   vectorU rates(s_^-1);
   vectorU X(m_,U1_ );
   double k = 1;
   if (beam.benum == BUNCHED) k = 2;

   if (itype > 0)
   {
      X[4] = 0;
      X[5] = ((2.*beam.Emit[2])^0.5);
      X[0] = X[5]*Lattice.Dx;
      X[1] = ((beam.Emit[0]*2./Lattice.betax)^0.5)+ (X[5]*Lattice.Dpx);
      X[2] = X[5]*Lattice.Dy;
      X[3] = ((beam.Emit[1] *2./Lattice.betay)^0.5)+ (X[5]*Lattice.Dpy);
      if (itype == 1)
         probability = CylinderProbability(time,X,ring);
      else
         probability = FiberProbability(time,X,ring);
   }else
      probability = 1;

   rates[0] = dEhor/(beam.Emit[0]*ring.Trev) * probability;
   rates[1] = dEver/(beam.Emit[1]*ring.Trev) * probability;
   rates[2] = dPP2 /(k*beam.Emit[2]*ring.Trev) * probability;
   //20.11.06
   if (Use_dE&&(beam.benum == COASTING))
      rates[2]+= dPB2 /(k*beam.Emit[2]*ring.Trev) * probability;

   if (Loss)
      rates[3] = - Lifetime(ring, beam, Lattice) * probability;
   iColl.Luminosity = beam.N_b * beam.Emit[3] * Nx * probability/ ring.Trev;

   return rates;
}

//---------------------------------------------------
doubleU xTarget::PelletProbability(xBeam&beam)
{
   doubleU probability(1, U1_);
   doubleU BeamRadius2(m_^2);
   doubleU BeamRadius(m_);
   doubleU Horda2(m_^2);
   doubleU Horda(m_);
   doubleU MeanHorda(0,m_);
   int i;
   int Nstep = 100;



   BeamRadius2 = ((((beam.Emit[0]*Lattice.betax) + (beam.Emit[2]*Lattice.Dx*Lattice.Dx))*
                     beam.Emit[1]*Lattice.betay)^0.5) * 4;
   BeamRadius = BeamRadius2^0.5;

   if(XWidth() == 0.0)
   {
      Horda2 = PShiftX^2;
      if(Horda2 < BeamRadius2)
         Horda2 = BeamRadius2 - Horda2;
      else
         Horda2 = 0;
      probability = (PelletX * PelletY * 2 * (Horda2^0.5)) /
                    (Interval * U_pi * BeamRadius2);
   }
   else
   {
      MeanHorda = 0;
      for (i = 1; i < Nstep; i++)
      {
      Horda = (Nstep-2*i)*BeamRadius/Nstep;
      MeanHorda += U_Exp(-0.5*((Horda - PShiftX)/XWidth)*((Horda - PShiftX)/XWidth))*
                        ((BeamRadius2-Horda*Horda)^0.5);
      }
      probability = (PelletX * PelletY * 2 * (2 * MeanHorda/(Nstep*((2.0*U_pi)^0.5)*XWidth))) /
                    (Interval * U_pi * BeamRadius);
   }
   if(probability > U_1)
      probability = U_1;

   return probability;
}
//-------------------------------------------------------------

doubleU xTarget::PelletCross(xTime&t, vectorU ion,xRing&ring)
{
   doubleU Nc(U1_);
   doubleU Ihor(m_);
   doubleU Ivert(m_);
   doubleU Ah(m_);
   doubleU Av(m_);
   doubleU Xb(m_);
   doubleU Xs(U1_);
   doubleU Yb(m_);
   doubleU Ys(U1_);
   int Nstep = 100;
   int i;
   doubleU Sum(U1_);
   doubleU Arg1(U1_);
   doubleU Arg2(U1_);
   doubleU Xsh(m_);

   Xsh = PShiftX + XWidth * xDistributor::Gaussian();

   doubleU gammax(m_^-1);
   gammax = (1. + (Lattice.alfax^2))/Lattice.betax;
   doubleU gammay(m_^-1);
   gammay = (1. + (Lattice.alfay^2))/Lattice.betay;

   Xb = ion[0] - Lattice.Dx * ion[5];
   Xs = ion[1] - Lattice.Dpx* ion[5];
   Yb = ion[2] - Lattice.Dy * ion[5];
   Ys = ion[3] - Lattice.Dpy* ion[5];

   Ihor = ((Lattice.betax*Xs*Xs)+(2.*Lattice.alfax*Xb*Xs)+(gammax*Xb*Xb));
   Ivert = ((Lattice.betay*Ys*Ys)+(2.*Lattice.alfay*Yb*Ys)+(gammay*Yb*Yb));
   Ah = (Ihor*Lattice.betax)^0.5;
   Av = (Ivert*Lattice.betay)^0.5;

   if((Ah - Xsh + (PelletX/2.0) + (Lattice.Dx * ion[5]))() < 0)
   {
      Nc = 0.0;
      return Nc;
   }

   Sum = 0.0;
   for(i = 0; i < Nstep; i++)
   {
     Arg1 = (((Av - (PelletY/2.0))*i/Nstep)+(PelletY/2.0))/Av;
     Arg2 = (((Av - (PelletY/2.0))*i/Nstep)-(PelletY/2.0))/Av;
     if (U_Abs(Arg1)() >= 1) Arg1 = 0.99999;
     if (U_Abs(Arg2)() >= 1) Arg2 = 0.99999;
     Sum += U_Acos(Arg1)- U_Acos(Arg2);
   }

   Nc = t.dt/ring.Trev;
   Nc *= (Av - (PelletY/2.0)) * Sum / (Nstep*Interval);
   Arg1 = (Xsh + (PelletX/2.0)- (Lattice.Dx * ion[5])) / Ah;
   Arg2 = (Xsh - (PelletX/2.0)- (Lattice.Dx * ion[5])) / Ah;
   if (U_Abs(Arg1)() >= 1) Arg1 = 0.99999;
   if (U_Abs(Arg2)() >= 1) Arg2 = 0.99999;

   Nc *= U_Acos(Arg1)- U_Acos(Arg2);
   Nc *= 2.0/(U_pi*U_pi);
   if (Nc < 0.0) Nc = 0.0;

   return Nc;
}
//-------------------------------------------------------------
//-------------------------------------------------------------

doubleU xTarget::FiberProbability(xTime&t, vectorU ion,xRing&ring)
{
   doubleU probability(1, U1_);
   doubleU Horda2(m_^2);
   doubleU Horda(0,m_);
   doubleU MeanHorda(0,m_);
   int i;
   int Nstep = 100;

   doubleU Ihor(m_);
   doubleU Ah(m_);
   doubleU Xb(m_);
   doubleU Xs(U1_);
   doubleU Xmin(m_);
   doubleU Xmax(m_);
   doubleU Arg1(U1_);
   doubleU Arg2(U1_);
   doubleU Phymin(U1_);
   doubleU Phymax(U1_);
   doubleU Phystep(U1_);
   doubleU Xstep(m_);
   doubleU gammax(m_^-1);
   gammax = (1. + (Lattice.alfax^2))/Lattice.betax;

   Xb = ion[0] - Lattice.Dx * ion[5];
   Xs = ion[1] - Lattice.Dpx* ion[5];

   Ihor = ((Lattice.betax*Xs*Xs)+(2.*Lattice.alfax*Xb*Xs)+(gammax*Xb*Xb));
   Ah = (Ihor*Lattice.betax)^0.5;

   if(Ah() == 0)
   {
      if(((Lattice.Dx * ion[5]) > (PShiftX - XWidth))&&
         ((Lattice.Dx * ion[5]) < (PShiftX + XWidth)))
           probability = PelletY / Interval;
      else probability = 0.0;
      return probability;
   }

   Xmin = -Ah + Lattice.Dx * ion[5];
   Xmax =  Ah + Lattice.Dx * ion[5];

   if(Xmin > (PShiftX - XWidth)) Phymax = U_pi;
   else
   {
      Arg2 = (PShiftX - XWidth - (Lattice.Dx * ion[5]))/Ah;
      if (U_Abs(Arg2)() <= 1.)
         Phymax = U_Acos(Arg2);
   }
   if(Xmax < (PShiftX + XWidth)) Phymin = 0.0;
   else
   {
      Arg1 = (PShiftX + XWidth - (Lattice.Dx * ion[5]))/Ah;
      if (U_Abs(Arg1)() <= 1.)
         Phymin = U_Acos(Arg1);
   }

   MeanHorda = 0;
   if (U_Abs(Arg1)() <= 1. && U_Abs(Arg2)() <= 1.)
   {
      Phystep = (Phymax - Phymin) / Nstep;
      for (i = 1; i < Nstep; i++)
      {
          Horda = (Ah*U_Cos(Phymin + i*Phystep)) +(Lattice.Dx * ion[5])-PShiftX;
          MeanHorda += 2.0*((XWidth*XWidth-Horda*Horda)^0.5);
      }
      MeanHorda *=U_Abs(Phystep)/U_pi;
   }
   probability = PelletY / Interval;
   if(probability > U_1)
      probability = U_1;

   probability *= (PelletX * MeanHorda) / (XWidth*XWidth*U_pi);

   if(probability > U_1)
      probability = U_1;

   return probability;
}
//---------------------------------------------------------
doubleU xTarget::CylinderProbability(xTime&t, vectorU ion,xRing&ring)
{
   doubleU probability(1, U1_);
   doubleU Horda2(m_^2);
   doubleU Horda(0,m_);
   doubleU MeanHorda(0,m_);
   int i;
   int Nstep = 100;

   doubleU Ihor(m_);
   doubleU Ah(m_);
   doubleU Xb(m_);
   doubleU Xs(U1_);
   doubleU Xmin(m_);
   doubleU Xmax(m_);
   doubleU Arg1(U1_);
   doubleU Arg2(U1_);
   doubleU Phymin(U1_);
   doubleU Phymax(U1_);
   doubleU Phystep(U1_);
   doubleU Xstep(m_);
   doubleU gammax(m_^-1);
   gammax = (1. + (Lattice.alfax^2))/Lattice.betax;

   Xb = ion[0] - Lattice.Dx * ion[5];
   Xs = ion[1] - Lattice.Dpx* ion[5];

   Ihor = ((Lattice.betax*Xs*Xs)+(2.*Lattice.alfax*Xb*Xs)+(gammax*Xb*Xb));
   Ah = (Ihor*Lattice.betax)^0.5;

   if(Ah() == 0)
   {
      probability = 0.0;
      return probability;
   }

   Xmin = -Ah + Lattice.Dx * ion[5];
   Xmax =  Ah + Lattice.Dx * ion[5];

   if ((Xmin < l/2) && (Xmax >- l/2))                             //01.02.2011
   {
      if (Xmin > -l/2 ) Phymax = U_pi;
      else
      {
         Arg2 = (-l/2  - (Lattice.Dx * ion[5]))/Ah;
         Phymax = U_Acos(Arg2);
      }
      if(Xmax < l/2 ) Phymin = 0.0;
      else
      {
         Arg1 = ( l/2 - (Lattice.Dx * ion[5]))/Ah;
            Phymin = U_Acos(Arg1);
      }

      Phystep = (Phymax - Phymin) / Nstep;
      MeanHorda = 0;

      for (i = 1; i < Nstep; i++)
      {
          Horda = (Ah*U_Cos(Phymin + i*Phystep)) +(Lattice.Dx * ion[5]);
          MeanHorda += 2.0*((l*l/4 - Horda*Horda)^0.5);
      }
      MeanHorda *=U_Abs(Phystep)/U_pi;
      probability = MeanHorda / l;

      if(probability > U_1)
         probability = U_1;
   }else
      probability = 0;
   return probability;
}
//---------------------------------------------------------

void xTarget::Kick (xTime&time, xBeam&beam, xRing&ring)
{
   beam.RMS(beam.Emit,Lattice);
   doubleU NN(U1_);
   doubleU dE_urban(eV_);
   static long rancal = 121256;

   double loss = 0;
   if (Loss)
      loss = (time.dt * Lifetime(ring, beam, Lattice))(U1_);
   double turn =  (time.dt/ring.Trev)(U1_);
   double teta = ((Teta2*turn)^0.5)(U1_);
   double dpp2 = ((dPP2 *turn  )^0.5)(U1_);
   double dpb2 = ((dPB2 ^0.5)*turn)(U1_);
   double xx, xxx;
   doubleU Ntot(U1_);
   Ntot = 0;
   Ntar = 0;
   for (int j = 0; j < beam.Number(); j++)
   {
      if (itype == 0)
      {
         if (!beam.Loss(Lattice, j, loss,false))
         {
            beam[j][1] += teta * xDistributor::Gaussian();
            beam[j][3] += teta * xDistributor::Gaussian();
            beam[j][5] += dpp2 * xDistributor::Gaussian();
            if (Use_dE)
            beam[j][5] -= dpb2;
            CrossNumber += (int)turn;
            iColl.Luminosity = beam.N_b * beam.Emit[3] * Nx / ring.Trev;
         }
      }else
      {  if (itype == 1)
            NN = CylinderProbability(time, beam(j), ring) * turn;
         else
            NN = FiberProbability(time, beam(j), ring) * turn;

         if(!beam.Loss(Lattice, j, (loss*NN()/turn),false))
         {
            if (itype == 1) //!!!!!!!!!!!!!!!!!!!!!!!!!
            {
               Ntot +=NN;
               beam[j][1] += ((Teta2*NN/2)^0.5)(U1_) * xDistributor::Gaussian();
               beam[j][3] += ((Teta2*NN/2)^0.5)(U1_) * xDistributor::Gaussian();
               beam[j][5] += ((dPP2 *NN  )^0.5)(U1_) * xDistributor::Gaussian();
               if (Use_dE)
               beam[j][5] -= ((dPB2 ^0.5)*NN)(U1_);
            }else
            {  
               Ntot += NN();
               NN += fabs(xDistributor::Flatten());

            //Longitudinal momentum variation
               if (Urban)
               {
                  dE_urban = 0.0;
                  for (int kk = 0; kk < int(NN()); kk++)
                  {  Ntar++;
                     dE_urban += double(xDistributor::Poisson(Nex1()))*Eex1 +
                                 double(xDistributor::Poisson(Nex2()))*Eex2;
                     for (int ll = 0; ll < xDistributor::Poisson(Nion()); ll++)
                     {
                        dE_urban += I/(1.0 - (GGG()*xDistributor::ran1(&rancal)));
                     }
                  }
                  beam[j][5] -= (ring.Energy.Gamma/(1+ring.Energy.Gamma)*
                                dE_urban /(ring.Energy.A *ring.Energy.Kinetic))(U1_);
               }else
               {
                  beam[j][5] += ((dPP2 *NN  )^0.5)(U1_) * xDistributor::Gaussian();
                  if (Use_dE)
                  beam[j][5] -= ((dPB2 ^0.5)*NN)(U1_);
               }

            //Plural scattering in transvese plane
               if (Plural)
               {
                  for (int nn = 0; nn < int(NN()*Omega0()); nn++)
                  {
                     do
                     {
                        xx = double(rand()) / RAND_MAX;
                     } while(xx==0.0);
                     xxx = xDistributor::Flatten();
                     beam[j][1] += ChiA(U1_)*sqrt(1.0/xx-1.0)*cos(3.1415926*xxx);
                     beam[j][3] += ChiA(U1_)*sqrt(1.0/xx-1.0)*sin(3.1415926*xxx);
                  }
               }else
               {
                  beam[j][1] += ((Teta2*NN/2)^0.5)(U1_) * xDistributor::Gaussian();
                  beam[j][3] += ((Teta2*NN/2)^0.5)(U1_) * xDistributor::Gaussian();
               }
            }
         }
      }
   }
   if (itype > 0)
   {
      iColl.Luminosity = beam.N_b * beam.Emit[3] * Nx * Ntot/(turn*beam.Number()*ring.Trev);
      Ntar = Ntar / beam.Number() / time.dt();
   }
   iColl.Events += iColl.Luminosity * CrossSection * time.dt;

}

