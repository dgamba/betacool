//---------------------------------------------------------------------------
#include "stdafx.h"
#define _CRT_SECURE_NO_WARNINGS
#include "xDynamic.h"           
#include "xRing.h"
#include "xBeam.h"

//---------------------------------------------------------------------------
xRing::xRing(xBeam& beam)
{
   pBeam = &beam;
   Index = -1;
   Circ._(m_);
   Arc._(m_);
   AcceptH._(m_);
   AcceptV._(m_);
   AcceptDp._(U1_);
   AcceptS._(m_);
   AcceptHinj._(m_);
   AcceptVinj._(m_);
   R_m._(m_);
   Trev._(s_);
   TLife._(s_);
   BetaH._(m_);
   BetaV._(m_);
   Eta._(U1_);
   DispH._(m_);
   DispV._(m_);
   // for RF
   h._(U1_);  //harmonic number
   V._(k_3*V_);   //RF voltage
   NonLinear = false;
   h_low._(U1_);  //harmonic number
   V_low._(k_3*V_);   //RF voltage

   //Tsweep._(s_);
   Inductor = false;
   Vind._(V_);
   L_s._(m_);  //Separatrix length
   Q_s._(U1_); //Synchrotron tune
   //for Vacuum
   a_m._(cm_); // mean chamber radius
   extend._(cm_);
   //20.11.06 For barrier bucket
   V_B._(k_3*V_);                          //Voltage amplitude
   T1._(s_);
   T2._(s_);
   S1._(m_);
   S2._(m_);
   BarHeight._(U1_);                       //Barrier height

}

int xRing::OnGet()
{
   Energy.A = Data.OnGet(11,6,1);
   Energy.Z = Data.OnGet(11,7,1);
   TLife    = Data.OnGet(11,8,1e6);
   int uk   = Data.OnGetI(11,9,0);
   if (uk > 3)
     Energy.PerNucleon = false;

   switch(uk)
   {case 0: case 4: UKinetic =     eV_; break;
    case 1: case 5: UKinetic = k_3*eV_; break;
    case 2: case 6: UKinetic = M_6*eV_; break;
    case 3: case 7: UKinetic = G_9*eV_; break;
   }

   switch((int)Data[11][1])
   {case 0:
   	Energy.Gamma = Data.OnGet(11,2,1.01);
   	if (Energy.Set(U_GAMMA)) return 2;
   	break;
    case 1:
   	Energy.Beta = Data.OnGet(11,3,0.99);
   	if (Energy.Set(U_BETA)) return 2;
   	break;
    case 2:
   	if (Energy.Set(U_KINETIC, Data.OnGet(11,4,1), UKinetic)) return 2;
   	break;
    case 3:
   	Energy.Momentum = Data.OnGet(11,5,1);
   	if (Energy.Set(U_MOMENTUM)) return 2;
   	break;
	}

   Circ    = Data.OnGet(12,1,100);
   Arc     = Circ;
   GammaTr = Data.OnGet(12,2,2);
   TunesH  = Data.OnGet(12,3,3.33);
   TunesV  = Data.OnGet(12,4,3.31);
   HromatH = Data.OnGet(12,5,0);
   HromatV = Data.OnGet(12,6,0);

   h = Data.OnGet(17,1,1);
   V = Data.OnGet(17,2,1);
   Size_s = Data.OnGet(17,5,100);
   NonLinear = Data.OnGetB(17,6,false);
   h_low = Data.OnGet(17,7,1);
   V_low = Data.OnGet(17,8,1);
   //dPini = Data.OnGet(17,6,0);
   //dPfin = Data.OnGet(17,7,0);
   //Tsweep= Data.OnGet(17,8,0);

   if (Data.OnGetB(12,9,false))
     Imagenary = -1;
   else                     
     Imagenary = 1;

   LatticeFile = Data.OnGetI(13,1,0);
   ReduceExtend = Data.OnGetI(13,2,2);
   extend = Data.OnGet(13,3,10);
   Autoskip= Data.OnGetB(13,5, false);

	Trev  = Circ / Energy.Velocity;
   BetaH = Circ / (2 * U_pi * TunesH);
   BetaV = Circ / (2 * U_pi * TunesV);
   Eta = (1/Energy.Gamma2) - (Imagenary/(GammaTr*GammaTr));
   R_m = (Circ)/(2.0*U_pi);
   DispH = BetaH/TunesH;
   DispV = 0; //BetaV/TunesV;

   InjectionPoint = Data.OnGetI(18,1,0);
   if (InjectionPoint == 0)
   {
      LATTICE.betax = BetaH;
      LATTICE.betay = BetaV;
      LATTICE.Dx    = DispH;
      LATTICE.Dy    = DispV;
   }else
   {
      LATTICE.betax = Data.OnGet(18,2,1);
      LATTICE.betay = Data.OnGet(18,3,1);
      LATTICE.alfax = Data.OnGet(18,4,0);
      LATTICE.alfay = Data.OnGet(18,5,0);
      LATTICE.Dx    = Data.OnGet(18,6,1);
      LATTICE.Dy    = Data.OnGet(18,7,1);
      LATTICE.Dpx   = Data.OnGet(18,8,0);
      LATTICE.Dpy   = Data.OnGet(18,9,0);
   }
   xEffect::SetLattice(LATTICE);

// for Vacuum
   a_m = Data.OnGet(43,1,1e-10);

   doubleU Q;
   doubleU E_p(938279600, eV_);
   Q = (U_e * h * U_Abs(Eta) * Energy.Z * V) /(2* U_pi * E_p * Energy.A * Energy.Gamma);
   L_s = Circ/h;
   Q_s = (Q^0.5) / Energy.Beta;

   if(pBeam->benum == BUNCHED)
      xLattice::B_s = Circ*U_Abs(Eta)/(2*U_pi*Q_s);     // syncrotron function,
   else                                                 // used for initial
      xLattice::B_s = Circ;                             // beam distibution

   Matrix20 = Data.OnGetI(20,1,0);
   Index20  = Data.OnGetI(20,2,0);
   Value20  = Data.OnGetI(20,3,0);
   for (int i = 0; i < 6; i++)
   for (int j = 0; j < 6; j++)
      Matrix[i][j] = Data.OnGet(19,i*6+j+1,1);

   Bucket.Analytic = Data.OnGetB(26,8,true);
   Bucket.Stationary = Data.OnGetB(26,9,false);
   Bucket.Moving = Data.OnGetB(26,11, false);
   Bucket.Show3D = Data.OnGetB(26,17, false);
   Bucket.HourGlass = Data.OnGet(26,18, 0);
   Bucket.cut = Data.OnGetB(26,19, false);
   Bucket.vertex = fabs(Data.OnGet(26,20, 0))*2.;

   Bucket.Harmonic1  = Data.OnGetB(26,22, false);
   Bucket.Number1    = Data.OnGetI(26,23, 22);
   Bucket.Harmonic2 = Data.OnGetB(26, 24, false);
   Bucket.Number2 = Data.OnGetI(26, 25, 66);
   Bucket.Induction = Data.OnGetB(26, 26, false);
   Bucket.TimeStart = Data.OnGet(26, 27, 0);
   Bucket.TimeFinish = Data.OnGet(26,28, 1);
   Bucket.Impedance = complex<double>(Data.OnGet(26, 29, 1), Data.OnGet(26, 30, 0));
   Bucket.LongSpaceCharge = Data.OnGet(26, 31, false);

   // Remove Shen modification
   // Bucket.OverlapBucket = Data.OnGetI(26, 35, 0);
   // Bucket.Issin = Data.OnGetB(26, 36, false);
   
   if (pBeam->benum == BUCKET && Bucket.Moving)
   {
      h     = Bucket.Number1;
      h_low = Bucket.Number2;
      Bucket.Analytic = Bucket.Harmonic1;
   }
   return OnRun();
}

int xRing::OnRun()
{
   AcceptH = Data.OnGet(12,7,1e-5);
   AcceptV = Data.OnGet(12,8,1e-5);
   AcceptDp = Data.OnGet(12,10,1e-2);
   AcceptS = Data.OnGet(12,11,1e10);
   AcceptHinj = Data.OnGet(14,9,1e-5);
   AcceptVinj = Data.OnGet(14,10,1e-5);
   
   Inductor = Data.OnGetB(16,7,false);
   Vind = Data.OnGet(16,8,0.001);

   V_B= Data.OnGet(26,1,1);                             //Voltage amplitude
   T1 = Trev*Data.OnGet(26,2,0.1);                      //RF duration in T0
   T2 = Trev*Data.OnGet(26,3,0.7);                      //Gap duration in T0
   S1 = Circ * (Data.OnGet(26,15,0.5)/2);               //Kicker gap in m
   S2 = Circ * (T2/2) / Trev;                           //Gap width in meters
   //if (S1 < S2) S1 = S2;
   BarHeight = (2.0*T1*Energy.Z*V_B*U_e/(Trev*Energy.Velocity*Energy.Momentum* U_Abs(Eta)))^0.5;

   Bucket.bCoeff = true;
   Bucket.IBSnorma = Data.OnGetB(26,16, false);
   Bucket.Periodical = Data.OnGetB(26,21, false);

   if (iBeam.benum == BUCKET)
   {
      if (Bucket.Moving)
      {
         Bucket.BarData.SetSeparators("\t,");
         if (!Bucket.BarData.Load(Data.OnGetC(26,12,"NoFile")))
            return 1;
         doubleU t(0,s_);
         Bucket.SetMoving(t);
      }else
      if (Bucket.Stationary)
      {
         Bucket.BarData.SetSeparators("\t,");
         if (!Bucket.BarData.Load(Data.OnGetC(26,10,"NoFile")))
            return 1;
         Bucket.SetBucket(Bucket.BarData);
      }else
      {
         double s0 =  Circ(m_) / 2;
         double s1 = (Circ * (T1 + T2/2) / Trev )(m_);
         double s2 = S2(m_);
         Bucket.barrier.SetNumber(5);
         Bucket.barrier[0].SetBarrier(-s0, -s1, 0);
         Bucket.barrier[1].SetBarrier(-s1, -s2, V_B(V_));
         Bucket.barrier[2].SetBarrier(-s2,  s2, 0);
         Bucket.barrier[3].SetBarrier( s2,  s1,-V_B(V_));
         Bucket.barrier[4].SetBarrier( s1,  s0, 0);

         Bucket.barrier[0].UP=(Bucket.barrier[0].s1-Bucket.barrier[0].s2)*Bucket.barrier[0].Ub;
         for (int i = 1; i < 5; i++)
         Bucket.barrier[i].UP=(Bucket.barrier[i].s1-Bucket.barrier[i].s2)*Bucket.barrier[i].Ub
                             + Bucket.barrier[i-1].UP;
         Bucket.Ucoeff = V_B(k_3*V_);
         Bucket.UPcoeff = (V_B * (Bucket.barrier[1].s2 - Bucket.barrier[1].s1))(k_3*V_*m_);
      }
   }   
   return 0;
}
//------------------------------------------------------------

int xRing::OnSet()
{
   static bool once = true;
   if(once)
   {  once = false;

      if(LatticeFile < 2)
      {  if (LatticeFile == 0)
         {  if(OutputMAD(Data.OnGetC(13,4,"NoFile")))return 1;
         }else
         {  if(InputMAD(Data.OnGetC(13,6,"NoFile")))return 1;
         }
         switch (ReduceExtend)
         {case 0:
             ReduceLattice(Data.OnGetC(13,8,"NoFile"));
             break;
          case 1:
             ExtendLattice();
             break;
         }
         if (LatticeFile == 1)
            SetLattice();
         if (InjectionPoint == 1)
            LATTICE = GetLattice(0);
      }
   }

	Data.OnSet(11,2, Energy.Gamma);
	Data.OnSet(11,3, Energy.Beta);
	Data.OnSet(11,4, Energy.Kinetic(UKinetic));
	Data.OnSet(11,5, Energy.Momentum);
   Data.OnSet(12,1, Circ);

	Data.OnSet(16,1, R_m);
	Data.OnSet(16,2, BetaH);
	Data.OnSet(16,3, BetaV);
	Data.OnSet(16,4, DispH);
   Data.OnSet(16,5, Trev);
   Data.OnSet(16,6, Eta);

   Data.OnSet(17,3, L_s);
   Data.OnSet(17,4, Q_s);

   Data.OnSet(26,4,BarHeight);
   Data.OnSet(26,18,Bucket.HourGlass);

   if (LatticeFile == 1)
   {  switch (Matrix20)
      {case 0: Matrix =  GetMatrix(iTime, Index20);
         break;
       case 1: Matrix = RingMatrix(iTime, Index20);
         break;
      }
      double m2d;
      for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
      {  switch (Value20)
         {case 0: m2d = real(Matrix[i][j]);
            break;
          case 1: m2d = imag(Matrix[i][j]);
            break;
          case 2: m2d = abs(Matrix[i][j]);
            break;
          case 3: m2d = arg(Matrix[i][j]);
            break;
         }
         if ((fabs(m2d) < 1e-5) || (fabs(m2d) > 1e5))
            m2d = BData::char2double(BData::double2char(m2d,3));
         Data.OnSet(19,i*6+j+1, m2d);
      }
   }
   Det20 = 1. - real(Matrix.Det());
   Data.OnSet(20,4, Det20);


   Data.OnSet(20,5, LATTICE.Qx());
   Data.OnSet(20,6, LATTICE.Qy());
	return 0;
}
//---------------------------------------------------------------------------
void xRing::ReduceLattice(char* filename)
{  int r,r1;
   doubleU length(m_),length1(m_);
   matrixU matrix, matrix1;
   BData Reduce;
   xRing& ring = *this;
   xTime time(ring);
   Reduce.AddSeparators(",");
   if (!Reduce.Load(filename))
   {  Warning("Ring | Lattice Structure | Reduce file does not exist", PressEnter);
      return;
   }
   Reduce.Col(Reduce.Row(),2);

   length = 0;
   matrix.Diag(1,0);

   r = (int)Reduce[Reduce.Row()-1][0];
   for (int i = r; i < Number(); i++)
   {
      if (LatticeFile == 1)
         matrix = ring[i].GetMatrix(time, Energy) * matrix;
      length += ring[i].Length;
   }
   ring.Number(r+1);

   for (int j = Reduce.Row()-1; j > 0; j--)
   {
      length1 = 0;
      matrix1.Diag(1,0);
      r = (int)Reduce[j][0];
      r1= (int)Reduce[j-1][0];
      for (int i = r1; i < r; i++)
      {
         if (LatticeFile == 1)
            matrix1 = ring[i].GetMatrix(time, Energy) * matrix1;
         length1 += ring[i].Length;
      }
      for (int i1 = r1; i1 < r-1; i1++)
         ring.LRemove(r1+1);
      ring[r1].Matrix = matrix1;
      ring[r1].Length = length1;
   }

   r = (int)Reduce[0][0];
   for (int i3 = 0; i3 < r; i3++)
   {
      if (LatticeFile == 1)
         matrix = ring[i3].GetMatrix(time, Energy) * matrix;
      length += ring[i3].Length;
   }
   for (int i2 = 0; i2 < r-1; i2++)
      ring.LRemove(1);

   ring[Number()-1].Matrix = matrix;
   ring[Number()-1].Length = length;

   for (int i4 = 0; i4 < Number(); i4++)
   {  Reduce[i4](0).sr = " , ";
      if (ring[i4].LGetSize() && (LatticeFile == 1))
         Reduce[i4](1).ch = ring[i4][0]->OpticName;
      else
         Reduce[i4](1).ch = " ";
   }
   Reduce.Save(filename);

}

xOpticsLibrary* xRing::GetOptics(chars name)
{
   //if (name == "RFCAVITY")                        return new xRFcavity;
   //if (name == "XTARGET")                         return &iTarget;
   if (name == "XECOOL")                          return &iEcool;
   if (name == "XLCOOL")                          return &iLaserCooling;
   if((name == "SBEND") || (name == "RBEND"))     return new xBend;
   if (name == "EBEND")                           return new xEBend;
   if (name == "FFIELD")                          return new xFField;
   if((name == "QUADRUPOLE") || (name == "QUAD")) return new xQuadrupole;
   if (name == "SEXTUPOLE")                       return new xSextupole;
   if (name == "SOLENOID")                        return new xSolenoid;
   if (name == "RFCAVITY")                        return new xRFcavity;
   return NULL;
}

void xRing::ExtendLattice()
{
   for (int i = Number()-1; i >=0; i--)
   {  xRing& ring = *this;
      if (ring[i].Length > extend)
      {  int Steps = (int)((ring[i].Length / extend)(U1_) + 0.5);
         if (Steps < 1) Steps = 1;
         ring[i].Length /= Steps;
         for (int k = 0; k < ring[i].LGetSize(); k++)
         {  ring[i][k]->Length = ring[i].Length;
            ring[i][k]->exit = false;
            ring[i][k]->entr = true;
         }

         for (int j = 0; j < Steps-1; j++)
         {  LInsert(i+1);
            ring[i+1].EL_LATTICE = ring[i].EL_LATTICE;
            ring[i+1].RE_MATRIX  = ring[i].RE_MATRIX;
            ring[i+1].Lattice    = ring[i].Lattice;
            ring[i+1].Length     = ring[i].Length;
            ring[i+1].OpticName  = ring[i].OpticName;

            ring[i+1].LSetSize(ring[i].LGetSize());

            for (int k = 0; k < ring[i].LGetSize(); k++)
            {  ring[i+1][k] = GetOptics(ring[i][k]->OpticName);
               if (ring[i+1][k] == NULL) break;
               ring[i+1][k]->Length     = ring[i][k]->Length;
               ring[i+1][k]->OpticName  = ring[i][k]->OpticName;
               ring[i+1][k]->EL_FORCES  = ring[i][k]->EL_FORCES;
               ring[i+1][k]->EL_FIELDS  = ring[i][k]->EL_FIELDS;
               ring[i+1][k]->EL_MATRIX  = ring[i][k]->EL_MATRIX;
               ring[i+1][k]->RE_MATRIX  = ring[i][k]->RE_MATRIX;
               ring[i+1][k]->entr = false;
               if (j == Steps-2)
                  ring[i+1][k]->exit = true;
               else
                  ring[i+1][k]->exit = false;

               if (ring[i][k]->OpticName == "SBEND")
               {
                  xBend* bend1 = (xBend*)(ring[i+1][k]);
                  xBend* bend0 = (xBend*)(ring[i][k]);
                  if (j == 0)
                  {  bend1->E2     = bend0->E2;
                     bend0->E2     = 0;
                     bend1->HGAP2  = bend0->HGAP2;
                     bend0->HGAP2  = 0;
                     bend1->FINT2  = bend0->FINT2;
                     bend0->FINT2  = 0;
                     bend0->ANGLE /= Steps;
                  }
                  bend1->ANGLE = bend0->ANGLE;
                  bend1->K1 = bend0->K1;
                  bend1->K2 = bend0->K2;
               }else
               if (ring[i][k]->OpticName == "EBEND")
               {
                  xEBend* ebend1 = (xEBend*)(ring[i+1][k]);
                  xEBend* ebend0 = (xEBend*)(ring[i][k]);
                  ebend1->FB  = ebend0->FB;
                  ebend1->FE  = ebend0->FE;
                  ebend1->Ro  = ebend0->Ro;
                  ebend1->n   = ebend0->n;
               }else
               if (ring[i][k]->OpticName == "FFIELD")
               {
                  xFField* ebend1 = (xFField*)(ring[i+1][k]);
                  xFField* ebend0 = (xFField*)(ring[i][k]);
                  ebend1->FB  = ebend0->FB;
                  ebend1->FE  = ebend0->FE;
                  ebend1->Ro  = ebend0->Ro;
                  ebend1->n   = ebend0->n;
               }else
               if (ring[i][k]->OpticName == "QUADRUPOLE")
               {
                  xQuadrupole* quad1 = (xQuadrupole*)(ring[i+1][k]);
                  xQuadrupole* quad0 = (xQuadrupole*)(ring[i][k]);
                  quad1->K1   = quad0->K1;
                  quad1->TILT = quad0->TILT;
               }else
               if (ring[i][k]->OpticName == "SEXTUPOLE")
               {
                  xSextupole* sext1 = (xSextupole*)(ring[i+1][k]);
                  xSextupole* sext0 = (xSextupole*)(ring[i][k]);
                  sext1->K2   = sext0->K2;
               }else
               if (ring[i][k]->OpticName == "SOLENOID")
               {
                  xSolenoid* sol1 = (xSolenoid*)(ring[i+1][k]);
                  xSolenoid* sol0 = (xSolenoid*)(ring[i][k]);
                  sol1->Ks   = sol0->Ks;
               }else
               if (ring[i][k]->OpticName == "RFCAVITY")
               {
                  xRFcavity* rfc1 = (xRFcavity*)(ring[i+1][k]);
                  xRFcavity* rfc0 = (xRFcavity*)(ring[i][k]);
                  rfc1->H  = rfc0->H;
                  rfc1->V  = rfc0->V;
               }
            }
         }
      }
   }
}

bool xRing::ReadMADfile(char* filename)
{
   BData MadIn;
   MadIn.SetSeparators(":=,()*!");
   if (MadIn.Load(filename))
   {
      char fn[128];
      strcpy(fn, filename);
      fn[strlen(fn)-3] = 'u';
      fn[strlen(fn)-2] = 's';
      fn[strlen(fn)-1] = 'e';

// Remove comments "!"
      for (int i = 0; i < MadIn.Row(); i++)
      {  for (int j = 0; j < MadIn.Col(i); j++)
            if(MadIn[i](j).sr == "!")
            {  if (j == 0)
               {  MadIn.LRemove(i--);
                  break;
               }  else
               {  MadIn.Col(i+1,j+1,i);
                  break;
               }
            }
      }

// Remove spaces and UpCase
      for (int i1 = 0; i1 < MadIn.Row(); i1++)
         for (int j1 = 0; j1 < MadIn.Col(i1); j1++)
         {  MadIn[i1](j1).ch.RemoveChar();
         /*
            if (MadIn[i1](j1).ch.Length() == 0)
            {  for (int j3 = j1+1; j3 < MadIn.Col(i1); j3++)
                  MadIn[i1](j3-1) = MadIn[i1](j3);
               MadIn.Col(i1+1, MadIn.Col(i1)-1, i1);
            }
         */   
            MadIn[i1](j1).ch.UpCase();
         }
      //MadIn.Save(fn);

/*
// Remove variables :=
      for (int k1 = 0; k1 < MadIn.Row(); k1++)
      {  if ((MadIn[k1](0).sr == ":") && (MadIn[k1](1).sr == "="))
         {  for (int k2 = k1+1; k2 < MadIn.Row(); k2++)
               for (int j2 = 1; j2 < MadIn.Col(k2); j2++)
                  if (MadIn[k2](j2).ch == MadIn[k1](0).ch)
                      MadIn[k2][j2] = MadIn[k1][2];

         }
      }
*/

// Replace end of lines "&"
      for (int i3 = 0; i3 < MadIn.Row(); i3++)
      {
         int col = MadIn.Col(i3)-1;
         chars ch = MadIn[i3](col).ch;
         int len = ch.Length()-1;
         if (len && ch[len] == '&')
         {  MadIn[i3](col).ch[len]='\0';
            MadIn[i3](col).sr = ",";
            MadIn.Col(i3+1, col+2, i3);
            MadIn[i3](col+1).ch = "&";
            MadIn[i3](col+1).sr = "\n";
            for (int j3 = 1; j3 < MadIn.Col(i3+1); j3++)
               MadIn[i3+1](j3-1) = MadIn[i3+1](j3);
            MadIn.Col(i3+2, MadIn.Col(i3+1)-1, i3+1);
         }
      }
      //MadIn.Save(fn);

      for (int i3 = 0; i3 < MadIn.Row(); i3++)
      {  bool more;
         do
         {  more = false;
            int col = MadIn.Col(i3)-1;
            if (MadIn[i3](col).ch == "&")
            {  MadIn.Col(i3+1,col + MadIn.Col(i3+1), i3);
               for (int j3 = 0; j3 < MadIn.Col(i3+1); j3++)
                  MadIn[i3](j3+col) = MadIn[i3+1](j3);
               MadIn.LRemove(i3+1);
               more = true;
            }
         }while (more);
      }
      //MadIn.Save(fn);

// Replace multiplication "*" in LINEs
      for (int i4 = 0; i4 < MadIn.Row(); i4++)
      {  if (MadIn.Col(i4) > 2)
         if (MadIn[i4](1).ch == "LINE")
         {  for (int j4 = 0; j4 < MadIn.Col(i4); j4++)
            if ((MadIn[i4](j4).sr == "*") && (MadIn[i4][j4] != EmptyData))
            {  int m = (int)MadIn[i4][j4];
               MadIn[i4](j4) = MadIn[i4](j4+1);
               for (int k4 = 0; k4 < m-2; k4++)
               {  MadIn[i4].LInsert(j4+1);
                  MadIn[i4](j4+1) = MadIn[i4](j4);
               }
            }
         }
      }

// looking for USE
      for (int i5 = 0; i5 < MadIn.Row(); i5++)
         if (MadIn[i5](0).ch == "USE")
            MadLine[0](0).ch = MadIn[i5](1).ch;

//expand of MadLine
      bool more;
      do
      {
         more = false;
         for (int i6 = 0; i6 < MadLine.Row(); i6++)
         {  for (int j6 = 0; j6 < MadIn.Row(); j6++)
            {  if (MadIn.Col(j6) > 2)
               {  chars ch = MadLine[i6](0).ch;
                  bool NoReverse = true;
                  if (ch[0] == '-')
                  {  NoReverse = false;
                     ch.RemoveChar('-');
                  }
                  if ((MadIn[j6](0).ch == ch) &&
                      (MadIn[j6](1).ch == "LINE"))
                  {  MadLine[i6](0).ch = MadIn[j6](3).ch;
                     more = true;
                     for (int k6 = 4; k6 < MadIn.Col(j6)-1; k6++)
                     {  if (NoReverse) i6++;
                        MadLine.LInsert(i6);
                        MadLine[i6](0).ch = MadIn[j6](k6).ch;
                     }
                  }
               }
            }
         }
      }while (more);
      //MadIn.Save(fn);

//Replace nikes on optics
      for (int i7 = 0; i7 < MadLine.Row(); i7++)
         for (int j7 = 0; j7 < MadIn.Row(); j7++)
            if(MadLine[i7](0).ch == MadIn[j7](0).ch)
            {  MadLine[i7]       =  MadIn[j7];
               break;
            }

      MadLine.Save(fn);
      return false;
   }else
      return true;
}

bool xRing::InputMAD(char* filename)
{
   if (ReadMADfile(filename))
      return true;

   int index1 = 0;
   Circ = 0;

   for (int i = 0; i < MadLine.Row(); i++)
   {
      if(LGetSize() <= index1) LSetSize(index1*2);
      xOpticsLibrary* plib;
      xRingOptics& opt = *LList[index1];
      opt.LSetSize(0);
      opt.Length = 0;
      int index2 = 0;
      bool lib1 = false;

      for (int j = 2; j < MadLine.Col(i); j++)
      {  lib1 = true;
//         chars ch = MadLine[i](j).ch;
//         if (ch[ch.Length()-1] == 'L')
         if (MadLine[i](j).ch == "L")
         {  opt.Length = MadLine[i][j+1];
            Circ += opt.Length;
            break;
         }
      }

      plib = GetOptics(MadLine[i](0).ch);
      if (plib)
      {  if ((plib->OpticName == "XTARGET") ||
             //(plib->OpticName == "RFCAVITY")||
             (plib->OpticName == "XECOOL")  ||
             (plib->OpticName == "XLCOOL")  )
         {  opt.LAdd();
            opt[index2] = plib;
            opt[index2]->Length = opt.Length;
            index2++;
            if (plib->OpticName == "XTARGET")
               lib1 = true;
         }
      }
      if(lib1) index1++;


      plib = GetOptics(MadLine[i](1).ch);
      if (plib)
      {  opt.LAdd();
         opt[index2] = plib;
         opt[index2]->Length = opt.Length;

         if (plib->OpticName ==  "SBEND")
         {  xBend* bend = (xBend*) plib;
            for (int j = 2; j < MadLine.Col(i); j++)
            {  if (MadLine[i](j).ch == "ANGLE")
                  bend->ANGLE = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "E1")
                  bend->E1 = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "E2")
                  bend->E2  = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "K1")
                  bend->K1  = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "HGAP")
               {  bend->HGAP1  = MadLine[i][j+1];
                  bend->HGAP2  = MadLine[i][j+1];
               }else
               if (MadLine[i](j).ch == "FINT")
               {  bend->FINT1  = MadLine[i][j+1];
                  bend->FINT2  = MadLine[i][j+1];
               }
            }
         }else

         if (plib->OpticName ==  "EBEND")
         {  xEBend* ebend = (xEBend*) plib;
            for (int j = 2; j < MadLine.Col(i); j++)
            {  if (MadLine[i](j).ch == "RO")
                  ebend->Ro  = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "FB")
                  ebend->FB = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "FE")
                  ebend->FE = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "N")
                  ebend->n    = MadLine[i][j+1];
            }
         }else

         if (plib->OpticName ==  "FFIELD")
         {  xEBend* ebend = (xEBend*) plib;
            for (int j = 2; j < MadLine.Col(i); j++)
            {  if (MadLine[i](j).ch == "RO")
                  ebend->Ro  = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "FB")
                  ebend->FB = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "FE")
                  ebend->FE = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "N")
                  ebend->n    = MadLine[i][j+1];
            }
         }else

         if (plib->OpticName == "QUADRUPOLE")
         {  xQuadrupole* quad = (xQuadrupole*) plib;
            for (int j = 2; j < MadLine.Col(i); j++)
            {  if (MadLine[i](j).ch == "K1")
                  quad->K1 = MadLine[i][j+1];
               if (MadLine[i](j).ch == "TILT")
                  quad->TILT = MadLine[i][j+1];
            }
         }else

         if (plib->OpticName == "SEXTUPOLE")
         {  xSextupole* sext = (xSextupole*) plib;
            for (int j = 2; j < MadLine.Col(i); j++)
            {  if (MadLine[i](j).ch == "K2")
                  sext->K2 = MadLine[i][j+1];
            }
         }
         else
         if (plib->OpticName == "SOLENOID")
         {  xSolenoid* sol = (xSolenoid*) plib;
            for (int j = 2; j < MadLine.Col(i); j++)
            {  if (MadLine[i](j).ch == "KS")
                  sol->Ks = MadLine[i][j+1];
            }
         }
         else
         if (plib->OpticName == "RFCAVITY")
         {  xRFcavity* rf = (xRFcavity*) plib;
            for (int j = 2; j < MadLine.Col(i); j++)
            {  if (MadLine[i](j).ch == "VOLT")
                  rf->V = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "HARMON")
                  rf->H = MadLine[i][j+1];
            }
         }

//=============== from here elements for CM are =========
/*
         else
         if (MadLine[i](1).ch == "CM_TOROID")
         {  opt.LAdd();
            xCM_Toroid* cm_tor = new xCM_Toroid;
            opt[index2] = cm_tor;
            cm_tor->Length = opt.Length;
            for (int j = 2; j < MadLine.Col(i); j++)
            {
               if (MadLine[i](j).ch == "B")
                  cm_tor->B = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "R")
                  cm_tor->R = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "K1")
                  cm_tor->Kdis = MadLine[i][j+1];
            }
         }
         else
         if (MadLine[i](1).ch == "CM_SOLENOID")
         {  opt.LAdd();
            xCM_Solenoid* cm_sol = new xCM_Solenoid;
            opt[index2] = cm_sol;
            cm_sol->Length = opt.Length;
            for (int j = 2; j < MadLine.Col(i); j++)
            {
               if (MadLine[i](j).ch == "B")
                  cm_sol->B = MadLine[i][j+1];
            }
         }
         else
         if (MadLine[i](1).ch == "CM_SOL_QUAD")
         {  opt.LAdd();
            xCM_Sol_Quad* cm_sol_q = new xCM_Sol_Quad;
            opt[index2] = cm_sol_q;
            cm_sol_q->Length = opt.Length;
            for (int j = 2; j < MadLine.Col(i); j++)
            {
               if (MadLine[i](j).ch == "B")
                  cm_sol_q->B = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "W")
                  cm_sol_q->Wind = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "K1")
                  cm_sol_q->K1 = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "K2")
                  cm_sol_q->K2 = MadLine[i][j+1];
            }
         }
         else
         if (MadLine[i](1).ch == "CM_TOR_QUAD")
         {  opt.LAdd();
            xCM_Tor_Quad* cm_tor_q = new xCM_Tor_Quad;
            opt[index2] = cm_tor_q;
            cm_tor_q->Length = opt.Length;
            for (int j = 2; j < MadLine.Col(i); j++)
            {
               if (MadLine[i](j).ch == "B")
                  cm_tor_q->B = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "W")
                  cm_tor_q->Wind = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "K1")
                  cm_tor_q->K1 = MadLine[i][j+1];
               else
               if (MadLine[i](j).ch == "K2")
                  cm_tor_q->K2 = MadLine[i][j+1];
            }
         }
//========================================================
*/
      }
   }
   LSetSize(index1);
   return false;
}

bool xRing::OutputMAD(char* FileName)
{
	BData MADFile, ColFile;
   int i, tag = 15;
   ColFile.Size(2, 15);
   int tab = (int) Data[tag][34];

   switch(tab)
   {
    case 0: //132 columns
      ColFile[0][0] = 19; ColFile[1][0] = 28;
      ColFile[0][1] = 29; ColFile[1][1] = 37;
      ColFile[0][2] = 38; ColFile[1][2] = 44;
      ColFile[0][3] = 45; ColFile[1][3] = 52;
      ColFile[0][4] = 53; ColFile[1][4] = 59;
      ColFile[0][5] = 60; ColFile[1][5] = 66;
      ColFile[0][6] = 67; ColFile[1][6] = 73;
      ColFile[0][7] = 74; ColFile[1][7] = 79;
      for(i = 1; i < 8; i++)
      {  ColFile[0][i+7] = ColFile[0][i] + 51;
         ColFile[1][i+7] = ColFile[1][i] + 51;
      }
      break;
    case 1: //157 columns
      ColFile[0][0] = 19; ColFile[1][0] = 28;
      ColFile[0][1] = 52; ColFile[1][1] = 61;
      ColFile[0][2] = 62; ColFile[1][2] = 71;
      ColFile[0][3] = 72; ColFile[1][3] = 79;
      ColFile[0][4] = 80; ColFile[1][4] = 85;
      ColFile[0][5] = 86; ColFile[1][5] = 91;
      ColFile[0][6] = 92; ColFile[1][6] = 97;
      ColFile[0][7] = 98; ColFile[1][7] = 103;
      for(i = 1; i < 8; i++)
      {  ColFile[0][i+7] = ColFile[0][i] + 52;
         ColFile[1][i+7] = ColFile[1][i] + 52;
      }
      break;
    case 2: case 4: //user columns
      for (i = 0; i < 15; i++)
      {  ColFile[0][i] = Data[tag][i*2+1]-1;
         ColFile[1][i] = Data[tag][i*2+2]-1;
      }
      break;
    case 3: //MAD (134 columns)
      ColFile[0][0] = 19; ColFile[1][0] = 28;
      ColFile[0][1] = 29; ColFile[1][1] = 37;
      ColFile[0][2] = 38; ColFile[1][2] = 44;
      ColFile[0][3] = 45; ColFile[1][3] = 52;
      ColFile[0][4] = 53; ColFile[1][4] = 60;
      ColFile[0][5] = 61; ColFile[1][5] = 67;
      ColFile[0][6] = 68; ColFile[1][6] = 74;
      ColFile[0][7] = 75; ColFile[1][7] = 80;
      for(i = 1; i < 8; i++)
      {  ColFile[0][i+7] = ColFile[0][i] + 52;
         ColFile[1][i+7] = ColFile[1][i] + 52;
      }
      break;
   }
   /*
   for (i = 0; i < 15; i++)
   {  Data[tag][i*2+1] = ColFile[0][i]+1;
      Data[tag][i*2+2] = ColFile[1][i]+1;
   }
   */

   if (tab == 4)
   {  if (!MADFile.Load(FileName)) return true;
      int index1 = 0;
      Arc = 0;
      for (i = 0; i < MADFile.Row()-1; i++)
      {  if(LGetSize() <= index1)
            LSetSize(index1*2);
         xRingOptics& opt = *LList[index1];

         opt.Length        = MADFile[i+1][0]-MADFile[i][0];
         opt.Lattice.dist  = MADFile[i][0];
         opt.Lattice.betax = MADFile[i][1];
         opt.Lattice.alfax = MADFile[i][2];
         opt.Lattice.mux   = MADFile[i][3];
         opt.Lattice.x     = MADFile[i][4];
         opt.Lattice.px    = MADFile[i][5];
         opt.Lattice.Dx    = MADFile[i][6];
         opt.Lattice.Dpx   = MADFile[i][7];
         opt.Lattice.betay = MADFile[i][8];
         opt.Lattice.alfay = MADFile[i][9];
         opt.Lattice.muy   = MADFile[i][10];
         opt.Lattice.y     = MADFile[i][11];
         opt.Lattice.py    = MADFile[i][12];
         opt.Lattice.Dy    = MADFile[i][13];
         opt.Lattice.Dpy   = MADFile[i][14];
         opt.EL_LATTICE    = true;

         Arc += opt.Length;
         opt.LSetSize(0);
         index1++;
      }
   }else
   { if (!MADFile.Load(FileName, &ColFile)) return true;

   int index1 = 0;
   Arc = 0;
   bool errors1;
   for (i = 0; i < MADFile.Row()-1; i++)
   {  errors1 = false;
      for (int k = 0; k < MADFile.Col(i); k++)
      {  if (MADFile[i][k] == EmptyData)
         {  errors1 = true;
            break;
         }
      }
      if (errors1)
      {  MADFile.LRemove(i);
         i--;
      }else
      {
         for (int k = 0; k < MADFile.Col(i+1); k++)
         {  if (MADFile[i+1][k] == EmptyData)
            {  errors1 = true;
               break;
            }
         }
         if (errors1)
         {  MADFile.LRemove(i+1);
            i--;
         }
      }
      if (!errors1)
      {  if(LGetSize() <= index1)
            LSetSize(index1*2);
         xRingOptics& opt = *LList[index1];

         opt.Length        = MADFile[i+1][0]-MADFile[i][0];
         opt.Lattice.dist  = MADFile[i][0];
         opt.Lattice.betax = MADFile[i][1];
         opt.Lattice.alfax = MADFile[i][2];
         opt.Lattice.mux   = MADFile[i][3];
         opt.Lattice.x     = MADFile[i][4];
         opt.Lattice.px    = MADFile[i][5];
         opt.Lattice.Dx    = MADFile[i][6];
         opt.Lattice.Dpx   = MADFile[i][7];
         opt.Lattice.betay = MADFile[i][8];
         opt.Lattice.alfay = MADFile[i][9];
         opt.Lattice.muy   = MADFile[i][10];
         opt.Lattice.y     = MADFile[i][11];
         opt.Lattice.py    = MADFile[i][12];
         opt.Lattice.Dy    = MADFile[i][13];
         opt.Lattice.Dpy   = MADFile[i][14];
         opt.EL_LATTICE    = true;

         Arc += opt.Length;
         opt.LSetSize(0);
         index1++;
      }
   }
   LSetSize(index1);
   }
   //MADFile.Save ("Lattice.tmp");
   if (Arc() == 0)
      Warning("Error in Lattice Structure file :",PressEnter,FileName);
   return false;
}

vectorU xRing::GetForces(xTime&time, vectorU&ion)
{  vectorU forces(U1_, m_^-1);
	xRingOptics& opt = *LList[Index];
  	iDrift.Length = opt.Length;
   forces = iDrift.GetForces(time, Energy, ion);

   for (int j = 0; j < opt.LGetSize(); j++)
   {  vectorU f;
      if(opt[j]->EL_FORCES)
         f = opt[j]->GetForces(time, Energy, ion);
      else
      if(opt[j]->EL_FIELDS)
      {  f = opt[j]->GetFields(time, Energy, ion);
         f = Field2Force(ion, f);
      }
      else
         Warning("Optics Library Forces (index=",j,") was not defined");
      forces += f;
   }
   return forces;
}

vectorU xRing::Field2Force(vectorU&X, vectorU&F)
{  vectorU forces(U1_, m_^-1);
   doubleU Ze_Pc(Energy.Z * U_e / (Energy.Momentum * U_c));

	for (int i = 0; i < 3; i++)
	forces[i*2+1] = Ze_Pc * F[i*2] / Energy.Beta;
	forces[1]    += Ze_Pc * ( (X[3]   *F[5]) -((X[5]+1)*F[3]));
	forces[3]    += Ze_Pc * (((X[5]+1)*F[1]) - (X[1]   *F[5]));
	forces[5]    += Ze_Pc * ( (X[1]   *F[2]) - (X[3]   *F[1]));

   return forces;
}

xLattice& xRing::GetLattice(int n)
{	if ((n >= 0) && (n < LGetSize()))
   {
   	if (LList[n]->EL_LATTICE)
         return LList[n]->Lattice;
      else
         Warning("Lattice (index=",n,") was not defined");
	}
   Warning("Ring has not index : ", n);
   return (*LList[0]).Lattice;
}

matrixU& xRing::GetMatrix(xTime& time, int n)
{
   if (n >= LGetSize())
      Warning("Ring has not index : ", n);
   return (*LList[n]).GetMatrix(time, Energy);
}

matrixU xRing::RingMatrix(xTime&t, int n)
{
   Matrix.Diag(1,0);
   {  if ((n >= 0) && (n < Number()))     
      {  for (int i = n; i < Number(); i++)
            Matrix = (*LList[i]).GetMatrix(t, Energy) * Matrix;
         for (int j = 0; j < n; j++)
            Matrix = (*LList[j]).GetMatrix(t, Energy) * Matrix;
      }else
         Warning("Ring has not index : ", n);
   }
   return Matrix;
}

vectorU xRing::OnEnter(xTime&t,vectorU ion)
{
	for (int i = 0; i < LList[Index]->LGetSize(); i++)
		ion = (*LList[Index])[i]->OnEnter(t, Energy, ion);
   return ion;
}

vectorU xRing::OnExit (xTime&t,vectorU ion)
{
	for (int i = 0; i < LList[Index]->LGetSize(); i++)
		ion = (*LList[Index])[i]->OnExit(t, Energy, ion);
   return ion;      
}

//********** returns Lattices from Transformation Matrix ********

xLattice xRing::CalcLattice(matrixU& TrMatr)
{
   xLattice L;
   double alpha, beta, gamma, mu;
   double sin_mu, cos_mu;
   double sign = 1.;
   int size = 2;

   matrixU S(4,4);

  for (int i = 0; i<4; i++)
   for (int j = 0; j<4; j++)
    S[i][j] = TrMatr[i][j];

     double alfa_func[2];
     double beta_func[2];
     double gamma_func[2];
     double mu_func[2];
     double Q_func[2];

      for (int i1 = 0; i1 < size; i1++)
      {
       alfa_func[i1] = 0.;
       beta_func[i1] = 0.;
       gamma_func[i1] = 0.;
       mu_func[i1] = 0.;
       Q_func[i1] = 0.;
      }

   int a = 0;
   int b = 1;

	for (int i2 = 0; i2 < size; i2++)
		{
         a = 2*i2;
         b = (2*i2)+1;

         cos_mu = (S(a,a)()+ S(b,b)())/2.;

         if (S(a,b)() >= 0) sign = 1.;
//         else if (S(a,b).v == 0) sign = 0.;
         else sign = -1.;

//           sin_mu = sign*sqrt(fabs ((S(a,b)*S(b,a))()*(-1.)- ((S(a,a)-S(b,b))*(S(a,a)-S(b,b)))()/4. ) );
           double R1 = (-1.)*(S(a,b)*S(b,a))();
           double R2 = ((S(a,a)-S(b,b))*(S(a,a)-S(b,b)))()/4.;
           sin_mu = sign*sqrt(fabs(R1-R2));

         //double Q = atan(sin_mu/cos_mu) / (2*M_PI);
         if (sin_mu == 0)
            {
              mu = asin(sin_mu);
              alpha = 1;
              beta = 1;
              gamma = 1;
             }
         else
             {
              mu = asin(sin_mu);
              //mu = acos(cos_mu);
              alpha = (S(a,a)-S(b,b))()/(2.*sin_mu);
              beta = S(a,b)()/sin_mu;
              gamma = (-1.)* S(b,a)()/sin_mu;
             }
         mu_func[i2] = mu;
         alfa_func[i2] = alpha;
         beta_func[i2] = beta;
         gamma_func[i2] = gamma;
         Q_func[i2] = atan(sin_mu/cos_mu) / (2*M_PI);
         if(Q_func[i2] < 0) Q_func[i2] += 0.5;
	}

   L.alfax = alfa_func[0];
   L.alfay = alfa_func[1];
   L.betax = beta_func[0];
   L.betay = beta_func[1];
   L.gammax = gamma_func[0];
   L.gammay = gamma_func[1];
   L.mux = mu_func[0];
   L.muy = mu_func[1];
   L.Qx = Q_func[0];
   L.Qy = Q_func[1];

   L.Dx =  ( ( (1.-TrMatr(1,1))*TrMatr(0,5)) + (TrMatr(0,1)*TrMatr(1,5)))/(2.-TrMatr(0,0)-TrMatr(1,1));
   L.Dpx = ( ( (1.-TrMatr(0,0))*TrMatr(1,5)) + (TrMatr(1,0)*TrMatr(0,5)))/(2.-TrMatr(0,0)-TrMatr(1,1));
   L.Dy =  ( ( (1.-TrMatr(3,3))*TrMatr(2,5)) + (TrMatr(2,3)*TrMatr(3,5)))/(2.-TrMatr(2,2)-TrMatr(3,3));
   L.Dpy = ( ( (1.-TrMatr(2,2))*TrMatr(3,5)) + (TrMatr(3,2)*TrMatr(2,5)))/(2.-TrMatr(2,2)-TrMatr(3,3));

   return L;
}

xLattice xRing::CalcLattice(matrixU& TM, xLattice& L1)
{
   xLattice L2;

   L2.betax = (1./L1.betax)*((((TM(0,0)*L1.betax)-(TM(0,1)*L1.alfax))*((TM(0,0)*L1.betax)-(TM(0,1)*L1.alfax)))+(TM(0,1)*TM(0,1)));
   L2.alfax = (-1./L1.betax)*((((TM(0,0)*L1.betax)-(TM(0,1)*L1.alfax))*((TM(1,0)*L1.betax)-(TM(1,1)*L1.alfax)))+(TM(0,1)*TM(1,1)));
   L2.gammax = (L1.betax*(TM(1,0)*TM(1,0))) - (2.*L1.alfax*(TM(1,1)*TM(1,0))) + (L1.gammax*(TM(1,1)*TM(1,1)));
   L2.mux = L1.mux + U_Atan(TM(0,1)/((TM(0,0)*L1.betax)-(TM(0,1)*L1.alfax)) );

   L2.betay = (1./L1.betay)*((((TM(2,2)*L1.betay)-(TM(2,3)*L1.alfay))*((TM(2,2)*L1.betay)-(TM(2,3)*L1.alfay)))+(TM(2,3)*TM(2,3)));
   L2.alfay = (-1./L1.betay)*((((TM(2,2)*L1.betay)-(TM(2,3)*L1.alfay))*((TM(3,2)*L1.betay)-(TM(3,3)*L1.alfay)))+(TM(2,3)*TM(3,3)));
   L2.gammay = (L1.betay*(TM(3,2)*TM(3,2))) - (2.*L1.alfay*(TM(3,3)*TM(3,2))) + (L1.gammay*(TM(3,3)*TM(3,3)));
   L2.muy = L1.muy + U_Atan(TM(2,3)/((TM(2,2)*L1.betay)-(TM(2,3)*L1.alfay)) );

   L2.Dx = (TM(0,0)*L1.Dx) + (TM(0,1)*L1.Dpx) + TM(0,5);
   L2.Dpx = (TM(1,0)*L1.Dx) + (TM(1,1)*L1.Dpx) + TM(1,5);
   L2.Dy = (TM(2,2)*L1.Dy) + (TM(2,3)*L1.Dpy) + TM(2,5);
   L2.Dpy = (TM(3,2)*L1.Dy) + (TM(3,3)*L1.Dpy) + TM(3,5);

   return L2;
}

matrixU xRing::CalcTransformMatrix (xLattice& L1, xLattice& L2)
{
   matrixU TM;
   doubleU mu_x;
   doubleU mu_y;

   mu_x = L2.mux - L1.mux;
   mu_y = L2.muy - L1.muy;

   TM(0,0, ((L2.betax/L1.betax)^0.5)*(U_Cos(mu_x) + (L1.alfax*U_Sin(mu_x))));
   TM(0,1, ((L1.betax*L2.betax)^0.5)*U_Sin(mu_x));
   TM(1,0, (-1.)*((L2.alfax-L1.alfax)*U_Cos(mu_x) + ((1.+(L2.alfax*L1.alfax))*U_Sin(mu_x)) )/((L1.betax*L2.betax)^0.5));
   TM(1,1, ((L1.betax/L2.betax)^0.5)*(U_Cos(mu_x)-(L2.alfax*U_Sin(mu_x))));

   TM(2,2, ((L2.betay/L1.betay)^0.5)*(U_Cos(mu_y) + (L1.alfay*U_Sin(mu_y))));
   TM(2,3, ((L1.betay*L2.betay)^0.5)*U_Sin(mu_y));
   TM(3,2, (-1.)*( (L2.alfay - L1.alfay)*U_Cos(mu_y) + ((1.+(L2.alfay*L1.alfay))*U_Sin(mu_y)) )/((L1.betay*L2.betay)^0.5));
   TM(3,3, ((L1.betay/L2.betay)^0.5)*(U_Cos(mu_y)-(L2.alfay*U_Sin(mu_y))));

  return TM;
}

void xRing::SetLattice()
{  xTime t(*this);
   matrixU matrix;
   matrix = RingMatrix(t,0);
   (*LList[0]).Lattice = CalcLattice(matrix);
   (*LList[0]).EL_LATTICE = true;
   for (int i = 1; i < Number(); i++)
   {  xLattice& lat = (*LList[i]).Lattice;
      lat = CalcLattice((*LList[i-1]).GetMatrix(t, Energy),
                        (*LList[i-1]).Lattice);
      (*LList[i]).EL_LATTICE = true;
   }
}



