//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xDistributor.h"
#define PERCENTSIGMA 1
#include <iostream>
#include <omp.h>
#include <algorithm>

//---------------------------------------------------------------------------
void xBunch::GetMean(vectorU ion, Tray<double>&bunch)
{
   for (int i = 0; i < bunch.Number; i++)               // Square distances array
   {  Dist2[i] = 0;
      for (int j = 0; j < 6; j+=2)
         Dist2[i] +=(bunch[i][j]-ion[j]())*(bunch[i][j]-ion[j]());
   }

   double max, min = 0;                                 // Searching for min
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

   for (int i = 0; i < 6; i++)                          // Reset mean arrays
      Mean[i] = 0;
   for (int i = 0; i < 6; i++)
      Mean2[i] = 0;

   for (int j = 0; j < 6; j++)                          // Find mean array
   {  for (int i = 0; i < Minst.Number; i++)
         Mean[j] += bunch[Minst[i]][j];
      Mean[j] /= Minst.Number;
   }

   for (int j = 0; j < 6; j++)                          // Find mean square array
   {  for (int i = 0; i < Minst.Number; i++)
         Mean2[j] += (bunch[Minst[i]][j]-Mean[j]) * (bunch[Minst[i]][j]-Mean[j]);
      Mean2[j] /= Minst.Number;
   }
   /*
   Mean2[0] = 0;
   Mean2[2] = 0;
   Mean2[4] = 0;
   
   for (int i = 0; i < Minst.Number; i++)
       {
         Mean2[0] += (bunch[Minst[i]][0]-ion[0]()) * (bunch[Minst[i]][0]-ion[0]());
         Mean2[2] += (bunch[Minst[i]][2]-ion[2]()) * (bunch[Minst[i]][2]-ion[2]());
         Mean2[4] += (bunch[Minst[i]][4]-ion[4]()) * (bunch[Minst[i]][4]-ion[4]());
       }
      Mean2[0] /= Minst.Number;
      Mean2[2] /= Minst.Number;
      Mean2[4] /= Minst.Number;
   */
   /*
   for (int i = 0; i < Minst.Number; i++)
       {
         if(Mean2[0] < ((bunch[Minst[i]][0]-Mean[0]) * (bunch[Minst[i]][0]-Mean[0])))
            Mean2[0] = (bunch[Minst[i]][0]-Mean[0]) * (bunch[Minst[i]][0]-Mean[0]);
         if(Mean2[2] < ((bunch[Minst[i]][2]-Mean[2]) * (bunch[Minst[i]][2]-Mean[2])))
            Mean2[2] = (bunch[Minst[i]][2]-Mean[2]) * (bunch[Minst[i]][2]-Mean[2]);
         if(Mean2[4] < ((bunch[Minst[i]][4]-Mean[4]) * (bunch[Minst[i]][4]-Mean[4])))
            Mean2[4] = (bunch[Minst[i]][4]-Mean[4]) * (bunch[Minst[i]][4]-Mean[4]);
       }
  */
}

void xBunch::LRF_PRF(U_Energy& e, Tray<double>&local)
{
   for (int i = 0; i < local.Number; i++)
   {
      local[i][1] *= (U_c * e.Beta * e.Gamma)(m_/s_);
      local[i][3] *= (U_c * e.Beta * e.Gamma)(m_/s_);
      local[i][4] *= (e.Gamma)();
      local[i][5] *= (U_c * e.Beta)(m_/s_);
   }
}

void xBunch::PRF_LRF(U_Energy& e, Tray<double>&local)
{
   for (int i = 0; i < local.Number; i++)
   {
      local[i][1] /= (U_c * e.Beta * e.Gamma)(m_/s_);
      local[i][3] /= (U_c * e.Beta * e.Gamma)(m_/s_);
      local[i][4] /= (e.Gamma)();
      local[i][5] /= (U_c * e.Beta)(m_/s_);
   }
}


//---------------------------------------------------------------------------
xDistributor::xDistributor(int n)
{
   time_t t;
   srand((unsigned) time(&t)); //initializes the random number generator

   Bunch.SetNumber(n,6);
   Invs.SetNumber(n,6);

   PRF = false;

   benum = COASTING;
   initial = 0;
   dpshift = 0;
   horshift = 0;
   vershift = 0;
   longshift = 0;
   UseForIBS = true;;
   MeanPercents = false;
                                  
   for (int i = 0; i < 3; i++)
   {  Sigma  [i*2]._(m_);
      Sigma_2[i*2]._(m_^2);
   }
   Sigma_s._(m_);
   Emit[0]._(m_);
   Emit[1]._(m_);
   Emit[4]._(m_);
   //Emit[4]._(eV_ * s_);
   X_FWHM = 0;
   Y_FWHM = 0;
   Array_FWHM = 0;
   xi.SetNumber(5,5);
}

void xDistributor::Number(int n)
{
   Bunch.SetNumber(n,6);
   Invs.SetNumber(n,6);
}

bool xDistributor::index(int i)
{
        if ((i >= 0) && (i < Bunch.Number)) return true;
        Warning("xDistributor has not index [", i,"]");
        return false;
}

double* xDistributor::operator[] (int i)
{  if (index(i)) return Bunch[i];
   return &OutRangeData;
}

vectorU xDistributor::operator() (int i)
{  vectorU v(m_,U1_);
   if (index(i))
      for (int j = 0; j < 6; j++)
         v[j] = Bunch[i][j];
   return v;
}

void xDistributor::operator() (int i, vectorU v)
{  if (index(i))
      for (int j = 0; j < 6; j++)
      {  if (j%2) Bunch[i][j] = v[j](U1_);
         else     Bunch[i][j] = v[j](m_);
      }
}

void xDistributor::operator() (int i, matrixU& m)
{  if (index(i))
   {  double ion[6];
      for (int j0 = 0; j0 < m.Row; j0++)
      {  ion[j0] = 0;
         for (int k = 0; k < m.Col; k++)
            ion[j0] += Bunch[i][k] * real(m[j0][k]);
      }
      for (int j1 = 0; j1 < m.Row; j1++)
         Bunch[i][j1] = ion[j1];
   }
}

void xDistributor::Matrix_2x2(int i, matrixU& m)
{  if (index(i))
   {  double ion[4];
      ion[0] = real(m[0][0])*Bunch[i][0] + real(m[0][1])*Bunch[i][1];
      ion[1] = real(m[1][0])*Bunch[i][0] + real(m[1][1])*Bunch[i][1];
      ion[2] = real(m[2][2])*Bunch[i][2] + real(m[2][3])*Bunch[i][3];
      ion[3] = real(m[3][2])*Bunch[i][2] + real(m[3][3])*Bunch[i][3];
      for (int j = 0; j < 4; j++)
         Bunch[i][j] = ion[j];
   }
}

doubleU xDistributor::operator() (int i, int j)
{  doubleU d;
   if (j%2)
   {  if (PRF) d._(m_/s_);
      else     d._(U1_);
   }else       d._(m_);
   if (index(i))
      d = Bunch[i][j];
   return d;
}

void xDistributor::operator() (int i, int j, doubleU d)
{  if (index(i))
   {  if (j%2)
      {  if (PRF) Bunch[i][j] = d(m_/s_);
         else     Bunch[i][j] = d(U1_);
      }else       Bunch[i][j] = d(m_);
   }
}
/*
void xDistributor::PRF2LRF(U_Energy& e)
{  if (PRF)
   {  PRF = false;
      for (int i = 0; i < Bunch.Number; i++)
      {  Bunch[i][1] /= (U_c * e.Beta * e.Gamma)();
         Bunch[i][3] /= (U_c * e.Beta * e.Gamma)();
         Bunch[i][4] /= (e.Gamma)();
         Bunch[i][5] /= (U_c * e.Beta)();
      }
   }
}

void xDistributor::LRF2PRF(U_Energy& e)
{  if (!PRF)
   {  PRF = true;
      for (int i = 0; i < Bunch.Number; i++)
      {  Bunch[i][1] *= (U_c * e.Beta * e.Gamma)();
         Bunch[i][3] *= (U_c * e.Beta * e.Gamma)();
         Bunch[i][4] *= (e.Gamma)();
         Bunch[i][5] *= (U_c * e.Beta)();
      }
   }
}
*/
void xDistributor::Copy(int i, int j)
{  if (index(i) && index(j))
      for (int k = 0; k < 6; k++)
         Bunch[i][k] = Bunch[j][k];
}

void xDistributor::Add (int i, int j)
{  if (index(i) && index(j))
      for (int k = 0; k < 6; k++)
         Bunch[i][k] += Bunch[j][k];
}

void xDistributor::Add (int i, vectorU v)
{  if (index(i))
      for (int j = 0; j < 6; j++)
      {  if (j%2) Bunch[i][j] += v[j](U1_);
         else     Bunch[i][j] += v[j](m_);
      }
}

void xDistributor::Add (int i, int j, doubleU d)
{  if (index(i))
   {  if (j%2) Bunch[i][j] += d(U1_);
      else     Bunch[i][j] += d(m_);
   }
}

//------------------------------------------------------

double* xDistributor::InvI(int i)
{  if (index(i)) return Invs[i];
   return &OutRangeData;
}

vectorU xDistributor::Inv(int i)
{  vectorU v(m_);
   v[2]._(U1_);
   v[5]._(U1_);
   if (index(i))
      for (int j = 0; j < 6; j++)
         v[j] = Invs[i][j];
   return v;
}

void xDistributor::Inv(int i, vectorU v)
{  if (index(i))
      for (int j = 0; j < 6; j++)
      if ((j == 2) || (j == 5))
         Invs[i][j] = v[j](U1_);
      else
         Invs[i][j] = v[j](m_);
}

doubleU xDistributor::Inv(int i, int j)
{  doubleU d;
   if ((j != 2) && (j != 5))
      d._(m_);
   if (index(i))
      d = Invs[i][j];
   return d;
}

void    xDistributor::Inv(int i, int j, doubleU d)
{  if (index(i))
      if ((j == 2) || (j == 5))
         Invs[i][j] = d(U1_);
      else
         Invs[i][j] = d(m_);
}
//=============================================================================

double xDistributor::Flatten()
{
   return 2. * rand() / RAND_MAX - 1.;
}
//------------------------------------------------------
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)    

double xDistributor::ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
//------------------------------------------------------

double xDistributor::Gaussian()
{
   double rand1 = double(rand()) / RAND_MAX;
   double rand2 = double(rand()) / RAND_MAX;
        if (rand1 < 1e-9) rand1 = 1e-9;
   double g = sqrt(-2.0*log(rand1))*sin(2.0*3.1415926*rand2);
   return g;
}
//------------------------------------------------------
int xDistributor::Poisson(double x)
{
   int i;
   double expx;
   double pir;

   if(x >= 88)
   {
      i = int(x + Gaussian()* sqrt(x));
   }
   else
   {
      expx = exp(-1.*x);
      i = -1;
      pir = 1.0;
      do
      {
         i += 1;
         pir *= double(rand()) / RAND_MAX;
      }while (pir > expx);
   }
   return i;
}
//-----------------------------------------------------------------------
void xDistributor::RMS(vectorU emit, xLattice& lat1)
{
   lat1.gammax = (((lat1.alfax)^2) + 1)/lat1.betax;
   lat1.gammay = (((lat1.alfay)^2) + 1)/lat1.betay;

   Sigma[0] = (emit[0] * lat1.betax )^0.5;
   Sigma[1] = (emit[0] * lat1.gammax)^0.5;
   Sigma[2] = (emit[1] * lat1.betay )^0.5;
   Sigma[3] = (emit[1] * lat1.gammay)^0.5;
   Sigma[5] =  emit[2]               ^0.5;

   if (benum == BUCKET)
   {
      Sigma[4] = CalcBucket(emit[2]);
   }else
   {
      Sigma[4] =  lat1.B_s;
      if (benum == BUNCHED) Sigma[4] *= Sigma[5];
   }
   
   for (int i = 0; i < 6; i++)
      Sigma_2[i] = Sigma[i]^2;
}
//------------------------------------------------------

void xDistributor::SetIon(int i)
{  if (!index(i)) return;

////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07
   for (int j = 0; j < 6; j++)
   {  Bunch[i][j] =  Sigma[j]();
      if (j==4 && benum == COASTING)
         Bunch[i][j] *= Flatten()/2.;
      else
      if (initial == 1)
         Bunch[i][j] *= Flatten() * sqrt(3.);
      else
      if (initial == 3)
         Bunch[i][j] *= Flatten() * 2;
      else
         Bunch[i][j] *= Gaussian();
   }
   if (initial == 3) {
      while ( (Bunch[i][0]*Bunch[i][0]+Bunch[i][2]*Bunch[i][2]) > (Sigma[0]*Sigma[2])()*3 ) 
      {
         Bunch[i][0] = Sigma[0]() * Flatten() * sqrt(3.);
         Bunch[i][2] = Sigma[2]() * Flatten() * sqrt(3.);
      }
   }
   if (benum == BUCKET)
   {
      vectorU X(m_, U1_);
      SetBucket(X);
      Bunch[i][4] = X[4](m_);
      Bunch[i][5] = X[5](U1_);
   }
}
//------------------------------------------------------

void xDistributor::Matching(xLattice& lat2, int i)
{
   lat2.gammax = (((lat2.alfax)^2) + 1) / lat2.betax;
   lat2.gammay = (((lat2.alfay)^2) + 1) / lat2.betay;

   //for (int i = 0; i < number; i++)
   {
      Bunch[i][0] *= ((1/(lat2.betax*lat2.gammax))^0.5)();
      Bunch[i][2] *= ((1/(lat2.betay*lat2.gammay))^0.5)();

      Bunch[i][0] -= (lat2.alfax* lat2.betax)() * Bunch[i][1] /
         (1 + lat2.alfax * lat2.alfax)(); // rotate on x coord by thin lens
      Bunch[i][2] -= (lat2.alfay* lat2.betay)() * Bunch[i][3] /
         (1 + lat2.alfay * lat2.alfay)(); // rotate on y coord by thin lens
   }
}
//------------------------------------------------------

void xDistributor::PlusDisp(xLattice& lat2, int i)
{
   //for (int i = 0; i < number; i++)
   {  Bunch[i][0] += lat2.Dx()  * Bunch[i][5]; // moves beam on x coord accordingly to D
      Bunch[i][1] += lat2.Dpx() * Bunch[i][5]; // moves beam on x' coord accordingly to D'
      Bunch[i][2] += lat2.Dy()  * Bunch[i][5]; // moves beam on y coord accordingly to D
      Bunch[i][3] += lat2.Dpy() * Bunch[i][5]; // moves beam on y' coord accordingly to D'
   }
}
//------------------------------------------------------

void xDistributor::MinusDisp(xLattice& lat2, int i)
{
   //for (int i = 0; i < number; i++)
   {  Bunch[i][0] -= lat2.Dx()  * Bunch[i][5]; // moves beam on x coord accordingly to D
      Bunch[i][1] -= lat2.Dpx() * Bunch[i][5]; // moves beam on x' coord accordingly to D'
      Bunch[i][2] -= lat2.Dy()  * Bunch[i][5]; // moves beam on y coord accordingly to D
      Bunch[i][3] -= lat2.Dpy() * Bunch[i][5]; // moves beam on y' coord accordingly to D'
   }
}
//------------------------------------------------------

void xDistributor::Distribution(vectorU emit, xLattice& lat, int n0)
{
   RMS(emit, lat);
   for (int i = n0; i < Bunch.Number; i++)
   {  SetIon(i);
      Matching(lat,i);
      PlusDisp(lat,i);
   }
}
//------------------------------------------------------
vectorU xDistributor::Emit_Tail(xLattice&lat, Tine<bool>&InCore, int Ncore)
{
   for (int i0 = 0; i0 < Bunch.Number; i0++)
      MinusDisp(lat,i0);

   double MO[6];   // array of mean beam coordinates
   double SIG[6];  // array of RMS params for mean beam coordinates
   double Corr[3]; // array of correlations between coord and angle

   for (int j0 = 0; j0 < 6; j0++)
   {  MO[j0] = 0;
      SIG[j0] = 0;
   }
   for (int j1 = 0; j1 < 3; j1++)
      Corr[j1] = 0;

   for (int i1 = 0; i1 < Bunch.Number; i1++)
   if  (!InCore[i1])
   for (int j2 = 0; j2 < 6; j2++)
      MO[j2] += Bunch[i1][j2] ;
   for (int j3 = 0; j3 < 6; j3++)
      MO[j3] /= Bunch.Number - Ncore;

   for (int i2 = 0; i2 < Bunch.Number; i2++)
   if  (!InCore[i2])
   for (int j4 = 0; j4 < 6; j4++)
      SIG[j4] += pow(Bunch[i2][j4] - MO[j4], 2);
   for (int j5 = 0; j5 < 6; j5++)
      SIG[j5] /= Bunch.Number - Ncore;

   for (int i3 = 0; i3 < Bunch.Number; i3++)
   if  (!InCore[i3])
   for (int j6 = 0; j6 < 3; j6++)
      Corr[j6] += (Bunch[i3][j6*2]-MO[j6*2])*(Bunch[i3][j6*2+1]-MO[j6*2+1]);
   for (int j7 = 0; j7 < 3; j7++)
      Corr[j7] /= Bunch.Number - Ncore;

   vectorU emittance;
   emittance[0]._(m_);
   emittance[1]._(m_);
   emittance[4]._(m_);

   emittance[0] = sqrt((SIG[0] * SIG[1]) - (Corr[0] * Corr[0]));
   emittance[1] = sqrt((SIG[2] * SIG[3]) - (Corr[1] * Corr[1]));
   emittance[2] = SIG[5];
   emittance[3] = Emit[3];
   Sigma_s      = sqrt(SIG[4]);
   emittance[4] = sqrt((SIG[4] * SIG[5]) - (Corr[2] * Corr[2]));;

   for (int i4 = 0; i4 < Bunch.Number; i4++)
      PlusDisp(lat,i4);
   return emittance;
}


vectorU xDistributor::Emit_RMS(xLattice& lat)
{
   for (int i0 = 0; i0 < Bunch.Number; i0++)
      MinusDisp(lat,i0);

   double MO[6];   // array of mean beam coordinates
   double SIG[6];  // array of RMS params for mean beam coordinates
   double Corr[3]; // array of correlations between coord and angle

   for (int j0 = 0; j0 < 6; j0++)
   {  MO[j0] = 0;
      SIG[j0] = 0;
   }
   for (int j1 = 0; j1 < 3; j1++)
      Corr[j1] = 0;

   for (int i1 = 0; i1 < Bunch.Number; i1++)
   for (int j2 = 0; j2 < 6; j2++)
      MO[j2] += Bunch[i1][j2] ;
   for (int j3 = 0; j3 < 6; j3++)
      MO[j3] /= Bunch.Number;

   for (int i2 = 0; i2 < Bunch.Number; i2++)
   for (int j4 = 0; j4 < 6; j4++)
      SIG[j4] += pow(Bunch[i2][j4] - MO[j4], 2);
   for (int j5 = 0; j5 < 6; j5++)
      SIG[j5] /= Bunch.Number;

   for (int i3 = 0; i3 < Bunch.Number; i3++)
   for (int j6 = 0; j6 < 3; j6++)
      Corr[j6] += (Bunch[i3][j6*2]-MO[j6*2])*(Bunch[i3][j6*2+1]-MO[j6*2+1]);
   for (int j7 = 0; j7 < 3; j7++)
      Corr[j7] /= Bunch.Number;

   vectorU emittance;
   emittance[0]._(m_);
   emittance[1]._(m_);
   emittance[4]._(m_);

   emittance[0] = sqrt((SIG[0] * SIG[1]) - (Corr[0] * Corr[0]));
   emittance[1] = sqrt((SIG[2] * SIG[3]) - (Corr[1] * Corr[1]));
   emittance[2] = SIG[5];
   emittance[3] = Emit[3];
   Sigma_s      = sqrt(SIG[4]);
   emittance[4] = sqrt((SIG[4] * SIG[5]) - (Corr[2] * Corr[2]));

   for (int i4 = 0; i4 < Bunch.Number; i4++)
      PlusDisp(lat,i4);
   return emittance;
}
//=============================================================================

void xDistributor::Courant_Snyder(xLattice& lat, int i)
{
   doubleU gammax(m_^-1);
   gammax = (1. + (lat.alfax^2))/lat.betax;
   doubleU gammay(m_^-1);
   gammay = (1. + (lat.alfay^2))/lat.betay;
   double Xb, Xs, Yb, Ys;

   Xb = Bunch[i][0] - lat.Dx() * Bunch[i][5];
   Xs = Bunch[i][1] - lat.Dpx()* Bunch[i][5];
   Yb = Bunch[i][2] - lat.Dy() * Bunch[i][5];
   Ys = Bunch[i][3] - lat.Dpy()* Bunch[i][5];

   Invs[i][0] = ((lat.betax()*Xs*Xs)+(2.*lat.alfax()*Xb*Xs)+(gammax()*Xb*Xb));
   Invs[i][1] = ((lat.betay()*Ys*Ys)+(2.*lat.alfay()*Yb*Ys)+(gammay()*Yb*Yb));
   Invs[i][2] = Bunch[i][5] * Bunch[i][5];
   if (benum == BUNCHED)
      Invs[i][2]+= Bunch[i][4] * Bunch[i][4] / (lat.B_s() * lat.B_s());
   else
      Invs[i][2] *= PERCENTSIGMA;   // ?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*
}

void xDistributor::Courant_Snyder(xLattice& lat)
{
   doubleU gammax(m_^-1);
   gammax = (1. + (lat.alfax^2))/lat.betax;
   doubleU gammay(m_^-1);
   gammay = (1. + (lat.alfay^2))/lat.betay;
   double X[6];

   double dP = 0;
   if (MeanPercents)
   {
      for (int i = 0; i < Bunch.Number; i++)
         dP += Bunch[i][5];
      dP /= Bunch.Number;
   }

   for (int i = 0; i < Bunch.Number; i++)
   {
      X[0] = Bunch[i][0] - lat.Dx() * Bunch[i][5];
      X[1] = Bunch[i][1] - lat.Dpx()* Bunch[i][5];
      X[2] = Bunch[i][2] - lat.Dy() * Bunch[i][5];
      X[3] = Bunch[i][3] - lat.Dpy()* Bunch[i][5];
      X[4] = Bunch[i][4];
      X[5] = Bunch[i][5];

      if (MeanPercents)
         X[5] -= dP;

      Invs[i][0] = ((lat.betax()*X[1]*X[1])+(2.*lat.alfax()*X[0]*X[1])+(gammax()*X[0]*X[0]));
      Invs[i][1] = ((lat.betay()*X[3]*X[3])+(2.*lat.alfay()*X[2]*X[3])+(gammay()*X[2]*X[2]));
      Invs[i][2] = X[5] * X[5];
      if (benum == BUNCHED)
         Invs[i][2]+= X[4] * X[4] / (lat.B_s() * lat.B_s());
      else
         Invs[i][2] *= PERCENTSIGMA; // ?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*
   }
}
//------------------------------------------------------

void xDistributor::QuickSort(Tray<double>&array, int size, int col)
{
   int i, j, k;
   double x;
   for(i=0; i < size; i++)
   {  k = i;
      x = array[i][col];
      for( j = i+1; j < size; j++)
         if (  array[j][col] < x )
         {  k=j;
            x=array[j][col];
         }
      array[k][col] = array[i][col];
      array[i][col] = x;
   }
}

//------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   shen modified     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void xDistributor::SimpleSort(Tray<double>&array, int size, int col)
{	
	//long atime=clock();/////////////////////////////////////////////////////////////////////////////////// 	 
	double *temp = new double[size];
	for (int i=0; i<size; i++)
		temp[i] = array[i][col];  
	sort(temp,temp+size);
	for (int i=0; i<size; i++)
		array[i][col] = temp[i];
	delete []temp;
	//double tmp;
	// for (int i = size-1; i > 0; i--)
	//    for (int j = 0; j < i; j++)
	//       if (array[j][col] > array[j+1][col])
	//       {  tmp = array[j][col];
	//          array[j][col]   = array[j+1][col];
	//          array[j+1][col] = tmp;
	//       }
	//cout<<(double)(clock()-atime)<<"\n";//////////////////////////////////////////////////  
}

//------------------------------------------------------

vectorU xDistributor::EmitSort(double percent, double percentl, xLattice& lat)
{
   int ind = (int)(floor(percent*Bunch.Number/100.)-1);
   int indl = (int)(floor(percentl*Bunch.Number/100.)-1);
   vectorU R_emit = Emit_RMS(lat);
   Courant_Snyder(lat);   
   for (int j = 0; j < 2; j++)
   {
      //QuickSort(Invs, Bunch.Number, j);
      SimpleSort(Invs, Bunch.Number, j);
      R_emit[j] = Invs[ind][j];
   }
   SimpleSort(Invs, Bunch.Number, 2);
   R_emit[2] = Invs[indl][2];   
   return R_emit;
}

//------------------------------------------------------

void xDistributor::Hystogra(double*hyst, int div, int col)
{
   double max = Invs[Bunch.Number-1][col];
   for (int k = 0; k <= div; k++)
      hyst[k] = 0;
   
   for (int i = 0; i < Bunch.Number; i++)
   {
      int index = Round(Invs[i][col] * div / max);
      if (index < 0)   index = 0;
      if (index > div) index = div;
      hyst[index] += 1;
   }
}

//------------------------------------------------------
void xDistributor::Profile(double *hyst, int col, double emit, int div)
{
   double h = 2. * Sigmas / div;
   double hn = h * emit;
   for (int k = 0; k <= div; k++)
      hyst[k] = 0;

   if (col != 5 || benum == BUNCHED)
   {
//#pragma omp parallel for
      for(int i = 0; i < Bunch.Number; i++)
      {
         for(double step=-Invs[i][col]; step < Invs[i][col]+hn/2; step+=hn)
         {
            int index=Round( ( step / emit + Sigmas) / h );
            if (index >= 0 && index <= div)
            hyst[index] += hn / sqrt(pow(Invs[i][col] + hn, 2) - step*step);
         }
      }

      h = 0;
      for (int j = 0; j <= div; j++)
         h += hyst[j];
      for (int l = 0; l <= div; l++)
         hyst[l] *= Bunch.Number / h;
   }
}
//------------------------------------------------------
#include "xPowell.h"
double*hst;     // bi-gaussian distribution
double sgm;     // sigma for Powell method
double dvs;     // divisions Bunch.Number

double OneGaussian(double *p)
{
   double x, sum = 0;
   for (int i = 0; i <= dvs; i++)
   {  x = -sgm + (2 * sgm / dvs) * i;
      sum += pow(hst[i] - p[2] * exp(-0.5 * pow(x / p[1], 2.)), 2.);
   }
   return sum;
}

void xDistributor::Powel1(double *p, double *hys, double sig, double div)
{
   hst = hys;
   sgm = sig;
   dvs = div;

   for (int i = 0; i < 5; i++)
   {  for (int j = 0; j < 5; j++)
      if (i == j) xi[i][j] = 1;
         else     xi[i][j] = 0;
   }
   int n = 2;            // Bunch.Number of fitting parameters
   double ftol = 1e-6;   // tolerance
   int iter = 0;         // Bunch.Number of iteration
   double fret = 0;

   powell(p, xi(), n, ftol, &iter, &fret, OneGaussian);

}

double BiGaussian(double *p)
{
   double x, sum = 0;
   for (int i = 0; i <= dvs; i++)
   {  x = -sgm + (2 * sgm / dvs) * i;
      sum += pow(hst[i] - p[3] * exp(-0.5 * pow(x / p[1], 2.)) -
                          p[4] * exp(-0.5 * pow(x / p[2], 2.)) , 2.);
   }
   return sum;
}

void xDistributor::Powel2(double *p, double *hys, double sig, double div)
{
   hst = hys;
   sgm = sig;
   dvs = div;

   for (int i = 0; i < 5; i++)
   {  for (int j = 0; j < 5; j++)
      if (i == j) xi[i][j] = 1;
         else     xi[i][j] = 0;
   }
   int n = 4;            // Bunch.Number of fitting parameters
   double ftol = 1e-6;   // tolerance
   int iter = 0;         // Bunch.Number of iteration
   double fret = 0;

   powell(p, xi(), n, ftol, &iter, &fret, BiGaussian);

}

//------------------------------------------------------
void xDistributor::Amplitudes(xLattice& L)
{
   for (int i = 0; i < Bunch.Number; i++)
   {
      Invs[i][3] = sqrt(Bunch[i][0]*Bunch[i][0] +
                        Bunch[i][1]*Bunch[i][1]*L.betax()*L.betax());
      Invs[i][4] = sqrt(Bunch[i][2]*Bunch[i][2] +
                        Bunch[i][3]*Bunch[i][3]*L.betay()*L.betay());
      if (benum == BUNCHED)
      Invs[i][5] = sqrt(Bunch[i][5]*Bunch[i][5] +
                        Bunch[i][4]*Bunch[i][4]/(L.B_s()*L.B_s()));
      else
      Invs[i][5] = fabs(Bunch[i][5]);
   }

}

vectorU xDistributor::Emit_FWHM(xLattice& Lattice)
{
   vectorU E_RMS, Result;
   Result[0]._(m_);
   Result[1]._(m_);
   Result[4]._(m_);
   double h = 2 * Sigmas / Division;
   Amplitudes(Lattice);
   E_RMS = Emit_RMS(Lattice);

   double width[3], emit[3];
   emit[0] =((E_RMS[0]*Lattice.betax)^0.5)();
   emit[1] =((E_RMS[1]*Lattice.betay)^0.5)();
   emit[2] = (E_RMS[2]               ^0.5)();

   for (int n = 0; n < 3; n++)
   {
      Profile(Hyst3[n], n+3, emit[n], Division);
      width[n] = Width(Hyst3[n], Division);
   }
   Result[0] = pow(emit[0] * width[0] * h / 2.355, 2) / Lattice.betax();
   Result[1] = pow(emit[1] * width[1] * h / 2.355, 2) / Lattice.betay();
   Result[2] = pow(emit[2] * width[2] * h / 2.355, 2);
   Result[3] = E_RMS[3];
   Result[4] = E_RMS[4];
   return Result;
}


//------------------7------------------------------------
double xDistributor::Width(double* Massive, int div)
{
  int j;
  Array_FWHM = 0;
  double maxX = 0;             // Bunch.Number of particles in cell with maximum
  int Xindex = 0;              // the index of this cell
  int leftX = 0, rightX = 0;   // left&right levels for half maximum

   for (j = 0; j <= div; j++)  // seeking for maximum and the cell index
   if (Massive[j]>maxX)
   {
      maxX = Massive[j];
      Xindex = j;
   }
   maxX = maxX / 2;             // half maximum

   j = Xindex;                  // index of center cell
   while (j > 0)                // seeking for left cell corresponded to half maximum
   {  j--;
      if (Massive[j] < maxX)
      {  leftX = j;
         break;
      }
   }

   j = Xindex;
   while (j < div)              // seeking for right cell corresponded to half maximum
   {  j++;
      if (Massive[j] < maxX)
      { rightX = j;
        break;
      }
   }

   for (int i = leftX; i<= rightX; i++)
//   for (int i = 0; i< Bunch.Number; i++)
     Array_FWHM += (int)Massive[i];

   return rightX - leftX;
}
//------------------------------------------------------

vectorU xDistributor::Emittance_FWHM(xLattice& lat, int DivP )
{
   vectorU R_emit;    // result vector
   R_emit[0]._(m_);
   R_emit[1]._(m_);

  double Ex, Ey, dP_P;
  vectorU Em;

     Em = Emit_RMS(lat);
     Ex   = ((Em[0]*lat.betax)^0.5)();
     Ey   = ((Em[1]*lat.betay)^0.5)();
     dP_P = (Em[2]           ^0.5)();

//  static bool recalc = true;
  /*
  if (recalc)
  {
     Em = Emittance_RMS(lat);
     Ex   = ((Em[0]*lat.betax)^0.5)();
     Ey   = ((Em[1]*lat.betay)^0.5)();
     dP_P = (Em[2]           ^0.5)();
     recalc = false;
   }
   else
   {
     Ex   = ((Emit[0]*lat.betax)^0.5)();
     Ey   = ((Emit[1]*lat.betay)^0.5)();
     dP_P = (Emit[2]           ^0.5)();

     Ex   = Emit[0]();
     Ey   = Emit[1]();
     dP_P = Emit[2]();
*/

  int index;
  double difX, difY, difP;
  double hX, hY, hP;
  double MinX, MaxX, MinY, MaxY, MinP, MaxP;
  MinX = -4;
  MaxX = 4;
  MinY = -4;
  MaxY = 4;
  MinP = -4;
  MaxP = 4;
  X_FWHM = 0;
  Y_FWHM = 0;


  Tine<double>IntX(DivP+1);   // array containing Bunch.Number of particles per division (for X)
  Tine<double>IntY(DivP+1);   // array containing Bunch.Number of particles per division (for Y)
  Tine<double>IntP(DivP+1);   // array containing Bunch.Number of particles per division (for dP)

   for (int k = 0; k <= DivP; k++)
   {
      IntX[k] = 0;
      IntY[k] = 0;
      IntP[k] = 0;
   }

      for (int i = 0; i < Bunch.Number; i++)
        {
// for transverse coordinate X
         index = Round((Bunch[i][0]/Ex-MinX)*DivP/(MaxX-MinX));  // calculates the interval where particle is located (for every particle)
         if (index < 0) index = 0;
         if (index > DivP) index = DivP;
         IntX[index] += 1;                         // adds particle inside the interval
// for transverse coordinate Y
         index = Round((Bunch[i][2]/Ey-MinY)*DivP/(MaxY-MinY));  // calculates the interval where particle is located (for every particle)
         if (index < 0) index = 0;
         if (index > DivP) index = DivP;
         IntY[index] += 1;                         // adds particle inside the interval
// for momentum spread
         index = Round((Bunch[i][5]/dP_P-MinP)*DivP/(MaxP-MinP));  // calculates the interval where particle is located (for every particle)
         if (index < 0) index = 0;
         if (index > DivP) index = DivP;
         IntP[index] += 1;                         // adds particle inside the interval
      }

      hX = (MaxX-MinX)/DivP;
      hY = (MaxY-MinY)/DivP;
      hP = (MaxP-MinP)/DivP;

   difX = Width(IntX(), DivP);                        // width on X half maximum
   R_emit[0] = ((difX*Ex*hX/2.52)*(difX*Ex*hX/2.52))/lat.betax();
   X_FWHM = Array_FWHM;

   difY = Width(IntY(), DivP);                        // width on Y half maximum
   R_emit[1] = ((difY*Ey*hY/2.52)*(difY*Ey*hY/2.52))/lat.betay();
   Y_FWHM = Array_FWHM;

   difP = Width(IntP(), DivP);                        // width on dP half maximum
   R_emit[2] = (difP*hP*dP_P/2.52)*(difP*hP*dP_P/2.52);

   R_emit[3] = Emit[3];

   return R_emit;
}

//------------------------------------------------------

vectorU xDistributor::CS_Coord2Emit(xLattice& lat)
{
  vectorU emit(U1_); // result vector
  emit[0]._(m_);
  emit[1]._(m_);

  double inv0,inv1,inv2;
  inv0 = 0;
  inv1 = 0;
  inv2 = 0;

  Courant_Snyder(lat);

  for (int i = 0; i < Bunch.Number; i++)
  {
   inv0 += Invs[i][0];
   inv1 += Invs[i][1];
   inv2 += Invs[i][2];
  }

   emit[0] = inv0/Bunch.Number/2;
   emit[1] = inv1/Bunch.Number/2;
   emit[2] = inv2/Bunch.Number/2;
   emit[3] = Emit[3];

   return emit;
}
//**************************************************

vectorU xDistributor::Emit_Def(xLattice& lattice)
{
   vectorU Res_E(U1_); // result vector
   Res_E[0]._(m_);
   Res_E[1]._(m_);

   switch(EmitDef)
   {
    case 0: Res_E = Emit_RMS(lattice); break;
    case 1: Res_E = CS_Coord2Emit(lattice); break;
    case 2: Res_E = Emit_FWHM(lattice); break;
    case 3: Res_E = EmitSort(percentage, percentlong, lattice); break;
    default: Warning("Wrong Emittance Definition");
   }

   double dP_MO = 0;
   for (int i = 0; i < Bunch.Number; i++)
      dP_MO += Bunch[i][5];
   dP_MO /= Bunch.Number;
   Res_E[5] = dP_MO;

   return Res_E;
}

//------------------------------------------------------

doubleU xDistributor::LocalDensity(double* X, int Count, xLattice&lat)
{
   if (Count < 2) Count = 2;
   if (Count > Bunch.Number) Count = Bunch.Number;

   Tine<double>Radius(Bunch.Number);
   double hourx = 1;
   double houry = 1;
   for (int i0 = 0; i0 < Bunch.Number; i0++)
   {
      if (HourGlass)
      {
         hourx = (X[4] - Bunch[i0][4]) / lat.betax();
         hourx = sqrt(hourx * hourx + 1);
         houry = (X[4] - Bunch[i0][4]) / lat.betay();
         houry = sqrt(houry * houry + 1);
      }
     else
     {
      hourx = 1;
      houry = 1;
     }
      // ------ 05.06
      hourx  = X[0] + Bunch[i0][0]*hourx;
      hourx *= hourx;
      houry  = X[2] - Bunch[i0][2]*houry;
      houry *= houry;
      Radius[i0] = sqrt(hourx + houry);

   }
   //--------------------05.06
   double radius;

   for (int i = 0; i < Count ; i++)
    for (int j = Bunch.Number-1; j > i; j--)
     if(Radius[j]   < Radius[j-1])
     {  radius      = Radius[j];
        Radius[j]   = Radius[j-1];
        Radius[j-1] = radius;
     }

   radius = 0;

   for (int k = 0; k < Count; k++)
       radius += Radius[k]* Radius[k];
        radius /= Count;
   doubleU r( Count / (2.0*M_PI * radius), m_^-2);
   /*
    for (int k = 0; k < Count; k++)
       if(Radius[k])radius += 1.0/ (Radius[k]* Radius[k]);
    doubleU r( radius / M_PI , m_^-2);
   */
   /*
   for (int k = 0; k < Bunch.Number; k++)
       radius += Radius[k]* Radius[k];
   radius /= Bunch.Number;
   */
   //-----------------
   //doubleU r( Count / (M_PI * Radius[Count-1]* Radius[Count-1]), m_^-2);
   //doubleU r( Bunch.Number / (M_PI * radius) , m_^-2);

   return  r;
}
//Local density for Gaussian beam 05.06
doubleU xDistributor::LocalDensity_R(double* X, double part, xLattice&lat)
{
   doubleU r(m_^-2);

   doubleU sigma_x_G(m_);
   doubleU sigma_y_G(m_);


   int numberI = 0;


   sigma_x_G = (( Emit[0] * lat.betax )^0.5);
   sigma_y_G = ((Emit[1]*lat.betay)^0.5);
   r = 0;

   for (int i0 = 0; i0 < Bunch.Number; i0++)

   if((fabs(X[0] - Bunch[i0][0]) < part*sigma_x_G(m_))&&(fabs(X[2] - Bunch[i0][2]) < part*sigma_y_G(m_)))
   {
   numberI++;
   }

   r = numberI/(4.0*part*part*sigma_x_G*sigma_y_G);

   return  r;
}
//Local density for Gaussian beam 05.06
doubleU xDistributor::LocalDensity_G(double* X, xLattice&lat)
{
   doubleU r(m_^-2);

   doubleU sigma_x_G(m_);
   doubleU sigma_y_G(m_);
   doubleU sigma_s_G(m_);

   doubleU x(X[0], m_);
   doubleU y(X[2], m_);
   doubleU s(X[4], m_);



   sigma_x_G = (( Emit[0] * lat.betax )^0.5);
   sigma_y_G = ((Emit[1]*lat.betay)^0.5);
   r = 0;
   /*
   if((x < 6.0*sigma_x_G)&&(y < 6.0*sigma_y_G))
      {
        r = Bunch.Number/(2.*U_pi*sigma_x_G*sigma_y_G);
        r *= 1./U_Exp(((x/sigma_x_G)^2.)/2.);
        r *= 1./U_Exp(((y/sigma_y_G)^2.)/2.);
      }
      */

   if(benum == BUNCHED)
   {
      if((x < 6.0*sigma_x_G)&&(y < 6.0*sigma_y_G))
      {
        r = Bunch.Number/(2.*U_pi*sigma_x_G*sigma_y_G);
        r *= 1./U_Exp(((x/sigma_x_G)^2.)/2.);
        r *= 1./U_Exp(((y/sigma_y_G)^2.)/2.);
      }
   }else
   {
     sigma_s_G = lat.B_s*(Emit[3]^0.5);
     if((x < 6.0*sigma_x_G)&&(y < 6.0*sigma_y_G)&&(s < 6.0*sigma_s_G))
      {
        r = Bunch.Number/(2.*U_pi*sigma_x_G*sigma_y_G*sqrt(X[1]*X[1]+X[3]*X[3]+1.0));
        r *= 1./U_Exp(((x/sigma_x_G)^2.)/2.);
        r *= 1./U_Exp(((y/sigma_y_G)^2.)/2.);
        r *= 1./U_Exp(((s/sigma_s_G)^2.)/2.);
        r *= U_Exp((((x*X[1]+y*X[3]-s)/sigma_s_G)^2.)/
                  (2.0 *(X[1]*X[1]+X[3]*X[3]+1.0)));
      }
   }
   

   return  r;
}
//---- change not seriosly 05.06
void xDistributor::RadialDensity(double*density, int Div, double hx, double hy)
{
   for (int d = 0; d < Div; d++)
      density[d] = 0.0;

   for (int i = 0; i < Bunch.Number; i++)
    for (int j = 1; j <= Div; j++)
     if ((((Bunch[i][0]/(hx*j))*(Bunch[i][0]/(hx*j))) + ((Bunch[i][2]/(hy*j))*(Bunch[i][2]/(hy*j)))) < 1.)
     {
        density[j-1] += 1.0;
        break;
     }
   for (int n = 0; n < Div; n++)
      density[n] /= M_PI * hx * hy * (2.*n + 1.);
}
// Rebealded 05.06
void xDistributor::ProfileDensity(double* density, int Div, double**profile,
                                  double hx, double hy)
{
    /*
    for (int j = 0; j < Div; j++)
      density[j] = sqrt(profile[0][j+Div] * profile[0][j+Div]) /
                       (M_PI * hx * hy * (2.*j + 1.));
    */
    for (int j = 0; j < Div; j++)
      density[j] = sqrt(profile[0][j+Div] * profile[1][j+Div]) /
                       (M_PI * hx * hy * (2.*j + 1.));
}

int xDistributor::index(int j, int num)
{
   double circ;
   if (benum == BUNCHED)
      circ = iRing.L_s(m_);
   else
      circ = iRing.Circ(m_);
   return Round((Bunch[j][4] + circ/2) / circ * num);
}

void xDistributor::LongProfile()
{
   int ix;
   double Pmax = 0;

   for (int i = 0; i < LProfile.Number; i++)
   {  LProfile[i] = 0;
      VProfile[i] = 0;
      Average[i]  = 0;
      Momentum[i] = 0;
      //RatesDp[i]  = 0;
   }

   for (int j = 0; j < Bunch.Number; j++)
   {  ix = index(j, LProfile.Number);
      if (ix >= 0 && ix < LProfile.Number)
      {  LProfile[ix] += 1;
         Average[ix]  += Bunch[j][5];
      }
      if(Pmax < LProfile[ix])
         Pmax = LProfile[ix];
   }
   for (int i = 0; i < LProfile.Number; i++)
   {  if (LProfile[i])
         Average[i]  /= LProfile[i];
      if(i)
         Integral[i] = Integral[i-1] + LProfile[i];
      else
         Integral[0] = LProfile[0];
   }
   for (int j = 0; j < Bunch.Number; j++)
   {  ix = index(j, LProfile.Number);
      if (ix >= 0 && ix < LProfile.Number)
         Momentum[ix] += pow(Bunch[j][5] - Average[ix], 2);
   }
   Pmax /= Integral[LProfile.Number-1] / double(LProfile.Number) * Bunch.Number;
   for (int i = 0; i < LProfile.Number; i++)
   {  if (LProfile[i])
      {  Momentum[i]  /= LProfile[i];
         LProfile[i] *= double(LProfile.Number) / (double)Bunch.Number;
      }
      Integral[i] *= Pmax;
   }
   for (int i = 0; i < LProfile.Number-1; i++)
      VProfile[i] = LProfile[i] - LProfile[i+1];
}

