//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xRunge.h"
#include "xRing.h"

//---------------------------------------------------------------------------
xTime::xTime(xRing& ring)
{
   pVelocity = &ring.Energy.Velocity;
   pCirc     = &ring.Circ;

	t ._(0, s_);
   dt._(0, s_);
   so._(0, m_);
   sr._(0, m_);
   ds._(0, m_);
}

xTime::xTime(xTime& time)
{
	t ._(s_);
   dt._(s_);
   so._(m_);
   sr._(m_);
   ds._(m_);
   *this = time;
}

int xTime::OnGet()
{
   Turns = 0;
   Steps = 0;
   return 0;
}

void xTime::TimeStep(doubleU h)
{
  	dt = h;
   ds = dt * (*pVelocity);
}

void xTime::DistStep(doubleU h)
{
  	ds = h;
   dt = ds / (*pVelocity);
}

void xTime::operator*=(double h)
{
   dt = dt * h;
   ds = ds * h;
}

void xTime::operator+=(double h)
{
   *this *= h;
   t += dt;
   so+= ds;
   sr+= ds;
   int k = (int)floor(((sr+(ds*1e-9))/ *pCirc)(U1_));
   if (k)
   {  sr = 0;
      Turns += k;
   }
   Steps++;
}

void xTime::operator++()
{
   *this += 1;
}

void xTime::operator=(xTime& t1)
{
	ds = t1.ds;
   so = t1.so;
   sr = t1.sr;
   dt = t1.dt;
   t  = t1.t;
	pVelocity = t1.pVelocity;
   pCirc     = t1.pCirc;
   Turns     = t1.Turns;
   Steps     = t1.Steps;
}
//---------------------------------------------------------------------------

vectorU Kutta4(xTime&t,doubleU&h,vectorU X,f_Methods f1)
{
   vectorU result;
   static vectorU k1, k2, k3, k4;
   xTime t2(t), t3(t);
   t2 += 0.5;
   t3 += 1;
   k1 = f1(t,  X)  * h;
   vectorU koef1;
   koef1 = X + (k1/2.);
   k2 = f1(t2, koef1) * h;
   koef1 = X + (k2/2.);
   k3 = f1(t2, koef1) * h;
   koef1 = X + (k3);
   k4 = f1(t3, koef1) * h;
   result = X + ((k1 + (k2*2) + (k3*2) + k4) / 6);
   return result;
}

vectorU Euler(xTime&t,doubleU&h,vectorU X,f_Methods f1)
{
   vectorU result;
   result = f1(t, X) * h + X;
   return result;
}

vectorU Symple(xTime&t,doubleU&h,vectorU X,f_Methods f1)
{
	vectorU X1, X2;
   X1 = f1(t, X)  * h + X;
   X2 = f1(t, X1) * h + X;
   for (int i = 0; i < Size6 / 2; i++)
   	X1[i*2] = X2[i*2];
	return X1;
}
/*
//******************* from Betatron*************************

//************* Square equation's solution ****************

vectorU xRunge::Squareq(vectorU koeff)
 {
  vectorU Multy;
  complex<double> multy[2];

  complex<double> a, b, c;

  a = koeff[0]();
  b = koeff[1]();
  c = koeff[2]();

  multy[0] = ( -b + sqrt(b*b - 4.*a*c))/(2.*a);
  multy[1] = ( -b - sqrt(b*b - 4.*a*c))/(2.*a);

  for (int i=0; i < 2; i++)
   Multy[i] = multy[i];

  return Multy;
}

//************* Qubic square of complex number *********************************

vectorU xRunge::qubstep (complex<double> a)
{
  double x, y, z, phi;
  vectorU mass;

  x = real(a);
  y = imag(a);
  z = sqrt(x*x + y*y);

  phi = atan(y/x);
  if (x < 0)
  phi+= M_PI;

  mass[0] = exp((1./3.)*log (z)) * (cos(phi/3.) + J*sin(phi/3.) );
  mass[1] = exp((1./3.)*log (z)) * (cos((phi + 2.*M_PI)/3.) + J*sin((phi + 2.*M_PI)/3.));
  mass[2] = exp((1./3.)*log (z)) * (cos((phi + 2.*2.*M_PI)/3.) + J*sin((phi + 2.*2.*M_PI)/3.));

  return mass;
}

//************* Qubic equation's decision ****************

vectorU xRunge::Qubeq (vectorU koeff)
{
  vectorU Multy;
  complex<double> multy[3];

  complex<double> a, b, c;
  complex<double> p, q;
  complex<double> A, B, Q;
  complex<double> alpha;

  vectorU Am, Bm;

  a = koeff[0]();
  b = koeff[1]();
  c = koeff[2]();

  p = ((-1.)*a*a/3.) + b;
  q = 2.*((a*a*a)/27.) - (a*b/3.) + c;

  Q = (p/3.)*(p/3.)*(p/3.) + (q/2.)*(q/2.);
  if (real(Q) < 0)
  {
    alpha = acos(real( (-1.)* q /(2.*sqrt((-1.)*(p*p*p)/27. ) )));
    multy[0] = 2.*sqrt((-1.)*p/3.)*cos(alpha/3.) ;
    multy[1] = (-2.)*sqrt((-1.)*p/3.)*cos(alpha/3. + 2.*M_PI/3.);
    multy[2] = (-2.)*sqrt((-1.)*p/3.)*cos(alpha/3. - 2.*M_PI/3.);
  }
  else
  {
  Am = qubstep((-1.)*(q/2.) + sqrt(Q));
  Bm = qubstep((-1.)*(q/2.) - sqrt(Q));

  int i, j;
  for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        if (U_Abs(Am[i]*Bm[j] + p/3.)<0.000001)
	        goto label1;

 label1: A = Am[i]();
         B = Bm[j]();

  multy[0] = A + B;
  multy[1] = (-1.)*(A+B)/2. + J*((A-B)/2.)*sqrt(3.);
  multy[2] = (-1.)*(A+B)/2. - J*((A-B)/2.)*sqrt(3.);
 }
  multy[0] = multy[0] - a/3.;
  multy[1] = multy[1] - a/3.;
  multy[2] = multy[2] - a/3.;

  for (int i = 0; i < 3; i++)
   Multy[i] = multy[i];

  return Multy;
}

//************* Quadra equation's decision ****************
vectorU xRunge::Quateq (vectorU koeff)
{
  vectorU Multy;
  complex<double> multy[4];

  vectorU temp, t1_i, t2_i, t1, t2;

  complex<double> a, b, c, d;
  complex<double> f, g, h, j;
  complex<double> p, q, r, s, y1; //r = 2pq   s - test

  a = koeff[0]();
  b = koeff[1]();
  c = koeff[2]();
  d = koeff[3]();

    f = (-1.)*b;
    g = a*c - 4.*d;
    h = (-1.)*a*a*d + 4.*b*d - c*c;

    temp[0] = f;
    temp[1] = g;
    temp[2] = h;

    temp = Qubeq(temp);

    y1 = temp[0]();

        p = sqrt((a*a/4.) - b + y1);
        q = sqrt((y1*y1/4.)-d);
        r = (a/2.)*y1 - c;
          if (imag(r)==0 && real(r)<0)
          q=-q;

    t1_i[0] = 1.;
    t1_i[1] = (a/2.) - p;
    t1_i[2] = (y1/2.) - q;

    t2_i[0] = 1.;
    t2_i[1] = (a/2.) + p;
    t2_i[2] = (y1/2.) + q;

    t1 = Squareq(t1_i);
    t2 = Squareq(t2_i);

  multy[0] = t1[0]();
  multy[1] = t1[1]();
  multy[2] = t2[0]();
  multy[3] = t2[1]();

  for (int i = 0; i < 4; i++)
   Multy[i] = multy[i];

 return Multy;
}

//*************** Koefficients of characteristic equation **********************

vectorU xRunge::Koeff (matrixU M4x4)
{
  vectorU Multy;
  complex<double> multy[4];

  matrixU m3x3(3,3);

  matrixU m4x4(4,4);
  for (int i = 0; i<4; i++)
   for (int j = 0; j<4; j++)
    m4x4[i][j] = M4x4[i][j];

  complex<double> B = 0;


	multy[0] = (-1.)*m4x4.Trace();

      for (int i = 0; i < 3; i++)
      {  m3x3 = m4x4.Minor(i, i);
         for (int j = i; j < 3; j++)
            B += ( m3x3.Minor(j, j) ).Det();
      }
   multy[1] = B;
      B = 0;
      for (int i = 0; i < 4; i++)
      {
         B += ( m4x4.Minor(i, i) ).Det();
      }
   multy[2] = -B;

   multy[3] = m4x4.Det();

  for (int i = 0; i < 4; i++)
   Multy[i] = multy[i];

 return Multy;
}
//*********************************************************

bool xRunge::QuickEigen(matrixU M4x4, double accuracy)
{
	matrixU m3x3(3, 3);
   matrixU t4x4(4, 4);

   matrixU m4x4(4,4);
   for (int i = 0; i<4; i++)
    for (int j = 0; j<4; j++)
     m4x4[i][j] = M4x4[i][j];

	int i, j;
	complex<double> com, det, SpA, B, cos_mu;
   double Angle[2];
   const double four = 4;
   const double two = 2;

	SpA = B = 0;

	for (i = 0; i < 4; i++) SpA += m4x4[i][i];
	for (i = 0; i < 3; i++)
	{  m3x3 = m4x4.Minor(i, i);
		for (j = i; j < 3; j++)
			B += ( m3x3.Minor(j, j) ).Det();
	}
	B = sqrt( pow(SpA/four, 2) - (B-two)/four );

	Angle[0] = real(arccos(SpA/four + B));
	Angle[1] = real(arccos(SpA/four - B));

   for (i = 0; i < 2; i++)
   {
      Eigen[i*2] = cos(Angle[i]) + J*(complex<double>)sin(Angle[i]);
      Eigen[i*2+1] = conj((complex<double>)Eigen[i*2]());               // !!!!!!!!!!!!!!!!!!
      t4x4.Diag((complex<double>)Eigen[i*2](),0.0);                     //!!!!!!!!!!!!!!!!!!!
      t4x4 = m4x4 - t4x4;
      if (abs(t4x4.Det()) > accuracy) return false;
   }
   return true;
}

//*********************************************************

void xRunge::EigenValue(matrixU Space)
{
   matrixU m4x4(4,4);
   for (int i = 0; i<4; i++)
    for (int j = 0; j<4; j++)
     m4x4[i][j] = Space[i][j];
   if (!QuickEigen(m4x4)) Eigen = Quateq(Koeff(m4x4));
}

//____________EUGENEVECTOR_______________________________

//**********  EugeneVector mormalisation ****************

void xRunge::NormVector()
{
   complex<double> c[4];

   matrixU S(4,4), F(4,4);
   matrixU e1(1,4), e2(1,4), e3(1,4), e4(1,4);
   matrixU b1(4,1), b2(4,1), b3(4,1), b4(4,1);
   matrixU s1(1,1), s2(1,1), s3(1,1), s4(1,1);
   matrixU  r1;
	S.Simplectic();

   for (int i = 0; i<4; i++)
   {
      b1[i][0] = conj(EigenV[i][0]);
      b2[i][0] = conj(EigenV[i][1]);
      b3[i][0] = conj(EigenV[i][2]);
      b4[i][0] = conj(EigenV[i][3]);
   }

   for (int i = 0; i<4; i++)
   {
      e1[0][i] = EigenV[i][0];
      e2[0][i] = EigenV[i][1];
      e3[0][i] = EigenV[i][2];
      e4[0][i] = EigenV[i][3];
   }

   e1 = EigenV.Trans();

   s1 = e1*S*b1;
   s2 = e2*S*b2;
   s3 = e3*S*b3;
   s4 = e4*S*b4;

   c[0] = imag(s1[0][0]);
   if (real(c[0]) == 0) c[0]=1.;

   c[1] = imag(s2[0][0]);
   if (real(c[1]) == 0) c[1]=1;

   c[2] = imag(s3[0][0]);
   if (real(c[2]) == 0) c[2]=1;

   c[3] = imag(s4[0][0]);
   if (real(c[3]) == 0) c[3]=1;

   e1 = e1*(-1)*sqrt(2./c[0]);
   e2 = e2*(-1)*sqrt(2./c[1]);
   e3 = e3*(-1)*sqrt(2./c[2]);
   e4 = e4*(-1)*sqrt(2./c[3]);

   for (int i = 0; i<4; i++)
   {
      F[i][0] = e1[0][i];
      F[i][1] = e2[0][i];
      F[i][2] = e3[0][i];
      F[i][3] = e4[0][i];
   }

   for (int i = 0; i<4; i++)
   {
      b1[i][0] = conj(e1[0][i]);
      b2[i][0] = conj(e2[0][i]);
      b3[i][0] = conj(e3[0][i]);
      b4[i][0] = conj(e4[0][i]);
   }

   s1 = e1*S*b1;
   s2 = e2*S*b2;
   s3 = e3*S*b3;
   s4 = e4*S*b4;

   c[0] = s1[0][0];
   c[1] = s2[0][0];
   c[2] = s3[0][0];
   c[3] = s4[0][0];

   int h = 0;
   for (int i = 0; i < 4; i++)
   {
		if ( (imag(c[i]) + 2.0 ) < 1e-8)
      {
         for (int j = 0; j<4; j++)
         	r1[j][h] = F[j][i];

         r1[5][h] = Eigen[i]();
         h++;
      }
   }

   for (int i = 0; i < 4; i++)
   {
   	NormV[i][0] =      r1[i][0];
   	NormV[i][1] = conj(r1[i][0]);
   	NormV[i][2] =      r1[i][1];
   	NormV[i][3] = conj(r1[i][1]);
	}
   NormE[0] =      r1[5][0];
   NormE[1] = conj(r1[5][0]);
   NormE[2] =      r1[5][1];
   NormE[3] = conj(r1[5][1]);
}

//****** EigenVector estimation **********

void xRunge::EigenVector(matrixU Space, bool OnlyOne)
{
   matrixU m4x4(4,4), m3x3(3,3);
   matrixU X(3,1), B(3,1);

   for (int i = 0; i<4; i++)
    for (int j = 0; j<4; j++)
     m4x4[i][j] = Space[i][j];

   if (OnlyOne) EigenValue(m4x4);
	complex<double> norm = 1.-J;

   for (int b = 0; b < 3; b++)
   	B[b][0] = - m4x4[b][3] * norm;

   for (int i = 0; i < 4; i++)
   {
      m3x3.Diag(Eigen[i](), 0.);

   for (int i = 0; i<3; i++)
    for (int j = 0; j<3; j++)
      m3x3[i][j] = m4x4[i][j] - m3x3[i][j];

      X = m3x3.Back() * B;
      for (int b = 0; b < 3; b++)
         EigenV[b][i] = X[b][0];
      EigenV[3][i] = norm;
   }
}

void xRunge::SetCanonic(matrixU space, double ro)
{
   complex<double> RO;
   matrixU T1, T2, S, J_;
   S.Simplectic();
   J_.Diag(0,0);
   J_[2][0] = 1;
   J_[0][2] = -1;

   matrixU M;
   M = space;

   T1 = M.Trans() * S;
   T1 = T1 * M;
   T1 = S - T1;
   T2 = M.Trans() * J_;
   T2 = T2 * M;
   T2 = T2 - J_;

   RO = T2[1][0] / T1[1][0];
//   RO = ro;

   Canonic[1][2] = -1./(2.*RO);
   Canonic[3][0] =  1./(2.*RO);
}

//**********  Estimation of non-linear equations ********************
bool xRunge::SNEquation (SNUP function, matrixU& unknown, PVOID pvoid, double accuracy)
{
  int R = unknown.Row;
  matrixU M_new(R,1);
  matrixU M_current(R,1);
  matrixU tmp(R,1);
  matrixU tmp1(R,1);
  matrixU tmp_p(R,1);
  matrixU Jacob(R,R);

  complex<double> maxim, detJ;
  matrixU Dx(R,1);

  double datchik=0;
 do
 {
  M_new = unknown;
  function(M_current, M_new, pvoid);

   for (int j = 0; j < R; j++)
   {
      Dx[j][0] = M_new[j][0]*accuracy;
      tmp_p[j][0] = M_new[j][0];
   }

  for (int i = 0; i < R; i++)
  {
      tmp_p[i][0] = M_new[i][0] + Dx[i][0];
      if (i>0) tmp_p[i-1][0] = M_new[i-1][0];

      function(tmp1, tmp_p, pvoid);

   for (int k = 0; k < R; k++)
    {
        Jacob[k][i] = (tmp1[k][0] - M_current[k][0])/Dx[i][0];
    }
  }

  detJ = Jacob.Det();
  if (abs(detJ) == 0) return true;

  tmp = Jacob.Back()*M_current*(-1);

     maxim =  tmp[0][0]/M_new[0][0];
     for (int i = 1; i < R; i++ )
      {
      if ( (abs(tmp[i][0]/M_new[i][0])) > abs(maxim))
         maxim =  tmp[i][0];
      }
   unknown = M_new + tmp;
   datchik++;
   }  while (abs(maxim) > accuracy);

 return true;
}
*/
