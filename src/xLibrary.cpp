//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xLibrary.h"
#include "xDynamic.h"

//---------------------------------------------------------------------------
xDrift::xDrift()
{
   OpticName = "DRIFT";
   EL_FORCES = true;
   EL_MATRIX = true;
}

vectorU xDrift::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);
   forces[0] = ion[1];
   forces[1] = 0;
   forces[2] = ion[3];
   forces[3] = 0;
   forces[4] = ion[5] / e.Gamma2;
   forces[5] = 0;                    
   return forces;
}

matrixU xDrift::GetMatrix (xTime&t,U_Energy&e)
{
   matrixU matrix;
   matrix(0, 1, Length);
   matrix(2, 3, Length);
   matrix(4, 5, Length / e.Gamma2);
   return matrix;
}

//---------------------------------------------------------------------------
xBend::xBend()
{
   OpticName = "SBEND";
   EL_FORCES = true;
   EL_MATRIX = true;
   ANGLE._(0);
   E1._(0);
   E2._(0);
   K1._(0, m_^-2);
   K2._(0, m_^-3);
   HGAP1._(0);
   HGAP2._(0);
   FINT1._(0);
   FINT2._(0);
}

vectorU xBend::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);
   if (ANGLE())
   {  forces[1] += -(ion[0] /((Length/ANGLE)^2))
                   +(ion[5] / (Length/ANGLE));
      forces[4] += - ion[0] / (Length/ANGLE);
   }
   if (K1())
   {  forces[1] += - K1 * ion[0];
      forces[3] +=   K1 * ion[2];
   }
   if (K2())
   {  forces[1] += - K2/2 * ((ion[0]^2)-(ion[2]^2));
      forces[3] +=   K2   *   ion[0]   * ion[2];
   }
   return forces;
}

vectorU xBend::OnEnter(xTime&t,U_Energy&e,vectorU ion)
{
   if (ANGLE())
   {  ion[1] += ion[0] * U_Tan(E1) / (Length/ANGLE);
      ion[3] -= ion[2] * U_Tan(E1) / (Length/ANGLE);
   }
   return ion;
}

vectorU xBend::OnExit (xTime&t,U_Energy&e,vectorU ion)
{
   if (ANGLE())
   {	ion[1] += ion[0] * U_Tan(E2) / (Length/ANGLE);
	   ion[3] -= ion[2] * U_Tan(E2) / (Length/ANGLE);
   }
   return ion;
}

matrixU xBend::GetMatrix (xTime&t,U_Energy&e)
{
   matrixU matrix;
   iDrift.Length = Length;
   matrix = iDrift.GetMatrix(t,e);

   complexU h;
   h = (ANGLE/Length)();

   if (ANGLE() || K1())
   {
      complexU cx, sx, cy, sy, kx, ky, dx, J1;

      kx = sqrt(complexU((h*h) + K1()));
      ky = sqrt(complexU(-K1()));

      cx = cos(kx*Length());
      sx = sin(kx*Length())/kx;

      if (ANGLE == 0)
      {  dx = 0;
         J1 = 0;
      }else
      {  dx = (1. - cos(kx*Length()))/(kx*kx);
         J1 = (Length() - sx)/(kx*kx);
      }

      matrix[0][0] = cx;
      matrix[0][1] = sx;
      matrix[0][5] = dx*h;
      matrix[1][0] = -1.*kx*kx*sx;
      matrix[1][1] = cx;
      matrix[1][5] = sx*h;

      if (K1() == 0)
      {  matrix[2][2] = 1;
         matrix[2][3] = Length();
         matrix[3][2] = 0;
         matrix[3][3] = 1;
      }else
      {  cy = cos(ky*Length());
         sy = sin(ky*Length())/ky;

         matrix[2][2] = cy;
         matrix[2][3] = sy;
         matrix[3][2] = -1.*ky*ky*sy;
         matrix[3][3] = cy;
      }

      matrix[4][0] = -1.*h*sx;
      matrix[4][1] = -1.*h*dx;
      matrix[4][5] = (Length()/(e.Gamma2())) - (J1*h*h);
   }

   matrixU Entrance, Exit;
   complexU e1, e2, e1_, e2_;
   e1 = U_Tan(E1);
   e2 = U_Tan(E2);
   e1_ = e1 - (h * 2. * HGAP1() * FINT1() * (1. + sin(e1)*sin(e1)));
   e2_ = e2 - (h * 2. * HGAP2() * FINT2() * (1. + sin(e2)*sin(e2)));
   Entrance[1][0] =   h * e1;
   Entrance[3][2] = - h * e1_;
   Exit[1][0] =   h * e2;
   Exit[3][2] = - h * e2_;

   matrix.SetZeroImag();
   matrix = Exit * matrix * Entrance;
   return matrix;
}

//---------------------------------------------------------------------------
xEBend::xEBend()
{
   OpticName = "EBEND";
   EL_MATRIX = true;
   Ro._(1, m_);
   FB._(0, m_);
   FE._(0, m_);
   n._(0);
}

matrixU xEBend::GetMatrix (xTime&t,U_Energy&e)
{
   matrixU matrix;
   iDrift.Length = Length;
   matrix = iDrift.GetMatrix(t,e);

   doubleU kx(m_^-1), ky(m_^-1);
   kx = ( (1 +(n * e.Gamma2)) / (Ro * Ro) )^0.5;
   ky = ( (1 - n)* e.Gamma2   / (Ro * Ro) )^0.5;

   matrix(0, 0, U_Cos(kx*Length));
   matrix(1, 0, U_Sin(kx*Length) * (-kx));
   matrix(0, 1, U_Sin(kx*Length) /   kx );
   matrix(1, 1, U_Cos(kx*Length));

   matrix(2, 2, U_Cos(ky*Length));
   matrix(3, 2, U_Sin(ky*Length) * (-ky));
   matrix(2, 3, U_Sin(ky*Length) /   ky );
   matrix(3, 3, U_Cos(ky*Length));

   matrix(4, 5, Length / e.Gamma2);

   matrixU Entr, Exit;
    
   if (FB() > 0)
   {if (entr)
     Entr(3,2,(1+e.Gamma2)/(24*FB*Ro*Ro)*((FB*FB*(4+e.Gamma2))-(FE*FE*e.Gamma2)));
    if (exit)
     Exit(3,2,(1+e.Gamma2)/(24*FB*Ro*Ro)*((FB*FB*(4+e.Gamma2))-(FE*FE*e.Gamma2)));
   }
   matrix = Exit * matrix * Entr;

   return matrix;
}

//---------------------------------------------------------------------------
xFField::xFField()
{
   OpticName = "FFIELD";
   EL_MATRIX = true;
   Ro._(1, m_);
   FB._(0.1, m_);
   FE._(0, m_);
   n._(0);
}

matrixU xFField::GetMatrix (xTime&t,U_Energy&e)
{
   matrixU Entr;
   if (FB() > 0)
     Entr(3,2,(1+e.Gamma2)/(24*FB*Ro*Ro)*((FB*FB*(4+e.Gamma2))-(FE*FE*e.Gamma2)));
   return Entr;
}

//---------------------------------------------------------------------------
xQuadrupole::xQuadrupole()
{
   OpticName = "QUADRUPOLE";
   EL_FORCES = true;
   EL_MATRIX = true;
   K1._(0, m_^-2);
   TILT._(0);
}

vectorU xQuadrupole::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);
   forces[1] = - K1 * ion[0];
   forces[3] =   K1 * ion[2];
   return forces;
}

matrixU xQuadrupole::GetMatrix (xTime&t,U_Energy&e)
{
   matrixU matrix;
   iDrift.Length = Length;
   matrix = iDrift.GetMatrix(t,e);

   if (K1())
   {
      complexU cx, sx, cy, sy, kx, ky;

      kx = sqrt(complexU(K1()));
      ky = sqrt(complexU(-K1()));

      cx = cos(kx*Length());
      sx = sin(kx*Length())/kx;

      matrix[0][0] = cx;
      matrix[0][1] = sx;
      matrix[1][0] = -1.*kx*kx*sx;
      matrix[1][1] = cx;

      if (K1() == 0)
      {  matrix[2][2] = 1;
         matrix[2][3] = Length();
         matrix[3][2] = 0;
         matrix[3][3] = 1;
      }else
      {  cy = cos(ky*Length());
         sy = sin(ky*Length())/ky;;

         matrix[2][2] = cy;
         matrix[2][3] = sy;
         matrix[3][2] = -1.*ky*ky*sy;
         matrix[3][3] = cy;
      }

      matrix[4][5] = Length()/e.Gamma2();

      matrix.SetZeroImag();
   }

   matrixU Entrance, Exit;

   Entrance[0][0] = U_Cos(TILT);
   Entrance[1][1] = Entrance[0][0];
   Entrance[2][2] = Entrance[0][0];
   Entrance[3][3] = Entrance[0][0];

   Entrance[0][2] = U_Sin(TILT);
   Entrance[1][3] = Entrance[0][2];
   Entrance[2][0] =-Entrance[0][2];
   Entrance[3][1] =-Entrance[0][2];

   Exit = Entrance.Back();

   matrix = Exit * matrix * Entrance;
   return matrix;
}

//---------------------------------------------------------------------------
xSextupole::xSextupole()
{
   OpticName = "SEXTUPOLE";
   EL_FORCES = true;
   K2._(0, m_^-3);
}

vectorU xSextupole::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);
   forces[1] = - K2/2 * ((ion[0]^2)-(ion[2]^2));
   forces[3] =   K2   *   ion[0]   * ion[2];
   return forces;
}
//---------------------------------------------------------------------------
xSolenoid::xSolenoid()
{
   OpticName = "SOLENOID";
   EL_FORCES = true;
   EL_MATRIX = true;
   Ks._(0, m_^-1);
}

vectorU xSolenoid::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);

   forces[1] = ion[3]*Ks;
   forces[3] = (-1.)*ion[1]*Ks;

   return forces;
}
vectorU xSolenoid::OnEnter(xTime&t,U_Energy&e,vectorU ion)
{
   ion[1] += ion[2]*Ks;
   ion[3] -= ion[0]*Ks;
   return ion;
}

vectorU xSolenoid::OnExit (xTime&t,U_Energy&e,vectorU ion)
{
   ion[1] -= ion[2]*Ks;
   ion[3] += ion[0]*Ks;
   return ion;
}
matrixU xSolenoid::GetMatrix (xTime&t,U_Energy&e)
{
   matrixU matrix;

   doubleU C = U_Cos(Ks*Length);
   doubleU S = U_Sin(Ks*Length);

   matrix(0,0, C*C);
   matrix(0,1, (1./Ks)*S*C);
   matrix(0,2, S*C);
   matrix(0,3, (1./Ks)*S*S);

   matrix(1,0, -1*Ks*S*C);
   matrix(1,1, C*C);
   matrix(1,2, -1*Ks*S*S);
   matrix(1,3, S*C);

   matrix(2,0, -1*S*C);
   matrix(2,1, (-1./Ks)*S*S);
   matrix(2,2, C*C);
   matrix(2,3, (1./Ks)*S*C);

   matrix(3,0, Ks*S*S);
   matrix(3,1, -1*S*C);
   matrix(3,2, -1*Ks*S*C);
   matrix(3,3, C*C);

   matrix(4,5, Length / e.Gamma2);

  return matrix;
}
//---------------------------------------------------------------------------

xRFcavity::xRFcavity()
{
   OpticName = "RFCAVITY";
   //EL_FORCES = true;
   EL_MATRIX = true;
   RE_MATRIX = true;
   H._(0.);
   V._(0., M_6*V_);
}

matrixU xRFcavity::GetMatrix (xTime&t,U_Energy&e)
{
   matrixU matrix;
// doubleU omega(2. * U_pi * iRing.h / iRing.Trev);
   doubleU omega(2. * U_pi * H / iRing.Trev);
   doubleU fi, fi0(0);
   fi = fi0 - omega * t.t;
// matrix(5,4, -omega*U_e*iRing.V*U_Cos(fi)/(e.Beta2*U_c*U_c*e.Momentum));
   matrix(5,4, -omega*U_e*V*U_Cos(fi)/(e.Beta2*U_c*U_c*e.Momentum));
   return matrix;
}
//---------------------------------------------------------------------------
SimpleECool::SimpleECool()
{
   EL_FORCES = true;
   Tapered._(0, m_^-1);
   SwitchOn._(0, m_3*s_);
   /*
	Param.Col(1, 5);
   for (int i = 0; i < 3; i++)
   	Param[0](i) = &K[i].v;
 	Param[0](3) = &Tapered.v;
 	Param[0](4) = &SwitchOn.v;
   OpticName = "COOL";
   Param[0](0).ch = "KX";
   Param[0](1).ch = "KY";
   Param[0](2).ch = "KZ";
   Param[0](3).ch = "TAPERED";
   Param[0](4).ch = "OFF";
   */
}

vectorU SimpleECool::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);
   if (t.t < SwitchOn)
   {
      for (int i = 0; i < 3; i++)
         forces[i*2+1] = ion[i*2+1] * K[i] / Length;
      forces[5] += K[2] / Length * (e.Gamma^2) * Tapered * ion[0];
      forces[5] /= 2;
	}
	return forces;
}
//---------------------------------------------------------------------------

xPallas::xPallas()
{
   EL_FORCES = true;
   Udc._(0, V_);
   Urf._(0, V_);
   omega._(0, Hz_);
   radius._(0, m_);
   /*
	Param.Col(1, 5);
   Param[0](0) = &ANGLE.v;
   Param[0](1) = &Udc.v;
   Param[0](2) = &Urf.v;
   Param[0](3) = &omega.v;
   Param[0](4) = &radius.v;
   OpticName = "PALLAS";
   Param[0](0).ch = "ANGLE";
   Param[0](1).ch = "UDC";
   Param[0](2).ch = "URF";
   Param[0](3).ch = "OMEGA";
   Param[0](4).ch = "R0";
   */
}

vectorU xPallas::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);

	doubleU phi(U_pi * 2 * omega * t.t);
	doubleU Ze_pvr2(e.Z*U_e/(e.Momentum*e.Velocity*(radius^2)));

	forces[1]= Ze_pvr2 * ion[2] * (Udc - (Urf*U_Cos(phi)));
	forces[3]= Ze_pvr2 * ion[0] * (Udc - (Urf*U_Cos(phi)));

   //forces[1]-= ion[0] /((Length/ANGLE)^2);
   forces[1]+= ion[5] / (Length/BendAngle);
   forces[4]-= ion[0] / (Length/BendAngle);

   return forces;
}

//---------------------------------------------------------------------------

ConstFocusing::ConstFocusing()
{
   EL_FORCES = true;
  // Param[0](0) = &LambdaL_N;
}

vectorU ConstFocusing::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
   doubleU Kc( ((iBeam.Number() / iBeam.CellSize)^3) * 1.5 *
   (e.Z^2) * U_rp / (e.A * (e.Gamma^5) * (e.Beta^2) *
   pow(LambdaL_N, 3)) );

  	vectorU forces(U1_, m_^-1);
   forces[1] = - Kc * ion[0];
   forces[3] = - Kc * ion[2];
   return forces;
}
//---------------------------------------------------------------------------

RadialField::RadialField()
{
   EL_FIELDS = true;
   Gradient._(0, V_/(m_^2));
   //Param[0](0) = &Gradient.v;
}

vectorU RadialField::GetFields(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU fields(V_/m_, G_);
   fields[0] = Gradient * ion[0];
   fields[2] = Gradient * ion[2];
	return fields;
}
//---------------------------------------------------------------------------

ConstSolenoid::ConstSolenoid()
{
   EL_FIELDS = true;
   BFiled._(0, G_);
   //Param[0](0) = &BFiled.v;
}

vectorU ConstSolenoid::GetFields(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU fields(V_/m_, G_);
   fields[5] = BFiled;
	return fields;
}

//************************************************************************
// optics elements for coupled motion
//************************************************************************

//------------------------ Toroid -------------------------------------------

xCM_Toroid::xCM_Toroid()
{
   OpticName = "CM_TOR";
   EL_FORCES = true;
   EL_MATRIX = true;
   B._(0, G_);
   R._(0, m_);
   Kdis._(0, m_^-1);
}

vectorU xCM_Toroid::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);
  	 forces[1] -= ion[0] * U_Pow((U_pi*R/Length),2) - Kdis*ion[5]*U_pi*R/Length;
   return forces;
}

vectorU xCM_Toroid::OnEnter(xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

vectorU xCM_Toroid::OnExit (xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

matrixU xCM_Toroid::GetMatrix (xTime&t,U_Energy&en)
{
  matrixU matrix;

   double s = t.sr();                                           // current long coord
   double RT = Length()/(M_PI*R());                             // toroid curvature

   double q = 3e8;                                              // KOULON
   double ro = ((en.Velocity*en.Gamma*U_me*U_c)/(U_e*B))();         // Larmour radius
   double j = (iBeam.I/(M_PI*iBeam.a*iBeam.a*en.Velocity*q))();     // electron current density (J)
   double wp2 = (-2.*M_PI*U_e*U_e*j/(en.Gamma3*U_me))();        // plasma Lengmur frequency (Omega_p^2)
   double v = en.Velocity();

	complex<double> w2p = wp2*pow(1/v,2);
	complex<double> w1 = sqrt(w2p + 1 / (RT*RT));
	complex<double> Qtt2 = (w2p + w1*w1 +  1 / (ro*ro))/(double)2.0;
	complex<double> Qt11 = Qtt2 + sqrt( pow(Qtt2, 2) -  w1*w1 * w2p );
	complex<double> Qt22 = Qtt2 - sqrt( pow(Qtt2, 2) -  w1*w1 * w2p );
	complex<double> Qtd = Qt22 - Qt11;

	complex<double> Qt1 = sqrt(Qt11);
	complex<double> Qt2 = sqrt(Qt22);
	complex<double> vco = cos(Qt1 * s) - cos(Qt2 * s);
	complex<double> vsi = sin(Qt1 * s) / Qt1 - sin(Qt2 * s) / Qt2 ;

	matrix[0][0] = (1. / Qtd) * ( Qt22 * cos(Qt1 * s) - Qt11 * cos(Qt2 * s) - pow(w1, 2) * vco );
	matrix[0][2] = -1. * ( 1. / (Qtd*ro) ) * vsi * w2p;
	matrix[0][1] = (1. / Qtd) * (sin(Qt1 * s) * (w2p-Qt11) / Qt1 - sin(Qt2 * s) * (w2p-Qt22) / Qt2 );
	matrix[0][3] = (1./(Qtd * ro)) * (cos(Qt1 * s) - cos(Qt2 * s));

	matrix[2][0] = w1 * w1 * vsi / (Qtd * ro);
	matrix[2][2] = (1. / Qtd) * ( (w2p - Qt11) * cos(Qt2 * s) - (w2p - Qt22) * cos(Qt1 * s) ) ;
	matrix[2][1] = -1. * (1. / Qtd) *  vco / ro;
	matrix[2][3] = (1. / Qtd) * ( Qt2 * sin (Qt2 * s) - Qt1 * sin (Qt1 * s) - pow(w1, 2)*vsi*(-1.));

	matrix[1][0] =-pow(w1, 2) * matrix[0][1];
	matrix[1][2] = -1. * w2p * vco / (Qtd * ro);
	matrix[1][1] = (1. / Qtd) * ( (w2p - Qt11) * cos(Qt1 * s) - (w2p - Qt22) * cos(Qt2 * s) );
	matrix[1][3] = (Qt1 * sin(Qt1 * s) - Qt2 * sin(Qt2 * s) ) / (-1. * Qtd * ro );

	matrix[3][0] = (1. / Qtd) * pow(w1, 2) * vco / ro;
	matrix[3][2] = (1. / Qtd) * (Qt1 * sin(Qt1 * s) * (w2p - Qt22)- Qt2 * sin(Qt2 * s) * (w2p - Qt11));
	matrix[3][1] = (1. / (Qtd * ro)) * ( Qt1 * sin(Qt1 * s) - Qt2 * sin(Qt2 * s) );
	matrix[3][3] =-(1. / Qtd) * (Qt11 * cos(Qt1 * s) - Qt22 * cos(Qt2 * s) - pow(w1, 2) * vco);

	matrix[0][5] = (1.- matrix[0][0])* RT;
	matrix[1][5] =    - matrix[1][0] * RT;
	matrix[2][5] =    - matrix[2][0] * RT;
	matrix[3][5] =    - matrix[3][0] * RT;

  return matrix;
}

//----------------------- Solenoid -------------------------------------------

xCM_Solenoid::xCM_Solenoid()
{
   OpticName = "CM_SOL";
   EL_FORCES = true;
   EL_MATRIX = true;
   B._(0, G_);
}

vectorU xCM_Solenoid::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);
   double ro = ((e.Velocity*e.Gamma*U_me*U_c)/(U_e*B))();         // Larmour radius
  	 forces[1] += ion[3]/ro;
  	 forces[3] += -ion[1]/ro;
   return forces;
}

vectorU xCM_Solenoid::OnEnter(xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

vectorU xCM_Solenoid::OnExit (xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

matrixU xCM_Solenoid::GetMatrix (xTime&t,U_Energy&en)
{
   double s = t.sr();                                           // current long coord
   double q = 3e8;                                              // KOULON
   double v = en.Velocity();
   double ro = ( (en.Velocity*en.Gamma*U_me*U_c)/(U_e*B) )();     // Larmour radius
   double j = ( iBeam.I/(M_PI*iBeam.a*iBeam.a*en.Velocity*q) )(); // electron current density (J)
   double wp2 = ( -2.*M_PI*U_e*U_e*j/(en.Gamma3*U_me) )();        // plasma Lengmur frequency (Omega_p^2)
   double wd = sqrt(1./(ro*ro*4.)- wp2/(v*v));                  // Omega_d - frequency connected to Debye radius (?)

   matrixU matrix;

	matrix[0][0] = cos(s/(2*ro))*cos(wd*s)+sin(s/(2*ro))*sin(wd*s)/(2*ro*wd);
	matrix[0][2] = sin(s/(2*ro))*cos(wd*s)-cos(s/(2*ro))*sin(wd*s)/(2*ro*wd);
	matrix[0][1] = cos(s/(2*ro))*sin(wd*s)/wd;
	matrix[0][3] = sin(s/(2*ro))*sin(wd*s)/wd;

	matrix[2][0] =-matrix[0][2];
	matrix[2][2] = matrix[0][0];
	matrix[2][1] =-matrix[0][3];
	matrix[2][3] = matrix[0][1];

	matrix[1][0] =-cos(s/(2*ro))*sin(wd*s)*wd*(1.-1./(4.*ro*ro*wd*wd));
	matrix[1][2] =-sin(s/(2*ro))*sin(wd*s)*wd*(1.-1./(4.*ro*ro*wd*wd));
	matrix[1][1] = cos(s/(2*ro))*cos(wd*s)-sin(s/(2*ro))*sin(wd*s)/(2*ro*wd);
	matrix[1][3] = sin(s/(2*ro))*cos(wd*s)+cos(s/(2*ro))*sin(wd*s)/(2*ro*wd);

	matrix[3][0] =-matrix[1][2];
	matrix[3][2] = matrix[1][0];
	matrix[3][1] =-matrix[1][3];
	matrix[3][3] = matrix[1][1];


  return matrix;
}

//----------------------- Stright solenoid + spiral quadrupole ----------------

xCM_Sol_Quad::xCM_Sol_Quad()
{
   OpticName = "CM_SOL_Q";
   EL_FORCES = true;
   EL_MATRIX = true;
   B._(0, G_);
   Wind._(0, U1_);
   K1._(0, m_^-1);
   K2._(0, m_^-1);
   Quad = 1;                                                   // ?????????  (sign of the Quad? )
}

vectorU xCM_Sol_Quad::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);

   double s = t.sr();                                          // current long coord
   double ro = ( (e.Velocity*e.Gamma*U_me*U_c)/(U_e*B) )();     // Larmour radius
	double k = 2 * M_PI * Wind()/Length();                      // period of the quadrupole spiral field

  	 forces[1] += Quad*(s*(K2-K1)/Length+K1)*(ion[0]*cos(2*k*s)+ion[2]*sin(2*k*s))/(ro*B);;
  	 forces[3] += Quad*(s*(K2-K1)/Length+K1)*(ion[0]*sin(2*k*s)-ion[2]*cos(2*k*s))/(ro*B);

   return forces;
}

vectorU xCM_Sol_Quad::OnEnter(xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

vectorU xCM_Sol_Quad::OnExit (xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

matrixU xCM_Sol_Quad::GetMatrix (xTime&t,U_Energy&e)
{
  matrixU matrix;

   double q = 3e8;                                              // KOULON
   double v = e.Velocity();                                     // particle Velocity
   double ro = ( (e.Velocity*e.Gamma*U_me*U_c)/(U_e*B) )();     // Larmour radius
   double j = ( iBeam.I/(M_PI*iBeam.a*iBeam.a*e.Velocity*q) )(); // electron current density (J)
   double wp2 = ( -2.*M_PI*U_e*U_e*j/(e.Gamma3*U_me) )();        // plasma Lengmur frequency (Omega_p^2)

	double G  = Quad;
	double xk = 2. * M_PI * Wind()/Length();                      // period of the quadrupole spiral field
	double sq = t.sr();

	double etaq = 2.*xk +1./ro;
	double alfaq = sqrt(etaq*etaq*(1./(4.*ro*ro)-wp2/(v*v))+pow(G/(B()*ro),2));

	double Q22 = (etaq*etaq + 1/(ro*ro))/4 - wp2/(v*v) - alfaq;
	double Q12 = Q22 + 2 * alfaq;
	double Q1  = sqrt(Q12);
	double Q2  = sqrt(Q22);

	double Qr2 = xk*xk + xk/ro - G/(B()*ro) + wp2/(v*v);
	double Qr1 = Q12 + Qr2;
			 Qr2 = Q22 + Qr2;
	double fi1 = Q12 / Qr1;
	double fi2 = Q22 / Qr2;

	double sq2 = xk / (2*xk + 1/ro);
	double sq1 = sq2 - fi1;
			 sq2 = sq2 - fi2;
	double zq1 = Q12/etaq - xk * fi1;
	double zq2 = Q22/etaq - xk * fi2;
	double Qdif = Q12 * fi2 - Q22 * fi1;

	double psi1x = sin(Q1*sq)*cos(xk*sq) - etaq *Q1/Qr1 * cos(Q1*sq)*sin(xk*sq);
	double psi2x = sin(Q2*sq)*cos(xk*sq) - etaq *Q2/Qr2 * cos(Q2*sq)*sin(xk*sq);

	double g1x = cos(Q1*sq)*cos(xk*sq) + etaq *Q1/Qr1 * sin(Q1*sq)*sin(xk*sq);
	double g2x = cos(Q2*sq)*cos(xk*sq) + etaq *Q2/Qr2 * sin(Q2*sq)*sin(xk*sq);

	matrix[0][2] = (zq1*Q2*psi2x - zq2*Q1*psi1x) / Qdif;
	matrix[0][3] = (g2x - g1x) / (etaq*(fi1-fi2));
	matrix[0][0] = (g1x*sq2 - g2x*sq1) / (fi1 - fi2);
	matrix[0][1] = (fi2*Q1*psi1x - fi1*Q2*psi2x) / Qdif;

	double psi1y = sin(Q1*sq)*sin(xk*sq) + etaq*Q1/Qr1 * cos(Q1*sq)*cos(xk*sq);
	double psi2y = sin(Q2*sq)*sin(xk*sq) + etaq*Q2/Qr2 * cos(Q2*sq)*cos(xk*sq);
	double g1y = cos(Q1*sq)*sin(xk*sq) - etaq*Q1/Qr1 * sin(Q1*sq)*cos(xk*sq);
	double g2y = cos(Q2*sq)*sin(xk*sq) - etaq*Q2/Qr2 * sin(Q2*sq)*cos(xk*sq);

	matrix[2][2] = (zq1*Q2*psi2y - zq2*Q1*psi1y) / Qdif;
	matrix[2][0] = (g1y*sq2 - g2y*sq1) / (fi1 - fi2);
	matrix[2][1] = (fi2*Q1*psi1y - fi1*Q2*psi2y) / Qdif;
	matrix[2][3] = (g2y - g1y) / (etaq*(fi1-fi2));

	double psivx01 = Q1*cos(Q1*sq)*cos(xk*sq) + etaq*fi1*sin(Q1*sq)*sin(xk*sq);
	double psivx02 = xk*( sin(Q1*sq)*sin(xk*sq) + etaq*Q1/Qr1*cos(Q1*sq)*cos(xk*sq) );
	double psivx1 = psivx01 - psivx02;
	double psivx03 = Q2*cos(Q2*sq)*cos(xk*sq) + etaq*fi2*sin(Q2*sq)*sin(xk*sq);
	double psivx04 = xk*( sin(Q2*sq)*sin(xk*sq) + etaq*Q2/Qr2*cos(Q2*sq)*cos(xk*sq) );
	double psivx2 = psivx03 - psivx04;

	double gvx01 = - Q1*sin(Q1*sq)*cos(xk*sq) + etaq*fi1*cos(Q1*sq)*sin(xk*sq);
	double gvx02 = xk*(-cos(Q1*sq)*sin(xk*sq) + etaq*Q1/Qr1*sin(Q1*sq)*cos(xk*sq) );
	double gvx1 =   gvx01 +   gvx02;
	double gvx03 =   - Q2*sin(Q2*sq)*cos(xk*sq) + etaq*fi2*cos(Q2*sq)*sin(xk*sq);
	double gvx04 = xk*(-cos(Q2*sq)*sin(xk*sq) + etaq*Q2/Qr2*sin(Q2*sq)*cos(xk*sq) );
	double gvx2 =   gvx03 +   gvx04;

	matrix[1][2] = (zq1*Q2*psivx2 - zq2*Q1*psivx1) / Qdif;
	matrix[1][3] = (gvx2 - gvx1) / (etaq*(fi1-fi2));
	matrix[1][0] = (gvx1*sq2 - gvx2*sq1) / (fi1 - fi2);
	matrix[1][1] = (fi2*Q1*psivx1 - fi1*Q2*psivx2) / Qdif;

	double psivy01 = Q1*cos(Q1*sq)*sin(xk*sq) - etaq*fi1*sin(Q1*sq)*cos(xk*sq);
	double psivy02 = xk*( sin(Q1*sq)*cos(xk*sq) - etaq*Q1/Qr1*cos(Q1*sq)*sin(xk*sq) );
	double psivy1 = psivy01 + psivy02;
	double psivy03 = Q2*cos(Q2*sq)*sin(xk*sq) - etaq*fi2*sin(Q2*sq)*cos(xk*sq);
	double psivy04 = xk*( sin(Q2*sq)*cos(xk*sq) - etaq*fi2/Q2*cos(Q2*sq)*sin(xk*sq) );
	double psivy2 = psivy03 + psivy04;

	double gvy01 = - Q1*sin(Q1*sq)*sin(xk*sq) - etaq*fi1*cos(Q1*sq)*cos(xk*sq);
	double gvy02 = xk*( cos(Q1*sq)*cos(xk*sq) + etaq*Q1/Qr1*sin(Q1*sq)*sin(xk*sq) );
	double gvy1 =   gvy01 +   gvy02;
	double gvy03 = - Q2*sin(Q2*sq)*sin(xk*sq) - etaq*fi2*cos(Q2*sq)*cos(xk*sq);
	double gvy04 = xk*( cos(Q2*sq)*cos(xk*sq) + etaq*Q2/Qr2*sin(Q2*sq)*sin(xk*sq) );
	double gvy2 = gvy03 +   gvy04;

	matrix[3][2] = (zq1*Q2*psivy2 - zq2*Q1*psivy1) / Qdif;
	matrix[3][3] = (gvy2 - gvy1) / (etaq*(fi1-fi2));
	matrix[3][0] = (gvy1*sq2 - gvy2*sq1) / (fi1 - fi2);
	matrix[3][1] = (fi2*Q1*psivy1 - fi1*Q2*psivy2) / Qdif;

  return matrix;
}

//----------------------- Toroid + spiral quadrupole ---------------------

xCM_Tor_Quad::xCM_Tor_Quad()
{
   OpticName = "CM_TOR_Q";
   EL_FORCES = true;
   EL_MATRIX = true;
   B._(0, G_);
   Wind._(0, U1_);
   K1._(0, m_^-1);
   K2._(0, m_^-1);
   Quad = 1;                                                   // ?????????  (sign of the Quad? )
}

vectorU xCM_Tor_Quad::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
  	vectorU forces(U1_, m_^-1);

   double s = t.sr();                                          // current long coord
   double ro = ( (e.Velocity*e.Gamma*U_me*U_c)/(U_e*B) )();     // Larmour radius
	double k = 2 * M_PI * Wind()/Length();                      // period of the quadrupole spiral field

  	 forces[1] += Quad*(s*(K2-K1)/Length+K1)*(ion[0]*cos(2*k*s)+ion[2]*sin(2*k*s))/(ro*B);;
  	 forces[3] += Quad*(s*(K2-K1)/Length+K1)*(ion[0]*sin(2*k*s)-ion[2]*cos(2*k*s))/(ro*B);

   return forces;
}

vectorU xCM_Tor_Quad::OnEnter(xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

vectorU xCM_Tor_Quad::OnExit (xTime&t,U_Energy&e,vectorU ion)
{
   return ion;
}

matrixU xCM_Tor_Quad::GetMatrix (xTime&t,U_Energy&e)
{
  matrixU matrix;
  return matrix;
}

//---------------------------------------------------------------------------

