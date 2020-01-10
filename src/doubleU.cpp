//---------------------------------------------------------------------------
#include "stdafx.h"
#define _CRT_SECURE_NO_WARNINGS
#include "doubleU.h"
#include <stdlib.h>

int ExpKick = 0;

//---------------------------------------------------------------------------
double light = 2.99792458e10;
bool dini = true;
void Dini()
{  dini = false;

// Attachments
   a_18(ATTO);
   f_15(FEMTO);
   p_12(PICO);
   n_9 (NANO);
   u_6 (MICRO);
   m_3 (MILI);
   k_3 (KILO);
   M_6 (MEGA);
   G_9 (GIGA);
   T_12(TERA);
   P_15(PETA);
   E_18(EXA);

// Base units
   s_ (SECOND);
   Q_ (QULOUN);
   m_ (METER);
   U1_(UNIT);

// Some derivative units
   cm_  (m_/100);                                        // centimetre
   barn_((cm_^2)*1e-24);                                 // barn
   H_   (cm_*1e9);                                       // Henry
   Hz_  (s_^-1);                                         // Herz
   min_ (s_*60);                                         // minute
   hour_(min_*60);                                       // hour
   day_ (hour_*24);                                      // day
   year_(day_*365);                                      // year

// Electromagnetic, perveance
   e_  (Q_*1.6021893e-19);                               // electron charge
   q_  (Q_/3e9);                                         // SGSI charge
   F_  (cm_*9e11);                                       // Farade
   V_  (Q_/F_);                                          // Volt
   A_  (Q_/s_);                                          // Ampere
   Ohm_(V_/A_);                                          // Ohm
   uP_ (u_6*A_/(V_^1.5));                                // Perveance
   G_  (V_/cm_*300);                                     // Gauss
   T_  (G_*1e4);                                         // Tesla

// Power, temperature, momentum, massa
   eV_   (e_*V_);                                        // electron Volts
   K_    (eV_*8.61738573e-5);                            // Kelvin
   eV_c2_(eV_*((s_/(cm_*light))^2));                     // mass
   eV_c_ (eV_*  s_/(cm_*light));                         // momentum
   J_    (Q_*V_);                                        // Joule
   W_    (J_/s_);                                        // Watt
   kg_   (J_*(s_/m_^2));                                 // kilogram
   g_    (kg_/1000);                                     // gram

// Massa of particles
   amu_      (M_6*eV_c2_ * 931.5016);                    // atom mass unit
   neutron_  (M_6*eV_c2_ * 939.5731);                    // neutron mass
   proton_   (M_6*eV_c2_ * 938.2796);                    // proton mass
   electron_ (M_6*eV_c2_ * 0.5110034);                   // electron mass

// Force, pressure
   N_   (J_/m_);                                         // Newton
   Pa_  (N_/(m_^2));                                     // Pascal
   atm_ (Pa_*1.013e5);                                   // atmosphere
   Torr_(atm_/760);                                      // Torr

//Base constants
   U_0._      ( 0,                 U1_ );                // zero
   U_1._      ( 1,                 U1_ );                // unit
   U_hbar._   ( 0.658212202e-15,   eV_*s_ );             // Plank constant
   U_amu._    ( 1,                 amu_ );               // atom mass unit
   U_mn._     ( 1,                 neutron_ );           // neutron mass
   U_mp._     ( 1,                 proton_ );            // proton mass
   U_me._     ( 1,                 electron_ );          // electron mass
   U_k._      ( 8.61738573e-5,     eV_/K_ );             // Boltzman
   U_c._      ( light/100.,        m_/s_ );              // light speed
   U_e._      ( 1,                 e_ );                 // electron charge
   U_Grav._   ( 6.6725985e-8,     (cm_^3)/(g_*(s_^2)));  // gravitational
   U_grav._   ( 9.80665,           m_/(s_^2) );          // acceleration
   U_NA._     ( 6.022136736e23,    U1_ );                // Avagadro
   U_pi._     ( M_PI,              U1_ );                // pi
   U_exp._    ( 2.718281828459045235,U1_ );              // exponent
   U_eps0._   ( 1./(4.*U_pi()*9e+9), F_/m_ );            // permittinity
   U_mu0._    (     4.*U_pi()*1e-7,  H_/m_ );            // permeability

// Derivative constants
   U_fine._   ( (U_e^2) / (U_hbar * U_c) );              // fine structure
   U_re._     ( (U_e^2) / (U_me*(U_c^2)) );              // electron radius
   U_rp._     ( (U_e^2) / (U_mp*(U_c^2)) );              // proton radius
   U_r1._     ( (U_hbar^2) / (U_me*(U_e^2)) );           // first Borh radius
   U_R._      ( U_me * (U_e^4) / ( (U_hbar^3) * 2) );    // Ridberg
   U_Ecoup._  ( U_me * (U_e^4) / ( (U_hbar^2) * 2) );    // coupling energy
   U_lambdaC._( U_pi * U_hbar * 2 / (U_me * U_c) );      // Compthon length
   U_muBorn._ ( U_e  * U_hbar     / (U_me * U_c * 2) );  // Borh magnethon

}
//---------------------------------------------------------------------------

const double accuracy = 1e-10;

MicroMega a_18, f_15, p_12, n_9, u_6, m_3;
MicroMega k_3, M_6, G_9, T_12, P_15, E_18;
Units_ s_(SECOND),Q_(QULOUN),m_(METER),U1_(UNIT);
Units_ cm_(UNIT),barn_(UNIT),H_(UNIT),Hz_(UNIT);
Units_ min_(UNIT),hour_(UNIT),day_(UNIT),year_(UNIT);
Units_ e_(UNIT),q_(UNIT),F_(UNIT),V_(UNIT),A_(UNIT);
Units_ Ohm_(UNIT),uP_(UNIT),G_(UNIT),T_(UNIT);
Units_ eV_(UNIT),K_(UNIT),eV_c_(UNIT),eV_c2_(UNIT);
Units_ J_(UNIT),W_(UNIT),kg_(UNIT),g_(UNIT);
Units_ amu_(UNIT),neutron_(UNIT),proton_(UNIT),electron_(UNIT);
Units_ N_(UNIT),Pa_(UNIT),atm_(UNIT),Torr_(UNIT);

doubleU U_0(UNIT),U_1(UNIT),U_hbar(UNIT),U_amu(UNIT),U_mn(UNIT);
doubleU U_mp(UNIT),U_me(UNIT),U_k(UNIT),U_c(UNIT),U_e(UNIT);
doubleU U_Grav(UNIT),U_grav(UNIT),U_NA(UNIT),U_pi(UNIT);
doubleU U_exp(UNIT),U_eps0(UNIT),U_mu0(UNIT);
doubleU U_fine(UNIT),U_re(UNIT),U_rp(UNIT),U_r1(UNIT);
doubleU U_R(UNIT),U_Ecoup(UNIT),U_lambdaC(UNIT),U_muBorn(UNIT);

//--------------------------------------------------------------------------

Units_::Units_()
{
   if (dini) Dini();
   k = 1.;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      sQm[i] = 0;
#endif
}

Units_::Units_(double koef)
{
   k = koef;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      sQm[i] = 0;
#endif
}

void Units_::operator()(BaseEnum baseenum)
{
   k = 1.;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      sQm[i] = 0;
   if (baseenum != UNIT)
      sQm[baseenum] = 1;
#endif
}

void Units_::operator()(Units_ u)
{
   k = u.k;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      sQm[i] = u.sQm[i];
#endif
}

#ifdef PoWeRs
void Units_::UnitsWarning()
{
        Warning("  s^",sQm[0],", Q^",sQm[1],", m^",sQm[2]);
}
#endif

bool Units_::operator==(Units_ u)
{
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      if ( fabs(sQm[i] - u.sQm[i]) > accuracy )
        return false;
#else
   Warning("! --- Comparison of dimensions is not correct --- !");
#endif
   return true;
}

bool Units_::operator!=(Units_ u)
{
   return !(*this == u);
}

bool Units_::Zero()
{
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      if ( fabs(sQm[i]) > accuracy)
         return false;
#endif
   return true;
}

Units_  Units_::operator*(Units_ u)
{
   Units_  tUnit;
   tUnit.k = k * u.k;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      tUnit.sQm[i] = sQm[i] + u.sQm[i];
#endif
   return tUnit;
}

Units_  Units_::operator/(Units_ u)
{
   Units_  tUnit;
   tUnit.k = k / u.k;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      tUnit.sQm[i] = sQm[i] - u.sQm[i];
#endif
   return tUnit;
}

Units_  Units_::operator*(double koef)
{
   Units_  tUnit;
   tUnit.k = k * koef;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      tUnit.sQm[i] = sQm[i];
#endif
   return tUnit;
}

Units_  Units_::operator/(double koef)
{
   Units_  tUnit;
   tUnit.k = k / koef;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
      tUnit.sQm[i] = sQm[i];
#endif
	return tUnit;
}

Units_  Units_::operator^(double power)
{
	Units_  tUnit ;
   tUnit.k = pow(fabs(k), power);
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
   	tUnit.sQm[i] = sQm[i] * power;
#endif
	return tUnit;
}

Units_  Units_::operator*(MicroMega mega)
{
	Units_  tUnit ;
   tUnit.k = k * mega.MegaValue;
#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
   	tUnit.sQm[i] = sQm[i];
#endif
	return tUnit;
}

//---------------------------------------------------------------------------
/*
MicroMega::MicroMega(MegaEnum megaenum)
{	MegaValue = pow(10, megaenum * 3);
}
*/
void MicroMega::operator()(MegaEnum megaenum)
{	MegaValue = pow(10, megaenum * 3.);
}

Units_  MicroMega::operator*(Units_  u)
{
	Units_  tUnit;
   tUnit.k = u.k * MegaValue;
	#ifdef PoWeRs
   for (int i = 0; i < UNIT; i++)
   	tUnit.sQm[i] = u.sQm[i];
	#endif
	return tUnit;
}
//---------------------------------------------------------------------------

doubleU::doubleU()
{
	v = 1;
}

doubleU::doubleU(double v1)
{
	v = v1;
}

doubleU::doubleU(Units_ u1)
{
	v = 1;
   u = u1;
}

doubleU::doubleU(double v1, Units_ u1)
{
	v = v1;
   u = u1;
}

doubleU::doubleU(double v1, double u1)
{
	v = v1;
   u.k = u1;
}

doubleU::doubleU(BaseEnum unit) : u(unit) { ; }

#ifdef PoWeRs
bool doubleU::CheckUnits(doubleU du, char* sign)
{
	if (u != du.u)
   {	Warning("Can not execute operator", EmptyData,
              sign, EmptyData,
              " for different Units_ :");
   	u.UnitsWarning();
      du.u.UnitsWarning();
//**********************************************************
      return false;    // ******** Debuger point ***********
//**********************************************************
   }
	return true;
}

bool doubleU::CheckZero(char* sign)
{
	if (!u.Zero())
   {	Warning("Can not execute operator",EmptyData,
              sign, EmptyData,
              " for Non Zero Units_ :");
   	u.UnitsWarning();
      return false;
   }
	return true;
}
#endif

double doubleU::operator()(Units_ u1)
{
	#ifdef PoWeRs
	if (u != u1)
   {	Warning("Can not execute operator(Units_ ) for different Units_ :");
   	u.UnitsWarning();
   	u1.UnitsWarning();
      return false;
   }
	#endif
  	return v * u.k / u1.k;
}

void doubleU::_(double d)
{
	v = d / u.k;
}

void doubleU::_(double d, Units_ u1)
{
	v = d;
   u = u1;
}

void doubleU::_(double d, double du)
{
	v = d;
   u.k = du;
}

void doubleU::_(Units_ u1)
{
	u = u1;
}

void doubleU::_(doubleU zvars)
{
	v = zvars.v;
   u = zvars.u;
}

void doubleU::operator=(const doubleU& zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "=");
	#endif
   v = zvars.v * zvars.u.k / u.k;
}

void doubleU::operator=(double d)
{
	v = d;
}

void doubleU::operator+=(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "+=");
	#endif
   v += zvars.v * zvars.u.k / u.k;
}

void doubleU::operator-=(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "-=");
	#endif
   v -= zvars.v * zvars.u.k / u.k;
}

void doubleU::operator*=(doubleU zvars)
{
   u = u * zvars.u;
   v *= zvars.v;
}

void doubleU::operator/=(doubleU zvars)
{
   u = u / zvars.u;
   if (zvars.v)
      v /= zvars.v;
   else
      Warning("doubleU::operator/=(doubleU): division by zero");

}

doubleU doubleU::Saturation(doubleU zvars, doubleU level)
{
   doubleU temp;
   temp  = *this;
   temp -= level;
   temp %= zvars;
   temp += level;
   return temp;
}

void doubleU::operator%=(doubleU zvars)
{
   doubleU& d = *this;
   if (d.v)
   {  if (U_Abs(d) > U_Abs(zvars))
         d += zvars;
      else
      {  double p = (zvars/d)(U1_);
         if (p < 0.)
         {  if(p < -111.)
               p = -111.;
            d *= exp(p);
            ExpKick++;
         }else
            d += zvars;
      }
   }   
   if ( fabs(d.v) < 1e-111 )
      d.v = 1e-111;
}

doubleU doubleU::operator+(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "+");
	#endif
	doubleU  temp ;
   temp.u = u;
   temp.v = v + zvars.v * zvars.u.k / temp.u.k;
   return temp;
}

doubleU doubleU::operator-(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "-");
	#endif
	doubleU  temp ;
   temp.u = u;
   temp.v = v - zvars.v * zvars.u.k / temp.u.k;
   return temp;
}

doubleU doubleU::operator-()
{
	doubleU temp;
   temp.v = -v;
   temp.u = u;
   return temp;
}


doubleU doubleU::operator*(doubleU zvars)
{
	doubleU  temp ;
   temp.u = u * zvars.u;
   temp.v = v * zvars.v;
   return temp;
}

doubleU doubleU::operator/(doubleU zvars)
{
	doubleU  temp ;
   temp.u = u / zvars.u;
   if (zvars.v)
      temp.v = v / zvars.v;
   else
      Warning("doubleU::operator/(doubleU): division by zero");
   return temp;
}

doubleU doubleU::operator%(doubleU zvars)
{
	doubleU temp(*this);
   temp %= zvars;
   return temp;
}

doubleU doubleU::operator+(double d)
{
	#ifdef PoWeRs
	CheckZero("+(double)");
	#endif
	doubleU  temp ;
   temp.u = u;
   temp.v = v + d / u.k;
   return temp;
}

doubleU doubleU::operator-(double d)
{
	#ifdef PoWeRs
	CheckZero("-(double)");
	#endif
	doubleU  temp ;
   temp.u = u;
   temp.v = v - d / u.k;
   return temp;
}

doubleU doubleU::operator*(double d)
{
	doubleU  temp ;
   temp.u = u;
   temp.v = v * d;
   return temp;
}

doubleU doubleU::operator/(double d)
{
	doubleU  temp ;
   temp.u = u;
   if (d)
      temp.v = v / d;
   else
      Warning("doubleU::operator/(double): division by zero");
   return temp;
}

doubleU doubleU::operator^(double d)
{
	doubleU  temp ;
   temp.u = u ^ d;
   temp.v = pow(fabs(v), d);
   return temp;
}

bool doubleU::operator==(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "==");
	#endif
	return (v*u.k) == (zvars.v*zvars.u.k);
}

bool doubleU::operator!=(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "!=");
	#endif
	return (v*u.k) != (zvars.v*zvars.u.k);
}

bool doubleU::operator>(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, ">");
	#endif
	return (v*u.k) > (zvars.v*zvars.u.k);
}

bool doubleU::operator>=(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, ">=");
	#endif
	return (v*u.k) >= (zvars.v*zvars.u.k);
}

bool doubleU::operator<(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "<");
	#endif
	return (v*u.k) < (zvars.v*zvars.u.k);
}

bool doubleU::operator<=(doubleU zvars)
{
	#ifdef PoWeRs
	CheckUnits(zvars, "<=");
	#endif
	return (v*u.k) <= (zvars.v*zvars.u.k);
}

//---------------------------------------------------------------------------
doubleU operator+(double d, doubleU du)
{
	#ifdef PoWeRs
	du.CheckZero("(double)+");
	#endif
	doubleU temp;
   temp.u = du.u;
   temp.v = d / du.u.k + du.v;
   return temp;
}

doubleU operator-(double d, doubleU du)
{
	#ifdef PoWeRs
	du.CheckZero("(double)-");
	#endif
	doubleU temp;
   temp.u = du.u;
   temp.v = d / du.u.k - du.v;
   return temp;
}

doubleU operator*(double d, doubleU du)
{
	doubleU temp;
   temp.u = du.u;
   temp.v = d * du.v;
   return temp;
}

doubleU operator/(double d, doubleU du)
{
	doubleU temp;
   temp.u = U1_ / du.u;
   temp.v = d / du.v;
   return temp;
}

double U_Ln(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Ln(doubleU )");
	#endif
   return log(zvars(U1_));
}

double U_Log(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Log(doubleU )");
	#endif
   return log10(zvars(U1_));
}

double U_Exp(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Exp(doubleU )");
	#endif
   return exp(zvars(U1_));
}

double U_Pow(doubleU zvars, double power)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Pow(doubleU )");
	#endif
//.............................................................................. May 2006
   if (power * U_Log(zvars) > 300)
   {  Warning("U_Pow result lager than double range (1e300) :", zvars(),"^",power);
      return 1;
   }
   return pow(fabs(zvars(U1_)), power);
}

double U_Sin(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Sin(doubleU )");
	#endif
   return sin(zvars(U1_));
}

double U_Cos(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Cos(doubleU )");
	#endif
   return cos(zvars(U1_));
}

double U_Tan(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Tan(doubleU )");
	#endif
   return tan(zvars(U1_));
}

double U_Asin(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Asin(doubleU )");
	#endif
   if (U_Abs(zvars)() > 1)
   {  Warning("U_Asin argument larger than unit :", zvars());
      return M_PI / 2.;
   }
   return asin(zvars(U1_));
}

double U_Acos(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Acos(doubleU )");
	#endif
   if (U_Abs(zvars)() > 1)
   {  Warning("U_Acos argument larger than unit :", zvars());
      return 0;
   }
   return acos(zvars(U1_));
}

double U_Atan(doubleU zvars)
{
	#ifdef PoWeRs
	zvars.CheckZero(" U_Atan(doubleU )");
	#endif
   return atan(zvars(U1_));
}

double U_Atan2(doubleU zvars1, doubleU zvars2)
{
	#ifdef PoWeRs
	zvars1.CheckUnits(zvars2, " U_Atan2(doubleU, doubleU)");
	#endif
   return ArcTan2(zvars1(U1_), zvars2(U1_));
}

double ArcTan2(double y, double x)
{  double angle;
   if ( x )
   { 	angle =  atan(y/x);
      if (x < 0)
      {	if (y < 0) angle -=  U_pi();
         else       angle +=  U_pi();
      }
   }else
   {	if ( y )
   	{	if (y < 0) angle = -U_pi()/2.;
   		else       angle =  U_pi()/2.;
      }else
      	angle = 0;
   }
   return angle;
}

double SqrtSum(double a, double b)
{
	return sqrt(a*a + b*b);
}

doubleU U_SqrtSum(doubleU a, doubleU b)
{
	#ifdef PoWeRs
	a.CheckUnits(b, " U_SqrtSum(doubleU, doubleU)");
	#endif
	return (((a*a)+(b*b))^0.5);
}

doubleU  U_Abs(doubleU zvars)
{
	doubleU  temp;
   temp.v = fabs(zvars.v);
   temp.u = zvars.u;
   return temp;
}
//---------------------------------------------------------------------------
bool U_Energy::PerNucleon = true;

U_Energy::U_Energy()
{
	Velocity._(m_/s_);
   Energy.  _(G_9 * eV_);
   Kinetic. _(G_9 * eV_);
   Momentum._(G_9 * eV_c_);
   Regidity._(T_*m_);
   Massa.   _(amu_);
   Charge.  _(e_);
}

void U_Energy::operator=(const U_Energy& ue)
{
	Gamma    = ue.Gamma;
	Gamma2   = ue.Gamma2;
	Gamma3   = ue.Gamma3;
   Beta     = ue.Beta;
   Beta2    = ue.Beta2;
   Beta3    = ue.Beta3;
   Velocity = ue.Velocity;
   Energy   = ue.Energy;
   Kinetic  = ue.Kinetic;
   Momentum = ue.Momentum;
   Regidity = ue.Regidity;
   Massa    = ue.Massa;
   Charge   = ue.Charge;
   Z        = ue.Z;
   A        = ue.A;
}

bool U_Energy::Set(U_EnergyEnum ue)
{
	bool errors = false;
   doubleU amu;
   if (!PerNucleon) amu = A;

	switch(ue)
   {
	 case U_GAMMA:
   	if (Gamma() <= 1)
      {	Warning("GAMMA should be more than unit");
      	errors = true;
		}
		break;
    case U_BETA:
   	if ((Beta() >0) && (Beta() < 1))
	   {	Gamma = (1 - (Beta^2))^(-0.5);
      }else
     	{	Warning("BETA should be less than unit and positive");
      	errors = true;
      }
		break;
    case U_VELOCITY:
   	if (( Velocity() > 0 ) && ( Velocity < U_c ))
	   {	Beta = Velocity / U_c;
			Gamma = (1 - (Beta^2))^(-0.5);
      }else
     	{	Warning("VELOCITY should be less than light speed and positive");
      	errors = true;
      }
		break;
    case U_ENERGY:
		Gamma = Energy / (A * Massa * (U_c^2));
   	if (Gamma() <= 1)
      {	Warning("Total Energy should be more than rest mass");
      	errors = true;
      }
		break;
    case U_KINETIC:
   	if (Kinetic() > 0)
		{	Gamma = Kinetic / (amu * Massa * (U_c^2)) + 1;
      }else
      {	Warning("Kinetic Energy should be positive");
      	errors = true;
      }
 		break;
    case U_MOMENTUM:
   	if (Momentum() > 0)
		{  Energy = (( (Momentum *  U_c    )^2)+
      				 ( (A * Massa * (U_c^2) )^2))^0.5;
      	Gamma = Energy / (A * Massa * (U_c^2));
      }else
      {	Warning("MOMENTUM should be positive");
      	errors = true;
      }
		break;
    case U_REGIDITY:
   	if (Regidity() > 0)
		{  Momentum = Regidity * U_e / U_c;
		   Energy = (( (Momentum *  U_c    )^2)+
      		     ( (A * Massa  * (U_c^2) )^2))^0.5;
      	Gamma = Energy / (A * Massa * (U_c^2));
      }else
      {	Warning("REGIDITY should be positive");
      	errors = true;
      }
		break;
    default: Warning("Unknown zEnergyEnum parameter : ", ue);
   }

   if (!errors)
   {
      if (ue != U_BETA)
         Beta = (1 - (Gamma^-2)) ^ 0.5;
      if (ue != U_VELOCITY)
         Velocity = Beta * U_c;
      if (ue != U_ENERGY)
         Energy = Gamma * A * Massa * (U_c^2);
      if (ue != U_KINETIC)
         Kinetic = (Gamma - 1) * amu * Massa * (U_c^2);
      if (ue != U_MOMENTUM)
			Momentum = (((Energy^2) - (A*Massa*(U_c^2)^2))^0.5) / U_c;
      if (ue != U_REGIDITY)
         Regidity = Momentum * U_c / U_e;
		Gamma2 = Gamma  * Gamma;
      Gamma3 = Gamma2 * Gamma;
      Beta2 = Beta  * Beta;
      Beta3 = Beta2 * Beta;

	}
   return errors;
}

bool U_Energy::Set(U_EnergyEnum ue, doubleU du)
{
	switch(ue)
   {case U_GAMMA:	   Gamma    = du; break;
    case U_BETA:     Beta     = du; break;
    case U_VELOCITY: Velocity = du;	break;
    case U_ENERGY:   Energy   = du;	break;
    case U_KINETIC:  Kinetic  = du;	break;
    case U_MOMENTUM: Momentum = du;	break;
    case U_REGIDITY: Regidity = du;	break;
    default: Warning("Unknown zEnergyEnum parameter : ", ue);
   }
   return Set(ue);
}

bool U_Energy::Set(U_EnergyEnum ue, double d, Units_ u)
{
	doubleU  du(d, u);
   return Set(ue, du);
}
//---------------------------------------------------------------------------

