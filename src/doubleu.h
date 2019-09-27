#ifdef _DEBUG
#define PoWeRs // Comment this line for more fast calculation
#endif
//---------------------------------------------------------------------------
#ifndef doubleUH
#define doubleUH
#include <math.h>
#ifndef M_PI
#define M_PI (3.141592653589793)
#endif
//---------------------------------------------------------------------------
class Units_;
class MicroMega;
class doubleU;

extern MicroMega a_18, f_15, p_12, n_9, u_6, m_3;
extern MicroMega k_3, M_6, G_9, T_12, P_15, E_18;

extern Units_ s_, Q_, m_, U1_;
extern Units_ cm_, barn_, H_, Hz_, min_, hour_, day_, year_;
extern Units_ e_, q_, F_, V_, A_, Ohm_, uP_, G_, T_;
extern Units_ eV_, K_, eV_c_, eV_c2_, J_, W_, kg_, g_;
extern Units_ amu_, neutron_, proton_, electron_;
extern Units_ N_, Pa_, atm_, Torr_;

extern doubleU U_0,U_1,U_hbar,U_amu,U_mn,U_mp,U_me,U_k,U_c,U_e;
extern doubleU U_Grav,U_grav,U_NA,U_pi,U_exp,U_eps0,U_mu0;
extern doubleU U_fine,U_re,U_rp,U_r1,U_R,U_Ecoup,U_lambdaC,U_muBorn;

//---------------------------------------------------------------------------
enum  BaseEnum{SECOND, QULOUN, METER, UNIT};

class Units_
{
public:
#ifdef PoWeRs
   void UnitsWarning();
   double sQm[UNIT];
#endif
   double k;
   Units_ ();
   Units_ (double);
   Units_ (BaseEnum) {;}
   void operator()(BaseEnum);
   void operator()(Units_);
   bool operator==(Units_);
   bool operator!=(Units_);
   bool Zero();
	Units_  operator*(Units_);
	Units_  operator/(Units_);
	Units_  operator*(double);
	Units_  operator/(double);
	Units_  operator^(double);
	Units_  operator*(MicroMega);
};

enum MegaEnum{ATTO=-6, FEMTO,PICO, NANO, MICRO,MILI,
				  KILO= 1, MEGA, GIGA, TERA, PETA, EXA};
class MicroMega
{public:
	double MegaValue;
   //MicroMega(MegaEnum);
   void operator()(MegaEnum);
	Units_ operator*(Units_);
};

class doubleU
{
 public:
   double v;
	Units_ u;
   doubleU();
   doubleU(double);
   doubleU(Units_);
   doubleU(double, Units_);
   doubleU(double, double);
   doubleU(BaseEnum);
#ifdef PoWeRs
   bool CheckUnits(doubleU, char*);
   bool CheckZero (char*);
#endif
   double operator()() { return v; }
   double operator()(Units_);

   void _(double);
   void _(double, Units_);
   void _(double, double);
   void _(Units_);
   void _(doubleU);

	void operator=(const doubleU&);
   void operator=(double d);

   void operator+=(doubleU);
   void operator-=(doubleU);
   void operator*=(doubleU);
   void operator/=(doubleU);
   void operator%=(doubleU);

	doubleU  operator+(doubleU);
	doubleU  operator-(doubleU);
   doubleU  operator-();
	doubleU  operator*(doubleU);
	doubleU  operator/(doubleU);
	doubleU  operator%(doubleU);
   doubleU  Saturation(doubleU zvars, doubleU level);

	doubleU  operator+(double);
	doubleU  operator-(double);
	doubleU  operator*(double);
	doubleU  operator/(double);
	doubleU  operator^(double);

   bool operator==(doubleU);
   bool operator!=(doubleU);
   bool operator> (doubleU);
   bool operator>=(doubleU);
   bool operator< (doubleU);
   bool operator<=(doubleU);
};

//---------------------------------------------------------------------------

doubleU operator+(double, doubleU);
doubleU operator-(double, doubleU);
doubleU operator*(double, doubleU);
doubleU operator/(double, doubleU);

double U_Ln (doubleU);
double U_Log(doubleU);
double U_Exp(doubleU);
double U_Pow(doubleU, double);
double U_Sin(doubleU);
double U_Cos(doubleU);
double U_Tan(doubleU);
double U_Asin(doubleU);
double U_Acos(doubleU);
double U_Atan(doubleU);
double U_Atan2(doubleU, doubleU);
double ArcTan2(double, double);
double SqrtSum(double a, double b);
doubleU U_SqrtSum(doubleU, doubleU);
doubleU U_Abs(doubleU);

//---------------------------------------------------------------------------
enum U_EnergyEnum{U_GAMMA, U_BETA, U_VELOCITY, U_ENERGY,
					  	U_KINETIC, U_MOMENTUM, U_REGIDITY};
class U_Energy
{public:
	static bool PerNucleon;
	doubleU Gamma, Beta, Velocity, Energy;
   doubleU Kinetic, Momentum, Regidity, Massa;
   doubleU Charge, Z, A;
   doubleU Gamma2, Gamma3, Beta2, Beta3;
   U_Energy();
	void operator=(const U_Energy&);
   bool Set(U_EnergyEnum);
   bool Set(U_EnergyEnum, doubleU);
   bool Set(U_EnergyEnum, double, Units_);
};

extern int ExpKick;

//---------------------------------------------------------------------------
#endif
