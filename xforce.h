//---------------------------------------------------------------------------
#ifndef xForceH
#define xForceH
#include "xdynamic.h"

class xEcool;
//---------------------------------------------------------------------------
class xForce : public xData
{
 public:
   int type;             // type of FrForce formula
   doubleU delta;        // Friction force constant
   BData force_transv;   // array for transverse force data
   BData force_long;     // array for longitudinal force data

   xForce();

   doubleU Vtr;          // transverse ion velocities in PRF
   doubleU v[3];         // velocity components for 3D force
   doubleU Ftr;          // tarnsverse forces in PRF
   doubleU f[3];         // forces components for 3D
   doubleU D[3][3];      // difusion tenzor

   void Budker(xFrParam);            // Budker model for friction force of ecool
   void NonMag(xFrParam);            // Non-magnetized model for friction force of ecool
   void DerSkr(xFrParam);            // Derbenev-Skrinsky-Meshkov model for friction force of ecool
   void Parhom(xFrParam);            // Parkhomchuk model for friction force of ecool
   void Toepffer(xFrParam);          // Toepffer model for friction force
   void D3(xFrParam);                // 3D model
   // 06.07
   void D4(xFrParam);                // 3D analyrical model
   //
   void Table (xFrParam);            // tabulated values of calculated friction force (by another program)
   void Linear (xFrParam);           // Linear model for friction force of ecool
   double phi_f(double );            //Functions for Budker's formula calculation
   double Erf(double );              // parameter of error integral
   double phi(double );              // parameter of error integral


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// here changed by GT, 01.02.05

   xSurf FFtr;                               // 3D plot for transverse friction force of electron cooling effect
   xSurf FFlong;                             // 3D plot for longitudinal friction force of electron cooling effect

// for tr. component of tr. force
   doubleU TT_Range1, TT_Range2, TT_Range3;
   int TT_div1, TT_div2, TT_div3;

// for long. component of tr. force
   doubleU TL_Range1, TL_Range2, TL_Range3;
   int TL_div1, TL_div2, TL_div3;

// for tr. component of long. force
   doubleU LT_Range1, LT_Range2, LT_Range3;
   int LT_div1, LT_div2, LT_div3;

// for long. component of long. force
   doubleU LL_Range1, LL_Range2, LL_Range3;
   int LL_div1, LL_div2, LL_div3;

   int interpolation;
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

   //--------01.02.05-----------------Non magnetized friction force----
   int dl;                    //Nuber of integration steps over long velocity
   int dt;                    //Nuber of integration steps over transverse velocity
   int nfi;                   //Nuber of integration steps over azimuthal
   int asimptotic;            //Asimptotoc presentation
   int smooth;                //maximum impact parameter smoothing
   int rhomin;                //constant minimum impact parameter
   int rmsplus;               //min impat Vi + rms_e
   //---------------------------------
   //---------------------------28.11.05---------------------------------------
   int DelPopolo;
   doubleU q_step;
   double q_max;
   //---------------------------28.11.05---------------------------------------
   //--------------------------D-S friction force --------------------
   int N_M;
   int Mag_As;
   int Magnetized;
   int Adiabatic;
   int Fast;
   int Constant;
   int Pestrikov;              //Longitudinal force in accordance with Pestrikov
   int nfiP;                   //number of integration steps over alpha
   //-----------------------------------------------------------------
   // Toepffer formulae 05.06 ----------
   int TFast;
   int Tight;
   int Stretched;              //Types of collisions
   int Tdl;                    //Nuber of integration steps over long velocity
   int Tdt;                    //Nuber of integration steps over transverse velocity
   int Tnfi;                   //Nuber of integration steps over azimuthal
   // 06.07
   // 3D force
   int D3dl;                    //Nuber of integration steps over long velocity
   int D3dx;                    //Nuber of integration steps over horizontal velocity
   int D3dy;                   //Nuber of integration steps vertical velocity
   int From_array;             //Switch from analytical
   // --------------
   int OnGet();
   int OnSet();

   bool ReadForceFile(char* filename, bool);                                  // read tabulated friction force values from file
   double Simple_Interpol(BData& force, double X, double Y);                  // simple interpolation of the table
   double Poly_Interpol(BData& force, double X, double Y);                    // polynomial interpolation of the table
   void FF(xTime&, xEbeam&, xRing&);

   int Surf, Line, Circle; //draw 3-D force and line
   doubleU Vtr_min, Vtr_max, Vlong_min, Vlong_max, div;                      // minimum and maximum values for transverse and long. velocities (for fr.force plotting)
   doubleU V_min, V_max, Angle, div2;                                        // minimum and maximum velocity values, angle and div number
   doubleU Velocity, div3;
};

extern xForce iForce;

#endif
