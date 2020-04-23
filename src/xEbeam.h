//---------------------------------------------------------------------------
#ifndef xEbeamH
#define xEbeamH
#include "bolideU.h"
#include "xRing.h"
#include "xBeam.h"
#include "xOptics.h"

//---------------------------------------------------------------------------
class xFrParam
{
public:
   xFrParam();

   doubleU Z;            // ion Z
   doubleU A;            // ion A
   doubleU mfield;       //longitudinal field in cooler in kG
   doubleU Ttemp;        // local transverse temperature in meV
   doubleU Ltemp;        //Longitudinal and transverse electron temperatures in meV
   doubleU Ttemp_centre; // transverse temperature in the beam centre in meV
   doubleU TempEff;      // effective tenperature in meV (for Parkhomchuk formula)
   doubleU Smoos;        // smoosing coefficient for Derbenev-Skrinsky formula
   doubleU Theta_Eff;    // effective electron angular spread
   doubleU V_eff_e;      // effective electron velocity spread
   doubleU V_tr_e;       // electron velocity spread in transverse plane
   // ---- 05.06   for 3D force
   doubleU V_tr_x; // electron velocity spread in horizontal plane
   doubleU V_tr_y; // electron velocity spread in vertical plane
   // ----------------
   doubleU V_long_e; // electron velocity spread along the magnetic field
   doubleU n_e;      //electron density in PRF
   doubleU tau;      //time of flight in PRF
   doubleU T_plasma; //period of plasma oscillations in PRF
   //--------11.02.05----------------Undulator-------
   int undulator;   //Recombination suppression using undulator
   doubleU lambda;  //Period
   doubleU B_field; //Field at axis
   doubleU r_0;     //Rotation radius
   doubleU Theta_U; //Coherent rotation angle
   doubleU V_und;   //Coherent rotation velocity
};

class xEbeam : public xData, public xBunch
{
public:
   xBeam *pBeam;
   xEbeam(xRing &ring);

   doubleU LinearRates[3]; // Cooling Rates for Linear model of Forces
   U_Energy e_Energy;
   bool inside; // check if the particle is inside the cooled beam or not
   doubleU *pZ; // pointer on ion Z
   doubleU *pA; // pointer on ion A

   xFrParam F;     // set of parameters for friction force calculation
   doubleU *pBeta; //pointer on ion ring Beta

   //Parameters for uniform cilinder
   doubleU bradius;  //electron beam radius
   doubleU bcurrent; //electron current
   doubleU dVdr;     //Derivative of transverse temperature over radius
   bool LongPulse;   //Longitudinal pulsed ebeam
   double PulseFrom; //from longitudinal position, Circumference
   double PulseUpTo; //up to longitudinal position, Circumference

   double Neutralization; // neutralization factor
   doubleU CoLength;      // the length of the cooling section
   int model;             // switch of electron beam geometry model
   doubleU e_emit_tr;     // emittance of the electron beam
   doubleU e_dpp;         // momentum spread of the electron beam
   int emit;
   int temp;
   int velocity;    // emittance or beam temperature definition
   int spread_temp; // angular spread or effective temperature definition

   //Electrn beam alighment
   bool EnableShifts;
   bool EnableFinal;
   vectorU Initial;
   vectorU Final;
   vectorU Shift;
   bool UsePeriod;
   int PaintingPeriod;
   bool FromFile;
   BData Solenoid;
   bool PaintingTable;
   BData PaintingData;
   int ScallingFactor;

   // 04.12.06
   doubleU Bdistance; // Distance between electron bunches, m
   int BunchNumber;   // Number of electron bunches

   //parameters of Uniform Gaussian bunch
   doubleU size_x; //horizontal dimension
   doubleU size_y; //vertical dimension
   doubleU size_s; //RMS bunch length
   //doubleU dist_uni;                               //Distance between ion and electron bunch centers
   doubleU Ne_uni; //Electron number in the bunch
   doubleU Ie_uni; //Peak electron current

   //     06.07
   //parameters of Gaussian bunch
   bool From_model; //
   vectorU Gsigma;  //RMS dimension
   doubleU Ne;      //Electron number in the bunch
   doubleU Ieb;     //Peak electron current

   //parameters of Gaussian cilinder
   doubleU sigma_x_cil; //RMS horizontal dimension
   doubleU sigma_y_cil; //RMS vertical dimension
   doubleU Ne_cil;      //Linear electron density
   doubleU Ie_cil;      //Electron current

   // 05.06 parameters of electron array 06.07
   bool GaussN; //local density in accordance with Gaussian low
   bool AFromFile;
   bool BiGauss;      //Generation of bi-Gaussian distribution
   bool Uniform;      //Generation of uniform transverse distribution
   int N_core;        //Number of particles inaside core
   double sigma_core; //ratio between core and tail rms
   double test11;
   double AEnergy;
   double AGamma;  // gamma divided by gamma+1
   BData EArray;   // electron bunch parameters from file
   vectorU Asigma; //RMS dimension
   doubleU ANe;    //Electron number in the bunch
   bool Slice;
   double From;
   double Upto;
   double Box;
   int N_hor;
   int N_long;
   //doubleU m_x;
   //doubleU m_y;
   //doubleU m_p;                                   //components of mean electron velocity
   //  05.06

   //parameters of Hollow beam
   doubleU r_hole;   //the hole radius
   doubleU n_hole;   //electron density in lRF in the hole
   doubleU I_hole;   //electron current for hole
   doubleU r_circle; //the circle radius
   doubleU n_circle; //density in the circle
   doubleU I_circle; //electron current for circle
   doubleU I_hocle;  //current for hole with circle density
   bool SpaceCharge;

   //        06.07
   //parameters of parabolic cylinder
   doubleU pradius;  //electron beam radius
   doubleU pcurrent; //electron current
   doubleU n_centre; //density at the axis
   doubleU dVdrp;    //Derivative of transverse temperature over radius

   ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07
   // density and trnacverse velicity gradient distributions from file
   doubleU dradius;
   doubleU dVdrd;
   BData Density;

   //-----------------------------------------------------------

   void PRF_2_LRF(xEbeam, doubleU &, doubleU &);

   vectorU Ibeam2Ebeam(xTime &, vectorU); // function transforms ion coords to electron beam system
   vectorU Ebeam2Ibeam(vectorU &);        // function transforms force components to ion beam system

   vectorU CoordinateShifts(xTime &t, vectorU &ion); // solenoid errors
   vectorU UniBunch(xTime &t, vectorU &ion);         // models electron beam like a transversly uniform bunch
   vectorU UniCilinder(xTime &t, vectorU &ion);      // models electron beam like a uniformly charged cilinder
   vectorU GaussCilinder(xTime &t, vectorU &ion);    // models electron beam like a Gaussian charged cilinder
   vectorU GaussBunch(xTime &t, vectorU &ion);       // models electron beam like a bunch with Gaussian distribution
   vectorU HollowBeam(xTime &t, vectorU &ion);       // models hollow electron beam
   vectorU Array(xTime &t, vectorU &ion);            // elecron co-ordinates from file
   vectorU Array2(xTime &t, vectorU &ion);           // elecron co-ordinates from file
                                                     //   06.07
   vectorU Parabolic(xTime &t, vectorU &ion);        // parabolic cylinder
   vectorU DensityFile(xTime &t, vectorU &ion);

   doubleU UC_dp_P(doubleU r, doubleU Ie, doubleU Re);   // calculates momentum shift for Uniformely distr. cilinder
   doubleU UC_Vdrift(doubleU r, doubleU Ie, doubleU Re); // calculates drift angle (drift velocity over long. velocity) for Uniformely distr. cilinder
   doubleU GC_dp_P(doubleU r);                           // calculates momentum shift for Gaussian distr. cilinder
   doubleU GC_Vdrift(doubleU r);                         // calculates drift angle (drift velocity over long. velocity) for Gaussian distr. cilinder
   double FormFactor(double t);                          // calculates field distribution in radial direction, returns dimensionless result
   // 05.06
   int OnGet();
   int OnSet();
   int OnRun();
};

extern xEbeam iEbeam;

#endif
