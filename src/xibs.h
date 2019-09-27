//---------------------------------------------------------------------------
#ifndef xIBSH
#define xIBSH
#include "xeffect.h"
#include "bptune.h"
//---------------------------------------------------------------------------
class xIBS : public xEffect
{
 public:
   int Model;                    // 0 for Piwinski model, 1 for Jie Wei model, 2 for Martini model
   int Model2;                   // 0 for Bjorken model, 1 for Aanalytical model, 2 - numerical integration, 3 - Zenkevich (excluded from interface)
   doubleU JieWeiL;              // Jie Wei Coulomb logarithm
   double BjorkenL;              // Bjorken Coulomb coeff
   double stepmu, stepnu, stepz; // steps over Martini integral
   double B_M_log, lam_limit, step_lam, uplimit;    // for Bjorken-Mtingwa
   double nb;
   double BUCKETnb;              // nb for barrier bucket
   double AverageTrans;          // Coupled - uncoupled transverse motion
   bool coupledJie;              // Coupling in the simplified Jie Wei model
   bool NoDispersion;            // simple approximation for d
   bool DispersionCoeff;         // - koefficient for gas relaxation formula - to take into account <Disp^2/beta_x> or not
   int  IBSkickModel;            // IBS kick model
   bool Slices;                  // Longitudinal slices
   int bpStepSize;               // Integrate step size for JINR model
   int bpStepNumber;             // Integrate step number
   int bpSliceNumber;            // number of longitudinal slices
   bool BMHEA;                   // switch for Bjorken-Mtingwa high energy approximation
   doubleU BK;                   //Test of detail IBS
   doubleU Mx;
   doubleU Mpx;
   doubleU M;
   double CoreNumber;            // Particle number in core
   double CoreDefinition;        // sigma value for core
   int    InCoreDefine;          // number of degree of freedom in core
   // 15.12.06
   bool CoreNumCor;              //Core number definition in Bi-Gauss as in other models
   bool CoreFWHM;                //Core emittance definition in Bi-Gauss via FWHM
   
   bpTune iTune;
   bool TuneShift;
   double CellNum[3];
   double GridSize[3];
   int FieldComponent;
   double VerticalCell;
   bool Plottuneshift;
   BData Bumpfile;
   double TSminX;
   double TSminY;


   //Diffusion tensor and friction components for individual kick 04.06
   doubleU D_xx, D_yy, D_zz, D_zx, D_zy;
   doubleU F_x, F_y, F_z;
   //   06.07 parameters for local model
   int IBSloc;
   double Box_IBS;
   double TranCoup;
   //

   doubleU s_x, s_xb, s_yb, s_s;          //Variables for beam rms calculation
   doubleU s_y, s_x_, s_y_, s_p;          //Variables for beam rate calculation
   doubleU s_h, D_, ri;                   //Variables for koefficient A0 calculation
   doubleU q, a, b, c, d, d_, Ao;         //Variables for koefficient A0 calculation
   void Sigma(xLattice&, xBeam&, xRing&); //calculates rms beam parameters (sigmas, ', etc)

   xIBS();
   int OnGet();
   int OnRun();                                       // 05.06 - ibs kick
   int OnSet();
   vectorU Rates(xTime&, xBeam&, xRing&);             // calculates IBS

   double PiwinskiIntegral(double, double, double);   // calculates Piwinski integral
   vectorU Piwinski(xBeam&, xRing&);                  // calculates IBS rates for Piwinski model

   doubleU Hi_function(doubleU);                      // calculates HI- function (JieWei model)
   vectorU JieWei(xBeam&, xRing&);                    // calculates IBS rates for JieWei model

   double explogint(double D);                                         // approximation by Abramowitz & Stegun
   vectorU MartiniIntegral(double, double, double, double, double);    // calculates Martini integral
   vectorU Martini();                                                  // calculates IBS rates for Martini model

   doubleU Bessel_I0(doubleU);                        // calculates modified Bessel function of 0 order
   doubleU G_F(doubleU, doubleU);                     // calculates function with twice integral - for Burov's model
   vectorU Burov(xLattice&, vectorU, xRing&);         // calculates rates for detailed IBS model (Burov's)

   vectorU Gas_Relaxation(xBeam&, xRing&);           // calculates rates by Gas Relaxation Formula

   vectorU Bjorken_Mtingwa2(xLattice&, xBeam&, xRing&);   // calculates rates by Bjorken_Mtingwa theory
   vectorU Bjorken_Mtingwa(xLattice&, xBeam&, xRing&);   // calculates rates by Bjorken_Mtingwa theory
   double B_M_integral(int, int, matrixU);               // integral in Bjorken_Mtingwa IBS
   double B_M_integral_I(matrixU L, int nomer);          // Fast integral in Bjorken_Mtingwa IBS (simplified)
   vectorU Bjorken_Mtingwa_HE(xBeam&, xRing& );         // calculates rates by Bjorken_Mtingwa theory for H.E.

   vectorU GP_BiGauss(xLattice&, xBeam&, xRing&);        // George Parsen Bi-Gaussian IBS theory

   vectorU MD_Rates(xTime&, xBeam&, xRing&);          // calculates rates for Molecular Dynamics model

   void Kick (xTime&, xBeam&, xRing&);                                        // kick from IBS
   void IndividualKick(xTime&, xBeam&, int j, vectorU, vectorU&, double k);  // kick for j-th particle
   void SidorinKick(xTime&, xBeam&, xRing&);
   //05.06
   void LocalKick(xTime&, xBeam&, xRing&);              //Local longitudinal diffusion
   void ZenkevichKick(xTime&, xBeam&, xRing&);         //Zenkevich kinetic model with kick in each optic element
   void ZenkevichKick2(xTime&, xBeam&, xRing&);         //Zenkevich kinetic model with averaged diffusion and friction
   void B_M_diffusion(xLattice&, xBeam&, xRing&);      //Diffusion components
   void B_M_friction(xLattice&, xBeam&, xRing&);      //Friction coefficients
};

extern xIBS iIBS;

#endif
