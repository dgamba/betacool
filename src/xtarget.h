//---------------------------------------------------------------------------
#ifndef xTargetH
#define xTargetH
#include "xeffect.h"
//---------------------------------------------------------------------------
class xMaterial                                  // class for parameters of target material
{
   public:
   int dens;
   int mat;
   double Percent;                               // percent ratio of the material in complex target
   doubleU Ro;                                   // target density gramm/cm^3
   doubleU l;                                    // length of the target
   doubleU N;                                    // target density atom/cm^3
   doubleU gr;
   doubleU A;                                    // target material mass number, A
   doubleU Z;                                    // target material charge number, A
   doubleU Rox;                                  // target cross-section density gramm/cm^2
   doubleU Nx;                                   // target cross-section density atom/cm^2
   doubleU Omega0;                               // Number of scattering
   doubleU ChiA;                                 // Screening angle

   doubleU Emax;                                  // maximum energy loss in 1 collision
   doubleU I;                                    // ionization potential
   //Urban model
   doubleU R;                                    //fraction of ionization
   doubleU Nex1;                                 //number of first level excitations
   doubleU Eex1;                                 //potential of first level excitation
   doubleU Nex2;                                 //number of second level excitations
   doubleU Eex2;                                 //potential of second level excitation
   doubleU Nion;                                 //number of ionizations
   doubleU GGG;                                  //relative energy loss in ionization
   //-----------
   doubleU Teta2;                                // angle of scattering (stripping angle)
   doubleU dE;                                   // ionization loss on scattering
   bool Use_dE;                                  // switch
   doubleU Estr2;                                // energy spectrum changing (energy stranglung)

   doubleU dPP2, dPB2;                           // koefficients for calculation the emittances deviation
   doubleU Hhor, Hver, dEhor, dEver;             // emittances deviations due to momentum spread

   void Density(int);                            // calculates material density
   void Bethe(U_Energy&);                        // calculates scattering on target
   void RMSANGLE(U_Energy&);                     // calculates rms angles after scattering
   void Deviations(U_Energy&, xLattice&);         // calculation of emittance and momentum deviations

   bool ElectronCapture;
   bool SingleScattering;
   bool NuclearReaction;
   bool InteractionEvents;
   doubleU CrossSection;
   doubleU Lifetime(xRing&, xBeam&, xLattice&);             // calculation beam lifetime due to interaction with target
   xMaterial();
};
//---------------------------------------------------------------------------

class xTarget : public xMaterial, public xEffect//, public xOpticsLibrary
{public:
   xTarget(xBeam&);                                     // target
   xBeam* pBeam;                                        // pointer on beam
   int itype;                                           // type of the target (gas cell, pellet target or fibber)
   //int imethod;                                         // type of the target calculation method (Effective dencity or Monte Carlo)
   //int Gaussian;                                        // Scattering in the pellet in accordance with Gaussian distribution
   int Urban;                                            // Scattering in the pellet with plural transverse and with excitation in longitudinal distribution
   int Plural;
   doubleU Cross;				                             //
   doubleU probability;                                 // probability for every particle in beam (needed when beam is generated)
   int num_MC;                                          // number of particles for Monte Carlo
   int CrossNumber;

   doubleU PShiftX;                                     // pellet shift from the axis in horizontal plane
   doubleU XWidth;                                      // Horizontal width of the pellet flux
   doubleU PelletX;                                     // Pellet horizontal size
   doubleU PelletY;                                     // Pellet vertical size
   doubleU Velosity;                                    // Pellet vertical velocity
   doubleU Interval;                                    // Interval between dropping pellets
   double  Ntar;                                        // events nubmer per model particle per second
   int OnGet();
   int OnSet();
   vectorU Rates(xTime&, xBeam&, xRing&);               // calculates rates for target effect
	//vectorU OnEnter(xTime&t,U_Energy&e,vectorU);         // calculates ion angles on target enter
   //vectorU MonteCarlo(xTime&, xBeam&, xRing&);          // calculates target as Monte Carlo
   //bool GetPellet(xTime&t, vectorU&ion);                // calculates pellet position
   //bool GetPelletA(xTime&t, vectorU&ion);               // calculates pellet position averaged over pellet motion
   doubleU PelletProbability(xBeam&beam);               // generate probability for every particle in the beam
   doubleU FiberProbability(xTime&,vectorU,xRing&);     //probability for given particle
   doubleU CylinderProbability(xTime&,vectorU,xRing&);   //probability for given particle
   doubleU PelletCross(xTime&t, vectorU ion,xRing&ring); //calculates number of ion crossing the target
   void Kick (xTime&, xBeam&, xRing&);
};

extern xTarget iTarget;

#endif
