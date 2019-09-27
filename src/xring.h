//---------------------------------------------------------------------------
#ifndef xRingH
#define xRingH
#include "bolideu.h"
#include "xlibrary.h"
#include "xbucket.h"
//---------------------------------------------------------------------------
class xBeam;

class xRing : public LTemplate<xRingOptics>, public xData
{
   Units_ UKinetic;                         // keep Energy definition

 public:
   int Index;                               // index of optics element
   xLattice LATTICE;                        // ring Lattice structure
   int InjectionPoint;                      // definition of lattice functions
   U_Energy Energy;                         // structure of UEnergy
   doubleU TLife;                           // particle lifetime
   doubleU Circ;                            // ring circumeference
   doubleU Arc;                             //
   doubleU GammaTr;                         // gamma transition
   doubleU Imagenary;                       // imagenary of gamma transition
   doubleU TunesH;                          // horizontal tune
   doubleU TunesV;                          // vertical tune
   doubleU HromatH;                         // horizontal chromaticity
   doubleU HromatV;                         // vertical chromaticity
   doubleU AcceptH;                         // horizontal acceptance
   doubleU AcceptV;                         // vertical acceptance
   doubleU AcceptDp;
   doubleU AcceptS;
   doubleU AcceptHinj;                      // horizontal acceptance at injection point
   doubleU AcceptVinj;                      // verticalal acceptance at injection point
   doubleU BetaH;                           // mean horizontal beta-function
   doubleU BetaV;                           // mean vertical beta-function
   doubleU DispH;                           // mean horizontal ring dispersion
   doubleU DispV;                           // mean verticaltal ring dispersion
   doubleU Trev;                            // revolution period

   // mean ring parameters:
   doubleU Eta;                             // off momentum factor
   doubleU R_m;                             //mean radius
   // for RF
   doubleU h;                               //harmonic number
	doubleU V;                               //RF voltage
   bool NonLinear;                          // non linear RF
   doubleU h_low;                           //harmonic number
	doubleU V_low;                           //RF voltage
   //doubleU dPini;                           //intital momentum deviation
   //doubleU dPfin;                           //final momentum deviation
   //doubleU Tsweep;                          //deviation sweeping time
   bool Inductor;                           //Induction acceleration
   doubleU Vind;                            //Amplitude of inductor
	doubleU Q_s;                             //Synchrotron tune
   doubleU L_s;                             //Separatrix length
   doubleU Size_s;                          //Separatrix size
   //doubleU B_s;                           // synchrotron function
   BData MadLine;
   //20.11.06 For barrier bucket
   doubleU V_B;                             //Voltage amplitude
   doubleU T1;                              //RF duration in T0
   doubleU T2;                              //Gap duration in T0
   doubleU S1;                              //Kicker gap im meters
   doubleU S2;                              //Gap width in meters
   doubleU BarHeight;                       //Barrier height
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~12/07
   //bool Analytic;
   //bool BucketFile;
   xBucket Bucket;
   
   //Vacuum system parameters  
   doubleU a_m; 										//mean radius of vacuum chamber

   // Matrix Form
   matrixU Matrix;                          // (ring) transformation matrix
   int Matrix20;                            // switch to show ring or selected optics element matrix
   int Index20;                             // index of the element to visualize tr. matrix
   int Value20;                             // type of the matrix values appearence
   double Det20;                            // matrix determinant

   xBeam* pBeam;                                    // pointer to beam
	xRing(xBeam&);
   int OnGet();
   int OnRun();
   int OnSet();

   void Number(int n) { LSetSize(n); }               // set size of the ring elements array
   int  Number(){return LGetSize( ); }               // returns size of the ring elements array

   bool ReadMADfile(char* filename);                // if true - read lattices from the pointed MAD file
   bool InputMAD(char* filename);                   // if true - read lattices from the pointed input MAD file
   bool OutputMAD(char* FileName);                  // if true - read lattices from the pointed output MAD file
   //void InsertMAD(BData&);

   //int OpticsLib;
  // bool AutoSkip;                                   // flag to skip points in lattice appearence
   //bool RedLattice;                                 // flag to show reduced littice structure
   bool Autoskip;
   int LatticeFile;
   int ReduceExtend;
   doubleU extend;                          //step for Extend Lattice

   void ReduceLattice(char* filename);              // plots reduced lattice structure
   xOpticsLibrary* GetOptics(chars name);
   void ExtendLattice();

   vectorU GetForces(xTime&t, vectorU&i);           // returns vector with right parts (forces)
	vectorU Field2Force(vectorU&X, vectorU&f);       // calculates vector of right parts (forces) accordingly to fields

	xLattice& GetLattice(int n);                     // returns lattices of the selected element
	matrixU&  GetMatrix(xTime&t, int n);             // returns matrix of the selected element
	matrixU  RingMatrix(xTime&t, int n = 0);         // calculates ring transformation matrix
	matrixU  RingMatrix(xTime&t, int n1, int n2);    // calculates transformation matrix of the pointed part of the ring

	vectorU OnEnter(xTime&t, vectorU ion);
	vectorU OnExit (xTime&t, vectorU ion);

   xLattice CalcLattice(matrixU& TrMatr);               // calculates lattices from the Ring matrix
   xLattice CalcLattice(matrixU& TM, xLattice& L1);     // Twiss tracking
   matrixU CalcTransformMatrix (xLattice&, xLattice&);  // calculates matrix between pointed lattices

   void SetLattice();                                  // calculates lattices from the Ring matrix and after makes Twiss tracking for ring
};

extern xRing  iRing;

#endif

