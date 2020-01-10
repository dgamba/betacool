//---------------------------------------------------------------------------
#ifndef xLibraryH
#define xLibraryH
#include "xOptics.h"
//---------------------------------------------------------------------------
class xDrift : public xOpticsLibrary                         // optic element Drift
{public:
   xDrift();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU ion);        // returns rigth part of Drift
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of Drift
   virtual ~xDrift() { ; }
};
extern xDrift iDrift;

class xBend : public xOpticsLibrary                     // optic element Bend
{public:
   doubleU ANGLE;                                       // rotation angle
   doubleU E1;                                          // entrance koefficient
   doubleU E2;                                          // exit koefficient
   doubleU K1;                                          // quadrupole koefficient
   doubleU K2;                                          // sextupole koefficient
   doubleU HGAP1;                                       // gap
   doubleU HGAP2;                                       // gap
   doubleU FINT1;                                       // fringe integrale
   doubleU FINT2;                                       // fringe integrale
   xBend();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Bend
   vectorU OnEnter(xTime&t,U_Energy&e,vectorU);          // returns rigth part of Entrance of Bend
   vectorU OnExit (xTime&t,U_Energy&e,vectorU);          // returns rigth part of Exit of Bend
   matrixU GetMatrix (xTime&t,U_Energy&e);               // returns transformation matrix of Bend
   virtual ~xBend() { ; }
};

class xEBend : public xOpticsLibrary                      // optic element Electric Bend
{public:
   doubleU Ro;              // bending rsdius
   doubleU FB;              // fringe field length of magnetic field
   doubleU FE;              // fringe field length of electric field
   doubleU n;               // field index of dispestion-suppresser
   xEBend();
   //vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Bend
   //vectorU OnEnter(xTime&t,U_Energy&e,vectorU);          // returns rigth part of Entrance of Bend
   //vectorU OnExit (xTime&t,U_Energy&e,vectorU);          // returns rigth part of Exit of Bend
   matrixU GetMatrix (xTime&t,U_Energy&e);               // returns transformation matrix of Bend
   virtual ~xEBend() { ; }
};

class xFField : public xOpticsLibrary                      // optic element Electric Bend
{public:
   doubleU Ro;              // bending rsdius
   doubleU FB;              // fringe field length of magnetic field
   doubleU FE;              // fringe field length of electric field
   doubleU n;               // field index of dispestion-suppresser
   xFField();
   //vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Bend
   //vectorU OnEnter(xTime&t,U_Energy&e,vectorU);          // returns rigth part of Entrance of Bend
   //vectorU OnExit (xTime&t,U_Energy&e,vectorU);          // returns rigth part of Exit of Bend
   matrixU GetMatrix (xTime&t,U_Energy&e);               // returns transformation matrix of Bend
   virtual ~xFField() { ; }
};

class xQuadrupole : public xOpticsLibrary
{public:
   doubleU K1;                                               //quadrupole koefficient
   doubleU TILT;                                             // rotation angle of skew quadrupole
   xQuadrupole();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Quadrupole
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of Quadrupole
   virtual ~xQuadrupole() { ; }
};

class xSextupole : public xOpticsLibrary
{public:
   doubleU K2;                                               //sextupole koefficient
   xSextupole();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Sextupole
   virtual ~xSextupole() { ; }
};

class xSolenoid : public xOpticsLibrary
{public:
   doubleU Ks;                                               // solenoid strength
   xSolenoid();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Solenoid
   vectorU OnEnter(xTime&t,U_Energy&e,vectorU);             // returns rigth part of Entrance of Solenoid
   vectorU OnExit (xTime&t,U_Energy&e,vectorU);             // returns rigth part of Exit of Solenoid
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of Solenoid
   virtual ~xSolenoid() { ; }
};

class xRFcavity : public xOpticsLibrary
{public:
   doubleU V;                                                // RF voltage
   doubleU H;                                                // harmonic number
   xRFcavity();
   //vectorU GetForces(xTime&t,U_Energy&e,vectorU&ion);
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of RF Cavity
   virtual ~xRFcavity() { ; }
};

class SimpleECool : public xOpticsLibrary                    // Simple ECOOL effect
{public:
   doubleU K[3];                                             //vector of momenta variation per one turn
   doubleU Tapered;                                          //horizontal gradient of the electron velocity
   doubleU SwitchOn;                                         // switch for ecool effect
   SimpleECool();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Solenoid
   virtual ~SimpleECool() { ; }
};

class xPallas : public xOpticsLibrary                        // Pallas ring
{public:
   doubleU BendAngle;                                        // bending angle of one focusing period
   doubleU Udc;                                              // direct current voltage
   doubleU Urf;                                              // rf voltage
   doubleU omega;                                            // rf frequency
   doubleU radius;                                           // quadrupole aperture radius
   xPallas();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Pallas
   virtual ~xPallas() { ; }
};

class ConstFocusing : public xOpticsLibrary                  // Constant focusing element
{public:
   double LambdaL_N;                                         // lenear density
   ConstFocusing();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of the lens
   virtual ~ConstFocusing() { ; }
};

//---------------------------------------------------------------------------

class RadialField : public xOpticsLibrary
{public:
        doubleU Gradient;                                         // field gradient
   RadialField();
   vectorU GetFields(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Radial field
   virtual ~RadialField() { ; }
};

class ConstSolenoid : public xOpticsLibrary
{public:
        doubleU BFiled;                                          // solenoid field
   ConstSolenoid();
   vectorU GetFields(xTime&t,U_Energy&e,vectorU);       // returns rigth part of Const Solenoid
   virtual ~ConstSolenoid() { ; }
};

//************************************************************************
// optics elements for coupled motion
//************************************************************************

class xCM_Toroid : public xOpticsLibrary                     // coupled-motion toroid
{public:

   doubleU B;                                                // magnetic field
   doubleU R;                                                // toroid radius
   doubleU Kdis;                                             // dispersion koefficient in toroid
   xCM_Toroid();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Solenoid
   vectorU OnEnter(xTime&t,U_Energy&e,vectorU ion);             // returns rigth part of Entrance of Solenoid
   vectorU OnExit (xTime&t,U_Energy&e,vectorU ion);             // returns rigth part of Exit of Solenoid
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of Solenoid
   virtual ~xCM_Toroid() { ; }
};

class xCM_Solenoid : public xOpticsLibrary                   // coupled-motion stright solenoid
{public:

   doubleU B;                                                // magnetic field
   xCM_Solenoid();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Solenoid
   vectorU OnEnter(xTime&t,U_Energy&e,vectorU);             // returns rigth part of Entrance of Solenoid
   vectorU OnExit (xTime&t,U_Energy&e,vectorU);             // returns rigth part of Exit of Solenoid
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of Solenoid
   virtual ~xCM_Solenoid() { ; }
};

class xCM_Sol_Quad : public xOpticsLibrary                   // coupled-motion (stright solenoid + spiral quadrupole)
{public:
   doubleU B;                                                // magnetic field
   doubleU Wind;                                             // number of fotation steps
   doubleU K1;                                               // entrance koefficient
   doubleU K2;                                               // exit koefficient
   double Quad;                                            // ?????????  (sign of the Quad? )
   xCM_Sol_Quad();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Solenoid
   vectorU OnEnter(xTime&t,U_Energy&e,vectorU);             // returns rigth part of Entrance of Solenoid
   vectorU OnExit (xTime&t,U_Energy&e,vectorU);             // returns rigth part of Exit of Solenoid
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of Solenoid
   virtual ~xCM_Sol_Quad() { ; }
};

class xCM_Tor_Quad : public xOpticsLibrary                   // coupled-motion (toroid + spiral quadrupole)
{public:
   doubleU B;                                                // magnetic field
   doubleU Wind;                                             // number of fotation steps
   doubleU K1;                                               // entrance koefficient
   doubleU K2;                                               // exit koefficient
   double Quad;                                            // ?????????  (sign of the Quad? )
   xCM_Tor_Quad();
   vectorU GetForces(xTime&t,U_Energy&e,vectorU);        // returns rigth part of Solenoid
   vectorU OnEnter(xTime&t,U_Energy&e,vectorU);             // returns rigth part of Entrance of Solenoid
   vectorU OnExit (xTime&t,U_Energy&e,vectorU);             // returns rigth part of Exit of Solenoid
   matrixU GetMatrix (xTime&t,U_Energy&e);                   // returns transformation matrix of Solenoid
   virtual ~xCM_Tor_Quad() { ; }
};

//************************************************************************

//---------------------------------------------------------------------------

#endif

