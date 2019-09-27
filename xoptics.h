//---------------------------------------------------------------------------
#ifndef xOpticsH
#define xOpticsH
#include "xrunge.h"
#include "bolideu.h"
#include "doubleu.h"
//---------------------------------------------------------------------------

class xLattice
{
 public:
   doubleU dist;                                              // current long. coordinate of optic element (distance from the S0 point)
   doubleU betax, alfax, mux, x, px, Dx, Dpx, gammax, Qx;     // Lattice functions for x plane
   doubleU betay, alfay, muy, y, py, Dy, Dpy, gammay, Qy;     // Lattice functions for y plane
   static doubleU B_s;                                        // synchrotron function
   xLattice();
};

class xOptics
{
 public:
   doubleU Length;                                            // length of the optic element
   xOptics() { Length._(m_); }
   chars OpticName;                                           // label of the optic element
};

class xOpticsLibrary : public xOptics
{
 public:                    // below types of representatiom form of optic element
   bool EL_FORCES;          // optic element represented by vector of Forces (right parts)
   bool EL_FIELDS;          // optic element represented by vector of Electrical Field component
   bool EL_MATRIX;          // optic element represented by transformation matrix
   bool RE_MATRIX;          // optic element represented by transformation matrix (for elements, dependent on time)
   bool entr;               // entrance fringe field
   bool exit;               // exit fringe field
   //bool Symplectic;
   //int id;
   //BData Param;
   //chars Name;
   xOpticsLibrary();
   virtual vectorU GetForces (xTime&t,U_Energy&e,vectorU ion);       // template for function to return vector of forces
   virtual vectorU GetFields (xTime&t,U_Energy&e,vectorU ion);       // template for function to return vector of fields
   virtual matrixU GetMatrix (xTime&t,U_Energy&e);                   // template for function to return matrix of element
   virtual vectorU OnEnter(xTime&t,U_Energy&e,vectorU ion) { return ion; }        // template for function to return vector of forces at the element entrance
   virtual vectorU OnExit (xTime&t,U_Energy&e,vectorU ion) { return ion; }        // template for function to return vector of forces at the element exit
   virtual ~xOpticsLibrary() { ; }
};

class xRingOptics : public xOptics, public LTemplate<xOpticsLibrary*>
{public:
   bool EL_LATTICE;                                                  // if element of ring has a Lattice structure
   bool RE_MATRIX;                                                   // if element of ring has a transformation matrix

   //bool Symplectic;
   xLattice Lattice;                                                 // Latiice structure
   matrixU Matrix;                                                   // Transformation matrix
   xRingOptics();
   matrixU& GetMatrix (xTime&t,U_Energy&e);                          // returns transformation matrix
};
//---------------------------------------------------------------------------

#endif
