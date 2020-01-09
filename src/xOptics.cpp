//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xOptics.h"
#include "xLibrary.h"
//---------------------------------------------------------------------------
doubleU xLattice::B_s(m_);

xLattice::xLattice()
{  xLattice::B_s._(m_);
   dist._ (0, m_);
   betax._(0, m_);
   alfax._(0, 1);
   mux._  (0, 2*U_pi());
   gammax._(0, m_^-1);         
   x._    (0, m_3 * m_);
   px._   (0, 0.001);
   Dx._   (0, m_);
   Dpx._  (0, 1);
   Qx._   (0, 1);
   betay._(0, m_);
   alfay._(0, 1);
   muy._  (0, 2*U_pi());
   gammay._(0, m_^-1);
   y._    (0, m_3 * m_);
   py._   (0, 0.001);
   Dy._   (0, m_);
   Dpy._  (0, 1);
   Qy._   (0, 1);
}

//---------------------------------------------------------------------------
xOpticsLibrary::xOpticsLibrary()
{
   EL_FORCES = false;
   EL_FIELDS = false;
   EL_MATRIX = false;
   RE_MATRIX = false;
   entr = true;
   exit = true;
}

vectorU xOpticsLibrary::GetForces(xTime&t,U_Energy&e,vectorU ion)
{
	Warning("GetForces was not defined");
  	vectorU forces(U1_, m_^-1);
	return forces;
}

vectorU xOpticsLibrary::GetFields(xTime&t,U_Energy&e,vectorU ion)
{
	Warning("GetEfield was not defined");
  	vectorU fields(V_/m_, G_);
	return fields;
}

matrixU xOpticsLibrary::GetMatrix(xTime&t,U_Energy&e)
{
	Warning("GetMatrix was not defined");
  	matrixU matrix;
	return matrix;
}
//---------------------------------------------------------------------------

xRingOptics::xRingOptics()
{
   EL_LATTICE = false;
   RE_MATRIX = true;
}

matrixU& xRingOptics::GetMatrix(xTime&t,U_Energy&e)
{  if(RE_MATRIX)
   {  RE_MATRIX = false;
      bool drift = true;
      for (int i = 0; i < LGetSize(); i++)
      {  if ((*LList[i])->EL_MATRIX)
         {  Matrix = (*LList[i])->GetMatrix(t, e);
            drift = false;
            RE_MATRIX = (*LList[i])->RE_MATRIX;
            break;
         }
      }
      if (drift)
      {  iDrift.Length = Length;
         Matrix = iDrift.GetMatrix(t, e);
      }
   }
   return Matrix;
}


