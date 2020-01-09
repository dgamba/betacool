//---------------------------------------------------------------------------
//#include <vcl.h>
//#pragma hdrstop
#include "stdafx.h"
#include "vectorU.h"
//---------------------------------------------------------------------------

vectorU::vectorU()
{
   for (int i = 0; i < Size6; i++)
      ARRAY[i]._(0);
}
vectorU::vectorU(doubleU d1)
{  for (int i = 0; i < Size6; i++)
      ARRAY[i]._(d1);
}
vectorU::vectorU(double d1)
{  for (int i = 0; i < Size6; i++)
      ARRAY[i]._(d1, U1_);
}
vectorU::vectorU(double d1, Units_ u1)
{  for (int i = 0; i < Size6; i++)
      ARRAY[i]._(d1, u1);
}
vectorU::vectorU(Units_ u1)
{  for (int i = 0; i < Size6; i++)
      ARRAY[i]._(0, u1);
}
vectorU::vectorU(Units_ u1, Units_ u2)
{  for (int i = 0; i < Size6 / 2; i++)
   {	ARRAY[i*2]  ._(0, u1);
      ARRAY[i*2+1]._(0, u2);
   }
}
void vectorU::operator =(double d)
{
   for (int i = 0; i < Size6; i++)
      ARRAY[i] = d;
}
void vectorU::operator =(doubleU d)
{
   for (int i = 0; i < Size6; i++)
      ARRAY[i] = d;
}
void vectorU::operator =(const vectorU& v)
{
   for (int i = 0; i < Size6; i++)
      ARRAY[i]._(v.ARRAY[i]);
}

#undef  SIGN_i
#undef  SIGN_
#define SIGN_i +=
#define SIGN_  +
#include "vectorU.h"
#undef  SIGN_i
#undef  SIGN_
#define SIGN_i -=
#define SIGN_  -
#include "vectorU.h"
#undef  SIGN_i
#undef  SIGN_
#define SIGN_i *=
#define SIGN_  *
#include "vectorU.h"
#undef  SIGN_i
#undef  SIGN_
#define SIGN_i /=
#define SIGN_  /
#include "vectorU.h"  

//---------------------------------------------------------------------------

