//---------------------------------------------------------------------------

#ifndef MatrixUH
#define MatrixUH
//---------------------------------------------------------------------------
#define Size6 6
#define complexU complex<double>

#if __BORLANDC__ == 1360
#include <complex.h>
#else
#include <complex>
#endif
#define _CRT_SECURE_NO_WARNINGS
//#include stdafx.h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "doubleU.h"
#include "bolideU.h"
#include "vectorU.h"

extern complex<double> J;

//____________________ARCCOS__&__ARCTAN_______________________
complexU arccos(complexU&);
complexU arctan(complexU&);
//---------------------------------------------------------------------------

template <class tin>
class Tine
{  int NUMBER;
   void Destructor()
   {  if (NUMBER)
      {  delete[] Line;
         NUMBER = 0;
         Number = 0;
      }
   }
 protected:
   tin* Line;
 public:
   int Number;
  ~Tine() { Destructor(); }
   Tine(int number = 0)
   {  NUMBER = 0;
      Line = NULL;
      SetNumber(number);
   }
   tin* operator() () { return Line; }
   tin& operator[] (int number)
   {
      if (number<0 || number>=NUMBER)
      {  //Warning("Tine index [",number,"] out of range");
         return *Line;
      }
      return Line[number];
   }
   void SetNumber(int number)
   {  if (number > NUMBER)
      {  tin *l = new tin[number];
         for (int i = 0; i < NUMBER; i++)
         if (i < number)
            l[i] = Line[i];
         Destructor();
         Line = l;
         NUMBER = number;
      }
      Number = number;
   }
   void operator=(Tine& tine)
   {  SetNumber(tine.Number);
      for(int n = 0; n <tine.Number; n++)
         Line[n] = tine.Line[n];
   }
};
//---------------------------------------------------------------------------

template <class tip>
class Tray
{  int NUMBER;
   int INDEX;
   void Destructor()
   {  if (NUMBER)
      {  for(int i = 0; i < NUMBER; i++)
            delete[] Array[i];
         delete[] Array;
         NUMBER = 0;
         Number = 0;
      }
   }
 protected:
   tip** Array;
 public:
   int Number;
   int Index;
  ~Tray() { Destructor(); }
   Tray(int number = 0, int index = 0)
   {  NUMBER = 0;
      INDEX = 0;
      Array = NULL;
      SetNumber(number, index);
   }
   tip**operator() () { return Array; }
   tip* operator[] (int number)
   {
      if (number<0 || number>=NUMBER)
      {  //Warning("Tray index [",number,"] out of range");
         return *Array;
      }
      return Array[number];
   }
   void SetNumber(int number, int index)
   {  if ((number > NUMBER) || (index > INDEX))
      {  tip **a = new tip*[number];
         for (int i = 0; i < number; i++)
            a[i] = new tip[index];
         for (int j = 0; j < NUMBER; j++)
         if (j < number)
            for (int k = 0; k < INDEX; k++)
            if (k < index)
               a[j][k] = Array[j][k];
         Destructor();
         Array = a;
         NUMBER = number;
         INDEX = index;
      }
      Number = number;
      Index = index;
   }
   void operator=(Tray& tray)
   {  SetNumber(tray.Number, tray.Index);
      for(int n = 0; n <tray.Number; n++)
      for(int i = 0; i <tray.Index;  i++)
         Array[n][i] = tray.Array[n][i];
   }
};
//---------------------------------------------------------------------------

class Matrix6
{public:
   int Row, Col;
   complexU Array6[Size6][Size6];

 public:
   Matrix6();
   Matrix6(int, int);
   complexU* operator[](int);
   void Diag(complexU, complexU);

   Matrix6 operator+(Matrix6);
   Matrix6 operator+(complexU);
   Matrix6 operator+(double);

   Matrix6 operator-(Matrix6);
   Matrix6 operator-(complexU);
   Matrix6 operator-(double);

   Matrix6 operator*(Matrix6);
   Matrix6 operator*(complexU);
   Matrix6 operator*(double);

   Matrix6 operator/(Matrix6);
   Matrix6 operator/(complexU);
   Matrix6 operator/(double);

//   Matrix6 operator =(const Matrix6&);

   void operator+=(Matrix6);
   void operator+=(complexU);
   void operator+=(double);

   void operator-=(Matrix6);
   void operator-=(complexU);
   void operator-=(double);

   void operator*=(Matrix6);
   void operator*=(complexU);
   void operator*=(double);

   void operator/=(Matrix6);
   void operator/=(complexU);
   void operator/=(double);
//----------------------------------------------------------------------------
   //int Min(int, int);
   void Simplectic();
   void SetZeroImag();
   Matrix6 Hat();
   Matrix6 Range(int, int, int, int);
  	Matrix6 Minor(int, int);
	Matrix6 Back();
	Matrix6 Trans();
	complexU Det();
   complexU Trace();
};
//---------------------------------------------------------------------
class matrixU : public Matrix6
{  Units_ un;
   doubleU vu[3];
 public:
   matrixU(int r=Size6, int c=Size6, Units_ u=m_);
   void operator=(const Matrix6&);
   void operator=(const matrixU&);
   //vectorU operator%(vectorU);
   void operator()(int, int, doubleU);
   doubleU operator()(int, int);
};
#endif
