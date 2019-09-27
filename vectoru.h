//---------------------------------------------------------------------------
#ifndef vectorUH
#define vectorUH

#include "doubleu.h"
#include "matrixu.h"
//---------------------------------------------------------------------------
class vectorU
{
  doubleU ARRAY[Size6];

 public:
  doubleU& operator[](int index)
  {
    if((index<0)||(Size6<=index))
       printf("vectorU index = ",index," out of range 0 - ",Size6);
    return ARRAY[index];
  }
  public:
  vectorU();
  vectorU(doubleU);
  vectorU(double);
  vectorU(double, Units_);
  vectorU(Units_);
  vectorU(Units_, Units_);

  void operator=(double);
  void operator+=(double);
  void operator-=(double);
  void operator*=(double);
  void operator/=(double);

  vectorU operator+(double);
  vectorU operator-(double);
  vectorU operator*(double);
  vectorU operator/(double);

  void operator=(doubleU);
  void operator+=(doubleU);
  void operator-=(doubleU);
  void operator*=(doubleU);
  void operator/=(doubleU);

  vectorU operator+(doubleU);
  vectorU operator-(doubleU);
  vectorU operator*(doubleU);
  vectorU operator/(doubleU);

  void operator=(const vectorU&);
  void operator+=(vectorU);
  void operator-=(vectorU);
  void operator*=(vectorU);
  void operator/=(vectorU);

  vectorU operator+(vectorU);
  vectorU operator-(vectorU);
  vectorU operator*(vectorU);
  vectorU operator/(vectorU);
};
//-------------------------------------------------------
#endif

#ifdef SIGN_i
void vectorU::operator SIGN_i(double m)
{
	for (int i = 0; i < Size6; i++)
		ARRAY[i] SIGN_i m;
}
vectorU vectorU::operator SIGN_(double m)
{  vectorU tmp;
   tmp = *this;
   tmp SIGN_i m;
	return tmp;
}
//---------------------------------------------------------------------------
void vectorU::operator SIGN_i(doubleU m)
{
	for (int i = 0; i < Size6; i++)
		ARRAY[i] SIGN_i m;
}
vectorU vectorU::operator SIGN_(doubleU m)
{  vectorU tmp;
   tmp = *this;
   tmp SIGN_i m;
	return tmp;
}
//---------------------------------------------------------------------------
void vectorU::operator SIGN_i(vectorU m)
{
	for (int i = 0; i < Size6; i++)
		ARRAY[i] SIGN_i m.ARRAY[i];
}
vectorU vectorU::operator SIGN_(vectorU m)
{  vectorU tmp;
   tmp = *this;
   tmp SIGN_i m;
	return tmp;
}   

#endif
