//---------------------------------------------------------------------------
#include "stdafx.h"
#define _CRT_SECURE_NO_WARNINGS
#include "matrixU.h"
//---------------------------------------------------------------------------
complexU J(0, 1);

complexU arccos(complexU &rez1)
{
   return -J * log(rez1 + J * sqrt(1. - rez1 * rez1));
}

complexU arctan(complexU &rez1)
{
   return (-J / 2.) * log((1. + J * rez1) / (1. - J * rez1));
}
//----------------------------------------------------------------------------
Matrix6::Matrix6()
{
   Row = Size6;
   Col = Size6;
   Diag(1, 0);
}
Matrix6::Matrix6(int row, int col)
{
   if ((Size6 < row) || (Size6 < col) || (row < 0) || (col < 0))
      printf("Wrong Matrix6 size : ", row, " x ", col);
   Row = row;
   Col = col;
   Diag(1, 0);
}
complexU *Matrix6::operator[](int index)
{
   if ((index < 0) || (Size6 <= index))
      printf("Matrix6 index =", index, " out of range 0 - ", Size6);
   return Array6[index];
}
void Matrix6::Diag(complexU diag, complexU rest)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         if (i == j)
            Array6[i][j] = diag;
         else
            Array6[i][j] = rest;
}
//----------------------------------------------------------------------------
Matrix6 Matrix6::operator+(Matrix6 m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] + m1.Array6[i][j];

   return m2;
}
Matrix6 Matrix6::operator+(complexU m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] + m1;

   return m2;
}
Matrix6 Matrix6::operator+(double m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] + m1;

   return m2;
}
//----------------------------------------------------------------------------
Matrix6 Matrix6::operator-(Matrix6 m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] - m1.Array6[i][j];

   return m2;
}
Matrix6 Matrix6::operator-(complexU m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] - m1;

   return m2;
}
Matrix6 Matrix6::operator-(double m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] - m1;

   return m2;
}
//---------------------------------------------------------------------------
Matrix6 Matrix6::operator*(Matrix6 m1)
{
   Matrix6 tmp;
   tmp = *this;
   tmp *= m1;
   return tmp;
}
Matrix6 Matrix6::operator*(complexU m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] * m1;

   return m2;
}
Matrix6 Matrix6::operator*(double m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] * m1;

   return m2;
}
//-----------------------------------------------------------------------------
Matrix6 Matrix6::operator/(Matrix6 m1)
{
   Matrix6 tmp;
   tmp = *this;
   tmp *= m1.Back();
   return tmp;
}
Matrix6 Matrix6::operator/(complexU m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] / m1;

   return m2;
}
Matrix6 Matrix6::operator/(double m1)
{
   Matrix6 m2;

   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         m2.Array6[i][j] = Array6[i][j] / m1;

   return m2;
}
//-----------------------------------------------------------------------------
void Matrix6::operator+=(Matrix6 m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] += m1.Array6[i][j];
}
void Matrix6::operator+=(complexU m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] += m1;
}
void Matrix6::operator+=(double m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] += m1;
}
//-----------------------------------------------------------------------------
void Matrix6::operator-=(Matrix6 m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] -= m1.Array6[i][j];
}
void Matrix6::operator-=(complexU m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] -= m1;
}
void Matrix6::operator-=(double m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] -= m1;
}
//-----------------------------------------------------------------------------
void Matrix6::operator*=(Matrix6 m2)
{
   complexU tmp[6];
   for (int i = 0; i < Row; i++)
   {
      for (int j = 0; j < Col; j++)
      {
         tmp[j] = 0;
         for (int k = 0; k < Row; k++)
            tmp[j] += Array6[i][k] * m2[k][j];
      }
      for (int j1 = 0; j1 < Col; j1++)
         Array6[i][j1] = tmp[j1];
   }
}
void Matrix6::operator*=(complexU m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] *= m1;
}
void Matrix6::operator*=(double m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] *= m1;
}
//----------------------------------------------------------------------------
void Matrix6::operator/=(Matrix6 m1)
{
   Matrix6 tmp;
   tmp = *this;
   tmp *= m1.Back();
}
void Matrix6::operator/=(complexU m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] /= m1;
}
void Matrix6::operator/=(double m1)
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         Array6[i][j] /= m1;
}
//----------------------------------------------------------------------------
/*
int Matrix6::Min(int a, int b)
{
   if(a<b) return a;
   else return b;
}
*/
void Matrix6::Simplectic()
{
   Diag(0, 0);
   for (int i = 0; i < Row / 2; i++)
   {
      Array6[i * 2][i * 2 + 1] = 1;
      Array6[i * 2 + 1][i * 2] = -1;
   }
}
void Matrix6::SetZeroImag()
{
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         if (fabs(imag(Array6[i][j])) < 1e-10)
            Array6[i][j] = real(Array6[i][j]);
}
Matrix6 Matrix6::Hat()
{
   if ((Row != 2) || (Col != 2))
      printf("tMatrix can not to perform Hat()");
   Matrix6 tmp(2, 2);
   tmp[0][0] = Array6[1][1];
   tmp[0][1] = -Array6[0][1];
   tmp[1][0] = -Array6[1][0];
   tmp[1][1] = Array6[0][0];
   return tmp;
}
Matrix6 Matrix6::Range(int r0, int c0, int r1, int c1)
{
   Matrix6 tmp(r1 - r0, c1 - c0);
   for (int i = r0; i < r1; i++)
      for (int j = c0; j < c1; j++)
         tmp[i - r0][j - r0] = Array6[i][j];
   return tmp;
}
Matrix6 Matrix6::Minor(int r, int c)
{
   int j, i = 0;
   Matrix6 tmp(Row - 1, Col - 1);
   for (int l = 0; l < Row - 1; i++, l++)
   {
      if (i == r)
         i++;
      j = 0;
      for (int k = 0; k < Col - 1; j++, k++)
      {
         if (j == c)
            j++;
         tmp[l][k] = Array6[i][j];
      }
   }
   return tmp;
}
Matrix6 Matrix6::Back()
{
   int size = Min(Row, Col);
   Matrix6 m2(size, size);
   Matrix6 tmp(size, size);
   m2 = (*this).Trans();
   complexU det = Det();
   for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
         tmp[i][j] = (m2.Minor(i, j)).Det() * ((i + j) % 2 ? -1. : 1.) / det;
   return tmp;
}
Matrix6 Matrix6::Trans()
{
   Matrix6 tmp(Col, Row);
   for (int i = 0; i < Row; i++)
      for (int j = 0; j < Col; j++)
         tmp[j][i] = Array6[i][j];
   return tmp;
}
complexU Matrix6::Det()
{
   complexU D;
   if (Row != Col)
      return 0;
   switch (Row)
   {
   case 1:
      D = Array6[0][0];
      break;
   case 2:
      D = Array6[0][0] * Array6[1][1] -
          Array6[0][1] * Array6[1][0];
      break;
   case 3:
      D = Array6[0][0] * Array6[1][1] * Array6[2][2] +
          Array6[0][1] * Array6[1][2] * Array6[2][0] +
          Array6[0][2] * Array6[1][0] * Array6[2][1] -
          Array6[0][2] * Array6[1][1] * Array6[2][0] -
          Array6[0][0] * Array6[1][2] * Array6[2][1] -
          Array6[0][1] * Array6[1][0] * Array6[2][2];
      break;
   default:
      D = 0;
      for (int j = 0; j < Row; j++)
         D += Array6[0][j] * (Minor(0, j)).Det() * (j % 2 ? -1. : 1.);
      break;
   }
   return D;
}
complexU Matrix6::Trace()
{
   complexU trace = 0;
   int size = Min(Row, Col);
   for (int i = 0; i < size; i++)
      trace = trace + Array6[i][i];
   return trace;
}
//----------------------------------------------------------------------------
matrixU::matrixU(int r, int c, Units_ u) : Matrix6(r, c)
{
   vu[1]._(u);
   vu[2]._(u ^ -1);
}
void matrixU::operator=(const Matrix6 &m2)
{
   int row = Min(Row, m2.Row);
   int col = Min(Col, m2.Col);
   for (int i = 0; i < row; i++)
      for (int j = 0; j < col; j++)
         Array6[i][j] = m2.Array6[i][j];
}
void matrixU::operator=(const matrixU &m2)
{
   int row = Min(Row, m2.Row);
   int col = Min(Col, m2.Col);
   for (int i = 0; i < row; i++)
      for (int j = 0; j < col; j++)
         Array6[i][j] = m2.Array6[i][j];
}
void matrixU::operator()(int i, int j, doubleU du)
{
   int index = ((i + j) % 2 ? 1 : 0) * (i % 2 ? 2 : 1);
   vu[index] = du;
   (*this)[i][j] = vu[index]();
}
doubleU matrixU::operator()(int i, int j)
{
   int index = ((i + j) % 2 ? 1 : 0) * (i % 2 ? 2 : 1);
   vu[index].v = real((*this)[i][j]);
   return vu[index];
}

//----------------------------------------------------------------------------
void test(Tine<int> &t2)
{
   Tine<int> t1;
   t1 = t2;
   Tray<int> t3, t4;
   t3 = t4;
   t3[1][1] = t4[0][0];
}
