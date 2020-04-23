#include "stdafx.h"
#include "pellets.h"
#include <time.h>
#include "math.h"

//=====================TIME===================================================//

pTime::pTime()
{
   t = 0;
   dt = 0;
   s = 0;
   ds = 0;
   Steps = 0;
}

void pTime::TimeStep(double h, double v)
{
   Velocity = v;
   dt = h;
   ds = dt * v;
}

void pTime::DistStep(double h, double v)
{
   Velocity = v;
   ds = h;
   dt = ds / v;
}

void pTime::operator+=(double h)
{
   t += dt * h;
   s += ds * h;
   Steps += (int)h;
}

//==========================COORDINATES=======================================//
Coordinates::Coordinates()
{
   for (int i = 0; i < 6; i++)
      X[i] = 0;
}

double &Coordinates::operator[](int i)
{
   return X[i];
}

void Coordinates::operator+=(Coordinates x0)
{
   Coordinates x1;
   for (int i = 0; i < 6; i++)
      X[i] += x0.X[i];
}

Coordinates Coordinates::operator+(Coordinates f2)
{
   Coordinates x2;
   for (int i = 0; i < 6; i++)
      x2.X[i] = X[i] + f2.X[i];
   return x2;
}

Coordinates Coordinates::operator*(double dt)
{
   Coordinates x0;
   for (int i = 0; i < 6; i++)
      x0[i] = X[i] * dt;
   return x0;
}

Coordinates Coordinates::operator/(double h)
{
   Coordinates x0;
   for (int i = 0; i < 6; i++)
      x0[i] = X[i] / h;
   return x0;
}

//=====================PELLETS================================================//
void Pellet::setPosition(pTime t)
{
   position[0] += t.dt * position[1];
   position[2] += t.dt * position[3];
   position[4] += t.dt * position[5];
}

void Pellet::set(Coordinates C0, Coordinates C1)
{
   double X, Y;

   do
   {
      X = (1. - rand() % 20000 / 10000.) * C1[0];
      Y = (1. - rand() % 20000 / 10000.) * C1[2];
   } while (X * X / (C1[0] * C1[0]) + Y * Y / (C1[2] * C1[2]) > 1.0);
   //   }while (sqrt( X*X*C1[0]*C1[0] + Y*Y*C1[2]*C1[2] )  C1[0]*C1[2]);

   position[0] = C0[0] + X;
   position[1] = C0[1] + (1. - rand() % 20000 / 10000.) * C1[1];
   position[2] = C0[2] + Y;
   position[3] = C0[3] + (1. - rand() % 20000 / 10000.) * C1[3];
   position[4] = C0[4] + (1. - rand() % 20000 / 10000.) * C1[4];
   position[5] = C0[5] + (1. - rand() % 20000 / 10000.) * C1[5];
}

//===========FLUX=========================================================//
Flux::~Flux()
{
   delete[] pellets;
}

Pellet &Flux::operator[](int i)
{
   return pellets[i];
}

void Flux::init(int N)
{
   Nplts = N;
   pellets = new Pellet[Nplts];
}

void Flux::set()
{

   for (int i = 0; i < Nplts; i++)
   {
      C0[4] = L * i;
      pellets[i].set(C0, C1);
   }
}

double Pellet::getP(double px[], int num, double X, double R)
{
   int index = (int)(num * X / (2. * R) + num / 2.);
   if (index < 0)
      return 0;
   if (index >= num)
      return 0;
   return px[index];
}

double Flux::getPxPy(Pellet &P, double px[], double ps[], int num, double Rx, double Rs)
{
   double pxpy, px1, py1;

   px1 = P.getP(px, num, P.position[0], Rx);
   py1 = P.getP(ps, num, P.position[4], Rs);
   pxpy = px1 * py1;

   return pxpy;
}

//=====================DYNAMICS===============================================//
int xPellet::OnGet()
{
   X0[0] = Data.OnGet(35, 1, 5);
   X0[1] = Data.OnGet(35, 2, 0);
   X0[2] = Data.OnGet(35, 3, 0);
   X0[3] = Data.OnGet(35, 4, 0);
   X0[4] = Data.OnGet(35, 5, 0);
   X0[5] = Data.OnGet(35, 6, -60000);

   X1[0] = Data.OnGet(35, 7, 4);
   X1[1] = Data.OnGet(35, 8, 0);
   X1[2] = Data.OnGet(35, 9, 4);
   X1[3] = Data.OnGet(35, 10, 0);
   X1[4] = Data.OnGet(35, 11, 2);
   X1[5] = Data.OnGet(35, 12, 0);

   pellets.L = Data.OnGet(35, 13, 4);
   pellets.Detector = Data.OnGet(35, 14, 1);
   pellets.C0 = X0;
   pellets.C1 = X1;

   step = Data.OnGet(35, 15, 1e-5);
   tt = Data.OnGet(35, 16, 0.0051);
   Nplts = Data.OnGetI(35, 17, 1000);
   //  num   =  Data.OnGet(35,18,400);   //Krestnikov 14.01.2009
   Rx = Data.OnGet(35, 19, 8);
   Rs = Data.OnGet(35, 20, 4);

   pellets.init(Nplts);
   t.TimeStep(step, X0[5]);

   pellets.set();

   return 0;
}

void xPellet::Run(double px[], double ps[])
{
   int i = 0;
   Average = 0;
   N = (int)(tt / t.dt);
   do
   {
      PxPy = 0;

      for (int j = 0; j < pellets.Nplts; j++)
      {
         PxPy += pellets.getPxPy(pellets[j], px, ps, num, Rx, Rs);
         pellets[j].setPosition(t);
      }
      Average += PxPy;
      t += 1;
      i++;
   } while (i < N);
   Average = Average / N;
}
