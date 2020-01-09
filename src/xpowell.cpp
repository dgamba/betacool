// Powell's Method
// http://www.library.cornell.edu/nr/cbookcpdf.html
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xPowell.h"
#include <math.h>
//#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include "bolideU.h"
//---------------------------------------------------------------------------
double *vector(int nl, int nh)
{
   double *v;
   v = (double *) malloc ((unsigned) (nh-nl+1) * sizeof(double));
   if (!v) Warning("allocation failure in vector()");
   return v - nl;
}

void free_vector(double *v, int nl, int nh)
{
   free ((char*) (v+nl));
}
//---------------------------------------------------------------------------

#define ITMAX100 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
/*
Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
a small number that protects against trying to achieve fractional accuracy for a minimum that
happens to be exactly zero.
*/
#define FMAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
/*
Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
the minimum to a fractional precision of about tol using Brent's method. The abscissa of
the minimum is returned as xmin, and the minimum function value is returned as brent, the
returned function value.
*/
{
   int iter;
   double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double d=0.0;
   double e=0.0;           // This will be the distance moved on the step before last.
   a=((ax < cx) ? ax : cx); // a and b must be in ascending order, but input abscissas need not be.
   b=((ax > cx) ? ax : cx);
   x=bx;
   w=bx;
   v=bx;              // Initializations...
   fw=(*f)(x);
   fv=fw;
   fx=fw;
   for (iter=1;iter<=ITMAX100;iter++) {    // Main program loop.
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) { // Test for done here.
         *xmin=x;
         return fx;
      }
      if (fabs(e) > tol1) {                 // Construct a trial parabolic fit.
         r=(x-w)*(fx-fv);
         q=(x-v)*(fx-fw);
         p=(x-v)*q-(x-w)*r;
         q=2.0*(q-r);
         if (q > 0.0) p = -p;
         q=fabs(q);
         etemp=e;
         e=d;
         if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
         // The above conditions determine the acceptability of the parabolic fit. Here we
         // take the golden section step into the larger of the two segments.
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
         else {
            d=p/q; // Take the parabolic step.
            u=x+d;
            if (u-a < tol2 || b-u < tol2)
               d=SIGN(tol1,xm-x);
         }
      } else {
         d=CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=(*f)(u);                    // This is the one function evaluation per iteration.
      if (fu <= fx) {                // Now decide what to do with our function evaluation.
         if (u >= x) a=x; else b=x;  // Housekeeping follows:
         SHFT(v,w,x,u)
         SHFT(fv,fw,fx,fu)
      } else {
         if (u < x) a=u; else b=u;
         if (fu <= fw || w == x) {
            v=w;
            w=u;
            fv=fw;
            fw=fu;
         } else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
         }
      } // Done with housekeeping. Back for another iteration.
   }
   Warning("Too many iterations in brent");
   *xmin=x; // Never get here.
   return fx;
}
//---------------------------------------------------------------------------

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY20 1.0e-20
/*
Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the
maximum magnification allowed for a parabolic-fit step.
*/
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double))
/*
Given a function func, and given distinct initial points ax and bx, this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa, fb, and fc.
*/
{
   double ulim,u,r,q,fu,dum;
   *fa=(*func)(*ax);
   *fb=(*func)(*bx);
   if (*fb > *fa) { // Switch roles of a and b so that we can go downhill in the direction from a to b.
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+GOLD*(*bx-*ax); // First guess for c.
   *fc=(*func)(*cx);
   while (*fb > *fc) {       // Keep returning here until we bracket.
      r=(*bx-*ax)*(*fb-*fc); // Compute u by parabolic extrapolation from a, b, c.
      q=(*bx-*cx)*(*fb-*fa); // TINY is used to prevent any possible division by zero.
      u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r) / (2.0*SIGN(FMAX(fabs(q-r),TINY20),q-r));
      ulim=(*bx)+GLIMIT*(*cx-*bx);
      // We won't go farther than this. Test various possibilities:
      if ((*bx-u)*(u-*cx) > 0.0) { // Parabolic u is between b and c: try it.
         fu=(*func)(u);
         if (fu < *fc) {           // Got a minimum between b and c.
            *ax=(*bx);
            *bx=u;
            *fa=(*fb);
            *fb=fu;
            return;
         } else if (fu > *fb) {    // Got a minimum between between a and u.
            *cx=u;
            *fc=fu;
            return;
         }
         u=(*cx)+GOLD*(*cx-*bx);    // Parabolic fit was no use. Use default magnification.
         fu=(*func)(u);
      } else if ((*cx-u)*(u-ulim) > 0.0) {     // Parabolic fit is between c and its allowed limit.
         fu=(*func)(u);
         if (fu < *fc) {
            SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
            SHFT(*fb,*fc,fu,(*func)(u))
         }
      } else if ((u-ulim)*(ulim-*cx) >= 0.0) { // Limit parabolic u to maximum allowed value.
         u=ulim;
         fu=(*func)(u);
      } else {                                // Reject parabolic u, use default magnification.
         u=(*cx)+GOLD*(*cx-*bx);
         fu=(*func)(u);
      }
      SHFT(*ax,*bx,*cx,u)                     // Eliminate oldest point and continue.
      SHFT(*fa,*fb,*fc,fu)
   }
}
//---------------------------------------------------------------------------

#define TOL 2.0e-4 // Tolerance passed to brent.
int ncom = 0; // Global variables communicate with f1dim.
double *pcom = 0,*xicom = 0,(*nrfunc)(double []);


double f1dim(double x) // Must accompany linmin.
{
   int j;
   double f,*xt;
   xt=vector(1,ncom);
   for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
   f=(*nrfunc)(xt);
   free_vector(xt,1,ncom);
   return f;
}

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
/*
Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
resets p to where the function func(p) takes on a minimum along the direction xi from p,
and replaces xi by the actual vector displacement that p was moved. Also returns as fret
the value of func at the returned location p. This is actually all accomplished by calling the
routines mnbrak and brent.
*/
{
   int j;
   double xx,xmin,fx,fb,fa,bx,ax;
   ncom=n; // Define the global variables.
   pcom=vector(1,n);
   xicom=vector(1,n);
   nrfunc=func;
   for (j=1;j<=n;j++) {
      pcom[j]=p[j];
      xicom[j]=xi[j];
   }
   ax=0.0; // Initial guess for brackets.
   xx=1.0;
   mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
   *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
   for (j=1;j<=n;j++) { // Construct the vector results to return.
      xi[j] *= xmin;
      p[j] += xi[j];
   }
   free_vector(xicom,1,n);
   free_vector(pcom,1,n);
}
//---------------------------------------------------------------------------

#define TINY 1.0e-25 // A small number.
#define ITMAX 200 // Maximum allowed iterations.
static double sqrarg;
#define SQR(a) (sqrarg = (a), sqrarg*sqrarg)

void powell(double p[],double**xi,int n,double ftol,int*iter,double*fret,double(*func)(double[]))
/*
Minimization of a function func of n variables. Input consists of an initial starting point
p[1..n]; an initial matrix xi[1..n][1..n], whose columns contain the initial set of directions
(usually the n unit vectors); and ftol, the fractional tolerance in the function value
such that failure to decrease by more than this amount on one iteration signals doneness. On
output, p is set to the best point found, xi is the then-current direction set, fret is the returned
function value at p, and iter is the number of iterations taken. The routine linmin is used.
*/
{
   
   int i,ibig,j;
   double del,fp,fptt,t,*pt,*ptt,*xit;
   pt=vector(1,n);
   ptt=vector(1,n);
   xit=vector(1,n);
   *fret=(*func)(p);
   for (j=1;j<=n;j++) pt[j]=p[j]; // Save the initial point.
   for (*iter=1;;(*iter)++) {
      fp=(*fret);
      ibig=0;
      del=0.0; // Will be the biggest function decrease.
      for (i=1;i<=n;i++) { // In each iteration, loop over all directions in the set.
         for (j=1;j<=n;j++) xit[j]=xi[j][i]; // Copy the direction,
         fptt=(*fret);
         linmin(p,xit,n,fret,func); // minimize along it,
         if (fabs(fptt-(*fret)) > del) { // and record it if it is the largest decrease so far.
            del=fabs(fptt-(*fret));
            ibig=i;
         }
      }
      if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
         free_vector(xit,1,n); // Termination criterion.
         free_vector(ptt,1,n);
         free_vector(pt,1,n);
         return;
      }
      if (*iter == ITMAX) Warning("powell exceeding maximum iterations.");
      for (j=1;j<=n;j++) {    // Construct the extrapolated point and the average direction moved. Save the old starting point.
         ptt[j]=2.0*p[j]-pt[j];
         xit[j]=p[j]-pt[j];
         pt[j]=p[j];
      }
      fptt=(*func)(ptt);      // Function value at extrapolated point.
      if (fptt < fp) {
         t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
         if (t < 0.0) {
            linmin(p,xit,n,fret,func); //Move to the minimum of the new direction, and save the new direction.
            for (j=1;j<=n;j++) xi[j][ibig]=xit[j];
         }
      }
   } // Back for another iteration.
}

#define DIMENSION 4

const double sgm = 5, dvs = 100;
double y[101];

double functest(double *p)
{
   double x, sum = 0;
   for (int i = 0; i < dvs + 1; i++)
   {  x = -sgm + (2 * sgm / dvs) * i;
      sum += pow(y[i] - p[3] * exp(-0.5 * pow(x / p[1], 2.)) -
                        p[4] * exp(-0.5 * pow(x / p[2], 2.)) , 2.);
   }
   return sum;
}

void powelltest()
{  double x;
   for (int j = 0; j < 101; j++)
   {  x = -sgm + (2 * sgm / dvs) * j;
      y[j] = 0.1* exp(-0.5 * pow(x / 2.0, 2.)) +
             0.9* exp(-0.5 * pow(x / 0.2, 2.));
   }

   double p[DIMENSION+1];
   p[0] = 0;
   p[1] = 1;
   p[2] = 1;
   p[3] = 0.5;
   p[4] = 0.5;

   double **xi;
   xi = new double*[DIMENSION+1];
   for (int i = 0; i < DIMENSION+1; i++)
   {  xi[i] = new double[DIMENSION+1];
      for (int j = 0; j < DIMENSION+1; j++)
      if (i == j) xi[i][j] = 1;
         else     xi[i][j] = 0;
   }
   int n = DIMENSION;
   double ftol = 1e-6;
   int iter = 0;
   double fret = 0;

   powell(p, xi, n, ftol, &iter, &fret, functest);

   for (int i1 = 0; i1 < DIMENSION+1; i1++)
      delete []xi[i1];
   delete []xi;

}
