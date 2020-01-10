//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xIBS.h"
#include "xDraw.h"
#include "bpIBS.h"
#include "bpData.h"
#include <omp.h>
#include <math.h>

//---------------------------------------------------------------------------
xIBS::xIBS()
{
	ID = 6;
   EffectName = "XIBS";
   stepmu = 10;
   stepnu = 10;
   stepz  = 10;
   JieWeiL._(U1_);
   BjorkenL = 20;
   AverageTrans = 0;
   coupledJie = false;
   NoDispersion = false;
   DispersionCoeff = true;
   IBSkickModel = 0;
   Slices = false;
   CoreNumber = 0;
   CoreDefinition = 1;
   InCoreDefine = 3;
   BUCKETnb = 2;
   s_x._(m_);
   s_xb._(m_);
   s_yb._(m_);
   s_s._(m_);
   D_._(m_);
   ri._(m_);
   Ao._(s_^-1);
   Mx._(m_);
   Mpx._(m_);
   M._(1, m_);
   D_xx._(s_^-1);
   D_yy._(s_^-1);
   D_zz._(s_^-1);
   D_zx._(s_^-1);
   D_zy._(s_^-1);
   F_x._(s_^-1);
   F_y._(s_^-1);
   F_z._(s_^-1);
}

int xIBS::OnGet()
{
   CoreDefinition = Data.OnGet(23,6,1);
   if (CoreDefinition <= 0)
      CoreDefinition = 1;
   InCoreDefine = Data.OnGetI(23,9,3);
   if (InCoreDefine <= 0)
      InCoreDefine = 1;
   if (InCoreDefine > 3)
      InCoreDefine = 3;

   Model = Data.OnGetI(42,1,2);
   Model2= Data.OnGetI(42,2,0);
   AverageTrans = Data.OnGet(42,3,0);
   stepmu = Data.OnGetI(42,4,20);
   stepnu = Data.OnGetI(42,5,20);
   stepz  = Data.OnGetI(42,6,20);
   BjorkenL = Data.OnGet(42,7,40);
   JieWeiL = Data.OnGet(42,8,20);
   coupledJie = Data.OnGetB(42,9,false);
   NoDispersion = Data.OnGetB(42,10,false);
   B_M_log = Data.OnGet(42,14,20);
   lam_limit = Data.OnGet(42,15,3);
   step_lam = Data.OnGetI(42,16,1000);
   BMHEA = Data.OnGetB(42,17,false);
   Slices = Data.OnGetB(42,18,false);
   bpStepSize = Data.OnGetI(42,19,100);
   bpStepNumber = Data.OnGetI(42,20,1000);
   bpSliceNumber = Data.OnGetI(42,21,10);

   TuneShift = Data.OnGetB(44,1,false);
   CellNum[0] = Data.OnGet(44,2,100);
   CellNum[1] = Data.OnGet(44,3,100);
   CellNum[2] = Data.OnGet(44,4,100);
   GridSize[0] = Data.OnGet(44,5,100);
   GridSize[1] = Data.OnGet(44,6,100);
   GridSize[2] = Data.OnGet(44,7,100);
   GridSize[2] = Data.OnGet(44,7,100);
   FieldComponent = Data.OnGetI(44,8,0);
   VerticalCell = Data.OnGetI(44,9,0);
   Plottuneshift = Data.OnGetB(44,11,false);
   Bumpfile.SetSeparators("\t,");
   Bumpfile.Load(Data.OnGetC(44,10,"NoFile"));
   TSminX = Data.OnGet(44,12,0);
   TSminY = Data.OnGet(44,13,0);

   IBSkickModel = Data.OnGetI(23,2,0);
   CoreNumCor = Data.OnGetB(27,1,false);
   CoreFWHM = Data.OnGetB(27,2,false);
   IBSloc = Data.OnGetI(27,3,200);
   Box_IBS = Data.OnGet(27,4,0.5);
   TranCoup = Data.OnGet(27,5,0.0);

	return 0;
}

int xIBS::OnRun()
{
   IBSkickModel = Data.OnGetI(23,2,0);

   CoreDefinition = Data.OnGet(23,6,1);
   if (CoreDefinition <= 0)
      CoreDefinition = 1;
   InCoreDefine = Data.OnGetI(23,9,3);
   if (InCoreDefine <= 0)
      InCoreDefine = 1;
   if (InCoreDefine > 3)
      InCoreDefine = 3;
   CoreNumCor = Data.OnGetB(27,1,false);
   CoreFWHM = Data.OnGetB(27,2,false);

   bpStepSize = Data.OnGetI(42,19,100);
   bpStepNumber = Data.OnGetI(42,20,1000);
   bpSliceNumber = Data.OnGetI(42,21,10);
   TSminX = Data.OnGet(44,12,0);
   TSminY = Data.OnGet(44,13,0);

  	return 0;
}
//--------------
int xIBS::OnSet()
{
   Data.OnSet(23,6,CoreDefinition);
   Data.OnSet(23,9,InCoreDefine);
   Data.OnSet(42,12,Mx);
   Data.OnSet(42,13,Mpx);
	return 0;
}

void xIBS::Sigma(xLattice& lat, xBeam& beam, xRing& ring)
{
   s_xb = (beam.Emit[0] * lat.betax)^0.5;
   s_yb = (beam.Emit[1] * lat.betay)^0.5;

   s_x_=((1+(lat.alfax^2))*beam.Emit[0]/(lat.betax))^0.5;
   s_y_=((1+(lat.alfay^2))*beam.Emit[1]/(lat.betay))^0.5;
   s_p = beam.Emit[2]^0.5;

   s_x = ( (s_xb*s_xb) + (lat.Dx*lat.Dx*s_p*s_p) )^0.5;
   s_y = s_p * s_xb / (ring.Energy.Gamma * s_x);

   s_h = 1 / (((1/(s_p*s_p)) + ((lat.Dx /s_xb)^2))^0.5);
   D_  = (lat.alfax * lat.Dx) + (lat.betax * lat.Dpx);
   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;

   q = 2*ring.Energy.Beta*ring.Energy.Gamma*((s_yb/ri)^0.5);
   a = s_y / s_x_ * ((lat.alfax*lat.alfax+1)^0.5);
   b = s_y / s_y_ * ((lat.alfay*lat.alfay+1)^0.5);
   c = q * s_y;
   d = (lat.Dx * s_p) / (((s_xb*s_xb)+(lat.Dx*lat.Dx*s_p*s_p))^0.5);
   d_= s_p * D_ / s_x;

   if (beam.benum == BUNCHED)
   {  nb = 1;
      beam.CalcBunch();
      if (Slices)
         s_s = ring.L_s / (2 * (U_pi^0.5));
      else
         s_s = beam.s_s;
   }else
   {  s_s = ring.Circ / (2 * (U_pi^0.5));
      if (beam.benum == COASTING)
         nb = 2;
      else
         nb = BUCKETnb;
   }

   Ao = (((lat.alfax*lat.alfax+1)*(lat.alfay*lat.alfay+1))^0.5) *
        U_c * ri * ri * beam.Emit[3] /
        ( 32 * U_pi * U_pi * s_xb * s_x_ * s_yb * s_y_* s_s * s_p *
        ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);
}

// ============================== IBS Rates ====================================
vectorU xIBS::Rates(xTime& time, xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   vectorU average(m_/s_);
   vectorU Vect(m_, U1_);
   double ave0 = 0;
   double ave1 = 0;
   double ave2 = 0;
   Sigma(Lattice, beam, ring);
																					
	bpRates Growth;																									
	bpBeam bplusBeam;
	bpRing bplusRing;
	//int row = loader("nc2splt2.prn",'\t');
	int i;
   switch (Model)
   {case 0:
      rates = Piwinski(beam, ring);
      break;

    case 1:
    case 2:
      //#pragma omp parallel for reduction(+: ave0, ave1, ave2) - not ok
      for (i = 0; i < ring.Number(); i++)
      {
         vectorU ave(m_/s_);
         Sigma(ring.GetLattice(i), beam, ring);

         if (Model == 1)
            ave = JieWei(beam, ring) * ring[i].Length;
         if (Model == 2)
            ave = Martini() * ring[i].Length;
         ave0 += ave[0](m_/s_);
         ave1 += ave[1](m_/s_);
         ave2 += ave[2](m_/s_);

        }
      rates[0] = ave0 / ring.Arc(m_);
      rates[1] = ave1 / ring.Arc(m_);
      rates[2] = ave2 / ring.Arc(m_);
      break;

    case 5:
      #pragma omp parallel for reduction(+: ave0, ave1, ave2)
      for (i = 0; i < ring.Number(); i++) {
         vectorU ave(m_/s_);
         ave = Bjorken_Mtingwa2(ring.GetLattice(i), beam, ring) * ring[i].Length;
         ave0 += ave[0](m_/s_);
         ave1 += ave[1](m_/s_);
         ave2 += ave[2](m_/s_);
      }
      rates[0] = ave0 / ring.Arc(m_);
      rates[1] = ave1 / ring.Arc(m_);
      rates[2] = ave2 / ring.Arc(m_);
     break;

    case 3:
      beam.RMS(beam.Emit,Lattice);
	  
	  bplusBeam.emit.emittancex = beam.Emit[0](m_);
	  bplusBeam.emit.emittancey = beam.Emit[1](m_);
	  bplusBeam.emit.momentum2 = beam.Emit[2]();
	  bplusBeam.emit.sigmab = beam.s_s(m_);
	  bplusBeam.q0 = ring.Energy.Z();
	  bplusBeam.m0 = ring.Energy.A();
	  bplusBeam.energy = ring.Energy.Kinetic(M_6*eV_);
	  bplusBeam.nion = beam.Emit[3]();
	  bplusBeam.bunched = beam.benum;
     bplusBeam.Init();

	  bplusRing.Circ = ring.Circ(m_);
     bplusRing.GammaTr = ring.GammaTr();
     bplusRing.Imagenery = ring.Imagenary();                          
     bplusRing.Q_h = ring.TunesH();
     bplusRing.Q_v = ring.TunesV();
     bplusRing.Hrf = ring.h();
     bplusRing.Vrf = ring.V();
     bplusRing.set(ring.Number());
	  for (int l = 0; l < bplusRing.size; l++)
	  {
		  bplusRing[l].s = ring[l].Lattice.dist(m_);
		  bplusRing[l].ds = ring[l].Length(m_);
		  bplusRing[l].betax = ring[l].Lattice.betax(m_);
		  bplusRing[l].alphax = ring[l].Lattice.alfax();
		  bplusRing[l].betay = ring[l].Lattice.betay(m_);
		  bplusRing[l].alphay = ring[l].Lattice.alfay();
		  bplusRing[l].Dx = ring[l].Lattice.Dx(m_);
		  bplusRing[l].Dpx = ring[l].Lattice.Dpx();
		  bplusRing[l].Dy = ring[l].Lattice.Dy(m_);
		  bplusRing[l].Dpy = ring[l].Lattice.Dpy();
	  }
     bplusRing.Init(bplusBeam);
      bplusRing.lattice.s = ring.LATTICE.dist(m_);
      bplusRing.lattice.ds = ring.Circ(m_);
      bplusRing.lattice.betax = ring.LATTICE.betax(m_);
      bplusRing.lattice.alphax = ring.LATTICE.alfax();
      bplusRing.lattice.betay = ring.LATTICE.betay(m_);
      bplusRing.lattice.alphay = ring.LATTICE.alfay();
      bplusRing.lattice.Dx = ring.LATTICE.Dx(m_);
      bplusRing.lattice.Dpx = ring.LATTICE.Dpx();
      bplusRing.lattice.Dy = ring.LATTICE.Dy(m_);
      bplusRing.lattice.Dpy = ring.LATTICE.Dpy();

     bpIntegrateAcc(bpStepSize, bpStepNumber);
	  bpIBSrates(bplusRing, Growth, bplusBeam, bplusBeam.emit);
	  rates[0]=Growth.ratex;
	  rates[1]=Growth.ratey;
	  rates[2]=Growth.ratez;
	  
      break;

    case 4:
      rates = Gas_Relaxation(beam, ring);
      break;

    default:
      Warning("IBS Rates are not defined", PressEnter);
   }

   if (AverageTrans)
   {  vectorU r = rates;
      rates[0] = r[0]/2 * (2.-AverageTrans) +
                 r[1]/2 *     AverageTrans  * beam.Emit[1] / beam.Emit[0];
      rates[1] = r[1]/2 * (2.-AverageTrans) +
                 r[0]/2 *     AverageTrans  * beam.Emit[0] / beam.Emit[1];

      //rates[0] = (rates[0] + rates[1]) / 2;
      //rates[1] = rates[0];
   }
   return rates;
}

// ============================== Piwinski =====================================

double xIBS::PiwinskiIntegral(double a, double b, double c)
{
   double x = 0;
   double p, q;
   double f = 0;
   for (int i = 0; i < stepmu; i++)
   {  p = sqrt(a*a + x*x*(1. - a*a));
      q = sqrt(b*b + x*x*(1. - b*b));
      f += (log(0.5*c*c*(1./p +1./q)) - 0.577)*
           (1. - 3.*x*x)/(p*q) / stepmu;
      x += 1. / stepmu;
   }//while(x < 1.);
   f *= 8. * U_pi();
   return f;
}

vectorU xIBS::Piwinski(xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   rates[0] = Ao *(PiwinskiIntegral(1./a(),b()/a(),c()/a()) +
      (PiwinskiIntegral(a(),b(),c()) * ((ring.DispH*s_p/s_xb)^2)));
   rates[1] = Ao * PiwinskiIntegral(1./b(),a()/b(),c()/b());
   rates[2] = Ao * PiwinskiIntegral(a(),b(),c()) * ((s_h/s_p)^2) * nb;
   return rates;
}

// ============================== Jie Wei ======================================

doubleU xIBS::Hi_function(doubleU hi)
{
 doubleU result(0);
 doubleU I(0);
 doubleU h1(0);
 doubleU h2(0);
 doubleU arth(0);

 if (hi >= 1)
  {
   h1 = (((hi-1)/hi)^0.5);
   arth = 0.5*U_Ln((1+h1)/(1-h1));
   I = arth/((hi*(hi-1))^0.5);
  }
 if (hi < 1)
  {
   h2 = (((1-hi)/hi)^0.5);
   I = U_Atan(h2)/((hi*(1-hi))^0.5);
  }
 result = ((1 + (2*hi))*I- 3)/(1.-hi);

 return result;

}

vectorU xIBS::JieWei(xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);

   if (coupledJie)
   {  if (NoDispersion)
         d = nb * 2. / (1. + (nb * 2.));
      rates[0] = 2*4*U_pi*Ao*JieWeiL*(d / 2.)*(ring.Energy.Gamma * U_pi)/(2. * ring.GammaTr);
      rates[1] = rates[0];
      rates[2] = 2*4*U_pi*Ao*JieWeiL* nb*(1. - (d*d))*(ring.Energy.Gamma * U_pi)/(2. * ring.GammaTr * d);
   }else
   {
      doubleU hi;

      hi = Hi_function(((a*a) + (b*b))/2.);

      rates[0] = 2*4*U_pi*Ao*JieWeiL*hi*((-a*a/2.) + (d*d));
      rates[1] = 2*4*U_pi*Ao*JieWeiL*hi* (-b*b)/2.;
      rates[2] = 2*4*U_pi*Ao*JieWeiL*hi*nb*(1. - (d*d));
   }
   return rates;
}

// ============================== Martini ======================================

// integral_0^\infinity exp(-Dz) log(1+z^2) dz
// approximation by Abramowitz & Stegun 5.2.39
// d = 2.0 / D * ((0.5 * M_PI - gsl_sf_Si(D)) * sin(D) - gsl_sf_Ci(D) * cos(D));
// d = 40.0 / D; // This gives the value which is the same as one by B-M model

double xIBS::explogint(double D)
{
  double d;
  if(D>=1){
    const double a1 = 42.242855;
    const double a2 = 302.757865;
    const double a3 = 352.018498;
    const double a4 = 21.821899;
    const double b1 = 48.196927;
    const double b2 = 482.485984;
    const double b3 = 1114.978885;
    const double b4 = 449.690326;
    double D2I = 1.0 / (D * D);
    double num = (((a4 * D2I + a3) * D2I + a2) * D2I + a1) * D2I + 1.0;
    double den = (((b4 * D2I + b3) * D2I + b2) * D2I + b1) * D2I + 1.0;
    d = 2.0 * D2I * num / (D * den);
  }else{
    const double EulerGamma = 0.57721566490153286061;
    const double am1 = -2.0 * EulerGamma;
    const double bm1 = -2.0;
    const double a0 = U_pi();
    const double a1 = -1.5 + EulerGamma;
    const double b1 = 1.0;
    const double a2 = -U_pi() / 6.0;
    const double a3 = (25.0 - 12.0 * EulerGamma) / 144.0;
    const double b3 = -12.0 / 144.0;
    const double a4 = U_pi() / 120.0;
    const double a5 = (-49.0 + 20.0 * EulerGamma) / 7200.0;
    const double b5 = 20.0 / 7200.0;
    const double a6 = -U_pi() / 5040.0 * 0.818;
    double x =  (((((a6 * D + a5) * D + a4) * D + a3) * D + a2) * D + a1) * D
                + a0 + am1 / D;
    double D2 = D * D;
    double y = log(D) * (((b5 * D2 + b3) * D2 + b1) * D + bm1 / D);
    d = x + y;
  }
  return(d);
}

#define HIROSHI // Integration derived by Hiroshi Tsutsui (Japan)

vectorU xIBS::MartiniIntegral(double a, double b, double c, double d, double d_)
{
    double f1 = 0, f2 = 0, f3 = 0;

//#ifdef HIROSHI
   double hmu = 6.5 / stepmu;
   double hnu = 6.5 / stepnu;
//find minimum Dmn position
    double mu0 = 0;
    double nu0 = 0;
    double Dmin = b * b;
    for(int j = 0; j <= stepnu; j++){
      double nu = M_PI * j / stepnu;
      double sinnu = sin(nu);
      double cosnu = cos(nu);
      double Dmn = sinnu * sinnu + pow(a * cosnu - d_ * sinnu, 2.0);
      if(Dmn < Dmin){
        Dmin = Dmn;
        mu0 = 0.5 * M_PI;
        nu0 = nu;
      }
    }
//#pragma omp parallel for reduction(+:f1)// reduction(+: f2) reduction(+: f3)
   for(int i = 0; i < stepmu; i++)
   {  
      double coef = 1.;
      double sin_mu,cos_mu,sin_nu,cos_nu;
      double sin_mu_2,cos_mu_2,sin_nu_2,cos_nu_2;
      double g1, g2, g3, DCAP, int_over_z;
      int k1;
      
      double mu = 0.5 * M_PI * (i + 0.5) / stepmu;
      double sinhyp = sinh((i - stepmu * 0.5) * hmu);
      double coshyp = sqrt(1.0 + sinhyp * sinhyp);
      double sinhyp2 = sinh(M_PI * 0.5 * sinhyp);
      double coshyp22 = 1.0 + sinhyp2 * sinhyp2;
      double tanhyp2 = sinhyp2 / sqrt(coshyp22);
      mu = (tanhyp2 + 1.0) * M_PI * 0.25 + mu0;
      double coef1 = stepmu * hmu * 0.25 * M_PI * coshyp / (coshyp22);
      sin_mu = sin(mu);
      cos_mu = cos(mu);
      sin_mu_2 = sin_mu*sin_mu;
      cos_mu_2 = cos_mu*cos_mu;
      for(int j = 0; j < stepnu; j++)
      {  double nu = M_PI * (j + 0.5) / stepnu;
         double sinhyp = sinh((j - stepnu * 0.5) * hnu);
         double coshyp = sqrt(1.0 + sinhyp * sinhyp);
         double sinhyp2 = sinh(M_PI * 0.5 * sinhyp);
         double coshyp22 = 1.0 + sinhyp2 * sinhyp2;
         double tanhyp2 = sinhyp2 / sqrt(coshyp22);
         nu = (tanhyp2 + 1.0) * M_PI * 0.5 + nu0;
         double coef2 = stepnu * hnu * M_PI * 0.25 * coshyp / (coshyp22);
         coef = coef1 * coef2;
         sin_nu = sin(nu);
         cos_nu = cos(nu);
         sin_nu_2 = sin_nu * sin_nu;
         cos_nu_2 = cos_nu * cos_nu;
         DCAP = ((sin_mu_2*sin_nu_2)+(sin_mu_2*((a*cos_nu)-(d_*sin_nu))*
                ((a*cos_nu)-(d_*sin_nu)))+(b*b*cos_mu_2))/(c*c);

/* #else
   for (int i = 0; i < stepnu; i++)
   {
      double nu = 2. * U_pi() * i / stepnu;
      sin_nu = sin(nu);
      cos_nu = cos(nu);
      sin_nu_2 = sin_nu * sin_nu;
      cos_nu_2 = cos_nu * cos_nu;
      for (int j = 0; j < stepmu; j++)
      {
         double mu = U_pi() * j / stepmu;
         sin_mu = sin(mu);
         cos_mu = cos(mu);
         sin_mu_2 = sin_mu*sin_mu;
         cos_mu_2 = cos_mu*cos_mu;
         DCAP = ((sin_mu_2*cos_nu_2)+(sin_mu_2*((a*sin_nu)-(d_*cos_nu))*
                ((a*sin_nu)-(d_*cos_nu)))+(b*b*cos_mu_2))/(c*c);
#endif */
         switch(Model2)
         {case 0:
            int_over_z = BjorkenL / DCAP;
            break;
          case 1:
            int_over_z = explogint(DCAP);
            break;
          case 2:
            int_over_z = 0;
            double z, dz;
            dz = 20. / (DCAP * stepz);
            for (k1 = 0; k1 < stepz; k1++)
            {
               z = dz * k1;
               int_over_z += exp(-DCAP*z) * log(1.+z*z) * dz;
            }
            break;
          case 3:
            int_over_z = -2.0*((0.577 + log(DCAP))/DCAP + 1.0);
            break;
         }
//#ifdef HIROSHI
         g1 = 1.0 - (3.0*sin_mu_2*sin_nu_2);
         g2 = 1.0 - (3.0*sin_mu_2*cos_nu_2) + (6.0*d_*sin_mu*sin_nu*cos_nu/a);
         g3 = 1.0 - (3.0*cos_mu_2);
/* #else
         g1 = 1.0 - (3.0*sin_mu_2*cos_nu_2);
         g2 = 1.0 - (3.0*sin_mu_2*sin_nu_2) + (6.0*d_*sin_mu*sin_nu*cos_nu/a);
         g3 = 1.0 - (3.0*cos_mu_2);
#endif */
         int_over_z *= sin_mu * coef;
         f1 += int_over_z * g1;
         f2 += int_over_z * g2;
         f3 += int_over_z * g3;
      }
   }
   doubleU k (2*Ao*2*U_pi*U_pi/(stepmu*stepnu*c*c));
   vectorU f(s_^-1);
   f[2] = k * nb *         f1*(1    -(d*d))/ 2;
   f[0] = k *((f2*a*a) + f1*(d_*d_+(d*d))) / 2;
   f[1] = k *  f3*b*b                      / 2;
   return f;
}

vectorU xIBS::Martini()
{
   return MartiniIntegral(a(),b(),c(),d(),d_());
}

// ============================== Detailed =====================================

doubleU xIBS::Bessel_I0(doubleU x1) // Returns the modified Bessel function I0(x) for any x
{
 doubleU result(0);
 double x = x1();
 double ax,ans;
 double y;        // Accumulate polynomials in double precision.
  if ((ax=fabs(x)) < 3.75)
  {                         // Polynomial fit.
   y = x/3.75;
   y *= y;
   ans = 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
   +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  }
  else
  {
   y = 3.75/ax;
   ans = (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
   +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
   +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
   +y*0.392377e-2))))))));
  }
 result = ans;
 return result;
}

doubleU xIBS::G_F(doubleU x, doubleU y) // Here for integrals the same stepnu and stepmu from Martini integral are used
{
   doubleU result(0);
   doubleU K1;
   double infinity1 = 6.;
   double infinity2 = 6.;

   K1 = (2./U_pi)*U_Exp(-1.*((x*x) + (y*y))/2.);        ///  Pi was under 2 !

   for (int i = 0; i < stepmu; i++)
   {  double nu = infinity1 * i / stepmu;
      for (int j = 0; j < stepnu; j++)
      {  double mu = infinity2 * j / stepnu;
         result += exp(-1.*((nu*nu) + (mu*mu))/2.)*Bessel_I0(nu*x)*Bessel_I0(mu*y);
      }
   }

   result *= K1*(infinity1/stepmu)*(infinity2/stepnu);
   return result;
}


vectorU xIBS::Burov(xLattice& lat, vectorU X, xRing& ring)
{
   vectorU Inv;
   Inv[0]._(m_);
   Inv[1]._(m_);

   vectorU rates(s_^-1);    // result vector

   doubleU unit(U1_);
   doubleU ri(0,m_);
   doubleU xm(0,m_), zm(0,m_), x2m(0,m_^2),ym(0,m_);
   doubleU sigx(0,m_), sigz(0,m_), sig2x(0,m_^2), sig2z(0,m_^2), sigy(0,m_), sig2y(0,m_^2);
   doubleU Vxm, Vym, Vzm, Ux, Ux2, Uy;
   doubleU Ni, Li;
   doubleU K1(0,s_^-1), K2, Kx, Ky, K3, K4;


   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;

   //Beam parameters
   sigx = iBeam.Sigma[0];        // Horizontal beam rms dimension
   sigy = iBeam.Sigma[2];        // Horizontal beam rms dimension
   sigz = iBeam.Sigma[4];        // Longitudinal beam rms dimension
   sig2x = iBeam.Sigma_2[0];     // Horizontal beam rms dimension squared
   sig2y = iBeam.Sigma_2[2];
   sig2z = iBeam.Sigma_2[4];     // Longitudinal beam rms dimension squared
   Ux =  iBeam.Sigma[1];         // r.m.s. longitudinal ion velocity
   Ux2 = iBeam.Sigma_2[1];;
   Uy = iBeam.Sigma[3];          // r.m.s. longitudinal ion velocity squared

   Ni = iBeam.Emit[3];           // number of ions
   Li = 19.;                     // Coulomb logarithm (ion)


   //Particle parameters
   Inv = xEcool::Courant_Snyder(X,lat,lat.B_s,true);
// navernoe nado tak:  iBeam.Courant_Snyder(lat);
// no togda mnogo nado peredelyvat' v etoi funkcii.


   xm = (2.* Inv[0] * lat.betax)^0.5;   // Amplitude Xmax
   ym = (2.* Inv[1] * lat.betay)^0.5;
   zm = ((2.*Inv[2])^0.5)* lat.B_s;   // Ampltude Zmax
   x2m = xm*xm;                        // Amplitude Xmax squared
   Vxm = (2.*(1.+(lat.alfax^2))*Inv[0]/(lat.betax))^0.5;   //Amplitude Xmax prime (Vxm)
   Vym = (2.*(1.+(lat.alfay^2))*Inv[1]/(lat.betay))^0.5;   //Amplitude Ymax prime (Vym)
   Vzm = (2.*Inv[2])^0.5;                                  //Amplitude Zmax prime (Vzm)


   //Debug version
   Vzm = 1.0;

   K1 = ((2./U_pi)^0.5)*(Ni*ri*ri*U_c*Li)/
        (ring.Energy.Gamma2*ring.Energy.Gamma2*ring.Energy.Beta2*ring.Energy.Beta*Vzm*Vzm*sig2x*sigz*Ux);

   K2 = (zm*zm)/(4.*sig2z);
   Kx = Vxm/Ux;
   Ky = Vym/Uy;
   K3 = 1. - (( Vzm/(ring.Energy.Gamma*Ux*1.41))^0.5);
   K4 = ((lat.Dx*lat.Dx)+((lat.Dpx*lat.betax) + (lat.alfax*lat.Dx))*
                         ((lat.Dpx*lat.betax) + (lat.alfax*lat.Dx)))*Vzm*Vzm/x2m;
   /*
   rates[2] = K1*U_Exp(-K2)*Bessel_I0(K2)*G_F(Kx, Ky)*K3;
   rates[0] = K4*rates[2];
   rates[1] = 0.;
   */
   rates[0] = K1;
   rates[1] = K1;
   rates[2] = K1;

   BK = U_Exp(-K2)*Bessel_I0(K2)*G_F(Kx, Ky)*K3/2.343;

   return rates;
}

// ============================== Gas_Relaxation ===============================

vectorU xIBS::Gas_Relaxation(xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   doubleU ri(0,m_);
   doubleU MxC(0, m_^2);
   doubleU MpxC(0, m_^2);
   doubleU emitr(0,m_);
   doubleU coup;

   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;

   if(DispersionCoeff)
   {  DispersionCoeff = false;
      for (int i=0; i < ring.Number(); i++)
      {
         MxC +=  ring[i].Lattice.Dx  * ring[i].Lattice.Dx *
                 ring[i].Length      / ring[i].Lattice.betax;
         M =    (ring[i].Lattice.Dpx * ring[i].Lattice.betax) +
                (ring[i].Lattice.Dx  * ring[i].Lattice.alfax) ;
         MpxC+= M*M * ring[i].Length / ring[i].Lattice.betax;
      }

      Mx  = MxC  / ring.Arc;
      Mpx = MpxC / ring.Arc;
      M = Mx + Mpx;
   }

    emitr = ((2.*(beam.Emit[0]*beam.Emit[0]+beam.Emit[1]*beam.Emit[1]))^0.5);

   //
  //  Gas-relaxation formula for bunched beam (Alexei: August 2004):
  // It was rewritten in terms of emittance_x:
  //  rates[2] = ri*ri* beam.Emit[3]* U_c*20./((8.*((U_pi/2.0)^0.5))*ring.Energy.Gamma3*ring.Energy.Beta3*
  //             (beam.Emit[0]^1.5)*(Lattice.betax^0.5)*s_s*beam.Emit[2]);

  /* Old formula:
   rates[2] = 2.*ri*ri* beam.Emit[3]* U_c*20./((2*(U_pi)^0.5)*ring.Energy.Gamma3*ring.Energy.Beta3*
             (emitr^1.5)*(Lattice.betax^0.5)*s_s*beam.Emit[2]);
  */


   // High-energy approximation of Bjorken-Mtingwa (Alexei: August 2004).
   // This rates are already for emittances so that they are larger by factor
   // of 2 than those for sigmas.
   //
   //  rates[2] = ri*ri* beam.Emit[3]* U_c*20./(8.*ring.Energy.Gamma3*ring.Energy.Beta3*
   //             (beam.Emit[0]^1.5)*(Lattice.betax^0.5)*s_s*beam.Emit[2]);
   //


// Since no "Bolide" form is made for Bjorken-Mtingwa yet - we switch from
// B-M to G-R directly in the code - using "Gas Relaxation" in the interface.

  if (Model == 4)
// Use G-R:
  rates[2] = ri*ri* beam.Emit[3]* U_c*20./((8.*((U_pi/2.0)^0.5))* ring.Energy.Gamma3 *ring.Energy.Beta3*
            (beam.Emit[0]^1.5)*(Lattice.betax^0.5)* s_s * beam.Emit[2]);
//  else
// Use B-M:
//  rates[2] = ri*ri* beam.Emit[3]* U_c*20./(8.*ring.Energy.Gamma3*ring.Energy.Beta3*
//             (beam.Emit[0]^1.5)*(Lattice.betax^0.5)*s_s*beam.Emit[2]);

// coup is coupling coefficient. For full coupling coup=1/2

  coup=0.5;

  rates[0] = (beam.Emit[2]/beam.Emit[0])* M *rates[2]*(1.-coup);
 //  rates[1] = rates[0];
  rates[1] = (beam.Emit[2]/beam.Emit[1])*M*rates[2]*coup;


   return rates;
}

// ============================== Bjorken - Mtingwa high energy appr. ======================
vectorU xIBS::Bjorken_Mtingwa_HE(xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   doubleU ri(0,m_);
   doubleU MxC(0, m_^2);
   doubleU MpxC(0, m_^2);
   doubleU emitr(0,m_);
   doubleU coup;

   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;

   for (int i=0; i < ring.Number(); i++)
      {
         MxC +=  ring[i].Lattice.Dx  * ring[i].Lattice.Dx *
                 ring[i].Length      / ring[i].Lattice.betax;
         M =    (ring[i].Lattice.Dpx * ring[i].Lattice.betax) +
                (ring[i].Lattice.Dx  * ring[i].Lattice.alfax) ;
         MpxC+= M*M * ring[i].Length / ring[i].Lattice.betax;
      }

      Mx  = MxC  / ring.Arc;
      Mpx = MpxC / ring.Arc;
      M = Mx + Mpx;

    emitr = ((2.*(beam.Emit[0]*beam.Emit[0]+beam.Emit[1]*beam.Emit[1]))^0.5);

// Use B-M high energy approximation:
  rates[2] = ri*ri* beam.Emit[3]* U_c*20./(8.*ring.Energy.Gamma3*ring.Energy.Beta3*
             (beam.Emit[0]^1.5)*(Lattice.betax^0.5)*s_s*beam.Emit[2]);

// coup is coupling coefficient. For full coupling coup=1/2

  coup=0.5;

  rates[0] = (beam.Emit[2]/beam.Emit[0])* M *rates[2]*(1.-coup);
  rates[1] = (beam.Emit[2]/beam.Emit[1])*M*rates[2]*coup;

 return rates;
}

// ============================== Bjorken - Mtingwa ===============================

vectorU xIBS::Bjorken_Mtingwa2(xLattice& lat, xBeam& beam, xRing& ring)
{
   doubleU ri(m_);
   doubleU A(s_^-1);

   doubleU phi_Bx;
   doubleU phi_By;
   doubleU Hx(m_);
   doubleU Hy(m_);

   double uplimt;

   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;

   A =  U_c * ri * ri * beam.Emit[3]*B_M_log /
   (8. * U_pi * beam.Emit[0]*beam.Emit[1] * s_s * (beam.Emit[2]^0.5) *
   ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);

   phi_Bx = lat.Dpx + (lat.alfax*lat.Dx/lat.betax);
   phi_By = lat.Dpy + (lat.alfay*lat.Dy/lat.betay);
   Hx = lat.betax*( (lat.Dpx + ((lat.alfax*lat.Dx)/lat.betax) )^2) +(lat.Dx*lat.Dx/lat.betax);
   Hy = lat.betay*( (lat.Dpy + ((lat.alfay*lat.Dy)/lat.betay) )^2) +(lat.Dy*lat.Dy/lat.betay);

   double Lh[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
   double Lv[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
   double Lp[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
   double L [3][3] = {{1,0,0},{0,1,0},{0,0,1}};

   Lh[0][0] = lat.betax(m_)/beam.Emit[0](m_);
   Lh[0][1] = (-ring.Energy.Gamma*lat.betax*phi_Bx/beam.Emit[0])(U1_);
   Lh[1][0] = Lh[0][1];
   Lh[1][1] = (ring.Energy.Gamma2*Hx/beam.Emit[0])(U1_);
   Lh[2][2] = 0;

   Lv[0][0] = 0;
   Lv[1][1] = (ring.Energy.Gamma2*Hy/beam.Emit[1])(U1_);
   Lv[1][2] = (-ring.Energy.Gamma*lat.betay*phi_By/beam.Emit[1])(U1_);
   Lv[2][1] = Lv[1][2];
   Lv[2][2] = (lat.betay/beam.Emit[1])(U1_);

   Lp[0][0] = 0;
   Lp[1][1] = (ring.Energy.Gamma2/beam.Emit[2])(U1_);
   Lp[1][2] = 0;

   uplimt = Lh[0][0];
   for (int i = 0; i < 3; i++)
   for (int j = 0; j < 3; j++)
   {
      L[i][j] = Lh[i][j] + Lv[i][j] + Lp[i][j];
      if(uplimt < fabs(L[i][j]))
         uplimt = fabs(L[i][j]);
   }
   uplimt *= lam_limit;

   double lambda, det, Tr;
   double LL[3][3], LB[3][3];
   double K[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
   double I[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

   for (int k=0; k < step_lam; k++)
   {
      lambda = uplimt * k / step_lam;

      for (int i=0; i < 3; i++)
         for (int j=0; j < 3; j++)
            LL[i][j] = I[i][j]*lambda + L[i][j]; 

//det = LL.Det();
det = (LL[0][0]*LL[1][1]*LL[2][2]) + (LL[0][1]*LL[1][2]*LL[2][0]) + (LL[0][2]*LL[1][0]*LL[2][1]) -
      (LL[0][2]*LL[1][1]*LL[2][0]) - (LL[0][0]*LL[1][2]*LL[2][1]) - (LL[0][1]*LL[1][0]*LL[2][2]);

//LB = LL.Back();
LB[0][0]=(LL[1][1]*LL[2][2])-(LL[1][2]*LL[2][1]);    LB[0][1]=(LL[2][0]*LL[1][2])-(LL[1][0]*LL[2][2]);     LB[0][2]=0;
LB[1][0]=LB[0][1]; LB[1][1]=(LL[0][0]*LL[2][2])-(LL[0][2]*LL[2][0]); LB[1][2]=(LL[2][0]*LL[0][1])-(LL[0][0]*LL[2][1]);
LB[2][0]=0;        LB[2][1]= LB[1][2];                               LB[2][2]=(LL[0][0]*LL[1][1])-(LL[0][1]*LL[1][0]);
      for (int i=0; i < 3; i++)
         for ( int j=0; j < 3; j++)
            LB[i][j] /= det;

//Tr = LB.Trace();
Tr = LB[0][0] + LB[1][1] + LB[2][2];

      for (int i=0; i < 3; i++)
         for (int j=0; j < 3; j++)
            K[i][j] += pow(lambda/det,0.5)*(Tr*((i==j)?1.:0.) - LB[i][j]*3.);
   }
   for (int i=0; i < 3; i++)
      for (int j=0; j < 3; j++)
         K[i][j] *= uplimt / step_lam;

   vectorU rates(0, s_^-1);

   for (int i=0; i < 3; i++)
   for (int j=0; j < 3; j++)
   {
      rates[0] += A * (K[i][j] * Lh[i][j]);
      rates[1] += A * (K[i][j] * Lv[i][j]);
      rates[2] += A * (K[i][j] * Lp[i][j]);
   }

   return rates;
}

vectorU xIBS::Bjorken_Mtingwa(xLattice& lat, xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   doubleU ri(0,m_);
   doubleU A(0,s_^-1);

   double phi_Bx;
   double phi_By;
   doubleU Hx(0,m_);
   doubleU Hy(0,m_);

   doubleU Ixx(s_^-1);
   doubleU Ixz(s_^-1);
   doubleU Iyy(s_^-1);
   doubleU Iyz(s_^-1);
   doubleU Izz(s_^-1);

   matrixU L (3,3);

   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;        // ion classical radius

   A =  U_c * ri * ri * beam.Emit[3]*B_M_log /
   (8. * U_pi * beam.Emit[0]*beam.Emit[1] * s_s * (beam.Emit[2]^0.5) *
   ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);

   phi_Bx = (lat.Dpx + (lat.alfax*lat.Dx/lat.betax))(U1_);
   phi_By = (lat.Dpy + (lat.alfay*lat.Dy/lat.betay))(U1_);
   Hx = lat.betax*( (lat.Dpx + ((lat.alfax*lat.Dx)/lat.betax) )^2) +(lat.Dx*lat.Dx/lat.betax);
   Hy = lat.betay*( (lat.Dpy + ((lat.alfay*lat.Dy)/lat.betay) )^2) +(lat.Dy*lat.Dy/lat.betay);

   L[0][0] = lat.betax(m_)/beam.Emit[0](m_);
   L[0][1] = 0;
   L[0][2] = ((-1.)*lat.betax(m_)*phi_Bx*ring.Energy.Gamma/beam.Emit[0](m_))();

   L[1][0] = 0;
   L[1][1] = lat.betay(m_)/beam.Emit[1](m_);
   L[1][2] = ((-1.)*lat.betay(m_)*phi_By*ring.Energy.Gamma/beam.Emit[1](m_))();

   L[2][0] = ((-1.)*ring.Energy.Gamma*lat.betax(m_)*phi_Bx/beam.Emit[0](m_))();
   L[2][1] = ((-1.)*ring.Energy.Gamma*lat.betay(m_)*phi_By/beam.Emit[1](m_))();
   L[2][2] = (ring.Energy.Gamma2*((Hx(m_)/beam.Emit[0](m_)) + (Hy(m_)/beam.Emit[1](m_)) + (1./beam.Emit[2])))();

   uplimit = real(L[0][0]);
   if(uplimit < real(L[1][1]))
      uplimit = real(L[1][1]);
   if(uplimit < real(L[2][2]))
      uplimit = real(L[2][2]);
   if(uplimit < fabs(real(L[2][0])))
      uplimit = fabs(real(L[2][0]));
   uplimit *= lam_limit;

   Ixx = A * B_M_integral_I(L, 1);
   Iyy = A * B_M_integral_I(L, 2);
   Izz = A * B_M_integral_I(L, 3);

   // 05.06
   uplimit = real(L[1][1]);
   if(uplimit < fabs(real(L[2][0])))
      uplimit = fabs(real(L[2][0]));
   uplimit *= lam_limit;

   Ixz = A * B_M_integral_I(L, 4);
   Iyz = A * B_M_integral_I(L, 5);

   rates[0] = ( (Hx*ring.Energy.Gamma2*Izz) -
               (2.*lat.betax*phi_Bx*ring.Energy.Gamma*Ixz) +
               (lat.betax*Ixx) )/beam.Emit[0];

   rates[1] = ( (Hy*ring.Energy.Gamma2*Izz) -
               (2.*lat.betay*phi_By*ring.Energy.Gamma*Iyz) +
               (lat.betax*Iyy) )/beam.Emit[1];

   rates[2] = ring.Energy.Gamma2*Izz/beam.Emit[2];

   return rates;
}

// ============================== Bjorken - Mtingwa integral  XX ===================

double xIBS::B_M_integral_I(matrixU L, int nomer)
{
   double result = 0;
   double L1[3][3];
   double LL[3][3];
   double E[3][3];
   double besk = uplimit;
   int step_l = (int)step_lam;
   double lambda = 0;
   double Min00 = 0;
   double Min11 = 0;
   double Min22 = 0;
   double det = 0;
   double TLL = 0;
   int i,j,k;

   for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      L1[i][j] = real(L[i][j]);              // move matrix to double array

   for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
     if (i==j) E[i][j] = 1.;
     else E[i][j] = 0.;                      // making unit matrix

   for (k=0; k < step_l; k++)
    {
     lambda = besk * k / step_l;             // integration area on step

     for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
       LL[i][j] = E[i][j]*lambda + L1[i][j];     // new matrix Lambda building

      // determinator of original matrix
      det = (LL[0][0]*LL[1][1]*LL[2][2]) + (LL[0][1]*LL[1][2]*LL[2][0]) + (LL[0][2]*LL[1][0]*LL[2][1]) -
            (LL[0][2]*LL[1][1]*LL[2][0]) - (LL[0][0]*LL[1][2]*LL[2][1]) - (LL[0][1]*LL[1][0]*LL[2][2]);

      // Determinators of minors for diagonal elements:
      Min00 = (LL[1][1]*LL[2][2])-(LL[1][2]*LL[2][1]);
      Min11 = (LL[0][0]*LL[2][2])-(LL[0][2]*LL[2][0]);
      Min22 = (LL[0][0]*LL[1][1])-(LL[0][1]*LL[1][0]);

    if (nomer == 1)
     TLL = Min11 + Min22 - (Min00*2.);   // for xx nomer = 1
    else if (nomer == 2)
     TLL = Min00 + Min22 - (Min11*2.);   // for yy nomer = 2
    else if (nomer == 3)
     TLL = Min11 + Min00 - (Min22*2.);   // for zz nomer = 3
    else if (nomer == 4)
     TLL = 3.*LL[1][1]*LL[0][2];         // for xz nomer = 4
    else if (nomer == 5)
     TLL = 3.*LL[0][0]*LL[1][2];         // for yz nomer = 5

     result += (pow((lambda/det),0.5)/det)* TLL;
    }
   result = result*besk/step_l;

 return result;
}

//*****************************************************************************

double xIBS::B_M_integral(int i, int j, matrixU L)
{
 complex<double> result = 0;
 matrixU Lambda(3,3);
 matrixU Lambda_(3,3);
 matrixU E(3,3);
 matrixU T(3,3);
 double besk = uplimit;
 int step_l = (int)step_lam;
 double Kr = 0.;
 double lambda = 0;
 complex<double> det = 0;
 complex<double> tr = 0;

 if(i==j) Kr = 1.;

 E.Diag(1., 0.);              // unit matrix

   for (int k=0; k < step_l; k++)
    {
     lambda = besk * k / step_l;            // area for integral on every step

     Lambda = E*lambda + L;                 // new matrix Lambda
     Lambda_ = Lambda.Back();               // inverse matrix

     det = Lambda.Det();                    // matrix determinator

     result += pow((lambda/Lambda.Det()),0.5)*((Lambda_.Trace()*Kr) - (3.*Lambda_[i][j]));
    }
   result = result*besk/(complex<double>)step_l;

 return real(result);
}

// ============================== Molecular Dynamics ===========================

vectorU xIBS::MD_Rates(xTime& time, xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);

   vectorU emit0, emit1;
   doubleU keep_ds(time.ds);
   time.DistStep(iDynamics.ds);
   if (iDynamics.Algoritm == 0)
      beam.Distribution(beam.Emit,Lattice);
   emit0 = iBeam.Emit_Def(Lattice);
   beam.SpaceChargeMD(time);
   emit1 = iBeam.Emit_Def(Lattice);
   rates = (emit1 - emit0) / (emit0 * time.dt);
   time.DistStep(keep_ds);

   return rates;
}

// ============================== IBS Kick =====================================

void xIBS::Kick (xTime&time, xBeam&beam, xRing&ring)
{
   if (beam.UseForIBS)
      beam.Emit = beam.Emit_Def(Lattice);
   else
      beam.Emit = beam.Emit_RMS(Lattice);

   int i, j;
   double Nc[3], Nt[3];
   vectorU Ec, Et, Rc, Rt;
   vectorU rates;
   Tine<bool>InCore;
   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;
	bpRates Growth;																									
	bpBeam bplusBeam;
	bpRing bplusRing;
   switch (IBSkickModel)
   {case 0:
      if (Slices && Model == 3)
      {
	     bplusBeam.emit.emittancex = beam.Emit[0](m_);
	     bplusBeam.emit.emittancey = beam.Emit[1](m_);
	     bplusBeam.emit.momentum2 = beam.Emit[2]();
	     bplusBeam.emit.sigmab = beam.s_s(m_);
	     bplusBeam.q0 = ring.Energy.Z();
	     bplusBeam.m0 = ring.Energy.A();
	     bplusBeam.energy = ring.Energy.Kinetic(M_6*eV_);
	     bplusBeam.nion = beam.Emit[3]();
	     bplusBeam.bunched = beam.benum;
        bplusBeam.setsize(beam.Number());
        for (int j = 0; j < 6; j++)
          for (int i = 0; i < bplusBeam.row; i++)
            bplusBeam[j][i] = beam[i][j];
        bplusBeam.Init();

	     bplusRing.Circ = ring.Circ(m_);
        bplusRing.GammaTr = ring.GammaTr();
        bplusRing.Imagenery = ring.Imagenary();                          
        bplusRing.Q_h = ring.TunesH();
        bplusRing.Q_v = ring.TunesV();
        bplusRing.Hrf = ring.h();
        bplusRing.Vrf = ring.V();
        bplusRing.set(ring.Number());
	     for (int l = 0; l < ring.Number(); l++)
	     {
		     bplusRing[l].s = ring[l].Lattice.dist(m_);
		     bplusRing[l].ds = ring[l].Length(m_);
		     bplusRing[l].betax = ring[l].Lattice.betax(m_);
		     bplusRing[l].alphax = ring[l].Lattice.alfax();
		     bplusRing[l].betay = ring[l].Lattice.betay(m_);
		     bplusRing[l].alphay = ring[l].Lattice.alfay();
		     bplusRing[l].Dx = ring[l].Lattice.Dx(m_);
		     bplusRing[l].Dpx = ring[l].Lattice.Dpx();
		     bplusRing[l].Dy = ring[l].Lattice.Dy(m_);
		     bplusRing[l].Dpy = ring[l].Lattice.Dpy();
	     }
        bplusRing.Init(bplusBeam);
		  bplusRing.lattice.s = ring.LATTICE.dist(m_);
		  bplusRing.lattice.ds = ring.Circ(m_);
		  bplusRing.lattice.betax = ring.LATTICE.betax(m_);
		  bplusRing.lattice.alphax = ring.LATTICE.alfax();
		  bplusRing.lattice.betay = ring.LATTICE.betay(m_);
		  bplusRing.lattice.alphay = ring.LATTICE.alfay();
		  bplusRing.lattice.Dx = ring.LATTICE.Dx(m_);
		  bplusRing.lattice.Dpx = ring.LATTICE.Dpx();
		  bplusRing.lattice.Dy = ring.LATTICE.Dy(m_);
		  bplusRing.lattice.Dpy = ring.LATTICE.Dpy();
        bpIntegrateAcc(bpStepSize, bpStepNumber);
        bpSlices (bplusRing, bplusBeam, bpSliceNumber, time.dt(s_));
        for (int j = 0; j < 6; j++)
          for (int i = 0; i < bplusBeam.row; i++)
             beam[i][j] = bplusBeam[j][i];
      }else
      if (Slices && beam.benum != COASTING)
      {  if (beam.benum == BUNCHED || !ring.Bucket.IBSnorma)
         {  for (int i = 0; i < beam.LProfile.Number; i++)
            {  beam.Emit[2]  = beam.Momentum[i];
               if (beam.Emit[2]())
               {  rates = Rates(time, beam, ring);
                  beam.RatesDp[i] = rates[2]()*beam.LProfile[i];
                  for (j = 0; j < beam.Number(); j++)
                  if (i == beam.index(j, beam.LProfile.Number))
                     IndividualKick(time, beam, j, rates*beam.LProfile[i], beam.Emit, 1.);
               }else
                  beam.RatesDp[i] = 0;
            }
         }else
         {  for (int i = 0; i < ring.Bucket.barrier.Number; i++)
            {  beam.Emit[2]  = ring.Bucket.barrier[i].P;
               if (beam.Emit[2]())
               {  rates = Rates(time, beam, ring);
                  for (j = 0; j < beam.Number(); j++)
                  if (i == ring.Bucket.Pos(beam, j))
                  {  //beam.RatesDp[i] = rates[2]()*ring.Bucket.barrier[i].D;
                     IndividualKick(time, beam, j, rates*ring.Bucket.barrier[i].D, beam.Emit, 1.);
                  }
               }
            }
         }
      }else
      {
         rates = Rates(time, beam, ring);
         for (j = 0; j < beam.Number(); j++)
            IndividualKick(time, beam, j, rates, beam.Emit, 1.);
      }
      break;

   case 1: case 2: case 3:
      Et = beam.Emit;
      Ec = beam.Emit_FWHM(Lattice);
      beam.Courant_Snyder(Lattice);

      if (IBSkickModel == 1)
      for(i = 0; i < 3; i++)
      if ((i != 2) || (beam.benum == BUNCHED))
      {
         beam.Powel2(beam.Fitting[i], beam.Hyst3[i], beam.Sigmas, beam.Division);
         if (fabs(beam.Fitting[i][1]) > fabs(beam.Fitting[i][2]))
         {
            Ec[i] = beam.Emit[i] * beam.Fitting[i][2] * beam.Fitting[i][2];
            Et[i] = beam.Emit[i] * beam.Fitting[i][1] * beam.Fitting[i][1];
            Nc[i] = 1. / (1. + fabs
                  ( beam.Fitting[i][1] / beam.Fitting[i][2] *
                    beam.Fitting[i][3] / beam.Fitting[i][4] ) );
         }else
         {
            Ec[i] = beam.Emit[i] * beam.Fitting[i][1] * beam.Fitting[i][1];
            Et[i] = beam.Emit[i] * beam.Fitting[i][2] * beam.Fitting[i][2];
            Nc[i] = 1. / (1. + fabs
                  ( beam.Fitting[i][2] / beam.Fitting[i][1] *
                    beam.Fitting[i][4] / beam.Fitting[i][3] ) );
         }
      }
      if(CoreFWHM) Ec = beam.Emit_FWHM(Lattice);

      CoreNumber = 0;
      InCore.SetNumber(beam.Number());
      for (i = 0; i < beam.Number(); i++)
      {
         int InsideCore = 0;
         for (j = 0; j < 3; j++)
            if (Ec[j] * 2 * CoreDefinition > beam.Inv(i,j))
               InsideCore++;

         if (InsideCore >= InCoreDefine)
         {  CoreNumber ++;
            InCore[i] = true;
         }else
            InCore[i] = false;
      }
/////////////////////////////////////////////////////////////////////27/09/2007
      for (i = 0; i < 3; i++)
/////////if ((IBSkickModel != 1) || CoreNumCor || ((i==2)&&(beam.benum!=BUNCHED)) )
         Nc[i] = CoreNumber / beam.Number();

      beam.Emit = Ec;
      Rc = Rates(time, beam, ring);
      if (IBSkickModel == 3)
         Et = beam.Emit_Tail(Lattice, InCore, (int)CoreNumber);
      beam.Emit = Et;
      Rt = Rates(time, beam, ring);

      for (i = 0; i < 3; i++)
      {  Nt[i] = 1 - Nc[i];
         Nc[i] *= Nc[i];
         Nt[i] *= Nt[i];
      }
      for (i = 0; i < 3; i++)
         Rc[i] = (Rc[i] * Nc[i]) + (Rt[i] * Nt[i] * Et[i] / Ec[i]);

      for (int i = 0; i < beam.Number(); i++)
         if (InCore[i])
            //IndividualKick(time, beam, i, Rc*ring.Bucket.Density(beam,ring,j), Ec, 1.);
            IndividualKick(time, beam, i, Rc, Ec, 1.);
         else
            //IndividualKick(time, beam, i, Rt*ring.Bucket.Density(beam,ring,j), Et, 1.);
            IndividualKick(time, beam, i, Rt, Et, 1.);

      beam.Emit = beam.Emit_Def(Lattice);
      break;

   case 4:
      rates = GP_BiGauss(Lattice, beam, ring);
      for (j = 0; j < beam.Number(); j++)
         IndividualKick(time, beam, j, rates, beam.Emit, 1.);
      break;

    case 5:
      ZenkevichKick(time, beam, ring);
      break;

    case 6:
      beam.RMS(beam.Emit,Lattice);
      Model = 2;
      rates = Rates(time, beam, ring);
      Model = 3;
      for (j = 0; j < beam.Number(); j++)
      {  Burov(Lattice, beam(j), ring);
         IndividualKick(time, beam, j, rates, beam.Emit, BK());
      }
      break;

    case 7:
      LocalKick(time, beam, ring);
      break;
   }
}

void xIBS::IndividualKick(xTime&time, xBeam&beam, int j,
                          vectorU rates, vectorU&emit, double k)
{
   double kb = 1;
   if(beam.benum == BUNCHED) kb = 2;

   if (rates[0]() < 0)
      beam(j,1, beam(j,1) * U_Exp(rates[0] * time.dt * k));
   else
      beam.Add(j, 1, ( ((emit[0] / Lattice.betax) *
      rates[0] * time.dt * k * 2.0)^0.5)*xDistributor::Gaussian());

   if (rates[1]() < 0)
      beam(j,3, beam(j,3) * U_Exp(rates[1] * time.dt * k));
   else
      beam.Add(j, 3, ( ((emit[1] / Lattice.betay) *
      rates[1] * time.dt * k * 2.0)^0.5)*xDistributor::Gaussian());

   if (rates[2]() < 0)
      beam(j,5, beam(j,5) * U_Exp(rates[2] * time.dt * k));
   else
      beam.Add(j, 5, ( ( emit[2] *
      rates[2] * time.dt * k * kb)^0.5)*xDistributor::Gaussian());
}

//Calculation of diffusion and friction from distribution
/*
void xIBS::SidorinKick(xTime&time, xBeam&beam, xRing&ring)
{
   matrixU Matrix;

   for (int n = 0; n < ring.Number(); n++)
   {
      beam.LRF2PRF(ring.Energy);

// kick here

      beam.PRF2LRF(ring.Energy);


      Matrix = ring.GetMatrix(time, n);
      for (int i = 0; i < beam.Number(); i++)
         beam(i, Matrix);

   }
}
*/

//Zenkevich kinetic model 05.06

void xIBS::ZenkevichKick(xTime&time, xBeam&beam, xRing&ring)
{
matrixU Matrix;

double ct;

double longdrift = 0.0;

double kb = 1;
   if(beam.benum == BUNCHED) kb = 2;
   //kb = 1.0;


double ksi1, ksi2, ksi3;
doubleU Cz1(s_^-0.5), Cz2(s_^-0.5), Cz3(s_^-0.5), DetD(s_^-3), Cz4(s_^-0.5);
      /*
      for (int n = 0; n < ring.Number(); n++)
           {
           Matrix = ring.GetMatrix(time, n);
           longdrift += real(Matrix[4][5]);
           }
//longitudinal matching
      for (int i = 0; i < beam.Number(); i++)
         beam[i][4] -= 0.0*longdrift*beam[i][5];
      */

   for (int n = 0; n < ring.Number(); n++)
   {
//Diffusion and friction calculation
    B_M_diffusion(ring[n].Lattice, beam, ring);
    B_M_friction(ring[n].Lattice, beam, ring);

//Calculation of the kick coefficients (it should be independent procedure)
      DetD = D_xx * D_yy * D_zz - D_zx*D_zx*D_zy - D_zx*D_zy*D_zy;
      Cz1 = D_zx/(D_xx^0.5);
      Cz4 = ((U_Abs((DetD/(2.0*D_xx*D_yy))))^0.5);
      Cz2 = (D_zy/((D_yy/2.0)^0.5)) - Cz4;
      Cz3 = (D_zy/((D_yy*2.0)^0.5)) + Cz4;

      Cz4 = ((D_xx)^0.5);
// kick
        for (int j = 0; j < beam.Number(); j++)
      {
         ksi1 = xDistributor::Gaussian();
         ksi2 = xDistributor::Gaussian();
         ksi3 = xDistributor::Gaussian();


         beam[j][1] += (-1.0* F_x *(beam[j][1] - (ring[n].Lattice.Dpx*beam[j][5]))*
                       time.dt*ring[n].Length/ring.Arc)(U1_);
         beam[j][3] += (-1.0* F_y *(beam[j][3] - (ring[n].Lattice.Dpy*beam[j][5]))*
                       time.dt*ring[n].Length/ring.Arc)(U1_);
         beam[j][5] += (-1.0* F_z *kb*beam[j][5]*time.dt*ring[n].Length/ring.Arc)(U1_);

         //beam[j][1] += (((D_xx*time.dt*ring[n].Length/ring.Arc)^0.5))(U1_)*ksi1;

         beam[j][1] += (((time.dt*ring[n].Length/ring.Arc)^0.5)*Cz4)(U1_)*ksi1;

         beam[j][3] += (((0.5*D_yy*time.dt*ring[n].Length/ring.Arc)^0.5))(U1_)*ksi2;
         beam[j][3] += (((0.5*D_yy*time.dt*ring[n].Length/ring.Arc)^0.5))(U1_)*ksi3;


         beam[j][5] += (((time.dt*ring[n].Length*kb/ring.Arc)^0.5)*Cz1)(U1_)*ksi1;
         beam[j][5] += (((time.dt*ring[n].Length*kb/ring.Arc)^0.5)*Cz2)(U1_)*ksi2;
         beam[j][5] += (((time.dt*ring[n].Length*kb/ring.Arc)^0.5)*Cz3)(U1_)*ksi3;



      }

      Matrix = ring.GetMatrix(time, n);
      for (int i = 0; i < beam.Number(); i++)
         beam(i, Matrix);
      longdrift += real(Matrix[4][5]);
   }
   for (int i = 0; i < beam.Number(); i++)
         beam[i][4] -= 1.0*longdrift*beam[i][5];

   if (AverageTrans)
            {
             for (int j = 0; j < beam.Number(); j++)
               {
               ct = beam[j][0];
               beam[j][0] = -1.0*beam[j][2];
               beam[j][2] = ct;
               }
            }
}

//Local longitudinal diffusion 05.06

void xIBS::LocalKick(xTime&time, xBeam&beam, xRing&ring)
{
   matrixU Matrix;

   double ct;
   doubleU lengthU(1,m_);
   doubleU velocityU(1,m_/s_);
   doubleU tempF(1,m_/s_);
   int Loc_n;

   beam.Minst.Number = IBSloc;
   if(beam.Minst.Number < 1) beam.Minst.Number = 1;
   if(beam.Minst.Number > beam.Bunch.Number) beam.Minst.Number = beam.Bunch.Number;

   beam.Minst.SetNumber(beam.Minst.Number);
   beam.Dist2.SetNumber(beam.Bunch.Number);

   vectorU ion_local(m_,m_/s_);
   doubleU dV(m_/s_);
   doubleU v[3];
   doubleU f[3];                        // forces components for 3D
   doubleU D[3][3];                     // difusion tenzor

   doubleU DetD((eV_^6)*(s_^3)/(m_^6)); //determinant of the diffusion tenzot
   doubleU SpD((eV_^2)*s_/(m_^2));      //Spur of the diffusion tensor
   doubleU SumD((eV_^4)*(s_^2)/(m_^4)); //Sum of minors

   double ksi1, ksi2, ksi3;
   doubleU Cz1((eV_^2)*s_/(m_^2)), Cz2((eV_^2)*s_/(m_^2)), Cz3((eV_^2)*s_/(m_^2)), Cz4((eV_^2)*s_/(m_^2));

   for (int i = 0; i < 3; i++)
   {
      v[i]._(m_/s_);
      f[i]._(eV_/m_);
      for(int j=0; j < 3; j++)
      D[i][j]._((eV_^2)*s_/(m_^2));
   }

   double longdrift = 0.0;

   for (int n = 0; n < ring.Number(); n++)
   {
      beam.Local = beam.Bunch;
      beam.LRF_PRF(ring.Energy, beam.Local);

      beam.RMS(beam.Emit,ring.GetLattice(n));

      for (int j1 = 0; j1 < beam.Number(); j1++)
      {
         ion_local[0] = beam(j1)[0];
         ion_local[1] = beam(j1)[1] * U_c * ring.Energy.Beta * ring.Energy.Gamma;
         ion_local[2] = beam(j1)[2];
         ion_local[3] = beam(j1)[3] *U_c * ring.Energy.Beta * ring.Energy.Gamma;
         ion_local[4] = beam(j1)[4] *ring.Energy.Gamma;
         ion_local[5] = beam(j1)[5] *U_c * ring.Energy.Beta;
         beam.GetMean(ion_local, beam.Local);                 // get Mean in PRF

         //---------- 26.04.07 ------------
         //Mean density and maximum impact parameter

         v[0] = ion_local[1];
         v[1] = ion_local[3];
         v[2] = ion_local[5];

         doubleU n_loc(m_^-3);
         doubleU n_loc_an(m_^-3);
         doubleU deference(U1_);
         /*
         //2-D uniform + Gauss
         n_loc = beam.Minst.Number*beam.Emit[3] /(2.*sqrt(2.)*2.*M_PI*
                 sqrt(2.*M_PI*beam.Mean2[0]*beam.Mean2[2]*beam.Mean2[4])*
                 beam.Bunch.Number*lengthU*lengthU*lengthU);

         //ellipsoid
         n_loc = 2.0*beam.Minst.Number*beam.Emit[3] /(5.0*sqrt(5.0)*(4.0/3.0)*M_PI*
                 sqrt(beam.Mean2[0]*beam.Mean2[2]*beam.Mean2[4])*
                 beam.Bunch.Number*lengthU*lengthU*lengthU);
         
         //cube
         n_loc = beam.Minst.Number*beam.Emit[3] /(3.0*sqrt(3.0)*8.0*
                 sqrt(beam.Mean2[0]*beam.Mean2[2]*beam.Mean2[4])*
                 beam.Bunch.Number*lengthU*lengthU*lengthU);
         */

         //box
         Loc_n = 0;
         for (int i = 0; i < beam.Minst.Number; i++)
            if(((( beam.Local[beam.Minst[i]][4]-ion_local[4]()) * (beam.Local[beam.Minst[i]][4]-ion_local[4]()))<(1.0*beam.Mean2[4]))&&
               ((((beam.Local[beam.Minst[i]][0]-ion_local[0]()) * (beam.Local[beam.Minst[i]][0]-ion_local[0]())/(beam.Mean2[0]*Box_IBS*Box_IBS))+
                 ((beam.Local[beam.Minst[i]][2]-ion_local[2]()) * (beam.Local[beam.Minst[i]][2]-ion_local[2]())/(beam.Mean2[2]*Box_IBS*Box_IBS)))<1.0
               )) Loc_n++;

         n_loc = Loc_n*beam.Emit[3] /(2.0*M_PI* Box_IBS*Box_IBS*sqrt(1.0* beam.Mean2[0]* beam.Mean2[2]* beam.Mean2[4])*
                 beam.Bunch.Number*lengthU*lengthU*lengthU);

         /*
         //analytical density
         n_loc_an = beam.Emit[3]* U_Exp(((-1.0*ion_local[0]*ion_local[0])/(2.0*beam.Sigma_2[0])) +
                                        ((-1.0*ion_local[2]*ion_local[2])/(2.0*beam.Sigma_2[2])) +
                                        ((-1.0*ion_local[4]*ion_local[4])/(2.0*beam.Sigma_2[4]*ring.Energy.Gamma*ring.Energy.Gamma)) )/
                   (2.0*M_PI*sqrt(2.0*M_PI)*beam.Sigma[0]*beam.Sigma[2]*beam.Sigma[4]*ring.Energy.Gamma);

         deference = (n_loc - n_loc_an)/n_loc_an;
         n_loc = n_loc_an;
         */

         //Diffusion  and friction calculation
         if (n_loc() > 0.0)
         {
            doubleU delta = 8 * U_pi * n_loc * ring.Energy.Z * ring.Energy.Z * ring.Energy.Z *ring.Energy.Z *(U_e^4)/(ring.Energy.A *U_amu); //Friction constant
            doubleU theta(m_/s_);  // total ion velocity
            doubleU delta_v[3], delta_d[3][3], v_local[3];

            for (int i = 0; i < 3; i++)
            {
               delta_v[i]._(0, (s_^2)/(m_^2));
               v_local[i]._(m_/s_);
               for (int j = 0; j < 3; j++)
                  delta_d[i][j]._(0, s_/m_);
            }
            doubleU U_rel((m_^2)/(s_^2));  // square of relative velocity

            theta   = ((v[0]^2) + (v[1]^2) +(v[2]^2))^0.5;
            for (int j = 0; j < 3; j++)
            {
               f[j] = 0;             
               delta_v[j] =0;
               for(int jj=0; jj < 3; jj++)
               {
                  D[j][jj]= 0;
                  delta_d[j][jj] = 0;
               }
            }

            for (int i = 0; i < beam.Minst.Number; i++)
            {
               U_rel = 0;
               for (int j = 0; j < 3; j++)
               {
                  v_local[j] = beam.Local[beam.Minst[i]][2*j+1];
                  U_rel += (v[j] - v_local[j])^2;
               }
               if (U_rel() != 0.0)
               {
                  for (int a = 0; a < 3; a++)
                  {
                     if (theta() > 0.)delta_v[a] += (v[a] -  v_local[a]) * B_M_log / (U_rel^1.5);

                     for (int b = 0; b < 3; b++)
                     if (a == b)
                      delta_d[a][b]+= ((U_rel - ((v[a]-v_local[a])*(v[b]-v_local[b]))) * B_M_log)
                                    /(U_rel^1.5);
                     else
                      delta_d[a][b]+=        - ((v[a]-v_local[a])*(v[b]-v_local[b])) * B_M_log
                                    /(U_rel^1.5);
                  }
               }
            }
            for (int j = 0; j < 3; j++)
            {
               if (theta() > 0.) f[j] = -delta_v[j] * delta / beam.Minst.Number;

               for(int k = 0; k < 3; k++)
                  D[j][k] = 0.5 * delta_d[j][k] * delta * ring.Energy.A *U_amu/ beam.Minst.Number;
            }

            DetD=     D[0][0] * D[1][1] * D[2][2] +
                  0.0*D[0][1] * D[1][2] * D[2][0] +
                  0.0*D[0][2] * D[1][0] * D[2][1] -
                      D[0][2] * D[1][1] * D[2][0] -
                      D[0][0] * D[1][2] * D[2][1] -
                  0.0*D[0][1] * D[1][0] * D[2][2] ;
            SpD = D[0][0] + D[1][1] + D[2][2];
            SumD= D[0][0]*D[1][1] - D[1][0]*D[0][1] +
                  D[0][0]*D[2][2] - D[0][2]*D[2][0] +
                  D[1][1]*D[2][2] - D[1][2]*D[2][1];

            // kick
            // friction
            /*
            if (AverageTrans)
            {  f[0] = (f[0] + f[1]) / 2.;
               f[1] = f[0];
            }
            */
            tempF = ion_local[1];
            tempF -= beam.Mean[1]*velocityU;
            tempF %= f[0]*time.dt*ring[n].Length/(ring.Arc*ring.Energy.Gamma*ring.Energy.A *U_amu);
            tempF += beam.Mean[1]*velocityU;
            ion_local[1] = tempF;

            tempF = ion_local[3];
            tempF -= beam.Mean[3]*velocityU;
            tempF %= f[1]*time.dt*ring[n].Length/(ring.Arc*ring.Energy.Gamma*ring.Energy.A *U_amu);
            tempF += beam.Mean[3]*velocityU;
            ion_local[3] = tempF;

            tempF = ion_local[5];
            tempF -= beam.Mean[5]*velocityU;
            tempF %= f[2]*time.dt*ring[n].Length/(ring.Arc*ring.Energy.Gamma*ring.Energy.A *U_amu);
            tempF += beam.Mean[5]*velocityU;
            ion_local[5] = tempF;

            //diffusion
            //Calculation of the kick coefficients
            Cz1 = D[0][2]*D[0][2]/D[0][0];
            if (TranCoup)
               Cz1 /= 2.;
            Cz4 = U_Abs(DetD/(2.0*D[0][0]*D[1][1]));
            Cz2 = ((D[1][2]/((D[1][1]/2.0)^0.5)) - (Cz4^0.5))*((D[1][2]/((D[1][1]/2.0)^0.5)) - (Cz4^0.5));
            Cz3 = ((D[1][2]/((D[1][1]*2.0)^0.5)) + (Cz4^0.5))*((D[1][2]/((D[1][1]*2.0)^0.5)) + (Cz4^0.5));
            Cz4 = D[0][0];

            // kick
            ksi1 = xDistributor::Gaussian();
            ksi2 = xDistributor::Gaussian();
            ksi3 = xDistributor::Gaussian();

            ion_local[1] += ((Cz4*time.dt*ring[n].Length/ring.Arc/ring.Energy.Gamma)^0.5)
                            *ksi1/(ring.Energy.A *U_amu);
            ion_local[3] += ((0.5*D[1][1]*time.dt*ring[n].Length/ring.Arc/ring.Energy.Gamma)^0.5)
                            *ksi1/(ring.Energy.A *U_amu);
            ion_local[3] += ((0.5*D[1][1]*time.dt*ring[n].Length/ring.Arc/ring.Energy.Gamma)^0.5)
                            *ksi2/(ring.Energy.A *U_amu);
            ion_local[5] += ((Cz1*time.dt*ring[n].Length/ring.Arc/ring.Energy.Gamma)^0.5)
                            *ksi1/(ring.Energy.A *U_amu);
            ion_local[5] += ((Cz2*time.dt*ring[n].Length/ring.Arc/ring.Energy.Gamma)^0.5)
                            *ksi2/(ring.Energy.A *U_amu);
            ion_local[5] += ((Cz3*time.dt*ring[n].Length/ring.Arc/ring.Energy.Gamma)^0.5)
                            *ksi3/(ring.Energy.A *U_amu);
            beam[j1][1] = (ion_local[1] / (U_c * ring.Energy.Beta * ring.Energy.Gamma))(U1_);
            beam[j1][3] = (ion_local[3]/(U_c * ring.Energy.Beta * ring.Energy.Gamma))(U1_);
            beam[j1][5] = (ion_local[5] /(U_c * ring.Energy.Beta))(U1_);

         } // n_loc > 0.0
      } // cyrcle over particles

      Matrix = ring.GetMatrix(time, n);
      //for(int i=4;i<6;i++)for(int j=0;j<6;j++)if(i!=j)Matrix[i][j]=0;
      for (int i = 0; i < beam.Number(); i++)
         beam(i, Matrix);
      longdrift += real(Matrix[4][5]);

   } // cycle over matrixes

   //bunch matching
   for(int i=0;i<beam.Number();i++)beam[i][4]-=longdrift*beam[i][5];

   if (AverageTrans)
   {  for (int j = 0; j < beam.Number(); j++)
      {  ct = beam[j][0];
         beam[j][0] = -beam[j][2];
         beam[j][2] = ct;
      }
   }
}


void xIBS::ZenkevichKick2(xTime&time, xBeam&beam, xRing&ring)
{
matrixU Matrix;

double kb = 1;
   if(beam.benum == BUNCHED) kb = 2;

double ksi1, ksi2, ksi3;
doubleU Cz1(s_^-0.5), Cz2(s_^-0.5), Cz3(s_^-0.5), DetD(s_^-3), Cz4(s_^-0.5);
doubleU F_x1(0,m_/s_),F_y1(0,m_/s_),F_z1(0,m_/s_);
doubleU D_xx2(0,m_/s_), D_yy2(0,m_/s_),D_zz2(0,m_/s_),D_zx2(0,m_/s_),D_zy2(0,m_/s_);

   for (int n = 0; n < ring.Number(); n++)
   {
//Diffusion and friction calculation





    B_M_diffusion(ring[n].Lattice, beam, ring);
    B_M_friction(ring[n].Lattice, beam, ring);
    F_x1 += F_x*ring[n].Length;
    F_y1 += F_y*ring[n].Length;
    F_z1 += F_z*ring[n].Length;
    D_xx2 += D_xx*ring[n].Length;
    D_yy2 += D_yy*ring[n].Length;
    D_zz2 += D_zz*ring[n].Length;
    D_zx2 += D_zx*ring[n].Length;
    D_zy2 += D_zy*ring[n].Length;
   }
   D_xx = D_xx2/ring.Arc;
   D_yy = D_yy2/ring.Arc;
   D_zz = D_zz2/ring.Arc;
   D_zx = D_zx2/ring.Arc;
   D_zy = D_zy2/ring.Arc;
   F_x = F_x1/ring.Arc;
   F_y = F_y1/ring.Arc;
   F_z = F_z1/ring.Arc;
//Calculation of the kick coefficients (it should be independent procedure)
      DetD = D_xx * D_yy * D_zz - D_zx*D_zx*D_zy - D_zx*D_zy*D_zy;
      Cz1 = D_zx/(D_xx^0.5);
      Cz4 = ((DetD/(2.0*D_xx*D_yy))^0.5);
      Cz2 = (D_zy/((D_yy/2.0)^0.5)) - Cz4;
      Cz3 = (D_zy/((D_yy*2.0)^0.5)) + Cz4;

      Cz4 = ((D_xx)^0.5);
// kick
        for (int j = 0; j < beam.Number(); j++)
      {
         ksi1 = xDistributor::Gaussian();
         ksi2 = xDistributor::Gaussian();
         ksi3 = xDistributor::Gaussian();

         beam[j][1] += (-1.0* F_x * beam[j][1]* time.dt)(U1_);
         beam[j][3] += (-1.0* F_y * beam[j][3] * time.dt)(U1_);
         beam[j][5] += (-1.0* F_z *beam[j][5] * time.dt)(U1_);

         //beam[j][1] += (((D_xx*time.dt*ring[n].Length/ring.Arc)^0.5))(U1_)*ksi1;

         beam[j][1] += (((time.dt)^0.5)*Cz4)(U1_)*ksi1;

         beam[j][3] += (((0.5*D_yy*time.dt)^0.5))(U1_)*ksi1;
         beam[j][3] += (((0.5*D_yy*time.dt)^0.5))(U1_)*ksi2;

         beam[j][5] += (((time.dt*kb)^0.5)*Cz1)(U1_)*ksi1;
         beam[j][5] += (((time.dt*kb)^0.5)*Cz2)(U1_)*ksi2;
         beam[j][5] += (((time.dt*kb)^0.5)*Cz3)(U1_)*ksi3;
      }


}

void xIBS::B_M_friction(xLattice& lat, xBeam& beam, xRing& ring)
{
   doubleU ri(0,m_);
   doubleU A(0,s_^-1);

   double phi_Bx;
   double phi_By;
   doubleU Hx(0,m_);
   doubleU Hy(0,m_);

   matrixU L (3,3);

   double L1[3][3];
   double LL[3][3];
   double E[3][3];

   double result1, result2, result3;
   double besk;
   int step_l;

   double lambda = 0;
   double Min00 = 0;
   double Min11 = 0;
   double Min22 = 0;
   double det = 0;
   double TLL = 0;
   int i,j,k;

   if (beam.benum == BUNCHED)
   {  nb = 1;
      beam.CalcBunch();
      s_s = beam.s_s;
   }else
   {  if (beam.benum == COASTING)
         nb = 2;
      else
         nb = BUCKETnb;   
      s_s = ring.Circ / (2 * (U_pi^0.5));
   }

   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;        // ion classical radius

        A =  U_c * ri * ri * beam.Emit[3]*B_M_log /
        (8. * U_pi * beam.Emit[0]*beam.Emit[1] * s_s * (beam.Emit[2]^0.5) *
        ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);

       phi_Bx = (lat.Dpx + (lat.alfax*lat.Dx(m_)/lat.betax(m_)))();
       phi_By = (lat.Dpy + (lat.alfay*lat.Dy(m_)/lat.betay(m_)))();
       Hx = lat.betax*( (lat.Dpx + ((lat.alfax*lat.Dx)/lat.betax) )^2) +(lat.Dx*lat.Dx/lat.betax);
       Hy = lat.betay*( (lat.Dpy + ((lat.alfay*lat.Dy)/lat.betay) )^2) +(lat.Dy*lat.Dy/lat.betay);

       L[0][0] = lat.betax(m_)/beam.Emit[0](m_);
       L[0][1] = 0;
       L[0][2] = ((-1.)*lat.betax(m_)*phi_Bx*ring.Energy.Gamma/beam.Emit[0](m_))();

       L[1][0] = 0;
       L[1][1] = lat.betay(m_)/beam.Emit[1](m_);
       L[1][2] = ((-1.)*lat.betay(m_)*phi_By*ring.Energy.Gamma/beam.Emit[1](m_))();

       L[2][0] = ((-1.)*ring.Energy.Gamma*lat.betax(m_)*phi_Bx/beam.Emit[0](m_))();
       L[2][1] = ((-1.)*ring.Energy.Gamma*lat.betay(m_)*phi_By/beam.Emit[1](m_))();
       L[2][2] = (ring.Energy.Gamma2*((Hx(m_)/beam.Emit[0](m_)) + (Hy(m_)/beam.Emit[1](m_)) + (1./beam.Emit[2])))();

       uplimit = real(L[0][0]);
       if (uplimit < real(L[1][1])) uplimit = real(L[1][1]);
       if (uplimit < real(L[2][2])) uplimit = real(L[2][2]);

       uplimit *= lam_limit;
       besk = uplimit;
       step_l = (int)step_lam;



   for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      L1[i][j] = real(L[i][j]);              // move matrix to double array

   for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
     if (i==j) E[i][j] = 1.;
     else E[i][j] = 0.;                      // making unit matrix

   result1 = 0.;
   result2 = 0.;
   result3 = 0.;

   for (k=0; k < step_l; k++)
    {
     lambda = besk * k / step_l;             // integration area on step

     for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
       LL[i][j] = E[i][j]*lambda + L1[i][j];     // new matrix Lambda building

     det = (LL[0][0]*LL[1][1]*LL[2][2]) - (LL[0][2]*LL[1][1]*LL[0][2]) - (LL[0][0]*LL[1][2]*LL[1][2]);   // determinator
     Min00 = (LL[1][1]*LL[2][2])-(LL[1][2]*LL[1][2]);            // minors for diagonal elements:
     Min11 = (LL[0][0]*LL[2][2])-(LL[0][2]*LL[0][2]);
     Min22 =  LL[0][0]*LL[1][1];

     TLL = Min00;   // for x nomer = 1
     result1 += (pow((lambda/det),0.5)/det)* TLL;

     TLL = Min11;   // for y nomer = 2
     result2 += (pow((lambda/det),0.5)/det)* TLL;

     TLL = Min22;   // for z nomer = 3
     result3 += (pow((lambda/det),0.5)/det)* TLL;

    }
   result1 = result1*besk/step_l;
   result2 = result2*besk/step_l;
   result3 = result3*besk/step_l;


       F_x = A * result1*lat.betax/(beam.Emit[0]*(1.0+lat.alfax*lat.alfax));
       F_y = A * result2*lat.betay/(beam.Emit[1]*(1.0+lat.alfay*lat.alfay));
       F_z = A * ring.Energy.Gamma * result3/beam.Emit[2];
       //F_x = A * result1*lat.betax/beam.Emit[0];
       //F_y = A * result2*lat.betay/beam.Emit[1];
       //F_z = A * ring.Energy.Gamma2 * result3/beam.Emit[2];

}

void xIBS::B_M_diffusion(xLattice& lat, xBeam& beam, xRing& ring)
{
   doubleU ri(0,m_);
   doubleU A(0,s_^-1);

   double phi_Bx;
   double phi_By;
   doubleU Hx(0,m_);
   doubleU Hy(0,m_);

   matrixU L (3,3);

   double L1[3][3];
   double LL[3][3];
   double E[3][3];

   double result1, result2, result3, result4, result5;
   double besk;
   int step_l;

   double lambda = 0;
   double Min00 = 0;
   double Min11 = 0;
   double Min22 = 0;
   double det = 0;
   double TLL = 0;
   int i,j,k;

   if (beam.benum == BUNCHED)
   {  nb = 1;
      beam.CalcBunch();
      s_s = beam.s_s;
   }else
   {  if (beam.benum == COASTING)
         nb = 2;
      else
         nb = BUCKETnb;
      s_s = ring.Circ / (2 * (U_pi^0.5));
   }

   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;        // ion classical radius

        A =  U_c * ri * ri * beam.Emit[3]*B_M_log /
        (8. * U_pi * beam.Emit[0]*beam.Emit[1] * s_s * (beam.Emit[2]^0.5) *
        ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);

       phi_Bx = (lat.Dpx + (lat.alfax*lat.Dx(m_)/lat.betax(m_)))();
       phi_By = (lat.Dpy + (lat.alfay*lat.Dy(m_)/lat.betay(m_)))();
       Hx = lat.betax*( (lat.Dpx + ((lat.alfax*lat.Dx)/lat.betax) )^2) +(lat.Dx*lat.Dx/lat.betax);
       Hy = lat.betay*( (lat.Dpy + ((lat.alfay*lat.Dy)/lat.betay) )^2) +(lat.Dy*lat.Dy/lat.betay);

       L[0][0] = lat.betax(m_)/beam.Emit[0](m_);
       L[0][1] = 0;
       L[0][2] = ((-1.)*lat.betax(m_)*phi_Bx*ring.Energy.Gamma/beam.Emit[0](m_))();

       L[1][0] = 0;
       L[1][1] = lat.betay(m_)/beam.Emit[1](m_);
       L[1][2] = ((-1.)*lat.betay(m_)*phi_By*ring.Energy.Gamma/beam.Emit[1](m_))();

       L[2][0] = ((-1.)*ring.Energy.Gamma*lat.betax(m_)*phi_Bx/beam.Emit[0](m_))();
       L[2][1] = ((-1.)*ring.Energy.Gamma*lat.betay(m_)*phi_By/beam.Emit[1](m_))();
       L[2][2] = (ring.Energy.Gamma2*((Hx(m_)/beam.Emit[0](m_)) + (Hy(m_)/beam.Emit[1](m_)) + (1./beam.Emit[2])))();

       uplimit = real(L[0][0]);
       if (uplimit < real(L[1][1])) uplimit = real(L[1][1]);
       if (uplimit < real(L[2][2])) uplimit = real(L[2][2]);

       uplimit *= lam_limit;
       besk = uplimit;
       step_l = (int)step_lam;



   for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      L1[i][j] = real(L[i][j]);              // move matrix to double array

   for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
     if (i==j) E[i][j] = 1.;
     else E[i][j] = 0.;                      // making unit matrix

   result1 = 0.;
   result2 = 0.;
   result3 = 0.;
   result4 = 0.;
   result5 = 0.;

   for (k=0; k < step_l; k++)
    {
     lambda = besk * k / step_l;             // integration area on step

     for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
       LL[i][j] = E[i][j]*lambda + L1[i][j];     // new matrix Lambda building

     det = (LL[0][0]*LL[1][1]*LL[2][2]) - (LL[0][2]*LL[1][1]*LL[0][2]) - (LL[0][0]*LL[1][2]*LL[1][2]);   // determinator
     Min00 = (LL[1][1]*LL[2][2])-(LL[1][2]*LL[1][2]);            // minors for diagonal elements:
     Min11 = (LL[0][0]*LL[2][2])-(LL[0][2]*LL[0][2]);
     Min22 =  LL[0][0]*LL[1][1];

     TLL = Min11 + Min22;   // for xx nomer = 1
     result1 += (pow((lambda/det),0.5)/det)* TLL;

     TLL = Min00 + Min22;   // for yy nomer = 2
     result2 += (pow((lambda/det),0.5)/det)* TLL;

     TLL = Min11 + Min00;   // for zz nomer = 3
     result3 += (pow((lambda/det),0.5)/det)* TLL;

     TLL = LL[1][1]*LL[0][2];         // for xz nomer = 4
     result4 += (pow((lambda/det),0.5)/det)* TLL;

     TLL = LL[0][0]*LL[1][2];         // for yz nomer = 5
     result5 += (pow((lambda/det),0.5)/det)* TLL;
    }
   result1 = result1*besk/step_l;
   result2 = result2*besk/step_l;
   result3 = result3*besk/step_l;
   result4 = result4*besk/step_l;
   result5 = result5*besk/step_l;


       D_xx = A * result1;
       D_yy = A * result2;
       D_zz = A * ring.Energy.Gamma2 * result3;
       D_zx = A * ring.Energy.Gamma * result4;
       D_zy = A * ring.Energy.Gamma * result5;
}


// ========================== George Parzen Bi-gaussian theory =================

vectorU xIBS::GP_BiGauss(xLattice& lat, xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   doubleU ri(0,m_);

   doubleU Ac(0,s_^-1);
   doubleU At(0,s_^-1);
   doubleU Act(0,s_^-1);

   double phi_Bx;
   double phi_By;
   doubleU Hx(0,m_);
   doubleU Hy(0,m_);

   doubleU Ixx(s_^-1);
   doubleU Ixz(s_^-1);
   doubleU Iyy(s_^-1);
   doubleU Iyz(s_^-1);
   doubleU Izz(s_^-1);

   doubleU Ic_xx(s_^-1);
   doubleU Ic_xz(s_^-1);
   doubleU Ic_yy(s_^-1);
   doubleU Ic_yz(s_^-1);
   doubleU Ic_zz(s_^-1);

   doubleU It_xx(s_^-1);
   doubleU It_xz(s_^-1);
   doubleU It_yy(s_^-1);
   doubleU It_yz(s_^-1);
   doubleU It_zz(s_^-1);

   doubleU Ict_xx(s_^-1);
   doubleU Ict_xz(s_^-1);
   doubleU Ict_yy(s_^-1);
   doubleU Ict_yz(s_^-1);
   doubleU Ict_zz(s_^-1);

   matrixU Lc (3,3);
   matrixU Lt (3,3);
   matrixU Lct (3,3);

   double Nc[3];
   double Numc, Numt, Numct;
   vectorU Ec, Et, Ect;
   Et = beam.Emit;
   Ec = beam.Emit_FWHM(lat);
   beam.Courant_Snyder(lat);

      for(int i = 0; i < 3; i++)
      {  beam.Powel2(beam.Fitting[i], beam.Hyst3[i], beam.Sigmas, beam.Division);
         if (fabs(beam.Fitting[i][1]) > fabs(beam.Fitting[i][2]))
         {
            Ec[i] = beam.Emit[i] * beam.Fitting[i][2] * beam.Fitting[i][2];
            Et[i] = beam.Emit[i] * beam.Fitting[i][1] * beam.Fitting[i][1];
            Nc[i] = 1. / (1. + fabs
                  ( beam.Fitting[i][1] / beam.Fitting[i][2] *
                    beam.Fitting[i][3] / beam.Fitting[i][4] ) );
         }else
         {
            Ec[i] = beam.Emit[i] * beam.Fitting[i][1] * beam.Fitting[i][1];
            Et[i] = beam.Emit[i] * beam.Fitting[i][2] * beam.Fitting[i][2];
            Nc[i] = 1. / (1. + fabs
                  ( beam.Fitting[i][2] / beam.Fitting[i][1] *
                    beam.Fitting[i][4] / beam.Fitting[i][3] ) );
         }
      }
      double NN = beam.Emit[3]();
      Numc = NN*((Nc[0]+Nc[1]+Nc[2])/3.);
      Numt = NN - Numc;
      Numct = Numc*Numt/NN;

      Ect=Ec;
      for(int d = 0; d < 3; d++)
       Ect[d] = (Ec[d]* Et[d])/(Ec[d]+Et[d]);
/*
      double CoreNumber = 0;
      bool* InCore = new bool[beam.Number()];
      for (int i0 = 0; i0 < beam.Number(); i0++)
      {  bool InsideCore = true;
         for (int j = 0; j < 3; j++)
            if (Ec[j] * 2 < beam.Inv(i0,j))
               InsideCore = false;
         if (InsideCore)
         {  CoreNumber += 1;
            InCore[i0] = true;
         }else
            InCore[i0] = false;
      }

      for (int i1 = 0; i1 < beam.Number(); i1++)
         if (InCore[i1])
            IndividualKick(time, beam, i1, Rc, Ec, 1.);
         else
            IndividualKick(time, beam, i1, Rt, Et, 1.);

      delete []InCore;
*/


   ri  = ring.Energy.Z * ring.Energy.Z * U_rp / ring.Energy.A;

   Ac =  U_c * ri * ri * Numc*B_M_log /(8. * U_pi *Ec[0]*Ec[1] * s_s * (Ec[2]^0.5) *
        ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);

   At =  U_c * ri * ri * Numt*B_M_log /(8. * U_pi *Et[0]*Et[1] * s_s * (Et[2]^0.5) *
        ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);

   Act =  U_c * ri * ri * Numct*B_M_log /(8. * U_pi * Ect[0]*Ect[1] * s_s * (Ect[2]^0.5) *
        ring.Energy.Beta3*ring.Energy.Gamma2*ring.Energy.Gamma2);

       phi_Bx = (lat.Dpx + (lat.alfax*lat.Dx(m_)/lat.betax(m_)))();
       phi_By = (lat.Dpy + (lat.alfay*lat.Dy(m_)/lat.betay(m_)))();
       Hx = lat.betax*( (lat.Dpx + ((lat.alfax*lat.Dx)/lat.betax) )^2) +(lat.Dx*lat.Dx/lat.betax);
       Hy = lat.betay*( (lat.Dpy + ((lat.alfay*lat.Dy)/lat.betay) )^2) +(lat.Dy*lat.Dy/lat.betay);
//*************************
       Lc[0][0] = lat.betax(m_)/Ec[0](m_);
       Lc[0][1] = 0;
       Lc[0][2] = ((-1.)*lat.betax(m_)*phi_Bx*ring.Energy.Gamma/Ec[0](m_))();

       Lc[1][0] = 0;
       Lc[1][1] = lat.betay(m_)/Ec[1](m_);
       Lc[1][2] = ((-1.)*lat.betay(m_)*phi_By*ring.Energy.Gamma/Ec[1](m_))();

       Lc[2][0] = ((-1.)*ring.Energy.Gamma*lat.betax(m_)*phi_Bx/Ec[0](m_))();
       Lc[2][1] = ((-1.)*ring.Energy.Gamma*lat.betay(m_)*phi_By/Ec[1](m_))();
       Lc[2][2] = (ring.Energy.Gamma2*((Hx(m_)/Ec[0](m_)) + (Hy(m_)/Ec[1](m_)) + (1./Ec[2])))();
//*************************
       Lt[0][0] = lat.betax(m_)/Et[0](m_);
       Lt[0][1] = 0;
       Lt[0][2] = ((-1.)*lat.betax(m_)*phi_Bx*ring.Energy.Gamma/Et[0](m_))();

       Lt[1][0] = 0;
       Lt[1][1] = lat.betay(m_)/Et[1](m_);
       Lt[1][2] = ((-1.)*lat.betay(m_)*phi_By*ring.Energy.Gamma/Et[1](m_))();

       Lt[2][0] = ((-1.)*ring.Energy.Gamma*lat.betax(m_)*phi_Bx/Et[0](m_))();
       Lt[2][1] = ((-1.)*ring.Energy.Gamma*lat.betay(m_)*phi_By/Et[1](m_))();
       Lt[2][2] = (ring.Energy.Gamma2*((Hx(m_)/Et[0](m_)) + (Hy(m_)/Et[1](m_)) + (1./Et[2])))();
//*************************
       Lct[0][0] = lat.betax(m_)/Ect[0](m_);
       Lct[0][1] = 0;
       Lct[0][2] = ((-1.)*lat.betax(m_)*phi_Bx*ring.Energy.Gamma/Ect[0](m_))();

       Lct[1][0] = 0;
       Lct[1][1] = lat.betay(m_)/Ect[1](m_);
       Lct[1][2] = ((-1.)*lat.betay(m_)*phi_By*ring.Energy.Gamma/Ect[1](m_))();

       Lct[2][0] = ((-1.)*ring.Energy.Gamma*lat.betax(m_)*phi_Bx/Ect[0](m_))();
       Lct[2][1] = ((-1.)*ring.Energy.Gamma*lat.betay(m_)*phi_By/Ect[1](m_))();
       Lct[2][2] = (ring.Energy.Gamma2*((Hx(m_)/Ect[0](m_)) + (Hy(m_)/Ect[1](m_)) + (1./Ect[2])))();
//*************************
        uplimit = real(Lct[0][0]);
        if (uplimit < real(Lct[1][1])) uplimit = real(Lct[1][1]);
        if (uplimit < real(Lct[2][2])) uplimit = real(Lct[2][2]);
        uplimit *= lam_limit;

       Ic_xx = Ac * B_M_integral_I(Lc, 1);
       Ic_yy = Ac * B_M_integral_I(Lc, 2);
       Ic_zz = Ac * B_M_integral_I(Lc, 3);
       Ic_xz = Ac * B_M_integral_I(Lc, 4);
       Ic_yz = Ac * B_M_integral_I(Lc, 5);

       It_xx = At * B_M_integral_I(Lt, 1);
       It_yy = At * B_M_integral_I(Lt, 2);
       It_zz = At * B_M_integral_I(Lt, 3);
       It_xz = At * B_M_integral_I(Lt, 4);
       It_yz = At * B_M_integral_I(Lt, 5);

       Ict_xx = Act * B_M_integral_I(Lct, 1);
       Ict_yy = Act * B_M_integral_I(Lct, 2);
       Ict_zz = Act * B_M_integral_I(Lct, 3);
       Ict_xz = Act * B_M_integral_I(Lct, 4);
       Ict_yz = Act * B_M_integral_I(Lct, 5);

       Ixx = (Numc*Ic_xx/beam.Emit[3]) + (Numt*It_xx/beam.Emit[3]) + 2.*Ict_xx;
       Iyy = (Numc*Ic_yy/beam.Emit[3]) + (Numt*It_yy/beam.Emit[3]) + 2.*Ict_yy;
       Izz = (Numc*Ic_zz/beam.Emit[3]) + (Numt*It_zz/beam.Emit[3]) + 2.*Ict_zz;
       Ixz = (Numc*Ic_xz/beam.Emit[3]) + (Numt*It_xz/beam.Emit[3]) + 2.*Ict_xz;
       Iyz = (Numc*Ic_yz/beam.Emit[3]) + (Numt*It_yz/beam.Emit[3]) + 2.*Ict_yz;

   rates[0] = ( (Hx*ring.Energy.Gamma2*Izz) -
               (2.*lat.betax*phi_Bx*ring.Energy.Gamma*Ixz) +
               (lat.betax*Ixx) )/beam.Emit[0];

   rates[1] = ( (Hy*ring.Energy.Gamma2*Izz) -
               (2.*lat.betay*phi_By*ring.Energy.Gamma*Iyz) +
               (lat.betax*Iyy) )/beam.Emit[1];

   rates[2] = ring.Energy.Gamma2*Izz/beam.Emit[2];

   return rates;
}







