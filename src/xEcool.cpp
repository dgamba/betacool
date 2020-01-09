//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xEcool.h"
#include "xDynamic.h"
//---------------------------------------------------------------------------
xEcool::xEcool()
{  
   ID = 1;
   EL_FORCES = true;
   OpticName = "XECOOL";
   EffectName = "XECOOL";
   SectionNumber = 1;
};

int xEcool::OnGet()
{
   E_model       = Data.OnGetI(50,1,0);
   integr_step   = Data.OnGetI(50,2,2);
   SectionNumber = Data.OnGetI(50,5,1);
   FringeField   = Data.OnGetB(50,6,false);

   I_model   = Data.OnGetI(51,1,0);
   step_tr   = Data.OnGetI(51,2,8);
   step_long = Data.OnGetI(51,3,8);
   num_MC    = Data.OnGetI(51,4,1000);

   return 0;
}                                   

int xEcool::OnSet()
{
   Lattice.betax = Data.OnGet(52,1,1);
   Lattice.betay = Data.OnGet(52,2,1);
   Lattice.alfax = Data.OnGet(52,3,0);
   Lattice.alfay = Data.OnGet(52,4,0);
   Lattice.Dx    = Data.OnGet(52,5,0);
   Lattice.Dy    = Data.OnGet(52,6,0);
   Lattice.Dpx   = Data.OnGet(52,7,0);
   Lattice.Dpy   = Data.OnGet(52,8,0);

   return 0;
}
//------------------------------------------------------

vectorU xEcool::Rates(xTime&t, xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   doubleU Init(m_);
   doubleU Fin(m_);

   if (beam.benum == BUNCHED)
      beam.CalcBunch();
   Init = iEbeam.Initial[4];
   Fin = iEbeam.Final[4];

   for (int jj = 0; jj < iEbeam.BunchNumber; jj++)
   {
      switch(I_model)
      {
       case 0: rates += SingleParticle(t,beam,ring); break;
       case 1: rates += MonteCarlo(t,beam,ring); break;
      }

      iEbeam.Initial[4] += iEbeam.Bdistance;
      iEbeam.Final[4] += iEbeam.Bdistance;

   }
   iEbeam.Initial[4] = Init;
   iEbeam.Final[4] = Fin;

   return rates * SectionNumber;
}
//------------------------------------------------------
//--  here the vector of force components is organizing
//--  for usage in iMethods call (see below
//------------------------------------------------------

vectorU SumForces(xTime& t, vectorU& X)
{
 return iEcool.GetForces(t, iRing.Energy, X) +
        iDrift.GetForces(t, iRing.Energy, X);
}
//------------------------------------------------------
// these functions are now in xdistributor:

vectorU xEcool::Courant_Snyder(vectorU X, xLattice& lat, doubleU Bs, bool bunched)
{
   vectorU Inv;    // result vector
   Inv[0]._(m_);
   Inv[1]._(m_);

   doubleU gammax(m_^-1);
   gammax = (1. + (lat.alfax^2))/lat.betax;
   doubleU gammay(m_^-1);
   gammay = (1. + (lat.alfay^2))/lat.betay;

   doubleU Xb(m_);
   doubleU Xs(U1_);
   doubleU Yb(m_);
   doubleU Ys(U1_);

   Xb = X[0] - lat.Dx * X[5];
   Xs = X[1] - lat.Dpx* X[5];
   Yb = X[2] - lat.Dy * X[5];
   Ys = X[3] - lat.Dpy* X[5];

   Inv[0] = ((lat.betax*Xs*Xs)+(2.*lat.alfax*Xb*Xs)+(gammax*Xb*Xb));
   Inv[1] = ((lat.betay*Ys*Ys)+(2.*lat.alfay*Yb*Ys)+(gammay*Yb*Yb));
   Inv[2] = X[5]* X[5];

   if (bunched)
      Inv[2] += (X[4] * X[4]) / (Bs * Bs);

   //Inv[3] = U_Acos(X[0]*X[0] / (Inv[0]*lat.betax*2));
   //Inv[4] = U_Acos(X[2]*X[2] / (Inv[1]*lat.betay*2));
   return Inv;
}

//------------------------------------------------------

xLattice xEcool::FinalLattice(xLattice& lat)
{
   xLattice ecool_lat_n;
   doubleU gammax(m_^-1);
   gammax = (1 + (Lattice.alfax^2))/Lattice.betax;
   doubleU gammay(m_^-1);
   gammay = (1 + (Lattice.alfay^2))/Lattice.betay;

   ecool_lat_n = Lattice;
   if(E_model > 0)
   {
      ecool_lat_n.betax = Lattice.betax - (2*iEbeam.CoLength*Lattice.alfax) + (iEbeam.CoLength*iEbeam.CoLength*gammax);
      ecool_lat_n.betay = Lattice.betay - (2*iEbeam.CoLength*Lattice.alfay) + (iEbeam.CoLength*iEbeam.CoLength*gammay);

      ecool_lat_n.alfax = Lattice.alfax - (iEbeam.CoLength*gammax);
      ecool_lat_n.alfay = Lattice.alfay - (iEbeam.CoLength*gammay);
   }

   return ecool_lat_n;
}

//------------------------------------------------------
//#define CHECKCOOL // check ion distribution in cooler

// 20.11.06 Barrier bucket
vectorU xEcool::SingleParticle(xTime&t, xBeam& beam, xRing& ring)
{
  vectorU rates(s_^-1); // result rates sec^-1
  vectorU rateT(U1_);   // temporary rates (unitless)
  vectorU X(m_,U1_ );   // initial particle coord
  vectorU X1(m_,U1_ );  // initial temp particle coord
  vectorU Z(m_,U1_ );   // vector of Fr. force components for "Thin lens" model
  doubleU total_number = 0;   // total number of steps (over all phazes)
  vectorU EmitN;              // vector of new Emittances (after effect pass)
  EmitN[0]._(m_);
  EmitN[1]._(m_);
  doubleU T_int(s_);

  xTime tE(t);    // local time variable for Cooling section
  tE.DistStep(iEbeam.CoLength/integr_step); // sets step increment for integration methods

  xLattice ecool_lat_n = FinalLattice(Lattice);

  double step = 2.0*M_PI/step_tr;   // step over betatron phazes

  double step_l;
  int step_p;
  if (beam.benum == BUNCHED)
  {  step_l = 2.0*M_PI/step_long;   // step over long phazes
     step_p = step_long;
  }else
  {if (beam.benum == BUCKET)
  {  beam.CalcBucket(beam.Emit[2]);
     step_l = 1.0/step_long;   // step over long phazes
     step_p = step_long;
  }else
  {
     step_l = 2.0; // 2.0 /step_long;
     step_p = 2;   // step_long + 1;
  }
  }
#ifdef CHECKCOOL
      for (int n0 = 0; n0 < 4; n0++)
      if(iDraw.RealSpace[n0].Graf.Enabled)
         iDraw.RealSpace[n0].Graf.Reset(0);
#endif

  for(int i = 0; i < step_tr; i++)
  {
     X1[2] = ((beam.Emit[1] *2.*Lattice.betay)^0.5)*sin(i*step);
     X1[3] = ((beam.Emit[1] *2./Lattice.betay)^0.5)*(cos(i*step) - (Lattice.alfay * sin(i*step)));

     for(int j = 0; j < step_tr; j++)
     {
	     X1[0] = ((beam.Emit[0]*2.*Lattice.betax)^0.5)*sin(j*step);
   	  X1[1] = ((beam.Emit[0]*2./Lattice.betax)^0.5)*(cos(j*step) - (Lattice.alfax * sin(j*step)));

        for(int k=0; k < step_p; k++)
        {
            if (beam.benum == BUNCHED)
            {
                X[4] = ((2.*beam.Emit[2])^0.5)*sin(k*step_l)*xLattice::B_s;
                X[5] = ((2.*beam.Emit[2])^0.5)*cos(k*step_l);
            }else
            {if (beam.benum == BUCKET)
            {
                T_int = k*step_l*beam.T_s;
                beam.GetBucket(X, T_int);
            }else
            {
            X[5] = ((beam.Emit[2])^0.5)*(k*step_l - 1.);
            X[4] = 0;
            }
            }

            X[0] = X1[0] + (X[5]*Lattice.Dx);
            X[1] = X1[1] + (X[5]*Lattice.Dpx);
            X[2] = X1[2] + (X[5]*Lattice.Dy);
            X[3] = X1[3] + (X[5]*Lattice.Dpy);
#ifdef CHECKCOOL
            double x, s = X[4]();
            if (iBeam == COASTING)
               s -= iBeam.CellSize() * floor(s / iBeam.CellSize());
            for (int n3 = 0; n3 < 2; n3++)
            {  if (iDraw.RealSpace[n3].x_axis == 4)
                  x = s;
               else
                  x = X[iDraw.RealSpace[n3].x_axis]();
               iDraw.RealSpace[n3].Graf.Point(x, X[iDraw.RealSpace[n3].y_axis]());
            }
#endif
            rates[3] += New_Coord(tE, beam, ring, X);
            EmitN = Courant_Snyder(X, ecool_lat_n, xLattice::B_s, beam.benum)/2.;
            if (beam.benum == COASTING)
               EmitN[2] *= 2.;
            if (beam.benum == BUCKET)
            {
               EmitN[2] = beam.Bdp2(X[4],X[5]);
            }
            rateT[0] += (EmitN[0] - beam.Emit[0])/beam.Emit[0];
            rateT[1] += (EmitN[1] - beam.Emit[1])/beam.Emit[1];
            rateT[2] += (EmitN[2] - beam.Emit[2])/beam.Emit[2];

            total_number += 1;
#ifdef CHECKCOOL
            s = X[4]();
            if (iBeam.benum == COASTING)
               s -= iBeam.CellSize() * floor(s / iBeam.CellSize());
            for (int n3 = 2; n3 < 4; n3++)
            {  if (iDraw.RealSpace[n3].x_axis == 4)
                  x = s;
               else
                  x = X[iDraw.RealSpace[n3].x_axis]();
               iDraw.RealSpace[n3].Graf.Point(x, X[iDraw.RealSpace[n3].y_axis]());
            }
#endif
        }
     }
  }
  rates[0] = rateT[0]/(ring.Trev* total_number);
  rates[1] = rateT[1]/(ring.Trev*total_number);
  rates[2] = rateT[2]/(ring.Trev*total_number);
  rates[3] /= total_number;
#ifdef CHECKCOOL
      for (int n8 = 0; n8 < 4; n8++)
         iDraw.RealSpace[n8].Graf.Save();
#endif
  return rates;
}
//
//------------------------------------------------------

doubleU xEcool::LifeTime(xTime&t, xBeam& beam, xRing& ring)
{
 doubleU lifetime(s_^-1);
 static doubleU alfa(3.02e-13,(cm_^3)/s_);
 static doubleU A(2.11e-22,(cm_^2));
 static doubleU eV(1.0, eV_);
 doubleU v_l(m_/s_);  // step over longitudinal electron velocity
 doubleU v_t(m_/s_);  // step over transverse electron velocity
 int i, j, n;
 doubleU Maxw;
 doubleU d_x((cm_^4)/(s_^2));
 doubleU d_z(cm_/s_);
 doubleU sigma(cm_^2);
 doubleU hnu0(eV_);
 doubleU Ener(eV_);
 doubleU alfa_und((cm_^3)/s_);
 doubleU result(3.02e-13,(cm_^3)/s_);

  lifetime = iEbeam.CoLength* iEbeam.F.n_e/(iEbeam.e_Energy.Gamma*ring.Circ)*alfa;
  result = (ring.Energy.Z^2)*((eV/iEbeam.F.Ttemp)^0.5)*(U_Ln(11.32*ring.Energy.Z/((iEbeam.F.Ttemp/eV)^0.5))
           + (0.14*((iEbeam.F.Ttemp/(eV*(ring.Energy.Z^2)))^(1./3.))))*alfa;
  lifetime *= (ring.Energy.Z^2)*((eV/iEbeam.F.Ttemp)^0.5)*(U_Ln(11.32*ring.Energy.Z/((iEbeam.F.Ttemp/eV)^0.5))
                                                           + (0.14*((iEbeam.F.Ttemp/(eV*(ring.Energy.Z^2)))^(1./3.))));
  if(iEbeam.F.undulator)
    {
            hnu0=13.6*eV*ring.Energy.Z*ring.Energy.Z;
            v_t = (iEbeam.F.V_tr_e * 3)/iForce.dt;
            v_l = (iEbeam.F.V_long_e * 3)/iForce.dl;
            n = 0;
            d_x = 0;
            d_z = 0;

            if(iEbeam.F.V_und()>(6*iEbeam.F.V_tr_e()))
            {
            Ener = U_me*(iEbeam.F.V_und*iEbeam.F.V_und)/2;
            sigma = A*(hnu0/Ener)*(0.5*U_Ln(hnu0/Ener)+0.525*((Ener/hnu0)^(1/3))+0.1402);
            alfa_und = iEbeam.F.V_und*sigma;
            }
            else
            {
            for (i = 1-iForce.dl; i < iForce.dl; i++)
            for (j = 1; j < iForce.dt; j++)

                {

                Maxw = U_Exp(((v_t*j+iEbeam.F.V_und)*(v_t*j+iEbeam.F.V_und))/(iEbeam.F.V_tr_e*iEbeam.F.V_tr_e*2.)+
                             ((v_l*i)*(v_l*i))/(iEbeam.F.V_long_e*iEbeam.F.V_long_e*2.));

                Ener = U_me*((v_t*j+iEbeam.F.V_und)*(v_t*j+iEbeam.F.V_und)+(v_l*i)*(v_l*i))/2;

                sigma = A*(hnu0/Ener)*(0.5*U_Ln(hnu0/Ener)+0.525*((Ener/hnu0)^(1/3))+0.1402);

                d_x += sigma*(((v_t*j+iEbeam.F.V_und)*(v_t*j+iEbeam.F.V_und)+(v_l*i)*(v_l*i))^0.5)*(v_t*j)/Maxw;

                d_z += (v_t*j)/Maxw;

                n++;

                }
            //alfa_und = d_x*18/(iEbeam.F.V_tr_e*((2. * U_pi)^0.5)*n);
              alfa_und = d_x/d_z;
             }
            lifetime = iEbeam.CoLength* iEbeam.F.n_e/(iEbeam.e_Energy.Gamma*ring.Circ)*alfa_und;

    }

 return lifetime;
}
//------------------------------------------------------
doubleU xEcool::New_Coord(xTime&tE, xBeam& beam, xRing& ring, vectorU&X)
{
   doubleU liferate(0, s_^-1);
   vectorU Force;
   int m;

   doubleU ro(ring.Energy.Momentum*U_c/(ring.Energy.Z*U_e*iEbeam.F.mfield));
   //double cm = ro(cm_);
   if (FringeField)
   {  X[1] += X[2] / (ro * 2.);
      X[3] -= X[0] / (ro * 2.);
   }
   
   tE.so = 0;
   switch(E_model)
   {
    case 0:
      tE.so = iEbeam.CoLength/2.;
      Force =  GetForces (tE, ring.Energy, X);
      X[1] +=  Force[1] * iEbeam.CoLength;
      X[3] +=  Force[3] * iEbeam.CoLength;
      X[5] +=  Force[5] * iEbeam.CoLength;
      if (Loss && iEbeam.inside)
         liferate -= LifeTime(tE, beam, ring);
      break;

    case 1:
      for(m = 0; m < integr_step; m++)
      {
         X = Euler(tE, tE.ds, X, SumForces);
         if (Loss && iEbeam.inside)
            liferate -= LifeTime(tE, beam, ring);
         ++tE;
      }
      liferate /= integr_step;
      break;

    case 2:
      for(m = 0; m < integr_step; m++)
      {
         X = Kutta4(tE, tE.ds, X, SumForces);
         if (Loss && iEbeam.inside)
            liferate -= LifeTime(tE, beam, ring);
         ++tE;
      }
      liferate /= integr_step;
      break;
   }
   if (FringeField)
   {  X[1] -= X[2] / (ro * 2.);
      X[3] += X[0] / (ro * 2.);
   }
   return liferate;
}

vectorU xEcool::MonteCarlo(xTime&t, xBeam& beam, xRing& ring)
{
  vectorU em(U1_); // vector of emittances of generated beam
  em[0]._(m_);
  em[1]._(m_);
  vectorU em1(U1_); // vector of new emittances of generated beam
  em[0]._(m_);
  em[1]._(m_);
  // ----- 06.07
  beam.EmitDef = 0;

  xLattice ecool_lat_n = FinalLattice(Lattice);

  xTime tE(t);    // local time variable for Cooling section
  tE.DistStep(iEbeam.CoLength/integr_step); // sets step increment for integration methods

  vectorU rates(s_^-1); // result rates sec^-1
  vectorU Z(m_,U1_ );   // vector of Fr. force components for "Thin lens" model
  vectorU X;            // vector of coordinate;

   beam.Number(num_MC);
   if (beam.benum == BUNCHED)
      beam.CalcBunch();
   beam.Distribution(beam.Emit,Lattice);
   em = beam.Emit_Def(Lattice);

   for (int i = 0; i < num_MC; i++)
   {
      X = beam(i);
      rates[3] += New_Coord(tE, beam, ring, X);
      beam(i, X);
   }

   em1 = beam.Emit_Def(ecool_lat_n);

   for (int j = 0; j < 3; j++)
      rates[j] = (em1[j]-em[j])/(em[j]*ring.Trev);
   rates[3] /= num_MC;
   return rates;
}
//------------------------------------------------------

vectorU xEcool::GetForces (xTime&t, U_Energy& e, vectorU ion)
{
   doubleU E_e (M_6*eV_ * 0.5110034);
   vectorU forces(U1_, m_^-1);
   vectorU forces_new(U1_, m_^-1);
   vectorU ion_new(m_, U1_);
   int i = 0;

   //iEbeam.inside = false;
   ion_new = iEbeam.Ibeam2Ebeam(t, ion);  

   if (!iEbeam.inside || (Multi == 0))
   {
      for (i = 0; i < 6; i++)
         forces_new[i] = 0;
      return forces_new;
   }

   iForce.Vtr = (((ion_new[1]*ion_new[1]) + (ion_new[3]*ion_new[3]))^0.5)*
                  (iEbeam.e_Energy.Beta * U_c * iEbeam.e_Energy.Gamma);
   iForce.v[2] = ion_new[5]*(iEbeam.e_Energy.Beta * U_c);
   iForce.v[0] = ion_new[1]*(iEbeam.e_Energy.Beta * U_c * iEbeam.e_Energy.Gamma);
   iForce.v[1] = ion_new[3]*(iEbeam.e_Energy.Beta * U_c * iEbeam.e_Energy.Gamma);

   switch(iForce.type)
   {
    case 0: iForce.Budker(iEbeam.F); break;
    case 1: iForce.NonMag(iEbeam.F); break;
    case 2: iForce.DerSkr(iEbeam.F); break;
    case 3: iForce.Parhom(iEbeam.F); break;
    case 4: iForce.Toepffer(iEbeam.F); break;
    case 5: iForce.Table(iEbeam.F); break;
    case 6: iForce.D3(iEbeam.F); break;
   }
   if(iForce.type == 6)
   {
     iForce.f[0] /= iEbeam.e_Energy.Gamma;
     iForce.f[1] /= iEbeam.e_Energy.Gamma;
   }else
   {
      iForce.Ftr /= iEbeam.e_Energy.Gamma;
      if (iForce.Ftr() == 0)
      {
       iForce.f[0] = 0;
       iForce.f[1] = 0;
      }
      else
      {
       iForce.f[0] = iForce.Ftr*ion_new[1]/(((ion_new[1]*ion_new[1]) + (ion_new[3]*ion_new[3]))^0.5);
       iForce.f[1] = iForce.Ftr*ion_new[3]/(((ion_new[1]*ion_new[1]) + (ion_new[3]*ion_new[3]))^0.5);
      }
   }
   forces[1] = iForce.f[0]/(e.Velocity * e.Momentum);
   forces[3] = iForce.f[1]/(e.Velocity * e.Momentum);
   forces[5] = iForce.f[2]/(e.Velocity * e.Momentum);
/*
   forces[0] = ion[1];
   forces[2] = ion[3];
   forces[4] = ion[5] / (e.Gamma^2);
*/
   forces_new = iEbeam.Ebeam2Ibeam(forces);  // ????????????????????????????????????????????????????????

   return forces_new;
}

void xEcool::Kick (xTime&time, xBeam&beam, xRing&ring)
{
   beam.RMS(beam.Emit,Lattice);   
   iDraw.InEcool.Reset(0);
   xTime tE(iTime);   // local time variable for Cooling section
   tE.DistStep(iEbeam.CoLength/integr_step);
   double loss;
   vectorU X;
   doubleU Init(m_);
   doubleU Fin(m_);
   doubleU lifetime(s_^-1);

//#pragma omp parallel for - not ok
   for (int j = 0; j < beam.Number(); j++)
   if	(!iEbeam.LongPulse || ((iBeam(j,4)>iRing.Circ*iEbeam.PulseFrom)&&(iBeam(j,4)<iRing.Circ*iEbeam.PulseUpTo)))
   {
      Init = iEbeam.Initial[4];
      Fin = iEbeam.Final[4];
      for (int jj = 0; jj < iEbeam.BunchNumber; jj++)
      {
         X = beam(j);
         lifetime = New_Coord(tE, beam, ring, X);
         for (int k = 1; k < 6; k +=2)
         {  beam(j, k, beam(j,k).Saturation(
              ((X[k] - beam(j,k)) * time.dt * SectionNumber / ring.Trev),
              iEbeam.Shift[k]));
         }

         if (Loss)
         {  loss = (-time.dt * lifetime * SectionNumber)(U1_);
            beam.Loss(Lattice, j, loss,false);
         }
         iEbeam.Initial[4] += iEbeam.Bdistance;
         iEbeam.Final[4] += iEbeam.Bdistance;

      }
      iEbeam.Initial[4] = Init;
      iEbeam.Final[4] = Fin;
   }
}
