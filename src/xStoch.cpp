//---------------------------------------------------------------------------
#include "stdafx.h"
#include "xEffect.h"
#include "xStoch.h"
//---------------------------------------------------------------------------

xStochastic::xStochastic()
{
   ID = 8;
   //Tag = 5;
   EffectName = "XSTOCHASTIC";

   //Horizontal chain

   HF_min._(G_9*Hz_);
   HF_max._(G_9*Hz_);
   HF_0._(G_9*Hz_);
   HBand._(G_9*Hz_);
   HG._(U1_);
   HG_opt._(U1_);
   HG_log._(U1_);
   H_l._(cm_);
   //PU parameters

   HPU_N._(U1_);
   HPU_S._(U1_);
   HPU_al._(m_);
   HPU_beta._(m_);
   HPU_w._(cm_);
   HPU_h._(cm_);
   //Kicker parameters
   HK_N._(U1_);
   HK_S._(U1_);
   HK_al._(m_);
   HK_beta._(m_);
   HK_w._(cm_);
   HK_h._(cm_);

   HP_th._(W_);
   HP_Sch._(W_);
   HP_tot._(W_);

   HE_inf._(m_);
   Htau._(Hz_);
   Htau_opt._(Hz_);
   Hx._(U1_);

   HA1._(U1_);
   HA2._(U1_);
   HA3._(U1_);

   //Vertical chain
   VF_min._(G_9*Hz_);
   VF_max._(G_9*Hz_);
   VF_0._(G_9*Hz_);
   VBand._(G_9*Hz_);
   VG._(U1_);
   VG_opt._(U1_);
   VG_log._(U1_);
   V_l._(cm_);
   //PU parameters

   VPU_N._(U1_);
   VPU_S._(U1_);
   VPU_al._(m_);
   VPU_beta._(m_);
   VPU_w._(cm_);
   VPU_h._(cm_);
   //Kicker parameters
   VK_N._(U1_);
   VK_S._(U1_);
   VK_al._(m_);
   VK_beta._(m_);
   VK_w._(cm_);
   VK_h._(cm_);

   VP_th._(W_);
   VP_Sch._(W_);
   VP_tot._(W_);

   VE_inf._(m_);
   Vtau._(Hz_);
   Vtau_opt._(Hz_);
   Vx._(U1_);

   VA1._(U1_);
   VA2._(U1_);
   VA3._(U1_);

   //Longitudinal chain
   LF_min._(G_9*Hz_);
   LF_max._(G_9*Hz_);
   LF_0._(G_9*Hz_);
   LBand._(G_9*Hz_);
   LG._(U1_);
   LG_opt._(U1_);
   LG_log._(U1_);
   L_l._(cm_);
   //PU parameters
   LPU_N._(U1_);
   LPU_al._(m_);
   //Kicker parameters
   LK_N._(U1_);
   LK_al._(m_);

   LP_th._(W_);
   LP_Sch._(W_);
   LP_tot._(W_);

   LE_inf._(U1_);
   kappa._(U1_);
   Ltau._(Hz_);
   Ltau0._(Hz_);
   Ltau_opt._(Hz_);
   A._(Hz_);
   B._(eV_*Hz_);

   LA1._(U1_);
   LA2._(U1_);
   LA3._(U1_);

   //Loop geometry and other parameters
   T_PU._(K_);
   T_PA._(K_);
   Z._(Ohm_);
   Msigma._(U1_);
   CombLoss._(U1_);
   Losses._(U1_);

   
   LBand._(G_9*Hz_);

   //     06.07
   Friction_c._(Hz_);
   Diffusion_c._(U1_);
   Rgain._(Ohm_);
   C_tilda._(m_*(s_^-1));
   //12.12.2007

   RLF_min._(G_9*Hz_);  //minimum frequency
   RLF_max._(G_9*Hz_);  //maximum frequency

   T_0._(s_);    //revolution period
   T_n._(U1_);      //time delay in T_0
   eta._(U1_);

}
//---------------------------------------------------------------------------
int xStochastic::OnGet()
{
   int i;
   doubleU y(1, cm_);
   double x;

        //Horizontal chain
        H_use = Data.OnGetI(5,1,0);
        HF_min = Data.OnGet(5,2,0);
        HF_max = Data.OnGet(5,3,0);
        HF_0 = 0.5*(HF_max + HF_min);
        HBand = HF_max - HF_min;

        i=(int) Data.OnGet(5,4,0);
        if(i)
            {
            HG_log = Data.OnGet(5,7,0);
            HG = U_Pow(10.0,HG_log()/20.0);
            }
        else
            {
            HG = Data.OnGet(5,6,0);
            HG_log = 20.0*U_Log(HG);
            }
        H_l = Data.OnGet(5,8,0);
        //Pickup parameters
        HPU_w = Data.OnGet(5,9,0);
        HPU_h = Data.OnGet(5,10,0);
        HPU_N = Data.OnGet(5,11,0);
        HPU_beta = Data.OnGet(5,13,0);
        HPU_al = HPU_N * (y + H_l);
        x = M_PI * HPU_w()/(2.0 * HPU_h());
        HPU_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));
        //Kicker parameters
        HK_w = Data.OnGet(5,16,0);
        HK_h = Data.OnGet(5,17,0);
        HK_N = Data.OnGet(5,18,0);
        HK_beta = Data.OnGet(5,20,0);
        HK_al = HK_N * (y + H_l);
        x = M_PI * HK_w()/(2.0 * HK_h());
        HK_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));
        //other parameters

        Average(HF_min, HF_max, H_l);

        //Vertical chain
        V_use = Data.OnGetI(6,1,0);
        VF_min = Data.OnGet(6,2,0);
        VF_max = Data.OnGet(6,3,0);
        VF_0 = 0.5*(VF_max + VF_min);
        VBand = VF_max - VF_min;

        i=(int) Data.OnGet(6,4,0);
        if(i)
            {
            VG_log = Data.OnGet(6,7,0);
            VG = U_Pow(10.0,VG_log()/20.0);
            }
        else
            {
            VG = Data.OnGet(6,6,0);
            VG_log = 20.0*U_Log(VG);
            }
        V_l = Data.OnGet(6,8,0);
        //Pickup parameters
        VPU_w = Data.OnGet(6,9,0);
        VPU_h = Data.OnGet(6,10,0);
        VPU_N = Data.OnGet(6,11,0);
        VPU_beta = Data.OnGet(6,13,0);
        VPU_al = VPU_N * (y + V_l);
        x = M_PI * VPU_w()/(2.0 * VPU_h());
        VPU_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));
        //Kicker parameters
        VK_w = Data.OnGet(6,16,0);
        VK_h = Data.OnGet(6,17,0);
        VK_N = Data.OnGet(6,18,0);
        VK_beta = Data.OnGet(6,20,0);
        VK_al = VK_N * (y + V_l);
        x = M_PI * VK_w()/(2.0 * VK_h());
        VK_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));
        //other parameters

        Average2(VF_min, VF_max, V_l);

        //Longitudinal chain
        L_use = Data.OnGetI(7,1,0);
        LF_min = Data.OnGet(7,2,0);
        LF_max = Data.OnGet(7,3,0);
        LF_0 = 0.5*(LF_max + LF_min);
        LBand = LF_max - LF_min;

        i=(int) Data.OnGet(7,4,0);
        if(i)
            {
            LG_log = Data.OnGet(7,7,0);
            LG = U_Pow(10.0,LG_log()/20.0);
            }
        else
            {
            LG = Data.OnGet(7,6,0);
            LG_log = 20.0*U_Log(LG);
            }
        L_l = Data.OnGet(7,8,0);
        //Pickup parameters
        LPU_N = Data.OnGet(7,11,0);
        LPU_al = LPU_N * (y + L_l);
        //Kicker parameters
        LK_N = Data.OnGet(7,18,0);
        LK_al = LK_N * (y + L_l);
        //other parameters
        Average3(LF_min, LF_max, L_l);

        //Loop geometry and other parameters
        Msigma = Data.OnGet(8,1,0);
        T_PU = Data.OnGet(8,2,0);
        T_PA = Data.OnGet(8,3,0);
        Z = Data.OnGet(8,4,0);
        CombLoss = Data.OnGet(8,5,0);
        Losses = Data.OnGet(8,6,0);

        HP_th = (T_PU + T_PA)*HG*HG*HBand;
        VP_th = (T_PU + T_PA)*VG*VG*VBand;
        LP_th = (T_PU + T_PA)*LG*LG*LBand/3.0;


   //             06.07

       RHIC = Data.OnGetB(7,29,false);
       Friction_c = Data.OnGet(90,1,0.001);
       Diffusion_c = Data.OnGet(90,2,0.001);
       Rgain =  Data.OnGet(90,3,0.5);
       alpha =  Data.OnGet(90,4,0.1);
       beta =   Data.OnGet(90,5,0.9);
//12.12.2007
       First = true;
       RLF_min =  Data.OnGet(90,6,0.5);
       RLF_max =  Data.OnGet(90,7,0.1);
       N_sub =   Data.OnGetI(90,8,15);
       n_series = Data.OnGetI(90,9,1);
       if(n_series == 0) n_series = 2;
       T_n = Data.OnGet(90,10,0.66);
//12.12.2007
        return 0;
}
int xStochastic::OnRun()
{
   int i;
   //doubleU y(1, cm_);
   double x;

        //Horizontal chain
        H_use = Data.OnGetI(5,1,0);

        i=(int) Data.OnGet(5,4,0);
        if(i)
            {
            HG_log = Data.OnGet(5,7,0);
            HG = U_Pow(10.0,HG_log()/20.0);
            }
        else
            {
            HG = Data.OnGet(5,6,0);
            HG_log = 20.0*U_Log(HG);
            }
        HPU_h = Data.OnGet(5,10,0);
        x = M_PI * HPU_w()/(2.0 * HPU_h());
        HPU_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));

        HK_h = Data.OnGet(5,17,0);
        x = M_PI * HK_w()/(2.0 * HK_h());
        HK_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));

        Data.OnSet(5,6,HG(U1_));
        Data.OnSet(5,7,HG_log(U1_));

       //Vertical chain
        V_use = Data.OnGetI(6,1,0);

        i=(int) Data.OnGet(6,4,0);
        if(i)
            {
            VG_log = Data.OnGet(6,7,0);
            VG = U_Pow(10.0,VG_log()/20.0);
            }
        else
            {
            VG = Data.OnGet(6,6,0);
            VG_log = 20.0*U_Log(VG);
            }

        VPU_h = Data.OnGet(6,10,0);
        x = M_PI * VPU_w()/(2.0 * VPU_h());
        VPU_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));

        VK_h = Data.OnGet(6,17,0);
        x = M_PI * VK_w()/(2.0 * VK_h());
        VK_S = 2.0 * (exp(x) - exp(-1.0*x))/(exp(x) + exp(-1.0*x));

        Data.OnSet(6,6,VG(U1_));
        Data.OnSet(6,7,VG_log(U1_));

         //Longitudinal chain
        L_use = Data.OnGetI(7,1,0);

        i=(int) Data.OnGet(7,4,0);
        if(i)
            {
            LG_log = Data.OnGet(7,7,0);
            LG = U_Pow(10.0,LG_log()/20.0);
            }
        else
            {
            LG = Data.OnGet(7,6,0);
            LG_log = 20.0*U_Log(LG);

        Data.OnSet(7,6,LG(U1_));
        Data.OnSet(7,7,LG_log(U1_));

        //OnSet();
            }
   return 0;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int xStochastic::OnSet()
{
        //Output horizontal
        Data.OnSet(5,5,Hx(U1_));
        Data.OnSet(5,6,HG(U1_));
        Data.OnSet(5,7,HG_log(U1_));
        Data.OnSet(5,27,HG_opt(U1_));
        Data.OnSet(5,12,HPU_al(m_));
        Data.OnSet(5,19,HK_al(m_));
        Data.OnSet(5,14,0.5*HPU_S(U1_));
        Data.OnSet(5,21,0.5*HK_S(U1_));
        Data.OnSet(5,22,HP_th(W_));
        Data.OnSet(5,23,HP_Sch(W_));
        Data.OnSet(5,24,HE_inf(m_));
        Data.OnSet(5,25,Htau(s_^-1));
        Data.OnSet(5,28,Htau_opt(s_^-1));
        Data.OnSet(5,26,HP_tot(W_));

        //Output vertical
        Data.OnSet(6,5,Vx(U1_));
        Data.OnSet(6,6,VG(U1_));
        Data.OnSet(6,7,VG_log(U1_));
        Data.OnSet(6,27,VG_opt(U1_));
        Data.OnSet(6,12,VPU_al(m_));
        Data.OnSet(6,19,VK_al(m_));
        Data.OnSet(6,14,0.5*VPU_S(U1_));
        Data.OnSet(6,21,0.5*VK_S(U1_));
        Data.OnSet(6,22,VP_th(W_));
        Data.OnSet(6,23,VP_Sch(W_));
        Data.OnSet(6,24,VE_inf(m_));
        Data.OnSet(6,25,Vtau(s_^-1));
        Data.OnSet(6,28,Vtau_opt(s_^-1));
        Data.OnSet(6,26,VP_tot(W_));
        //Output longitudinal
        Data.OnSet(7,6,LG(U1_));
        Data.OnSet(7,7,LG_log(U1_));
        Data.OnSet(7,27,LG_opt(U1_));
        Data.OnSet(7,12,LPU_al(m_));
        Data.OnSet(7,19,LK_al(m_));
        Data.OnSet(7,22,LP_th(W_));
        Data.OnSet(7,23,LP_Sch(W_));
        Data.OnSet(7,24,LE_inf(U1_));
        Data.OnSet(7,25,Ltau(s_^-1));
        Data.OnSet(7,28,Ltau_opt(s_^-1));
        Data.OnSet(7,26,LP_tot(W_));

        return 0;
}
//---------------------------------------------------------------------------

vectorU xStochastic::Rates(xTime&t, xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);


   doubleU xx(1,m_);
   if(H_use)
   {
   HE_inf = HA3*U_Pow(((HK_N*HK_beta)/(HPU_N*HPU_beta)),0.5)*(T_PU+T_PA)*
            HK_S*HPU_h*H_l*(1.0+ring.Energy.Beta)*HG/
            (HPU_S*HK_h*ring.Energy.Beta*ring.Energy.Momentum*U_c*4.0*HA2);

   Hx = U_e*U_e*HG*3.0*HA2*beam.Emit[3]*xx*U_Pow((HK_N*HK_beta*HPU_N*HPU_beta/(xx*xx)),0.5)*
       HPU_S*HK_S*Z* (1.0+ring.Energy.Beta)* H_l/(ring.Trev*ring.Trev*
       16.0*U_Abs(ring.Eta)* Msigma * U_Pow(beam.Emit[2], 0.5)* HK_h*HPU_h * ring.Energy.Beta*
       ring.Energy.Momentum*U_c*HBand);

   HG_opt = 1.5*HG/Hx;

   Htau_opt = HBand*4.0* Msigma * U_Abs(ring.Eta)* U_Pow(beam.Emit[2], 0.5)* HBand*ring.Trev/
              (3.0 * beam.Emit[3]*U_Ln(2.0));

   Htau = 16.0*U_Abs(ring.Eta)* Msigma * U_Pow(beam.Emit[2], 0.5)* HBand *HBand * ring.Trev*
         Hx * (1.0 - 2.0*Hx*U_Ln((Hx+HF_max/HBand)/(Hx+HF_min/HBand)) + Hx*Hx/((Hx+HF_max/HBand)*(Hx+HF_min/HBand)))/
         (3.0 * beam.Emit[3]);

   rates[0] = -1.0*Htau * (beam.Emit[0]-HE_inf)/beam.Emit[0];

   HP_Sch = HA1*beam.Emit[3] * HPU_N * HPU_beta * Z * HPU_S * HPU_S*
            U_e*U_e*ring.Energy.Z*ring.Energy.Z* beam.Emit[0]*HG*HG*HBand/
            (ring.Trev*HPU_h*HPU_h);
   HP_tot = HP_th * pow(10.0,(CombLoss(U1_)/10.0))*pow(10.0,(Losses(U1_)/10.0));
   HP_tot += HP_Sch*pow(10.0,(Losses(U1_)/10.0));
   }

   if(V_use)
   {
   VE_inf = VA3*U_Pow(((VK_N*VK_beta)/(VPU_N*VPU_beta)),0.5)*(T_PU+T_PA)*
            VK_S*VPU_h*V_l*(1.0+ring.Energy.Beta)*VG/
            (VPU_S*VK_h*ring.Energy.Beta*ring.Energy.Momentum*U_c*4.0*VA2);

   Vx = U_e*U_e*VG*3.0*VA2*beam.Emit[3]*xx*U_Pow((VK_N*VK_beta*VPU_N*VPU_beta/(xx*xx)),0.5)*
       VPU_S*VK_S*Z* (1.0+ring.Energy.Beta)* V_l/(ring.Trev*ring.Trev*
       16.0*U_Abs(ring.Eta)* Msigma * U_Pow(beam.Emit[2], 0.5)* VK_h*VPU_h * ring.Energy.Beta*
       ring.Energy.Momentum*U_c*VBand);

   VG_opt = 1.5*VG/Vx;

   Vtau_opt = VBand*4.0* Msigma* U_Abs(ring.Eta)* U_Pow(beam.Emit[2], 0.5)* VBand*ring.Trev/
              (3.0 * beam.Emit[3]*U_Ln(2.0));

   Vtau = 16.0*U_Abs(ring.Eta)* Msigma * U_Pow(beam.Emit[2], 0.5)* VBand *VBand * ring.Trev*
         Vx * (1.0 - 2.0*Vx*U_Ln((Vx+VF_max/VBand)/(Vx+VF_min/VBand)) + Vx*Vx/((Vx+VF_max/VBand)*(Vx+VF_min/VBand)))/
         (3.0 * beam.Emit[3]);

   rates[1] = -1.0*Vtau * (beam.Emit[1]-VE_inf)/beam.Emit[1];

   VP_Sch = VA1*beam.Emit[3] * VPU_N * VPU_beta * Z * VPU_S * VPU_S*
            U_e*U_e*ring.Energy.Z*ring.Energy.Z* beam.Emit[0]*VG*VG*VBand/
            (ring.Trev*VPU_h*VPU_h);
   VP_tot = VP_th * pow(10.0,(CombLoss(U1_)/10.0))*pow(10.0,(Losses(U1_)/10.0));
   VP_tot += VP_Sch*pow(10.0,(Losses(U1_)/10.0));
   }
//13.12.2007

   if(L_use)
   {
   if(RHIC)
   {
      eta = ring.Eta;
   }
   else{
   kappa = ring.Eta * ring.Energy.Gamma / (ring.Energy.Gamma + U_1);

   Ltau0 = 4.0*LA1*U_e*U_e*U_Pow((LK_N*LPU_N),0.5)*Z*LG*LBand*LF_0*U_Abs(kappa)/ring.Energy.Kinetic;
   /*
   A = U_e*U_e*LA1*(T_PU+T_PA)*Z*LK_N*kappa*kappa*LG*LG*ring.Trev*LBand*
       (LF_0*LF_0 + LBand*LBand/12.0)/
       (ring.Energy.Kinetic*ring.Energy.Kinetic*U_pi*U_pi);
   */
   A = 4.0*U_e*U_e*LA1*(T_PU+T_PA)*Z*LK_N*kappa*kappa*LG*LG*ring.Trev*LBand*
       (LF_0*LF_0 + LBand*LBand/12.0)/
       (ring.Energy.Kinetic*ring.Energy.Kinetic);

   Ltau = Ltau0 - 3.0*A;

   LG_opt = (LA1/LA4)*8.0*U_Pow(U_pi/(LK_N*LPU_N), 0.5)* ring.Energy.Beta*
            ring.Energy.Momentum*U_c * U_Pow(beam.Emit[2], 0.5)* ring.Trev/
            (3.0*U_e*U_e* beam.Emit[3]* Z);

   Ltau_opt = 2.0*LBand*16.0* U_Abs(ring.Eta)* U_Pow((U_pi*beam.Emit[2]), 0.5)* ring.Trev*LF_0*LA1*LA1/
              (3.0 * beam.Emit[3]*LA4);
   //Ltau_opt = 2.0*LBand*2.0* U_Abs(ring.Eta)* U_Pow((2.0*U_pi*beam.Emit[2]), 0.5)* ring.Trev*LF_0/
   //           (beam.Emit[3]*LA1);

   B = 2.0*LA4 * U_e*U_e*U_e*U_e*LK_N*LPU_N*Z*Z*U_Abs(kappa)*LG*LG*LF_0*LBand/
       (ring.Trev*ring.Energy.Kinetic);

   rates[2] = (-2.0*Ltau) + ((0.75*B*beam.Emit[3]*ring.Energy.Gamma)/
              (U_Pow(U_pi*beam.Emit[2], 0.5)*ring.Energy.Kinetic*(1.0+ring.Energy.Gamma)));

   LE_inf = 3.0*B*beam.Emit[3]*ring.Energy.Gamma/
            (8.0*U_Pow(U_pi, 0.5)*Ltau*ring.Energy.Kinetic*(1.0+ring.Energy.Gamma));

   LP_Sch = 4.0 * LA1 * beam.Emit[3] * LPU_N * Z *
            U_e*U_e* ring.Eta* ring.Eta * beam.Emit[2]*
            LG*LG*LBand * ring.Trev * (LF_0*LF_0 + LBand*LBand/12.0);
   LP_tot = LP_th * pow(10.0,(CombLoss(U1_)/10.0))*pow(10.0,(Losses(U1_)/10.0));
   LP_tot += LP_Sch*pow(10.0,(Losses(U1_)/10.0));
   }
   }

   return rates;
}
//---------------------------------------------------------------------------
void xStochastic::Kick (xTime&time, xBeam&beam, xRing&ring)
{
   vectorU rates;
   beam.Emit = beam.Emit_Def(Lattice);
   rates = Rates(time, beam, ring);
   double r, k = 1;
   if (beam.benum == BUNCHED) k = 2;

   for (int j = 0; j < beam.Number(); j++)
   {
      if(H_use)
      {
          //Drift calculation
          r = (-1.0*Htau* time.dt)(U1_);
          if(r > -1.)
            r += 1.;
          else
            r = exp(r);
          beam(j,1, beam(j,1) * r);
          //Diffusion calculation

          r = ((((2.0*HE_inf*Htau)*time.dt/Lattice.betax)^0.5))(U1_);
          beam.Add(j,1,r*xDistributor::Gaussian());
      }

      if(V_use)
      {
          //Drift calculation
          r = (-1.0*Vtau* time.dt)(U1_);
          if(r > -1.)
            r += 1.;
          else
            r = exp(r);
          beam(j,3, beam(j,3) * r);
          //Diffusion calculation
          r = ((((2.0*VE_inf*Vtau)*time.dt/Lattice.betay)^0.5))(U1_);
          beam.Add(j,3,r*xDistributor::Gaussian());
      }

      if(L_use)
      {
          //         06.07
          //13.12.2007
          if(RHIC)
          {

             vectorU forces(U1_, m_^-1);
             doubleU diffusion(m_^-1);
             eta = ring.Eta;
             forces = GetForces (time, ring.Energy, beam(j));

             //Friction calculation
             //r = (-1.0*Friction_c* Rgain*time.dt * (alpha*beam(j,5)+beta*beam(j,5)*beam(j,5)*beam(j,5)))(U1_);
             r =  (forces[5]*ring.Circ*time.dt/ring.Trev)(U1_);
             beam(j, 5, beam(j,5).Saturation(r,0.0));
             //beam(j, 5, beam(j,5).Saturation(forces[5]*ring.Circ,0.0));
             //Diffusion calculation
             diffusion = GetDiffusion(time, ring.Energy, beam, beam(j));

             r = (( diffusion*time.ds)^0.5)(U1_);
             //r = 0.0;
             beam.Add(j,5,r*xDistributor::Gaussian());

          }else
          {

             //Drift calculation
             r = (-1.0*Ltau0* time.dt*k)(U1_);
             if(r > -1.)
             r += 1.;
             else
             r = exp(r);
             beam(j,5, beam(j,5) * r);
             //Diffusion calculation
             r = ((B*0.5*(beam.Emit[2]^0.5)*beam.Emit[3]*ring.Energy.Gamma/
                 (ring.Energy.Kinetic*(1.0+ring.Energy.Gamma))*time.dt*k)^0.5)(U1_);
             beam.Add(j,5,r*xDistributor::Gaussian());
             r = ((3.0*A*beam.Emit[2]*time.dt*k)^0.5)(U1_);
             beam.Add(j,5,r*xDistributor::Gaussian());
          }
      }
   }
}
//---------------------------------------------------------------------------
void xStochastic::Average(doubleU f1, doubleU f2, doubleU l)
{
 int stepnumber = 401;
 int i;
 doubleU arg(U1_);
 doubleU Sin(U1_);
 doubleU step = (f2 - f1)/stepnumber;
 doubleU f = f1-step;
 doubleU Sum1(0,U1_);
 doubleU Sum2(0,U1_);
 doubleU Sum3(0,U1_);

 for(i = 0; i < stepnumber; i++)
     {
     f += step;
     arg = 2*U_pi*f*l/(iRing.Energy.Beta*U_c);
     Sin = U_Sin(arg);
     Sum1 += Sin*Sin;
     Sum2 += Sin*Sin/arg;
     //Sum3 += 1.0;
     Sum3 += Sin*Sin/(arg*arg);
     }
     HA1 = Sum1/stepnumber;
     HA2 = Sum2/stepnumber;
     HA3 = Sum3/stepnumber;
}
//---------------------------------------------------------------------------
void xStochastic::Average2(doubleU f1, doubleU f2, doubleU l)
{
 int stepnumber = 401;
 int i;
 doubleU arg(U1_);
 doubleU Sin(U1_);
 doubleU step = (f2 - f1)/stepnumber;
 doubleU f = f1-step;
 doubleU Sum1(0,U1_);
 doubleU Sum2(0,U1_);
 doubleU Sum3(0,U1_);

 for(i = 0; i < stepnumber; i++)
     {
     f += step;
     arg = 2*U_pi*f*l/(iRing.Energy.Beta*U_c);
     Sin = U_Sin(arg);
     Sum1 += Sin*Sin;
     Sum2 += Sin*Sin/arg;
     //Sum3 += 1.0;
     Sum3 += Sin*Sin/(arg*arg);
     }
     VA1 = Sum1/stepnumber;
     VA2 = Sum2/stepnumber;
     VA3 = Sum3/stepnumber;
}
//---------------------------------------------------------------------------
void xStochastic::Average3(doubleU f1, doubleU f2, doubleU l)
{
 int stepnumber = 401;
 int i;
 doubleU arg(U1_);
 doubleU Sin(U1_);
 doubleU step = (f2 - f1)/stepnumber;
 doubleU f = f1-step;
 doubleU Sum1(0,U1_);
 doubleU Sum2(0,U1_);
 doubleU Sum3(0,U1_);
 doubleU Sum4(0,U1_);
 
 for(i = 0; i < stepnumber; i++)
     {
     f += step;
     arg = 2*U_pi*f*l/(iRing.Energy.Beta*U_c);
     Sin = U_Sin(arg);
     Sum1 += Sin*Sin;
     Sum2 += Sin*Sin/arg;
     //Sum3 += 1.0;
     Sum3 += Sin*Sin/(arg*arg);
     Sum4 += Sin*Sin*Sin*Sin;
     }
     LA1 = Sum1/stepnumber;
     LA2 = Sum2/stepnumber;
     LA3 = Sum3/stepnumber;
     LA4 = Sum4/stepnumber;
}
//---------------------------------------------------------------------------

//12.12.2007
vectorU xStochastic::GetForces (xTime&t, U_Energy& e, vectorU ion)
{
   vectorU forces(U1_, m_^-1);
   doubleU domega(Hz_);
   int i;
   if(First)
   {
      T_0 = (*t.pCirc)/(*t.pVelocity);
      h_min = int((RLF_min*T_0)(U1_));
      h_max = int((RLF_max*T_0)(U1_));
      h_sub = (h_max - h_min)/N_sub;
      n_min = h_min/h_sub;
      n_max = h_max/h_sub;
      C_tilda = (e.Z*e.Z * U_rp*eta)/(T_0*e.A*e.Beta*e.Beta*e.Gamma);
      First = false;
   }

   domega = ion[5]* eta *2.0*U_pi/T_0;


      for(i = n_min; i <= n_max; i++)
      {
       if(n_series == 2)
       forces[5] += (2.0*U_Sin(domega*T_0*(U_1 + T_n)*h_sub*i)*
                    (U_Cos(domega*T_0*h_sub*i) - U_1))*(C_tilda*Rgain)*2*U_pi/(T_0*e.Velocity);
       if(n_series == 1)
       forces[5] += (U_Sin(domega*T_0* T_n*h_sub*i)-
                    U_Sin(domega*T_0*(U_1 + T_n)*h_sub*i))*(C_tilda*Rgain)*2*U_pi/(T_0*e.Velocity);
      }



    return forces;
}

doubleU xStochastic::GetDiffusion(xTime& t, U_Energy& e, xBeam& beam, vectorU ion)
{
doubleU diffusion(m_^-1);
   diffusion = 0.0;
   doubleU mult(U_1);
   doubleU domega(Hz_);
   int i;
   if(First)
   {
      T_0 = (*t.pCirc)/(*t.pVelocity);
      h_min = int((RLF_min*T_0)(U1_));
      h_max = int((RLF_max*T_0)(U1_));
      h_sub = (h_max - h_min)/N_sub;
      n_min = h_min/h_sub;
      n_max = h_max/h_sub;
      C_tilda = (e.Z*e.Z * U_rp*eta)/(T_0*e.A*e.Beta*e.Beta*e.Gamma);
      First = false;
   }

   domega = ion[5]* eta *2.0*U_pi/T_0;
   beam.RMS(beam.Emit, Lattice);

   mult = beam.Emit[3]/(U_Exp((ion[4]*ion[4])/(2.0*beam.Sigma_2[4]))*((2.0*U_pi)^0.5)*beam.Sigma[5]);

      for(i = n_min; i <= n_max; i++)
      {
       if(n_series == 2)
       diffusion += U_pi*4.0*(1.0 - U_Cos(domega*T_0*h_sub*i))*
                    (1.0 - U_Cos(domega*T_0*h_sub*i))* mult * (C_tilda*Rgain)*(C_tilda*Rgain)/(T_0*e.Velocity);
       if(n_series == 1)
       diffusion += U_pi*2.0*(1.0 - U_Cos(domega*T_0*h_sub*i))* mult* (C_tilda*Rgain)*(C_tilda * Rgain)/(T_0*e.Velocity);
      }

    diffusion *= Diffusion_c;

    return diffusion;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

xStochCool::xStochCool()
{
   use = false;
   type = 0;
   band._(G_9*Hz_);
   noise._(0);
   T_pk._(s_);
   T_kp._(s_);

}

doubleU xStochCool::Rates(xBeam&beam, xRing&ring, double P, double N, bool simple, bool optimal)
{
   doubleU rate(s_^-1);

   if (simple)
   {
      doubleU Mkp, Mpk2;
      Mkp  = 1. / (2.* band * eta_kp * T_kp * sqrt(P));
      Mpk2 =      (2.* band * eta_pk * T_pk * sqrt(P)) ^ 2.;
      if (optimal)
         gain = (1. - Mpk2) / (Mkp + noise);
      rate = - band / N * ((2.*gain*(1. - Mpk2)) - (gain*gain*(Mkp + noise)));
   }else
      rate = - 2. * band * band * U_Abs(ring.Eta) * ring.Trev * sqrt(P) / N;
   return rate;
}

//---------------------------------------------------------------------------
xGated::xGated(int N)
{
   ID = 9;
   Number = N;
   simple = true;
   optimal = false;
   System.SetNumber(N);
   F._(1/s_);
   D._(1/s_);
}

vectorU xGated::Rates(xTime& time, xBeam& beam, xRing& ring)
{
   vectorU rates(s_^-1);
   for(int i = 0; i < Number; i++)
   if (System[i].use)
      if (Numerical && i==0)
        {
         rates[2] = 2 * F * getForce((sqrt(beam.Emit[2]())/dpF)()) / sqrt(beam.Emit[2]())+
         D * getDifusion( (sqrt(beam.Emit[2]())/dpD)(), 1 ) / beam.Emit[2]();
        }else
        {
         rates[System[i].type] += System[i].Rates(beam, ring, beam.Emit[2](),
         beam.Emit[3](), simple, optimal);
        }
   return rates;
}

int xGated::OnGet()
{
   simple  = Data.OnGetB(91,41, true);
   optimal = Data.OnGetB(91,42, false);
   for (int i = 0; i < Number; i++)
   {
      System[i].use     = Data.OnGetB(91,1 + i*10, false);
      System[i].type    = Data.OnGetI(91,2 + i*10, 0);
      System[i].band    = Data.OnGet (91,3 + i*10, 1);
      System[i].initial = Data.OnGet (91,4 + i*10, 1);
      System[i].final   = Data.OnGet (91,5 + i*10, 2);
      System[i].T_pk    = Data.OnGet (91,6 + i*10, 0);
      if (simple && System[i].use && (System[i].T_pk > iRing.Trev))
         Warning("Time flight from pick-up to kicker more than revolution time", PressEnter);
      System[i].T_kp    = iRing.Trev - System[i].T_pk;
      System[i].eta_pk  = Data.OnGet (91,7 + i*10, 0);
      System[i].eta_kp  = Data.OnGet (91,8 + i*10, 0);
      System[i].gain    = Data.OnGet (91,9 + i*10, 1);
   }
      Numerical         = Data.OnGetB (92,1, 0);
      F                 = Data.OnGet (92,2, 0.0001);
      dpF               = Data.OnGet (92,3, -0.002);
      D                 = Data.OnGet (92,4, 1E-7);
      dpD               = Data.OnGet (92,5, -0.008);
      sigma             = Data.OnGet (92,6, 1e-3);
   return 0;
}

int xGated::OnSet()
{
   if (simple && optimal)
      for (int i = 0; i < Number; i++)
      if (System[i].use)
      {
             doubleU Mkp, Mpk2;
             Mkp  = 1. / (2.* System[i].band * System[i].eta_kp * System[i].T_kp * (iBeam.Emit[2]^0.5));
             Mpk2 =      (2.* System[i].band * System[i].eta_pk * System[i].T_pk * (iBeam.Emit[2]^0.5)) ^ 2.;
             System[i].gain = (1. - Mpk2) / (Mkp + System[i].noise);
             Data.OnSet(91,10+i*10, System[i].gain());
      }
   return 0;   
}

void xGated::Kick (xTime&time, xBeam&beam, xRing&ring)
{
   doubleU r;
   int n, p;

   for (int i = 0; i < Number; i++)
   if (System[i].use)
   {  n = System[i].type * 2 + 1;
      for (int j = 0; j < beam.Number(); j++)
      {
         if (System[i].initial * ring.Circ <= beam(j,4) &&
             System[i].final   * ring.Circ >= beam(j,4) )
         {
            if (Numerical && i==0)
            {
                r = F * getForce((beam(j,5)/dpF)()) * time.dt / beam(j,5) ;
            }else
            {
              if (beam.benum == BUCKET)
              {
                 p = ring.Bucket.Pos(beam, j);
                 if (ring.Bucket.barrier[p].P)
                 r = System[i].Rates(beam, ring, ring.Bucket.barrier[p].P,
                     ring.Bucket.barrier[p].D * beam.Emit[3](), simple, optimal) * time.dt;
                 else
                    r = 0;
              }else
                 r = System[i].Rates(beam, ring, beam.Emit[2](),
                     beam.Emit[3](), simple, optimal) * time.dt;
            }

            if(r > -1.)
               r += 1.;
            else
               r = U_Exp(r);

            beam(j, n, beam(j,n) * r);

            if (Numerical && i==0)
            {
               doubleU K = (beam.Emit[2]^0.5)/sigma;
               double k0 = getDifusion( (beam(j,5)/dpD)(), K() );
               double k1 = sqrt(((D *  k0 * time.dt)()))  * xDistributor::Gaussian() / beam(j,5)();
               double k2 = 1. + k1;
               beam(j, n, beam(j,n) * k2 );
            }
         }
      }
   }
}

double xGated::getForce (double dp)
{
      if (-1. <= dp && dp <= 1.)
            return - dp;
      else
      if (-1. > dp)
            return dp + 2.;
      else
            return dp - 2.;
}

double xGated::getDifusion(double dp, double sigma)
{
         if (-sigma  <= dp && dp <= 0)
            return -dp;
         else
         if (0 <= dp && dp <= sigma)
            return dp;
         else
         if (-2.* sigma <= dp && dp<= -sigma)
            return 2*sigma + dp;
         else
         if (sigma <= dp && dp <= 2. * sigma  )
            return -(dp - 2*sigma);
         else return 0;


}
