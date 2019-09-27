//---------------------------------------------------------------------------
#ifndef xstochH     
#define xstochH
//---------------------------------------------------------------------------
class xStochastic : public xEffect
{
 public:

   //Horizontal degree of freedom
   int H_use;     // usage of horizontal cooling
   doubleU HF_min; //lower frequency
   doubleU HF_max; //upper frequency
   doubleU HF_0; //medium frequency
   doubleU HBand; //Horizontal bandwidth
   doubleU HG;    //Linear Electronic gain
   doubleU HG_opt; // Optimum linear gain
   doubleU HG_log;  //logarithmic gain
   doubleU H_l;     //electrode length
   //PU parameters
   doubleU HPU_N;  //Number of horizontal pickup loop numbers
   doubleU HPU_S;  //Horizontal pickup Sensitivity
   doubleU HPU_al;  //Horizontal pickup app. length
   doubleU HPU_beta;  //Horizontal pickup beta function

   doubleU HPU_w;     //electrode width
   doubleU HPU_h;     //gap height
   //Kicker parameters
   doubleU HK_N;  //Number of horizontal kicker loop numbers
   doubleU HK_S;  //Horizontal kicker Sensitivity
   doubleU HK_al;  //Horizontal kicker app. length
   doubleU HK_beta;  //Horizontal kicker beta function
   doubleU HK_w;     //electrode width
   doubleU HK_h;     //gap height
   //Output parameters
   doubleU HE_inf;   //emittance at infinit time
   doubleU Htau;     //cooling rate
   doubleU Htau_opt; //optimum cooling rate
   doubleU HP_th;    //Thermal noise power
   doubleU HP_Sch;   //Schottky power in the cooling band
   doubleU HP_tot;   //Total power consumption
   doubleU Hx;        //horizontal x parameter

   doubleU HA1;    //Averaged over horizontal bandwidth
   doubleU HA2;
   doubleU HA3;

   //Vertical degree of freedom
   int V_use;     // usage of vertical cooling
   doubleU VF_min; //lower frequency
   doubleU VF_max; //upper frequency
   doubleU VF_0; //medium frequency
   doubleU VBand; //Vertical bandwidth
   doubleU VG;    //Electronic gain
   doubleU VG_opt; //Optimum electronic gain
   doubleU VG_log; // logarithmic gain
   doubleU V_l;     //electrode length
   //PU parameters
   doubleU VPU_N;  //Number of vertical pickup loop numbers
   doubleU VPU_S;  //Vertical pickup Sensitivity
   doubleU VPU_al;  //Vertical pickup app. length
   doubleU VPU_beta;  //Vertical pickup beta function
   doubleU VPU_w;     //electrode width
   doubleU VPU_h;     //gap height
   //Kicker parameters
   doubleU VK_N;  //Number of vertical kicker loop numbers
   doubleU VK_S;  //Vertical kicker Sensitivity
   doubleU VK_al;  //Vertical kicker app. length
   doubleU VK_beta;  //Vertical kicker beta function
   doubleU VK_w;     //electrode width
   doubleU VK_h;     //gap height
   //Output parameters
   doubleU VE_inf;   //emittance at infinit time
   doubleU Vtau;     //cooling rate
   doubleU Vtau_opt; //optimum cooling rate
   doubleU VP_th;    //Thermal noise power
   doubleU VP_Sch;   //Schottky power in the cooling band
   doubleU VP_tot;   //Total consumption power
   doubleU Vx;        //horizontal x parameter

   doubleU VA1;    //Averaged over vertical bandwidth
   doubleU VA2;
   doubleU VA3;

   //Longitudinal degree of freedom
   int L_use;     // usage of longitudinal cooling
   doubleU LF_min; //lower frequency
   doubleU LF_max; //upper frequency
   doubleU LF_0; //medium frequency
   doubleU LBand; //Vertical bandwidth
   doubleU LG;    //Electronic gain
   doubleU LG_opt; //Optimum electronic gain
   doubleU LG_log; // logarithmic gain
   doubleU L_l;     //electrode length
   //PU parameters
   doubleU LPU_N;  //Number of longitudinal pickup loop numbers
   doubleU LPU_al;  //Longitudinal pickup app. length
   //Kicker parameters
   doubleU LK_N;  //Number of longitudinal kicker loop numbers
   doubleU LK_al;  //Longitudinal kicker app. length
   //Output parameters
   doubleU LE_inf;   //emittance at infinit time
   doubleU kappa;    //ring parameter
   doubleU Ltau0;    //singleparticle cooling rate
   doubleU Ltau_opt; //optimum cooling rate
   doubleU A;        //constant for cooling time calculation
   doubleU Ltau;     //cooling rate
   doubleU B;        //constant for rate calculation
   doubleU LP_th;    //Thermal noise power
   doubleU LP_Sch;   //Schottky power in the cooling band
   doubleU LP_tot;   //Total consumption power

   doubleU LA1;    //Averaged over longitudinal bandwidth
   doubleU LA2;
   doubleU LA3;
   doubleU LA4;

   //Loop geometry and other parameters
   doubleU T_PU;  //Pickup effective temperature
   doubleU T_PA;  //Preamplifier effective temperature
   doubleU Msigma;//Total width of momentum spread
   doubleU Z;     //characteristic inpedance
   doubleU CombLoss; //Losses in combiner
   doubleU Losses;  // Other losses

   //            06.07
   //Simplified model for RHIC longitudinal cooling
   bool RHIC;
   doubleU Friction_c;
   doubleU Diffusion_c;
   doubleU Rgain;
   doubleU C_tilda;
   doubleU alpha;
   doubleU beta;
   //12.12.2007
   doubleU RLF_min;  //minimum frequency
   doubleU RLF_max;  //maximum frequency
   bool First;
   doubleU T_0;    //revolution period
   int n_series;     //nuber of the filters in series
   int N_sub;        //number of the subbands
   doubleU T_n;      //time delay in T_0
   int h_min;        //minimum harmonics of the revolution frecuency
   int h_max;        //maximum harmonics of the revolution frecuency
   int h_sub;        //number of harmonics per subband
   int n_min;        //minimum number of the subband
   int n_max;        //maximum number of the subband
   doubleU eta;      //the ring off-momentum factor
   vectorU xStochastic::GetForces (xTime&, U_Energy&, vectorU);
   doubleU xStochastic::GetDiffusion(xTime&, U_Energy&, xBeam&, vectorU);
   doubleU Diff;
   //12.12.2007
   //            06.07
   xStochastic();
   int OnGet();
   int OnRun();
   int OnSet();
   vectorU Rates(xTime&, xBeam&, xRing&);    // calculates effect's rates
   void Kick (xTime&, xBeam&, xRing&);
   void Average(doubleU, doubleU, doubleU);
   void Average2(doubleU, doubleU, doubleU);
   void Average3(doubleU, doubleU, doubleU);
};
//---------------------------------------------------------------------------

class xStochCool
{public:

   bool use;
   int type;
   doubleU band;
   doubleU initial;
   doubleU final;
   doubleU noise;                    
   doubleU gain;
   doubleU eta_pk;
   doubleU eta_kp;
   doubleU T_pk;
   doubleU T_kp;

   xStochCool();
   doubleU Rates(xBeam&, xRing&, double P, double N, bool simple, bool optimal);
};

class xGated : public xEffect
{public:

   int Number;
   bool simple;
   bool optimal;
   bool Numerical;
   doubleU F;
   doubleU D;
   doubleU dpF;
   doubleU dpD;
   doubleU sigma;

   Tine<xStochCool> System;

   xGated(int N);
   int OnGet();
   int OnSet();
   int OnRun() { return OnGet(); }

   vectorU Rates(xTime&, xBeam&, xRing&);
   void Kick (xTime&, xBeam&, xRing&);
//----------------------------------------- Krestnikov. Stoch cooling numerical model 15.01.09
   double getForce (double);
   double getDifusion(double, double);
};

extern xGated iGated;
#endif
