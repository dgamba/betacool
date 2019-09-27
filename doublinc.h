// Unit Variables
//-----------------------------------------------------------------
// Should be included in main file of project
// before any global variables

// Attachments
MicroMega a_18(ATTO),f_15(FEMTO),p_12(PICO),n_9(NANO),u_6(MICRO),m_3(MILI);
MicroMega k_3(KILO),M_6(MEGA),G_9(GIGA),T_12(TERA),P_15(PETA),E_18(EXA);

// Base u
Units_ s_ (SECOND);
Units_ Q_ (QULOUN);
Units_ m_ (METER);
Units_ U1_(BASESIZE);

// Some derivative u
Units_ cm_  (m_/100);
Units_ barn_((cm_^2)*1e-24);
Units_ H_   (cm_*1e9);
Units_ Hz_  (s_^-1);
Units_ min_ (s_*60);
Units_ hour_(min_*60);
Units_ day_ (hour_*24);
Units_ year_(day_*365);

// Electromagnetic, perveance
Units_ e_  (Q_*1.6021893e-19);
Units_ q_  (Q_/3e9);
Units_ F_  (cm_*9e11);
Units_ V_  (Q_/F_);
Units_ A_  (Q_/s_);
Units_ Ohm_(V_/A_);
Units_ uP_ (u_6*A_/(V_^1.5));
Units_ G_  (V_/cm_*300);
Units_ T_  (G_*1e4);

// Power, temperature, momentum, massa
Units_ eV_   (e_*V_);
Units_ K_    (eV_*8.61738573e-5);
Units_ eV_c2_(eV_*((s_/(cm_*light))^2));
Units_ eV_c_ (eV_*  s_/(cm_*light));
Units_ J_    (Q_*V_);
Units_ W_    (J_/s_);
Units_ kg_   (J_*(s_/m_^2));
Units_ g_    (kg_/1000);

// Massa of particles
Units_ amu_      (M_6*eV_c2_ * 931.5016);
Units_ neutron_  (M_6*eV_c2_ * 939.5731);
Units_ proton_   (M_6*eV_c2_ * 938.2796);
Units_ electron_ (M_6*eV_c2_ * 0.5110034);

// Force, pressure
Units_ N_   (J_/m_);
Units_ Pa_  (N_/(m_^2));
Units_ atm_ (Pa_*1.013e5);
Units_ Torr_(atm_/760);

//Base constants
doubleU U_0		  ( 0 );
doubleU U_1      ( 1 );
doubleU U_hbar   ( 0.658212202e-15,   eV_*s_ );
doubleU U_amu    ( 1,                 amu_ );
doubleU U_mn     ( 1,                 neutron_ );
doubleU U_mp     ( 1,                 proton_ );
doubleU U_me     ( 1,                 electron_ );
doubleU U_k      ( 8.61738573e-5,     eV_/K_ );
doubleU U_c      ( light,             cm_/s_ );
doubleU U_e      ( 1,                 e_ );
doubleU U_Grav   ( 6.6725985e-8,     (cm_^3)/(g_*(s_^2)) );
doubleU U_grav   ( 9.80665,           m_/(s_^2) );
doubleU U_NA     ( 6.022136736e23 );
doubleU U_pi     ( 3.141592653589793238 );
doubleU U_exp    ( 2.718281828459045235 );
doubleU U_eps0   ( 1./(4.*U_pi()*9e+9), F_/m_ );
doubleU U_mu0    (     4.*U_pi()*1e-7,  H_/m_ );

// Derivative constants
doubleU U_fine   ( (U_e^2) / (U_hbar * U_c) );
doubleU U_re     ( (U_e^2) / (U_me*(U_c^2)) );
doubleU U_rp     ( (U_e^2) / (U_mp*(U_c^2)) );
doubleU U_r1     ( (U_hbar^2) / (U_me*(U_e^2)) );
doubleU U_R      ( U_me * (U_e^4) / ( (U_hbar^3) * 2) );
doubleU U_Ecoup  ( U_me * (U_e^4) / ( (U_hbar^2) * 2) );
doubleU U_lambdaC( U_pi * U_hbar * 2 / (U_me * U_c) );
doubleU U_muBorn ( U_e  * U_hbar     / (U_me * U_c * 2) );

