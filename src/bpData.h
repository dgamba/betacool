#ifndef BPDATA // prevent variable redefinition    
#define BPDATA
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#ifndef M_PI
#define M_PI (3.141592653589793)
#endif
#define M_2PI (2.*M_PI)
extern const double NULLdata;

double Flatten();
double Gaussian();
class bpRing;

template <class L>
class Line {
protected:
	L*line;
public:
	int size;
   Line(int s = 1) { line = new L[size=s]; }
   virtual ~Line() { delete []line; }
   L& operator[](int i) { return line[i]; }
   void set(int s) { delete []line; line = new L[size=s]; }
   void operator=(Line& Line1){ set(Line1.size); for(int i=0;i<size;i++)line[i]=Line1.line[i]; }
};

template <class M>
class Mass {
protected:
   M**mass;
public:
   int row;
   int col;
   void NewMass(int r=1,int c=1){mass=new M*[row=r];for(int i=0;i<row;i++)mass[i]=new M[col=c];}
   void DelMass(){for(int i=0;i<row;i++)delete mass[i];delete[]mass;}
   Mass(int r=1,int c=1){NewMass(r,c);}
   virtual~Mass(){DelMass();}
   void set(int r=1,int c=1){DelMass();NewMass(r,c);}
   M*operator[](int i){return mass[i];}
   void operator=(Mass& Mass1){set(Mass1.row,Mass1.col);for(int i=0;i<row;i++)
      for(int j=0;j<col;j++)mass[i][j]=Mass1.mass[i][j];}
};

class bpData: public Mass<double> {
public:
   bpData(int r=1,int c=1) : Mass(r,c) {
    for (int i = 0; i < row; i++) 
      for (int j = 0; j < col; j++) 
         mass[i][j] = NULLdata;
   }
   int LoadFileSep(char* fn, char sep);
   void SaveFileSep(char* fn, char* sep, int RN = 0);
};

class bpRates {
public:
	double ratex;                                // horizontal rate, 1/s
	double ratey;                                // vertical rate, 1/s
	double ratez;                                // longitudinal rate, 1/s
   double raten;                                // particle number rate, 1/s
   bpRates() { ratex = 0; ratey = 0; ratez = 0; raten = 0; }
};

class bpEmittance {
public:
	double emittancex;                           // hirzontal emittance, m
	double emittancey;                           // vertical emittance, m
	double emittancez;                           // longitudinal emittance, m
	double momentum2;                            // momentum spread
   double sigmab;                               // beam length, m. If (bunched != 1) sigmabb = ring.Circ; 
   double sigma_s;                              // bunch size, m
}; 

class bpLattice {
public:
	double s;
	double ds;
	double betax;
	double alphax;
	double betay;	
	double alphay;
   double gammax;
   double gammay;
	double Dx;
	double Dpx;
	double Dy;
	double Dpy;
   bpLattice() { s = 0; ds = 1; gammax = 0; gammay =0;
	   betax = 1; alphax = 0; betay = 1; alphay = 0;
	   Dx = 0;    Dpx = 0;    Dy = 0;    Dpy = 0;
   }
};

class bpBeam : public Mass<double> {
public:
   bpBeam(int s=1) : Mass(6,s) { ; }
   void setsize(int s=1) { set(6,s); }
   void operator+=(bpBeam& beam1);
   void operator-=(int i);

  	double q0;                                   // charge
	double m0;                                   // mass
	double energy;                               // kinetic energy, MeV
   double beta;                                 // relativistic beta
   double gamma;                                // relativistic gamma
   double velocity;                             // absolute velocity, m
	double nion;                                 // real ion number
	int bunched;                                 // 0 - coasting, 1 - bunched
   bpEmittance emit;                            // beam emittance
   double Sigma[6];                             // RMS sigma values

   void Init();
   void MinusD(bpLattice& lat);
   void PlusD(bpLattice& lat);
   void RmsEmittance(bpLattice&);
   void Distribution(bpLattice& lat);
   void Invariants(bpLattice& lat, double B_s, int i, double Invs[3]);
   void TransRotation (bpLattice& L1, bpLattice& L2, double mu_x, double mu_y, int i);
   void LongRotation (bpRing& ring, double mu, double dt);
};

class bpRing : public Line<bpLattice> {
public:
   bpRing(int s=1) : Line(s) { ; }
   double Circ;                                 // Circumference, m
   double Imagenery;                            // -1 - imagenery, +1 - real
   double eta;                                  // of momemtum factor
   double GammaTr;                              // gamma trnsition
   double Q_h;                                  // horizontal tune
   double Q_v;                                  // vertical tune
   double Q_s;                                  // synchrotron tune
   double B_s;                                  // synchrotron function
   double L_s;                                  // separatrix length, m
   double Trev;                                 // revolution period, s
   double Hrf;                                  // RF harmonic number 
   double Vrf;                                  // RF voltage, kV
   bpLattice lattice;                           // lattice at injection point
   double Acceptancex;                          // horizontal acceptance, m
   double Acceptancey;                          // vertical acceptance, m

   void Init(bpBeam& beam, bool Vrf_OR_sigmab = false);
};

#endif
