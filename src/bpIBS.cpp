#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <time.h>
#include "bpIBS.h"
#include <stdlib.h>
//#define M_PI (3.141592653589793)
#define M_2PI (2.*M_PI)
using namespace std;

double q = 79;																										
double m0 = 197;																									
double delta = 4500;																								
double nb = pow(10, 9);																								
double gammaRel = 1 + (delta / 931.5);																					
double beta = sqrt(((pow(gammaRel, 2)) - 1) / pow(gammaRel, 2));															
double emittancex = pow(10, (-6));																					
double emittancey = pow(10, (-6));																					
double sigmab = 0.6;																								
double momentum = pow(10, (-3));																					
double c0 = 3 * pow(10, 8);																							
double r0 = 2.8179403267 * pow(10, (-15)) * q * q * 9.10938 * pow(10, -31) / (m0 * 1.66054 * pow(10, -27));			
double e0 = m0 * 1.66054 * pow(10, -27) * 3 * 3 * pow(10, 8) * pow(10, 8) / (1.60218 * pow(10, -19));				
double AIBS = (c0*nb*pow(r0, 2)) / (8.*M_PI*pow(beta, 3)*pow(gammaRel, 4)*emittancex*emittancey*sigmab*momentum);			
int bem = 1;
double maxy1 = 200;
double maxy2 = 1000;

void bpIBSinit(bpBeam& beamParam, bpEmittance& emit){
	q = beamParam.q0;																											//charge of ion
	m0 = beamParam.m0;																											//mass of ion
	delta = beamParam.energy;																									//energy, MeV
	nb = beamParam.nion;																											//number of particles
	gammaRel = 1 + (delta / 931.5);																								//relativistic factor
   beta = sqrt(((pow(gammaRel, 2)) - 1) / pow(gammaRel, 2));																//longditudial speed 
	emittancex = emit.emittancex;																								//emittance-x
	emittancey = emit.emittancey;																								//emittance-y
	sigmab = emit.sigmab;																										//size of bunch
	momentum = sqrt(emit.momentum2);																									//impulse scatter
	r0 = 2.8179403267 * pow(10, (-15)) * q * q * 9.10938 * pow(10, -31) / (m0 * 1.66054 * pow(10, -27));						//radius
	e0 = m0 * 1.66054 * pow(10, -27) * 3 * 3 * pow(10, 8) * pow(10, 8) / (1.60218 * pow(10, -19));								//charge
	AIBS = (c0*nb*pow(r0, 2)) / (8.*M_PI*pow(beta, 3)*pow(gammaRel, 4)*emittancex*emittancey*sigmab*momentum);	
	bem = beamParam.bunched;																									//bunched(1) or not(0)
   if (bem != 1) AIBS *= 2. * sqrt(M_PI);
}

void bpIntegrateAcc(double mmaxy1, double mmaxy2){
	maxy1 = mmaxy1;
	maxy2 = mmaxy2;
}

double Cubic(double a, double b, double c) {
	double q, r, r2, q3, x;
	q = (a*a - 3.*b) / 9.; r = (a*(2.*a*a - 9.*b) + 27.*c) / 54.;
	r2 = r*r; q3 = q*q*q;
	if (r2<q3) {
		double t = acos(r / sqrt(q3));
		a /= 3.; q = -2.*sqrt(q);
		return (q*cos((t + M_2PI) / 3.) - a);
	}
	else {
		double aa, bb;
		if (r <= 0.) r = -r;
		aa = -pow(r + sqrt(r2 - q3), 1. / 3.);
		if (aa != 0.) bb = q / aa;
		else bb = 0.;
		a /= 3.; q = aa + bb; r = aa - bb;
		x = (-0.5)*q - a;
		if (x == 0.) return(0);
		return(x);
	}
}

double coulombLog (bpLattice& Latt) {
	double v;
	if (bem == 1) {
		v = 8 * sqrt(pow(M_PI, 3) * (Latt.betax * emittancex + pow(Latt.Dx * momentum, 2)) * 
         (Latt.betay * emittancey + pow(Latt.Dy * momentum, 2))) * sigmab;
	}
	else {
		v = 4 * M_PI * sqrt((Latt.betax * emittancex + pow(Latt.Dx * momentum, 2)) * 
         (Latt.betay * emittancey + pow(Latt.Dy * momentum, 2))) * sigmab;
	}
	double rho = nb / (pow(10,6) * v);
	double tt = ((gammaRel * gammaRel) - 1) * e0 * emittancex / Latt.betax;
	double lDebye = 7.434 * sqrt (tt / rho) / q;
	double rmcl = 1.44 * q * q / (pow(10,9) * tt);
	double rmqm = 1.973 / (pow(10,13) * (2 * sqrt (tt * e0)));
	double rmax;
	double rmin;
	if (sqrt(Latt.betax * emittancex + pow((Latt.Dx * momentum), 2)) < lDebye) {
		rmax = sqrt(Latt.betax * emittancex + pow((Latt.Dx * momentum), 2));
	} else { rmax = lDebye; }
	if (rmcl > rmqm) {
	  rmin = rmcl;
	} else { rmin = rmqm; }
	return (log(rmax / rmin));
}


double integrate(double a1, double b1, double c1, int ex)
{
	double max = 0;
	double df1 = 0;
	double df2 = 0;
	double f = 0;
	double lyambda = 0;
	double a = (2 * ex*a1 - 5 * a1) / (2 * ex - 8);
	double b = (2 * b1*ex - 2 * b1) / (2 * ex - 8);
	double c = c1*(2 * ex + 1) / (2 * ex - 8);	
	max = Cubic(a, b, c);
	for (double m = 0; m < max * 10; m = m + max/maxy1) {
		df1 = pow(m, ex) * sqrt(m) * pow(10, 25) / pow((pow(m, 3) + a1 * pow(m, 2) + b1 * m + c1), 1.5);
		df2 = pow(m + max / maxy1, ex) * sqrt(m + max / maxy1) * pow(10, 25) / pow((pow(m + max / maxy1, 3) + a1 * pow(m + max / maxy1, 2) + b1 * (m + max / maxy1) + c1), 1.5);
		f = f + (((df1 + df2) / 2) * max / maxy1);
	}	
	for (double m = max * 10; m < max * maxy2; m = m + max) {
		df1 = pow(m, ex) * sqrt(m) * pow(10, 25) / pow((pow(m, 3) + a1 * pow(m, 2) + b1 * m + c1), 1.5);
		df2 = pow(m + max, ex) * sqrt(m + max) * pow(10, 25) / pow((pow(m + max, 3) + a1 * pow(m + max, 2) + b1 * (m + max) + c1), 1.5);
		f = f + (((df1 + df2) / 2) * max);
	}
	return f;
}

void bpIBSrates (bpRing& ring, bpRates& Growth, bpBeam& beam, bpEmittance& emittance) {

   if (emittance.emittancex < 1e-111 || emittance.emittancey < 1e-111 || emittance.momentum2 < 1e-222)
   {
	   Growth.ratex = 0;
	   Growth.ratey = 0;
	   Growth.ratez = 0;
      return;
   }

	bpIBSinit(beam, emittance);
	double ratex = 0;
	double ratey = 0;
	double ratez = 0; 
#pragma omp parallel for reduction (+: ratex) reduction(+: ratey) reduction(+: ratez)  
	for (int i = 0; i < ring.size; i++) {
		double phix, phiy, Hx, Hy, deltax, deltay, deltaz, a1, b1, c1, ax, bx, az, bz, ay, by, part;
		bpLattice Latt = ring[i];

      if ((abs(Latt.Dx) < 1e-50) && (abs(Latt.Dpx) < 1e-50)) { Latt.Dx = pow(10, -50); Latt.Dpx = pow(10, -50); }
		if ((abs(Latt.Dy) < 1e-50) && (abs(Latt.Dpy) < 1e-50)) { Latt.Dy = pow(10, -50); Latt.Dpy = pow(10, -50); }
		phix = (Latt.betax*Latt.Dpx + Latt.alphax * Latt.Dx) / Latt.betax;
		phiy = (Latt.betay*Latt.Dpy + Latt.alphay*Latt.Dy) / Latt.betay;
		Hx = (pow(Latt.Dx, 2) + (pow(Latt.betax, 2) * pow(phix, 2))) / Latt.betax;
		Hy = (pow(Latt.Dy, 2) + (pow(Latt.betay, 2) * pow(phiy, 2))) / Latt.betay;
		deltax = (pow(gammaRel, 2) * Hx) / emittancex;
		deltay = Latt.betay / emittancey;
		deltaz = pow(gammaRel, 2) / pow(momentum, 2);
		a1 = Latt.betax / emittancex + Latt.betay / emittancey + pow(gammaRel, 2)*(Hx / emittancex + Hy / emittancey + pow(momentum, -2));
		b1 = (Latt.betax*Latt.betay) / (emittancex*emittancey) + pow(gammaRel, 2)*((Latt.betax*Hx) / pow(emittancex, 2) + (Latt.betay*Hx) / (emittancex*emittancey) + (Latt.betay*Hy) / pow(emittancey, 2) + (Latt.betax*Hy) / (emittancex*emittancey) - (pow(Latt.betax, 2)*pow(phix, 2)) / pow(emittancex, 2) - (pow(Latt.betay, 2)*pow(phiy, 2)) / pow(emittancey, 2) + Latt.betax / (emittancex*pow(momentum, 2)) + Latt.betay / (emittancey*pow(momentum, 2)));
		c1 = pow(gammaRel, 2)*((Latt.betax*Latt.betay*Hx) / (pow(emittancex, 2)*emittancey) + (Latt.betax*Latt.betay*Hy) / (emittancex*pow(emittancey, 2)) - (pow(Latt.betax, 2)*Latt.betay*pow(phix, 2)) / (pow(emittancex, 2)*emittancey) - (Latt.betax*pow(Latt.betay, 2)*pow(phiy, 2)) / (emittancex*pow(emittancey, 2)) + (Latt.betax*Latt.betay) / (emittancex*emittancey*pow(momentum, 2)));
		ax = (-2 * Latt.betax) / emittancex - Latt.betay / emittancey + (2 * pow(Latt.betax, 2)) / (emittancex*pow(gammaRel, 2)*Hx) - (Latt.betax*Latt.betay) / (emittancey*pow(gammaRel, 2)*Hx) + (2 * pow(gammaRel, 2)*Hx) / emittancex + (2 * pow(gammaRel, 2)*Hy) / emittancey - (Latt.betax*Hy) / (emittancey*Hx) + (6 * pow(Latt.betax, 2)*pow(phix, 2)) / (emittancex*Hx) + (2 * pow(gammaRel, 2)) / pow(momentum, 2) - Latt.betax / (Hx*pow(momentum, 2));
		bx = pow(Latt.betax, 2) / pow(emittancex, 2) - (4 * Latt.betax*Latt.betay) / (emittancex*emittancey) + (pow(Latt.betax, 2)*Latt.betay) / (emittancex*emittancey*pow(gammaRel, 2)*Hx) + (Latt.betax*pow(gammaRel, 2)*Hx) / pow(emittancex, 2) + (Latt.betay*pow(gammaRel, 2)*Hx) / (emittancex*emittancey) + (Latt.betay*pow(gammaRel, 2)*Hy) / pow(emittancey, 2) + (Latt.betax*pow(gammaRel, 2)*Hy) / (emittancex*emittancey) - (2 * Latt.betax*Latt.betay*Hy) / (pow(emittancey, 2)*Hx) + (pow(Latt.betax, 2)*Hy) / (emittancex*emittancey*Hx) - (pow(Latt.betax, 2)*pow(gammaRel, 2)*pow(phix, 2)) / pow(emittancex, 2) - (pow(Latt.betax, 3)*pow(phix, 2)) / (pow(emittancex, 2)*Hx) + (6 * pow(Latt.betax, 2)*Latt.betay*pow(phix, 2)) / (emittancex*emittancey*Hx) - (pow(Latt.betay, 2)*pow(gammaRel, 2)*pow(phiy, 2)) / pow(emittancey, 2) + (2 * Latt.betax*pow(Latt.betay, 2)*pow(phiy, 2)) / (pow(emittancey, 2)*Hx) + (Latt.betax*pow(gammaRel, 2)) / (emittancex*pow(momentum, 2)) + (Latt.betay*pow(gammaRel, 2)) / (emittancey*pow(momentum, 2)) + pow(Latt.betax, 2) / (emittancex*Hx*pow(momentum, 2)) - (2 * Latt.betax*Latt.betay) / (emittancey*Hx*pow(momentum, 2));
		az = -(Latt.betax / emittancex) - Latt.betay / emittancey + (2 * pow(gammaRel, 2)*Hx) / emittancex + (2 * pow(gammaRel, 2)*Hy) / emittancey + (2 * pow(gammaRel, 2)) / pow(momentum, 2);
		bz = (-2 * Latt.betax*Latt.betay) / (emittancex*emittancey) + (Latt.betax*pow(gammaRel, 2)*Hx) / pow(emittancex, 2) + (Latt.betay*pow(gammaRel, 2)*Hx) / (emittancex*emittancey) + (Latt.betay*pow(gammaRel, 2)*Hy) / pow(emittancey, 2) + (Latt.betax*pow(gammaRel, 2)*Hy) / (emittancex*emittancey) - (pow(Latt.betax, 2)*pow(gammaRel, 2)*pow(phix, 2)) / pow(emittancex, 2) - (pow(Latt.betay, 2)*pow(gammaRel, 2)*pow(phiy, 2)) / pow(emittancey, 2) + (Latt.betax*pow(gammaRel, 2)) / (emittancex*pow(momentum, 2)) + (Latt.betay*pow(gammaRel, 2)) / (emittancey*pow(momentum, 2));
		ay = -(Latt.betax / emittancex) + (2 * Latt.betay) / emittancey - (pow(gammaRel, 2)*Hx) / emittancex - (Latt.betax*pow(gammaRel, 2)*Hy) / (Latt.betay*emittancex) - (2 * pow(gammaRel, 2)*Hy) / emittancey + (2 * pow(gammaRel, 4)*Hx*Hy) / (Latt.betay*emittancex) + (2 * pow(gammaRel, 4)*pow(Hy, 2)) / (Latt.betay*emittancey) + (6 * Latt.betay*pow(gammaRel, 2)*pow(phiy, 2)) / emittancey - pow(gammaRel, 2) / pow(momentum, 2) + (2 * pow(gammaRel, 4)*Hy) / (Latt.betay*pow(momentum, 2));
		by = (Latt.betax*Latt.betay) / (emittancex*emittancey) - (2 * Latt.betax*pow(gammaRel, 2)*Hx) / pow(emittancex, 2) + (Latt.betay*pow(gammaRel, 2)*Hx) / (emittancex*emittancey) + (Latt.betay*pow(gammaRel, 2)*Hy) / pow(emittancey, 2) - (4 * Latt.betax*pow(gammaRel, 2)*Hy) / (emittancex*emittancey) + (Latt.betax*pow(gammaRel, 4)*Hx*Hy) / (Latt.betay*pow(emittancex, 2)) + (pow(gammaRel, 4)*Hx*Hy) / (emittancex*emittancey) + (pow(gammaRel, 4)*pow(Hy, 2)) / pow(emittancey, 2) + (Latt.betax*pow(gammaRel, 4)*pow(Hy, 2)) / (Latt.betay*emittancex*emittancey) + (2 * pow(Latt.betax, 2)*pow(gammaRel, 2)*pow(phix, 2)) / pow(emittancex, 2) - (pow(Latt.betax, 2)*pow(gammaRel, 4)*Hy*pow(phix, 2)) / (Latt.betay*pow(emittancex, 2)) - (pow(Latt.betay, 2)*pow(gammaRel, 2)*pow(phiy, 2)) / pow(emittancey, 2) + (6 * Latt.betax*Latt.betay*pow(gammaRel, 2)*pow(phiy, 2)) / (emittancex*emittancey) - (Latt.betay*pow(gammaRel, 4)*Hy*pow(phiy, 2)) / pow(emittancey, 2) - (2 * Latt.betax*pow(gammaRel, 2)) / (emittancex*pow(momentum, 2)) + (Latt.betay*pow(gammaRel, 2)) / (emittancey*pow(momentum, 2)) + (Latt.betax*pow(gammaRel, 4)*Hy) / (Latt.betay*emittancex*pow(momentum, 2)) + (pow(gammaRel, 4)*Hy) / (emittancey*pow(momentum, 2));
		//part = Latt.ds / (ring[ring.size - 1].s);
      part = Latt.ds / ring.lattice.ds;
		double lg = coulombLog(Latt);
		double f1 = integrate(a1, b1, c1, 1);
		double f2 = integrate(a1, b1, c1, 0);
		ratex += deltax * (ax * f1 / pow(10, 25) + bx * f2 / pow(10, 25)) * lg * part;
		ratey += deltay * (ay * f1 / pow(10, 25) + by * f2 / pow(10, 25)) * lg * part;
		ratez += deltaz * (az * f1 / pow(10, 25) + bz * f2 / pow(10, 25)) * lg * part;
	}

	Growth.ratex = ratex * AIBS;
	Growth.ratey = ratey * AIBS;
	Growth.ratez = ratez * AIBS;
}

void bpIbsKick(double dt, bpRates& growth, bpBeam& beam, bpEmittance& emit, bpLattice& latt) {

   double kb = 1; 
   if (beam.bunched == 1) kb = 2;
   
   for (int i = 0; i < beam.row; i++)
   {
      if (growth.ratex < 1e-111)
            beam[1][i] *= exp(growth.ratex * dt);
      else  beam[1][i] += sqrt(emit.emittancex/latt.betax*growth.ratex*dt*2.)*Gaussian();

      if (growth.ratey < 1e-111)
            beam[3][i] *= exp(growth.ratey * dt);
      else  beam[3][i] += sqrt(emit.emittancey/latt.betay*growth.ratey*dt*2.)*Gaussian();

      if (growth.ratez < 1e-111)
            beam[5][i] *= exp(growth.ratez * dt);
      else  beam[5][i] += sqrt(emit.momentum2*growth.ratez*dt*kb)*Gaussian();
   }
}

void bpSlices (bpRing& ring, bpBeam& beam, int nSlices, double dt) {
		
		int* n = new int[nSlices];
		bpBeam* beamSlices = new bpBeam[nSlices];
		bpRates Growth2;
		for (int j=0; j<nSlices; j++)
			n[j] = 0;
		double k = -ring.Circ/2;
      double dk = ring.Circ/nSlices;
      int fillSlices = 0;
		for (int j=0; j<nSlices; j++) {
			for (int i = 0; i<beam.row; i++) {
				if (beam[4][i]>k && beam[4][i]<k+dk)
					n[j] ++;
			}
			k += dk;
		}
		for (int j = 0; j<nSlices; j++) { 
			beamSlices[j].bunched = 0;
			beamSlices[j].emit = beam.emit;
			beamSlices[j].energy = beam.energy;
			beamSlices[j].m0 = beam.m0;
			beamSlices[j].nion = beam.nion;
			beamSlices[j].q0 = beam.q0;
         beamSlices[j].Init();
			beamSlices[j].setsize(n[j]);
		}		
		
		k = -ring.Circ/2;			
		for (int j = 0; j<nSlices; j++) {
			int t=0;
			for (int i = 0; i<beam.row; i++) {
				if (beam[4][i]>k && beam[4][i]<k+dk) {
					for (int r = 0; r<6; r++)
						beamSlices[j][r][t]=beam[r][i];
					t++;
				}				
			}
			k += dk;
		}
		int particleS = 0;
		int particleB = 0;
      double dr;
		for (int i = 0; i < nSlices; i++) {
			beamSlices[i].RmsEmittance(ring.lattice);
			if (n[i] > 2) {	
			   bpIBSrates(ring, Growth2, beamSlices[i], beamSlices[i].emit);
            dr = 1. * n[i] / beam.row * nSlices;
            Growth2.ratex *= dr;
            Growth2.ratey *= dr;
            Growth2.ratez *= dr;
			   bpIbsKick(dt, Growth2, beamSlices[i], beamSlices[i].emit, ring.lattice);			
			} 
			while (particleS < n[i]) {
				for (int j = 0; j<6; j++) {
					beam[j][particleB] = beamSlices[i][j][particleS];
				}
				particleB++;
				particleS++;
			}
			particleS = 0;
		}
		delete [] n;
		delete [] beamSlices;
}
