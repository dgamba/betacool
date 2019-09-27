#ifndef BPIBS // prevent class redefenition     
#define BPIBS

#include "bpData.h"
void bpIbsKick(double dt, bpRates& growth, bpBeam& beam, bpEmittance& emit, bpLattice& latt);
void bpIBSrates (bpRing&, bpRates&, bpBeam&, bpEmittance&);
void bpIBSinit(bpBeam&, bpEmittance&);
void bpIntegrateAcc(double, double);
void bpSlices (bpRing&, bpBeam&, int, double);

#endif