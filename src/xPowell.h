//---------------------------------------------------------------------------
#ifndef xPowellH
#define xPowellH
//---------------------------------------------------------------------------

void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret, double (*func)(double[]));
void powelltest();
#endif
