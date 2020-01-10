//---------------------------------------------------------------------------
#ifndef bolideUH
#define bolideUH
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
// #include <io.h>

#if __BORLANDC__ == 1360 // C++Builder 5.0
#include <complex.h>
#else
#include <complex>
using namespace std;
#endif

//#include "stdafx.h"
#include "dataU.h"
#include "matrixU.h"
#include "warning.h"
//---------------------------------------------------------------------------
extern double* BParam[5];

class xData : public BTemplate<xData>
{public:
	//int Tag;
   static  char Beta[128];
	static  char File[128];
	static  BData Data;
   static  int Get(char* filename = NULL);
   virtual int OnGet() { return 0; }
   static  int Set(char* filename = NULL);
   virtual int OnSet() { return 0; }
   static  int Run(char* filename = NULL);
   virtual int OnRun() { return 0; }
};

class xGraf : public BTemplate<xGraf>
{  int Header;
   int SkipPoint;
   int SkipSave;
 public:
   BData Data;
   bool AutoSkip;
   bool Enabled;
   int Skip;
 	char FileName[256];
   xGraf(int size = 1000, char* filename = NULL);
   void SetName(char* filename);
   void Reset(int skip);
   void Size(int s, int c = 2) { Data.Size(s+1, c); Reset(0); }
   int  Size(){return Data.Row(); }
   void Point(double x = EmptyData, double y = EmptyData);
   void Points(double *X, int n);
   void Save(int period = 0);
   static void SaveAllGraf(int period = 0);
   static void ResetAll(int, bool);

   void Plot();
};

class xSurf : public BTemplate<xSurf>
{
   bool LogX;
   bool LogY;
	int SkipSave;
 public:
   bool Enabled;
 	char FileName[256];
   BData Data;
	xSurf(char* filename = NULL);
   void SetName(char* filename);
   void Reset(double reset = EmptyData);

   void SetXaxis(int size,double min=0,double max=10,bool log=false);
   void SetYaxis(int size,double min=0,double max=10,bool log=false);

   void   SetPoint(int x, int y, double z) { Data[x][y] = z; }
   void   AddPoint(int x, int y, double z) { Data[x][y]+= z; }
   double GetPoint(int x, int y)    { return Data[x][y]; }

   void   SetValue(double x, double y, double z);
   void   AddValue(double x, double y, double z);
   double GetValue(double x, double y);
	double*   Value(double x, double y);

   void Save(int period = 0);
   static void SaveAll(int period = 0);
};

void ShowTime(const char*,bool show = false);

class xLoader
{public:
   FILE *FilePointer;
   bool Check(char*,bool removefile = true);
};

class xTimer
{
	time_t time1;
	time_t time2;
   bool stop;

 public:
   static int year;
   static int month;
   static int day;
   static bool finish;
   xLoader Loader;
   long hours;
   long minutes;
   long seconds;
   xTimer(bool Stop = true);
   bool Timer(int interval);
   static void Converter();
};


#endif
