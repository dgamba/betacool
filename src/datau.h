//---------------------------------------------------------------------------
#ifndef dataUH
#define dataUH
#define _CRT_SECURE_NO_WARNINGS
#define M_PHI U1_ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "doubleU.h"
#include "btemplat.h"
#include "warning.h"

extern double EmptyData;
extern double PressEnter;
extern double OutRangeData;
extern double SpaceData;
long Round(double d);

//---------------------------------------------------------------------------
// array of chars
class chars
{  char* pchar;                                   // string
 public:

	chars()        { pchar = NULL; }               // default constructor
   chars(char*  c) { pchar = NULL; set(c);}       // constructor with string
   chars(chars& c) { pchar = NULL; set(c.pchar);} // constructor with chars
  ~chars();                                       // destructor
   char& operator[](unsigned int index);          // indexing of string
   char* get()    { if (pchar) return   pchar; return (char*)""; } // pointer on string
   void  set(char* c);                                      // set of string
   void  operator =  (char* c) { set(c); }                  // set of string
   void  operator =  (chars c) { *this = c.pchar; }         // set of string
   bool  operator == (char* c);                             // comparison of strings
   bool  operator != (char* c) { return !(*this == c); }    // comparison of strings
   bool  operator == (chars c) { return  (*this == c.pchar); } // comparison
   bool  operator != (chars c) { return !(*this == c.pchar); } // comparison
   void  operator += (char* c);                                // adding of strings
   void  operator += (chars c) { *this += c.pchar; }           // adding of strings
   unsigned int Length();                                      // get length
   void  Length(unsigned int i);                               // set length
   int   Position(char* c);                              // position of string
   int   Position(chars c) { return Position(c.pchar); } // position of string
   //bool  Empty() { return pchar; }                       // check empty
   chars RemoveChar(char c = ' ');                       // remove char
   chars LeadsChar (char c = ' ');                       // remove leads char
   chars UpCase();                                       // upcase of string
   chars LowCase();                                      // lowcase of string
};

// Unit cell for reading (writing) information from (to) file
class BCell
{	bool bd;
 public:
	double* pd;    // double value
	chars ch;      // chars value
   chars sr;      // separator
   bool fl;       // flag
	BCell();
  ~BCell()             { if (bd) delete pd; }
	void set(double* p) { if (bd) delete pd; pd = p; bd = false; }
	void operator=(BCell&);
   void operator=(double* p) { set(p); }
   void operator=(char* c) { ch = c; }
   void operator=(chars c) { ch = c; }
};

// Array of BCell
class BData;
class BLine : public LTemplate<BCell>
{  friend class BData;
	bool index(int);
   char* sr(int i)
   { if (index(i)) return LList[i]->sr.get( ); return (char*)"";}
   void  sr(int i, char* c)
   { if (index(i)) LList[i]->sr.set(c); }

 public:
   bool flag;    // sometimes useful
	void Size(int s) {  LSetSize(s); }
	int  Size(){ return LGetSize( ); }

// reloading operator[] which comes from LTemplate
	double& operator[](int i)
	{ if (index(i)) return *LList[i]->pd; return OutRangeData; }

	BCell& operator()(int i)
	{ if (index(i)) return *LList[i]; return *LList[0]; }

   void operator=(double*p) { LList[0]->set(p); }
   void operator=(BLine&);

};

// Array of BLine
class BData : public LTemplate<BLine>
{
   chars separators;
 public:
   static char*  double2char(double d, int digits = 10);
   static double char2double(char* c);
   static double char2double(chars c);
   void AddSeparators(char* c) { separators += c; }
   void SetSeparators(char* c)
   	   { separators[1] = '\0'; separators += c; }
	int  ROWSHIFT;
   int  COLSHIFT;
   void SetShift(int i = 1) { ROWSHIFT = i; COLSHIFT = i; }

 	BData();
 	BData(int r, int c);
   char* operator()(int r, int c) { return (*this)[r](c).ch.get(); }

   int  Row() { return LGetSize( ); }
   void Row(int r)   { LSetSize(r); }

   int  Col (int r = 0) { return LList[r]->Size( ); }
   void Col (int r, int c, int r0 = 0);
   void Set (int r, int c);
   void Size(int r, int c);

   void Empty(double d = EmptyData, int r = 0, int r0 = 0);

   int Load(char* FileName, BData* pmad = NULL);
	int Save(char* FileName, int Digit = 10);

   bool Translate(bool, char*, char* savename = NULL);
   void ResetFlag(bool);

   bool check;
   char*   OnGetC(int r, int c, char*   def);
   bool    OnGetB(int r, int c, bool    def);
   int     OnGetI(int r, int c, int     def);
   double  OnGet(int r, int c, double  def);
   doubleU OnGet(int r, int c, doubleU def);

   void OnSet(int r, int c, bool    d);
   void OnSet(int r, int c, int     d);
   void OnSet(int r, int c, double  d);
   void OnSet(int r, int c, doubleU d);

};

#endif
