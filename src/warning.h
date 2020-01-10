/*
 * Warning function moved to this independent file
 *
 * Davide - Jan 2020
 */
#ifndef xWARNING // prevent class redefenition     
#define xWARNING
#include <iostream>

extern double    EmptyData;
extern double    SpaceData;
extern int       WarningCount;
extern char      WarningFile[];
extern const int MaxWarning;
extern double    PressEnter;
void Warning(const char* str1, double ld1 = EmptyData,
			 const char* str2 = NULL, double ld2 = EmptyData,
			 const char* str3 = NULL, double ld3 = EmptyData,
                         bool endline = true);

#endif
