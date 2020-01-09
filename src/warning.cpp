/*
 * Warning function moved to this independent file
 *
 * Davide - Jan 2020
 */
#include "warning.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 

#include <cstddef>

using namespace std;

int WarningCount = 0;
char WarningFile[128] = "xbolide.war";
const int MaxWarning = 100;

void Warning(const char* str1, double ld1,
				 const char* str2, double ld2,
				 const char* str3, double ld3, bool endline)
{
	if (WarningCount > MaxWarning) return;
   WarningCount++;

   char Buffer[100], OutString[200];
   strcpy(OutString, str1);
   if ((ld1 != EmptyData) && (ld1 != PressEnter) && (ld1!= SpaceData))
   {	gcvt(ld1, 10, Buffer);
		if (Buffer[strlen(Buffer)-1] == '.') Buffer[strlen(Buffer)-1] = '\0';
		strcat(OutString, Buffer);
   }
   if (ld1 == SpaceData)
      strcat(OutString, " ");
   if (str2)
   	strcat(OutString, str2);
   if ((ld2 != EmptyData) && (ld2 != PressEnter) && (ld2!= SpaceData))
   {	gcvt(ld2, 10, Buffer);
		if (Buffer[strlen(Buffer)-1] == '.') Buffer[strlen(Buffer)-1] = '\0';
		strcat(OutString, Buffer);
   }
   if (ld2 == SpaceData)
      strcat(OutString, " ");
   if (str3)
   	strcat(OutString, str3);
   if ((ld3 != EmptyData) && (ld3 != PressEnter) && (ld3!= SpaceData))
   {	gcvt(ld3, 10, Buffer);
		if (Buffer[strlen(Buffer)-1] == '.') Buffer[strlen(Buffer)-1] = '\0';
		strcat(OutString, Buffer);
   }
   if (ld3 == SpaceData)
      strcat(OutString, " ");
   if (endline)
      strcat(OutString, "\n");
	FILE *FilePointer;
   FilePointer = fopen(WarningFile, "at");
	fputs(OutString, FilePointer);
	fclose(FilePointer);
   cout << OutString;
   if ((ld1 == PressEnter) || (ld2 == PressEnter) || (ld3 == PressEnter))
   {  cout << "Press Enter to continue...";
      scanf(OutString);
   }
}
