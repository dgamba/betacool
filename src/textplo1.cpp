//---------------------------------------------------------------------------
#include <iostream>
#include <stdio.h>
#include <dos.h>
#include <io.h>
#include "dataU.h"
//#pragma hdrstop
using namespace std;
//---------------------------------------------------------------------------

const int pSize = 24;
chars plot[pSize];
char* space80 = "                                                                          ";
char* caxis;
int gLeft  = 7;
int gWidth = 64;
int gTop   = 1;
int gHeight= pSize - 4;
int Xcol = 0;
int Ycol = 1;
int timer = 0;
int stoptime = 0;
int curtime = 0;
int counter = 60;
char TitleX[40] = "X_Axis";
char TitleY[40] = "Y_Axis";
char pname[80] = "";
bool remonce = true;
double Xmin, Xmax, Ymin, Ymax;
BData graf;
BData param;

void Warning(char* str1,        double ld1 = EmptyData,
				 char* str2 = NULL, double ld2 = EmptyData,
				 char* str3 = NULL, double ld3 = EmptyData,
             bool endline = true);
void Warning(char* str1, double ld1,
				 char* str2, double ld2,
				 char* str3, double ld3, bool endline)
{
   char Buffer[100], OutString[200];
   strcpy(OutString, str1);
   if ((ld1 != EmptyData) && (ld1 != PressEnter) && (ld1!= SpaceData))
   {	gcvt(ld1, 10, Buffer);
      strcat(OutString, Buffer);
   }
   if (ld1 == SpaceData)
      strcat(OutString, " ");
   if (str2)
   	strcat(OutString, str2);
   if ((ld2 != EmptyData) && (ld2 != PressEnter) && (ld2!= SpaceData))
   {	gcvt(ld2, 10, Buffer);
      strcat(OutString, Buffer);
   }
   if (ld2 == SpaceData)
      strcat(OutString, " ");
   if (str3)
   	strcat(OutString, str3);
   if ((ld3 != EmptyData) && (ld3 != PressEnter) && (ld3!= SpaceData))
   {	gcvt(ld3, 10, Buffer);
      strcat(OutString, Buffer);
   }
   if (ld3 == SpaceData)
      strcat(OutString, " ");
   if (endline)
      strcat(OutString, "\n");
   cout << OutString;
   if (strlen(pname))
   {
      FILE *FilePointer;
      FilePointer = fopen(pname, "at");
      fputs(OutString, FilePointer);
      fclose(FilePointer);
   }
   if ((ld1 == PressEnter) || (ld2 == PressEnter) || (ld3 == PressEnter))
   {  cout << "Press Enter to continue...";
      scanf(OutString);
   }
}

void axis(double&min, double&max, BData&data, int col)
{
   double mini = 0;
   double maxi = 0;
   for (int r = 0; r < graf.Row(); r++)
   if (data[r][col] != EmptyData)
   {
      if (mini > data[r][col]) mini = data[r][col];
      if (maxi < data[r][col]) maxi = data[r][col];
   }

   double maximum;
   if (fabs(maxi) > fabs(mini)) maximum = fabs(maxi);
      else                      maximum = fabs(mini);
   double maxlog;
   if (maximum > 0)
      maxlog = log10(maximum);
   else
      maxlog = 0.1;   
   int quotient  = maxlog;
   if (maxlog < 0) quotient -= 1;
   double remainder = maxlog - quotient;
   if      (remainder < log10(2)) maximum = 2;
   else if (remainder < log10(4)) maximum = 4;
   else if (remainder < log10(8)) maximum = 8;
   else                           maximum = 10;
   maximum *= pow(10, quotient);

   if (maxi > 0)
   {  max = maximum;
      if (mini >= 0) min = 0;
      else           min = -max;
   }else
   {  min = -maximum;
      max = 0;
   }
}

int main(int argc, char* argv[])
{
// Read param file
   if (argc == 1)
   {  cout << "Use: textplot filename1 [filename2] [...]\n";
      cout << "See example in textplot.txt file :\n";
      cout << "----------------------------------\n";
      strcpy(pname, "textplot.txt");
      remove(pname);
      Warning("space.cur; input file name");
      Warning("; output file name");
      Warning("5; saving interval [sec]");
      Warning("; finish time [sec]");
      Warning("luminosity; y-axis name");
      Warning("time [sec]; x-axis name");
      Warning("2; input file column of y-axis");
      Warning("0; input file column of x-axis");
      Warning("; y-axis max (nothing means auto)");
      Warning("; y-axis min");
      Warning("; x-axis max");
      Warning("; x-axis min");
      Warning("223; lower plot symbol");
      Warning("220; upper plot symbol");
      Warning("219; double plot symbol");
      return 1;
   }

while(true)
{
   for (int p = 1; p < argc; p++)
   {
      if (param.Load(argv[p]) == 0) return 1;
      if (graf.Load(param[0](0).ch.get()) == 0) return 1;
      strcpy(pname, param[1](0).ch.get());
      if(remonce) remove(pname);
      if (p == 1)
      {
         if (param[2][0] != EmptyData) timer = param[2][0];
         if (param[3][0] != EmptyData) stoptime = param[3][0];
      }
      strcpy(TitleY, param[4](0).ch.get());
      strcpy(TitleX, param[5](0).ch.get());
      Ycol = param[6][0];
      Xcol = param[7][0];

   // Axis gridlines and titles

      int i;
      for (i = 0; i < pSize; i++)
         plot[i] = space80;
      for (i = gLeft; i < gLeft+gWidth+1; i++)
      {  int d = gHeight / 4;
         plot[gTop]    [i] = '-';
         plot[gTop+d]  [i] = '-';
         plot[gTop+d*2][i] = '-';
         plot[gTop+d*3][i] = '-';
         plot[gTop+d*4][i] = '-';
      }
      for (i = gTop+1; i < gTop+gHeight; i++)
      {  int d = gWidth / 4;
         plot[i][gLeft]     = '|';
         plot[i][gLeft+d]   = '|';
         plot[i][gLeft+d*2] = '|';
         plot[i][gLeft+d*3] = '|';
         plot[i][gLeft+d*4] = '|';
      }

      for (unsigned j = 0; j < strlen(TitleY); j++)
         plot[0][gLeft+j+1] = TitleY[j];
      for (unsigned j = 0; j < strlen(TitleX); j++)
         plot[gTop+gHeight+2][j+gLeft+(gWidth-strlen(TitleX))/2] = TitleX[j];


// Axis min and max
   /*
      graf.Size(1000,2);
      for (int r = 0; r < graf.Row(); r++)
      {
         graf[r][0] = r*3;
         graf[r][1] = r-1000;
         //graf[r][0] = random(10000) - 2000;
         //graf[r][1] = random(400);
      }
   */

      axis(Ymin, Ymax, graf, Ycol);
      if (param[8][0] != EmptyData) Ymax = param[8][0];
      if (param[9][0] != EmptyData) Ymin = param[9][0];

      caxis = BData::double2char(Ymin,gLeft-2);
      int slen = strlen(caxis);
      for (unsigned int u = 0; u < strlen(caxis); u++)
         plot[gTop+gHeight-1][u+gLeft-slen] = caxis[u];

      caxis = BData::double2char((Ymax+Ymin)/2,gLeft-2);
      slen = strlen(caxis);
      for (unsigned int u = 0; u < strlen(caxis); u++)
         plot[gTop+gHeight/2][u+gLeft-slen] = caxis[u];

      caxis = BData::double2char(Ymax,gLeft-2);
      slen = strlen(caxis);
      for (unsigned int u = 0; u < strlen(caxis); u++)
         plot[gTop+1][u+gLeft-slen] = caxis[u];

      axis(Xmin, Xmax, graf, Xcol);
      if (param[10][0] != EmptyData) Xmax = param[10][0];
      if (param[11][0] != EmptyData) Xmin = param[11][0];

      caxis = BData::double2char(Xmin,gLeft-2);
      for (unsigned int u = 0; u < strlen(caxis); u++)
         plot[gTop+gHeight+1][u+gLeft+1] = caxis[u];

      caxis = BData::double2char((Xmax+Xmin)/2,gLeft-2);
      slen = strlen(caxis);
      for (unsigned int u = 0; u < strlen(caxis); u++)
         plot[gTop+gHeight+1][u+gLeft+1+(gWidth-slen)/2] = caxis[u];

      caxis = BData::double2char(Xmax,gLeft-2);
      slen = strlen(caxis);
      for (unsigned int u = 0; u < strlen(caxis); u++)
         plot[gTop+gHeight+1][u+gLeft+gWidth-slen] = caxis[u];

   // Curve

   /*
      for (i = 5; i < 34; i++)
      {
         plot[i-2][i*2+1] = 220;
         plot[i-2][i*2] = 223;
      }
   */
      int x,y;
      char up, lo, jo;
      if (param[12][0]!= EmptyData) up = char(param[12][0]);
      if (param[13][0]!= EmptyData) lo = char(param[13][0]);
      if (param[14][0]!= EmptyData) jo = char(param[14][0]);

      for (int r = 1; r < graf.Row(); r++)
      if ((graf[r][Xcol] != EmptyData) && (graf[r][Ycol] != EmptyData))
      {
         if (Ymin > graf[r][Ycol])  graf[r][Ycol] = Ymin;
         if (Ymax < graf[r][Ycol])  graf[r][Ycol] = Ymax;
         if (Xmin > graf[r][Xcol])  graf[r][Xcol] = Xmin;
         if (Xmax < graf[r][Xcol])  graf[r][Xcol] = Xmax;
      }


      for (int r = 0; r < graf.Row(); r++)
      if ((graf[r][Xcol] != EmptyData) && (graf[r][Ycol] != EmptyData))
      {
         x = Round(gWidth   *(graf[r][Xcol]-Xmin)/(Xmax-Xmin));
         y = Round(gHeight*2*(graf[r][Ycol]-Ymin)/(Ymax-Ymin));
         if (y % 2)
         {
            if(plot[gTop+gHeight-y/2][gLeft+x] == lo)
               plot[gTop+gHeight-y/2][gLeft+x] =  jo;
            if(plot[gTop+gHeight-y/2][gLeft+x] != jo)
               plot[gTop+gHeight-y/2][gLeft+x] =  up;
         }else
         {
            if(plot[gTop+gHeight-y/2][gLeft+x] == up)
               plot[gTop+gHeight-y/2][gLeft+x] =  jo;
            if(plot[gTop+gHeight-y/2][gLeft+x] != jo)
               plot[gTop+gHeight-y/2][gLeft+x] =  lo;
         }
      }
   // Output
      time_t t;
      struct tm *gmt;
      t = time(NULL);
      gmt = localtime(&t);
      Warning("       --------------------------------------------- ",
              gmt->tm_year+1900,"/", gmt->tm_mon+1,"/",gmt->tm_mday, false);
      Warning("-",gmt->tm_hour,":",gmt->tm_min,":",gmt->tm_sec);
      for (i = 0; i < pSize; i++)
         Warning(plot[i].get());
   }

   if (stoptime)
   {  curtime += timer;
      if (curtime > stoptime)
         return 1;
   }

   if (timer)
//      sleep(timer);

   {  //cout << "         0-------------1/4------------1/2------------3/4-------------1\n";
      cout << " wait ";
      cout << timer;
      cout << " sec : ";
      if (timer < counter)
      {  for (int i=0; i<timer; i++)
         {  sleep(1);
            for (int j=0; j < counter/timer; j++)
               cout << ".";
         }
      }else
      {  for (int i=0; i<counter; i++)
         {  sleep(timer/counter);
            cout << ".";
         }
      }
      cout << "\n";
   }

   else return 0;
   remonce = false;
}
}
//---------------------------------------------------------------------------
