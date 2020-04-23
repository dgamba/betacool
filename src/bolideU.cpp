//---------------------------------------------------------------------------
#include "stdafx.h"
#define _CRT_SECURE_NO_WARNINGS
#include "bolideU.h"
#include "xEffect.h"
#include "xDraw.h"
#include <iostream>
#include <stdio.h>
//#include <io.h>

//---------------------------------------------------------------------------
double InRangeData = 124e7;
double Interval1 = M_PI / 2. * 1e7;
double *BParam[5];
char *betacool_bld = "betacool.bld";

//---------------------------------------------------------------------------
char xData::Beta[];
char xData::File[];
BData xData::Data;

int xData::Get(char *FileName)
{
   if (FileName)
   {
      Data.AddSeparators(",");
      Data.SetShift();
      if (!Data.Translate(true, FileName))
         return 1;
   }
   for (int i = 0; i < BCount; i++)
      if (BItems[i]->OnGet())
         return 2;
   return 0;
}

int xData::Set(char *FileName)
{
   for (int i = 0; i < BCount; i++)
      if (BItems[i]->OnSet())
         return 2;
   if (FileName)
   {
      if (!Data.Translate(false, FileName))
         return 1;
   }
   return 0;
}

int xData::Run(char *FileName)
{
   if (FileName)
   {
      Data.ResetFlag(false);
      if (!Data.Translate(true, FileName))
         return 1;
   }
   for (int i = 0; i < BCount; i++)
      if (BItems[i]->OnRun())
         return 2;
   return 0;
}

//---------------------------------------------------------------------------

xGraf::xGraf(int s, char *filename)
{
   Data.Col(s, 2);
   SkipSave = 0;
   SkipPoint = 0;
   Skip = 0;
   AutoSkip = false;
   Enabled = true;
   Reset(0);
   if (filename == NULL)
   {
      strcpy(FileName, "xGraf");
      strcat(FileName, BData::double2char(BIndex));
      strcat(FileName, ".cur");
   }
   else
   {
      strcpy(FileName, filename);
   }
}

void xGraf::SetName(char *filename)
{
   strcpy(FileName, filename);
}

void xGraf::Reset(int skip)
{
   if (!Enabled)
      return;
   Header = 0;
   Skip = skip;
   SkipPoint = 0;
   Data.Empty(EmptyData);
}

void xGraf::Point(double x, double y)
{
   if (!Enabled)
      return;
   if (!SkipPoint || (SkipPoint > Skip))
   {
      SkipPoint = 0;
      Data[Header][0] = x;
      Data[Header][1] = y;
      Header++;
      if (Header > Data.Row() - 2)
      {
         if (AutoSkip)
         {
            for (int i = 0; i <= Header; i++)
            {
               if (i <= Header / 2)
               {
                  Data[i][0] = Data[i * 2][0];
                  Data[i][1] = Data[i * 2][1];
               }
               else
               {
                  Data[i][0] = EmptyData;
                  Data[i][1] = EmptyData;
               }
            }
            Header /= 2;
            if (Skip)
               Skip *= 2;
            else
               Skip = 1;
         }
         else
            Header = 0;
      }
      Data[Header][0] = EmptyData;
      Data[Header][1] = EmptyData;
   }
   SkipPoint++;
}

void xGraf::Points(double *X, int n)
{
   if (!Enabled)
      return;
   if (!SkipPoint || (SkipPoint > Skip))
   {
      SkipPoint = 0;
      for (int j = 0; j < n; j++)
         Data[Header][j] = X[j];
      Header++;
      if (Header > Data.Row() - 2)
      {
         if (AutoSkip)
         {
            for (int i = 0; i <= Header; i++)
            {
               if (i <= Header / 2)
               {
                  for (int j = 0; j < n; j++)
                     Data[i][j] = Data[i * 2][j];
               }
               else
               {
                  for (int j = 0; j < n; j++)
                     Data[i][j] = EmptyData;
               }
            }
            Header /= 2;
            if (Skip)
               Skip *= 2;
            else
               Skip = 1;
         }
         else
            Header = 0;
      }
      for (int j = 0; j < n; j++)
         Data[Header][j] = EmptyData;
   }
   SkipPoint++;
}

void xGraf::Save(int period)
{
   if (!Enabled)
      return;
   if (!SkipSave || (SkipSave >= period))
   {
      SkipSave = 0;
      Data.Save(FileName);
   }
   SkipSave++;
}

void xGraf::SaveAllGraf(int period)
{
   for (int i = 0; i < BCount; i++)
      if (!BItems[i]->Enabled)
      {
         //BItems[i]->Data.Empty(EmptyData);
         BItems[i]->Data.Save(BItems[i]->FileName);
      }
}

void xGraf::ResetAll(int csize, bool askip)
{
   for (int i = 0; i < BCount; i++)
   { //BItems[i]->Size(csize);
      BItems[i]->Data.Size(csize + 2, BItems[i]->Data.Col());
      BItems[i]->AutoSkip = askip;
   }
}

const int pSize = 40;
chars plot[pSize];
char *space80 = "                                                                             \n";
int gLeft = 5;
int gWidth = 71;
int gTop = 2;
int gHeight = 36;

void xGraf::Plot()
{
   int i;
   for (i = 0; i < pSize; i++)
      plot[i] = space80;
   for (i = gLeft; i < gLeft + gWidth; i++)
   {
      int d = (gHeight - 1) / 5;
      plot[gTop][i] = '-';
      plot[gTop + d][i] = '.';
      plot[gTop + d * 2][i] = '.';
      plot[gTop + d * 3][i] = '.';
      plot[gTop + d * 4][i] = '.';
      plot[gTop + d * 5][i] = '-';
   }
   for (i = gTop + 1; i < gTop + gHeight - 1; i++)
   {
      int d = (gWidth - 1) / 5;
      plot[i][gLeft] = '|';
      plot[i][gLeft + d] = ':';
      plot[i][gLeft + d * 2] = ':';
      plot[i][gLeft + d * 3] = ':';
      plot[i][gLeft + d * 4] = ':';
      plot[i][gLeft + d * 5] = '|';
   }
   for (i = 0; i < 5; i++)
   {
      plot[1][i] = '2';
      plot[38][i] = '3';
   }
   for (i = 5; i < 10; i++)
   {
      plot[39][i] = '4';
      plot[39][i + 66] = '5';
   }
   for (i = 10; i < 30; i++)
   {
      plot[i][10] = '@';
   }
   for (i = 0; i < 40; i++)
   {
      cout << plot[i].get();
   }
}

//---------------------------------------------------------------------------

xSurf::xSurf(char *filename)
{
   Enabled = true;
   SetXaxis(10);
   SetYaxis(10);
   SkipSave = 0;
   Reset();
   if (filename == NULL)
   {
      strcpy(FileName, "xSurf");
      strcat(FileName, BData::double2char(BIndex));
      strcat(FileName, ".sur");
   }
   else
   {
      strcpy(FileName, filename);
   }
}

void xSurf::SetName(char *filename)
{
   strcpy(FileName, filename);
}

void xSurf::Reset(double reset)
{
   for (int i = 1; i < Data.Row(); i++)
      for (int j = 1; j < Data.Col(); j++)
         Data[i][j] = reset;
}

void xSurf::SetXaxis(int size, double min, double max, bool log)
{
   if (!Enabled)
      return;
   Data.Size(size + 2, Data.Col());
   LogX = log;
   double step;
   if (max <= min)
   {
      Warning(FileName, EmptyData,
              " - Max of x-axis should be more than min");
      return;
   }
   if (log)
   {
      if ((max <= 0) || (min <= 0))
      {
         Warning(FileName, EmptyData,
                 " - Logarithmic x-axis should be positive");
         return;
      }
      step = log10(max / min) / size;
      for (int i = 1; i < size + 2; i++)
         Data[i][0] = pow(10, step * (i - 1)) * min;
   }
   else
   {
      step = (max - min) / size;
      for (int i = 1; i < size + 2; i++)
         Data[i][0] = step * (i - 1) + min;
   }
}

void xSurf::SetYaxis(int size, double min, double max, bool log)
{
   if (!Enabled)
      return;
   Data.Size(Data.Row(), size + 2);
   LogY = log;
   double step;
   if (max <= min)
   {
      Warning(FileName, EmptyData,
              " - Max of y-axis should be more than min");
      return;
   }
   if (log)
   {
      if ((max <= 0) || (min <= 0))
      {
         Warning(FileName, EmptyData,
                 " - Logarithmic y-axis should be positive");
         return;
      }
      step = log10(max / min) / size;
      for (int i = 1; i < size + 2; i++)
         Data[0][i] = pow(10, step * (i - 1)) * min;
   }
   else
   {
      step = (max - min) / size;
      for (int i = 1; i < size + 2; i++)
         Data[0][i] = step * (i - 1) + min;
   }
}

double *xSurf::Value(double x, double y)
{
   if (!Enabled)
      return &OutRangeData;
   if ((x >= Data[1][0]) && (y >= Data[0][1]) &&
       (x <= Data[Data.Row() - 1][0]) &&
       (y <= Data[0][Data.Col() - 1]))
   {
      int i, j;
      if (LogX)
         i = Round(log10(x / Data[1][0]) * (Data.Row() - 2) /
                   log10(Data[Data.Row() - 1][0] / Data[1][0]));
      else
         i = Round((x - Data[1][0]) * (Data.Row() - 2) /
                   (Data[Data.Row() - 1][0] - Data[1][0]));
      if (LogY)
         j = Round(log10(y / Data[0][1]) * (Data.Col() - 2) /
                   log10(Data[0][Data.Col() - 1] / Data[0][1]));
      else
         j = Round((y - Data[1][0]) * (Data.Col() - 2) /
                   (Data[0][Data.Row() - 1] - Data[0][1]));
      return &Data[i + 1][j + 1];
   }
   else
      return &OutRangeData;
}

void xSurf::SetValue(double x, double y, double z)
{
   if (!Enabled)
      return;
   double *pv = Value(x, y);
   if (pv == NULL)
      Warning(FileName, EmptyData,
              "- xSurf cannot SetValue : out of range x =", x, ", y =", y);
   else
      *pv = z;
}

void xSurf::AddValue(double x, double y, double z)
{
   if (!Enabled)
      return;
   double *pv = Value(x, y);
   if (pv == NULL)
      Warning(FileName, EmptyData,
              "- xSurf cannot AddValue : out of range x =", x, ", y =", y);
   else
      *pv += z;
}

double xSurf::GetValue(double x, double y)
{
   if (!Enabled)
      EmptyData;
   double *pv = Value(x, y);
   if (pv == NULL)
   {
      Warning(FileName, EmptyData,
              "- xSurf cannot GetValue : out of range x =", x, ", y =", y);
      return EmptyData;
   }
   return *pv;
}

void xSurf::Save(int period)
{
   if (!Enabled)
      return;
   if (!SkipSave || (SkipSave >= period))
   {
      SkipSave = 0;
      Data.Save(FileName);
   }
   SkipSave++;
}

void xSurf::SaveAll(int period)
{
   for (int i = 0; i < BCount; i++)
      BItems[i]->Save(period);
}

//---------------------------------------------------------------------------
#if defined(__BORLANDC__)
#include <dos.h>
#endif

void ShowTime(const char *ch, bool show)
{
   time_t t;
   struct tm *gmt;
   t = time(NULL);
   gmt = localtime(&t);
   Warning(ch, gmt->tm_year + 1900, "/", gmt->tm_mon + 1, "/", gmt->tm_mday, false);
   Warning("-", gmt->tm_hour, ":", gmt->tm_min, ":", gmt->tm_sec, !show);
   if (show)
      Warning(" : ", (t - InRangeData - Interval1) / Interval1);
}

int xTimer::year = 2009;
int xTimer::month = 2;
int xTimer::day = 17;
bool xTimer::finish = false;

xTimer::xTimer(bool Stop)
{
   time1 = 0;
   hours = 0;
   minutes = 0;
   seconds = 0;
   stop = Stop;
}

bool xTimer::Timer(int interval)
{
   time2 = time(NULL);
   if (time1 <= time2 - interval)
   {
      WarningCount = 0;
      seconds += long(time2 - time1);
      if (seconds > 1e9)
         seconds = 0;
      minutes = seconds / 60;
      hours = minutes / 60;
      time1 = time2;
      time2 -= long(InRangeData + Interval1);
      if ((time1 < InRangeData || time2 > 0) ||
          (stop && Loader.Check("bolide.stp")))
         finish = true;
      if (Loader.Check("bolide.pau", false))
      {
         Warning("BETACOOL is staying on PAUSE");
         while (Loader.Check("bolide.pau", false))
         {
            ;
         }
      }
      return true;
   }
   return false;
}

void xTimer::Converter()
{
#if defined(__DATE__)
   chars dd;
   dd = "00";
   chars ss;
   ss = __DATE__;
   ss.UpCase();
   dd[0] = ss[4];
   dd[1] = ss[5];
   dd.RemoveChar();
   day = (int)BData::char2double(dd);

   if (ss.Position("JAN") >= 0)
      month = 1;
   else if (ss.Position("FEB") >= 0)
      month = 2;
   else if (ss.Position("MAR") >= 0)
      month = 3;
   else if (ss.Position("APR") >= 0)
      month = 4;
   else if (ss.Position("MAY") >= 0)
      month = 5;
   else if (ss.Position("JUN") >= 0)
      month = 6;
   else if (ss.Position("JUL") >= 0)
      month = 7;
   else if (ss.Position("AUG") >= 0)
      month = 8;
   else if (ss.Position("SEP") >= 0)
      month = 9;
   else if (ss.Position("OCT") >= 0)
      month = 10;
   else if (ss.Position("NOV") >= 0)
      month = 11;
   else if (ss.Position("DEC") >= 0)
      month = 12;

   dd = "2000";
   dd[2] = ss[9];
   dd[3] = ss[10];
   year = (int)BData::char2double(dd);
#endif
   InRangeData = (day + (month - 1) * 30.5 + (year - 1970) * 365) * 86400;
}

bool xLoader::Check(char *filename, bool removefile)
{
   FilePointer = fopen(filename, "rt");
   if (FilePointer)
   {
      fclose(FilePointer);
      if (removefile)
         remove(filename);
      return true;
   }
   return false;
}
