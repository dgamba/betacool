//---------------------------------------------------------------------------
#include "stdafx.h"
#if __BORLANDC__// == 1360 // C++Builder 5.0
#include <vcl.h>
#pragma hdrstop
#include <condefs.h>
#pragma hdrstop
//---------------------------------------------------------------------------
USEUNIT("doubleu.cpp");
USEUNIT("datau.cpp");
USEUNIT("xring.cpp");
USEUNIT("xbeam.cpp");
USEUNIT("xoptics.cpp");
USEUNIT("xeffect.cpp");
USEUNIT("xrunge.cpp");
USEUNIT("xdynamic.cpp");
USEUNIT("xlibrary.cpp");
USEUNIT("xdraw.cpp");
USEUNIT("xibs.cpp");
USEUNIT("xdistributor.cpp");
USEUNIT("xebeam.cpp");
USEUNIT("xecool.cpp");
USEUNIT("xforce.cpp");
USEUNIT("xtarget.cpp");
USEUNIT("xrestgas.cpp");
USEUNIT("xpowell.cpp");
USEUNIT("xstoch.cpp");
USEUNIT("xhiroshi.cpp");
USEUNIT("matrixu.cpp");
USEUNIT("vectoru.cpp");
USEUNIT("xbucket.cpp");
USEUNIT("pellets.cpp");
USEUNIT("bolideu.cpp");
//---------------------------------------------------------------------------
#pragma argsused
#endif

#define _CRT_SECURE_NO_WARNINGS
#include "xdynamic.h"
//#pragma link "bolideu.obj"
//---------------------------------------------------------------------------

// Global variables
xDrift     iDrift;
xDraw      iDraw;
xBeam1     iBeam(iRing);
xRing      iRing(iBeam);
xTime      iTime(iRing);
xTaskRates iTaskRates;
xDynamics  iDynamics;
CollidingBeam cBeam;
//xDynamics  linj_loss;

// ECOOL objects
xEcool iEcool;
xForce iForce;
xEbeam iEbeam(iRing);

// xEffects
xIBS       iIBS;
xLosses    iLosses;
xHeat      iHeat;
xTarget    iTarget(iBeam);
xRestGas   iRestGas;
xColl      iColl;
xStochastic iStochastic;
xGated     iGated (4);
xLaserCooling iLaserCooling;
xLastEffect iLastEffect;
xComponent iComponent;

bool linj_loss;

int main(int argc, char* argv[])
{
   xGraf test;

   strcpy(WarningFile, argv[0]);
   int len = strlen(WarningFile);
   if(WarningFile[len-4] == '.')
      WarningFile[len-4] = '\0';
   strcat(WarningFile, ".war");
   xTimer::Converter();
#ifdef _DEBUG
	Warning("BETACOOL----------debug ver.6.2 (omp)----------",xTimer::year,"/",xTimer::month,"/",xTimer::day);
#else
	Warning("BETACOOL----------release ver.6.2 (omp)----------",xTimer::year,"/",xTimer::month,"/",xTimer::day);
#endif

   switch(argc)
   {case 2: 
    case 3: xData::Data.check = true;
      strcpy(xData::Beta, argv[0]);
      strcpy(xData::File, argv[1]);
      Warning(argv[0],SpaceData,argv[1],SpaceData,argv[2]);
      break;
    default:
      Warning("Use:",SpaceData,argv[0],SpaceData,"inputfile /parameters",PressEnter);
      return 3;
   }
   if (xData::Get(xData::File))
   	return 1;
   if (xData::Set(xData::File))
   	return 2;

   for (int i = 2; i < argc; i++)
   {  if (argv[i][0] == '/')
      {  switch (argv[i][1])
         {case 'd':
            iDynamics.Dynamics();
            break;
          case 'm':
				iDynamics.ModelBeam();
         	break;
          case 't':
				iDynamics.Tracking();
         	break;
          case 'l':
				iDynamics.Lattice();
         	break;
          case 'r':
				iDynamics.Rates();
         	break;
          case 'b':
				iDynamics.BeamTest();
         	break;
          case 'f':
				iDynamics.FFTest();
         	break;
          case 'c':
				iDynamics.LuminosityCalculation();
         	break;
          case 'g':
				iForce.FF(iTime, iEbeam, iRing);
         	break;
          case 's':
				iDraw.SpaceChargeDraw(iBeam,iEbeam);
         	break;
          case '9':  
          case 'i':
            iDynamics.Distribution(iBeam.InitialEmit,iRing.LATTICE, 0, false);
            if (iBeam.initial == 2)
            {  iBeam.Emit = iBeam.Emit_Def(iRing.LATTICE);
               iBeam.InitialEmit = iBeam.Emit;
            }
            iDraw.evolution = true;
            if(iBeam.benum == BUNCHED)
               iBeam.CalcBunch();
            if (iBeam.benum == BUCKET)
            {  iRing.Bucket.CalcParticle(iBeam, iRing);
               if (iRing.Bucket.Analytic && (iRing.Bucket.Stationary || iRing.Bucket.Moving))
                  for (int j4 = 0; j4 < iBeam.Number(); j4++)
                  {  while (iBeam(j4,4) > iRing.Circ/(2*iRing.h()))
                        iBeam[j4][4] -= iRing.Circ(m_)/iRing.h();
                     while (iBeam(j4,4) <-iRing.Circ/(2*iRing.h()))
                        iBeam[j4][4] += iRing.Circ(m_)/iRing.h();
                  }
            }      
            iBeam.LongProfile();
            iDraw.DrawRealSpace (iBeam, iRing, iRing.LATTICE, 0);
            if (iBeam.benum == BUCKET && iRing.Bucket.Show3D)
               iRing.Bucket.LuminosityTest();
            if (argv[i][1] == '9')
               iDraw.TuneShift();
         	break;

          case '3':
            iDraw.Rates(iTaskRates, iTime, iBeam, iRing);
            break;
          case '1':
            iDraw.DrawLaserForce(iTime, iRing, iLaserCooling);
            break;
          case '5':
            iDraw.DrawStochForce(iTime, iRing, iBeam, iStochastic);
            break;
          case 'e':
            iDraw.DrawEBunch(iEbeam);
            break;
          case 'n':
            iDraw.DrawEDensity(iEbeam);
            break;
          case 'p':
            iDynamics.Distribution(iBeam.InitialEmit,iRing.LATTICE, 0, false);
            if (iBeam.initial == 2)
            {  iBeam.Emit = iBeam.Emit_Def(iRing.LATTICE);
               iBeam.InitialEmit = iBeam.Emit;
            }
            iDraw.evolution = true;
            if (iBeam.benum == BUCKET) iRing.Bucket.CalcParticle(iBeam, iRing);
            iDraw.DrawPellet(0);
            break;
          case 'o':
            iDraw.Oscilograph();
            break;
          case '2':
            iDraw.DrawGain();
            break;
          case '4':
            iDraw.DrawLStoch();
            break;
          case 'x':
            iRing.Bucket.LuminosityTest();
            break;     
          case 'z':
            test.Plot();
            break;
          case 'v':
            iComponent.Plot();
            break;
          case 'w':
            iRing.Bucket.Generate();
            break; 
          default:
          	Warning("Undefined parameter : ",EmptyData,argv[i]);
				break;
         }
      }
   }

   xData::Set(xData::File);
   ShowTime("FINISH-");

	return 0;
}

