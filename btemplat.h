//---------------------------------------------------------------------------
#ifndef BTemplateH
#define BTemplateH
//---------------------------------------------------------------------------
//BOLIDE Self Pointer Template(s)

#define DEFINEbTemplate ATemplate
#define DEFINEbItems AItems
#define DEFINEbCount ACount
#define DEFINEbIndex AIndex
#include "btempinc.h"

#define DEFINEbTemplate BTemplate
#define DEFINEbItems BItems
#define DEFINEbCount BCount
#define DEFINEbIndex BIndex
#include "btempinc.h"

#define DEFINEbTemplate CTemplate
#define DEFINEbItems CItems
#define DEFINEbCount CCount
#define DEFINEbIndex CIndex
#include "btempinc.h"

#define Max(a,b) (a > b ? a : b)
#define Min(a,b) (a < b ? a : b)
//---------------------------------------------------------------------------
//BOLIDE List Pointer Template

template <class L>
class LTemplate
{
 //protected:
	int LSize;
 public:
	L** LList;
   LTemplate(int size = 1);
   virtual ~LTemplate();
	L& operator[](int index);
	void operator= (LTemplate&);
   inline int LGetSize() { return LSize; }
   void LSetSize(int size);
   void LInsert (int index);
   void LRemove (int index);
   void LAdd() { LSetSize(LSize+1); }
   void LReset(L);
};

template <class L>
LTemplate<L>::LTemplate(int size)
{
	LList = new L*[LSize=size];
   for (int i = 0; i < LSize; i++)
  		LList[i] = new L;
}

template <class L>
LTemplate<L>::~LTemplate()
{	for (int i = 0; i < LSize; i++)
		delete LList[i];
	delete []LList;
}

template <class L>
L& LTemplate<L>::operator[](int index)
{
	if ((index >= 0) && (index < LSize))
		return *LList[index];
#if __BORLANDC__ == 1360
	Warning("LTemplate has not [index] : ", index);
#endif
	return **LList;
}

template <class L>
void LTemplate<L>:: operator= (LTemplate<L>& LT)
{
	LSetSize(LT.LGetSize());
   for (int i = 0; i < LGetSize(); i++)
   	*LList[i] = LT[i];
}

template <class L>
void LTemplate<L>::LSetSize(int size)
{
#if __BORLANDC__ == 1360
	if (size < 0)
   { 	Warning("LTemplate size should be positive");
   	return;
   }
#endif
   if (LSize != size)
   {  L** NewList = new L*[size];
      for (int i = 0; i < Max(LSize,size); i++)
      {  if (i < Min(LSize,size))
         	NewList[i] = LList[i];
         else
         	if (LSize < size)
         		NewList[i] = new L;
            else
               delete LList[i];
      }
      delete []LList;
      LList = NewList;
      LSize = size;
	}
}

template <class L>
void LTemplate<L>::LInsert(int index)
{  if (index < 0)
   {
#if __BORLANDC__ == 1360
      Warning("Can not Insert [index] : ", index);
#endif
      return;
   }
   int j = 0;
   L** NewList = new L*[LSize+1];
   for (int i = 0; i < LSize+1; i++)
   {  if (i == index)
      {  NewList[j] = new L;
      	j++;
      }
      if (i < LSize)
         NewList[j] = LList[i];
      j++;
   }
   delete []LList;
   LList = NewList;
   LSize++;
}

template <class L>
void LTemplate<L>::LRemove(int index)
{

   if ((index < 0) || (LSize == 0))
   {
#if __BORLANDC__ == 1360
      Warning("Can not Remove [index] :", index,"for LSize :", LSize);
#endif
      return;
   }

   int j = 0;
   L** NewList = new L*[LSize-1];
   for (int i = 0; i < LSize; i++)
   {
      if (i == index)
      {  j--;
      	delete LList[i];
      }else
      {	NewList[j] = LList[i];
      }
      j++;
   }
   delete []LList;
   LList = NewList;
   LSize--;
}

template <class L>
void LTemplate<L>::LReset(L l)
{
	for (int i = 0; i < LSize; i++)
   	*LList[i] = l;
}

//---------------------------------------------------------------------------
//BOLIDE Double Reference Template

template <class W>
class WTemplate
{	static bool bW;
	W* pW;
 protected:
   WTemplate();
 public:
   ~WTemplate();
   W& WRef();
};

template <class W>
bool WTemplate<W>::bW = false;

template <class W>
WTemplate<W>::WTemplate()
{  bW = !bW;
	if (bW) pW = new W;
	else    pW = NULL;
}

template <class W>
WTemplate<W>::~WTemplate()
{	if (pW) delete pW;
}

template <class W>
W& WTemplate<W>::WRef()
{	if (pW) return *pW;
   else	  return *((W*)this);
}

//---------------------------------------------------------------------------
//BOLIDE Static Pointer Array

template <class Z>
class ZTemplate
{
	static int iZ;
	static Z** pZ;
 protected:
   bool bZ;
   ZTemplate();
 public:
   static Z& newZ(bool forward = true);
   void deleteZ() { bZ = true; }
};

template <class Z>
int  ZTemplate<Z>::iZ = 8;
template <class Z>
Z**  ZTemplate<Z>::pZ = NULL;

template <class Z>
ZTemplate<Z>::ZTemplate()
{
   bZ = true;
	if(!pZ)
   {  pZ = new Z*[iZ];
   	for (int i = 0; i < iZ; i++)
      	pZ[i] = new Z;
   }
}

template <class Z>
Z& ZTemplate<Z>::newZ(bool forward)
{
	if (forward)
   {
      for (int i = 0; i < iZ; i++)
      {	if(pZ[i]->bZ)
         {	pZ[i]->bZ = false;
            return *pZ[i];
         }
      }
   }else
   {  for (int i = iZ-1; i >= 0; i--)
      {	if(pZ[i]->bZ)
         {	pZ[i]->bZ = false;
            return *pZ[i];
         }
      }
   }
   Warning("No more temprorary ZTemplate");
   return *pZ[0];
}

//---------------------------------------------------------------------------
#endif

