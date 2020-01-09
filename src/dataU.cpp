//---------------------------------------------------------------------------
#include "stdafx.h"
#include "datau.h"     
//#include <stdio.h>
//#include <io.h>
//---------------------------------------------------------------------------
double EmptyData    = -1e301;
double PressEnter   = -2e301;
double OutRangeData = -3e301;
double SpaceData    = -4e301;

long Round(double d)
{  long l = (long) d;
	if ( fabs(d-l) >= 0.5-1e-9 )
	{  if (d > 0) l++;
		else       l--;
	}
	return l;
}

//---------------------------------------------------------------------------
// Class hierarchy for reading (writing) information from (to) file
char char0 = '\0';

chars::~chars()
{  if (pchar)
   {  delete []pchar;
      pchar = NULL;
   }
}

char& chars::operator[](unsigned int index)
{
   if (pchar && (index < strlen(pchar)))
      return pchar[index];
   else
      return char0;
}

void  chars::set(char* c)
{
	if (pchar) delete []pchar;
   if (c)
   {  pchar = new char[strlen(c) + 1];
      //strcpy(pchar,c);
      for (unsigned int j = 0; j <= strlen(c); j++)
         pchar[j] = c[j];
	}else
   	pchar = NULL;
}

bool chars::operator==(char* c)
{
   if ((!pchar) || (!c))
      return false;
   if (strlen(pchar) != strlen(c))
      return false;
   for (unsigned int i = 0; i < strlen(pchar); i++)
      if (pchar[i] != c[i])
         return false;
   return true;
}

void chars::operator += (char* c)
{
   if (!pchar)
      *this = c;
   else
      if (c)
      {  char* newchar = new char[strlen(pchar) + strlen(c) + 1];
         for (unsigned int i = 0; i <  strlen(pchar); i++)
            newchar[i] = pchar[i];
         for (unsigned int j = 0; j <= strlen(c); j++)
            newchar[j + strlen(pchar)] = c[j];
         delete []pchar;
         pchar = newchar;
      }
}

unsigned int chars::Length()
{  if (!pchar) return 0;
   return strlen(pchar);
}

void chars::Length(unsigned int i)
{
	if (pchar) delete []pchar;
   pchar = new char[i + 1];
   for (unsigned int j = 0; j <= i; j++)
      pchar[j] = char0;
}

int chars::Position(char* c)
{
   if ((!pchar) || (!c))
      return -1;
   if (strlen(pchar) < strlen(c))
      return -1;
   for (unsigned int i = 0; i < strlen(pchar)-strlen(c)+1; i++)
   {  bool equal = true;
      for (unsigned int j = 0; j < strlen(c); j++)
         if (pchar[i+j] != c[j])
            equal = false;
      if (equal)
         return i;
   }
   return -1;
}

chars chars::RemoveChar(char c)
{  bool more;
   do
   {  more = false;
      for (unsigned int i = 0; i < strlen(pchar); i++)
         if (pchar[i] == c)
         {  for (unsigned int j = i; j < strlen(pchar); j++)
               pchar[j] = pchar[j+1];
            more = true;
         }
   }while (more);
   return *this;
}

chars chars::LeadsChar(char c)
{
   while (pchar[0] == c)
      for(unsigned int j = 0; j < strlen(pchar); j++)
          pchar[j] = pchar[j+1];
   while (pchar[strlen(pchar)-1] == c)
          pchar[strlen(pchar)-1] =  '\0';
   return *this;
}

chars chars::UpCase()
{
   for (unsigned int i = 0; i < strlen(pchar); i++)
      if ((pchar[i] >= 'a') && (pchar[i] <= 'z'))
         pchar[i] += char('A' - 'a');
   return *this;
}

chars chars::LowCase()
{
   for (unsigned int i = 0; i < strlen(pchar); i++)
      if ((pchar[i] >= 'A') && (pchar[i] <= 'Z'))
         pchar[i] -= char('A' - 'a');
   return *this;
}

//---------------------------------------------------------------------------

BCell::BCell()
{  sr.set("\t");
   pd = new double;
   *pd = EmptyData;
   bd = true;              
   fl = false;
}

void BCell::operator=(BCell& bcell)
{
	if (bcell.bd)
   	*pd = *bcell.pd;
	else	set(bcell.pd);
   ch = bcell.ch;
   sr = bcell.sr;
}

bool BLine::index(int i)
{
	if ((i >= 0) && (i < LGetSize()))
   return true;
	Warning("Row has not [index] : ", i);
	return false;
}
/*
void BLine::rf(int i, double*p)
{
	if(LGetSize( ) <= i)
   	LSetSize(i+1);
	LList[i]->set(p);
}

*/
void BLine::operator=(BLine& bline)
{
	LSetSize(bline.LGetSize());
   for (int i = 0; i < LGetSize(); i++)
   	*LList[i]  = *bline.LList[i];
}
//---------------------------------------------------------------------------
BData::BData()
{
	separators = "\n\t;\\ ";
	SetShift(0);
	Size(1, 1);
   //onget_r = 0;
   //onget_c = 0;
   check = false;
}

BData::BData(int r, int c)
{
	separators = "\n\t;\\ ";
	SetShift(0);
	Col(r, c);
   //onget_r = 0;
   //onget_c = 0;
   check = false;
}

char* BData::double2char(double d, int digits)
{
	static char  d2c[32];
	if (d == EmptyData) return "";
	gcvt(d, digits, d2c);
   if(d2c[strlen(d2c)-1] == '.')
      d2c[strlen(d2c)-1] = '\0';
   return d2c;
}

double BData::char2double(char* c)
{
	if (!strlen(c))
   	return EmptyData;
	char *endptr;
	double d = strtod(c, &endptr);
   if (*endptr != NULL)
   	return EmptyData;
	return d;
}

double BData::char2double(chars c)
{
	//if (!strlen(c.get()))
	if (!c.Length())
   	return EmptyData;
   c.LeadsChar();
	char *endptr;
	double d = strtod(c.get(), &endptr);
   if (*endptr != NULL)
   	return EmptyData;
	return d;
}

void BData::Col(int r, int c, int r0)
{
	if(LGetSize() < r)
		LSetSize(r);
   for (int i = r0; i < r; i++)
  		LList[i]->Size(c);
}

void BData::Set(int r, int c)
{
	if(LGetSize() < r+1)
		LSetSize(r+1);
	if(LList[r]->LGetSize() < c+1)
		LList[r]->LSetSize(c+1);
}

void BData::Size(int r, int c)
{
	LSetSize(r);
   for (int i = 0; i < r; i++)
  		LList[i]->Size(c);
}

void BData::Empty(double d, int r, int r0)
{  if (!r) r = Row();
   for (int i = r0; i < r; i++)
   	for (int j = 0; j < Col(i); j++)
   		(*this)[i][j] = d;
}

#define LINELENGTH 4096 // maximum length of file row

int BData::Load(char* FileName, BData* pmad)
{
	FILE *FilePointer;
	FilePointer = fopen(FileName, "rt");
	if (!FilePointer)
   {
    	Warning("Cannot open file : ", EmptyData, FileName);
      return 0;
   }

   BData& data = *this;
   static char line[LINELENGTH+1];
   static char buff[LINELENGTH+1];
   static char seps[100];
   int  seplen = separators.Length();
   int  row = ROWSHIFT;
   int col, head, tail;

while (fgets(line, LINELENGTH, FilePointer))
{
	if(data.Row() <= row)
   	data.Row(10 * row);
   if(pmad)
   {
      BData& bmad  = *pmad;
      col = bmad.Col();
      data.Col(row+1, col, row);  //col+1 ??

      for (int i = 0; i < col; i++)
      {
         tail = (int) bmad[0][i];
         head = (int) bmad[1][i] + 1;
         if(head >= (int) strlen(line))
            head  =       strlen(line) - 2;

         if ((tail >= 0) && (tail < head))
         {
         	for (int j = tail; j < head; j++)
            	buff[j-tail] = line[j];
            buff[head-tail] = '\0';

            //while (buff[strlen(buff)-1] == ' ')
            //       buff[strlen(buff)-1] = '\0';
            data[row](i) = buff;
         }else
         {  if (tail == -1 && head == 0)
               data[row](i) = "0";
            else
               data[row](i) = "";
         }
         data[row].sr(i, "\t");
      }
      row++;
   }else
	{
      while(line[0] == ' ')
         for (unsigned int l = 0; l < strlen(line); l++)
            line[l] = line[l+1];
      col  = COLSHIFT;
      tail = 0;
   	head = -1;
   	do
      {  bool OnFindSeps = true;
      	do
			{ 	head++;
            for (int i = 0; i < seplen; i++)
            {  if ((separators[i] == '\\') && (i < seplen-1))
               {	if (head)
                  {  i++;
                     if((line[head-1] != separators[i]) &&
                        (line[head]   == separators[i]) &&
                        (line[head+1] == separators[i]))
                     {
                     	int j = 0;
                        for (; line[head+j]==separators[i]; j++)
                        	seps[j] = separators[i];
                        seps[j] = '\0';
                     	OnFindSeps = false;
								break;
                     }
						}
               }else
               {	if (separators[i] == line[head])
                  { 	seps[0] = line[head];
                  	if(seps[0] == '\n')
                     	seps[0] =  '\t';
                     seps[1] = '\0';
                     OnFindSeps = false;
         }  }  }  }
         while (OnFindSeps);

         for (int i = tail; i < head; i++)
         	buff[i-tail] = line[i];
         buff[head-tail] = '\0';

         if(data.Col(row) <= col)
      		data.Col(row+1, 2*col, row);

         if(buff[strlen(buff)-1] == '\r')
            buff[strlen(buff)-1] =  '\0';
         data[row](col) = buff;
         data[row].sr(col, seps);

        	head += strlen(seps) - 1;
         tail = head + 1;
         col++;
      }while ((line[head]!='\n')&&(head<(int)strlen(line)));
  		data.Col(row+1, col, row);
      row++;
   }
}
	fclose(FilePointer);
	if (row > 0)
   {  data.Row(row);
      for (int i = 0; i < Row(); i++)
         for (int j = 0; j < Col(i); j++)
            data[i][j] = char2double(data[i](j).ch);
	}
	return 1;
}

int BData::Save(char* FileName, int Digit)
{
	static char buf[LINELENGTH*5];
   strcpy(buf, "");
	FILE *FilePointer;
	FilePointer = fopen(FileName, "wt");
	if (!FilePointer)
	{	Warning("Cannot save to file : ", EmptyData, FileName);
      return 0;
   }
   BData& data = *this;
  int row = Row();
	for (int i = ROWSHIFT; i < row; i++)
	{  int col = Col(i);
   	for (int j = COLSHIFT; j < col; j++)
   	{ 	if (data[i][j] == EmptyData)
            strcat(buf, data(i,j));
         else
            strcat(buf, double2char(data[i][j],Digit));
			if (j < col-1)
            strcat(buf, data[i].sr(j));
         if (strlen(buf) > LINELENGTH*4)
         {  fputs(buf, FilePointer);
            strcpy(buf, "");
         }
		}
      strcat(buf, "\n");
	}
	fputs(buf, FilePointer);
	fclose(FilePointer);
   return 1;
}

bool BData::Translate(bool OnGet, char* loadname, char* savename)
{
   if (!OnGet && (loadname[0] == '*'))
   {  Save(savename);
      return true;
   }

   BData This;
   This.AddSeparators(",=;[]");
   if (!This.Load(loadname))
      return false;

   BData& That = *this;
   bool incrow;
   int col = That.COLSHIFT;
   int row = That.ROWSHIFT;

   for (int i = 0; i < This.Row(); i++)
   {  incrow = true;
      if (This[i](0).sr == "[")
      {  incrow = false;
         for (int j = 1; j < This.Col(i); j++)
         {  if ((This[i](j).ch == "row") && (This[i](j).sr == "="))
               row = (int)This[i][++j];
            if ((This[i](j).ch == "col") && (This[i](j).sr == "="))
               col = (int)This[i][++j];
            if ((This[i](j).ch == "file") && (This[i](j).sr == "="))
            {  j++;
               if (!Translate(OnGet, This(i,j), This(i,j)))
                  return false;
            }
            if (This[i](j).sr == "]")
               break;
         }
      }else
      {  for (int j = 0; j < This.Col(i); j++)
         {  if((This.Col(i) == 1) && (This[i](j).ch == ""))
               break;
            if (OnGet)
            {  That.Set(row, col);
               That[row][col]    = This[i][j];
               That[row](col).ch = This[i](j).ch;
               if (That[row](col).fl)
               {  Warning("Error : redefinition of Data [",
                  row,"] [",col,"] =", That[row][col]);
                  return false;
               }
               That[row](col).fl = true;
            }else
            {
               This[i][j]    = That[row][col];
               This[i](j).ch = That[row](col).ch;
               That[row](col).fl = false;
            }
            col++;
            if (This[i](j).sr == ";")
            {  break;
            }
            if (This[i](j).sr == "=")
            {  incrow = false;
               break;
            }
         }
      }
      if (incrow)
      {  row++;
         col = That.COLSHIFT;
      }
   }
   if (!OnGet)
   {  if (savename)
         This.Save(savename);
      else
         This.Save(loadname);
   }
   return true;
}

void BData::ResetFlag(bool flag)
{
	for (int i = 0; i < Row(); i++)
 	for (int j = 0; j < Col(i); j++)
      (*this)[i](j).fl = flag;
}
//---------------------------------------------------------------------------

char* BData::OnGetC(int r, int c, char* def)
{
	if((LGetSize()>r)&&(LList[r]->LGetSize()>c))
	return (*this)(r, c);
   if (check) Warning("Input parameter [",r,",",c,"] was initialised with default value");
   return def;
}

bool BData::OnGetB(int r, int c, bool def)
{
	if((LGetSize()>r)&&(LList[r]->LGetSize()>c)&&((*this)[r][c]!=EmptyData))
	return (bool)((*this)[r][c]);
   if (check) Warning("Input parameter [",r,",",c,"] was initialised with default value:", (int) def);
   return def;
}

int BData::OnGetI(int r, int c, int def)
{  
	if((LGetSize()>r)&&(LList[r]->LGetSize()>c)&&((*this)[r][c]!=EmptyData))
	return (int)(*this)[r][c];
   if (check) Warning("Input parameter [",r,",",c,"] was initialised with default value:", def);
   return def;
}

double BData::OnGet(int r, int c, double def)
{
	if((LGetSize()>r)&&(LList[r]->LGetSize()>c)&&((*this)[r][c]!=EmptyData))
	return (*this)[r][c];
   if (check) Warning("Input parameter [",r,",",c,"] was initialised with default value:", def);
   return def;
}

doubleU BData::OnGet(int r, int c, doubleU def)
{  
	if((LGetSize()>r)&&(LList[r]->LGetSize()>c)&&((*this)[r][c]!=EmptyData))
	return (doubleU)(*this)[r][c];
   if (check) Warning("Input parameter [",r,"][",c,"] was initialised with default value");
   return def;
}

void BData::OnSet(int r, int c, bool d)
{
	if(LGetSize() > r)
	   if(LList[r]->LGetSize() > c)
         if((*this)[r][c]!=EmptyData)
		   (*this)[r][c] = (double)d;
}

void BData::OnSet(int r, int c, int d)
{
	if(LGetSize() > r)
	   if(LList[r]->LGetSize() > c)
         if((*this)[r][c]!=EmptyData)
		   (*this)[r][c] = (double)d;
}

void BData::OnSet(int r, int c, double d)
{
	if(LGetSize() > r)
	   if(LList[r]->LGetSize() > c)
         if((*this)[r][c]!=EmptyData)
		   (*this)[r][c] = d;
}

void BData::OnSet(int r, int c, doubleU d)
{
	if(LGetSize() > r)
	   if(LList[r]->LGetSize() > c)
         if((*this)[r][c]!=EmptyData)
		   (*this)[r][c] = d();
}

