//BOLIDE Self Pointer Template(s)

template <class B>
class DEFINEbTemplate
{
 public:
	static B **DEFINEbItems;
   static int DEFINEbCount;
   int DEFINEbIndex;
   DEFINEbTemplate();
   virtual ~DEFINEbTemplate();
};

template <class B>
B** DEFINEbTemplate<B>::DEFINEbItems = NULL;
template <class B>
int DEFINEbTemplate<B>::DEFINEbCount = 0;
                        
template <class B>
DEFINEbTemplate<B>::DEFINEbTemplate()
{
	B** newItems = new B*[DEFINEbCount+1];
   for (int i = 0; i < DEFINEbCount; i++)
		newItems[i] = DEFINEbItems[i];
	newItems[DEFINEbCount] = (B*)this;
   if (DEFINEbCount) delete []DEFINEbItems;
   DEFINEbItems = newItems;
   DEFINEbIndex = DEFINEbCount;
	DEFINEbCount++;
}

template <class B>
DEFINEbTemplate<B>::~DEFINEbTemplate()
{
	DEFINEbCount--;
   if (DEFINEbCount)
   {  B** newItems = new B*[DEFINEbCount];
      for (int i = 0; i < DEFINEbCount; i++)
         newItems[i] = DEFINEbItems[i];
      delete []DEFINEbItems;
      DEFINEbItems = newItems;
	}else
      delete []DEFINEbItems;
}

#undef DEFINEbTemplate
#undef DEFINEbItems
#undef DEFINEbCount
#undef DEFINEbIndex


