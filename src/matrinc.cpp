// Icnlude file for vectorU arithmetic

void vectorU::operator SIGN_i(double m)
{
	for (int i = 0; i < SIZE_; i++)
		ARRAY[i] SIGN_i m;
}
vectorU vectorU::operator SIGN_(double m)
{  vectorU tmp;
   tmp = *this;
   tmp SIGN_i m;
	return tmp;
}
//---------------------------------------------------------------------------
void vectorU::operator SIGN_i(doubleU m)
{
	for (int i = 0; i < SIZE_; i++)
		ARRAY[i] SIGN_i m;
}
vectorU vectorU::operator SIGN_(doubleU m)
{  vectorU tmp;
   tmp = *this;
   tmp SIGN_i m;
	return tmp;
}
//---------------------------------------------------------------------------
void vectorU::operator SIGN_i(vectorU m)
{
	for (int i = 0; i < SIZE_; i++)
		ARRAY[i] SIGN_i m.ARRAY[i];
}
vectorU vectorU::operator SIGN_(vectorU m)
{  vectorU tmp;
   tmp = *this;
   tmp SIGN_i m;
	return tmp;
}

