// Icnlude file for tMatrix arithmetic

#ifdef PLUS_
template <class tM, int rM, int cM>
void tMatrix<tM, rM, cM>::operator SIGN_i(double m)
{
	for (int i = 0; i < Row; i++)
		(*this)[i][i] SIGN_i m;
}
#else
template <class tM, int rM, int cM>
void tMatrix<tM, rM, cM>::operator SIGN_i(double m)
{
	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			(*this)[i][j] SIGN_i m;
}
#endif
template <class tM, int rM, int cM>
tMatrix<tM, rM, cM> tMatrix<tM, rM, cM>::operator SIGN_(double m)
{
	tMatrix<tM, rM, cM> tmp(Row, Col);
	tmp = *this;
	tmp SIGN_i m;
	return tmp;
}
//---------------------------------------------------------------------------
#ifdef PLUS_
template <class tM, int rM, int cM>
void tMatrix<tM, rM, cM>::operator SIGN_i(tM m)
{
	for (int i = 0; i < Row; i++)
		(*this)[i][i] SIGN_i m;
}
#else
template <class tM, int rM, int cM>
void tMatrix<tM, rM, cM>::operator SIGN_i(tM m)
{
	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			(*this)[i][j] SIGN_i m;
}
#endif
template <class tM, int rM, int cM>
tMatrix<tM, rM, cM> tMatrix<tM, rM, cM>::operator SIGN_(tM m)
{
	tMatrix<tM, rM, cM> tmp(Row, Col);
	tmp = *this;
	tmp SIGN_i m;
	return tmp;
}
//---------------------------------------------------------------------------
#ifdef PLUS_
template <class tM, int rM, int cM>
void tMatrix<tM, rM, cM>::operator SIGN_i(tMatrix<tM, rM, cM> m)
{
	int row = Min(m.Row, Row);
	int col = Min(m.Col, Col);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			(*this)[i][j] SIGN_i m[i][j];
}
#endif
template <class tM, int rM, int cM>
tMatrix<tM, rM, cM> tMatrix<tM, rM, cM>::operator SIGN_(tMatrix<tM, rM, cM> m)
{
	int row = Min(m.Row, Row);
	int col = Min(m.Col, Col);
	tMatrix<tM, rM, cM> tmp(row, col);
	tmp = *this;
	tmp SIGN_i m;
	return tmp;
}
