#include "Q_UStar_massU_U_QStar.hpp"

namespace heat{

Q_UStar_massU_U_QStar::
Q_UStar_massU_U_QStar
(const dealii::SparseMatrixEZ<heat::real> &Q_In,
 const dealii::SparseMatrix<heat::real> &U_In,
 const dealii::SparseMatrix<heat::real> &massU_In)
  :
  Q(&Q_In),
  U(&U_In),  
  massU(&massU_In),
  temp1( Q_In.n() ),
  temp2( U_In.m() ),
  temp3( massU_In.m() ),
  temp4( U_In.n()  )
{}

unsigned int
Q_UStar_massU_U_QStar::m(void) const
{
  return Q->m();
}

void
Q_UStar_massU_U_QStar
::vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const

{  
  
  Q->Tvmult(temp1, source);  
  U->vmult(temp2,temp1);
  massU->vmult(temp3,temp2);
  U->Tvmult(temp4, temp3);  
  Q->vmult(destination, temp4);  
  
}

void
Q_UStar_massU_U_QStar
::Tvmult(dealii::Vector<heat::real> & destination,
	 const dealii::Vector<heat::real> & source) const
{
  vmult(destination,source);
}


} //end namespace heat
