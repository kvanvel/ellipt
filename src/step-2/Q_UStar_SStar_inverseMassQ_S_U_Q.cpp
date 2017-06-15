#include "Q_UStar_SStar_inverseMassQ_S_U_Q.hpp"

namespace heat{

Q_UStar_SStar_inverseMassQ_S_U_Q::
Q_UStar_SStar_inverseMassQ_S_U_Q
(const dealii::SparseMatrixEZ<heat::real> &Q_In,
 const dealii::SparseMatrix<heat::real> &U_In,
 const dealii::SparseMatrix<heat::real> &totalQFromU_In,
 const dealii::SparseMatrix<heat::real> &inverseMassQ_In,
 const dealii::SparseMatrix<heat::real> &totalUFromQ_In)
  :
  Q(&Q_In),
  U(&U_In),
  totalQFromU(&totalQFromU_In),
  inverseMassQ(&inverseMassQ_In),
  totalUFromQ(&totalUFromQ_In),
  temp1( U_In.m() ),
  temp2( U_In.m() ),
  temp3( totalQFromU_In.m() ),
  temp4( inverseMassQ_In.m() ),
  temp5( totalUFromQ_In.m() ),
  temp6( U_In.m() )  
{}

unsigned int
Q_UStar_SStar_inverseMassQ_S_U_Q::m(void) const
{
  return Q->m();
}

void
Q_UStar_SStar_inverseMassQ_S_U_Q
::vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const

{  
  
  Q->Tvmult(temp1, source);  
  U->vmult(temp2,temp1);  
  totalQFromU->vmult(temp3,temp2);  
  inverseMassQ->vmult(temp4,temp3);  
  totalQFromU->Tvmult(temp5,temp4);  
  U->Tvmult(temp6, temp5);  
  Q->vmult(destination, temp6);  
  
}

void
Q_UStar_SStar_inverseMassQ_S_U_Q
::Tvmult(dealii::Vector<heat::real> & destination,
	 const dealii::Vector<heat::real> & source) const
{
  vmult(destination,source);
}


} //end namespace heat
