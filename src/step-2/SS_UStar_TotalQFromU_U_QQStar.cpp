#include "SS_UStar_TotalQFromU_U_QQStar.hpp"

namespace heat{

SS_UStar_TotalQFromU_U_QQStar::
SS_UStar_TotalQFromU_U_QQStar
(const dealii::SparseMatrixEZ<heat::real> &SS_In,
 const dealii::SparseMatrix<heat::real> &U_from_Neu_In,
 const dealii::SparseMatrix<heat::real> &TotalQFromU_In,
 const dealii::SparseMatrix<heat::real> &U_from_Dir_In,
 const dealii::SparseMatrixEZ<heat::real> &QQ_In)
  :
  SS(&SS_In),
  U_from_Neu(&U_from_Neu_In),
  TotalQFromU(&TotalQFromU_In),
  U_from_Dir(&U_from_Dir_In),
  QQ(&QQ_In),  
  temp1( QQ_In.n() ),
  temp2( U_from_Dir_In.m() ),
  temp3( TotalQFromU_In.m() ),  
  temp4( U_from_Neu_In.n() )
  //temp5( SS.m() )
{}

unsigned int
SS_UStar_TotalQFromU_U_QQStar::m(void) const
{
  return SS->m();
}

unsigned int
SS_UStar_TotalQFromU_U_QQStar::n(void) const
{
  return QQ->m();
}

void
SS_UStar_TotalQFromU_U_QQStar
::vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const

{
  QQ->Tvmult(temp1, source);
  U_from_Dir->vmult(temp2, temp1);
  TotalQFromU->vmult(temp3, temp2);
  U_from_Neu->Tvmult(temp4, temp3);
  SS->vmult(destination, temp4);  
}

void
SS_UStar_TotalQFromU_U_QQStar
::Tvmult(dealii::Vector<heat::real> & destination,
	 const dealii::Vector<heat::real> & source) const
{
  SS->Tvmult(temp4, source);
  U_from_Neu->vmult(temp3, temp4);
  TotalQFromU->Tvmult(temp2, temp3);
  U_from_Dir->Tvmult(temp1,temp2);
  QQ->vmult(destination, temp1);  
}


} //end namespace heat
