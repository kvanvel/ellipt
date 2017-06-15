#include "SS_UStar_massQ_U_SSStar.hpp"

namespace heat{

SS_UStar_massQ_U_SSStar::
SS_UStar_massQ_U_SSStar
(const dealii::SparseMatrixEZ<heat::real> &SS_In,
 const dealii::SparseMatrix<heat::real> &U_In,
 const dealii::SparseMatrix<heat::real> &massQ_In)
  :
  SS(&SS_In),
  U(&U_In),
  massQ(&massQ_In),
  temp1( SS_In.n() ),
  temp2( U_In.m() ),
  temp3( massQ_In.m() ),
  temp4( U_In.n() )
{}

unsigned int
SS_UStar_massQ_U_SSStar::m(void) const
{
  return SS->m();
}

void
SS_UStar_massQ_U_SSStar
::vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const

{  
  
  SS->Tvmult(temp1, source);  
  U->vmult(temp2,temp1);  
  massQ->vmult(temp3,temp2);
  U->Tvmult(temp4, temp3);
  SS->vmult(destination, temp4);
  
}

void
SS_UStar_massQ_U_SSStar
::Tvmult(dealii::Vector<heat::real> & destination,
	 const dealii::Vector<heat::real> & source) const
{
  vmult(destination,source);
}


} //end namespace heat
