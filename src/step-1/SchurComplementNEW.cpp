#include "SchurComplementNEW.hpp"
#include "SS_UStar_TotalQFromU_U_QQStar.hpp"




namespace heat{

template<class MATRIX1, class MATRIX2, class MATRIX3>
SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3>::SchurComplementNEW
(const MATRIX1 &Ainv_In,
 const MATRIX2 &B_In,
 const MATRIX3 &C_In)
  :
  Ainv(&Ainv_In),
  B(&B_In),
  C(&C_In),
  temp1(B_In.m() ),
  temp2(C_In.m() )
  
{}

template<class MATRIX1, class MATRIX2, class MATRIX3>
void
SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3>
::vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const
{
  
  B->vmult(temp1,source);
  Ainv->vmult(temp2, temp1);
  C->Tvmult(destination, temp2);
    
}

template<class MATRIX1, class MATRIX2, class MATRIX3>
void
SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3>
::Tvmult(dealii::Vector<heat::real> & destination,
	 const dealii::Vector<heat::real> & source) const

{
  C->vmult(temp2, source);
  Ainv->vmult(temp1, temp2);
  B->Tvmult(destination, temp1);
}

template<class MATRIX1, class MATRIX2, class MATRIX3>
unsigned int
SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3>
::m() const
{
 return  C->n();

}

template<class MATRIX1, class MATRIX2, class MATRIX3>
unsigned int
SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3>
::n() const
{
  return B->n();
}





} //end namespace heat

//SchurComplementNEW
template class heat::SchurComplementNEW< class dealii::IterativeInverse<dealii::Vector<heat::real> >,
					 class heat::SS_UStar_TotalQFromU_U_QQStar,
					 class heat::SS_UStar_TotalQFromU_U_QQStar >;
