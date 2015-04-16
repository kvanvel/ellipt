#include "ApproximateSchurComplement.hpp"

namespace heat{

ApproximateSchurComplement::ApproximateSchurComplement
(const dealii::SparseMatrix<heat::real> &AIN,
 const dealii::SparseMatrixEZ<heat::real> &ConstraintMatrixIN)
  :
  A( &AIN ),
  ConstraintMatrix( &ConstraintMatrixIN ),
  temp1( AIN.m() ),  
  temp2( AIN.m() )
{}

void 
ApproximateSchurComplement
::vmult(dealii::Vector<heat::real> & destination,
	const dealii::Vector<heat::real> & source) const
{
  ConstraintMatrix->Tvmult(temp1,source);
  A->precondition_Jacobi(temp2,temp1);
  ConstraintMatrix->vmult(destination,temp2);
}

void
ApproximateSchurComplement
::Tvmult(dealii::Vector<heat::real> & destination,
	 const dealii::Vector<heat::real> &source) const
{
  vmult(destination, source);
}
  
} // end namespace heat
