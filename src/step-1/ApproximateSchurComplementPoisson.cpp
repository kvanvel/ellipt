#include "ApproximateSchurComplementPoisson.hpp"

namespace heat{

ApproximateSchurComplementPoisson::ApproximateSchurComplementPoisson
(const dealii::SparseMatrixEZ<heat::real> & ConstraintMatrixIN,
 const dealii::SparseMatrix<heat::real> & StiffnessMatrixIN,
 const dealii::SparseMatrix<heat::real> & MassMatrixIN)
  :
  constraintMatrix(&ConstraintMatrixIN),
  stiffnessMatrix(&StiffnessMatrixIN),
  massMatrix(&MassMatrixIN),
  temp1(ConstraintMatrixIN.n()),
  temp2(ConstraintMatrixIN.n())
{
}


void 
ApproximateSchurComplementPoisson
::vmult(dealii::BlockVector<heat::real> & destination,
	const dealii::BlockVector<heat::real> & source) const
{
  stiffnessMatrix->vmult(temp1,source.block(0));
  constraintMatrix->Tvmult_add(temp1,source.block(1));
  massMatrix->precondition_Jacobi(temp2,temp1);
  stiffnessMatrix->Tvmult(destination.block(0), temp2);
  constraintMatrix->vmult(destination.block(1), temp2);
}

void
ApproximateSchurComplementPoisson
::Tvmult(dealii::BlockVector<heat::real> & destination,
	 const dealii::BlockVector<heat::real> &source) const
{
  vmult(destination, source);
}
  
} // end namespace heat
