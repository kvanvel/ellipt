#include "SchurComplementPoisson.hpp"

namespace heat{

SchurComplementPoisson::SchurComplementPoisson
(const dealii::SparseMatrixEZ<heat::real> & ConstraintMatrixIN,
 const dealii::SparseMatrix<heat::real> & StiffnessMatrixIN,
 const dealii::IterativeInverse<dealii::Vector<heat::real> > & MassInverseIN)
:
  constraintMatrix(&ConstraintMatrixIN),
  stiffnessMatrix(&StiffnessMatrixIN),
  massInverse(&MassInverseIN),
  temp1(ConstraintMatrixIN.n()),
  temp2(ConstraintMatrixIN.n())
{}



void
SchurComplementPoisson::vmult
(dealii::BlockVector<heat::real> & destination,
 const dealii::BlockVector<heat::real> &source) const
{
  stiffnessMatrix->vmult(temp1,source.block(0));
  constraintMatrix->Tvmult_add(temp1,source.block(1));
  massInverse->vmult(temp2,temp1);
  stiffnessMatrix->Tvmult(destination.block(0), temp2);
  constraintMatrix->vmult(destination.block(1), temp2);  
}

void
SchurComplementPoisson::Tvmult
(dealii::BlockVector<heat::real> & destination,
 const dealii::BlockVector<heat::real> & source) const
{
  vmult(destination, source);
}
 
} // end namespace heat
