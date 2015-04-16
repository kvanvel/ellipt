#include "SchurComplement.hpp"

namespace heat{

SchurComplement::SchurComplement
(const dealii::SparseMatrixEZ<heat::real> &constraintMatrixIN,
 const dealii::IterativeInverse<dealii::Vector<heat::real> > & Minverse)
:
constraintMatrix(&constraintMatrixIN),
m_inverse(&Minverse),
temp1(constraintMatrixIN.n()),
temp2(constraintMatrixIN.n())
{}

void
SchurComplement::vmult(dealii::Vector<heat::real> &destination,
		       const dealii::Vector<heat::real> &source) const

{
  constraintMatrix->Tvmult(temp1,source);  
  dealii::deallog.push("InnerSolve");
  try {
    m_inverse->vmult(temp2, temp1);    
  } catch (...) {
    std::cerr << "Failure in the inner solve of Schur Complement" << std::endl;
    throw;
  }
  
  dealii::deallog.pop();
  constraintMatrix->vmult(destination, temp2);
}

void
SchurComplement::Tvmult(dealii::Vector<heat::real> & destination,
			const dealii::Vector<heat::real> & source) const
{
  vmult(destination,source);
}


} //end namespace heat
