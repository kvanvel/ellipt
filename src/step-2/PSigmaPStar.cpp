#include "PSigmaPStar.hpp"

namespace heat{

PSigmaPStar::
PSigmaPStar(const dealii::SparseMatrixEZ<heat::real>  &P_In,
	    const dealii::SparseMatrix<heat::real> &Sigma_In)
  :
  P(&P_In),
  Sigma(&Sigma_In),
  temp1(P_In.n()),
  temp2(Sigma_In.m())
{}
  

void
PSigmaPStar::vmult(dealii::Vector<heat::real> &destination,
		   const dealii::Vector<heat::real> &source) const

{  
  P->Tvmult(temp1,source);
  Sigma->vmult(temp2,temp1);
  P->vmult(destination,temp2);
}

void
PSigmaPStar::Tvmult(dealii::Vector<heat::real> & destination,
		    const dealii::Vector<heat::real> & source) const
{
  vmult(destination,source);
}


} //end namespace heat
