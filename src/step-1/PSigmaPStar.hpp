#ifndef H_PSIGMAPSTAR__
#define H_PSIGMAPSTAR__

#include "globals.hpp"
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/iterative_inverse.h>

namespace heat{

class PSigmaPStar : public dealii::Subscriptor
{
public: 
  PSigmaPStar(const dealii::SparseMatrixEZ<heat::real>  &P_In,
	      const dealii::SparseMatrix<heat::real> &Sigma_In);
    

  void 
  vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const;

  void 
  Tvmult(dealii::Vector<heat::real> &destination,
	 const dealii::Vector<heat::real> &source) const;
  

private:
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > P;
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> > Sigma;

  
  mutable dealii::Vector<heat::real> temp1, temp2;
    
};

} // End Namespace heat

#endif //End header Guard

