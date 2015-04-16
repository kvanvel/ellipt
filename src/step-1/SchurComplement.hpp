#ifndef H_SCHURCOMPLEMENT__
#define H_SCHURCOMPLEMENT__

#include "globals.hpp"
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/iterative_inverse.h>

namespace heat{

class SchurComplement : public dealii::Subscriptor
{
public:
  SchurComplement(const dealii::SparseMatrixEZ<heat::real>  &constraintMatrix,
		  const dealii::IterativeInverse<dealii::Vector<heat::real> > &Minv);
  void 
  vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const;

  void 
  Tvmult(dealii::Vector<heat::real> &destination,
	 const dealii::Vector<heat::real> &source) const;
  

private:
  //const dealii::SmartPointer<const dealii::BlockSparseMatrix<heat::real> > system_matrix;
  //  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> > massMatrix;
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > constraintMatrix;
  const dealii::SmartPointer<const dealii::IterativeInverse<dealii::Vector<heat::real> > > 
  m_inverse;
  
  mutable dealii::Vector<heat::real> temp1, temp2;
    
};

} // End Namespace heat

#endif //End header Guard

