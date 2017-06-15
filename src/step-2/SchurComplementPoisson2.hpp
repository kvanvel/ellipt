#ifndef H_SCHURCOMPLEMENTPOISSON2__
#define H_SCHURCOMPLEMENTPOISSON2__

#include "globals.hpp"
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/iterative_inverse.h>
#include <deal.II/lac/block_vector.h>

namespace heat{

class SchurComplementPoisson : public dealii::Subscriptor
{
public:
  SchurComplementPoisson
  (const dealii::SparseMatrixEZ<heat::real> & constraintMatrix,
   const dealii::SparseMatrix<heat::real> & StiffnessMatrix,
   const dealii::IterativeInverse<dealii::Vector<heat::real> > & MassInverse);
  
  void
  vmult(dealii::BlockVector<heat::real> & destination,
	const dealii::BlockVector<heat::real> & source) const;
  
  void
  Tvmult(dealii::BlockVector<heat::real> & destination,
	 const dealii::BlockVector<heat::real> & source) const;

private:
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > constraintMatrix;
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> > stiffnessMatrix;
  const dealii::SmartPointer<const dealii::IterativeInverse<dealii::Vector<heat::real> > > massInverse;
  mutable dealii::Vector<heat::real> temp1, temp2;
};

}

#endif // End header Guard
