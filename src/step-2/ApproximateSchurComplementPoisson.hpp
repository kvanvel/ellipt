#ifndef H_APPROXIMATESCHURCOMPLEMENT_POISSON__
#define H_APPROXIMATESCHURCOMPLEMENT_POISSON__

#include "globals.hpp"
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/block_vector.h>


namespace heat{

class ApproximateSchurComplementPoisson : public dealii::Subscriptor
{  
public:
  ApproximateSchurComplementPoisson
  (const dealii::SparseMatrixEZ<heat::real> & ConstraintMatrixIN,
   const dealii::SparseMatrix<heat::real> & StiffnessMatrixIN,
   const dealii::SparseMatrix<heat::real> & MassMatrixIN);


  void vmult(dealii::BlockVector<heat::real> & destination,
	     const dealii::BlockVector<heat::real> & source) const;
  
  void Tvmult(dealii::BlockVector<heat::real> & destination,
	      const dealii::BlockVector<heat::real> & source) const;
  
private:
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > constraintMatrix;
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> > stiffnessMatrix;
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> > massMatrix;

  mutable dealii::Vector<heat::real> temp1, temp2;
};


} // end namespace heat

#endif // End Include Guard
