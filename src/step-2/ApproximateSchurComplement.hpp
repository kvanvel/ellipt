#ifndef H_APPROXIMATESCHURCOMPLEMENT__
#define H_APPROXIMATESCHURCOMPLEMENT__

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include "globals.hpp"


namespace heat{

class ApproximateSchurComplement : public dealii::Subscriptor
{
public:
  ApproximateSchurComplement
  (const dealii::SparseMatrix<heat::real> &AIn,
   const dealii::SparseMatrixEZ<heat::real> &ConstraintMatrixIN);

  void vmult(dealii::Vector<heat::real> & destination,
	     const dealii::Vector<heat::real> & source) const;
  
  void Tvmult(dealii::Vector<heat::real> & destination,
	      const dealii::Vector<heat::real> & source) const;
  
private:
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> > A;
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > ConstraintMatrix;
  mutable dealii::Vector<heat::real> temp1, temp2;
};


} // end namespace heat

#endif // End Include Guard
