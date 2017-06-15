#ifndef H_SS_UStar_MASSQ_U_SSSTAR__
#define H_SS_UStar_MASSQ_U_SSSTAR__

#include "globals.hpp"
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/iterative_inverse.h>

namespace heat{

class SS_UStar_massQ_U_SSStar : public dealii::Subscriptor
{
public: 
  SS_UStar_massQ_U_SSStar
  (const dealii::SparseMatrixEZ<heat::real> &SS_In,
   const dealii::SparseMatrix<heat::real> &U_In,   
   const dealii::SparseMatrix<heat::real> &massQ_In);


  unsigned int
  m(void) const;
  
  void 
  vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const;

  void 
  Tvmult(dealii::Vector<heat::real> &destination,
	 const dealii::Vector<heat::real> &source) const;
  

private:
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > SS;
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> >
  U,
    massQ;
    

  
  mutable dealii::Vector<heat::real> temp1, temp2, temp3, temp4;
    
};

} // End Namespace heat

#endif //End header Guard

