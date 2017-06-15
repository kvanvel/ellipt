#ifndef H_SS_UStar_TOTALQFROMU_U_QQSTAR__
#define H_SS_UStar_TOTALQFROMU_U_QQSTAR__

#include "globals.hpp"
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/iterative_inverse.h>

namespace heat{

class SS_UStar_TotalQFromU_U_QQStar : public dealii::Subscriptor
{
public: 
  SS_UStar_TotalQFromU_U_QQStar
  (const dealii::SparseMatrixEZ<heat::real> &SS_In,
   const dealii::SparseMatrix<heat::real> &U_from_Neu_In,   
   const dealii::SparseMatrix<heat::real> &TotalQFromU_In,
   const dealii::SparseMatrix<heat::real> &U_from_Dir_In,
   const dealii::SparseMatrixEZ<heat::real> &QQ_IN);

  unsigned int
  m(void) const;

  unsigned int
  n(void) const;
  
  void 
  vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const;

  void 
  Tvmult(dealii::Vector<heat::real> &destination,
	 const dealii::Vector<heat::real> &source) const;
  

private:
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > SS;
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> >
  U_from_Neu, TotalQFromU, U_from_Dir;
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > QQ;
    

  
  mutable dealii::Vector<heat::real> temp1, temp2, temp3, temp4;
    
};

} // End Namespace heat

#endif //End header Guard

