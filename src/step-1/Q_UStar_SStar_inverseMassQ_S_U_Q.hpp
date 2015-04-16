#ifndef H_Q_USTAR_SStar_inverseMassQ_S_U_Q__
#define H_Q_USTAR_SStar_inverseMassQ_S_U_Q__

#include "globals.hpp"
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/iterative_inverse.h>

namespace heat{

class Q_UStar_SStar_inverseMassQ_S_U_Q : public dealii::Subscriptor
{
public: 
  Q_UStar_SStar_inverseMassQ_S_U_Q
  (const dealii::SparseMatrixEZ<heat::real> &Q_In,
   const dealii::SparseMatrix<heat::real> &U_In,
   const dealii::SparseMatrix<heat::real> &totalQFromU_In,
   const dealii::SparseMatrix<heat::real> &inverseMassQ_In,
   const dealii::SparseMatrix<heat::real> &totalUFromQ_In);


  void 
  vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const;

  void 
  Tvmult(dealii::Vector<heat::real> &destination,
	 const dealii::Vector<heat::real> &source) const;
  

private:
  const dealii::SmartPointer<const dealii::SparseMatrixEZ<heat::real> > Q;
  const dealii::SmartPointer<const dealii::SparseMatrix<heat::real> >
    U,
    totalQFromU,
    inverseMassQ,
    totalUFromQ;

  
  mutable dealii::Vector<heat::real> temp1, temp2, temp3, temp4, temp5, temp6;
    
};

} // End Namespace heat

#endif //End header Guard

