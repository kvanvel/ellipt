#ifndef H_Q_USTAR_MASSU_U_QSTAR__
#define H_Q_USTAR_MASSU_U_QSTAR__


#include "globals.hpp"
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/iterative_inverse.h>

namespace heat{

class Q_UStar_massU_U_QStar : public dealii::Subscriptor
{
public: 
  Q_UStar_massU_U_QStar
  (const dealii::SparseMatrixEZ<heat::real> &Q_In,
   const dealii::SparseMatrix<heat::real> &U_In,
   const dealii::SparseMatrix<heat::real> &massU_In);
   

  unsigned int
  m(void) const;
  
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
    massU;


  
  mutable dealii::Vector<heat::real> temp1, temp2, temp3, temp4;
    
};

} // End Namespace heat

#endif //End header Guard

