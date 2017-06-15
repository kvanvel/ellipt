#ifndef H_SCHURCOMPLEMENT_NEW__
#define H_SCHURCOMPLEMENT_NEW__

#include "globals.hpp"
#include <deal.II/lac/iterative_inverse.h>

namespace heat{

template<class MATRIX1, class MATRIX2, class MATRIX3>
class SchurComplementNEW : public dealii::Subscriptor
{
public:

  typedef dealii::types::global_dof_index size_type;
  
  SchurComplementNEW(const MATRIX1 &Ainv,
		     const MATRIX2 &B,
		     const MATRIX3 &C);
  
  void 
  vmult(dealii::Vector<heat::real> &destination,
	const dealii::Vector<heat::real> &source) const;

  void 
  Tvmult(dealii::Vector<heat::real> &destination,
	 const dealii::Vector<heat::real> &source) const;

  unsigned int m() const;

  unsigned int n() const;
  

private:
  const dealii::SmartPointer<const MATRIX1, SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3> > Ainv;
  const dealii::SmartPointer<const MATRIX2, SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3> > B;
  const dealii::SmartPointer<const MATRIX3, SchurComplementNEW<MATRIX1, MATRIX2, MATRIX3> > C;

  mutable dealii::Vector<heat::real> temp1, temp2;
  
};



} // End Namespace heat

#endif //End header Guard

