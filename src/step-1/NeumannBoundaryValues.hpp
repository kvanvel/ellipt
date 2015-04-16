#ifndef H_NEUMANN_BOUNDARY_VALUES__
#define H_NEUMANN_BOUNDARY_VALUES__

#include "globals.hpp"
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>

namespace heat{

template <int dim>
class NeumannBoundaryValues : public dealii::TensorFunction<1,dim>
{
public:
  NeumannBoundaryValues() : dealii::TensorFunction<1,dim>(){}

  //using dealii::Function<dim>::vector_value;
  //using dealii::Function<dim>::vector_value_list;

  virtual
  dealii::Tensor<1,dim>
  value(const dealii::Point<dim> & point) const;
	

  void
  value_list(const std::vector< dealii::Point<dim> > & point,
	     std::vector< dealii::Tensor<1,dim> > & values) const;
  
};

} //namespace heat end

#endif
