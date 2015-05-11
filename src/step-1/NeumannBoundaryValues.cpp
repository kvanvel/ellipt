#include "NeumannBoundaryValues.hpp"

namespace heat {

// template <int dim>
// NeumannBoundaryValuesNEW<dim>::NeumannBoundaryValuesNEW()
//   :
//   dealii::TensorFunction<1,dim>() {}



template<int dim>
dealii::Tensor<1,dim> 
NeumannBoundaryValues<dim>::value(const dealii::Point<dim> & point) const

{
  dealii::Tensor<1,dim> return_value;
  return_value[0] = -cos(point[0]) * cos(point[1]);
  return_value[1] =  sin(point[0]) * sin(point[1]);

  return return_value;    
    
}

template<int dim>
void
NeumannBoundaryValues<dim>::value_list(const std::vector< dealii::Point<dim> > & points,
				             std::vector< dealii::Tensor<1,dim> > & values) const
{
  Assert(points.size() == values.size(),
	 dealii::ExcDimensionMismatch(points.size(), values.size() ) );
  for(unsigned int p = 0; p < points.size(); ++p){
    values[p] = ( this->value(points[p]) ) ;
  }
}
  
} // end namespace heat

template class heat::NeumannBoundaryValues<2>;
template class heat::NeumannBoundaryValues<3>;
