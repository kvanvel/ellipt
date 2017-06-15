#include "NeumannBoundaryValues.hpp"

namespace heat {

// template <int dim>
// NeumannBoundaryValuesNEW<dim>::NeumannBoundaryValuesNEW()
//   :
//   dealii::TensorFunction<1,dim>() {}



template<int dim>
dealii::Tensor<1,dim,heat::real> 
NeumannBoundaryValues<dim>::value(const dealii::Point<dim> & p) const

{
  assert(2 == dim);

  double x = p(0);
  double y = p(1);
  double phi = std::atan2(y,-x)+M_PI;
  double r43 = std::pow(x*x+y*y,2./3.);

  dealii::Tensor<1,dim, heat::real> result;
  result[0] = -2./3.*(std::sin(2./3.*phi)*x + std::cos(2./3.*phi)*y)/r43;
  result[1] = -2./3.*(std::sin(2./3.*phi)*y - std::cos(2./3.*phi)*x)/r43;
  return result;

  //return_value[0] = -cos(point[0]) * cos(point[1]);
  //return_value[1] =  sin(point[0]) * sin(point[1]);    
}

template<int dim>
void
NeumannBoundaryValues<dim>::value_list(const std::vector< dealii::Point<dim> > & points,
				       std::vector< dealii::Tensor<1,dim,heat::real> > & values) const
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
