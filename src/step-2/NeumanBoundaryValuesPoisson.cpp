#include "NeumanBoundaryValuesPoisson.hpp"

namespace heat {
  
template <int dim>
double
NeumanBoundaryValuesPoisson<dim>::value(const dealii::Point<dim> & p,
					 const unsigned int component) const
{
  Assert(component == 0, dealii::ExcInternalError() );
  return 0.0;
}


template <int dim>
void
NeumanBoundaryValuesPoisson<dim>
::value_list(const std::vector<dealii::Point<dim> > & points,
	     std::vector<heat::real > & values ) const
{
  Assert(points.size() == values.size(),
	 dealii::ExcDimensionMismatch(points.size(), values.size() ) );
  for(unsigned int p = 0; p < points.size(); ++p){
    this->value(points[p],values[p]);
  }
}
	     
} // end namespace heat

template class heat::NeumanBoundaryValuesPoisson<2>;
template class heat::NeumanBoundaryValuesPoisson<3>;
