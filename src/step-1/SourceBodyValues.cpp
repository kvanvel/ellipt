#include "SourceBodyValues.hpp"

namespace heat {
  
template <int dim>
heat::real
SourceBodyValues<dim>::value(const dealii::Point<dim> & /*p*/,
			     const unsigned int component) const
{
  Assert(component == 0, dealii::ExcInternalError() );
  return 0.0;
  //return 2.0 * sin(p[0]) * cos(p[1]);
  
}


template <int dim>
void
SourceBodyValues<dim>
::value_list(const std::vector<dealii::Point<dim> > & points,
	     std::vector<heat::real> & values,
	     const unsigned int component ) const
{
  Assert(points.size() == values.size(),
	 dealii::ExcDimensionMismatch(points.size(), values.size() ) );
  for(unsigned int p = 0; p < points.size(); ++p){    
    values[p] = this->value(points[p],component);
  }
}
	     
} // end namespace heat


template class heat::SourceBodyValues<2>;
template class heat::SourceBodyValues<3>;
