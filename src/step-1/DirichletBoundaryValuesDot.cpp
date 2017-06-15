#include "DirichletBoundaryValuesDot.hpp"
#include <iostream>
namespace heat {
  
template <int dim>
double
DirichletBoundaryValuesDot<dim>::value(const dealii::Point<dim> & /*p*/,
				       const unsigned int component) const
{
  Assert(component == 0, dealii::ExcInternalError() );  
  //std::cout << "time = " << this->get_time() << std::endl;
  const heat::real time = this->get_time();
  //return time;
  return sin(3.14 * time);
  //return 0.0;
}


template <int dim>
void
DirichletBoundaryValuesDot<dim>
::value_list(const std::vector<dealii::Point<dim> > & points,
	     std::vector<heat::real > & values ) const
{
  Assert(points.size() == values.size(),
	 dealii::ExcDimensionMismatch(points.size(), values.size() ) );
  for(unsigned int p = 0; p < points.size(); ++p)
    {
      this->value(points[p],values[p]);
    }
}
	     
} // end namespace heat


template class heat::DirichletBoundaryValuesDot<2>;
template class heat::DirichletBoundaryValuesDot<3>;
