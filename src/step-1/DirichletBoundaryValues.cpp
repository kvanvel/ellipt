#include "DirichletBoundaryValues.hpp"

using namespace dealii;  

namespace heat {
  
template <int dim>
double
DirichletBoundaryValues<dim>::value(const dealii::Point<dim> & p,
				    const unsigned int component) const
{
  Assert(component == 0, dealii::ExcInternalError() );
  
  double x = p(0);
  double y = p(1);

  if ((x>=0) && (y>=0))
    return 0.;
  
  double phi = std::atan2(y,-x)+M_PI;
  double r2 = x*x+y*y;

  double something = std::pow(r2,1./3.) * std::sin(2./3.*phi);
  AssertIsFinite(something);

  return something;

  // auto something = typename dealii::Functions::LSingularityFunction.value(p, component);
  // return something;
  //return 1.0;
  
  //return sin( p(0) ) * cos( p(1) );
  // return pow(pow(p(0),2) + pow(p(1),2),2.0/3.0)*
  //   sin((2*atan2(-p(1),-p(0)))/3.);
  
}


template <int dim>
void
DirichletBoundaryValues<dim>
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

template class heat::DirichletBoundaryValues<2>;
template class heat::DirichletBoundaryValues<3>;
