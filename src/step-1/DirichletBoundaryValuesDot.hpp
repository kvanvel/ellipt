#ifndef H_DIRICHELET_BOUNDARY_VALUES_DOT__
#define H_DIRICHELET_BOUNDARY_VALUES_DOT__

#include "globals.hpp"
#include <deal.II/base/function.h>

namespace heat{

template <int dim>
class DirichletBoundaryValuesDot : public dealii::Function<dim>
{
public:
  DirichletBoundaryValuesDot() : dealii::Function<dim>(1){}
  
  using dealii::Function<dim>::value;
  using dealii::Function<dim>::value_list;
  
  virtual double
  value(const dealii::Point<dim> & /*p */,
	const unsigned int componenet = 0) const;
  
  virtual void
  value_list(const std::vector<dealii::Point<dim > > & points,
	     std::vector<heat::real > & values ) const;
  
};

} //namespace heat end
#endif
