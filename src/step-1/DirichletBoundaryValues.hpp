#ifndef H_DIRICHELET_BOUNDARY_VALUES__
#define H_DIRICHELET_BOUNDARY_VALUES__

#include "globals.hpp"
#include <deal.II/base/function.h>

namespace heat{

template <int dim>
class DirichletBoundaryValues : public dealii::Function<dim>
{
public:
  DirichletBoundaryValues() : dealii::Function<dim>(1){}

  using dealii::Function<dim>::value;
  virtual double
  value(const dealii::Point<dim> & /*p */,
	const unsigned int componenet = 0) const;
  
  using dealii::Function<dim>::value_list;
  virtual void
  value_list(const std::vector<dealii::Point<dim > > & points,
	     std::vector<heat::real > & values ) const;
  
};

} //namespace heat end
#endif
