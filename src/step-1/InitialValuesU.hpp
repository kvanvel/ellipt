#ifndef H_INITIAL_VALUES_U__
#define H_INITIAL_VALUES_U__

#include "globals.hpp"
#include <deal.II/base/function.h>

namespace heat{ 

template <int dim>
class InitialValuesU : public dealii::Function<dim>
{
public:
  InitialValuesU() : dealii::Function<dim>(1) {}

  using dealii::Function<dim>::value;
  using dealii::Function<dim>::value_list;
  
  virtual double
  value(const dealii::Point<dim> & /*p */,
	const unsigned int component = 0) const;
    
  virtual void
  value_list(const std::vector<dealii::Point< dim > > & points,
	     std::vector< heat::real > & values ) const;

  
};

} // heat namespace end



#endif

