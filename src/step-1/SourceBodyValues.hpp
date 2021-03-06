#ifndef H_SOURCE_BODY_VALUES__
#define H_SOURCE_BODY_VALUES__

#include "globals.hpp"
#include <deal.II/base/function.h>

namespace heat{

template <int dim>
class SourceBodyValues : public dealii::Function<dim,heat::real>
{
public:
  SourceBodyValues() : dealii::Function<dim,heat::real>(1){}

  //using dealii::Function<dim>::value;

  heat::real
  value(const dealii::Point<dim> & /*p */,
	const unsigned int component = 0) const;
  
  //using dealii::Function<dim>::value_list;

  void
  value_list(const std::vector<dealii::Point<dim > > & points,
	     std::vector<heat::real> & values,
	     const unsigned int component = 0) const;
  
};

} //namespace heat end
#endif
