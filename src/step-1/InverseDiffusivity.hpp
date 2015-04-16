#ifndef H_INVERSE_DIFFUSIVITY__
#define H_INVERSE_DIFFUSIVITY__

#include "globals.hpp"
#include <deal.II/base/tensor_function.h>

namespace heat{

template <int dim>
class InverseDiffusivity : public dealii::TensorFunction<2,dim>
{
public:
  InverseDiffusivity() : dealii::TensorFunction<2,dim>() {}

  using dealii::TensorFunction<2,dim>::value;
  void value (const dealii::Point<dim> & point,
	      dealii::Tensor<2,dim> & value) const;
  
  using dealii::TensorFunction<2,dim>::value_list;
  void value_list(const std::vector< dealii::Point< dim > > & points,
		  std::vector< dealii::Tensor<2,dim> > & values) const;
};




} // heat namespace end



#endif
