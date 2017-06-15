// Mostly Copied from 
//http://www.dealii.org/8.0.0/doxygen/deal.II/step_20.html#Theinversepermeabilitytensor
#include "InverseDiffusivity.hpp"

template <int dim>
void
heat
::InverseDiffusivity<dim>::value (const dealii::Point<dim> & point,
				  dealii::Tensor<2,dim> & value) const
{
  
  value.clear();  
  for( unsigned int d = 0; d<dim; ++d){
    value[d][d] = 1.0;
  }
  
}

template <int dim>
void
heat
::InverseDiffusivity<dim>
::value_list(const std::vector< dealii::Point< dim > > & points,
	     std::vector< dealii::Tensor<2,dim> >& values) const
{
  Assert( points.size() == values.size(),
	  dealii::ExcDimensionMismatch(points.size(), values.size() ) );
    		  
  for (unsigned int p = 0; p< points.size(); ++p){
    values[p].clear();
    this->value(points[p],values[p]);
  }
}

template class heat::InverseDiffusivity<2>;
template class heat::InverseDiffusivity<3>;
