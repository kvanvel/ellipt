//Taken Mainly from 
//http://www.dealii.org/8.0.0/doxygen/deal.II/step_23.html#Equationdata

#include "InitialValuesU.hpp"

namespace heat {
  
template <int dim>
double
InitialValuesU<dim>::value(const dealii::Point<dim> & p /*p */,
			   const unsigned int component) const
{
  Assert (component ==0, dealii::ExcInternalError() );  
  dealii::Point<dim> zero;
  //double r = p.distance(zero);  
  //return 0;
  //return sin(r);  
  //return sin(3.1415 * p(0)) * sin(3.1415 * p(1));
  // return 1.0;
  
  if((fabs(p[0]) < 1) && (fabs(p[1]) < 1)){
    double temp =  cos(.5* 3.1415 * p[0] ) * cos(.5 * 3.1415 * p[1] );
    return temp * temp;
  }
  else{
    return 0.0;
  }
  

  
}

template <int dim>
void
InitialValuesU<dim>::value_list(const std::vector<dealii::Point<dim> > & points,
				std::vector<heat::real > & values) const
{
  Assert(points.size() == values.size(),
	 dealii::ExcDimensionMismatch(points.size(), values.size() ) );
  
  for(unsigned int p = 0; p < points.size(); ++p)
    {
      this->value(points[p],values[p]);
    }
}

} //heat namespace end


template class heat::InitialValuesU<2>;
template class heat::InitialValuesU<3>;
