#include "globals.hpp"
#include "EllipticProblem.hpp"
#include <deal.II/base/logstream.h>

int 
main ()
{
  dealii::deallog.depth_console(0);
  
  for(unsigned int order = 1; order <= 8; ++order){
    std::cout << "Order = " << order << std::endl;
    for(unsigned int refinements = 2; refinements <= 2; ++ refinements){
      std::cout << "refinements = " << refinements << std::endl;
      heat::EllipticProblem<2> EllipticProblem(order,refinements);
      EllipticProblem.run();
    }
  }
  return 0;
}

