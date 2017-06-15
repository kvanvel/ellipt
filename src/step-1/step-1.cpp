#include "globals.hpp"
#include "EllipticProblem.hpp"
#include <deal.II/base/logstream.h>

//#include <grvy.h>

int 
main ()
{
  //grvy_timer_init("Something");
  //grvy_timer_begin("Main Program");
  dealii::deallog.depth_console(0);  
  
  for(unsigned int order = 3; order <= 3;  ++order){
    //dealii::ConvergenceTable convergenceTable;

    std::cout << "Order = " << order << std::endl;
    
    for(unsigned int refinements = 2; refinements <= 2;	++refinements){
      if( (0 == order) && ( 0 == refinements)){
	continue;
      }

      //std::cout << "Refinements = " << refinements << std::endl;
      
      heat::EllipticProblem<2> EllipticProblem(order,refinements);
      EllipticProblem.run();

      
    }
    //grvy_timer_end("Main Program");
    //grvy_timer_finalize();
    //grvy_timer_summarize();
  }
  return 0;
}

