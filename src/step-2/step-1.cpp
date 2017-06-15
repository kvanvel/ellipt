#include "globals.hpp"
#include "EllipticProblem.hpp"
#include <deal.II/base/logstream.h>
#include <deal.II/base/convergence_table.h>
#include <grvy.h>

int 
main ()  
{

  //grvy_timer_init("Something");
  //grvy_timer_begin("Main Program");
  dealii::deallog.depth_console(-1);
  
  for(unsigned int order = 1; order <= 1; ++order){
    dealii::ConvergenceTable convergenceTable;

    std::cout << "Order = " << order << std::endl;
    
    for(unsigned int refinements = 1; refinements <= 1; ++refinements){
      if( (0 == order) && ( 0 == refinements)){
	continue;
      }
      std::cout << "Refinements = " << refinements << std::endl;
      convergenceTable.add_value("refinements",refinements);

      heat::EllipticProblem<2> EllipticProblem(order,refinements);

      EllipticProblem.run();
      
      convergenceTable.add_value("L2errorU", EllipticProblem.L2errorU);
      convergenceTable.add_value("L2errorQ", EllipticProblem.L2errorQ);
    }

    std::cout << "CONVERGENCE TABLE" << std::endl;
    
    convergenceTable.set_precision("L2errorU", 3);
    convergenceTable.set_precision("L2errorQ", 3);
    
    convergenceTable.set_scientific("L2errorU", true);
    convergenceTable.set_scientific("L2errorQ", true);
    
    convergenceTable.evaluate_convergence_rates
      ("L2errorU", dealii::ConvergenceTable::reduction_rate_log2);
    convergenceTable.evaluate_convergence_rates
      ("L2errorQ", dealii::ConvergenceTable::reduction_rate_log2);
    
    convergenceTable.write_text(std::cout);
    //grvy_timer_end("Main Program");
    //grvy_timer_finalize();
    //grvy_timer_summarize();
  }
  return 0;
}

