#include "LDGIntegrator.hpp"
#include "LDG.hpp"

namespace heat{
namespace LDGIntegrator{

template<int dim>
LDGIntegrator<dim>::LDGIntegrator(dealii::Point<dim> referenceDirection_In)
  :
  dealii::MeshWorker::LocalIntegrator<dim>(true, true, true),
  referenceDirection(referenceDirection_In)
  
  
{}

template<int dim>
void
LDGIntegrator<dim>::cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
			 dealii::MeshWorker::IntegrationInfo<dim> &info ) const
{  
  
  const heat::real one = 1.0;
  const heat::real minusOne = -1.0;

  const dealii::types::boundary_id dir_plus =  1;
  const dealii::types::boundary_id dir_minus = 0;
  const dealii::types::boundary_id neu_plus  = 3;
  const dealii::types::boundary_id neu_minus = 2;
  
  dealii::LocalIntegrators::L2::mass_matrix
    (dinfo.matrix(0,false), info.fe_values(0), one);
  {
    dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;
    tempLAPACKMatrix.copy_from(dinfo.matrix(0,false).matrix);
    tempLAPACKMatrix.invert();
    dinfo.matrix(4,false).matrix = tempLAPACKMatrix;
  }

  

  dealii::LocalIntegrators::L2::mass_matrix
    (dinfo.matrix(3,false), info.fe_values(1), one);
  {
    dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;
    tempLAPACKMatrix.copy_from(dinfo.matrix(3,false).matrix);
    tempLAPACKMatrix.invert();
    dinfo.matrix(7,false).matrix = tempLAPACKMatrix;
  }
  
  dealii::LocalIntegrators::Divergence::gradient_matrix(dinfo.matrix(2, false).matrix,
							info.fe_values(0),
							info.fe_values(1),
							minusOne);

  dealii::LocalIntegrators::Divergence::cell_matrix(dinfo.matrix(1,false).matrix,
						    info.fe_values(1),
						    info.fe_values(0),
						    minusOne);
      
}



  




template <int dim>
void
LDGIntegrator<dim>::face(dealii::MeshWorker::DoFInfo<dim> & dinfoSELF,
			 dealii::MeshWorker::DoFInfo<dim> & dinfoNEIG,
			 dealii::MeshWorker::IntegrationInfo<dim> &infoSELF,
			 dealii::MeshWorker::IntegrationInfo<dim> &infoNEIG ) const {


  AssertDimension(dinfoSELF.n_matrices(), 16);
  AssertDimension(dinfoNEIG.n_matrices(), 16);

  const heat::real One = 1.0;

  heat::LocalIntegrators::LDG::numericalFluxQFromU
    (
     referenceDirection,
     dinfoSELF.matrix(2,false),
     dinfoSELF.matrix(2,true),
     dinfoNEIG.matrix(2,true),
     dinfoNEIG.matrix(2,false),
     infoSELF.fe_values(1),
     infoNEIG.fe_values(1),
     infoSELF.fe_values(0),
     infoNEIG.fe_values(0),
     One);


  heat::LocalIntegrators::LDG::numericalFluxUFromQ
    (
     referenceDirection,
     dinfoSELF.matrix(1,false),
     dinfoSELF.matrix(1,true),
     dinfoNEIG.matrix(1,true),
     dinfoNEIG.matrix(1,false),
     infoSELF.fe_values(0),
     infoNEIG.fe_values(0),
     infoSELF.fe_values(1),
     infoNEIG.fe_values(1),
     One);
  
}
  
template<int dim>
void
LDGIntegrator<dim>::boundary(dealii::MeshWorker::DoFInfo<dim> & dinfo,
			     dealii::MeshWorker::IntegrationInfo<dim> & info) const
{
  const dealii::types::boundary_id dir_plus =  1;
  const dealii::types::boundary_id dir_minus = 0;
  const dealii::types::boundary_id neu_plus  = 3;
  const dealii::types::boundary_id neu_minus = 2;

  
  AssertDimension(dinfo.n_matrices(), 16);

  const heat::real One = 1.0;



  switch( dinfo.face->boundary_id()){
    
  case dir_plus:
    std::cout << "on dir plus" << std::endl;
    heat::LocalIntegrators::LDG::numericalFluxUFromQBoundary
      (referenceDirection,
       dinfo.matrix(1,false),
       info.fe_values(0),
       info.fe_values(1),
       One
       );

    break;
  case dir_minus:
    std::cout << "on dir minus" << std::endl;
    heat::LocalIntegrators::LDG::BoundaryMassU(referenceDirection,
    					       dinfo.matrix(8,false),
					       info.fe_values(0),
    					       One);    

    break;
  case neu_plus:
    std::cout << "on neu plus" << std::endl;
    heat::LocalIntegrators::LDG::BoundaryMassQ(referenceDirection,
					       dinfo.matrix(11,false),
					       info.fe_values(1),
					       One);
    break;
  case neu_minus:
    std::cout << "on nue minus" << std::endl;
    heat::LocalIntegrators::LDG::numericalFluxQFromUBoundary
    (referenceDirection,
     dinfo.matrix(2,false),
     info.fe_values(1),
     info.fe_values(0),
     One);    
    break;
  default:
    std::cout << "Something bad happened" << std::endl;
    assert(false);
    break;

    
  }
}



// MassQ
template <int dim>
MassQ<dim>::MassQ()
  :
  dealii::MeshWorker::LocalIntegrator<dim>(true, false, false) {}
  
template <int dim>
void
MassQ<dim>::cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
		 dealii::MeshWorker::IntegrationInfo<dim> &info ) const
{
  AssertDimension(dinfo.n_matrices(), 1);

  const double coefficient = 1.0;
  
  dealii::LocalIntegrators::L2::mass_matrix
    (dinfo.matrix(0,false).matrix, info.fe_values(0), coefficient);
}

// MassU
template <int dim>
MassU<dim>::MassU()
  :
  dealii::MeshWorker::LocalIntegrator<dim>(true, false, false) {}
  
template <int dim>
void
MassU<dim>::cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
		 dealii::MeshWorker::IntegrationInfo<dim> &info ) const
{
  AssertDimension(dinfo.n_matrices(), 1);

  const double coefficient = 1.0;
  
  dealii::LocalIntegrators::L2::mass_matrix
    (dinfo.matrix(0,false).matrix, info.fe_values(0), coefficient);
}

template<int dim>
StiffUFromQ<dim>::StiffUFromQ()
  :
  dealii::MeshWorker::LocalIntegrator<dim>(true, false, false) {}



template<int dim>
void
StiffUFromQ<dim>::cell
(dealii::MeshWorker::DoFInfo<dim> & dinfo,
 dealii::MeshWorker::IntegrationInfo<dim> &info ) const
{

  AssertDimension(dinfo.n_matrices(), 4);

  heat::real one = 1.0;
  
  
  dealii
    ::LocalIntegrators
    ::Divergence
    ::gradient_matrix
    (
     dinfo.matrix(2, false).matrix,
     info.fe_values(0),
     info.fe_values(1),
     one);

  
}



template<int dim>
JumpQ<dim>::JumpQ()
  :
  dealii::MeshWorker::LocalIntegrator<dim>(false, false, true) {}

template<int dim>
void
JumpQ<dim>::face( dealii::MeshWorker::DoFInfo<dim> & dinfo0,
		  dealii::MeshWorker::DoFInfo<dim> & dinfo1,
		  dealii::MeshWorker::IntegrationInfo<dim> &info0,
		  dealii::MeshWorker::IntegrationInfo<dim> &info1) const
{
  AssertDimension(dinfo0.n_matrices(), 4);
  AssertDimension(dinfo1.n_matrices(), 4);

  heat::real one = 1.0;

  dealii
    ::LocalIntegrators
    ::L2::jump_matrix(dinfo0.matrix(3,false),
		      dinfo0.matrix(3,true),
		      dinfo1.matrix(3,false),
		      dinfo1.matrix(3,true),
		      info0.fe_values(1),
		      info1.fe_values(1),
		      one,
		      one
		      );
}



} // End Namespace LDGIntegrator
} // End Namespace heat

template class heat::LDGIntegrator::MassQ<2>;
template class heat::LDGIntegrator::MassQ<3>;
template class heat::LDGIntegrator::MassU<2>;
template class heat::LDGIntegrator::MassU<3>;
template class heat::LDGIntegrator::StiffUFromQ<2>;
template class heat::LDGIntegrator::StiffUFromQ<3>;
template class heat::LDGIntegrator::JumpQ<2>;
template class heat::LDGIntegrator::JumpQ<3>;
template class heat::LDGIntegrator::LDGIntegrator<2>;
template class heat::LDGIntegrator::LDGIntegrator<3>;
