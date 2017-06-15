#include "LDGIntegratorDrift.hpp"
#include "LDG.hpp"

namespace heat{
namespace LDGIntegratorDrift{

template<int dim>
LDGIntegratorDrift<dim>::LDGIntegratorDrift(
					    const dealii::Point<dim> referenceDirection_In,
					    const dealii::Quadrature<dim-1> & face_quadrature_In,
					    const dealii::UpdateFlags & face_update_flags_In,
					    const dealii::MappingQ1<dim,dim>& mapping_In,
					    const dealii::Vector<heat::real> & ElecState_In)
					    
  :
  dealii::MeshWorker::LocalIntegrator<dim>(true, true, true),
  referenceDirection( referenceDirection_In),
  face_quadrature(& face_quadrature_In),
  face_update_flags( face_update_flags_In ),
  mapping(mapping_In),
  ElecState(& ElecState_In)
  
{}

template<int dim>
void
LDGIntegratorDrift<dim>::cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
			      dealii::MeshWorker::IntegrationInfo<dim> &info ) const
{  
  const heat::real minusOne = -1.0;
  
  heat::LocalIntegrators::LDG::StiffUFromUMat(info.fe_values(1),
					      info.fe_values(0),
					      *ElecState,					      
					      dinfo.matrix(0, false),
					      minusOne);
}


template <int dim>
void
LDGIntegratorDrift<dim>::face(dealii::MeshWorker::DoFInfo<dim> & dinfoSELF,
			     dealii::MeshWorker::DoFInfo<dim> & dinfoNEIG,
			     dealii::MeshWorker::IntegrationInfo<dim> &infoSELF,
			     dealii::MeshWorker::IntegrationInfo<dim> &infoNEIG ) const {

  const heat::real One = 1.0;
  heat::LocalIntegrators::LDG::FluxUFromUMat(infoSELF.fe_values(1),
					     infoNEIG.fe_values(1),
					     infoSELF.fe_values(0),
					     infoNEIG.fe_values(0),
					     *ElecState,
					     dinfoSELF.matrix(0, false),
					     dinfoSELF.matrix(0, true),
					     dinfoNEIG.matrix(0, true),
					     dinfoNEIG.matrix(0, false),
					     One);
}
  
template<int dim>
void
LDGIntegratorDrift<dim>::boundary(dealii::MeshWorker::DoFInfo<dim> & dinfo,
				  dealii::MeshWorker::IntegrationInfo<dim> & info) const
{
  
  const dealii::types::boundary_id dir_plus =  1;
  const dealii::types::boundary_id dir_minus = 0;
  const dealii::types::boundary_id neu_plus  = 3;
  const dealii::types::boundary_id neu_minus = 2;
  
  switch (dinfo.face->boundary_id() ){
  
  case dir_minus:
    heat::
      LocalIntegrators::
      LDG::
      U_MinusRHS_FromDirBC(referenceDirection,
			   info.fe_values(1),
			   info.fe_values(0),
			   *ElecState,
			   dinfo.vector(0).block(0),
			   1.0);
    break;
  case dir_plus:
    heat
      ::LocalIntegrators
      ::LDG
      ::FluxUFromUBoundaryMat( referenceDirection,
			       info.fe_values(1),
			       info.fe_values(0),
			       *ElecState,
			       dinfo.matrix(0,false),
			       1.0);
    break;
  case neu_plus:
    break;
  case neu_minus:
    break;  
  }
  
  // //AssertDimension(dinfo.n_matrices(), 16);

  // const heat::real One = 1.0;

  // switch( dinfo.face->boundary_id()){
    
  // case dir_plus:
  //   //std::cout << "on dir plus" << std::endl;
  //   heat::LocalIntegrators::LDG::numericalFluxUFromQBoundary
  //     (referenceDirection,
  //      dinfo.matrix(8,false),
  //      info.fe_values(0),
  //      info.fe_values(1),
  //      One
  //      );

  //   heat::LocalIntegrators::LDG::Q_MinusRHS_FromDirPlus
  //     (referenceDirection,
  //      info.fe_values(1),
  //      dinfo.vector(0).block(1),
  //      One);


  //   break;
  // case dir_minus:
  //   heat::LocalIntegrators::LDG::numericalFluxQFromUBoundary
  //     (referenceDirection,
  //      dinfo.matrix(11,false),
  //      info.fe_values(1),
  //      info.fe_values(0),
  //      One);

  //   heat::LocalIntegrators::LDG::ConstraintUDirMinus(referenceDirection,
  // 						     info.fe_values(0),
  // 						     dinfo.vector(1).block(0),
  // 						     One);

  //   break;
    
  // case neu_plus:
  //    heat::LocalIntegrators::LDG::numericalFluxUFromQBoundary
  //     (referenceDirection,
  //      dinfo.matrix(8,false),
  //      info.fe_values(0),
  //      info.fe_values(1),
  //      One
  //      );

  //    heat::LocalIntegrators::LDG::ConstraintQNeuPlus(referenceDirection,
  // 						     info.fe_values(1),
  // 						     dinfo.vector(1).block(1),
  // 						     One);

  //   //std::cout << "on neu plus" << std::endl;
  //   // heat::LocalIntegrators::LDG::BoundaryMassQ(referenceDirection,
  //   // 					       dinfo.matrix(10,false),
  //   // 					       info.fe_values(1),
  //   // 					       One);
  //   break;
  // case neu_minus:
  //   //std::cout << "on nue minus" << std::endl;
  //   heat::LocalIntegrators::LDG::numericalFluxQFromUBoundary
  //     (referenceDirection,
  //      dinfo.matrix(11,false),
  //      info.fe_values(1),
  //      info.fe_values(0),
  //      One);


  //   heat::LocalIntegrators::LDG::U_MinusRHS_FromNeuMinus(referenceDirection,
  // 							 info.fe_values(0),
  // 							 dinfo.vector(0).block(0),
  // 							 One);
    
  //   break;
  // default:
  //   std::cout << "Something bad happened" << std::endl;
  //   assert(false);
  //   break;    
  // }
}






// template<int dim>
// JumpQ<dim>::JumpQ()
//   :
//   dealii::MeshWorker::LocalIntegrator<dim>(false, false, true) {}

// template<int dim>
// void
// JumpQ<dim>::face( dealii::MeshWorker::DoFInfo<dim> & dinfo0,
// 		  dealii::MeshWorker::DoFInfo<dim> & dinfo1,
// 		  dealii::MeshWorker::IntegrationInfo<dim> &info0,
// 		  dealii::MeshWorker::IntegrationInfo<dim> &info1) const
// {
//   AssertDimension(dinfo0.n_matrices(), 4);
//   AssertDimension(dinfo1.n_matrices(), 4);

//   heat::real one = 1.0;

//   dealii
//     ::LocalIntegrators
//     ::L2::jump_matrix(dinfo0.matrix(3,false),
// 		      dinfo0.matrix(3,true),
// 		      dinfo1.matrix(3,false),
// 		      dinfo1.matrix(3,true),
// 		      info0.fe_values(1),
// 		      info1.fe_values(1),
// 		      one,
// 		      one
// 		      );
// }



} // End Namespace LDGIntegratorDrift
} // End Namespace heat

template class heat::LDGIntegratorDrift::LDGIntegratorDrift<2>;
template class heat::LDGIntegratorDrift::LDGIntegratorDrift<3>;
