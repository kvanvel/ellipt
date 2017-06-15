#include "LDGIntegratorKENT.hpp"
#include "LDG.hpp"

namespace heat{
namespace LDGIntegratorKENT{

template<int dim>
LDGIntegratorKENT<dim>::LDGIntegratorKENT(dealii::Point<dim> referenceDirection_In,
					  const dealii::Quadrature<dim-1> & face_quadrature_In,
					  const dealii::UpdateFlags & face_update_flags_In,
					  const dealii::MappingQ1<dim,dim>& mapping_In)
  :
  dealii::MeshWorker::LocalIntegrator<dim>(true, true, true),
  referenceDirection(referenceDirection_In),
  face_quadrature(& face_quadrature_In),
  face_update_flags( face_update_flags_In ),
  mapping(mapping_In)
  
{}

template<int dim>
void
LDGIntegratorKENT<dim>::cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
			     dealii::MeshWorker::IntegrationInfo<dim> &info ) const
{  
  const heat::real one = 1.0;
  const heat::real minusOne = -1.0;

  const dealii::types::boundary_id dir_plus =  1;
  const dealii::types::boundary_id dir_minus = 0;
  const dealii::types::boundary_id neu_plus  = 3;
  const dealii::types::boundary_id neu_minus = 2;

  // std::cout << "CELL dinfo.matrix(9,false).matrix = " << std::endl;
  
  // dinfo.matrix(9,false).matrix.print(std::cout);
  
  
  //massU  
  dealii::LocalIntegrators::L2::mass_matrix
    (dinfo.matrix(0,false), info.fe_values(0), one);  
  {
    dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;
    tempLAPACKMatrix.copy_from(dinfo.matrix(0,false).matrix);
    tempLAPACKMatrix.invert();
    dinfo.matrix(1,false).matrix = tempLAPACKMatrix;
  }  

  dealii::LocalIntegrators::L2::mass_matrix
    (dinfo.matrix(2,false), info.fe_values(1), one);
  {
    dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;
    tempLAPACKMatrix.copy_from(dinfo.matrix(2,false).matrix);
    tempLAPACKMatrix.invert();
    dinfo.matrix(3,false).matrix = tempLAPACKMatrix;
  }
  
  // dealii::LocalIntegrators::Divergence::gradient_matrix(dinfo.matrix(4, false).matrix,
  //  							info.fe_values(0),
  //  							info.fe_values(1),
  //  							minusOne);

  heat::LocalIntegrators::LDG::StiffQFromU(dinfo.matrix(4,false).matrix,
					   info.fe_values(1),
					   info.fe_values(0),
					   minusOne);


  

  // dealii::LocalIntegrators::Divergence::cell_matrix(dinfo.matrix(5,false).matrix,
  // 						    info.fe_values(1),
  // 						    info.fe_values(0),
  // 						    minusOne);

  heat::LocalIntegrators::LDG::StiffUFromQ(dinfo.matrix(5,false).matrix,
					   info.fe_values(0),
					   info.fe_values(1),
					   minusOne);



  heat::LocalIntegrators::LDG::U_MinusRHS_FromSource(referenceDirection,
						     info.fe_values(0),
						     dinfo.vector(0).block(0),
						     minusOne);
  
  
  //TODO This needs to be moved out.
  //dealii::MappingQ1<dim,dim> Mapping;
  
   
  for(unsigned int face_no = 0;
      face_no < dealii::GeometryInfo<dim>::faces_per_cell;
      ++face_no){

    typename dealii::Triangulation<dim>::face_iterator
      face_sys = dinfo.cell->face(face_no);
    if( face_sys->at_boundary() ){      

      dealii::FEFaceValues<dim>
  	fe_v_face_Sys_U(mapping,
			info.finite_element().base_element(0),
			*face_quadrature,
			dealii::update_values			
			| dealii::update_quadrature_points
			| dealii::update_JxW_values
			| dealii::update_normal_vectors);
      
      fe_v_face_Sys_U.reinit(dinfo.cell, face_no);
      

      //face_update_flags);

      dealii::FEFaceValues<dim>
	fe_v_face_Sys_Q(mapping,
			info.finite_element().base_element(1),
			*face_quadrature,
			dealii::update_values
			| dealii::update_quadrature_points
			| dealii::update_normal_vectors
			| dealii::update_JxW_values);

      fe_v_face_Sys_Q.reinit(dinfo.cell, face_no);
      
      
      //face_update_flags);
      

      //auto something= fe_v_face_Sys.get_fe();      
      
      switch( face_sys->boundary_id() ){
	

      case dir_plus:
  	//std::cout << "On dir_plus CELL" << std::endl;
	 
  	break;
      case dir_minus:
  	//std::cout << "On dir_minus CELL" << std::endl;
  	heat::LocalIntegrators::LDG::BoundaryMassU(referenceDirection,
						   dinfo.matrix(9,false),
						   fe_v_face_Sys_U);
	
  	break;
      case neu_plus:
  	//std::cout << "On neu_plus CELL" << std::endl;
	heat::LocalIntegrators::LDG::BoundaryMassQ(referenceDirection,
						   dinfo.matrix(10,false),
						   fe_v_face_Sys_Q);
  	break;
      case neu_minus:
  	//std::cout << "On neu_minus CELL" << std::endl;
  	break;
      default:
  	std::cout << "Something bad happened" << std::endl;
  	assert(false);
  	break;
      }            
    }
  } //End of if(face_sys-at_boundary() )



  { //Work For Boundary Mass U
       
    dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;    
   
    tempLAPACKMatrix.copy_from(dinfo.matrix(9,false).matrix);

    tempLAPACKMatrix.compute_svd();    
    dinfo.matrix(12,false).matrix = *(tempLAPACKMatrix.svd_u);    
    
    for(unsigned int i = 0; i < dinfo.matrix(9,false).matrix.m(); ++i){
      heat::real valueToInsert = 0.0;
      // std::cout << "tempLAPACKMatrix.singular_value( " << i << ") = "
      // 		<<  tempLAPACKMatrix.singular_value(i)  << std::endl;
      if(fabs(tempLAPACKMatrix.singular_value(i) ) > 1e-10){
	valueToInsert = tempLAPACKMatrix.singular_value(i);
      }
      dinfo.matrix(13,false).matrix(i,i) = valueToInsert;      
    }    
  }
    
  
  { //Work For Boundary Mass Q      
    
    dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;    
    
    tempLAPACKMatrix.copy_from(dinfo.matrix(10,false).matrix);

    tempLAPACKMatrix.compute_svd();    
    dinfo.matrix(14,false).matrix = *(tempLAPACKMatrix.svd_u);    
    
    for(unsigned int i = 0; i < dinfo.matrix(10,false).matrix.m(); ++i){
      heat::real valueToInsert = 0.0;
      if(fabs(tempLAPACKMatrix.singular_value(i) ) > 1e-10){
	valueToInsert = tempLAPACKMatrix.singular_value(i);
      }
      dinfo.matrix(15,false).matrix(i,i) = valueToInsert;      
    }    
  }  
}







template <int dim>
void
LDGIntegratorKENT<dim>::face(dealii::MeshWorker::DoFInfo<dim> & dinfoSELF,
			     dealii::MeshWorker::DoFInfo<dim> & dinfoNEIG,
			     dealii::MeshWorker::IntegrationInfo<dim> &infoSELF,
			     dealii::MeshWorker::IntegrationInfo<dim> &infoNEIG ) const {

  
  //  AssertDimension(dinfoSELF.n_matrices(), 16);
  //  AssertDimension(dinfoNEIG.n_matrices(), 16);

  const heat::real One = 1.0;

  heat::LocalIntegrators::LDG::numericalFluxQFromU
    (referenceDirection,
     dinfoSELF.matrix(6,false),
     dinfoSELF.matrix(6,true),
     dinfoNEIG.matrix(6,true),
     dinfoNEIG.matrix(6,false),
     infoSELF.fe_values(1),
     infoNEIG.fe_values(1),
     infoSELF.fe_values(0),
     infoNEIG.fe_values(0),
     One);

  heat::LocalIntegrators::LDG::numericalFluxUFromQ
    (referenceDirection,
     dinfoSELF.matrix(7,false),
     dinfoSELF.matrix(7,true),
     dinfoNEIG.matrix(7,true),
     dinfoNEIG.matrix(7,false),
     infoSELF.fe_values(0),
     infoNEIG.fe_values(0),
     infoSELF.fe_values(1),
     infoNEIG.fe_values(1),
     One);

    // heat::LocalIntegrators::LDG::numericalFluxUFromQ
    // (referenceDirection,
    //  dinfoSELF.matrix(7,false),
    //  dinfoNEIG.matrix(7,true),
    //  dinfoSELF.matrix(7,true),
    //  dinfoNEIG.matrix(7,false),
    //  infoSELF.fe_values(0),
    //  infoNEIG.fe_values(0),
    //  infoSELF.fe_values(1),
    //  infoNEIG.fe_values(1),
    //  One);
  
}
  
template<int dim>
void
LDGIntegratorKENT<dim>::boundary(dealii::MeshWorker::DoFInfo<dim> & dinfo,
				 dealii::MeshWorker::IntegrationInfo<dim> & info) const
{
  
  const dealii::types::boundary_id dir_plus =  1;
  const dealii::types::boundary_id dir_minus = 0;
  const dealii::types::boundary_id neu_plus  = 3;
  const dealii::types::boundary_id neu_minus = 2;

  
  //AssertDimension(dinfo.n_matrices(), 16);

  const heat::real One = 1.0;

  switch( dinfo.face->boundary_id()){
    
  case dir_plus:
    //std::cout << "on dir plus" << std::endl;
    heat::LocalIntegrators::LDG::numericalFluxUFromQBoundary
      (referenceDirection,
       dinfo.matrix(8,false),
       info.fe_values(0),
       info.fe_values(1),
       One
       );

    heat::LocalIntegrators::LDG::Q_MinusRHS_FromDirPlus
      (referenceDirection,
       info.fe_values(1),
       dinfo.vector(0).block(1),
       One);


    break;
  case dir_minus:
    heat::LocalIntegrators::LDG::numericalFluxQFromUBoundary
      (referenceDirection,
       dinfo.matrix(11,false),
       info.fe_values(1),
       info.fe_values(0),
       One);

    heat::LocalIntegrators::LDG::ConstraintUDirMinus(referenceDirection,
						     info.fe_values(0),
						     dinfo.vector(1).block(0),
						     One);

    break;
    
  case neu_plus:
     heat::LocalIntegrators::LDG::numericalFluxUFromQBoundary
      (referenceDirection,
       dinfo.matrix(8,false),
       info.fe_values(0),
       info.fe_values(1),
       One
       );

     heat::LocalIntegrators::LDG::ConstraintQNeuPlus(referenceDirection,
						     info.fe_values(1),
						     dinfo.vector(1).block(1),
						     One);

    //std::cout << "on neu plus" << std::endl;
    // heat::LocalIntegrators::LDG::BoundaryMassQ(referenceDirection,
    // 					       dinfo.matrix(10,false),
    // 					       info.fe_values(1),
    // 					       One);
    break;
  case neu_minus:
    //std::cout << "on nue minus" << std::endl;
    heat::LocalIntegrators::LDG::numericalFluxQFromUBoundary
      (referenceDirection,
       dinfo.matrix(11,false),
       info.fe_values(1),
       info.fe_values(0),
       One);


    heat::LocalIntegrators::LDG::U_MinusRHS_FromNeuMinus(referenceDirection,
							 info.fe_values(0),
							 dinfo.vector(0).block(0),
							 One);
    

    
    break;
  default:
    std::cout << "Something bad happened" << std::endl;
    assert(false);
    break;    
  }
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



} // End Namespace LDGIntegratorKENT
} // End Namespace heat

// template class heat::LDGIntegratorKENT::MassQ<2>;
// template class heat::LDGIntegratorKENT::MassQ<3>;
// template class heat::LDGIntegratorKENT::MassU<2>;
// template class heat::LDGIntegratorKENT::MassU<3>;
// template class heat::LDGIntegratorKENT::StiffUFromQ<2>;
// template class heat::LDGIntegratorKENT::StiffUFromQ<3>;
// template class heat::LDGIntegratorKENT::JumpQ<2>;
// template class heat::LDGIntegratorKENT::JumpQ<3>;
template class heat::LDGIntegratorKENT::LDGIntegratorKENT<2>;
template class heat::LDGIntegratorKENT::LDGIntegratorKENT<3>;
