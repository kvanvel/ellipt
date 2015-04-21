#include "EllipticProblem.hpp"
#include "InitialValuesU.hpp"
#include "DirichletBoundaryValues.hpp"
#include "DirichletBoundaryValuesDot.hpp"
#include "DirichletBoundaryValuesPoisson.hpp"
#include "NeumanBoundaryValuesPoisson.hpp"
#include "NeumannBoundaryValues.hpp"
#include "SourceBodyValues.hpp"
#include "PSigmaPStar.hpp"
#include "Q_UStar_SStar_inverseMassQ_S_U_Q.hpp"
#include "Q_UStar_massU_U_QStar.hpp"
#include "ApproximateSchurComplement.hpp"
#include "SchurComplementPoisson.hpp"
#include "ApproximateSchurComplementPoisson.hpp"
#include "InverseDiffusivity.hpp"

#include "myMatrixTools.hpp"


namespace heat{   

const heat::real pi = boost::math::constants::pi<heat::real>();

template<int dim>
EllipticProblem<dim>::EllipticProblem(const unsigned int degreeIn,
				      const unsigned int refinementsIn)
  :  
  degree(degreeIn),
  refinements(refinementsIn),
  n_eigenvalues(5),
  feSystem
  (
   dealii::FE_DGQ<dim>(degree), 1,    // Density
   dealii::FESystem<dim>(dealii::FE_DGQ<dim>(degree), dim),1 //Q
   ),
  feU(degree),
  feQ(dealii::FE_DGQ<dim>(degree),dim),
  feTraceU_Dir(degree),
  feTraceU_Neu(degree),
  feTraceU_Rob(degree),
  dof_handler(triangulation),
  dof_handlerU(triangulation),
  dof_handlerQ(triangulation),
  dof_handlerTraceU_Dir(triangulation_Dir),
  dof_handlerTraceU_Neu(triangulation_Neu),
  dof_handlerTraceU_Rob(triangulation_Rob),
  UeigenStates(n_eigenvalues),
  QeigenStates(n_eigenvalues)
{ 
}

template<int dim>
void EllipticProblem<dim>::run()
{ 

  //grvy_timer_begin("make_grid");
  make_grid();
  //grvy_timer_end("make_grid");

  //grvy_timer_begin("make_dofs"); 
  make_dofs();
  //grvy_timer_end("make_dofs");
  
  //grvy_timer_begin("report_dofs");
  //report_dofs();  
  //grvy_timer_end("report_dofs");

  //grvy_timer_begin("assemble_system");
  assemble_system();
  //grvy_timer_end("assemble_system");
  
  //DistillMatrices
  //grvy_timer_begin("Distill Matrices");
  heat::myMatrixTools::DistillMatrix(MassU,MassUSparsityPattern);
  heat::myMatrixTools::DistillMatrix(InverseMassU,InverseMassUSparsityPattern);
  heat::myMatrixTools::DistillMatrix(TotalUFromQ,TotalUFromQSparsityPattern);  
  heat::myMatrixTools::DistillMatrix(MassQ, MassQSparsityPattern);
  heat::myMatrixTools::DistillMatrix(InverseMassQ,
				     InverseMassQSparsityPattern);
  heat::myMatrixTools::DistillMatrix(TotalQFromU,
				     TotalQFromUSparsityPattern);  
  //grvy_timer_end("Distill Matrices");

  

  //grvy_timer_begin("ElliptSolve");
  ElliptSolve();
  
  //grvy_timer_end("ElliptSolve");
  //print_system();
  PrintState(0);
  
  //EigenSolve();

  //output_results_eigensystem();
  
  {
    dealii::Vector<heat::real> difference(triangulation.n_active_cells() );
    //const heat::DirichletBoundaryValues<dim> dirBC;
    
    dealii::VectorTools::integrate_difference
      (dealii::MappingQ1<dim>(),
       dof_handlerU,
       Ustate,
       heat::DirichletBoundaryValues<dim>(),
       difference,
       dealii::QGauss<dim>(degree+4),
       dealii::VectorTools::NormType::L2_norm);

    //std::cout << "difference.l2_norm() = " << difference.l2_norm() << std::endl;

    this->L2errorU = difference.l2_norm();
   
    //std::cout << "UState = " << Ustate << std::endl;    
  }
  {
    dealii::Vector<heat::real> differenceQ(triangulation.n_active_cells() );
    const heat::NeumannBoundaryValues<dim> neuBC;
    dealii::VectorFunctionFromTensorFunction<dim> neuBC_vector_function(neuBC,0,dim);
    
    dealii::VectorTools::integrate_difference
      (
       dealii::MappingQ1<dim>(),
       dof_handlerQ,
       Qstate,
       neuBC_vector_function,
       differenceQ,
       dealii::QGauss<dim>(degree+4),
       dealii::VectorTools::NormType::L2_norm
       );
    
    //std::cout << "differenceQ.l2_norm() = " << differenceQ.l2_norm() << std::endl;
    this->L2errorQ = differenceQ.l2_norm();
    
  }
  //Make initial condition  


  
  //makeU_RHS_FromInitialCondition();
  //makeConstraintDirU(0.0);
  
  //USolve();
  // //This next line is here because USolve actually puts its solution
  // //in UStateDot;
  // //TODO:  This is a hack and should be fixed.

  // Taking care of the initial condition is extremely confusing.
  // TODO:  Fix this.
  //Ustate = UStateDot;
  //UStateDot = 0;
  
  //TODO:  This is assuming that the start time is 0,
  // Something that may not be true.
  //makeConstraintNeuQ(0.0);
  //makeQ_RHS_FromU();
  
  //QSolve_NEW(Qstate,lagrangeQ,Q_RHS, constraintNeuQ);
  
  //makeConstraintDirU(0.0);

  //makeU_RHS_FromQ();
  //makeU_RHS_From_U();
  //U_RHS += U_RHS_From_U;  
  //USolve();

  // {
  //   dealii::Vector<heat::real> tempUStateDot(UStateDot.size() );
  //   USolve_NEW(tempUStateDot,lagrangeU,U_RHS,constraintDirU);
  //   tempUStateDot -= UStateDot;
  //   std::cout << "tempUStateDot.norm_sqr() = " << tempUStateDot.norm_sqr()
  // 	      << std::endl;
  // }
  
  // Print Initial condition
  
  //makeElec_RHS_From_DirichletBoundary(0.0);
  //makeConstraint_RHS_From_NeumanBoundary();  
  //PoissonSolve();
  
  //PrintState(0);
  //simulate();
  
}

//23456789012345678901345678901234567890123456789012345678901234567890123456789

template<int dim>
void
EllipticProblem<dim>::make_grid()
{
  //Initialize ReferenceDirection;
  
  for(unsigned int i = 0; i < dim; ++i){
    referenceDirection[i] = 1.0;
  }
  
  //const unsigned int n_refinements = 2;
  dealii::Point<dim> 
    bottom_left,
    top_right;


  if(2==dim){
    bottom_left = dealii::Point<dim>( (-1.0/2.0)*pi, (-1.0/2.0)*pi);
    top_right = dealii::Point<dim>( pi, pi);
  } else if(3 == dim){
    bottom_left
      = dealii::Point<dim>( (-1.0/2.0)*pi, (-1.0/2.0)*pi,(-1.0/2.0)*pi);
    top_right = dealii::Point<dim>( pi, pi, pi);    
  } else {
    assert(false); // Not implemented  //TODO: Change to dealii assert
  }

  std::vector<unsigned int> repetitions;
  for(unsigned int i = 0; i < dim; ++i){
    repetitions.emplace_back(1);
  }  
  
  //const bool colorize = true;
  const dealii::types::boundary_id dir_plus = 3;
  const dealii::types::boundary_id dir_minus = 5;
  // dealii::GridGenerator::subdivided_hyper_rectangle
  //   (triangulation, repetitions,bottom_left,top_right, colorize);
  
  dealii::GridGenerator::hyper_L(triangulation,-1,1);

  //We have to do some gymnastics to assign the boundary indicator based
  //on the normal vector.  We do that here.
  {
    dealii::QGauss<dim - 1> face_quadrature_formula(1);
    const dealii::UpdateFlags face_update_flags = 
      dealii::update_values
      | dealii::update_quadrature_points
      | dealii::update_normal_vectors;
  
    for(auto cell_sys = triangulation.begin_active();
	cell_sys != triangulation.end();
	++cell_sys){
      //grvy_timer_begin("loop over faces");
      for(unsigned int face_no = 0; 
	  face_no < dealii::GeometryInfo<dim>::faces_per_cell;
	  ++face_no){

	auto face_sys = cell_sys->face(face_no);

	if( face_sys->at_boundary() ){
	  //determine what boundary and assign the correct boundary indicator.
	  dealii::FEFaceValues<dim>
	    fe_v_face_U(mapping,
			feU,
			face_quadrature_formula,
			face_update_flags);
	  fe_v_face_U.reinit(cell_sys,face_no);
	  const auto normal = fe_v_face_U.normal_vector(0);
	  if(normal * referenceDirection > 0){
	    face_sys->set_boundary_indicator(dir_plus);
	  } else {
	    face_sys->set_boundary_indicator(dir_minus);
	  }	  
	}	
      }
    }	
  }

  
  
  // dealii::GridGenerator::hyper_ball
  //   (triangulation, p1, /*radius = */ 1);

  // static const dealii::HyperBallBoundary<dim> boundary;
  // triangulation.set_boundary(0,boundary);
  
  triangulation.refine_global(refinements);
  

  //Interior Refinement
  unsigned int nInteriorRefinements =0;
  for(unsigned int i = 0; i < nInteriorRefinements; ++i){
    const dealii::Point<dim> center = (bottom_left+top_right)* 0.5;
    for(auto cellT = triangulation.begin_active();
	cellT != triangulation.end();
	++cellT){
      for(unsigned int v = 0;
	  v < dealii::GeometryInfo<dim>::vertices_per_cell;
	  ++v){
	const double DistanceFromCenter = center.distance(cellT->vertex(v));
	if(std::fabs(DistanceFromCenter) < 0.25){
	  cellT->set_refine_flag();
	}
      }
    }
    triangulation.execute_coarsening_and_refinement();
  }

  //Border Refinement
  unsigned int nBorderRefinements = 0;
  for(unsigned int i = 0; i < nBorderRefinements; ++i){
    for(auto cell = triangulation.begin_active();
	cell != triangulation.end();
	++cell){
      for(unsigned int face_no = 0; 
	  face_no < dealii::GeometryInfo<dim>::faces_per_cell;
	  ++face_no){
	auto face = cell->face(face_no);	
	if(face->at_boundary() ){
	  cell->set_refine_flag();
	}
      }      
    }     
    triangulation.execute_coarsening_and_refinement();
  }
  
  
  //dealii::GridTools::distort_random(0.1,triangulation,true);
  
  {
  std::ofstream out ("grid-1.eps");
  dealii::GridOut grid_out;
  grid_out.write_eps(triangulation,out);
  //std::cout << "Grid written to grid-1.eps" << std::endl;
  }
  
  
  // Dirichelet Boundaries
  boundary_ids_Dir_plus.insert(dir_plus);
  boundary_ids_Dir_minus.insert(dir_minus);

  
  
  //boundary_ids_Dir_minus.insert(2);
  //boundary_ids_Dir_plus.insert(3);
  
  //boundary_ids_Dir.insert(2);
  //boundary_ids_Dir.insert(3); 
   if(3 == dim){
     boundary_ids_Dir_minus.insert(4);
     boundary_ids_Dir_plus.insert(5);
   }
  
   
   //Neuman Boundaries   
   //boundary_ids_Neu.insert(0);
   //boundary_ids_Neu.insert(1);
   //boundary_ids_Neu_minus.insert(2);
   //boundary_ids_Neu_plus.insert(3);   

   // if(3 == dim){
  //   boundary_ids_Neu_minus.insert(4);
  //   boundary_ids_Neu_plus.insert(5);
  // }
   

  //std::set<dealii::types::boundary_id> boundary_id_Rob;
  //boundary_id_Rob.insert(2);
  //boundary_ids_Rob.insert(1);
  
  // Robin Boundaries
  //boundary_ids_Rob.insert();
  //boundary_ids_Rob.insert(3);
  //assert(boundary_ids_Rob.empty() );
}

template<int dim>
void
EllipticProblem<dim>::make_dofs()
{
  //This is pretty tricky
  // https://groups.google.com/forum/#!searchin/dealii/2-d$203-d/dealii/0cwQcfGc2zQ/Pc2G3zIoH-YJ

  // There is a good reason why we use the dof_handlerU for the trace spaces
  // of Q.  It has to do with the fact that we will be dotting with the unit
  // normal vector.
   
  dof_handler.distribute_dofs(feSystem);
  dof_handlerU.distribute_dofs(feU);
  dof_handlerQ.distribute_dofs(feQ);

  if( !boundary_ids_Dir_minus.empty() ){
    DirSurfaceToVolumeTriangulationMap = 
      dealii::
      GridTools::
      extract_boundary_mesh(dof_handlerU, 
			    dof_handlerTraceU_Dir,
			    boundary_ids_Dir_minus);    
    
    for(auto it = DirSurfaceToVolumeTriangulationMap.begin();
	it != DirSurfaceToVolumeTriangulationMap.end(); 
	++it){    
      DirVolumeToSurfaceTriangulationMap.emplace(std::make_pair(it->second, it->first));
      
    }
    dof_handlerTraceU_Dir.distribute_dofs(feTraceU_Dir);
  }
 

  if ( !boundary_ids_Neu_plus.empty() ){
    this->NeuSurfaceToVolumeTriangulationMap = 
      dealii::GridTools::extract_boundary_mesh(dof_handlerU,
					       dof_handlerTraceU_Neu,
					       boundary_ids_Neu_plus);  
  
    for(auto it = NeuSurfaceToVolumeTriangulationMap.begin(); 
	it != NeuSurfaceToVolumeTriangulationMap.end(); 
	++it){
      NeuVolumeToSurfaceTriangulationMap.emplace(std::make_pair(it->second, it->first) );
    }
    dof_handlerTraceU_Neu.distribute_dofs(feTraceU_Neu);
  }

  //TODO:  Robin stuff has never been tested.
  if( !boundary_ids_Rob_plus.empty() ){
    RobSurfaceToVolumeTriangulationMap = 
      dealii::GridTools::extract_boundary_mesh(dof_handlerU,
					       dof_handlerTraceU_Rob,
					       boundary_ids_Rob_plus);
    
    for(auto it = RobSurfaceToVolumeTriangulationMap.begin();      
	it != RobSurfaceToVolumeTriangulationMap.end(); 
	++it){
      RobVolumeToSurfaceTriangulationMap.emplace(std::make_pair(it->second, it->first) );
    }
    dof_handlerTraceU_Rob.distribute_dofs(feTraceU_Rob);
  }
  
  

  
  

  dealii::DoFRenumbering::block_wise(dof_handler);
  std::vector<dealii::types::global_dof_index> dofs_per_block(2);
  dealii::DoFTools::count_dofs_per_block(dof_handler,dofs_per_block);

  // dealii::DoFTools::make_sparsity_pattern(dof_handlerElec,
  // 					  c_sparsityElec,
  // 					  constraintElec,
  // 					  /* keep_constrained_dofs = */ false);
  
  const unsigned int 
    n_u = dofs_per_block[0],
    n_q = dofs_per_block[1];
    
    
  const unsigned int
    n_U = dof_handlerU.n_dofs(),
    n_Q = dof_handlerQ.n_dofs(),
    n_TraceU_Dir = dof_handlerTraceU_Dir.n_dofs(),
    n_TraceU_Neu = dof_handlerTraceU_Neu.n_dofs();
    //n_TraceU_Rob = dof_handlerTraceU_Rob.n_dofs();

  //TODO: Make this a dealii style assert.
  assert(n_U == n_u);   
  assert(n_Q == n_q);


  dealii::
    Table<2,dealii::DoFTools::Coupling>
    interiorMask(feSystem.n_components(),feSystem.n_components());

  dealii::
    Table<2,dealii::DoFTools::Coupling> 
    fluxMask(feSystem.n_components(), feSystem.n_components());

  for(unsigned int c = 0; c< feSystem.n_components(); ++c){
    for(unsigned int d = 0; d < feSystem.n_components(); ++d){
      interiorMask[c][d]= dealii::DoFTools::none;
      fluxMask[c][d]= dealii::DoFTools::none;
    }
  }
  
  // U and U itself
  interiorMask[0][0] = dealii::DoFTools::always; 
  
  // U and Q
  for(unsigned int j = 1; j < 1 + dim; ++j){
    interiorMask[0][j] = dealii::DoFTools::always;
    fluxMask[0][j] = dealii::DoFTools::nonzero;
  }
   
  //Q and U;
  for(unsigned int i = 1; i < 1 + dim; ++i){
    interiorMask[i][0] = dealii::DoFTools::always;
    fluxMask[i][0] = dealii::DoFTools::nonzero;
  }
  
  //Q and Q;   
  for(unsigned int i = 1; i < 1 + dim; ++i){
    for(unsigned int j = 1; j < 1 + dim; ++j){
      if(i == j){
	interiorMask[i][j] = dealii::DoFTools::always;
      }
    }
  }

   
  
  //dealii::DoFTools::count_dofs_per_block(dof_handler,
  dealii::BlockIndices bi(dofs_per_block);
  dealii::BlockCompressedSparsityPattern csp(bi,bi);
  
   
  sparsity_pattern.reinit(2,2);
   
  const unsigned int
    n_u_dpc = feU.n_dofs_per_cell(),
    n_q_dpc = feQ.n_dofs_per_cell();

  sparsity_pattern.block(0,0).reinit(n_u,n_u,n_u_dpc);
  sparsity_pattern.block(0,1).reinit(n_u,n_q, 2 * (2 * dim + 1) * n_q_dpc);

   
  sparsity_pattern.block(1,0).reinit(n_q,n_u,2 * (2 * dim + 1) * n_u_dpc);
  sparsity_pattern.block(1,1).reinit(n_q,n_q, n_q_dpc);  

    

  // std::cout << "sparsity_pattern.memory_consumption() = "
  // 	    <<  sparsity_pattern.memory_consumption() << std::endl;
  
  sparsity_pattern.collect_sizes();
      
  dealii::
    DoFTools::
    make_flux_sparsity_pattern(dof_handler,
			       sparsity_pattern,
			       interiorMask,
			       fluxMask);  
   
  
  // std::cout << "After make flux sparsity pattern" << std::endl;
  // std::cout << "sparsity_pattern.memory_consumption() = "
  // 	    <<  sparsity_pattern.memory_consumption() << std::endl;
  
  sparsity_pattern.compress();   
  // std::cout << "After compress" << std::endl;
  // std::cout << "sparsity_pattern.memory_consumption() = "
  // 	    <<  sparsity_pattern.memory_consumption() << std::endl;
  
  system_matrix.reinit(sparsity_pattern);
  
  state.reinit(2);
  state.block(0).reinit(n_u);
  state.block(1).reinit(n_q);
  
  state.collect_sizes();

  Ustate = state.block(0);
  lambdaU = state.block(0);
  for(unsigned int i = 0; i < n_eigenvalues; ++i){
    (UeigenStates[i]).reinit(n_u);
    (QeigenStates[i]).reinit(n_q);
  }
  
  Qstate = state.block(1);
  //U_RHS = state.block(0);
  U_MinusRHS = state.block(0);
  U_RHS_From_U = state.block(0);
  UdirConstraint_RHS = state.block(0);
  //Q_RHS = state.block(1);
  Q_MinusRHS = state.block(1);

  
  lagrangeU.reinit(n_TraceU_Dir);
  lagrangeQ.reinit(n_TraceU_Neu);
  
  constraintDirU.reinit(n_TraceU_Dir);
  constraintNeuQ.reinit(n_TraceU_Neu);

}


// template<int dim>
// void
// EllipticProblem<dim>::makeU_RHS_From_U()
// {
//   U_RHS_From_U = 0;
//   //TODO:  Check that these are appropriate quadrature formulae.
//   dealii::QGauss<dim> body_quadrature_formula(2* degree+2);
//   dealii::QGauss<dim - 1> face_quadrature_formula(2 * degree+2);  
  
//   const dealii::UpdateFlags 
//     body_update_flags = 
//     dealii::update_values
//     | dealii::update_gradients
//     | dealii::update_quadrature_points
//     | dealii::update_JxW_values,
    
//     face_update_flags = 
//     dealii::update_values
//     | dealii::update_quadrature_points
//     | dealii::update_JxW_values
//     | dealii::update_normal_vectors;

//   const unsigned int
//     dofs_per_cell_U = dof_handlerU.get_fe().dofs_per_cell,
//     dofs_per_cell_Elec = dof_handlerElec.get_fe().dofs_per_cell;
  
//   std::vector<dealii::types::global_dof_index> 
//     dofs_self_U(dofs_per_cell_U),
//     dofs_neig_U(dofs_per_cell_U),
//     dofs_self_Elec(dofs_per_cell_Elec),
//     dofs_neig_Elec(dofs_per_cell_Elec);
  
//   dealii::FEValues<dim>
//     fe_v_body_U(mapping, feU,  body_quadrature_formula,
// 		body_update_flags);
  
//   dealii::FEFaceValues<dim>    
//     fe_v_face_U(mapping,feU,face_quadrature_formula,
// 		face_update_flags);
  
//   dealii::FESubfaceValues<dim> 
//     fe_v_subface_U(mapping, feU, face_quadrature_formula,
// 		   face_update_flags);
  
//   dealii::FEFaceValues<dim> 
//     fe_v_face_neig_U(mapping,
// 		     feU,
// 		     face_quadrature_formula,
// 		     face_update_flags);
  
//   dealii::FESubfaceValues<dim>
//     fe_v_subface_neig_U(mapping,
// 			feU,
// 			face_quadrature_formula,
// 			face_update_flags);
  
//   dealii::FEValues<dim>
//     fe_v_body_Elec(mapping, 
// 		   feElec,
// 		   body_quadrature_formula,
// 		   body_update_flags);

//   dealii::FEFaceValues<dim>
//     fe_v_face_Elec(mapping, 
// 		   feElec,
// 		   face_quadrature_formula, 
// 		   face_update_flags);

//   dealii::FESubfaceValues<dim> 
//     fe_v_subface_Elec(mapping, feElec, face_quadrature_formula,
// 		      face_update_flags);
  
//   dealii::FEFaceValues<dim> 
//     fe_v_face_neig_Elec(mapping,
// 			feElec,
// 			face_quadrature_formula,
// 			face_update_flags);
  
//   dealii::FESubfaceValues<dim>
//     fe_v_subface_neig_Elec(mapping,
// 			   feElec,
// 			   face_quadrature_formula,
// 			   face_update_flags);

//   dealii::Vector<heat::real>
//     U_Vector(dofs_per_cell_U),
//     U_Vector_ForFaces(dofs_per_cell_U);    
  
//   typename dealii::DoFHandler<dim>::active_cell_iterator
//     cell_sys = dof_handler.begin_active(),
//     endc_sys = dof_handler.end(),
//     cell_U = dof_handlerU.begin_active(),
//     cell_Q = dof_handlerQ. begin_active(),//  Don't think we need the Q.
//     cell_Elec = dof_handlerElec.begin_active();
  
//   for(; 
//       cell_sys != endc_sys;
//       ++cell_sys, ++cell_U, ++cell_Q, ++cell_Elec){
    
//     cell_U->get_dof_indices(dofs_self_U);    
//     fe_v_body_U.reinit(cell_U);
//     fe_v_body_Elec.reinit(cell_Elec);    
//     // Compute body term
//     U_Vector = 0;
//     assembleStiffUFromULocal(fe_v_body_U,fe_v_body_Elec,U_Vector);
//     //Insert body term
//     U_RHS_From_U.add(dofs_self_U,U_Vector);
//     U_Vector = 0;
    
//     // Now compute the face terms
//     for(unsigned int face_no = 0;
// 	face_no < dealii::GeometryInfo<dim>::faces_per_cell;
// 	++face_no){
      
//       typename dealii::DoFHandler<dim>::face_iterator
// 	face_sys = cell_sys->face(face_no);      
//       if(face_sys->at_boundary() ){
// 	continue;  
// 	// We only do this for the faces that are in the interior.
//       }
//       else {
// 	fe_v_face_U.reinit(cell_U,face_no);
// 	fe_v_face_Elec.reinit(cell_Elec,face_no);
	
// 	Assert( cell_sys-> neighbor(face_no).state() 
// 		== dealii::IteratorState::valid,
// 		dealii::ExcInternalError() );
	
// 	if(face_sys->has_children() ){
// 	  Assert(false, dealii::ExcNotImplemented() );	  
// 	  // Implemented but not tested.
// 	  // We first must worry about the mesh refinement and the electric
// 	  // field.
	  
// 	  const unsigned int neighbor2 = 
// 	    cell_sys->neighbor_face_no(face_no);
// 	  for( unsigned int subface_no = 0;
// 	       subface_no < face_sys->number_of_children();
// 	       ++subface_no){
// 	    typename dealii::DoFHandler<dim>::cell_iterator
// 	      neighbor_child_U
// 	      = cell_U->neighbor_child_on_subface(face_no,subface_no),
// 	      neighbor_child_Elec
// 	      = cell_Elec->neighbor_child_on_subface(face_no,subface_no);
	    
// 	    Assert(!(neighbor_child_U->has_children())
// 		   , dealii::ExcInternalError() );
	    
// 	    fe_v_subface_U.reinit(cell_U,face_no,subface_no);
// 	    fe_v_subface_Elec.reinit(cell_Elec,face_no,subface_no);
// 	    fe_v_face_neig_U.reinit(neighbor_child_U,neighbor2);
// 	    fe_v_face_neig_Elec.reinit(neighbor_child_Elec,neighbor2);
	    
// 	    U_Vector = 0;
// 	    assembleFluxUFromULocal(fe_v_subface_U,
// 				    fe_v_subface_Elec,
// 				    fe_v_face_neig_U,
// 				    fe_v_face_neig_Elec,
// 				    U_Vector);
// 	    U_Vector *= 1.0;
// 	    U_RHS_From_U.add(dofs_self_U,U_Vector);
// 	    U_Vector = 0;
// 	  }
// 	}
	
// 	else if( !(cell_sys->neighbor_is_coarser( face_no ) ) ){	  
// 	  typename dealii::DoFHandler<dim>::cell_iterator
// 	    neighbor_U = cell_U->neighbor(face_no),
// 	    neighbor_Elec = cell_Elec->neighbor(face_no);
		
// 	  // Now we know the neighbors face is neither finer nor coarser
// 	  const unsigned int neighbor2
// 	    = cell_sys->neighbor_of_neighbor(face_no);
// 	  // TODO:  Put the assert in here that checks the functionality 
// 	  // of the previous assignment
	  
// 	  fe_v_face_U.reinit(cell_U,face_no);
// 	  fe_v_face_Elec.reinit(cell_Elec,face_no);
// 	  fe_v_face_neig_U.reinit(neighbor_U, neighbor2);
// 	  fe_v_face_neig_Elec.reinit(neighbor_Elec, neighbor2);
	  
// 	  cell_U->get_dof_indices(dofs_self_U);
// 	  neighbor_U->get_dof_indices(dofs_neig_U);
// 	  neighbor_Elec->get_dof_indices(dofs_neig_Elec);
// 	  U_Vector = 0;
// 	  assembleFluxUFromULocal(fe_v_face_U,
// 				  fe_v_face_Elec,
// 				  fe_v_face_neig_U,
// 				  fe_v_face_neig_Elec,
// 				  U_Vector);
// 	  U_Vector *= -1.0;
// 	  U_RHS_From_U.add(dofs_self_U,U_Vector);
// 	  U_Vector = 0;
// 	}
// 	else{
// 	  Assert(false, dealii::ExcNotImplemented() );
// 	  // Implemented but not tested;
	  
// 	  const std::pair<unsigned int, unsigned int> neighbor_face_subface
// 	    = cell_sys->neighbor_of_coarser_neighbor(face_no);
// 	  const unsigned int nFace_no = neighbor_face_subface.first;
// 	  const unsigned int nSubface_no = neighbor_face_subface.second;

// 	  assert
// 	    (cell_sys
// 	     ->neighbor(face_no)
// 	     ->face(nFace_no)
// 	     ->child(nSubface_no) 
// 	     == cell_sys->face(face_no)
// 	     );

// 	  typename dealii::DoFHandler<dim>::cell_iterator
// 	    neighbor_U
// 	    = cell_U->neighbor(face_no),
// 	    neighbor_Elec
// 	    = cell_Elec->neighbor(face_no);

// 	  fe_v_face_U.reinit(cell_U,face_no);
// 	  fe_v_face_Elec.reinit(cell_Elec,face_no);
// 	  fe_v_subface_neig_U.reinit(neighbor_U,nFace_no,nSubface_no);
// 	  fe_v_subface_neig_Elec.reinit(neighbor_Elec, nFace_no, nSubface_no);

// 	  U_Vector = 0;
// 	  assembleFluxUFromULocal(fe_v_face_U,
// 				  fe_v_face_Elec,
// 				  fe_v_subface_neig_U,
// 				  fe_v_subface_neig_Elec,
// 				  U_Vector);
// 	  U_Vector *= -1.0;
// 	  U_RHS_From_U.add(dofs_self_U,U_Vector);
// 	  U_Vector = 0;
// 	}
//       }
//     }
//   }
// }

/*
template<int dim>
void
EllipticProblem<dim>::makeElec_RHS_From_DirichletBoundary(const heat::real time)
{
  heat::DirichletBoundaryValuesPoisson<dim> dirBC;
  dirBC.set_time(time);
  // These are boundary values of the Potentional.
  
  dealii::Vector<heat::real> temp(TraceElecDir.m());
  dealii::QGauss<dim-1> face_quadrature_formula(degree + 2);
  dealii
    ::VectorTools
    ::create_right_hand_side
    (faceMapping,
     dof_handlerTraceU_Dir,
     face_quadrature_formula,
     dirBC,
     temp
     );
  
  //TODO: Somthing might be wrong here.  I would have expected that we 
  // would multiply by// TraceElecDir and not TraceElecDir Transposed.
  TraceElecDir.Tvmult(Elec_RHS, temp);
}
*/

/*
template<int dim>
void
EllipticProblem<dim>::makeConstraint_RHS_From_NeumanBoundary()
{
heat::NeumanBoundaryValuesPoisson<dim> neuBC;
  //dealii::Vector<heat::real> temp(TraceElecNeu.m() );
  dealii::QGauss<dim-1> face_quadrature_formula(degree + 2);
  dealii
    ::VectorTools
    ::create_right_hand_side(faceMapping,
			     dof_handlerTraceU_Neu,
			     face_quadrature_formula,
			     neuBC,
			     PoissonNeuBC_RHS);
}

*/

/*
template <int dim>
void 
EllipticProblem<dim>::makeU_RHS_FromInitialCondition()
{
  heat::InitialValuesU<dim> initU;
  dealii::QGauss<dim> quadrature_formula(degree + 2);
  dealii
    ::VectorTools
    ::create_right_hand_side(mapping,
			     dof_handlerU,
			     quadrature_formula,
			     initU,
			     U_RHS);
}

*/

/*
template <int dim>
void
EllipticProblem<dim>::makeQ_RHS_FromU()
{
  TotalQFromU.vmult(Q_RHS, Ustate);
}
*/

/*
template <int dim>
void
EllipticProblem<dim>::makeU_RHS_FromQ() 
{
  TotalUFromQ.vmult(U_RHS,Qstate);
}
*/


template <int dim>
void EllipticProblem<dim>::report_dofs()
{
  std::cout 
    << "Number of active cells: "
    << triangulation.n_active_cells() << std::endl
    << "Total number of cells: " 
    << triangulation.n_cells() << std::endl
    << "Number of degrees of freedom: "
    << dof_handler.n_dofs() << std::endl
    << "U Number of degrees of freedom: " 
    << dof_handlerU.n_dofs() << std::endl
    << "Q Number of degrees of freedom: " 
    << dof_handlerQ.n_dofs() << std::endl

    << dof_handler.n_dofs() << " = " 
    << dof_handlerU.n_dofs() << " + " << dof_handlerQ.n_dofs() 
    << std::endl;

  std::cout << std::endl
	    << "Dirichelet Boundary"
	    << std::endl
	    << "Number of active cells:"
	    << triangulation_Dir.n_active_cells() << std::endl
	    << "Total number of cells: " 
	    << triangulation_Dir.n_cells() << std::endl
	    << "Number of degrees of freedom: "
	    << dof_handlerTraceU_Dir.n_dofs()
	    << std::endl;
  
  std::cout << std::endl
	    << "Neuman Boundary" << std::endl
	    << "Number of active cells:"
	    << triangulation_Neu.n_active_cells() << std::endl
	    << "Total number of cells: "
	    << triangulation_Neu.n_cells() << std::endl
	    << "Number of degrees of freedom: "
	    << dof_handlerTraceU_Neu.n_dofs()
	    << std::endl;

  std::cout << std::endl
	    << "Robin Boundary" << std::endl
	    << "Number of active cells:"
	    << triangulation_Rob.n_active_cells() << std::endl
	    << "Total number of cells: "
	    << triangulation_Rob.n_cells() << std::endl
	    << "Number of degrees of freedom: "
	    << dof_handlerTraceU_Rob.n_dofs()
	    << std::endl;  
}

/*
template <int dim>
void
EllipticProblem<dim>::makeConstraintDirU(const heat::real time)
{
heat::DirichletBoundaryValuesDot<dim> DirDot;
  DirDot.set_time(time);
  dealii::QGauss<dim - 1> face_quadrature_formula(degree+2);
  dealii
    ::VectorTools
    ::create_right_hand_side(faceMapping,
			     dof_handlerTraceU_Dir,
			     face_quadrature_formula,
			     DirDot,
			     constraintDirU);
}
*/

/*
template <int dim>
void
EllipticProblem<dim>::makeConstraintNeuQ(const heat::real time)
{
  heat::DirichletBoundaryValues<dim> Neu;
  Neu.set_time(time);

  if(dof_handlerTraceU_Neu.has_active_dofs() ){  
    dealii::QGauss<dim - 1> face_quadrature_formula(degree+2);
    dealii
      ::VectorTools
      ::create_right_hand_side(faceMapping,
			     dof_handlerTraceU_Neu,
			     face_quadrature_formula,
			     Neu,
			     constraintNeuQ);
  }
}
*/

template <int dim>
void
EllipticProblem<dim>::assemble_system()
{

  //grvy_timer_begin("assemble_system_inside");
  heat::real tolerance = 1e-10;
  dealii::QGauss<dim> body_quadrature_formula(degree+4);
  dealii::QGauss<dim - 1> face_quadrature_formula(degree+4);
  
  dealii
    ::FEValues<dim> 
    fe_values(feSystem,body_quadrature_formula,
	      dealii::update_values 
	      | dealii::update_gradients
	      | dealii::update_quadrature_points
	      | dealii::update_JxW_values);  
  
  dealii
    ::FEValues<dim> 
    fe_face_values(feSystem,body_quadrature_formula,
		   dealii::update_values 
		   | dealii::update_gradients
		   | dealii::update_quadrature_points
		   | dealii::update_JxW_values);

  const unsigned int n_q_points = body_quadrature_formula.size();    
  const InverseDiffusivity<dim> DiffusivityInverse;
    
  std::vector<dealii::Tensor<2,dim> > DiffusivityInverseValues(n_q_points);
  
  const dealii::FEValuesExtractors::Scalar density(0);
  const dealii::FEValuesExtractors::Vector FLUX(1);
 
  
  // For now, just build everything separately.  Don't worry about saving
  
  dealii::DoFTools::map_dof_to_boundary_indices(dof_handlerQ,
						boundary_ids_Dir_minus,
						uDofsWithSupportOnBoundary);
  
						
  /*std::cout << "uDofsWithSupportOnBoundary.size() = " << uDofsWithSupportOnBoundary.size() << std::endl;

  for(unsigned int i = 0;
      i < uDofsWithSupportOnBoundary.size();
      ++i){
    const dealii::types::global_dof_index j = uDofsWithSupportOnBoundary[i];
    if(j != dealii::numbers::invalid_dof_index){
      std::cout
	<< "uDofsWithSupportOnBoundary["<< i << "]= "
	<< uDofsWithSupportOnBoundary[i] << std::endl;
    }
  }
  */
  
  //  dealii
  //     ::DoFTools
  // ::make_boundary_sparsity_pattern
  // (dof_handlerU,
  //  uDofsWithSupportOnBoundary,
  //  TRACE_U_DIR_NEW_SparsityPattern);

  //   typename dealii::FunctionMap<dim>::type dirichlet_boundary_functions;
  //   dealii::ZeroFunction<dim> homogeneous_dir_bc(1);
  //   dirichlet_boundary_functions[0] = &homogeneous_dir_bc;

    
  //   dealii::
  //     MatrixCreator::
  //     create_boundary_mass_matrix
  //     (mapping,
  //      dof_handlerU,
  //      face_quadrature_formula,
  //      TRACE_U_DIR_NEW,
  //      dirichlet_boundary_functions,
  //      TRACE_RHS_VECTOR,
  //      uDofsWithSupportOnBoundary       
  //      ); 
  
  
  // We are following the step-30 tutorial 
  const unsigned int n_dofs_U = dof_handlerU.n_dofs();
  const unsigned int n_dofs_Q = dof_handlerQ.n_dofs();
  const unsigned int n_dofs_TrU_Dir = dof_handlerTraceU_Dir.n_dofs();
  const unsigned int n_dofs_TrU_Neu = dof_handlerTraceU_Neu.n_dofs();  
  
  
  MassU.reinit(sparsity_pattern.block(0,0) );
  Svd_U.reinit(sparsity_pattern.block(0,0) );
  Svd_Sigma.reinit(sparsity_pattern.block(0,0) );
  //Svd_Sigma_Inverse.reinit(sparsity_pattern.block(0,0));
  InverseMassU.reinit(sparsity_pattern.block(0,0));
  BoundaryMassU.reinit(sparsity_pattern.block(0,0) );
  
  StiffUFromQ.reinit(sparsity_pattern.block(0,1));
  FluxBoundaryUFromQ.reinit(sparsity_pattern.block(0,1));
  FluxCenterUFromQ.reinit(sparsity_pattern.block(0,1));
  FluxPosUFromQ.reinit(sparsity_pattern.block(0,1));
  FluxNegUFromQ.reinit(sparsity_pattern.block(0,1));
  TotalUFromQ.reinit(sparsity_pattern.block(0,1));

  MassQ.reinit(sparsity_pattern.block(1,1) );
  InverseMassQ.reinit(sparsity_pattern.block(1,1));
  
  StiffQFromU.reinit(sparsity_pattern.block(1,0));  
  FluxBoundaryQFromU.reinit(sparsity_pattern.block(1,0));
  FluxCenterQFromU.reinit(sparsity_pattern.block(1,0));
  FluxPosQFromU.reinit(sparsity_pattern.block(1,0));
  FluxNegQFromU.reinit(sparsity_pattern.block(1,0));
  TotalQFromU.reinit(sparsity_pattern.block(1,0));
  
  
  TraceUDir.reinit(n_dofs_TrU_Dir,n_dofs_U,degree);
  TraceQNeu.reinit(n_dofs_TrU_Neu,n_dofs_Q,degree); 

  
  const unsigned int
    dofs_per_cell_U = dof_handlerU.get_fe().dofs_per_cell,
    dofs_per_cell_Q = dof_handlerQ.get_fe().dofs_per_cell,
    dofs_per_cell_TraceU_Dir = feTraceU_Dir.n_dofs_per_cell(),  
    dofs_per_cell_TraceU_Neu = feTraceU_Neu.n_dofs_per_cell();
  
  
  std::vector<dealii::types::global_dof_index> 
    dofs_self_U(dofs_per_cell_U),
    dofs_neig_U(dofs_per_cell_U),
    dofs_self_Q(dofs_per_cell_Q),
    dofs_neig_Q(dofs_per_cell_Q),
    dofs_self_TraceU_Dir(dofs_per_cell_TraceU_Dir),
    dofs_self_TraceU_Neu(dofs_per_cell_TraceU_Neu);

  
  const dealii::UpdateFlags 
    body_update_flags = 
    dealii::update_values
    | dealii::update_gradients
    | dealii::update_quadrature_points
    | dealii::update_JxW_values,    
    face_update_flags = 
    dealii::update_values
    | dealii::update_quadrature_points
    | dealii::update_JxW_values
    | dealii::update_normal_vectors;  
  
  dealii::FEValues<dim-1,dim>
    fe_v_TraceU_Dir(faceMapping,
  		    feTraceU_Dir,
  		    face_quadrature_formula,
  		    body_update_flags);
  
  dealii::FEValues<dim-1,dim>
    fe_v_TraceU_Neu(faceMapping,
  		    feTraceU_Neu,
  		    face_quadrature_formula,
  		    body_update_flags);
  
  dealii::Vector<heat::real> 
    U_Vector(dofs_per_cell_U),
    Q_Vector(dofs_per_cell_Q),
    TrU_Dir_U_Vector(dofs_per_cell_TraceU_Dir),
    TrU_Neu_Q_Vector(dofs_per_cell_TraceU_Neu);  
  
  
  typedef typename dealii::FullMatrix<heat::real> FullMatrix_t;
  
  FullMatrix_t 
    U_U_Matrix         (dofs_per_cell_U,dofs_per_cell_U),
    U_U_Matrix_Boundary(dofs_per_cell_U,dofs_per_cell_U),
    U_U_Matrix_Boundary_temp(dofs_per_cell_U,dofs_per_cell_U),
    U_Q_Matrix         (dofs_per_cell_U,dofs_per_cell_Q),
    U_Q_Matrix_2       (dofs_per_cell_U,dofs_per_cell_Q),
    Q_U_Matrix         (dofs_per_cell_Q,dofs_per_cell_U),
    Q_U_Matrix_2       (dofs_per_cell_Q,dofs_per_cell_U),
    Q_Q_Matrix         (dofs_per_cell_Q,dofs_per_cell_Q),
    TrU_Dir_U_Matrix   (dofs_per_cell_TraceU_Dir,dofs_per_cell_U),
    TrU_Neu_Q_Matrix   (dofs_per_cell_TraceU_Neu,dofs_per_cell_Q);

  //grvy_timer_end("assemble_system_inside");
  // http://www.dealii.org/developer/doxygen/deal.II/group__Iterators.html
  {
    for(auto 
	  cell_sys = dof_handler.begin_active(),
	  cell_U = dof_handlerU.begin_active(),
	  cell_Q = dof_handlerQ.begin_active();
	
	cell_sys != dof_handler.end();
	
	++cell_sys, 
	++cell_U,
	++cell_Q
	){

      dealii::FEValues<dim>
	fe_v_body_U(mapping, feU,  body_quadrature_formula,
		    body_update_flags);
     
      dealii
	::FEValues<dim>       
	fe_v_body_Q(mapping,
		    feQ,
		    body_quadrature_formula,
		    body_update_flags);
    

      fe_v_body_U.reinit(cell_U);
      fe_v_body_Q.reinit(cell_Q);

      //grvy_timer_begin("assemble MassU");
      // Assemble MassU and InverseMassU Local:
      U_U_Matrix = 0;
      U_U_Matrix_Boundary = 0;
      U_U_Matrix_Boundary_temp = 0;
      assembleMassULocal(fe_v_body_U,U_U_Matrix);    
      InsertGlobalFromLocal(cell_U,cell_U,U_U_Matrix,MassU);
      //grvy_timer_end("assemble MassU");
      
      //grvy_timer_begin("assemble MassUInverse");
      heat::myMatrixTools::ChopMatrix(U_U_Matrix,tolerance);
      {
	dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;
	tempLAPACKMatrix.copy_from(U_U_Matrix);
	tempLAPACKMatrix.invert();      
	U_U_Matrix = tempLAPACKMatrix;
      }
      
      InsertGlobalFromLocal(cell_U,cell_U,U_U_Matrix,InverseMassU);
      //grvy_timer_end("assemble MassUInverse");

      
      // Assemble StiffUfromQ Local;
      //grvy_timer_begin("assemble stiffUFromQ local");
      U_Q_Matrix = 0;
      assembleStiffUFromQLocal(fe_v_body_U,fe_v_body_Q,U_Q_Matrix);
      InsertGlobalFromLocal(cell_U,cell_Q,U_Q_Matrix,StiffUFromQ);
      //grvy_timer_end("assemble stiffUFromQ local");
    
      //Assemble MassQ Local;
      //grvy_timer_begin("assemble massQ local");
      Q_Q_Matrix = 0;
      assembleMassQLocal(fe_v_body_Q,Q_Q_Matrix);
      InsertGlobalFromLocal(cell_Q,cell_Q,Q_Q_Matrix,MassQ);
      //grvy_timer_end("assemble massQ local");

      //Assemble InverseMassQ Local
      //grvy_timer_begin("assemble inverseMassQ local");
      heat::myMatrixTools::ChopMatrix(Q_Q_Matrix,tolerance);
      {
	dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix;
	tempLAPACKMatrix.copy_from(Q_Q_Matrix);
	tempLAPACKMatrix.invert();      
	Q_Q_Matrix = tempLAPACKMatrix;
      }
      InsertGlobalFromLocal(cell_Q,cell_Q,Q_Q_Matrix,InverseMassQ);
      //grvy_timer_end("assemble inverseMassQ local");
      
      // Assemble StiffQFromU Local;
      //grvy_timer_begin("assemble stiffQFromU local");
      Q_U_Matrix = 0;
      assembleStiffQFromULocal(fe_v_body_Q,fe_v_body_U,Q_U_Matrix);
      InsertGlobalFromLocal(cell_Q,cell_U,Q_U_Matrix,StiffQFromU);
      //grvy_timer_end("assemble stiffQFromU local");

      // Assemble Source Vector
      //grvy_timer_begin("assemble source vector");
      U_Vector = 0;
      assembleU_RHS_FromSourceLocal(fe_v_body_U,U_Vector);      
      InsertGlobalFromLocalVector(cell_U,U_Vector,U_MinusRHS);

      //U_Vector *= -1.0;
      //InsertGlobalFromLocalVector(cell_U,U_Vector,U_RHS);
      
      //grvy_timer_end("assemble source vector");
    
      //Now loop over the faces;
      //grvy_timer_begin("loop over faces");
      for(unsigned int face_no = 0; 
	face_no < dealii::GeometryInfo<dim>::faces_per_cell;
	++face_no){

      typename 
	dealii::DoFHandler<dim>::face_iterator 
	face_sys = cell_sys->face(face_no),
	face_U = cell_U->face(face_no);
	//face_Q = cell_Q->face(face_no);
      
      if( face_sys->at_boundary() ){
	
	dealii::FEFaceValues<dim>    
	  fe_v_face_U(mapping, feU, face_quadrature_formula,
		      face_update_flags);

	dealii::FEFaceValues<dim>    
	  fe_v_face_Q(mapping,
		      feQ,
		      face_quadrature_formula,
		      face_update_flags);	

	fe_v_face_U.reinit(cell_U,face_no);
	fe_v_face_Q.reinit(cell_Q,face_no);	
	
	
	if( boundary_ids_Dir_plus.count(face_sys->boundary_indicator()) != 0){
	  //Dir plus
	  // FluxUFromQBoundary
	  U_Q_Matrix = 0;
	  assembleFluxUFromQBoundaryLocal(fe_v_face_U,fe_v_face_Q,U_Q_Matrix);
	  InsertGlobalFromLocal(cell_U,cell_Q,U_Q_Matrix,FluxBoundaryUFromQ);	  

	  Q_Vector = 0;
	  assembleQ_RHS_FromDirPlusLocal(fe_v_face_Q,Q_Vector);
	  InsertGlobalFromLocalVector(cell_Q,Q_Vector,Q_MinusRHS);
	  //Q_Vector *= -1.0;	  
	  //InsertGlobalFromLocalVector(cell_Q,Q_Vector,Q_RHS);		  
	}

	else if(boundary_ids_Dir_minus.count(face_sys->boundary_indicator()) != 0){	  

	  // FluxQFromUBoundary
	  Q_U_Matrix = 0;
	  assembleFluxQFromUBoundaryLocal(fe_v_face_Q,fe_v_face_U,Q_U_Matrix);
	  InsertGlobalFromLocal(cell_Q,cell_U,Q_U_Matrix,FluxBoundaryQFromU);

	  // Boundary massMatrix;
	  U_U_Matrix_Boundary_temp = 0;
	  assembleBoundaryMassULocal(fe_v_face_U, U_U_Matrix_Boundary_temp);
	  InsertGlobalFromLocal(cell_U,
				cell_U,
				U_U_Matrix_Boundary_temp,
				BoundaryMassU);

	  U_U_Matrix_Boundary.add(1.0, U_U_Matrix_Boundary_temp);

	  // Dirichelet Boundary vector;
	  U_Vector = 0;	  
	  assembleConstraintUDirMinusLocalNOTRACE(fe_v_face_U,U_Vector);
	  InsertGlobalFromLocalVector(cell_U,U_Vector,UdirConstraint_RHS);

	  // Lagrange Multiplier for U
	  auto cell_TrU_Dir
	    = DirVolumeToSurfaceTriangulationMap.at(face_U);

	  TrU_Dir_U_Matrix = 0;
	  fe_v_TraceU_Dir.reinit(cell_TrU_Dir);
	  assembleTraceULocalDir(fe_v_TraceU_Dir,
				 fe_v_face_U,
				 TrU_Dir_U_Matrix);	  

	  InsertGlobalFromLocal(cell_TrU_Dir,
				cell_U,
				TrU_Dir_U_Matrix,
				TraceUDir);

	  // Constraint Vector
	  TrU_Dir_U_Vector = 0;
	  assembleConstraintUDirMinusLocal(fe_v_TraceU_Dir,TrU_Dir_U_Vector);
	  InsertGlobalFromLocalVectorNEW(cell_TrU_Dir,TrU_Dir_U_Vector,constraintDirU);


	}
	else if(boundary_ids_Neu_plus.count(face_sys->boundary_indicator()) != 0){
	  std::cout << "On Neuman Boundary Plus" << std::endl;

	  // FluxUFromQBoundary
	  U_Q_Matrix = 0;
	  assembleFluxUFromQBoundaryLocal(fe_v_face_U,fe_v_face_Q,U_Q_Matrix);
	  InsertGlobalFromLocal(cell_U,cell_Q,U_Q_Matrix,FluxBoundaryUFromQ);
	  
	  auto cell_TrU_Neu
	    = NeuVolumeToSurfaceTriangulationMap.at(face_U);
	  
	  fe_v_TraceU_Neu.reinit(cell_TrU_Neu);	  
	  fe_v_face_Q.reinit(cell_Q,face_no);
	  

	  TrU_Neu_Q_Matrix = 0;	  
	  assembleTraceQLocalNeu
	    (fe_v_TraceU_Neu,
	     fe_v_face_Q,
	     TrU_Neu_Q_Matrix);	  
	  InsertGlobalFromLocal
	    (cell_TrU_Neu,
	     cell_Q,
	     TrU_Neu_Q_Matrix,
	     TraceQNeu);

	  //Constraint for Q;
	  TrU_Neu_Q_Vector = 0;

	  assembleConstraintQNeuPlusLocal(fe_v_TraceU_Neu, fe_v_face_Q, TrU_Neu_Q_Vector);
	  InsertGlobalFromLocalVectorNEW(cell_TrU_Neu, TrU_Neu_Q_Vector, constraintNeuQ);
	  
	  
	}
	else if(boundary_ids_Neu_minus.count(face_sys->boundary_indicator() ) != 0){
	  std::cout << "On Neuman Boundary minus" << std::endl;

	  // FluxQFromUBoundary
	  Q_U_Matrix = 0;	  
	  assembleFluxQFromUBoundaryLocal(fe_v_face_Q,fe_v_face_U,Q_U_Matrix);
	  InsertGlobalFromLocal(cell_Q,cell_U,Q_U_Matrix,FluxBoundaryQFromU);

	  U_Vector = 0;
	  assembleQ_RHS_FromNeuMinusLocal(fe_v_face_U, U_Vector);
	  InsertGlobalFromLocalVector(cell_U,U_Vector,U_MinusRHS);	

	  
	  
	}
	else if(boundary_ids_Rob_plus.count(face_sys->boundary_indicator() ) != 0){
	  std::cout << "On Robin Boundary plus" << std::endl;
	  // TODO:  Turn this into a dealii assertion.
	  std::cerr << "Robin BC's are not implemented yet" << std::endl;
	  std::abort();
	  //std::cout << "Robin Boundary" << std::endl;
	}

	else if(boundary_ids_Rob_minus.count(face_sys->boundary_indicator() ) != 0){
	  std::cout << "On Robin Boundary minus" << std::endl;
	  // TODO:  Turn this into a dealii assertion.
	  std::cerr << "Robin BC's are not implemented yet" << std::endl;
	  std::abort();
	  //std::cout << "Robin Boundary" << std::endl;
	}
	

	else {
	  //TODO:  Change this to a dealii assert
	  std::cerr
	    << "There is a face on the boundary" 
	    << "that isn't assigned a boundary condition."
	    << std::endl;
	  std::abort();
	}	
      }
      
      else{ 
	//now we know that the current face is not at a boundary
	//grvy_timer_begin("something");
	Assert( 
	       cell_sys-> neighbor(face_no).state() 
	       == dealii::IteratorState::valid,
	       dealii::ExcInternalError() );
	
	typename 
	  dealii::DoFHandler<dim>::cell_iterator 	  
	  neighbor_U = cell_U->neighbor(face_no),
	  neighbor_Q = cell_Q->neighbor(face_no);
	if(face_sys->has_children()){
	  //grvy_timer_begin("face_sys->has_children()");
	  const unsigned int neighbor2 = 
	    cell_sys->neighbor_face_no(face_no);
	  const unsigned int neighbor2_U = 
	    cell_U->neighbor_face_no(face_no);
	  const unsigned int neighbor2_Q = 
	    cell_Q->neighbor_face_no(face_no);
	  
	  assert(neighbor2 == neighbor2_U);
	  assert(neighbor2 == neighbor2_Q);
	  for( unsigned int subface_no = 0;
	       subface_no < face_sys->number_of_children();
	       ++subface_no){
	    const typename dealii::DoFHandler<dim>::cell_iterator 
	      neighbor_child_sys
	      = cell_sys->neighbor_child_on_subface(face_no, subface_no),
	      neighbor_child_U
	      = cell_U->neighbor_child_on_subface(face_no, subface_no),
	      neighbor_child_Q
	      = cell_Q->neighbor_child_on_subface(face_no, subface_no);

	    const auto 
	      neighbor_U_for_assembly = neighbor_child_U,
	      neighbor_Q_for_assembly = neighbor_child_Q;
	    
	    
	    Assert(!(neighbor_child_sys->has_children())
		   , dealii::ExcInternalError() );

	    dealii::FESubfaceValues<dim> 
	      fe_v_self_U(mapping, feU, face_quadrature_formula,
			  face_update_flags);
	    
	    dealii::FESubfaceValues<dim> 
	      fe_v_self_Q(mapping,
			  feQ,
			  face_quadrature_formula,
			  face_update_flags);
	    
	    dealii::FEFaceValues<dim> 
	      fe_v_neig_U(mapping,
			  feU,
			  face_quadrature_formula,
			  face_update_flags);

	    dealii::FEFaceValues<dim> 
	      fe_v_neig_Q(mapping,
			  feQ, 
			  face_quadrature_formula,
			  face_update_flags);
	    
	    fe_v_self_U.reinit(cell_U,face_no,subface_no);
	    fe_v_self_Q.reinit(cell_Q,face_no,subface_no);
	    fe_v_neig_U.reinit(neighbor_child_U,neighbor2);
	    fe_v_neig_Q.reinit(neighbor_child_Q,neighbor2);
	    
	    assembleFluxLocals(fe_v_self_U,
	    		       fe_v_neig_U,
	    		       fe_v_self_Q,
	    		       fe_v_neig_Q,
	    		       cell_U,
	    		       cell_Q,
	    		       neighbor_U_for_assembly,
	    		       neighbor_Q_for_assembly,
	    		       U_Q_Matrix,
	    		       U_Q_Matrix_2,
	    		       Q_U_Matrix,
	    		       Q_U_Matrix_2);	    
	    
	  }
	  //grvy_timer_end("face_sys->has_children()");
	}
	else{ // !(face->has_children())
	  if(!(cell_sys->neighbor_is_coarser( face_no ))){
	    //grvy_timer_begin("!cell_sys->neighbor_is_coarser");
	    //Now we know that the neighbor's face is neither finer nor coarser

	    //grvy_timer_begin("line1");
	    const unsigned int neighbor2 
	      = cell_sys->neighbor_of_neighbor(face_no);
	    //grvy_timer_end("line1");
	    
	    //grvy_timer_begin("line2");
	    dealii::FEFaceValues<dim>    
	      fe_v_self_U(mapping,feU,face_quadrature_formula,
			  face_update_flags);
	    //grvy_timer_end("line2");

	    
	    //grvy_timer_begin("line3");
	    dealii::FEFaceValues<dim>    
	      fe_v_self_Q(mapping,
			  feQ,
			  face_quadrature_formula,
			  face_update_flags);
	    //grvy_timer_end("line3");
	      

	    //grvy_timer_begin("line4");
	    dealii::FEFaceValues<dim> 
	      fe_v_neig_U(mapping,
			  feU,
			  face_quadrature_formula,
			  face_update_flags);
	    //grvy_timer_end("line4");

	    //grvy_timer_begin("line5");
	    dealii::FEFaceValues<dim> 
	      fe_v_neig_Q(mapping,
			  feQ, 
			  face_quadrature_formula,
			  face_update_flags);
	    //grvy_timer_end("line5");
	    

	    //grvy_timer_begin("line6");
	    fe_v_self_U.reinit(cell_U,face_no);
	    fe_v_self_Q.reinit(cell_Q,face_no);
	    fe_v_neig_U.reinit(neighbor_U, neighbor2);
	    fe_v_neig_Q.reinit(neighbor_Q, neighbor2);
	    //grvy_timer_end("line6");
	    
	    //grvy_timer_begin("line7");
	    const auto 
	      neighbor_U_for_assembly = neighbor_U,
	      neighbor_Q_for_assembly = neighbor_Q;
	    //grvy_timer_end("line7");
	    
	    //grvy_timer_begin("assemble flux locals");
	    assembleFluxLocals
	      (fe_v_self_U,
	       fe_v_neig_U,
	       fe_v_self_Q,
	       fe_v_neig_Q,
	       cell_U,
	       cell_Q,
	       neighbor_U_for_assembly,
	       neighbor_Q_for_assembly,
	       U_Q_Matrix,
	       U_Q_Matrix_2,
	       Q_U_Matrix,
	       Q_U_Matrix_2);
	    //grvy_timer_end("assemble flux locals");
	    //grvy_timer_end("!cell_sys->neighbor_is_coarser"); 
	  } else {
	    // Ok, now we now that we are in the interior, on a face whose
	    // neighbor is coarser.  This is the opposite of the situation
	    // of the first thing we checked.
	    
	    
	    const std::pair<unsigned int, unsigned int> neighbor_face_subface
	      = cell_sys->neighbor_of_coarser_neighbor(face_no);
	    
	    const unsigned int nFace_no = neighbor_face_subface.first;
	    const unsigned int nSubface_no = neighbor_face_subface.second;
	    
	    //TODO:  Change to a dealii style assert
	    assert
	      (cell_sys
	       ->neighbor(face_no)
	       ->face(nFace_no)
	       ->child(nSubface_no) 
	       == cell_sys->face(face_no)
	       );
	    
	    typename dealii::DoFHandler<dim>::cell_iterator
	      neighbor_U = cell_U->neighbor(face_no),
	      neighbor_Q = cell_Q->neighbor(face_no);

	    dealii::FEFaceValues<dim>    
	      fe_v_self_U(mapping,feU,face_quadrature_formula,
			  face_update_flags);
	    
	    dealii::FEFaceValues<dim>    
	      fe_v_self_Q(mapping,
			  feQ,
			  face_quadrature_formula,
			  face_update_flags);

	    dealii::FESubfaceValues<dim>
	      fe_v_neig_U(mapping,
			  feU,
			  face_quadrature_formula,
			  face_update_flags);
	      
	    dealii::FESubfaceValues<dim>
	      fe_v_neig_Q(mapping,
			  feQ, 
			  face_quadrature_formula,
			  face_update_flags);	    
	    
	    fe_v_self_U.reinit(cell_U,face_no);
	    fe_v_self_Q.reinit(cell_Q,face_no);
	    fe_v_neig_U.reinit(neighbor_U,nFace_no,nSubface_no);
	    fe_v_neig_Q.reinit(neighbor_Q,nFace_no,nSubface_no);

	    const auto 
	      neighbor_U_for_assembly = neighbor_U,
	      neighbor_Q_for_assembly = neighbor_Q;
	    
	    assembleFluxLocals
	      (fe_v_self_U,
	       fe_v_neig_U,
	       fe_v_self_Q,
	       fe_v_neig_Q,
	       cell_U,
	       cell_Q,
	       neighbor_U_for_assembly,
	       neighbor_Q_for_assembly,
	       U_Q_Matrix,
	       U_Q_Matrix_2,
	       Q_U_Matrix,
	       Q_U_Matrix_2);	    	    
	  }
	}
	//grvy_timer_end("something");
      } // End of face not at boundary	
      
      } // End of loops through faces

      //grvy_timer_end("loop over faces");
    // Now we build the Singular value decompositions for the boundary mass matrix;
      //grvy_timer_begin("compute svds");
    heat::myMatrixTools::ChopMatrix(U_U_Matrix_Boundary, tolerance);
    {
      dealii::LAPACKFullMatrix<heat::real> tempLAPACKMatrix,
	tempLAPACKMatrix2;

      dealii::FullMatrix<heat::real> myEigenvectors;

      dealii::Vector<heat::real> myEigenvalues;
      tempLAPACKMatrix.copy_from(U_U_Matrix_Boundary);
      tempLAPACKMatrix2.copy_from(U_U_Matrix);

      //tempLAPACKMatrix.compute_generalized_eigenvalues_symmetric(tempLAPACKM


      tempLAPACKMatrix.compute_eigenvalues_symmetric
    	(-2,
	 100,
    	 1e-38,
    	 myEigenvalues,
    	 myEigenvectors);

      tempLAPACKMatrix = myEigenvectors;
      InsertGlobalFromLocal(cell_U,cell_U,myEigenvectors,Svd_U);

      U_U_Matrix_Boundary = 0;

      assert(myEigenvalues.size() == U_U_Matrix_Boundary.m() );

      for(unsigned int i = 0; i < myEigenvalues.size(); ++i){
    	 if(fabs(myEigenvalues[i]) < tolerance)  {
    	   myEigenvalues[i] = 0.0;
	 }
	 U_U_Matrix_Boundary(i,i) = myEigenvalues[i];

      }	
      InsertGlobalFromLocal(cell_U,cell_U,U_U_Matrix_Boundary,Svd_Sigma);

    }
    //grvy_timer_end("compute svds");
    } // End of loops through cells

  }

  //grvy_timer_begin("PP and QQ");
  rankBoundaryMassU = 0;
  nullityBoundaryMassU = 0;
  std::vector<unsigned int> Prows,Qrows;
  
  
  for(unsigned int i = 0; i < Svd_Sigma.m(); ++i){
    if( fabs( Svd_Sigma(i,i) ) > tolerance ){
      ++rankBoundaryMassU;
      Prows.emplace_back(i);
    }
    else{
      ++nullityBoundaryMassU;
      Qrows.emplace_back(i);
    }    
  }

  PP.reinit(rankBoundaryMassU,n_dofs_U);
  QQ.reinit(nullityBoundaryMassU,n_dofs_U);

  for(unsigned int i = 0; i < rankBoundaryMassU; ++i){
    PP.add(i,Prows[i],1.0);
  }

  for(unsigned int i = 0; i < nullityBoundaryMassU; ++i){
    QQ.add(i,Qrows[i],1.0);
  }
  //grvy_timer_end("PP and QQ");
  //lambdaU.reinit(rankBoundaryMassU);

  //MassUInverseDirect.initialize(MassU);
  //MassQInverseDirect.initialize(MassQ);

  
  

  //TODO:  Not the correct spot for this.

  //grvy_timer_begin("TotalUFromQ summing");
  const heat::real theta = 0.0;
  TotalUFromQ.copy_from(StiffUFromQ);
  TotalUFromQ.add(-1.0 * theta,FluxCenterUFromQ);
  TotalUFromQ.add(-1.0 * (1-theta),FluxPosUFromQ);
  TotalUFromQ.add(-1.0 * (1-theta),FluxNegUFromQ);
  TotalUFromQ.add( -1.0 , FluxBoundaryUFromQ);
  //grvy_timer_end("TotalUFromQ summing");

  //grvy_timer_begin("TotalQFromU summing");
  TotalQFromU.copy_from(StiffQFromU);
  TotalQFromU.add(-1.0 * theta,FluxCenterQFromU);
  TotalQFromU.add(-1.0 * (1.0 - theta), FluxPosQFromU);
  TotalQFromU.add(-1.0 * (1.0 - theta), FluxNegQFromU);
  TotalQFromU.add(-1.0, FluxBoundaryQFromU);
  //grvy_timer_end("TotalQFromU summing");
}

template<int dim>
void
EllipticProblem<dim>
::print_system()
{
  
  std::vector<std::pair< dealii::SparseMatrix< heat::real > * , std::string > > list;

  list.emplace_back(&MassU,              "MassU.mtx");
  list.emplace_back(&BoundaryMassU,      "BoundaryMassU.mtx");
  list.emplace_back(&Svd_U,              "Svd_U.mtx");
  list.emplace_back(&Svd_Sigma,           "Svd_Sigma.mtx");
  list.emplace_back(&InverseMassU,       "InverseMassU.mtx");
  list.emplace_back(&StiffUFromQ,        "StiffUFromQ.mtx");
  list.emplace_back(&FluxPosUFromQ,      "FluxPosUFromQ.mtx");
  list.emplace_back(&FluxNegUFromQ,      "FluxNegUFromQ.mtx");
  list.emplace_back(&FluxCenterUFromQ,   "FluxCenterUFromQ.mtx");
  list.emplace_back(&FluxBoundaryUFromQ, "FluxBoundaryUFromQ.mtx");
  list.emplace_back(&TotalUFromQ,        "TotalUFromQ.mtx");  
  list.emplace_back(&MassQ,              "MassQ.mtx");
  list.emplace_back(&InverseMassQ,       "InverseMassQ.mtx");
  list.emplace_back(&StiffQFromU,        "StiffQFromU.mtx");
  list.emplace_back(&FluxPosQFromU,      "FluxPosQFromU.mtx");
  list.emplace_back(&FluxNegQFromU,      "FluxNegQFromU.mtx");
  list.emplace_back(&FluxCenterQFromU,   "FluxCenterQFromU.mtx");
  list.emplace_back(&FluxBoundaryQFromU, "FluxBoundaryQFromU.mtx");
  list.emplace_back(&TotalQFromU,        "TotalQFromU.mtx");
  
  for( auto it = list.begin(); it != list.end(); it++){
    std::cout 
      << "\n" << it->second << ".n_nonzero_elements() = " 
      << it->first->n_nonzero_elements() << "\n"
      << it->second << ".n_actual_nonzero_elements() = "
      << it->first->n_actually_nonzero_elements() << "\n";    

    std::ofstream ostr(it->second);
    ostr << std::setprecision(15);
    heat::myMatrixTools::PrintMatrixMarket(*(it->first), ostr);
  }
	
  std::pair< dealii::SparseMatrixEZ<heat::real > *, std::string > temp2;
  std::vector < decltype(temp2) > list2;

  list2.emplace_back(&TraceUDir, "TraceUDir.mtx");
  list2.emplace_back(&TraceQNeu, "TraceQNeu.mtx");
  list2.emplace_back(&PP,        "PP.mtx");
  list2.emplace_back(&QQ,        "QQ.mtx");
  
  for( auto it = list2.begin();
       it != list2.end();
       it++){

    std::string something = it->second;

    std::ofstream ostr(it->second);

    ostr << std::setprecision(15);
    heat::myMatrixTools::PrintMatrixMarket(*(it->first), ostr);
  }

  std::pair< dealii::Vector<heat::real>*, std::string > temp3;
  std::vector <decltype(temp3) > list3;

  //list3.emplace_back(&U_RHS, "U_RHS.dat");
  list3.emplace_back(&U_MinusRHS, "U_MinusRHS.dat");
  list3.emplace_back(&Ustate, "Ustate.dat");
  list3.emplace_back(&Qstate, "Qstate.dat");
  list3.emplace_back(&UperpNEW, "Uperp.dat");

  list3.emplace_back(&UparNEW, "Upar.dat");
  //list3.emplace_back(&Q_RHS, "Q_RHS.dat");
  list3.emplace_back(&Q_MinusRHS, "Q_MinusRHS.dat");
  list3.emplace_back(&constraintDirU, "constraintDirU.dat");
  list3.emplace_back(&UdirConstraint_RHS, "uDirConstraint_RHS.dat");
  if(constraintNeuQ.size() != 0){
    list3.emplace_back(&constraintNeuQ, "constraintNeuQ.dat");
  }
  for(auto it = list3.begin();
      it != list3.end();
      ++it){
    
    std::ofstream ostr(it->second);
    it->first->print(ostr,
		     12, // Precision
		     true, // Scientific
		     true  // across
		     );

  }
  
}


template<int dim>
void
EllipticProblem<dim>
::assembleFluxLocals
(const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> & fe_v_neig_U,
 const dealii::FEFaceValuesBase<dim> & fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> & fe_v_neig_Q,
 const typename dealii::DoFHandler<dim>::active_cell_iterator & cell_U,
 const typename dealii::DoFHandler<dim>::active_cell_iterator & cell_Q,
 const typename dealii::DoFHandler<dim>::active_cell_iterator & neighbor_U,
 const typename dealii::DoFHandler<dim>::active_cell_iterator & neighbor_Q,
 dealii::FullMatrix<heat::real> U_Q_Matrix,
 dealii::FullMatrix<heat::real> U_Q_Matrix_2,
 dealii::FullMatrix<heat::real> Q_U_Matrix,
 dealii::FullMatrix<heat::real> Q_U_Matrix_2)
{
  Q_U_Matrix = 0;
  Q_U_Matrix_2 = 0;
  assembleFluxCenterQFromULocal	      
    (fe_v_self_Q, fe_v_self_U,
     fe_v_neig_U, Q_U_Matrix, Q_U_Matrix_2);	    
  InsertGlobalFromLocal
    (cell_Q, cell_U,
     Q_U_Matrix,   FluxCenterQFromU);	    
	    
  InsertGlobalFromLocal
    (cell_Q, neighbor_U,
     Q_U_Matrix_2, FluxCenterQFromU);
	    
  Q_U_Matrix = 0;
  Q_U_Matrix_2 = 0;
  assembleFluxPosQFromULocal
    (fe_v_self_Q, fe_v_self_U,
     fe_v_neig_U, Q_U_Matrix, Q_U_Matrix_2);	    
  InsertGlobalFromLocal
    (cell_Q, cell_U,
     Q_U_Matrix, FluxPosQFromU);	    
  InsertGlobalFromLocal
    (cell_Q, neighbor_U,
     Q_U_Matrix_2, FluxPosQFromU);
	    
  Q_U_Matrix = 0;
  Q_U_Matrix_2 = 0;
  assembleFluxNegQFromULocal
    (fe_v_self_Q, fe_v_self_U,
     fe_v_neig_U, Q_U_Matrix, Q_U_Matrix_2);	    
  InsertGlobalFromLocal
    (cell_Q, cell_U,
     Q_U_Matrix, FluxNegQFromU);	    
  InsertGlobalFromLocal
    (cell_Q, neighbor_U,
     Q_U_Matrix_2, FluxNegQFromU);
	    
  U_Q_Matrix = 0;
  U_Q_Matrix_2 = 0;	    
  assembleFluxCenterUFromQLocal
    (fe_v_self_U, fe_v_self_Q,
     fe_v_neig_Q, U_Q_Matrix, U_Q_Matrix_2);	    
  InsertGlobalFromLocal
    (cell_U,cell_Q,
     U_Q_Matrix,FluxCenterUFromQ);
  InsertGlobalFromLocal
    (cell_U, neighbor_Q,
     U_Q_Matrix_2, FluxCenterUFromQ);
	    
  U_Q_Matrix = 0;
  U_Q_Matrix_2 = 0;
  assembleFluxPosUFromQLocal
    (fe_v_self_U, fe_v_self_Q,
     fe_v_neig_Q, U_Q_Matrix, U_Q_Matrix_2);	    
  InsertGlobalFromLocal
    (cell_U,cell_Q,
     U_Q_Matrix,  FluxPosUFromQ);
  InsertGlobalFromLocal
    (cell_U,neighbor_Q,
     U_Q_Matrix_2,FluxPosUFromQ);
	    
  U_Q_Matrix = 0;
  U_Q_Matrix_2 = 0;
  assembleFluxNegUFromQLocal
    (fe_v_self_U, fe_v_self_Q,
     fe_v_neig_Q, U_Q_Matrix, U_Q_Matrix_2);	    
  InsertGlobalFromLocal
    (cell_U,cell_Q,
     U_Q_Matrix,FluxNegUFromQ);
  InsertGlobalFromLocal
    (cell_U,neighbor_Q,
     U_Q_Matrix_2,FluxNegUFromQ);	    
}
 
 

 
		      

template<int dim>
void
EllipticProblem<dim>
::assembleMassULocal(const dealii::FEValues<dim> & fe_v,
		     dealii::FullMatrix<heat::real> & ui_vi_matrix)
{
  const dealii::FEValuesExtractors::Scalar density(0);
  const std::vector<double> & JxW = fe_v.get_JxW_values();
  for(unsigned int point = 0; point < fe_v.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v.dofs_per_cell; ++j){	
	ui_vi_matrix(i,j) += 
	  fe_v[density].value(i,point)
	  * JxW[point]
	  * fe_v[density].value(j,point);
      }
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleBoundaryMassULocal(const dealii::FEFaceValuesBase<dim> & fe_v,
			     dealii::FullMatrix<heat::real> & ui_vi_matrix)
{
  const std::vector<double> & JxW = fe_v.get_JxW_values();
  for(unsigned int point = 0; point < fe_v.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v.dofs_per_cell; ++j){
	ui_vi_matrix(i,j) +=
	  JxW[point]
	  * fe_v.shape_value(i,point)
	  * fe_v.shape_value(j,point);
      }
    }
  }
  
  // const dealii::FEValuesExtractors::Scalar density(0);
  // const std::vector<double> & JxW = fe_v.get_JxW_values();
  // for(unsigned int point = 0; point < fe_v.n_quadrature_points; ++point){
  //   for(unsigned int i = 0; i < fe_v.dofs_per_cell; ++i){
  //     for(unsigned int j = 0; j < fe_v.dofs_per_cell; ++j){	
  // 	ui_vi_matrix(i,j) += 
  // 	  fe_v[density].value(i,point)
  // 	  * JxW[point]
  // 	  * fe_v[density].value(j,point);
  //     }
  //   }
  // }
}


template<int dim>
void
EllipticProblem<dim>
::assembleMassQLocal(const dealii::FEValues<dim> & fe_v,		    
		     dealii::FullMatrix<heat::real> & Qi_Qi_matrix)
{
  //We are writing this in a more "speed optimum" way because profiling has shown this to be a bottleneck.
  const dealii::FEValuesExtractors::Vector flux(0);
  const std::vector<double> & JxW = fe_v.get_JxW_values();
  for(unsigned int point = 0; point < fe_v.n_quadrature_points; ++point){
    const auto JxWAtPoint = JxW[point];
    for(unsigned int i = 0; i < fe_v.dofs_per_cell; ++i){
      const auto value1 = fe_v[flux].value(i,point);
      //const auto valueProd = JxWAtPoint * value1;
      for(unsigned int j = 0; j < fe_v.dofs_per_cell; ++j){
	const auto value2 = fe_v[flux].value(j, point);
	
	// Qi_Qi_matrix(i,j) += 
	//   JxW[point] 
	//   * fe_v[flux].value(i,point)
	//   * fe_v[flux].value(j,point);

	Qi_Qi_matrix(i,j) += JxWAtPoint * value1 * value2;	 
      }
    }
  }
}





template<int dim>
void
EllipticProblem<dim>
::assembleStiffUFromQLocal(const dealii::FEValues<dim> & fe_v_U,
			   const dealii::FEValues<dim> & fe_v_Q,
			   dealii::FullMatrix<heat::real> & ui_qj_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);
  const std::vector<double> & JxW = fe_v_U.get_JxW_values();
  for(unsigned int point = 0; point < fe_v_U.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_U.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_Q.dofs_per_cell; ++j){
	ui_qj_matrix(i,j) += 
	  JxW[point] 
	  * fe_v_U[density].gradient(i,point)
	  * fe_v_Q[flux].value(j,point);	
      }
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleStiffQFromULocal(const dealii::FEValues<dim> & fe_v_Q,
			   const dealii::FEValues<dim> & fe_v_U,
			   dealii::FullMatrix<heat::real> & qi_uj_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);
  const std::vector<double> & JxW = fe_v_Q.get_JxW_values();
  for(unsigned int point = 0; point < fe_v_Q.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_Q.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_U.dofs_per_cell; ++j){
	qi_uj_matrix(i,j) += 
	  JxW[point]
	  * fe_v_Q[flux].divergence(i,point)
	  * fe_v_U[density].value(j,point);
      }
    }
  }  
}


template<int dim>
void
EllipticProblem<dim>
::assembleFluxCenterUFromQLocal
(const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> &fe_v_neig_Q,
 dealii::FullMatrix<heat::real> & Uin_Qin_matrix,
 dealii::FullMatrix<heat::real> & Uin_Qout_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const std::vector<double> & JxW = fe_v_self_U.get_JxW_values();
  const std::vector<dealii::Point<dim> > & normals 
    = fe_v_self_U.get_normal_vectors();
  heat::real onehalf = 0.5;
  for(unsigned int point = 0; point < fe_v_self_U.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_self_U.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_self_Q.dofs_per_cell; ++j){
	Uin_Qin_matrix(i,j) += 
	  onehalf * 
	  JxW[point] 
	  * fe_v_self_U.shape_value(i,point) 
	  * fe_v_self_Q[flux].value(j,point)
	  * normals[point];
      }
    }
  }
  for(unsigned int point = 0; point < fe_v_self_U.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_self_U.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_neig_Q.dofs_per_cell; ++j){
	Uin_Qout_matrix(i,j) += 
	  onehalf
	  * JxW[point] 
	  * fe_v_self_U.shape_value(i,point) 
	  * fe_v_neig_Q[flux].value(j,point) 
	  * normals[point];
      }
    }
  }
}

// template<int dim>
// void
// EllipticProblem<dim>
// ::assembleStiffUFromULocal
// (const dealii::FEValuesBase<dim> & fe_v_self_U,
//  const dealii::FEValuesBase<dim> & fe_v_self_Elec,
//  dealii::Vector<heat::real> & Uin_vector)
// {
//   const dealii::FEValuesExtractors::Vector elec(0);
//   const dealii::FEValuesExtractors::Scalar density(0);  

//   const std::vector<double> & JxW = fe_v_self_U.get_JxW_values();
  
//   std::vector < typename dealii::Tensor<1,dim> >
//     ElecSelfValues(fe_v_self_U.n_quadrature_points);
//   fe_v_self_Elec[elec].get_function_values(ElecState,ElecSelfValues);
  
//   std::vector < heat::real>
//     USelfValues(fe_v_self_U.n_quadrature_points);
//   fe_v_self_U[density].get_function_values(Ustate,USelfValues);

//   for(unsigned int point = 0; point < fe_v_self_U.n_quadrature_points;
//       ++point){
//     for(unsigned int i = 0; i < fe_v_self_U.dofs_per_cell; ++i){
//       Uin_vector(i) += JxW[point] 
// 	* ElecSelfValues[point] 
// 	* fe_v_self_U[density].gradient(i,point)
// 	* USelfValues[point];
//     }
//   }
// }

// template<int dim>
// void
// EllipticProblem<dim>
// ::assembleFluxUFromULocal
// (const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
//  const dealii::FEFaceValuesBase<dim> & fe_v_self_Elec,
//  const dealii::FEFaceValuesBase<dim> & fe_v_neig_U,
//  const dealii::FEFaceValuesBase<dim> & fe_v_neig_Elec,
//  dealii::Vector<heat::real> & Uin_vector)
// {
//   const dealii::FEValuesExtractors::Vector elec(0);
//   const dealii::FEValuesExtractors::Scalar density(0);
  
//   const auto & JxW = fe_v_self_U.get_JxW_values();
//   const auto & normals = fe_v_self_U.get_normal_vectors();
  
//   std::vector < typename dealii::Tensor<1,dim> >
//     ElecSelfValues(fe_v_self_U.n_quadrature_points),
//     ElecNeigValues(fe_v_neig_U.n_quadrature_points);
  
//   fe_v_self_Elec[elec].get_function_values(ElecState,ElecSelfValues);
//   fe_v_neig_Elec[elec].get_function_values(ElecState,ElecNeigValues);  

//   std::vector < heat::real> 
//     USelfValues(fe_v_self_U.n_quadrature_points),
//     UNeigValues(fe_v_neig_U.n_quadrature_points);
//   fe_v_self_U[density].get_function_values(Ustate, USelfValues);
//   fe_v_neig_U[density].get_function_values(Ustate, UNeigValues);
  
//   const heat::real zero = 0.0;
//   //heat::real temp1,temp2;
//   for(unsigned int point = 0;
//       point < fe_v_self_U.n_quadrature_points;
//       ++point){
//     for(unsigned int i = 0; i < fe_v_self_U.dofs_per_cell; ++i){
//       if(ElecSelfValues[point] * normals[point] > zero){
// 	Uin_vector(i) 
// 	  += JxW[point] 
// 	  * fe_v_self_U[density].value(i,point) * USelfValues[point]
// 	  * ElecSelfValues[point] * normals[point];
//       }
//       else {
// 	Uin_vector(i) 
// 	  += JxW[point] 
// 	  * fe_v_self_U[density].value(i,point) * UNeigValues[point]
// 	  * ElecSelfValues[point] * normals[point];      
//       }
//     }
//   }
// }

template<int dim>
void
EllipticProblem<dim>
::assembleFluxPosUFromQLocal
(const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> &fe_v_neig_Q,
 dealii::FullMatrix<heat::real> & Uin_Qin_matrix,
 dealii::FullMatrix<heat::real> & Uin_Qout_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);
  
  const std::vector<double> & JxW = fe_v_self_U.get_JxW_values();
  const std::vector<dealii::Point<dim> > & normals 
    = fe_v_self_U.get_normal_vectors();
  
  for(unsigned int point = 0; point < fe_v_self_U.n_quadrature_points; ++point){
    if(referenceDirection * normals[point] > 0){
      for(unsigned int i = 0; i < fe_v_self_U.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_v_neig_Q.dofs_per_cell; ++j){
	  Uin_Qin_matrix(i,j) += 
	    JxW[point] 
	    * fe_v_self_U[density].value(i,point) 
	    * fe_v_self_Q[flux].value(j,point)
	    * normals[point];
	}
      }
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleFluxNegUFromQLocal
(const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> &fe_v_neig_Q,
 dealii::FullMatrix<heat::real> & Uin_Qin_matrix,
 dealii::FullMatrix<heat::real> & Uin_Qout_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);  
  
  const auto & JxW = fe_v_self_U.get_JxW_values();
  const auto & normals = fe_v_self_U.get_normal_vectors();
  
  for(unsigned int point = 0; 
      point < fe_v_self_U.n_quadrature_points; 
      ++point){
    if(referenceDirection * normals[point] < 0){
      for(unsigned int i = 0; i < fe_v_self_U.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_v_self_Q.dofs_per_cell; ++j){
	  Uin_Qout_matrix(i,j) += 
	    JxW[point]
	    * fe_v_self_U[density].value(i,point) 
	    * fe_v_neig_Q[flux].value(j,point)
	    * normals[point];
	}
      }
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleFluxUFromQBoundaryLocal
(const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> & fe_v_self_Q,
 dealii::FullMatrix<heat::real> & Uin_Qin_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const auto & JxW = fe_v_self_U.get_JxW_values();
  const auto & normals 
    = fe_v_self_U.get_normal_vectors();
  
  for(unsigned int point = 0;
      point < fe_v_self_U.n_quadrature_points;
      ++point){    
    assert( referenceDirection * normals[point] > 0);
    for(unsigned int i = 0; i < fe_v_self_U.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_self_Q.dofs_per_cell; ++j){
	Uin_Qin_matrix(i,j) +=
	  JxW[point]
	  * fe_v_self_U.shape_value(i,point)
	  * fe_v_self_Q[flux].value(j,point)
	  * normals[point];
      }
    }
  }
}

  
  
template<int dim>
void
EllipticProblem<dim>
::assembleFluxPosQFromULocal
(const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> &fe_v_neig_U,
 dealii::FullMatrix<heat::real> & Qin_Uin_matrix,
 dealii::FullMatrix<heat::real> & Qin_Uout_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);

  const auto & JxW = fe_v_self_Q.get_JxW_values();
  const auto & normals = fe_v_self_Q.get_normal_vectors();
  
  for(unsigned int point = 0; 
      point < fe_v_self_Q.n_quadrature_points; 
      ++point){
    if(referenceDirection * normals[point] > 0){
      for(unsigned int i = 0; i < fe_v_self_Q.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_v_self_U.dofs_per_cell; ++j){
	  Qin_Uout_matrix(i,j) += 
	    JxW[point]
	    * fe_v_self_Q[flux].value(i,point)
	    * fe_v_neig_U[density].value(j,point)
	    * normals[point];				  
	}
      }
    }
  }  
}

template<int dim>
void
EllipticProblem<dim>
::assembleFluxNegQFromULocal
(const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> &fe_v_neig_U,
 dealii::FullMatrix<heat::real> & Qin_Uin_matrix,
 dealii::FullMatrix<heat::real> & Qin_Uout_matrix){
  
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);
  const auto & JxW = fe_v_self_Q.get_JxW_values();
  const auto & normals = fe_v_self_Q.get_normal_vectors();
  
  for(unsigned int point = 0; point < fe_v_self_Q.n_quadrature_points; ++point){
    if(referenceDirection * normals[point] < 0){
      for(unsigned int i = 0; i < fe_v_self_Q.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_v_neig_U.dofs_per_cell; ++j){
	  Qin_Uin_matrix(i,j) += 
	    JxW[point]
	    * fe_v_self_Q[flux].value(i,point)
	    * fe_v_self_U[density].value(j,point)
	    * normals[point];
	}
      }
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleFluxCenterQFromULocal
(const dealii::FEFaceValuesBase<dim> & fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
 const dealii::FEFaceValuesBase<dim> & fe_v_neig_U,
 dealii::FullMatrix<heat::real> & Qin_Uin_matrix,
 dealii::FullMatrix<heat::real> & Qin_Uout_matrix){
  
  const heat::real onehalf = 0.5;
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);
  const auto & JxW = fe_v_self_Q.get_JxW_values();
  const auto & normals 
    = fe_v_self_Q.get_normal_vectors();

  for(unsigned int point = 0; point < fe_v_self_Q.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_self_Q.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_self_U.dofs_per_cell; ++j){
	Qin_Uin_matrix(i,j) += 
	  onehalf
	  * JxW[point] 
	  * fe_v_self_Q[flux].value(i,point) 
	  * fe_v_self_U[density].value(j,point) 
	  * normals[point];
      }
    }
  }
  for(unsigned int point = 0; point < fe_v_self_Q.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_self_Q.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_neig_U.dofs_per_cell; ++j){
	Qin_Uout_matrix(i,j) += 
	  onehalf 
	  * JxW[point] 
	  * fe_v_self_Q[flux].value(i,point)
	  * fe_v_neig_U[density].value(j,point) 
	  * normals[point];
      }
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleFluxQFromUBoundaryLocal
(const dealii::FEFaceValuesBase<dim> & fe_v_self_Q,
 const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
 dealii::FullMatrix<heat::real> & Qin_Uin_matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const auto & JxW = fe_v_self_Q.get_JxW_values();
  const auto & normals = fe_v_self_Q.get_normal_vectors();
  for(unsigned int point = 0; point < fe_v_self_Q.n_quadrature_points; ++point){
      assert( referenceDirection * normals[point] < 0);
    for(unsigned int i = 0; i < fe_v_self_Q.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_self_U.dofs_per_cell; ++j){
	Qin_Uin_matrix(i,j) += 
	  JxW[point] 
	  * fe_v_self_Q[flux].value(i,point) 
	  * fe_v_self_U.shape_value(j,point) 
	  * normals[point];
      }
    }
  }
}
  



template<int dim>
void
EllipticProblem<dim>
::assembleTraceQLocalNeu
(const dealii::FEValuesBase<dim-1, dim> & fe_v_TraceU_Neu,
 const dealii::FEFaceValuesBase<dim> & fe_v_face_Q,
 dealii::FullMatrix<heat::real> & TrU_Neu_Q_Matrix)
{
  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);
  
  const std::vector<dealii::Point<dim> > & normals = fe_v_face_Q.get_normal_vectors();
  const std::vector<double> & JxW = fe_v_TraceU_Neu.get_JxW_values();

  for(unsigned int point = 0; point < fe_v_TraceU_Neu.n_quadrature_points;++point){
    //TODO:  Dealii assert.
    assert( referenceDirection * normals[point] > 0 );
    for(unsigned int i = 0; i < fe_v_TraceU_Neu.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_face_Q.dofs_per_cell; ++j){
	TrU_Neu_Q_Matrix(i,j) +=
	  JxW[point] 
	  * fe_v_TraceU_Neu[density].value(i,point)
	  * normals[point]
	  * fe_v_face_Q[flux].value(j,point);
      } 
    }
  }
}

// template<int dim>
// void
// EllipticProblem<dim>
// ::assembleTraceElecLocalNeu
// (
//  const dealii::FEValuesBase<dim-1,dim> & fe_v_TraceU_Neu,
//  const dealii::FEFaceValuesBase<dim> & fe_v_face_Elec,
//  dealii::FullMatrix<heat::real> & TrU_Neu_Elec_Matrix)
// {
//   const dealii::FEValuesExtractors::Scalar density(0);
//   const dealii::FEValuesExtractors::Vector elec(0);
//   const std::vector<dealii::Point<dim> > & normals = fe_v_face_Elec.get_normal_vectors();
//   const std::vector<double> & JxW = fe_v_face_Elec.get_JxW_values();
  
//   for(unsigned int point = 0; point < fe_v_TraceU_Neu.n_quadrature_points; ++point){
//     for(unsigned int i = 0; i < fe_v_TraceU_Neu.dofs_per_cell; ++i){
//       for(unsigned int j = 0; j < fe_v_face_Elec.dofs_per_cell; ++j){
// 	TrU_Neu_Elec_Matrix(i,j) +=
// 	  JxW[point]
// 	  * fe_v_TraceU_Neu[density].value(i,point)
// 	  * normals[point]
// 	  * fe_v_face_Elec[elec].value(j,point);	  
//       }
//     }
//   }
// }
  

template<int dim>
void
EllipticProblem<dim>
::assembleTraceULocalDir
(const dealii::FEValuesBase<dim-1,dim> & fe_v_self_TrU,
 const dealii::FEFaceValuesBase<dim> & fe_v_self_U,	      
 dealii::FullMatrix<heat::real> & TrU_U_matrix){
  const auto & JxW = fe_v_self_U.get_JxW_values();
  const auto & normals = fe_v_self_U.get_normal_vectors();
  for(unsigned int point = 0; point < fe_v_self_U.n_quadrature_points; ++point){
    assert( referenceDirection * normals[point] < 0);
    for(unsigned int i = 0; i < fe_v_self_TrU.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_v_self_U.dofs_per_cell; ++j){
	TrU_U_matrix(i,j) +=
	  JxW[point]
	  * fe_v_self_TrU.shape_value(i,point)
	  * fe_v_self_U.shape_value(j,point);
      }
    }
  }
}


template<int dim>
void
EllipticProblem<dim>
::assembleConstraintUDirMinusLocal
(const dealii::FEValuesBase<dim-1,dim> & fe_v_self_TrU,
 dealii::Vector<heat::real> & TrU_U_Vector){
  const heat::DirichletBoundaryValues<dim> dirBC;
  const auto & JxW = fe_v_self_TrU.get_JxW_values();  
  for(unsigned int point = 0; point < fe_v_self_TrU.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_self_TrU.dofs_per_cell; ++i){
      TrU_U_Vector(i) +=
	JxW[point]
	* fe_v_self_TrU.shape_value(i,point)
	* dirBC.value(fe_v_self_TrU.quadrature_point(point));	
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleConstraintUDirMinusLocalNOTRACE
(const dealii::FEFaceValuesBase<dim> & fe_v_face_U,
 dealii::Vector<heat::real> & U_Vector){
  const heat::DirichletBoundaryValues<dim> dirBC;
  const auto & JxW = fe_v_face_U.get_JxW_values();  
  for(unsigned int point = 0; point < fe_v_face_U.n_quadrature_points; ++point){
      for(unsigned int i = 0; i < fe_v_face_U.dofs_per_cell; ++i){
      U_Vector(i) +=
	JxW[point]
	* fe_v_face_U.shape_value(i,point)
	* dirBC.value(fe_v_face_U.quadrature_point(point));	
    }
  }
}


template<int dim>
void
EllipticProblem<dim>
::assembleConstraintQNeuPlusLocal
(const dealii::FEValuesBase<dim-1,dim> & fe_v_self_TrU,
 const dealii::FEFaceValuesBase<dim> & fe_v_face_Q,
 dealii::Vector<heat::real> & TrU_Q_Vector){
  const heat::NeumannBoundaryValues<dim> nueBC;
  const auto & JxW = fe_v_self_TrU.get_JxW_values();
  const auto & normals = fe_v_face_Q.get_normal_vectors();
  for(unsigned int point = 0; point < fe_v_self_TrU.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_v_self_TrU.dofs_per_cell; ++i){
      //dealii::Vector<heat::real> temp(dim);
      //nueBC.vector_value(fe_v_self_TrU.quadrature_point(point), temp);
  
      //dealii::Tensor<1,dim,heat::real> tempTensor;
      auto something = nueBC.value(fe_v_self_TrU.quadrature_point(point));
      //tempTensor = temp;
      
      TrU_Q_Vector(i) +=
	JxW[point]
	* fe_v_self_TrU.shape_value(i,point)
	* (something * normals[point]);
    }
  }
}


template<int dim>
void
EllipticProblem<dim>
  ::assembleQ_RHS_FromDirPlusLocal
(
      const dealii::FEFaceValuesBase<dim> & fe_v_face_Q,
	dealii::Vector<heat::real> & Q_Vector)
{

      const heat::DirichletBoundaryValues<dim> dirBC;


      const dealii::FEValuesExtractors::Vector flux(0);
      const auto & normals = fe_v_face_Q.get_normal_vectors();
      const auto & JxW = fe_v_face_Q.get_JxW_values();
      for(unsigned int point = 0; point < fe_v_face_Q.n_quadrature_points; ++point){
      //TODO:  Dealii assert.
      assert( referenceDirection * normals[point] > 0.0);
      for(unsigned int i = 0; i < fe_v_face_Q.dofs_per_cell; ++i){
      Q_Vector(i) -= JxW[point]
	* fe_v_face_Q[flux].value(i,point)
	* normals[point]
	* dirBC.value(fe_v_face_Q.quadrature_point(point));
    }    
  }  
}

template<int dim>
void
EllipticProblem<dim>
::assembleQ_RHS_FromNeuMinusLocal
(
 const dealii::FEFaceValuesBase<dim> & fe_v_face_U,
 dealii::Vector<heat::real> & U_Vector)
{
  const heat::NeumannBoundaryValues<dim> neuBC;
  const dealii::FEValuesExtractors::Scalar density(0);
  const auto & normals = fe_v_face_U.get_normal_vectors();
  const auto & JxW = fe_v_face_U.get_JxW_values();
  for(unsigned int point = 0; point < fe_v_face_U.n_quadrature_points; ++point){
    //TODO:  Dealii assert
    assert( referenceDirection * normals[point] < 0.0);
    for(unsigned int i = 0; i < fe_v_face_U.dofs_per_cell; ++i){
      U_Vector(i) +=
	JxW[point]
	* fe_v_face_U[density].value(i,point)
	* neuBC.value(fe_v_face_U.quadrature_point(point)) * normals[point];
    }
  }
}

template<int dim>
void
EllipticProblem<dim>
::assembleU_RHS_FromSourceLocal
( const dealii::FEValuesBase<dim> & fe_v_u,
  dealii::Vector<heat::real> & U_Vector
 )
{
  const heat::SourceBodyValues<dim> sourceBV;
  const dealii::FEValuesExtractors::Scalar density(0);
  const auto & JxW = fe_v_u.get_JxW_values();
  for(unsigned int point = 0; point < fe_v_u.n_quadrature_points; ++point){
      for(unsigned int i = 0; i < fe_v_u.dofs_per_cell; ++i){
      U_Vector(i) += JxW[point] * fe_v_u[density].value(i,point)
	* sourceBV.value(fe_v_u.quadrature_point(point));
    }
  }
}


// template<int dim>
// void
// EllipticProblem<dim>
// ::assembleTraceUElecDir
// (const dealii::FEValuesBase<dim-1,dim> & fe_v_self_TrU,
//  const dealii::FEFaceValuesBase<dim> & fe_v_face_E,
//  dealii::FullMatrix<heat::real> & TrU_Dir_Elec_Matrix)
// {
//   const dealii::FEValuesExtractors::Vector elec(0);
//   const dealii::FEValuesExtractors::Scalar density(0);  
//   const std::vector<dealii::Point<dim> > & normals = fe_v_face_E.get_normal_vectors();  
//   const std::vector<double> & JxW = fe_v_self_TrU.get_JxW_values();
  
//   for(unsigned int point = 0; point < fe_v_self_TrU.n_quadrature_points; ++point){
//     for(unsigned int i = 0; i < fe_v_self_TrU.dofs_per_cell; ++i){
//       for(unsigned int j = 0; j < fe_v_face_E.dofs_per_cell; ++j){
// 	TrU_Dir_Elec_Matrix(i,j) += 
// 	  JxW[point]
// 	  * fe_v_self_TrU[density].value(i,point)
// 	  * fe_v_face_E[elec].value(j,point)
// 	  * normals[point];	      
//       }
//     }
//   }
// }
 
 

// template<int dim>
// void
// EllipticProblem<dim>
// ::assembleTraceQLocal
// (const dealii::FEFaceValuesBase<dim> & fe_v_self_Q,
//  const dealii::FEFaceValuesBase<dim> & fe_v_self_Phi,
//  dealii::FullMatrix<heat::real> & Q_Phi_matrix)
// {
//   const dealii::FEValuesExtractors::Vector flux(0);
//   const auto & JxW = fe_v_self_Q.get_JxW_values();
//   const auto & normals = fe_v_self_Q.get_normal_vectors();  

//   for(unsigned int point = 0; 
//       point < fe_v_self_Q.n_quadrature_points;
//       ++point){
//     for(unsigned int i = 0; i < fe_v_self_Q.dofs_per_cell; ++i){
//       for(unsigned int j = 0; j < fe_v_self_Phi.dofs_per_cell; ++j){
// 	Q_Phi_matrix(i,j) +=
// 	  JxW[point]
// 	  * fe_v_self_Q[flux].value(i,point)
// 	  * fe_v_self_Phi.shape_value(j,point)
// 	  * normals[point];
//       }
//     }
//   }
// }

//  TODO:  Need to find the correct spot for this.  
//  Probably shouldn't be part of the heat problem class.

template<int dim>
void
EllipticProblem<dim>
::InsertGlobalFromLocal
(const typename dealii::DoFHandler<dim-1,dim>::active_cell_iterator & test_iterator,
 const typename dealii::DoFHandler<dim, dim>::active_cell_iterator & tria_iterator,
 const dealii::FullMatrix<heat::real> & local_matrix,
 dealii::SparseMatrixEZ<heat::real> & global_matrix)
{

  heat::real tolerance = 1e-10;
  const unsigned int dofs_per_cell_test = test_iterator->get_fe().dofs_per_cell;
  const unsigned int dofs_per_cell_tria = tria_iterator->get_fe().dofs_per_cell;
  
  std::vector<dealii::types::global_dof_index> dofs_test(dofs_per_cell_test);
  std::vector<dealii::types::global_dof_index> dofs_tria(dofs_per_cell_tria);

  test_iterator->get_dof_indices(dofs_test);
  tria_iterator->get_dof_indices(dofs_tria);

  for(unsigned int i = 0; i < dofs_per_cell_test; ++i){
    for(unsigned int j = 0; j < dofs_per_cell_tria; ++j){
      const heat::real temp = local_matrix(i,j);
      if(fabs(temp) > tolerance){
	global_matrix.add(dofs_test[i],dofs_tria[j],temp);
      }
    }
  }  
}

template<int dim>
void
EllipticProblem<dim>
::InsertGlobalFromLocal
(const typename dealii::DoFHandler<dim>::active_cell_iterator & test_iterator,
 const typename dealii::DoFHandler<dim>::active_cell_iterator & tria_iterator,
 const dealii::FullMatrix<heat::real> & local_matrix,
 dealii::SparseMatrix<heat::real> & global_matrix)
{
  heat::real tolerance = 1e-10;;
  const unsigned int dofs_per_cell_test = test_iterator->get_fe().dofs_per_cell;
  const unsigned int dofs_per_cell_tria = tria_iterator->get_fe().dofs_per_cell;

  std::vector<dealii::types::global_dof_index> dofs_test(dofs_per_cell_test);
  std::vector<dealii::types::global_dof_index> dofs_tria(dofs_per_cell_tria);

  test_iterator->get_dof_indices(dofs_test);
  tria_iterator->get_dof_indices(dofs_tria);

  for(unsigned int i = 0; i < dofs_per_cell_test; ++i){
    for(unsigned int j = 0; j < dofs_per_cell_tria; ++j){
      const heat::real temp = local_matrix(i,j);
      if(fabs(temp) > tolerance){
	global_matrix.add(dofs_test[i],dofs_tria[j],temp);
      }
    }
  }
}


template<int dim>
void
EllipticProblem<dim>
::InsertGlobalFromLocalVector
(const typename dealii::DoFHandler<dim>::active_cell_iterator & test_iterator,
 const dealii::Vector<heat::real> & local_vector,
 dealii::Vector<heat::real> & global_vector)
{
  heat::real tolerance = 1e-10;;
  const unsigned int dofs_per_cell_test = test_iterator->get_fe().dofs_per_cell;

  std::vector<dealii::types::global_dof_index> dofs_test(dofs_per_cell_test);

  test_iterator->get_dof_indices(dofs_test);

  for(unsigned int i = 0; i < dofs_per_cell_test; ++i){
    const heat::real temp = local_vector(i);
    if(fabs(temp) > tolerance){
      global_vector(dofs_test[i]) += temp;
    }
  }
}


//TODO:  This is almost exactly the same as the previous function.  Find out how to combine them.
// My guess is that they shouldn't be defined inside of the elliptic problem class.
template<int dim>
void
EllipticProblem<dim>
::InsertGlobalFromLocalVectorNEW
(const typename dealii::DoFHandler<dim-1, dim>::active_cell_iterator & test_iterator,
 const dealii::Vector<heat::real> & local_vector,
 dealii::Vector<heat::real> & global_vector)
{
  heat::real tolerance = 1e-10;;
  const unsigned int dofs_per_cell_test = test_iterator->get_fe().dofs_per_cell;

  std::vector<dealii::types::global_dof_index> dofs_test(dofs_per_cell_test);

  test_iterator->get_dof_indices(dofs_test);

  for(unsigned int i = 0; i < dofs_per_cell_test; ++i){
    const heat::real temp = local_vector(i);
    if(fabs(temp) > tolerance){
      global_vector(dofs_test[i]) += temp;
    }
  }
}

template<int dim>
void
EllipticProblem<dim>::EigenSolve(){  
  
  //TODO:  Assert here to make sure this is ok.  That we actually can get this number of
  // eigenvalues
  

  
  std::vector<dealii::Vector<heat::real> > reducedUEigenfunctions(n_eigenvalues);

  for(auto&& i : reducedUEigenfunctions){
    i.reinit(QQ.m());
  }
  
  std::vector<std::complex<double> > eigenvalues(n_eigenvalues);
  Q_UStar_SStar_inverseMassQ_S_U_Q 
    bigProduct
    (QQ,
     Svd_U,
     TotalQFromU,
     InverseMassQ,
     TotalUFromQ
     );

  const Q_UStar_massU_U_QStar smallProduct(QQ,Svd_U, MassU);
  

  auto Identity = dealii::PreconditionIdentity();
  dealii::SolverControl innerControl(1000,1e-14);

  dealii::IterativeInverse<dealii::Vector<heat::real>> bigProductInverse;
  bigProductInverse.initialize(bigProduct, Identity);
  bigProductInverse.solver.select("cg");
  bigProductInverse.solver.set_control(innerControl);

  
  
  
  const unsigned int num_arnoldi_vectors = 2 * n_eigenvalues + 2;
  dealii::ArpackSolver::AdditionalData additional_data
    (num_arnoldi_vectors,
     dealii::ArpackSolver::WhichEigenvalues::largest_real_part
     );

  dealii::SolverControl solver_control(dof_handler.n_dofs(), 1e-14);
  
  dealii::ArpackSolver eigenSolver(solver_control, additional_data);
  eigenSolver.solve(
		    bigProduct,
		    smallProduct,
		    bigProductInverse,
		    eigenvalues,
		    reducedUEigenfunctions,
		    eigenvalues.size() );
  
  
  for(unsigned int i = 0; i < eigenvalues.size(); ++i){
    dealii::Vector<heat::real> temp1( QQ.n() );
    dealii::Vector<heat::real> temp2( InverseMassQ.m() );
    QQ.Tvmult(temp1, reducedUEigenfunctions[i]);
    Svd_U.vmult(UeigenStates[i], temp1);
    
    TotalQFromU.vmult(temp2, UeigenStates[i]);
    InverseMassQ.vmult(QeigenStates[i], temp2);    
    QeigenStates[i] *= -1.0;
  }  
}

template <int dim>
void EllipticProblem<dim>::output_results_eigensystem ()
{
  typedef dealii::DataComponentInterpretation::DataComponentInterpretation DCI_type;
    
  dealii::DataOut<dim> dataOutSys;      
  dataOutSys.attach_dof_handler (dof_handler);

  std::vector<dealii::BlockVector<heat::real> > EigenStates(n_eigenvalues);

  for (unsigned int i=0; i<UeigenStates.size(); ++i){

    EigenStates[i].reinit(2);
    EigenStates[i].block(0).reinit(UeigenStates[i].size());
    EigenStates[i].block(1).reinit(QeigenStates[i].size());
    EigenStates[i].collect_sizes();
    
    EigenStates[i].block(0) = UeigenStates[i];
    EigenStates[i].block(1) = QeigenStates[i];    
    

    
    std::vector<std::string> SolutionNames;
    std::vector<DCI_type> DataComponentInterpretation;

    SolutionNames.emplace_back("Density_" + dealii::Utilities::int_to_string(i,2) );	
    DataComponentInterpretation
      .push_back(dealii
		 ::DataComponentInterpretation
		 ::component_is_scalar);

      
    for(unsigned int j = 0; j < dim; ++j){
      SolutionNames.emplace_back("DiffusiveCurrent" + dealii::Utilities::int_to_string(j,1)
				 + "_" + dealii::Utilities::int_to_string(i,2) );
      DataComponentInterpretation
	.push_back(dealii::DataComponentInterpretation::component_is_scalar);
    }

    
    dataOutSys.add_data_vector(EigenStates[i],
			       SolutionNames,
			       dealii::DataOut<dim>::type_dof_data,
			       DataComponentInterpretation);
      
  }    
  dataOutSys.build_patches (degree*4);
  std::ofstream output ("eigenvectors.vtk");
  dataOutSys.write_vtk (output);
}



template<int dim>
void
EllipticProblem<dim>::ElliptSolve(){
  //grvy_timer_begin("find Uperp");
  
      //Find Uperp;
      dealii::Vector<heat::real> uPerp(rankBoundaryMassU);
      UperpNEW.reinit(rankBoundaryMassU);
  {
      //Ustate = 0;
      PSigmaPStar PSigmaPStar_Instance(PP,Svd_Sigma);
    
    dealii::Vector<heat::real>
      temp1(Svd_Sigma.n()),
      temp2(PP.m()),
      temp3(PP.n());
   
    Svd_U.Tvmult(temp1, UdirConstraint_RHS);
    PP.vmult(temp2,temp1);

    dealii::SolverControl SolverControl(10000,1e-14);
    dealii::SolverCG<> cg(SolverControl);


    cg.solve(PSigmaPStar_Instance, UperpNEW, temp2,
	     dealii::PreconditionIdentity() );

    
    PP.Tvmult(temp3,UperpNEW);
    Svd_U.vmult(Ustate, temp3);    
    
  }

  
  //Find Upar;
  dealii::Vector<heat::real> uPar(nullityBoundaryMassU);
  UparNEW.reinit(nullityBoundaryMassU);
  {
    dealii::Vector<heat::real>
      temp1( rankBoundaryMassU + nullityBoundaryMassU ),
      temp2( Svd_U.m() ),
      temp3( TotalQFromU.m() ),
      temp4( TotalQFromU.m() ),
      temp5( TotalUFromQ.m() ),
      temp6( Svd_U.n() ),
      temp7( nullityBoundaryMassU ),
      temp8( QQ.n() ),
      temp9( Svd_U.m() );      
      
    //grvy_timer_begin("uPar RHS");
    PP.Tvmult(temp1, UperpNEW);    
    Svd_U.vmult(temp2,temp1);    
    TotalQFromU.vmult(temp3, temp2);
    temp3 += Q_MinusRHS;
    InverseMassQ.vmult(temp4, temp3);
    TotalQFromU.Tvmult(temp5,temp4);
    temp5 += U_MinusRHS;
    Svd_U.Tvmult(temp6, temp5);
    QQ.vmult(temp7, temp6);
    temp7 *= -1.0;     
    //grvy_timer_end("uPar RHS");


    
    Q_UStar_SStar_inverseMassQ_S_U_Q 
      bigProduct
      (QQ,
	Svd_U,
	TotalQFromU,
	InverseMassQ,
	TotalUFromQ
	); 

    dealii::SolverControl SolverControl(10000, 1e-16);
    dealii::SolverCG<>    cg(SolverControl);

    //grvy_timer_begin("Upar Solve");
    cg.solve(bigProduct, UparNEW, temp7, dealii::PreconditionIdentity() );
    //grvy_timer_end("Upar Solve");

    
    
    QQ.Tvmult(temp8, UparNEW);
    Svd_U.vmult(temp9,temp8);
    Ustate += temp9;    

    
  }  
  //grvy_timer_end("find Upar");
  //Compute Q;
  //grvy_timer_begin("compute Q");
  {
    dealii::Vector<heat::real>
      temp1(TotalQFromU.m() );
    temp1 = 0;
    TotalQFromU.vmult(temp1,Ustate);
    temp1 += Q_MinusRHS;
    temp1 *= -1.0;
    InverseMassQ.vmult(Qstate,temp1);
      
  }
  //grvy_timer_end("compute Q");

  
  //Compute lambda (don't really care though);
  //grvy_timer_begin("compute lambda");
  {
    dealii::Vector<heat::real>
      temp1( TotalUFromQ.m() ),
      temp2( Svd_U.m() ),
      temp3( PP.m() ),
      lambdaHat( PP.m() ),
      temp4( PP.n() );
    

    TotalUFromQ.vmult( temp1, Qstate );
    temp1 += U_MinusRHS;
    Svd_U.Tvmult( temp2, temp1 );
    PP.vmult( temp3, temp2);


    PSigmaPStar PSigmaPStar_Instance(PP,Svd_Sigma);

    dealii::SolverControl SolverControl(10000,1e-16);
    dealii::SolverCG<> cg(SolverControl);

    cg.solve(PSigmaPStar_Instance, lambdaHat, temp3,
	     dealii::PreconditionIdentity() );

    PP.Tvmult(temp4, lambdaHat);
    Svd_U.vmult(lambdaU,temp4);
    
  }
  //grvy_timer_end("compute lambda");
}

// template<int dim>
// void
// EllipticProblem<dim>::USolve()
// {
//   dealii::deallog.push("USolve");
//   // QSolve and USolve are almost identical,
//   // We should refactor this.
//   dealii::PreconditionIdentity identity;
//   // dealii::PreconditionJacobi<dealii::SparseMatrix< heat::real > > jacobiPre;
//   // jacobiPre.initialize(MassU);
//   dealii::IterativeInverse<dealii::Vector< heat::real > > MassUInverse;
//   MassUInverse.initialize(MassU, InverseMassU);
//   MassUInverse.solver.select("cg");
  
  
//   //TODO:  Magic NumberU  
//   static dealii::SolverControl innerControl(5, 1e-14);
//   MassUInverse.solver.set_control(innerControl);  
//   dealii::Vector<heat::real> temp(Ustate.size() );
//   {    
//     dealii::Vector<heat::real> schur_rhs(constraintDirU.size() );
//     try {
//       MassUInverse.vmult(temp, U_RHS);
//     } 
//     catch (...) {
//       std::cerr << "Failure in the first Usolve inversion" << std::endl;
//       std::abort();
//     }
    
//     TraceUDir.vmult(schur_rhs, temp);
//     schur_rhs -= constraintDirU;    
//     heat::SchurComplement
//       schurComplement(TraceUDir,MassUInverse);    
//     heat::ApproximateSchurComplement
//       approx_schurcomplement(MassU,TraceUDir);    
//     dealii::IterativeInverse<dealii::Vector<heat::real> >
//       preconditionInner;
    
//     preconditionInner.initialize(approx_schurcomplement, identity);
//     preconditionInner.solver.select("cg");
//     preconditionInner.solver.set_control(innerControl);

//     //TODO:  Magic Numbers
//     dealii::SolverControl solver_control( 40 , 1e-8);
//     dealii::SolverCG<> cg(solver_control); 
    
//     try {
//       cg.solve(schurComplement, lagrangeU, schur_rhs, identity);
//     } 
//     catch (...) {
//       std::cerr << "Failure in the USolve Schur Complement Solver" << std::endl;
//       std::abort();
//     }
//   }
  
//   TraceUDir.Tvmult(temp,lagrangeU);
//   temp *= -1;
//   temp += U_RHS;
//   try {
//     MassUInverse.vmult(UStateDot, temp);
//   }
//   catch (...){
//     std::cerr << "Failure in the second Usolve inversion" << std::endl;
//     std::abort();
//   }  
//   dealii::deallog.pop();
// }

/*
template<int dim>
void
EllipticProblem<dim>::USolve_NEW(dealii::Vector<heat::real> & Udest,
			     dealii::Vector<heat::real> & lambdadest,
			     const dealii::Vector<heat::real> & source1,
			     const dealii::Vector<heat::real> & source2)
{
  dealii::deallog.push("USolve");

  dealii::PreconditionIdentity identity;

  dealii::IterativeInverse<dealii::Vector< heat::real > > MassUInverse;
  MassUInverse.initialize(MassU, InverseMassU);
  MassUInverse.solver.select("cg");
  
  
  //TODO:  Magic NumberU  
  static dealii::SolverControl innerControl(5, 1e-14);
  MassUInverse.solver.set_control(innerControl);  
  dealii::Vector<heat::real> temp(Ustate.size() );
  {    
    dealii::Vector<heat::real> schur_rhs(constraintDirU.size() );
    try {
      MassUInverse.vmult(temp, source1);
    } 
    catch (...) {
      std::cerr << "Failure in the first Usolve inversion" << std::endl;
      std::abort();
    }
    
    TraceUDir.vmult(schur_rhs, temp);
    schur_rhs -= source2;    
    heat::SchurComplement
      schurComplement(TraceUDir,MassUInverse);    
    heat::ApproximateSchurComplement
      approx_schurcomplement(MassU,TraceUDir);    
    dealii::IterativeInverse<dealii::Vector<heat::real> >
      preconditionInner;
    
    preconditionInner.initialize(approx_schurcomplement, identity);
    preconditionInner.solver.select("cg");
    preconditionInner.solver.set_control(innerControl);

    //TODO:  Magic Numbers
    dealii::SolverControl solver_control( 40 , 1e-8);
    dealii::SolverCG<> cg(solver_control); 
    
    try {
      cg.solve(schurComplement, lambdadest, schur_rhs, identity);
    } 
    catch (...) {
      std::cerr << "Failure in the USolve Schur Complement Solver" << std::endl;
      std::abort();
    }
  }
  
  TraceUDir.Tvmult(temp,lambdadest);
  temp *= -1;
  temp += U_RHS;
  try {
    MassUInverse.vmult(Udest, temp);
  }
  catch (...){
    std::cerr << "Failure in the second Usolve inversion" << std::endl;
    std::abort();
  }  
  dealii::deallog.pop();
}

template<int dim>
void
EllipticProblem<dim>::PoissonSolve()
{
  //TODO:  This is confusing, write up a TeX file explaining this.  
  dealii::PreconditionIdentity identity;
  
  dealii::PreconditionJacobi<dealii::SparseMatrix<heat::real> > jacobiPre;
  jacobiPre.initialize(MassElec);
  dealii::IterativeInverse<dealii::Vector< heat::real> > MassElecInverse;
  MassElecInverse.initialize(MassElec, MassElecInverseDirect);
  MassElecInverse.solver.select("cg");
  
  //TODO:  Magic Numbers  
  static dealii::ReductionControl innerControl(2000,0, 1.0e-10);
  MassElecInverse.solver.set_control(innerControl);
  
  dealii::Vector<heat::real> temp1( ElecState.size() );
  dealii::Vector<heat::real> temp2( PotState.size() );
  dealii::Vector<heat::real> temp3( TraceElecNeu.m() );

  std::vector<dealii::BlockVector<heat::real>::size_type> block_sizes;  
  block_sizes.push_back( PotState.size() );
  block_sizes.push_back( TraceElecNeu.m() );

  
  
  dealii::BlockVector<heat::real> schur_rhs(block_sizes);
  schur_rhs.collect_sizes();
  dealii::BlockVector<heat::real> answer = schur_rhs;
  
  //TODO:  Figure out a way to use the previous state as an initial guess.
  {
    //Build the RHS's
    // First we will concentrate on the first block
    
    try {
      MassElecInverse.vmult(temp1, Elec_RHS);
    } catch (...){
      std::cerr << "Failure in the first PoissonSolver Inversion" << std::endl;
      std::abort();
    }
    StiffPotFromElec.vmult(temp2,temp1);


    PotFromU.Tvmult(schur_rhs.block(0), Ustate);
    //std::cout << "Turn off Coupling of poisson depending on U" << std::endl;
    
    schur_rhs.block(0) = 0;
    schur_rhs.block(0) -= temp2;
    
    // Now the second block.
    // We reuse temp1,    
    TraceElecNeu.vmult(temp3,temp1);
    schur_rhs.block(1)  = PoissonNeuBC_RHS;
    schur_rhs.block(1) -= temp3;    
    
    //Now the schur_rhs is properly populated    
    heat::SchurComplementPoisson
      schurComplement(TraceElecNeu,
		      StiffElecFromPot,
		      MassElecInverse);
    
    //schurComplement.vmult(answer, schur_rhs);    
    heat::ApproximateSchurComplementPoisson
      approximateSchurComplementPoisson(TraceElecNeu,
					StiffElecFromPot,
					MassElec);

    dealii::IterativeInverse<dealii::BlockVector<heat::real> > precondition;
    precondition.initialize(approximateSchurComplementPoisson, identity);
    precondition.solver.set_control(innerControl);
    precondition.solver.select("cg");
    
    //TODO:  Magic Numbers;
    dealii::SolverControl solver_control(2000,1.0e-7);
    dealii::SolverCG<dealii::BlockVector<heat::real> > cg(solver_control);
    
    try {
      cg.solve(schurComplement, answer, schur_rhs, identity);    
    } catch (...) {
      std::cerr << "Failure in the SchurComplementPoisson Solver" << std::endl;
      std::abort();
    }
    
    {
      dealii::BlockVector<heat::real> residual = schur_rhs;
      schurComplement.vmult(residual, answer);
      residual -= schur_rhs;
      schurComplement.vmult(schur_rhs, residual);
    }
  }
  
  // dealii::Vector<heat::real> 
  //   temp4;
  StiffElecFromPot.vmult(temp1, answer.block(0));
  TraceElecNeu.Tvmult_add(temp1, answer.block(1));
  temp1 += Elec_RHS;    
  PotState = answer.block(0);  
  try {
    MassElecInverse.vmult(ElecState, temp1);  
  } 
  catch (...) {
    std::cout << "Failure in PoissonSolve second inversion" << std::endl;
    std::abort(); 
  }
}

template<int dim>
void
EllipticProblem<dim>::simulate()
{
  //It seems like this should be in a different spot.
  const double h 
    = dealii::GridTools::minimal_cell_diameter(triangulation);
  
  std::cout << "h = " << h << std::endl;

  //Todo: Magic Number 0.1  
  heat::real deltaT = .005 * h * h / ((2 * degree +1) *  (2 * degree + 1));

  //Todo: Magic Number numTimeStamps.
  unsigned int numTimeStamps = 100;
  
  double 
    timeInitial = 0.0,
    timeFinal   = 0.1;
  
  std::vector<heat::real> timeStampVector(numTimeStamps);
  unsigned int numTimeSteps = numTimeStamps - 1;
  heat::real timePerTimeStep = (timeFinal - timeInitial) / numTimeSteps;
  
  for(unsigned int i = 0; i < timeStampVector.size(); ++i){
    timeStampVector[i] = timeInitial + i * timePerTimeStep;
  }
  
  //assert(timeStampVector[timeStampVector.size() - 1] == timeFinal);

  heat::real timeNow;
  heat::real timeNext;

  // We start at one and not zero because we printed the initial condition
  // already.
  for(unsigned int timeStampLabel = 1;  
      timeStampLabel < numTimeStamps;
      ++timeStampLabel){
    
    timeNow = timeStampVector[timeStampLabel - 1];
    timeNext = timeStampVector[timeStampLabel];
    
    assert(timeNow < timeNext);
    ComputeUnext(timeNow,timeNext,deltaT);

    std::cout << "\rPrinting timeStampLabel " << timeStampLabel << " of " 
	      << numTimeStamps - 1 << std::endl;    
    
    ComputeDensityDot(timeNext);
    //assert(timeStampLabel < 3);
    PrintState  (timeStampLabel);
    
  }
  std::cout << "\nFinished" << std::endl;
}

template<int dim>
void
EllipticProblem<dim>::ComputeDensityDot(const heat::real time)
{
  makeElec_RHS_From_DirichletBoundary(time);
  makeConstraint_RHS_From_NeumanBoundary();
  
  PoissonSolve();
  makeConstraintNeuQ(time);
  
  makeQ_RHS_FromU();
  
  
  QSolve_NEW(Qstate,lagrangeQ,Q_RHS, constraintNeuQ);
  


  //Compute Udot  
  makeConstraintDirU(time);
  
  makeU_RHS_FromQ();
  //makeU_RHS_From_U();
  //U_RHS += U_RHS_From_U;
  //std::cout << "U_RHS.l2_norm() = " << U_RHS.l2_norm() << std::endl;
  
  USolve_NEW(UStateDot,lagrangeU,U_RHS,constraintDirU);

  
}

template<int dim>
void
EllipticProblem<dim>::ComputeUnext
(heat::real timeNow,
 heat::real timeNext,
 heat::real deltaT)
{
  heat::real deltaTRecommended = deltaT;
  assert(timeNow < timeNext);
  
  bool acceptFlag = false;

  //TODO Magic Number
  assert(timeNext - timeNow > 0);    
  
  while (timeNext - timeNow > 0.0000000001){  
    //std::cout << "timeNow = " << timeNow << "\n";
    if(timeNow + deltaTRecommended > timeNext){
      deltaT = timeNext - timeNow;
    }
    else{
      deltaT = deltaTRecommended;
    }
    
    //std::cout << "timeNow = " << timeNow << std::endl;

    //TODO: Magic number
    //assert(fabs(timeNow + deltaT-timeNext) < 0.000001);
    
    ComputeDensityDot(timeNow);
    ForwardEuler(timeNow,deltaT,acceptFlag,deltaTRecommended);
    
    // AdaptiveExplicitRungeKutta54
    //       (timeNow,
    // 	   deltaT,
    // 	   acceptFlag,
    // 	   deltaTRecommended);

    if(acceptFlag){
      Ustate += DeltaUState;
      timeNow += deltaT;
    }
  }  
}

template<int dim>
void
EllipticProblem<dim>::ForwardEuler
(const heat::real timeNow,
 const heat::real deltaT,
 bool & acceptFlag,
 heat::real & deltaTRecommended)
{
  Assert( deltaT > 0, dealii::ExcMessage("Negative DeltaT") );
  
  DeltaUState.sadd(0.0, deltaT, UStateDot);
  acceptFlag = true;
  deltaTRecommended = deltaT;
}

				
template<int dim>
void
EllipticProblem<dim>::AdaptiveExplicitRungeKutta54
(const heat::real timeNow,
 const heat::real deltaT,
 bool & acceptFlag,
 heat::real & deltaTRecommended)
{
  //  This isn't really ideal, but compute density dot does its work
  //  using the vector that resides in Ustate;  And we need to feed different,
  //  stage vectors in there.  So, for now, we just store the state in this temp
  //  Storage vector.
  
  dealii::Vector<heat::real> TempStorage = Ustate;
  const unsigned int numDof = Ustate.size();
  const int stages = 7;

  dealii::FullMatrix<heat::real> ButcherA(stages, stages);
  dealii::FullMatrix<heat::real> ButcherB(2,stages);
  dealii::Vector<heat::real> ButcherC(stages);

  ButcherA = 0;
  ButcherB = 0;
  ButcherC = 0;

  ButcherA(1,0) = 1.0 / 5.0;
  
  ButcherA(2,0) = 3.0 / 40.0;
  ButcherA(2,1) = 9.0 / 40.0;
  
  ButcherA(3,0) =  44.0 / 45.0;
  ButcherA(3,1) = -56.0 / 15.0;
  ButcherA(3,2) =  32.0 / 9.0;
  
  ButcherA(4,0) =  19372.0 / 6561.0;
  ButcherA(4,1) = -25360.0 / 2187.0;
  ButcherA(4,2) =  64448.0 / 6561.0;
  ButcherA(4,3) = -212.0   / 729.0;
  
  ButcherA(5,0) =  9017.0 / 3168.0;
  ButcherA(5,1) = -355.0 / 33.0;
  ButcherA(5,2) =  46732.0 / 5247.0;
  ButcherA(5,3) =  49.0 / 176.0;
  ButcherA(5,4) = -5103.0 / 18656.0;
    
  ButcherA(6,0) =  35.0 / 384.0;
  ButcherA(6,1) =  0.0;
  ButcherA(6,2) =  500.0 / 1113;
  ButcherA(6,3) =  125.0 / 192.0;
  ButcherA(6,4) = -2187.0 / 6784.0;
  ButcherA(6,5) =  11.0 / 84.0;
  
  ButcherB(0,0) =  35.0 / 384.0;
  ButcherB(0,1) =  0.0;
  ButcherB(0,2) =  500.0 / 1113.0;
  ButcherB(0,3) =  125.0 / 192.0;
  ButcherB(0,4) = -2187.0 / 6784.0;
  ButcherB(0,5) =  11.0 / 84.0;
  ButcherB(0,6) =  0.0;

  ButcherB(1,0) =  5179.0 / 57600.0;
  ButcherB(1,1) =  0.0;
  ButcherB(1,2) =  7571.0 / 16695.0;
  ButcherB(1,3) =  393.0 / 640.0;
  ButcherB(1,4) = -92097.0 / 339200.0;
  ButcherB(1,5) =  187.0 / 2100.0;
  ButcherB(1,6) =  1.0 / 40.0;

  ButcherC(0) = 0.0;
  ButcherC(1) = 1.0 / 5.0;
  ButcherC(2) = 3.0 / 10.0;
  ButcherC(3) = 4.0 / 5.0;
  ButcherC(4) = 8.0 / 9.0;
  ButcherC(5) = 1.0;
  ButcherC(6) = 1.0;

  typedef std::vector< dealii::Vector<heat::real> > workspaceVectors;  
  dealii::Vector<heat::real> fourthOrder(numDof);
  dealii::Vector<heat::real> fifthOrder(numDof);
  fourthOrder = 0;
  fifthOrder = 0;
  workspaceVectors k(stages, dealii::Vector<heat::real>(numDof ));

  for(int stageNumber = 0; stageNumber < stages; ++stageNumber){
      //compute k[j]
      //  compute stage time
      const heat::real stageTime = timeNow + deltaT * ButcherC[stageNumber];
      //  compute stage Density;
      //dealii::Vector<heat::real> stageState( numDof );
      //stageState = Ustate;
      Ustate = TempStorage;
      for ( int j = 0; j < stageNumber-1; ++j){
	Ustate.add( ButcherA(stageNumber, j), k[j]);
      }

    // dealii::Vector<heat::real> stageDensityDot( numDof );
    // stageDensityDot = 0;
      
    ComputeDensityDot(stageTime);      			
      k[stageNumber] = UStateDot;
      k[stageNumber] *= deltaT;            
  }
  
  for(int j = 0; j < stages; ++j){
   
    heat::real b = ButcherB(0,j);
    fifthOrder.add( b , k[j] );
  }

  for(int j = 0; j < stages; ++j){
    heat::real b = ButcherB(1,j);
    fourthOrder.add( b, k[j] );
  }
  
  dealii::Vector<heat::real> difference;
  difference = fourthOrder;
  difference -=fifthOrder;

  const heat::real errorSquared = MassU.matrix_norm_square(difference);
  const heat::real error = sqrt(errorSquared);  
  const heat::real toleranceMultiplierSquared = MassU.matrix_norm_square(Ustate);
  const heat::real toleranceMultiplier = sqrt(toleranceMultiplierSquared);
  heat::real tolerance = toleranceMultiplier*.1;
  // tolerance = std::max(tolerance, 4.0);
  //  const heat::real beta = 1.0;
  heat::real factor = NAN;

  std::cout << "error = " << error << std::endl;  
  
  if(error < tolerance){
    factor = std::pow(tolerance / error, 0.25);
    acceptFlag = true;
    std::cout << "accept!" << std::endl;
    assert(factor > 1.0);
    DeltaUState = fourthOrder;
  }
  else{
    factor = std::pow(tolerance / error, 1.20);
    acceptFlag = false;
    std::cout << "reject!" << std::endl;
    
  }
  deltaTRecommended = deltaT * factor;
  Ustate = TempStorage;
}




template<int dim>
void
EllipticProblem<dim>::QSolve_NEW(dealii::Vector<heat::real> & Qdest,
			     dealii::Vector<heat::real> & mudest,
			     const dealii::Vector<heat::real> & source1,
			     const dealii::Vector<heat::real> & source2)
  
{
  // TODO:  Need to figure out a way to handle the case when there is no
  // Neumann boundary conditions.

  // Assert that the sizes are correct?
  
  dealii::deallog.push("QSolve_NEW");
  // The following is almost exactly the same as Usolve
  // TODO:  Figure out a way to combine these.
  
  dealii::PreconditionIdentity identity;
  dealii::IterativeInverse<dealii::Vector< heat::real > > MassQInverse;
  MassQInverse.initialize(MassQ, InverseMassQ);
  MassQInverse.solver.select("cg");
  //TODO:  Magic Numbers

  static dealii::SolverControl innerControl(5, 1e-10);
  MassQInverse.solver.set_control(innerControl);

  dealii::Vector<heat::real> temp(Qstate.size() );
  {
    dealii::Vector<heat::real> schur_rhs(constraintNeuQ.size() );
    try {
      MassQInverse.vmult(temp, source1);
    } catch (...) {
      std::cerr << "Failure in the first QSolve inversion." << std::endl;
      std::abort();
    }
    
    TraceQNeu.vmult(schur_rhs, temp);
    
    
    schur_rhs -= source2;

    heat::SchurComplement
      schurComplement(TraceQNeu,MassQInverse);
    heat::ApproximateSchurComplement
      approx_schurcomplement(MassQ,TraceQNeu);
    dealii::IterativeInverse<dealii::Vector<heat::real> >
      preconditionInner;
    preconditionInner.initialize(approx_schurcomplement, identity);
    preconditionInner.solver.select("cg");
    preconditionInner.solver.set_control(innerControl);
    //TODO:  Magic Numbers
    dealii::SolverControl solver_control( 40, 1e-8);        
    
    dealii::SolverCG<> cg(solver_control);

    //std::cout << "Q inner Solve" << std::endl;
    try {
      cg.solve(schurComplement, mudest, schur_rhs, identity);
    } catch (...) {
      std::cerr << "Failure in the QSolve SchurComplement" << std::endl;
      
      std::abort();
    }    
  }
  TraceQNeu.Tvmult(temp,mudest);
  temp *= -1;
  temp += source1;
  
  try {
    MassQInverse.vmult(Qdest,temp);  
  } catch (...) {
    std::cerr << "Failure in the QSolve Second Solve" << std::endl;
    std::abort();
  }
  dealii::deallog.pop();  
}
*/
template<int dim>
void
EllipticProblem<dim>::PrintState(const unsigned int timeStampNumber)
{
  // We have found that printing thing in the naive 
  // way was the biggest bottleneck in the program.
  // So we are printing things this way;
  {
    typedef dealii::DataComponentInterpretation::DataComponentInterpretation DCI_type;
    std::vector<std::string> SolutionNames;
    std::vector<DCI_type>
      DataComponentInterpretation;    

	
    SolutionNames.emplace_back("Density");
    DataComponentInterpretation
      .push_back(dealii
		 ::DataComponentInterpretation
		 ::component_is_scalar);

    for(unsigned int i = 0; i < dim; ++i){
      std::string name = "DiffusiveCurrent";
      name += dealii::Utilities::int_to_string(i,1);
      SolutionNames.push_back(name);
      DataComponentInterpretation
	.push_back(dealii
	::DataComponentInterpretation
	::component_is_scalar);
	//::component_is_part_of_vector);      
    }
    
    
    dealii::DataOut<dim> dataOutSys;
    dataOutSys.attach_dof_handler(dof_handler);
    
    //TODO:  Don't really like having a printing function change the
    //state.
    //We should just use a vector that is local to this function
    //store the state, i think.  Either that, or we should have a
    //different function set the state.

    
    state.block(0) = Ustate;
    state.block(1) = Qstate;

    
    dataOutSys.add_data_vector(state,
			       SolutionNames,
			       dealii::DataOut<dim>::type_dof_data,
			       DataComponentInterpretation);
    dataOutSys.build_patches(degree*4);
    std::ostringstream filename;
    filename 
      << "State-" 
      << dealii::Utilities::int_to_string(timeStampNumber, 4) 
      << ".vtu";
    
    std::ofstream output (filename.str().c_str() );
    dataOutSys.write_vtu(output);
  }
}

 
} // end namespace heat

template class heat::EllipticProblem<2>;
template class heat::EllipticProblem<3>;

