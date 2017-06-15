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
#include "SS_UStar_massQ_U_SSStar.hpp"
#include "SS_UStar_TotalQFromU_U_QQStar.hpp"
#include "ApproximateSchurComplement.hpp"
#include "SchurComplementPoisson.hpp"
#include "SchurComplementNEW.hpp"
#include "ApproximateSchurComplementPoisson.hpp"
#include "LDGIntegrator.hpp"
#include "LDGIntegratorKENT.hpp"
#include "LDGErrorIntegrator.hpp"
#include "InverseDiffusivity.hpp"
#include "myMatrixTools.hpp"
#include "zero_matrix.hpp"
#include "number_chain.hpp"

namespace heat{

const auto pi = boost::math::constants::pi<heat::real>();

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
  EstimatedError(1),
  UeigenStates(n_eigenvalues),
  QeigenStates(n_eigenvalues)
  
  //mapping(degree+2),
  //  faceMapping(degree+2)
{}

template<int dim>
void EllipticProblem<dim>::run()
{

  make_grid();
  make_dofs();
  
  Ustate = 0;
  Qstate = 10;

  state.block(0) = Ustate;
  state.block(1) = Qstate;

  PrintState(0);  

  auto StateOLD = state;
  auto UStateOLD = state.block(0);
  auto QStateOLD = state.block(1);  

  // check projection
  {
    dealii::Vector<double> difference(triangulation.n_active_cells() );
    dealii::VectorTools::integrate_difference(mapping,
					      dof_handlerU,
					      UStateOLD,
					      heat::DirichletBoundaryValues<dim>(),
					      difference,
					      dealii::QGauss<dim>(degree+2),
					      dealii::VectorTools::NormType::L2_norm
					     );    

    const auto error = difference.l2_norm();

    std::cout << "error = " << error << std::endl;
  }

  
    
  
  dealii::SolutionTransfer<dim,dealii::Vector<heat::real> >
    U_soltrans(dof_handlerU);

  dealii::SolutionTransfer<dim,dealii::Vector<heat::real> >
    Q_soltrans(dof_handlerQ);  
  
  dealii::GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
  							    EstimatedError.block(0),
  							    0.3,
  							    0);
  

  triangulation.prepare_coarsening_and_refinement();
  U_soltrans.prepare_for_coarsening_and_refinement(UStateOLD);
  Q_soltrans.prepare_for_coarsening_and_refinement(QStateOLD);

  triangulation.execute_coarsening_and_refinement();

  //dof_handler.clear();
  //dof_handler.initialize(triangulation, feSystem);  
  dof_handler.distribute_dofs(feSystem);
  dof_handler.initialize_local_block_info();
  
  dof_handlerU.distribute_dofs(feU);
  //dof_handlerU.initialize_local_block_info();  
  dof_handlerQ.distribute_dofs(feQ);
  //dof_handlerQ.initialize_local_block_info();
  
  Ustate.reinit(dof_handlerU.n_dofs());
  Qstate.reinit(dof_handlerQ.n_dofs());
  
  state.reinit(2);
  state.block(0).reinit(dof_handlerU.n_dofs());
  state.block(1).reinit(dof_handlerQ.n_dofs());
  state.collect_sizes();

  U_soltrans.interpolate(UStateOLD, Ustate);
  Q_soltrans.interpolate(QStateOLD, Qstate);

  state.block(0) = Ustate;
  state.block(1) = Qstate;



  
  
  {
    dealii::Vector<double> difference(triangulation.n_active_cells() );
    dealii::VectorTools::integrate_difference(mapping,
					      dof_handlerU,
					      Ustate,
					      heat::DirichletBoundaryValues<dim>(),
					      difference,
					      dealii::QGauss<dim>(degree+2),
					      dealii::VectorTools::NormType::L2_norm
					     );    
    
    const auto error = difference.l2_norm();
    std::cout << "error = " << error << std::endl;
  }
  
  
  PrintState(1);

  std::cout << "Ustate" << std::endl;
  std::cout <<  Ustate << std::endl;

  std::cout << "Qstate" << std::endl;
  std::cout <<  Qstate << std::endl;

  assert(false);
  
}
//23456789012345678901345678901234567890123456789012345678901234567890123456789

template<int dim>
void
EllipticProblem<dim>::make_grid()
{
  const bool colorize = false;
  dealii::GridGenerator::hyper_cube(triangulation, false);
  triangulation.refine_global(refinements);
}

template<int dim>
void
EllipticProblem<dim>::make_dofs()
{
  //This is pretty tricky
  // https://groups.google.com/forum/#!searchin/dealii/2-d$203-d/dealii/0cwQcfG
  //         2zQ/Pc2G3zIoH-YJ

  // There is a good reason why we use the dof_handlerU for the trace spaces
  // of Q.  It has to do with the fact that we will be dotting with the unit
  // normal vector.
  dof_handler.distribute_dofs(feSystem);
  dof_handler.initialize_local_block_info();
  dof_handlerU.distribute_dofs(feU);
  dof_handlerU.initialize_local_block_info();  
  dof_handlerQ.distribute_dofs(feQ);
  dof_handlerQ.initialize_local_block_info();
				    

  if( !boundary_ids_Dir_minus.empty() ){
    DirSurfaceToVolumeTriangulationMap =
      dealii::
      GridGenerator::
      extract_boundary_mesh(dof_handlerU,
			    dof_handlerTraceU_Dir,
			    boundary_ids_Dir_minus);

    for(auto it = DirSurfaceToVolumeTriangulationMap.begin();
	it != DirSurfaceToVolumeTriangulationMap.end();
	++it){
      DirVolumeToSurfaceTriangulationMap
	.emplace(std::make_pair(it->second, it->first));
    }
    dof_handlerTraceU_Dir.distribute_dofs(feTraceU_Dir);
  }


  if ( !boundary_ids_Neu_plus.empty() ){
    this->NeuSurfaceToVolumeTriangulationMap =
      dealii::GridGenerator::extract_boundary_mesh(dof_handlerU,
						   dof_handlerTraceU_Neu,
						   boundary_ids_Neu_plus);

    for(auto it = NeuSurfaceToVolumeTriangulationMap.begin();
	it != NeuSurfaceToVolumeTriangulationMap.end();
	++it){
      NeuVolumeToSurfaceTriangulationMap
	.emplace(std::make_pair(it->second, it->first) );
    }
    dof_handlerTraceU_Neu.distribute_dofs(feTraceU_Neu);
  }

  //TODO:  Robin stuff has never been tested.
  if( !boundary_ids_Rob_plus.empty() ){
    RobSurfaceToVolumeTriangulationMap =
      dealii::GridGenerator::extract_boundary_mesh(dof_handlerU,
						   dof_handlerTraceU_Rob,
						   boundary_ids_Rob_plus);

    for(auto it = RobSurfaceToVolumeTriangulationMap.begin();
	it != RobSurfaceToVolumeTriangulationMap.end();
	++it){
      RobVolumeToSurfaceTriangulationMap
	.emplace(std::make_pair(it->second, it->first) );
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
      if((i == j) || true){
	interiorMask[i][j] = dealii::DoFTools::always;
	fluxMask[i][j] = dealii::DoFTools::always;
      }
    }
  }



  //dealii::DoFTools::count_dofs_per_block(dof_handler,
  dealii::BlockIndices bi(dofs_per_block);
  dealii::BlockDynamicSparsityPattern csp(bi,bi);


  sparsity_pattern.reinit(2,2);

  const unsigned int
    n_u_dpc = feU.n_dofs_per_cell(),
    n_q_dpc = feQ.n_dofs_per_cell();

  sparsity_pattern.block(0,0).reinit(n_u,n_u,n_u_dpc);
  sparsity_pattern.block(0,1).reinit(n_u,n_q, 2 * (2 * dim + 1) * n_q_dpc);


  sparsity_pattern.block(1,0).reinit(n_q,n_u,2 * (2 * dim + 1) * n_u_dpc);
  sparsity_pattern.block(1,1).reinit(n_q,n_q, 4 * dim * n_q_dpc);



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

  BlockMinusRHS.reinit(2);
  BlockMinusRHS.block(0).reinit(n_u);
  BlockMinusRHS.block(1).reinit(n_q);
  BlockMinusRHS.collect_sizes();

  BlockConstraints.reinit(2);
  BlockConstraints.block(0).reinit(n_u);
  BlockConstraints.block(1).reinit(n_q);
  BlockConstraints.collect_sizes();
  

  

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
  QneuConstraint_RHS = state.block(1);
  //Q_RHS = state.block(1);
  Q_MinusRHS = state.block(1);


  lagrangeU.reinit(n_TraceU_Dir);
  lagrangeQ.reinit(n_TraceU_Neu);

  constraintDirU.reinit(n_TraceU_Dir);
  constraintNeuQ.reinit(n_TraceU_Neu);

  SystemMatricesCollection.add(0,0,"MassU");  //0 
  SystemMatricesCollection.add(0,0,"InverseMassU"); //1
  SystemMatricesCollection.add(1,1,"MassQ"); //2
  SystemMatricesCollection.add(1,1,"InverseMassQ");    //3
  SystemMatricesCollection.add(1,0,"StiffQFromU");//4 
  SystemMatricesCollection.add(0,1,"StiffUFromQ");//5
  SystemMatricesCollection.add(1,0,"NumericalFluxQFromU");//6
  SystemMatricesCollection.add(0,1,"NumericalFluxUFromQ");//7
  SystemMatricesCollection.add(0,1,"NumericalFluxUFromQBoundary");//8
  SystemMatricesCollection.add(0,0,"BoundaryMassU");//9
  SystemMatricesCollection.add(1,1,"BoundaryMassQ");//10
  SystemMatricesCollection.add(1,0,"NumericalFluxQFromUBoundary");//11
  SystemMatricesCollection.add(0,0,"Svd_U_for_Dir"); //12
  SystemMatricesCollection.add(0,0,"Svd_Sigma_for_Dir"); //13
  SystemMatricesCollection.add(1,1,"Svd_U_for_Neu"); //14
  SystemMatricesCollection.add(1,1,"Svd_Sigma_for_Neu"); //15

  SystemMatricesCollection.reinit(sparsity_pattern);

  SystemRHSsAndConstraints.add(&BlockMinusRHS, "BlockMinusRHS");
  SystemRHSsAndConstraints.add(&BlockConstraints, "BlockConstraints");

  

}


 







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
		   //::component_is_scalar);
		   ::component_is_part_of_vector);
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

    std::cout << "Ustate Inside Print State" << std::endl;
    std::cout << Ustate << std::endl;

    //std::cout << "state" << std::endl;
    

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
