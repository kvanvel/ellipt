#ifndef H_ELLIPTIC_PROBLEM__
#define H_ELLIPTIC_PROBLEM__

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
#include <limits>


#include <boost/math/constants/constants.hpp>


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/schur_matrix.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/transpose_matrix.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h> 
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>


#include "globals.hpp"

namespace heat {

template <int dim>
class EllipticProblem
{
public:
  EllipticProblem(const unsigned int degree, const unsigned int refinements);
  void run();

  heat::real L2errorU, L2errorQ;

private:
  void 
  assemble_system();

  void 
  assemble_system_UGLY();
  

  void
  print_system();

  void
  output_results_eigensystem();
  
  void 
  make_grid();
  
  void 
  make_dofs();
  
  void
  report_dofs();
  
  void
  PrintState(const unsigned int timeStampNumber);
  
  void
  makeConstraintDirU(const heat::real time);
  
  void
  makeConstraintNeuQ(const heat::real time);
  
  void
  makeU_RHS_From_U();

  void
  assembleFluxLocals
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
   dealii::FullMatrix<heat::real> Q_U_Matrix_2);  
  
  void 
  assembleMassULocal(const dealii::FEValues<dim> &fe_v,
		     dealii::FullMatrix<heat::real> &ui_vi_matrix
		     );

  void  
  assembleBoundaryMassULocal(const dealii::FEFaceValuesBase<dim> & fe_v,
			     dealii::FullMatrix<heat::real> & ui_vi_matrix);

  
  void 
  assembleMassQLocal(const dealii::FEValues<dim> &fe_v,
		     dealii::FullMatrix<heat::real> &ui_vi_matrix);

  void
  assembleBoundaryMassQLocal(const dealii::FEFaceValuesBase<dim> & fe_v,
			     dealii::FullMatrix<heat::real> & pi_qj_matrix);
  
  void 
  assembleStiffUFromQLocal(const dealii::FEValues<dim> & fe_v_U,
			   const dealii::FEValues<dim> & fe_v_Q,
			   dealii::FullMatrix<heat::real> & ui_qj_matrix);

  void
  assembleStiffQFromULocal(const dealii::FEValues<dim> & fe_v_Q,
			   const dealii::FEValues<dim> & fe_v_U,
			   dealii::FullMatrix<heat::real> & qi_uj_matrix);
  
  void
  assembleFluxUFromULocal
  (const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> & fe_v_self_Elec,
   const dealii::FEFaceValuesBase<dim> & fe_v_neig_U,
   const dealii::FEFaceValuesBase<dim> & fe_v_neig_Elec,
   dealii::Vector<heat::real> & Uin_Uin_vector);

  void
  assembleStiffUFromULocal
  (const dealii::FEValuesBase<dim> & fe_v_self_U,
   const dealii::FEValuesBase<dim> & fe_v_self_Elec,
   dealii::Vector<heat::real> & Uin_vector);
  
  void
  assembleFluxNegUFromULocal
  (const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> & fe_v_self_Elec,
   const dealii::FEFaceValuesBase<dim> & fe_v_neig_U,
   const dealii::FEFaceValuesBase<dim> & fe_v_neig_Elec,
   dealii::Vector<heat::real> & Uin_Uin_vector,
   dealii::Vector<heat::real> & Uin_Uout_vector);
    
  void
  assembleFluxCenterUFromQLocal
  (const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
   const dealii::FEFaceValuesBase<dim> &fe_v_neig_Q,
   dealii::FullMatrix<heat::real> & Uin_Qin_matrix,
   dealii::FullMatrix<heat::real> & Uin_Qout_matrix);

  void
  assembleFluxPosUFromQLocal
  (const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
   const dealii::FEFaceValuesBase<dim> &fe_v_neig_Q,
   dealii::FullMatrix<heat::real> & Uin_Qin_matrix,
   dealii::FullMatrix<heat::real> & Uin_Qout_matrix);
  
  void
  assembleFluxNegUFromQLocal
  (const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
   const dealii::FEFaceValuesBase<dim> &fe_v_neig_Q,
   dealii::FullMatrix<heat::real> & Uin_Qin_matrix,
   dealii::FullMatrix<heat::real> & Uin_Qout_matrix);
  
  void
  assembleFluxCenterQFromULocal
  (const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
   const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> &fe_v_neig_U,
   dealii::FullMatrix<heat::real> & Qin_Uin_matrix,
   dealii::FullMatrix<heat::real> & Qin_Uout_matrix);

  void
  assembleFluxPosQFromULocal
  (const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
   const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> &fe_v_neig_U,
   dealii::FullMatrix<heat::real> & Qin_Uin_matrix,
   dealii::FullMatrix<heat::real> & Qin_Uout_matrix);

  void
  assembleFluxNegQFromULocal
  (const dealii::FEFaceValuesBase<dim> &fe_v_self_Q,
   const dealii::FEFaceValuesBase<dim> &fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> &fe_v_neig_U,
   dealii::FullMatrix<heat::real> & Qin_Uin_matrix,
   dealii::FullMatrix<heat::real> & Qin_Uout_matrix);  

  void
  assembleFluxUFromQBoundaryLocal
  (const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
   const dealii::FEFaceValuesBase<dim> & fe_v_self_Q,
   dealii::FullMatrix<heat::real> & Uin_Qin_matrix);
  
  void
  assembleFluxQFromUBoundaryLocal
  (const dealii::FEFaceValuesBase<dim> & fe_v_self_Q,
   const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
   dealii::FullMatrix<heat::real> & Qin_Uin_matrix);

  void
  assembleMassElecLocal
  (const dealii::FEValuesBase<dim> & fe_v_self_Elec,
   dealii::FullMatrix<heat::real> & Elec_Elec_matrix);
  
  void
  assembleStiffElecFromPotLocal
  (const dealii::FEValuesBase<dim> & fe_v_Elec,
   const dealii::FEValuesBase<dim> & fe_v_Pot,
   dealii::FullMatrix<heat::real> & E_P_matrix);
  
  void
  assembleStiffPotFromElecLocal
  (const dealii::FEValuesBase<dim> & fe_v_Pot,
   const dealii::FEValuesBase<dim> & fe_v_Elec,
   dealii::FullMatrix<heat::real> & P_E_Matrix);

  
  void
  assemblePotFromULocal
  (const dealii::FEValuesBase<dim> & fe_v_Pot,
   const dealii::FEValuesBase<dim> & fe_v_U,
   dealii::FullMatrix<heat::real> & P_U_Matrix);
    
  void
  assembleTraceULocalDir
  (const dealii::FEValuesBase<dim-1,dim> & fe_v_self_TrU,
   const dealii::FEFaceValuesBase<dim> & fe_v_self_U,
   dealii::FullMatrix<heat::real> & TrU_U_matrix);

  void
  assembleTraceUElecDir
  (const dealii::FEValuesBase<dim-1,dim> & fe_v_self_TrU,
   const dealii::FEFaceValuesBase<dim> & fe_v_self_E,
   dealii::FullMatrix<heat::real> & TrU_Dir_Elec_Matrix);
  
  void
  assembleTraceQLocalNeu
  (const dealii::FEValuesBase<dim-1, dim> & fe_v_TraceU_Neu,
   const dealii::FEFaceValuesBase<dim> & fe_v_face_Q,
   dealii::FullMatrix<heat::real> & TrU_Neu_Q_Matrix);

  
  void
  assembleQ_RHS_FromDirPlusLocal
  (const dealii::FEFaceValuesBase<dim> & fe_v_face_Q,
   dealii::Vector<heat::real> & Q_Vector);

  
  void
  assembleConstraintUDirMinusLocal
  (const dealii::FEValuesBase<dim-1,dim> & fe_v_self_TrU,
   dealii::Vector<heat::real> & TrU_U_Vector);

  void
  assembleConstraintUDirMinusLocalNOTRACE
  (const dealii::FEFaceValuesBase<dim> & fe_v_face_U,
   dealii::Vector<heat::real> & U_Vector);


  void
  assembleConstraintQNeuPlusLocal
  (const dealii::FEValuesBase<dim-1, dim> & fe_v_self_TrU,
   const dealii::FEFaceValuesBase<dim> & fe_v_face_Q,
   dealii::Vector<heat::real> & TrU_Q_Vector);

  void
  assembleConstraintQNeuPlusLocalNOTRACE
  (const dealii::FEFaceValuesBase<dim> & fe_v_face_Q,
   dealii::Vector<heat::real> & Q_Vector);

  void
  assembleQ_RHS_FromNeuMinusLocal
  (
   const dealii::FEFaceValuesBase<dim> & fe_v_face_U,
   dealii::Vector<heat::real> & U_Vector);

  void
  assembleU_RHS_FromSourceLocal
  ( const dealii::FEValuesBase<dim> & fe_v_u,
    dealii::Vector<heat::real> & U_Vector);


  // void
  // assembleTraceElecLocalNeu
  // (const dealii::FEValuesBase<dim-1,dim> & fe_v_Trace_U_Neu,
  //  const dealii::FEFaceValuesBase<dim> & fe_v_face_Elec,
  //  dealii::FullMatrix<heat::real> & TrU_Neu_Elec_Matrix);

  // void
  // makeU_RHS_FromInitialCondition();
  
  // void
  // makeElec_RHS_From_DirichletBoundary(const heat::real time);
  
  // void
  // makeConstraint_RHS_From_NeumanBoundary();

  void
  makeU_RHS_FromQ();

  void
  makeQ_RHS_FromU();

  void
  ElliptSolve();

  void
  ElliptSolveNEW();

  void
  EigenSolve();

  void
  EigenSolveNEW();
  
  void
  USolve();

   
  void
  QSolve_NEW(dealii::Vector<heat::real> & Qdest,
	     dealii::Vector<heat::real> & mudest,
	     const dealii::Vector<heat::real> & source1,
	     const dealii::Vector<heat::real> & source2);

  void
  USolve_NEW(dealii::Vector<heat::real> & Udest,
	     dealii::Vector<heat::real> & lambdadest,
	     const dealii::Vector<heat::real> & source1,
	     const dealii::Vector<heat::real> & source2);

 
  // void
  // PoissonSolve();
  
  // void
  // ComputeUnext(heat::real timeNow,
  // 	       heat::real timeNext,
  // 	       heat::real deltaT);
  
  // void
  // ForwardEuler(const heat::real timeNow,
  // 	       const heat::real deltaT,
  // 	       bool & acceptFlag,
  // 	       heat::real & deltaTRecommended);
  
  // void
  // AdaptiveExplicitRungeKutta54(const heat::real timeComputed,
  // 			       const heat::real deltaT,
  // 			       bool & acceptFlag,
  // 			       heat::real & deltaTRecommended);
  
  // void
  // ComputeDensityDot(heat::real time);

  // void
  // simulate();    
  

  void  
  InsertGlobalFromLocal
  (const typename dealii::DoFHandler<dim>::active_cell_iterator& test_iterator,
   const typename dealii::DoFHandler<dim>::active_cell_iterator& tria_iterator,
   const dealii::FullMatrix<heat::real>& local_matrix,
   dealii::SparseMatrix<heat::real>& global_matrix);  
  
  void
  InsertGlobalFromLocal
  (
   const typename
   dealii::DoFHandler<dim-1, dim >::active_cell_iterator & test_iterator,

   const typename
   dealii::DoFHandler<dim,dim>::active_cell_iterator & tria_iterator,
   const dealii::FullMatrix<heat::real> & local_matrix,
   dealii::SparseMatrixEZ<heat::real> & global_matrix);


  void
  InsertGlobalFromLocalVector
  (const typename dealii::DoFHandler<dim>::active_cell_iterator & test_iterator,
   const dealii::Vector<heat::real> & local_vector,
   dealii::Vector<heat::real> & global_vector);

  void
  InsertGlobalFromLocalVectorNEW
  (const typename dealii::DoFHandler<dim-1,dim>::active_cell_iterator & test_iterator,
   const dealii::Vector<heat::real> & local_vector,
   dealii::Vector<heat::real> & global_vector);


  std::map<typename dealii::DoFHandler<dim-1, dim>::cell_iterator,
	   typename dealii::DoFHandler<dim  , dim>::face_iterator>
    DirSurfaceToVolumeTriangulationMap,
    NeuSurfaceToVolumeTriangulationMap,
    RobSurfaceToVolumeTriangulationMap;

  std::map<typename dealii::DoFHandler<dim  ,dim>::face_iterator,
	   typename dealii::DoFHandler<dim-1,dim>::cell_iterator>
    DirVolumeToSurfaceTriangulationMap,
    NeuVolumeToSurfaceTriangulationMap,
    RobVolumeToSurfaceTriangulationMap;
  
  dealii::Point<dim> referenceDirection;
  
  
  
  const unsigned int degree;
  const unsigned int refinements;
  const unsigned int n_eigenvalues;
  unsigned int rankBoundaryMassU;
  unsigned int nullityBoundaryMassU;
  unsigned int rankBoundaryMassQ;
  unsigned int nullityBoundaryMassQ;
  
  dealii::Triangulation<dim> triangulation;
  dealii::Triangulation<dim-1,dim> triangulation_Dir;
  dealii::Triangulation<dim-1,dim> triangulation_Neu;
  dealii::Triangulation<dim-1,dim> triangulation_Rob;
  
  dealii::FESystem<dim> feSystem;
  dealii::FE_DGQ<dim> feU;
  dealii::FESystem<dim> feQ;  
  dealii::FE_DGQ<dim-1,dim> feTraceU_Dir;
  dealii::FE_DGQ<dim-1,dim> feTraceU_Neu;
  dealii::FE_DGQ<dim-1,dim> feTraceU_Rob;
  
  dealii::ConstraintMatrix constraintU;
  dealii::ConstraintMatrix constraintElec;
  dealii::ConstraintMatrix constraintSys;
  
  dealii::DoFHandler<dim> dof_handler;
  dealii::DoFHandler<dim> dof_handlerU;
  dealii::DoFHandler<dim> dof_handlerQ;
  
  dealii::DoFHandler<dim-1,dim> dof_handlerTraceU_Dir;
  dealii::DoFHandler<dim-1,dim> dof_handlerTraceU_Neu;
  dealii::DoFHandler<dim-1,dim> dof_handlerTraceU_Rob;

  //dealii::DoFHandler<dim> dof_handlerPot;
  //dealii::DoFHandler<dim> dof_handlerElec;
  
  std::set<dealii::types::boundary_id> boundary_ids_Dir_plus, boundary_ids_Dir_minus;
  std::set<dealii::types::boundary_id> boundary_ids_Neu_plus, boundary_ids_Neu_minus;
  std::set<dealii::types::boundary_id> boundary_ids_Rob_plus, boundary_ids_Rob_minus;
  
  dealii::BlockSparsityPattern sparsity_pattern;
  dealii::BlockSparseMatrix<heat::real> system_matrix;
  dealii::BlockSparseMatrix<heat::real> system_matrix2;
  dealii::BlockSparseMatrix<heat::real> system_matrix3;
  dealii::BlockSparseMatrix<heat::real> system_matrix4;
  

  dealii::SparsityPattern 
    MassUSparsityPattern,
    InverseMassUSparsityPattern,
    TotalUFromQSparsityPattern,
    StiffUFromQSparsityPattern,
    FluxPosUFromQSparsityPattern,
    FluxNegUFromQSparsityPattern,
    FluxCenterUFromQSparsityPattern,
    FluxBoundaryUFromQSparsityPattern,
    MassQSparsityPattern,
    InverseMassQSparsityPattern,
    TotalQFromUSparsityPattern,
    StiffQFromUSparsityPattern,
    FluxPosQFromUSparsityPattern,
    FluxNegQFromUSparsityPattern,
    FluxCenterQFromQSparsityPattern,
    FluxBoundaryQFromUSparsityPattern,
    TraceUDirSparsityPattern,
    TraceQNeuSparsityPattern,
    TraceQRobSparsityPattern,
  // MassElecSparsityPattern,
  // StiffElecFromPotSparsityPattern,
  // StiffPotFromElecSparsityPattern,
  // PotFromUSparsityPattern,
   TRACE_U_DIR_NEW_SparsityPattern;

  
  
  // std::vector<std::vector<unsigned int> > MassUCols;
  // std::vector<std::vector<std::pair<unsigned int, heat::real> > >  
  // MassUColsVals;
  
  dealii::SparseMatrix<heat::real>
  InverseMassU,
    MassU,
    BoundaryMassU,    
    Svd_U_for_Dir,
    Svd_Sigma_for_Dir,
    StiffUFromQ,
    FluxPosUFromQ,
    FluxNegUFromQ,
    FluxCenterUFromQ,
    FluxBoundaryUFromQ,
    TotalUFromQ,
    InverseMassQ,
    InverseMassQCholesky,
    MassQ,
    BoundaryMassQ,
    Svd_U_for_Neu,
    Svd_Sigma_for_Neu,
    StiffQFromU,
    FluxPosQFromU,
    FluxNegQFromU,
    FluxCenterQFromU,
    FluxBoundaryQFromU,
    TotalQFromU,
  //MassElec,
  //    StiffElecFromPot,
  // StiffPotFromElec,
  //PotFromU,
  TRACE_U_DIR_NEW
    ;
  
  dealii::SparseDirectUMFPACK MassUInverseDirect;
  dealii::SparseDirectUMFPACK MassQInverseDirect;
  //dealii::SparseDirectUMFPACK MassElecInverseDirect;
  
  dealii::SparseMatrixEZ<heat::real>
    TraceUDir,
    TraceQNeu,
    PP,
    QQ,
    RR,
    SS;
  
  
  dealii::BlockVector<heat::real> state;
  dealii::BlockVector<heat::real> BlockMinusRHS;
  dealii::BlockVector<heat::real> BlockConstraints;
  
  dealii::Vector<heat::real> Ustate;
  dealii::Vector<heat::real> UperpNEW;
  dealii::Vector<heat::real> UparNEW;
  //dealii::Vector<heat::real> DeltaUState;
  //dealii::Vector<heat::real> UStateDot;
  dealii::Vector<heat::real> lagrangeU;
  dealii::Vector<heat::real> lambdaU;
  dealii::Vector<heat::real> constraintDirU;
  dealii::Vector<heat::real> Qstate;
  dealii::Vector<heat::real> QperpNEW;
  dealii::Vector<heat::real> QparNEW;
  dealii::Vector<heat::real> lagrangeQ;
  dealii::Vector<heat::real> constraintNeuQ; //And Rob?
  //dealii::Vector<heat::real> U_RHS;
  dealii::Vector<heat::real> U_MinusRHS;
  dealii::Vector<heat::real> U_RHS_From_U;
  dealii::Vector<heat::real> UdirConstraint_RHS;
  dealii::Vector<heat::real> QneuConstraint_RHS;
  //dealii::Vector<heat::real> Q_RHS;
  dealii::Vector<heat::real> Q_MinusRHS;
  //dealii::Vector<heat::real> PotState;
  //dealii::Vector<heat::real> ElecState;
  //dealii::Vector<heat::real> Elec_RHS;
  //dealii::Vector<heat::real> PoissonNeuBC_RHS;
  //dealii::Vector<heat::real> lagrangeElec;
  std::vector< dealii::Vector<heat::real> > UeigenStates;
  std::vector< dealii::Vector<heat::real> > QeigenStates;
  
  dealii::MappingQ<dim> mapping;
  dealii::MappingQ<dim-1,dim> faceMapping;

  std::vector<dealii::types::global_dof_index >
    uDofsWithSupportOnBoundary,
    qDotsWithOrthongalSupportOnBoundary;

  
};  

} //end namespace heat

#endif
