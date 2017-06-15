#include "LDGErrorIntegrator.hpp"
#include "LDG.hpp"

namespace heat{
namespace LDGErrorIntegrator{

template<int dim>
LDGErrorIntegrator<dim>::LDGErrorIntegrator(dealii::Point<dim> referenceDirection_In)
  :
  dealii::MeshWorker::LocalIntegrator<dim>(true, true, true)  
{}

template<int dim>
void
LDGErrorIntegrator<dim>::cell(dealii::MeshWorker::DoFInfo<dim> & dinfo_Self,
			      dealii::MeshWorker::IntegrationInfo<dim> &infoSelf ) const{

  heat::real normSquaredOfCellResidualQ = 0;
  heat::real normSquaredOfCellResidualU = 0;
  //heat::real normSquaredOfCellResidualQ = 0;
  
  

  //const heat::SourceBodyValues<dim> source;

  //std::vector<heat::real> SourceTerms(infoSelf.fe_values(0).n_quadrature_points);

  //source.value_list(infoSelf.fe_values(0).get_quadrature_points(), SourceTerms);
  
  //auto U_Vector_self = infoSelf.gradients;
  
  // heat::LocalIntegrators::LDG::CellResidualQ(infoSelf.fe_values(0),
  // 					     U_Vector_self[0],
  // 					     SourceTerms,
  // 					     normSquaredOfCellResidual);  

  // std::vector<dealii::types::global_dof_index>
  //   local_dof_indices_Q(infoSelf.fe_values(1).dofs_per_cell),
  //   local_dof_indices_U(infoSelf.fe_values(0).dofs_per_cell);

  dealii::Vector<heat::real>
    local_dofs_Q(infoSelf.fe_values(1).dofs_per_cell),
    local_dofs_U(infoSelf.fe_values(0).dofs_per_cell);

  const auto & Indices_By_Block = dinfo_Self.indices_by_block;

  const auto & Q_Indices = Indices_By_Block[1];
  const auto & U_Indices = Indices_By_Block[0];

  

  std::shared_ptr<dealii::MeshWorker::VectorDataBase<dim,dim> >
    GlobalDataPtr = infoSelf.global_data;

  const auto & globalDataCollection = GlobalDataPtr->data;

  const auto & SolutionVectorPtr
    = globalDataCollection.template read_ptr<dealii::BlockVector<heat::real> > ("state");


 SolutionVectorPtr->extract_subvector_to(Q_Indices.begin(),
					 Q_Indices.end(),
					 local_dofs_Q.begin() );
  

 SolutionVectorPtr->extract_subvector_to(U_Indices.begin(),
					 U_Indices.end(),
					 local_dofs_U.begin() );
 
 heat
   ::LocalIntegrators
   ::LDG
   ::CellResidualU(infoSelf.fe_values(1),
		   local_dofs_Q,
		   infoSelf.fe_values(0), 
		   local_dofs_U,
		   normSquaredOfCellResidualU);

  
 heat
   ::LocalIntegrators
   ::LDG
   ::CellResidualQ(infoSelf.fe_values(1),
		   local_dofs_Q,						
		   normSquaredOfCellResidualQ);
 
 // auto termCELL
 //   = normSquaredOfCellResidualQ
 //   * dinfo_Self.cell->diameter()
 //   * dinfo_Self.cell->diameter();  

 
 
  dinfo_Self.value(0) =  (normSquaredOfCellResidualQ + 0 * normSquaredOfCellResidualU)
    * dinfo_Self.cell->diameter()
    * dinfo_Self.cell->diameter();
}


template <int dim>
void
LDGErrorIntegrator<dim>::face(dealii::MeshWorker::DoFInfo<dim> & dinfoSELF,
			      dealii::MeshWorker::DoFInfo<dim> & dinfoNEIG,
			      dealii::MeshWorker::IntegrationInfo<dim> &infoSELF,
			      dealii::MeshWorker::IntegrationInfo<dim> &infoNEIG ) const {

  heat::real jumpSquared = 0.0;

  auto Q_Vector_self = infoSELF.values;
  auto Q_Vector_neig = infoNEIG.values;

  heat::LocalIntegrators::LDG::QJump(infoSELF.fe_values(1),
  				     infoNEIG.fe_values(1),
  				     Q_Vector_self[0],
				     Q_Vector_neig[0],  				     
				     jumpSquared);
  
  dinfoSELF.value(0) = jumpSquared * dinfoSELF.cell->diameter();
  dinfoNEIG.value(0) = jumpSquared * dinfoSELF.cell->diameter();
    // sqrt(
    // 	 dinfoSELF.cell->diameter() * dinfoNEIG.cell->diameter()
    // 	 );
  
}
  
template<int dim>
void
LDGErrorIntegrator<dim>::boundary
(dealii::MeshWorker::DoFInfo<dim> & dinfo,
 dealii::MeshWorker::IntegrationInfo<dim> & info) const
{
  
}

} // End Namespace LDGErrorIntegrator
} // End Namespace heat

template class heat::LDGErrorIntegrator::LDGErrorIntegrator<2>;
template class heat::LDGErrorIntegrator::LDGErrorIntegrator<3>;
