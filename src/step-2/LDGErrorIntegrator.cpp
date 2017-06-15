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

  heat::real normSquaredOfCellResidual = 0;

  const heat::SourceBodyValues<dim> source;

  std::vector<heat::real> SourceTerms(infoSelf.fe_values(0).n_quadrature_points);

  source.value_list(infoSelf.fe_values(0).get_quadrature_points(), SourceTerms);
  

  
  auto U_Vector_self = infoSelf.gradients;
  
  heat::LocalIntegrators::LDG::CellResidualQ(infoSelf.fe_values(0),
  					     U_Vector_self[0],
					     SourceTerms,
					     normSquaredOfCellResidual);
  auto termCELL
    = normSquaredOfCellResidual
    * dinfo_Self.cell->diameter()
    * dinfo_Self.cell->diameter();  

  dinfo_Self.value(0) =    normSquaredOfCellResidual
    * dinfo_Self.cell->diameter()
    * dinfo_Self.cell->diameter();  
  
  PRINT(termCELL)
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

  auto termFACE = jumpSquared * dinfoSELF.cell->diameter();

  PRINT(termFACE);
}
  
template<int dim>
void
LDGErrorIntegrator<dim>::boundary(dealii::MeshWorker::DoFInfo<dim> & dinfo,
				  dealii::MeshWorker::IntegrationInfo<dim> & info) const
{}

} // End Namespace LDGErrorIntegrator
} // End Namespace heat

template class heat::LDGErrorIntegrator::LDGErrorIntegrator<2>;
template class heat::LDGErrorIntegrator::LDGErrorIntegrator<3>;
