#ifndef H_LDG_INTEGRATOR_DRIFT__
#define H_LDG_INTEGRATOR_DRIFT__

#include "globals.hpp"
#include "myMatrixTools.hpp"


#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/output.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>
#include <deal.II/integrators/l2.h>
#include <deal.II/integrators/divergence.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/fe/mapping_q.h>


namespace heat{
namespace LDGIntegratorDrift{

template<int dim>
class LDGIntegratorDrift : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
  LDGIntegratorDrift(
		     const dealii::Point<dim> referenceDirection,
		     const dealii::Quadrature<dim-1> & face_quadrature_In,
		     const dealii::UpdateFlags & update_flags,
		     const dealii::MappingQ1<dim,dim>& mapping_In,
		     const dealii::Vector<heat::real> & ElecState_In);
		     
    

  void
  cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
       dealii::MeshWorker::IntegrationInfo<dim> &info) const; 

  void
  face(dealii::MeshWorker::DoFInfo<dim> & dinfoSELF,
       dealii::MeshWorker::DoFInfo<dim> & dinfoNEIG,
       dealii::MeshWorker::IntegrationInfo<dim> &infoSELF,
       dealii::MeshWorker::IntegrationInfo<dim> &infoNEIG ) const;

  void
  boundary(dealii::MeshWorker::DoFInfo<dim> & dinfo,
	   dealii::MeshWorker::IntegrationInfo<dim> & info) const;
  
  
private:
  const dealii::Point<dim> referenceDirection;
  const dealii::SmartPointer<const dealii::Quadrature<dim-1> > face_quadrature;
  const dealii::UpdateFlags face_update_flags;
  const dealii::MappingQ1<dim,dim> mapping;
  const dealii::SmartPointer<const dealii::Vector<heat::real> > ElecState;

};

} // End namespace LDGIntegratorDrift
} // End namepsace heat

#endif
