#ifndef H_LDG_INTEGRATOR_KENT__
#define H_LDG_INTEGRATOR_KENT__

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
namespace LDGIntegratorKENT{

template<int dim>
class LDGIntegratorKENT : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
  LDGIntegratorKENT(dealii::Point<dim> referenceDirection_In,
		    const dealii::Quadrature<dim-1> & face_quadrature_In,
		    const dealii::UpdateFlags & update_flags,
		    const dealii::MappingQ1<dim,dim>& mapping_In);
    

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

};

// template<int dim>
// class MassQ  : public dealii::MeshWorker::LocalIntegrator<dim>
// {
// public:
// MassQ();

// virtual
// void
// cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
//      dealii::MeshWorker::IntegrationInfo<dim> &info) const;

// };


// template<int dim>
// class MassU  : public dealii::MeshWorker::LocalIntegrator<dim>
// {
// public:
// MassU();

// virtual
// void
// cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
//      dealii::MeshWorker::IntegrationInfo<dim> &info) const;

// };


// template<int dim>
// class StiffUFromQ  : public dealii::MeshWorker::LocalIntegrator<dim>
// {
// public:
//   StiffUFromQ();

// virtual
// void
// cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
//      dealii::MeshWorker::IntegrationInfo<dim> &info) const;
// };


// template<int dim>
// class JumpQ  : public dealii::MeshWorker::LocalIntegrator<dim>
// {
// public:
//   JumpQ();
  

// virtual
// void
// face(dealii::MeshWorker::DoFInfo<dim> & dinfo0,
//      dealii::MeshWorker::DoFInfo<dim> & dinfo1,
//      dealii::MeshWorker::IntegrationInfo<dim> &info0,
//      dealii::MeshWorker::IntegrationInfo<dim> &info1) const;
// };


} // End namespace LDGIntegratorKENT
} // End namepsace heat

#endif
