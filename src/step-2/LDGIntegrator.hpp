#ifndef H_LDG_INTEGRATOR__
#define H_LDG_INTEGRATOR__

#include "globals.hpp"


#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/output.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>
#include <deal.II/integrators/l2.h>
#include <deal.II/integrators/divergence.h>
#include <deal.II/lac/lapack_full_matrix.h>


namespace heat{
namespace LDGIntegrator{

template<int dim>
class LDGIntegrator : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
  LDGIntegrator(dealii::Point<dim> referenceDirection_In);  

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
  dealii::Point<dim> referenceDirection;

};

template<int dim>
class MassQ  : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
MassQ();

virtual
void
cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
     dealii::MeshWorker::IntegrationInfo<dim> &info) const;

};


template<int dim>
class MassU  : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
MassU();

virtual
void
cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
     dealii::MeshWorker::IntegrationInfo<dim> &info) const;

};


template<int dim>
class StiffUFromQ  : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
  StiffUFromQ();

virtual
void
cell(dealii::MeshWorker::DoFInfo<dim> & dinfo,
     dealii::MeshWorker::IntegrationInfo<dim> &info) const;
};


template<int dim>
class JumpQ  : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
  JumpQ();

virtual
void
face(dealii::MeshWorker::DoFInfo<dim> & dinfo0,
     dealii::MeshWorker::DoFInfo<dim> & dinfo1,
     dealii::MeshWorker::IntegrationInfo<dim> &info0,
     dealii::MeshWorker::IntegrationInfo<dim> &info1) const;
};


} // End namespace LDGIntegrator
} // End namepsace heat

#endif
