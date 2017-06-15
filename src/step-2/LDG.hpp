#ifndef __H_LDG__
#define __H_LDG__


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/meshworker/dof_info.h>


#include "SourceBodyValues.hpp"
#include "DirichletBoundaryValues.hpp"
#include "NeumannBoundaryValues.hpp"
#include <typeinfo>


//#include <deal.II/

namespace heat{

namespace LocalIntegrators{
    
namespace LDG{

template <int dim>
void StiffQFromU(
		 dealii::FullMatrix<double> & Q_in_U_in,
		 const dealii::FEValuesBase<dim> & fe_Q_Self,
		 const dealii::FEValuesBase<dim> & fe_U_Self,
		 const double factor = 1.0)
{
  const unsigned int n_Q_dofs = fe_Q_Self.dofs_per_cell;
  const unsigned int n_U_dofs = fe_U_Self.dofs_per_cell;

  AssertDimension(Q_in_U_in.m(), n_Q_dofs);
  AssertDimension(Q_in_U_in.n(), n_U_dofs);

  AssertDimension(fe_U_Self.get_fe().n_components(), 1);
  AssertDimension(fe_Q_Self.get_fe().n_components(), dim);

  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);

  const auto & JxW = fe_Q_Self.get_JxW_values();

  for(unsigned int point = 0;
      point < fe_Q_Self.n_quadrature_points;
      ++point){
    for(unsigned int i = 0; i < fe_Q_Self.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_U_Self.dofs_per_cell; ++j){
	Q_in_U_in(i,j) += factor * JxW[point] * fe_Q_Self[flux].divergence(i,point)
	  * fe_U_Self[density].value(j,point);
      }
    }
  }
}
  
template <int dim>
void StiffUFromQ(dealii::FullMatrix<double> & U_in_Q_in,
		 const dealii::FEValuesBase<dim> & fe_U_Self,
		 const dealii::FEValuesBase<dim> & fe_Q_Self,
		 const double factor = 1.0)
{
  const unsigned int n_Q_dofs = fe_Q_Self.dofs_per_cell;
  const unsigned int n_U_dofs = fe_U_Self.dofs_per_cell;

  AssertDimension(U_in_Q_in.m(), n_U_dofs);
  AssertDimension(U_in_Q_in.n(), n_Q_dofs);

  AssertDimension(fe_U_Self.get_fe().n_components(), 1);
  AssertDimension(fe_Q_Self.get_fe().n_components(), dim);

  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);

  const auto & JxW = fe_Q_Self.get_JxW_values();
  for(unsigned int point = 0;
      point < fe_U_Self.n_quadrature_points;
      ++point){
    for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_Q_Self.dofs_per_cell; ++j){
	U_in_Q_in(i,j)
	  += factor
	  * JxW[point]
	  * fe_U_Self[density].gradient(i,point)
	  * fe_Q_Self[flux].value(j,point);
      }
    }
  }
}


  



template <int dim>
void numericalFluxQFromU(
			 const dealii::Point<dim> referenceDirection,
			 dealii::FullMatrix<double> & Q_in_U_in,
			 dealii::FullMatrix<double> & Q_in_U_out,
			 dealii::FullMatrix<double> & Q_out_U_in,
			 dealii::FullMatrix<double> & Q_out_U_out,
			 const dealii::FEValuesBase<dim> & fe_Q_Self,
			 const dealii::FEValuesBase<dim> & fe_Q_Neig,
			 const dealii::FEValuesBase<dim> & fe_U_Self,
			 const dealii::FEValuesBase<dim> & fe_U_Neig,
			 const double factor = 1.0){

  const unsigned int n_Q_dofs = fe_Q_Self.dofs_per_cell;
  const unsigned int n_U_dofs = fe_U_Self.dofs_per_cell;

  AssertDimension(n_Q_dofs, fe_Q_Neig.dofs_per_cell);
  AssertDimension(n_U_dofs, fe_U_Neig.dofs_per_cell );
  
  AssertDimension(Q_in_U_in.m(), n_Q_dofs);
  AssertDimension(Q_in_U_in.n(), n_U_dofs);
  AssertDimension(Q_in_U_out.m(), n_Q_dofs);
  AssertDimension(Q_in_U_out.n(), n_U_dofs);
  AssertDimension(Q_out_U_in.m(), n_Q_dofs);
  AssertDimension(Q_out_U_in.n(), n_U_dofs);
  AssertDimension(Q_out_U_out.m(), n_Q_dofs);
  AssertDimension(Q_out_U_out.n(), n_U_dofs);

  AssertDimension(fe_U_Self.get_fe().n_components(), 1);
  AssertDimension(fe_U_Neig.get_fe().n_components(), 1);
  AssertDimension(fe_Q_Self.get_fe().n_components(), dim);
  AssertDimension(fe_Q_Neig.get_fe().n_components(), dim);

  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);

  const auto & JxW = fe_Q_Self.get_JxW_values();
  const auto & normals = fe_Q_Self.get_normal_vectors();

  

  for(unsigned int point = 0;
      point < fe_Q_Self.n_quadrature_points;
      ++point){
    const auto & normalVector = fe_Q_Self.normal_vector(point);
    if(referenceDirection * normalVector > 0){
      for(unsigned int i = 0; i < fe_Q_Self.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_U_Neig.dofs_per_cell; ++j){
	  Q_in_U_out(i,j)
	    += factor
	    * JxW[point]
	    * fe_Q_Self[flux].value(i,point)
	    * fe_U_Neig[density].value(j,point)	    
	    * normals[point];
	}
      }
    }
    else {
      for(unsigned int i = 0; i < fe_Q_Self.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_U_Self.dofs_per_cell; ++j){
	  Q_in_U_in(i,j)
	    += factor
	    * JxW[point]
	    * fe_Q_Self[flux].value(i,point)
	    * fe_U_Self[density].value(j,point)
	    * normals[point];
	}
      }
    }
  }
}

template <int dim>
void numericalFluxQFromUBoundary(
			 const dealii::Point<dim> referenceDirection,
			 dealii::FullMatrix<double> & Q_in_U_in,
			 const dealii::FEValuesBase<dim> & fe_Q_Self,			 
			 const dealii::FEValuesBase<dim> & fe_U_Self,			 
			 const double factor = 1.0){

  const unsigned int n_Q_dofs = fe_Q_Self.dofs_per_cell;
  const unsigned int n_U_dofs = fe_U_Self.dofs_per_cell;
    
  AssertDimension(Q_in_U_in.m(), n_Q_dofs);
  AssertDimension(Q_in_U_in.n(), n_U_dofs);
  
  AssertDimension(fe_U_Self.get_fe().n_components(), 1);
  AssertDimension(fe_Q_Self.get_fe().n_components(), dim);
  

  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);
  
  const auto & JxW = fe_Q_Self.get_JxW_values();
  const auto & normals = fe_Q_Self.get_normal_vectors();

  for(unsigned int point = 0;
      point < fe_Q_Self.n_quadrature_points;
      ++point){
    assert(referenceDirection * normals[point] < 0);
    for(unsigned int i = 0; i < fe_Q_Self.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_U_Self.dofs_per_cell; ++j){
	Q_in_U_in(i,j)	    
	  += factor
	  * JxW[point]
	  * fe_Q_Self[flux].value(i,point)
	  * fe_U_Self[density].value(j,point)
	  * normals[point];
      }
    }
  }
}



template <int dim>
void numericalFluxUFromQ
(
 const dealii::Point<dim> referenceDirection,
 dealii::FullMatrix<double> & U_in_Q_in,
 dealii::FullMatrix<double> & U_in_Q_out,
 dealii::FullMatrix<double> & U_out_Q_in,
 dealii::FullMatrix<double> & U_out_Q_out,
 const dealii::FEValuesBase<dim> & fe_U_Self,
 const dealii::FEValuesBase<dim> & fe_U_Neig,
 const dealii::FEValuesBase<dim> & fe_Q_Self,
 const dealii::FEValuesBase<dim> & fe_Q_Neig,
 const double factor = 1.0){

  const unsigned int n_Q_dofs = fe_Q_Self.dofs_per_cell;
  const unsigned int n_U_dofs = fe_U_Self.dofs_per_cell;

  AssertDimension(n_Q_dofs ,fe_Q_Neig.dofs_per_cell);
  AssertDimension(n_U_dofs , fe_U_Neig.dofs_per_cell );  
  
  AssertDimension(U_in_Q_in.m(), n_U_dofs);
  AssertDimension(U_in_Q_in.n(), n_Q_dofs);
  AssertDimension(U_in_Q_out.m(), n_U_dofs);
  AssertDimension(U_in_Q_out.n(), n_Q_dofs);
  AssertDimension(U_out_Q_in.m(), n_U_dofs);
  AssertDimension(U_out_Q_in.n(), n_Q_dofs);
  AssertDimension(U_out_Q_out.m(), n_U_dofs);
  AssertDimension(U_out_Q_out.n(), n_Q_dofs);

  AssertDimension(fe_U_Self.get_fe().n_components(), 1);
  AssertDimension(fe_U_Neig.get_fe().n_components(), 1);
  AssertDimension(fe_Q_Self.get_fe().n_components(), dim);
  AssertDimension(fe_Q_Neig.get_fe().n_components(), dim);

  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);

  const auto & JxW = fe_U_Self.get_JxW_values();
  const auto & normals = fe_U_Self.get_normal_vectors(); 
  
  for(unsigned int point = 0;
      point < fe_Q_Self.n_quadrature_points;
      ++point){
    if(referenceDirection * normals[point] > 0.0){
      for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_Q_Self.dofs_per_cell; ++j){
	  U_in_Q_in(i,j)	    
	    += factor
	    * JxW[point]
	    * fe_U_Self[density].value(i,point)
	    * fe_Q_Self[flux].value(j,point)
	    * normals[point];
	}
      }
    }
    else {
      for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
	for(unsigned int j = 0; j < fe_Q_Neig.dofs_per_cell; ++j){
	  U_in_Q_out(i,j)
	    += factor
	    * JxW[point]
	    * fe_U_Self[density].value(i,point)
	    * fe_Q_Neig[flux].value(j,point)
	    * normals[point];
	}
      }
    }
  }
}

template <int dim>
void numericalFluxUFromQBoundary
(
 const dealii::Point<dim> referenceDirection,
 dealii::FullMatrix<double> & U_in_Q_in,
 const dealii::FEValuesBase<dim> & fe_U_Self,
 const dealii::FEValuesBase<dim> & fe_Q_Self,
 const double factor = 1.0){
 const unsigned int n_Q_dofs = fe_Q_Self.dofs_per_cell;
 const unsigned int n_U_dofs = fe_U_Self.dofs_per_cell;
  
   
  AssertDimension(U_in_Q_in.m(), n_U_dofs);
  AssertDimension(U_in_Q_in.n(), n_Q_dofs);

  AssertDimension(fe_U_Self.get_fe().n_components(), 1);
  AssertDimension(fe_Q_Self.get_fe().n_components(), dim);

  const dealii::FEValuesExtractors::Vector flux(0);
  const dealii::FEValuesExtractors::Scalar density(0);

  const auto & JxW = fe_Q_Self.get_JxW_values();
  const auto & normals = fe_Q_Self.get_normal_vectors();

  for(unsigned int point = 0;
      point < fe_Q_Self.n_quadrature_points;
      ++point){
    assert(referenceDirection * normals[point] > 0);
    for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_Q_Self.dofs_per_cell; ++j){
	U_in_Q_in(i,j)	    
	  += factor
	  * JxW[point]
	  * fe_U_Self[density].value(i,point)
	  * fe_Q_Self[flux].value(j,point)
	  * normals[point];
      }
    }
  }
}

template <int dim>
void BoundaryMassU
(
 const dealii::Point<dim> referenceDirection,
 dealii::FullMatrix<double>& U_in_U_in,
 const dealii::FEValuesBase<dim> & fe_U_Self,
 const double factor = 1.0){

  const dealii::FEValuesExtractors::Scalar density(0);
  const unsigned int n_U_dofs = fe_U_Self.dofs_per_cell;

  //AssertDimension(U_in_U_in.m(), n_U_dofs);
  //AssertDimension(U_in_U_in.n(), n_U_dofs);

  //AssertDimension(fe_U_Self.get_fe().n_components(), 1); 

  const auto & JxW = fe_U_Self.get_JxW_values();
  const auto & normals = fe_U_Self.get_normal_vectors();


  
  for(unsigned int point = 0;
      point < fe_U_Self.n_quadrature_points;
      ++point){
    assert(referenceDirection * normals[point] < 0);
    for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_U_Self.dofs_per_cell; ++j){
	U_in_U_in(i,j)
	  += factor
	  * JxW[point]
	  * fe_U_Self[density].value(i,point)
	  * fe_U_Self[density].value(j,point);
      }
    }
  }
}



template <int dim>
void BoundaryMassQ
(
 const dealii::Point<dim> referenceDirection,
 dealii::FullMatrix<double>& Q_in_Q_in,
 const dealii::FEValuesBase<dim> & fe_Q_Self,
 const double factor = 1.0){

  const unsigned int n_Q_dofs = fe_Q_Self.dofs_per_cell;

  AssertDimension(Q_in_Q_in.m(), n_Q_dofs);
  AssertDimension(Q_in_Q_in.n(), n_Q_dofs);

  AssertDimension(fe_Q_Self.get_fe().n_components(), dim);

  const dealii::FEValuesExtractors::Vector flux(0);
  
  const auto & JxW = fe_Q_Self.get_JxW_values();
  const auto & normals = fe_Q_Self.get_normal_vectors();

  for(unsigned int point = 0;
      point < fe_Q_Self.n_quadrature_points;
      ++point){
    assert(referenceDirection * normals[point] > 0);
    for(unsigned int i = 0; i < fe_Q_Self.dofs_per_cell; ++i){
      for(unsigned int j = 0; j < fe_Q_Self.dofs_per_cell; ++j){
	Q_in_Q_in(i,j)
	  += factor
	  * JxW[point]
	  * (fe_Q_Self[flux].value(i,point) * normals[point])
	  * (fe_Q_Self[flux].value(j,point) * normals[point]);
      }
    }
  }
}

template <int dim>
void U_MinusRHS_FromSource
(const dealii::Point<dim> referenceDirection,
 const dealii::FEValuesBase<dim> & fe_U_Self,
 dealii::Vector<heat::real> & U_Vector,
 const double factor = 1.0){

  const heat::SourceBodyValues<dim> sourceBV;
  const dealii::FEValuesExtractors::Scalar density(0);
  const auto & JxW = fe_U_Self.get_JxW_values();
  for(unsigned int point = 0; point < fe_U_Self.n_quadrature_points; ++point){
    for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
      U_Vector(i) += factor * JxW[point] * fe_U_Self[density].value(i,point)
	* sourceBV.value(fe_U_Self.quadrature_point(point));
    }
  }
}


template <int dim>
void Q_MinusRHS_FromDirPlus
(const dealii::Point<dim> referenceDirection,
 const dealii::FEValuesBase<dim> & fe_Q_Self,
 dealii::Vector<heat::real> & Q_Vector,
 const double factor = 1.0){

  const heat::DirichletBoundaryValues<dim> dirBC;

  const dealii::FEValuesExtractors::Vector flux(0);
  const auto & normals = fe_Q_Self.get_normal_vectors();
  const auto & JxW = fe_Q_Self.get_JxW_values();
  for(unsigned int point = 0; point < fe_Q_Self.n_quadrature_points; ++point){
    assert(referenceDirection * normals[point] > 0.0);
    for(unsigned int i = 0; i < fe_Q_Self.dofs_per_cell; ++i){
      Q_Vector(i) += factor * JxW[point] * fe_Q_Self[flux].value(i,point)
      * normals[point] * dirBC.value(fe_Q_Self.quadrature_point(point));
    }
  }  
}

template <int dim>
void U_MinusRHS_FromNeuMinus
(const dealii::Point<dim> referenceDirection,
 const dealii::FEValuesBase<dim> & fe_U_Self,
 dealii::Vector<heat::real> & U_Vector,
 const double factor = 1.0){
  const heat::NeumannBoundaryValues<dim> neuBC;
  const dealii::FEValuesExtractors::Scalar density(0);
  const auto & normals = fe_U_Self.get_normal_vectors();
  const auto & JxW = fe_U_Self.get_JxW_values();
  for(unsigned int point = 0; point < fe_U_Self.n_quadrature_points; ++point){
    assert(referenceDirection * normals[point] < 0.0);
    for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
      U_Vector(i) += factor * JxW[point] * fe_U_Self[density].value(i,point)
	* neuBC.value(fe_U_Self.quadrature_point(point)) * normals[point];
    }
  }
}

template<int dim>
void ConstraintUDirMinus
(const dealii::Point<dim> referenceDirection,
 const dealii::FEValuesBase<dim> & fe_U_Self,
 dealii::Vector<heat::real> & U_Vector,
 const double factor = 1.0){

  const heat::DirichletBoundaryValues<dim> dirBC;
  const auto & JxW = fe_U_Self.get_JxW_values();
  const auto & normals = fe_U_Self.get_normal_vectors();
  for(unsigned int point = 0; point < fe_U_Self.n_quadrature_points; ++point){
    assert(referenceDirection * normals[point] < 0.0);
    for(unsigned int i = 0; i < fe_U_Self.dofs_per_cell; ++i){
      U_Vector(i) += factor * JxW[point] * fe_U_Self.shape_value(i,point)
	* dirBC.value(fe_U_Self.quadrature_point(point));
    }
  }
}


template<int dim>
void ConstraintQNeuPlus
(
 const dealii::Point<dim> & referenceDirection,
 const dealii::FEValuesBase<dim> & fe_Q_Self,
 dealii::Vector<heat::real> & Q_Vector,
 const double factor = 1.0){

  const heat::NeumannBoundaryValues<dim> neuBC;
  const auto & JxW = fe_Q_Self.get_JxW_values();
  const auto & normals = fe_Q_Self.get_normal_vectors();
  const dealii::FEValuesExtractors::Vector flux(0);

  for(unsigned int point = 0; point < fe_Q_Self.n_quadrature_points; ++point){
    assert(referenceDirection * normals[point] > 0.0);
    for(unsigned int i = 0; i < fe_Q_Self.dofs_per_cell; ++i){
      Q_Vector(i) += (factor * JxW[point] * ( fe_Q_Self[flux].value(i,point) * normals[point] ) 
		      * ( neuBC.value(fe_Q_Self.quadrature_point(point)) * normals[point]));
    }
  }
}

template<int dim>
void
QJump( const dealii::FEValuesBase<dim> & fe_Q_self,
       const dealii::FEValuesBase<dim> & fe_Q_neig,
       const std::vector<std::vector<heat::real> > & Q_Vector_self,
       const std::vector<std::vector<heat::real> > & Q_Vector_neig,
       heat::real & jumpSquared)
{
  const auto & JxW = fe_Q_self.get_JxW_values();
  heat::real partialsum = 0;
  for(unsigned int dimension = 0; dimension < dim; ++dimension){
    for(unsigned int point = 0; point < fe_Q_self.n_quadrature_points; ++point){
      const heat::real diff = Q_Vector_self[1+dimension][point] - Q_Vector_neig[1+dimension][point];
      partialsum += JxW[point] * diff * diff;
    }
  }

  jumpSquared = partialsum;

}

template <int dim>
void
CellResidualQ( const dealii::FEValuesBase<dim> & fe_Q_self,
	       const std::vector<std::vector<dealii::Tensor<1,dim,heat::real >> > & Q_grad_Self,
	       const std::vector<heat::real> SourceTerms,
	       heat::real & normSquaredOfCellResidual){  
  
  //std::cout << "inside cell" << std::endl;

  const auto & JxW = fe_Q_self.get_JxW_values();  
  
  heat::real partialsum = 0;

  for(unsigned int point = 0; point < fe_Q_self.n_quadrature_points; ++point){
    heat::real DivQ = 0.0;
    for(unsigned int dimension = 0; dimension < dim; ++dimension){
      DivQ += Q_grad_Self[1+dimension][point][dimension];      
    }
    const heat::real diff = DivQ - SourceTerms[point];
    partialsum += JxW[point] * diff * diff;
  }

  normSquaredOfCellResidual = partialsum;

}




} //End Namepsace LDG
} //End Namespace LocalIntegrator
} //End Namespace heat


#endif
