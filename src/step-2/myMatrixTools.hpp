#ifndef _MYMATRIXTOOLS_H__
#define _MYMATRIXTOOLS_H__

#include "globals.hpp"
#include "EllipticProblem.hpp"
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/sparse_matrix_ez.h>


namespace heat{
namespace myMatrixTools{

template<class FULLMATRIX>
void
ChopMatrix(FULLMATRIX & In,
	   heat::real threshold = 1e-10);


template<class SPARSEMATRIX,class STREAM>
void
PrintMatrixMarket(const SPARSEMATRIX & In,
		  STREAM & out,
		  heat::real threshold = 1e-10
		  );

// Implementation

template<class FULLMATRIX>
void
ChopMatrix(FULLMATRIX & In,
	   heat::real threshold)
{
  for(unsigned int i = 0; i < In.m(); ++i){
    for(unsigned int j = 0; j < In.n(); ++j){
      if(fabs( In(i,j) ) < threshold){
	In(i,j) = 0;
      }
    }
  }
}


template<class SPARSEMATRIX,class STREAM>
void
PrintMatrixMarket(const SPARSEMATRIX & In,
		  STREAM & out,
		  heat::real threshold
		  )
{
  Assert( In.m() != 0,  dealii::ExcNotInitialized() );
  Assert( threshold >0, dealii::ExcMessage("Negative threshold!") );
  
  // Print the header
  out << "%%MatrixMarket matrix coordinate real general\n";
  const auto nnz = In.n_actually_nonzero_elements(threshold);
  out << In.m() << ' ' << In.n() << ' ' << nnz << '\n';
  
  // Print the body
  for(unsigned int i = 0; i < In.m(); ++i){
    for(auto it = In.begin(i); it != In.end(i); ++it){
      const auto value = it->value();
      if(fabs(value) >  threshold){
	const unsigned int j = it->column();
	out << i + 1 << ' ' << j + 1 << ' ' << value << '\n';
      }
    }
  }  
}


template<class STREAM>
void
PrintMatrixMarket(const dealii::SparseMatrixEZ<heat::real> & In,
		  STREAM & out,
		  heat::real threshold = 1e-10
		  )
{
  if(In.m() == 0){
    return;
  }
  
  //Assert( In.m() != 0, dealii::ExcNotInitialized() );
  Assert( threshold >0, dealii::ExcMessage("Negative threshold!") );
  
  out << "%%MatrixMarket matrix coordinate real general\n";  
  const auto nnz = In.n_nonzero_elements();
  out << In.m() << ' ' << In.n() << ' ' << nnz << std::endl;
  
  for(unsigned int i = 0; i < In.m(); ++i){

    for(auto it = In.begin(i); it != In.end(i); ++it){
      const auto value = it->value();
      if(fabs(value) >  threshold){
	const unsigned int j = it->column();
	out << i + 1 << ' ' << j + 1 << ' ' << value << std::endl;
      }
    }
  }  
}


template<class SPARSEMATRIX>
void
DistillMatrix(SPARSEMATRIX & matrixIN,
	      typename dealii::SparsityPattern & sparsityPattern,	      
	      heat::real threshold = 1e-10)
{
  
  Assert( matrixIN.m() != 0, dealii::ExcNotInitialized() );
  Assert( threshold > 0, dealii::ExcMessage("Negative threshold!") );  

  
  const dealii::SparseMatrix<heat::real>::size_type n_rows = matrixIN.m();
  const dealii::SparseMatrix<heat::real>::size_type n_cols = matrixIN.n();
  
  std::vector<std::vector<dealii::SparseMatrix<heat::real>::size_type > > col_indices(n_rows );
  std::vector<std::vector<std::pair<dealii::SparseMatrix<heat::real>::size_type, heat::real> > >
    colVal_indices(n_rows );

  for(unsigned int i = 0; i < n_rows; ++i){    
    for(auto it = matrixIN.begin(i); it != matrixIN.end(i); ++it){
      const auto value = it->value();
      if(fabs(value) > threshold){
	const auto col = it->column();
	const std::pair<dealii::SparseMatrix<heat::real>::size_type, heat::real> tempPair {col, value};
	col_indices[i].push_back(col);
	colVal_indices[i].push_back(tempPair);
      }
    }
  }
  //const bool optimizeDiagonal = true;
  sparsityPattern.copy_from(n_rows,
			    n_cols,
			    colVal_indices.begin(), 
			    colVal_indices.end()
			    //optimizeDiagonal
			    );
  
  matrixIN.clear();
  matrixIN.reinit(sparsityPattern);
  matrixIN.copy_from(colVal_indices.begin(),
		     colVal_indices.end() );
  
  sparsityPattern.compress();  
}



} // End namespace heat
} // End namespace myMatrixTools

#endif
