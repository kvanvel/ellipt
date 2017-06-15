#ifndef __deal2__matrix_chain_h
#define __deal2__matrix_chain_h

#include <deal.II/base/config.h>
#include <deal.II/lac/exceptions.h>


namespace heat{

template<typename... Args> class MatrixChain;

template<>
class MatrixChain<> : public dealii::Subscriptor {};

template<typename T, typename... Trest>
class MatrixChain<T, Trest...> : private MatrixChain<Trest...>
{
public:

  MatrixChain(const T & A, const Trest & ... B);
  
  const dealii::SmartPointer<const T, MatrixChain<T, Trest...> > A;
  MatrixChain<Trest...> B;  

};



template<typename T, typename... Trest>
MatrixChain<T, Trest...>::MatrixChain(const T &A_in, const Trest & ... rest)
  :
  A(&A_in),
  B(rest...)
{
}

// const dealii::SmartPointer<const dealii::SparseMatrix<double> >
// const dealii::SparseMatrix<double>*
// const dealii::SmartPointer<const dealii::SparseMatrix<double>, heat::MatrixChain<dealii::SparseMatrix<double> > >’
// const dealii::SparseMatrix<double>*

// template<typename... Args> class myTest;

// template<> class myTest<>
// {
// };

// template<typename T1, typename... TRest> class myTest<T1, TRest...>
//   : public myTest<TRest...>
// {
// public:
//   T1 Data;

//   myTest<TRest...>& Rest()
//   {
//     return *this;    
//   }
// };



// template<typ

// template<typename NUMBER, typename... Args>
// class NumberChain
// {
// public:
//   NumberChain(Args... most, NUMBER &last);

// private:
//   NUMBER* Last;
//   //const dealii::SmartPointer<const MATRIXA, MatrixChain<MATRIXA, Args...> A;
//   NumberChain<Args...>&  Most;
//   //const dealii::SmartPointer<const MatrixChain<Args...>, MatrixChain<Args...,MatrixA> > most;
// };




// template<typename MATRIXA, typename...Args>
// NumberChain<MATRIXA, Args...>
// ::NumberChain(Args... most, MATRIXA& A_In)
//   :
//   Last(&A_In),
//   Most(most...)
// {}
                             


}
#endif

