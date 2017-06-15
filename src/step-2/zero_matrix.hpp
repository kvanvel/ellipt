#ifndef __deal2__zero_matrix_h
#define __deal2__zero_matrix_h

#include <deal.II/base/config.h>
#include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN

class ZeroMatrix : public Subscriptor
{
 public:
  typedef types::global_dof_index size_type;

  ZeroMatrix();

  ZeroMatrix(const size_type m, const size_type n);
  void reinit (const size_type m, const size_type n);

  size_type m() const;

  size_type n() const;

  template<class VECTOR1, class VECTOR2>
  void vmult(VECTOR1 &out,
	     const VECTOR2 &in) const;

  template<class VECTOR1, class VECTOR2>
  void vmult_add(VECTOR1 & out,
		 const VECTOR2 & in) const;

  template <class VECTOR1, class VECTOR2>
  void Tvmult(VECTOR1 &out,
	      const VECTOR2 &in) const;

  template <class VECTOR1, class VECTOR2>
  void Tvmult_add(VECTOR1 &out,
		  const VECTOR2 &in) const;

private:
  size_type n_Rows;
  size_type n_Cols;
  
  
};

// ------------------------- inline and template functions -------------
#ifndef DOXYGEN

inline
ZeroMatrix::ZeroMatrix()
  :
  n_Rows(0),
  n_Cols(0)
{}

inline
ZeroMatrix::ZeroMatrix(const size_type m, const size_type n)
  :
  n_Rows(m),
  n_Cols(n)
{}

inline
void
ZeroMatrix::reinit(const size_type m, const size_type n)
{
  n_Rows = m;
  n_Cols = n;
}

inline
ZeroMatrix::size_type
ZeroMatrix::m() const
{
  return n_Rows;
}

inline
ZeroMatrix::size_type
ZeroMatrix::n() const
{
  return n_Rows;
}


template <class VECTOR1, class VECTOR2>
inline
void
ZeroMatrix::vmult(VECTOR1       &out,
		  const VECTOR2 &in) const
{
  Assert(out.size() == n_Rows, ExcDimensionsMismatch(out.size(), n_Rows));
  Assert(in.size() == n_Cols, ExcDimensionsMismatch(in.size(), n_Cols));

  out = 0;
}

template <class VECTOR1, class VECTOR2>
inline
void
ZeroMatrix::vmult_add(VECTOR1       &out,
		      const VECTOR2 &in) const
{
  Assert(out.size() == n_Rows, ExcDimensionsMismatch(out.size(), n_Rows));
  Assert(in.size() == n_Cols, ExcDimensionsMismatch(in.size(), n_Cols));

  //Don't do anything.
}

template<class VECTOR1, class VECTOR2>
inline
void
ZeroMatrix::Tvmult(VECTOR1        &out,
		   const VECTOR2  &in) const
{
  Assert(out.size() == n_Cols, ExcDimensionsMismatch(out.size(), n_Cols));
  Assert(in.size() == n_Rows, ExcDimesnionsMismatch(in.size(), n_Rows));

  out = 0;
}

template<class VECTOR1, class VECTOR2>
inline
void
ZeroMatrix::Tvmult_add(VECTOR1        &out,
		       const VECTOR2  &in) const
{
  Assert(out.size() == n_Cols, ExcDimensionsMismatch(out.size(), n_Cols));
  Assert(in.size() == n_Rows, ExcDimesnionsMismatch(in.size(), n_Rows));
  // Don't do anything.
}

                             
#endif

DEAL_II_NAMESPACE_CLOSE

#endif

