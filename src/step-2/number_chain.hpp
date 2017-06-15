#ifndef __number_chain_
#define __number_chain_

#include <deal.II/base/config.h>
#include <deal.II/lac/exceptions.h>

namespace heat{

template<typename... Ts> class NumberChain;

template<>
class NumberChain<> : public dealii::Subscriptor
{
public:
  unsigned int m() const
  {return 0;}
  unsigned int n() const
  {return 0;}

  template<typename VECTOR1, typename VECTOR2>
  void Tvmult(VECTOR1 & Destination,
	      const VECTOR2 &Source) const
  {}

  template<typename VECTOR1, typename VECTOR2>
  void vmult(VECTOR1 & Destination,
	     const VECTOR2 &Source) const
  {}
  

};



template<typename First, typename ... Rest>
class NumberChain<First, Rest...> : public NumberChain<Rest...>  
{
  const dealii::SmartPointer<const First> member;
  const unsigned int length;
public:
  NumberChain(First const& f, Rest const& ... rest)
    :
    NumberChain<Rest...>(rest...),
    member(&f),
    length(1 + sizeof...(Rest))
  {}

 
  const dealii::SmartPointer<const First>
  head() const
  {
    return member;
  }

  NumberChain<Rest...> const& rest() const
  {
    return *this;
  }

  template<typename VECTOR1, typename VECTOR2>
  void Tvmult(VECTOR1 & Destination,
	      const VECTOR2 &Source) const;

  template<typename VECTOR1, typename VECTOR2>
  void vmult(VECTOR1 & Destination,
	     const VECTOR2 &Source) const;

  unsigned int m() const;
  unsigned int n() const;

};


template<typename First, typename... Rest>
unsigned int NumberChain<First, Rest...>::m() const
{
  return head()->m();
}


template<typename First, typename... Rest>
unsigned int NumberChain<First, Rest...>::n() const
{
  if(1 != length){
    return
      this->rest().n();
  } else {
    return head()->n();
  }  
}


template<typename First, typename... Rest>
template<typename VECTOR1, typename VECTOR2>
void
NumberChain<First, Rest...>::Tvmult(VECTOR1 & Destination,
				   const VECTOR2 &Source) const
{
  
  if(1 != length){
    VECTOR1 temp1;

    temp1.reinit( this->head()->n() );

    this->head()->Tvmult(temp1, Source);
    rest().Tvmult(Destination, temp1);
  }
  else{
    this->head()->Tvmult(Destination, Source);
  }


}


template<typename First, typename... Rest>
template<typename VECTOR1, typename VECTOR2>
void
NumberChain<First, Rest...>::vmult(VECTOR1 & Destination,
				   const VECTOR2 &Source) const
{
  
  if(1 != length){
    VECTOR1 temp1;

    temp1.reinit( this->rest().m() );

    this->rest().vmult(temp1, Source);
    this->head()->vmult(Destination, temp1);
  }
  else{
    this->head()->vmult(Destination, Source);
  }


}





}
#endif
