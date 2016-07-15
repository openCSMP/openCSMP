#ifndef CSMP_CONCATENATE_H
#define CSMP_CONCATENATE_H

#include "Interrelation.h"

namespace csmp {

template<size_t dim, template <typename> class operation>
class Concatenate : public Interrelation< dim> {
  public:
    Concatenate( const PropertyDatabase<dim>& p, const char* resultProperty,
    			       const char* argumentProperty1, const char* argumentProperty2 );
                 
    Concatenate( const PropertyDatabase<dim>& p, const char* resultProperty,
                 std::vector<std::string>& arguments);	
                 		 
    ~Concatenate() {};
    
    void Calculate();

  private:
    Operand<dim>&  res_;
    operation<Operand< dim> > operate_;   
    bool result_is_argument;
    
    void Init( const char* resultProperty, std::vector<std::string>& arguments );
};

/**
 
@class Concatenate Concatenate "interrelations/Concatenate.h"
@author Adrian Burri
@date 2002

Concatenate is an interrelation which allows users to combine an
arbitrary number of variables. An STL binary function object is used to
define the operation used to combine the variables. For most 
applications, std::multiplies (for multiplication) and std::plus (for
summation) of the STL base library will be used, but in principle one
is free to define an arbitrary binary function. 

@attention for this interrelation to work all variables must have the 
same placement.

 
@section applicability Applicability
 
Use this interrelation when you need to combine different variables.   

 
@section implementation Implementation

The implementation makes use of the STL binary function classes
and the Operand class. Therefore, combining element and nodal variables
should be possible (not tested). The first operand is used to
initialize the result operand, which results in the following behaviour
when using non-commutative operands (here shown for std::minus):

@code
	result_property = first_property - (\sum other_properties)   
@endcode

If the result property is used as argument as well, first_property
will be result_property, no matter at which position it was named in the
list.

*/



/** Combines the arguments and stores result into resultProperty.
*/
template<size_t dim, template <typename> class operation>
inline void Concatenate<dim,operation>::Calculate()
 {
  typename std::map<std::string,Operand<dim> >::iterator 
  it(Interrelation<dim>::operand_list_.begin());

  if (it == Interrelation< dim>::result_) ++it;
 	
 	if (!result_is_argument) {
 	  res_ = static_cast<double64>(1.);
 	  res_ *= it->second;
 	  ++it;
 	}
 	
 	for (; it != Interrelation<dim>::operand_list_.end(); ++it) {
 		if (it != Interrelation<dim>::result_) {
 			res_ = operate_(res_, it->second);
 		}
 	}
    //res_.PrintValue();
     
 } // end specification of interrelation

} // csmp

#endif

