#ifndef CONSTANT_FACTOR_H
#define CONSTANT_FACTOR_H

#include "Interrelation.h"

namespace csmp {

template<size_t dim, template <typename> class operation>
class ConstantFactor : public Interrelation< dim> {
  private:
    Operand<dim>&  res_;
    Operand<dim>&  arg_;
    ScalarVariable factor_, val_;
    
    operation<ScalarVariable > operate_;
    bool  res_prop_equal_to_arg_prop_;
    
  public:
    ConstantFactor( const PropertyDatabase<dim>& p, const char* resultProperty,
                    const char* argumentProperty, const double factor);
    			          
    ~ConstantFactor() {};
    
    void ChangeFactor( double factor );
    
    void Calculate();
};

/**

@class ConstantFactor ConstantFactor "interrelations/ConstantFactor.h"
@author S.K. Matthaei
@author S. Geiger
@author S. Roberts
@date 2002


ConstantFactor is an interrelation which allows users to combine an
arbitrary number of variables. An STL binary function object is used to
define the operation used to combine the variables. For most 
applications, std::multiplies (for multiplication) and std::plus (for
summation) of the STL base library will be used, but in principle one
is free to define an arbitrary binary function.

 
@section applicability Applicability

Use this interrelation when you need to combine different variables.
 
 
@section implementation Implementation
 
The implementation makes use of the STL binary function classes
and the Operand class. Therefore, combining element and nodal variables
should be possible (not tested). The first operand is used to
initialize the result operand, which results in the following behaviour
when using non-commutative operands (here shown for std::minus):

@code
result_property = first_property - (\sum other_properties).
@endcode
*/



} // csmp

#endif

