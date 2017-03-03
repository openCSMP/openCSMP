#include "ArithmeticMean.h"

using namespace std;

namespace csmp {

/**
    Returns the average value of the supplied property into the result variable.
    For example, harmonic means of the permeability are returned into a region
    variable.
    
    For all positive data sets containing at least one pair of nonequal values, 
    the harmonic mean is always the least of the three means, 
    while the arithmetic mean is always the greatest of the three,
    and the geometric mean is always in between. 
*/
template<size_t dim, typename var>
ArithmeticMean<dim,var>::ArithmeticMean( const PropertyDatabase<dim>& p, const char* res_prop, const char* prop_to_avg ) 
      : Interrelation<dim>(p),
        I_( Interrelation<dim>::GlobalProperty(prop_to_avg) ),
        O_( Interrelation<dim>::GlobalProperty(res_prop) )
 {
    Interrelation<dim>::Name("ArithmeticMean");
    Interrelation<dim>::OutputCondition( O_, PLAIN );
    Interrelation<dim>::ResultProperty(res_prop);

    ///  @todo (2-D) (2-P) Use template parameter here!
    // making sure that the type of the result variable is consistent with the template parameter. 
    string  var_type(typeid(var).name());
    // ugly but robust
    if ( p.Type(res_prop) == SCALAR and var_type.find("Scalar") == 0 )
     throw Exception( ERROR, "ArithmeticArithmeticMean", "Result variable must match 'ScalarVariable' template parameter: ", res_prop );
    else if ( p.Type(res_prop) == VECTOR and var_type.find("Vector") == 0 )
     throw Exception( ERROR, "ArithmeticArithmeticMean", "Result variable must match 'VectorVariable' template parameter: ", res_prop );
    else if ( p.Type(res_prop) == TENSOR and var_type.find("Tensor") == 0 )
     throw Exception( ERROR, "ArithmeticArithmeticMean", "Result variable must match 'TensorVariable' template parameter: ", res_prop );
 }
 
 
 
template<size_t dim, typename var>
inline void ArithmeticMean<dim,var>::Calculate()
 {
    I_.AssignTo( prop_value_ );
    O_ = prop_value_;
 } 

template class ArithmeticMean<1U,ScalarVariable>;
template class ArithmeticMean<2U,ScalarVariable>;
template class ArithmeticMean<3U,ScalarVariable>;

} // csmp

