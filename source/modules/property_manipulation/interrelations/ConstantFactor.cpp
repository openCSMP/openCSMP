#include "ConstantFactor.h"
#include "ErrorHandler.h"

// Includes for function objects
#include <functional>

using namespace std;

namespace csmp {


/**

Constructor which is there for convenience, since most of the time to 
variables are being combined.
*/
template<size_t dim, template <typename> class operation>
ConstantFactor<dim,operation >::ConstantFactor( const PropertyDatabase<dim>& p, 
                                                const char* resultProperty,
                                                const char* argumentProperty,
                                                const double64 factor )
      : Interrelation<dim>(p),
        res_(Interrelation<dim>::GlobalProperty(resultProperty)),
        arg_(Interrelation<dim>::GlobalProperty(argumentProperty)),
        factor_(PLAIN,factor),
        res_prop_equal_to_arg_prop_(false)
 {   
    Interrelation<dim>::Name("ConstantFactor");
     // setting up the variables for the formula
    Interrelation<dim>::OutputCondition(res_, PLAIN);
    Interrelation<dim>::ResultProperty(resultProperty);

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( res_.Key() == arg_.Key() ) res_prop_equal_to_arg_prop_ = true;

    if ( res_.Type() != SCALAR or arg_.Type() != SCALAR )
      csmp_error.notice( CSMP_FATAL_ERROR, "ConstantFactor<dim,operation >(constructor)",
                                      "Interrelation only works for scalar variables" );

    if ( res_.Placement() != arg_.Placement() )
      throw Exception( CSMP_ERROR, "ConstantFactor", "in this interrelation, all variables must have the same placement");
 }
 
 
 
 /** Changes the factor.
 */
template<size_t dim, template <typename> class operation>
void ConstantFactor<dim,operation>::ChangeFactor(double64 factor) {
  factor_ = factor;
}



/** Combines the arguments and stores result into resultProperty
 */
template<size_t dim, template <typename> class operation>
void ConstantFactor<dim,operation>::Calculate()
 {
    if ( !res_prop_equal_to_arg_prop_ ) {
        res_ = static_cast<double64>(1.);
        res_ *= arg_;
      }
    res_.AssignTo(val_);
    res_ = operate_(val_, factor_);
     
 } // end specification of interrelation

 

template class ConstantFactor< 1, multiplies>;
template class ConstantFactor< 2, multiplies>;
template class ConstantFactor< 3, multiplies>;
template class ConstantFactor< 1, divides>;
template class ConstantFactor< 2, divides>;
template class ConstantFactor< 3, divides>;
template class ConstantFactor< 1, minus>;
template class ConstantFactor< 2, minus>;
template class ConstantFactor< 3, minus>;
template class ConstantFactor< 1, plus>;
template class ConstantFactor< 2, plus>;
template class ConstantFactor< 3, plus>;

}
