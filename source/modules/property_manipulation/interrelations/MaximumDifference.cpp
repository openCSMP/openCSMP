#include "MaximumDifference.h"

using namespace std;

namespace csmp {

template<size_t dim>
MaximumDifference<dim>::MaximumDifference( const PropertyDatabase<dim>& p,
                                           const char* variable1, const char* variable2 ) 
      : Interrelation<dim>(p),
        Var1( Interrelation<dim>::GlobalProperty(variable1) ),
        Var2( Interrelation<dim>::GlobalProperty(variable2) ),
        max_difference(0.)
 {
    Interrelation<dim>::Name("MaximumDifference");
     // setting up the variables for the formula
    Interrelation<dim>::OutputCondition( Var2, PLAIN );
    Interrelation<dim>::ResultProperty(variable2);
 }



template<size_t dim>
void MaximumDifference<dim>::Reset()
 {
    max_difference = 0.;
     
 } // end 



template<size_t dim>
double  MaximumDifference<dim>::Value() const
 {
    return max_difference;
     
 } // end 



template class MaximumDifference<1U>;
template class MaximumDifference<2U>;
template class MaximumDifference<3U>;

} // csmp
