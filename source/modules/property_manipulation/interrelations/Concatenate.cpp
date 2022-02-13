#include "Concatenate.h"

// Includes for function objects
#include <numeric>
#include <functional>
#include "Exception.h"

using namespace std;

namespace csmp {


/** Constructor which is there for convenience, since most of the time to
variables are being combined.
*/
template<uint32_t dim, template <typename> class operation>
Concatenate<dim,operation >::Concatenate( const PropertyDatabase<dim>& p,
                                          const char* resultProperty,
                                          const char* argumentProperty1,
                                          const char* argumentProperty2)  			 									                        
    : Interrelation<dim>(p),
      res_(Interrelation<dim>::GlobalProperty(resultProperty)),
      result_is_argument(false)
 {
    Interrelation<dim>::Name(" Concatenate");
    Interrelation<dim>::OutputCondition(res_, PLAIN);
    Interrelation<dim>::ResultProperty(resultProperty);
 
    vector<string> args;
    args.push_back(argumentProperty1);
    args.push_back(argumentProperty2);
    Init(resultProperty, args);
    
    // safety
    if ( res_.Placement() != p.Placement(argumentProperty1) or 
         p.Placement(argumentProperty1) != p.Placement(argumentProperty2) )
      throw Exception( ERROR, "Concatenate", "in this interrelation, all variables must have the same placement");
 }
 
 
 
 
 
/** Constructor. The argument properties are passed as a vector of strings and
processed in the order given therein.
*/
template<uint32_t dim, template <typename> class operation>
Concatenate<dim,operation>::Concatenate( const PropertyDatabase<dim>& p,  
                                         const char* resultProperty,
									                       vector<string>& arguments )
	: Interrelation<dim>(p),
	  res_(Interrelation<dim>::GlobalProperty(resultProperty)),
	  result_is_argument(false)
 {
 	  Interrelation<dim>::Name(" Concatenate");
    Interrelation<dim>::OutputCondition(res_, PLAIN);
    Interrelation<dim>::ResultProperty(resultProperty);
 	  Init(resultProperty, arguments);
   
    for (  vector<string>::const_iterator
           it=arguments.begin(); it!=arguments.end(); it++ )
      if ( res_.Placement() != p.Placement( (*it).c_str() )  )
        throw Exception( ERROR, "Concatenate", "in this interrelation, all variables must have the same placement");
 }






template<uint32_t dim, template <typename> class operation>
void Concatenate<dim,operation>::Init( const char* resultProperty, vector<string>& arguments ) 
{
  for (vector<string>::iterator it = arguments.begin(); it != arguments.end(); ++it) {
 		Interrelation<dim>::GlobalProperty(it->c_str());
 		if (!strcmp(resultProperty, it->c_str())) {
 		  result_is_argument = true;
 		}
 	}
}

template class Concatenate< 1, multiplies>;
template class Concatenate< 2, multiplies>;
template class Concatenate< 3, multiplies>;
template class Concatenate< 1, divides>;
template class Concatenate< 2, divides>;
template class Concatenate< 3, divides>;
template class Concatenate< 1, plus>;
template class Concatenate< 2, plus>;
template class Concatenate< 3, plus>;
template class Concatenate< 1, minus>;
template class Concatenate< 2, minus>;
template class Concatenate< 3, minus>;

}
