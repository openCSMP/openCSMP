#include "ExtractTensorVariableComponent.h"
#include "FiniteElementManager.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ExtractTensorVariableComponent<dim>::ExtractTensorVariableComponent( const PropertyDatabase<dim>& p,
                                                                     const char* tens_var, 
                                                                     const char* scalar_var, 
                                                                     uint32_t i, uint32_t j )
      : Interrelation<dim>(p),
        T( Interrelation<dim>::GlobalProperty(tens_var) ),
        S( Interrelation<dim>::GlobalProperty(scalar_var) ),
        comp_i(i),
        comp_j(j)
 {
    Interrelation<dim>::Name("ExtractTensorVariableComponent");
    Interrelation<dim>::OutputCondition( S, PLAIN );
    Interrelation<dim>::ResultProperty( scalar_var );
    
    if ( T.Type() != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableComponent::(constructor)", 
                                   "Extraction variable must be of TensorVariable type");
    if ( S.Type() != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableComponent::(constructor)", 
                                   "Extracted variable must be of ScalarVariable type");
 
    if ( T.Placement() == ELEMENT && S.Placement() != ELEMENT )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableComponent::(constructor)", 
                                   "Element variable cannot be extrapolated to elsewhere");
   
    if ( T.Placement() != ELEMENT && T.Placement() != S.Placement() )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableComponent::(constructor)", 
                                   "Incompatible placement of input variables");
   
    if ( comp_i > dim-1 )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableComponent::(constructor)", 
                                   "Desired vector component 'i' does not exist");

    if ( comp_j > dim-1 )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableComponent::(constructor)", 
                                   "Desired vector component 'j' does not exist");
 }



template class ExtractTensorVariableComponent<1U>;
template class ExtractTensorVariableComponent<2U>;
template class ExtractTensorVariableComponent<3U>;

} // csp














