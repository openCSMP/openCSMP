#include "ExtractVectorVariableComponent.h"
#include "FiniteElementManager.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ExtractVectorVariableComponent<dim>::ExtractVectorVariableComponent( const PropertyDatabase<dim>& p,
                                                                     const char* vec_var, 
                                                                     const char* scalar_var, 
                                                                     size_t comp )
      : Interrelation<dim>(p),
        V( Interrelation<dim>::GlobalProperty(vec_var) ),
        S( Interrelation<dim>::GlobalProperty(scalar_var) ),
        component(comp)
 {
    Interrelation<dim>::Name("ExtractVectorVariableComponent");
    Interrelation<dim>::OutputCondition( S, PLAIN );
    Interrelation<dim>::ResultProperty( scalar_var );
    
    if ( V.Type() != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableComponent::(constructor)", 
                                   "Extraction variable must be of VectorVariable type");
    if ( S.Type() != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableComponent::(constructor)", 
                                   "Extracted variable must be of ScalarVariable type");
 
    if ( V.Placement() == ELEMENT && S.Placement() != ELEMENT )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableComponent::(constructor)", 
                                   "Element variable cannot be extrapolated to elsewhere");
   
    if ( V.Placement() != ELEMENT && V.Placement() != S.Placement() )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableComponent::(constructor)", 
                                   "Incompatible placement of input variables");
   
    if ( comp > dim-1 )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableComponent::(constructor)", 
                                   "Desired vector component does not exist");
 }



template<uint32_t dim>
ExtractVectorVariableComponent<dim>::~ExtractVectorVariableComponent() {}


template class ExtractVectorVariableComponent<1U>;
template class ExtractVectorVariableComponent<2U>;
template class ExtractVectorVariableComponent<3U>;

} // csmp
