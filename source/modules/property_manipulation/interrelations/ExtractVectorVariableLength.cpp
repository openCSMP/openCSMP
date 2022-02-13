#include "ExtractVectorVariableLength.h"
#include "FiniteElementManager.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ExtractVectorVariableLength<dim>::ExtractVectorVariableLength( const PropertyDatabase<dim>& p,
                                                          const char* vec_var, 
                                                          const char* to_scalar_var )
      : Interrelation<dim>(p),
        V( Interrelation<dim>::GlobalProperty(vec_var) ),
        S( Interrelation<dim>::GlobalProperty(to_scalar_var) )
 {
    Interrelation<dim>::Name("ExtractVectorVariableLength");
    Interrelation<dim>::OutputCondition( S, PLAIN );
    Interrelation<dim>::ResultProperty( to_scalar_var );
    
    if ( V.Type() != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableLength::(constructor)", 
                                   "Extraction variable must be of VectorVariable type");
    if ( S.Type() != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableLength::(constructor)", 
                                   "Extracted variable must be of ScalarVariable type");
 
    if ( V.Placement() == ELEMENT && S.Placement() != ELEMENT )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableLength::(constructor)", 
                                   "Element variable cannot be extrapolated to elsewhere");
   
    if ( V.Placement() != ELEMENT && V.Placement() != S.Placement() )
      throw csmp::Exception( FATAL_ERROR, "ExtractVectorVariableLength::(constructor)", 
                                   "Incompatible placement of input variables");
 }



template<uint32_t dim>
ExtractVectorVariableLength<dim>::~ExtractVectorVariableLength() {}


template class ExtractVectorVariableLength<1U>;
template class ExtractVectorVariableLength<2U>;
template class ExtractVectorVariableLength<3U>;

} // csmp
