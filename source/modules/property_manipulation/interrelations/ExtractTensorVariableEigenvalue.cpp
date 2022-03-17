#include "ExtractTensorVariableEigenvalue.h"
#include "FiniteElementManager.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ExtractTensorVariableEigenvalue<dim>::ExtractTensorVariableEigenvalue( const PropertyDatabase<dim>& p,
                                                                       const char* tens_var, 
                                                                       const char* scalar_var, 
                                                                       size_t e )
      : Interrelation<dim>(p),
        t( Interrelation<dim>::GlobalProperty(tens_var) ),
        s( Interrelation<dim>::GlobalProperty(scalar_var) ),
        eval(e)
 {
    Interrelation<dim>::Name("ExtractTensorVariableEigenvalue");
    Interrelation<dim>::OutputCondition( s, PLAIN );
    Interrelation<dim>::ResultProperty( scalar_var );
    
    if ( t.Type() != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Extraction variable must be of TensorVariable type");
    if ( s.Type() != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Extracted variable must be of ScalarVariable type");
 
    if ( t.Placement() != s.Placement() )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Incompatible placement of input variables");
   
    if ( eval > dim-1 )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Eigenvalue number out of range.");
 }


template class ExtractTensorVariableEigenvalue<2U>;
template class ExtractTensorVariableEigenvalue<3U>;

} // csmp
















