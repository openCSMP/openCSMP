#include "ExtractTensorVariableRow.h"
#include "FiniteElementManager.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ExtractTensorVariableRow<dim>::ExtractTensorVariableRow( const PropertyDatabase<dim>& p,
                                                         const char* tens_var, 
                                                         const char* vector_var, 
                                                         uint32_t r )
      : Interrelation<dim>(p),
        t( Interrelation<dim>::GlobalProperty(tens_var) ),
        v( Interrelation<dim>::GlobalProperty(vector_var) ),
        row(r)
 {
    Interrelation<dim>::Name("ExtractTensorVariableRow");
    Interrelation<dim>::OutputCondition( v, PLAIN );
    Interrelation<dim>::ResultProperty( vector_var );
    
    if ( t.Type() != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Extraction variable must be of TensorVariable type");
    if ( v.Type() != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Extracted variable must be of VectorVariable type");
 
    if ( t.Placement() != v.Placement() )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Incompatible placement of input variables");
   
    if ( row > dim-1 )
      throw csmp::Exception( FATAL_ERROR, "ExtractTensorVariableRow::(constructor)", 
                                   "Row out of range.");
 }


template class ExtractTensorVariableRow<1U>;
template class ExtractTensorVariableRow<2U>;
template class ExtractTensorVariableRow<3U>;

} // csmp
















