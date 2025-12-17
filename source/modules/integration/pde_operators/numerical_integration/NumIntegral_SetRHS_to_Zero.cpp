#include "NumIntegral_SetRHS_to_Zero.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_SetRHS_to_Zero<dim,CELL>::NumIntegral_SetRHS_to_Zero( const PropertyDatabase<dim>& pref, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,test)
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_SetRHS_to_Zero", test );

    if ( MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR ||
         MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE )
      throw csmp::Exception( ERROR, "MathOperatorRHS->NumIntegral_SetRHS_to_Zero<dim>::(constructor):",
                             test, "Dependent-variable must be a scalar variable placed on the nodes." );
                      
 } // end constructor


/**

Generats a RHS vector of zeros that is added to the global solution matrix.

The result is returned into the MathOperatorRHS vector<fT> V.

@test O.K.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_SetRHS_to_Zero<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // create a RHS vector of zeros
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );

} // end ComputeContribution


template class NumIntegral_SetRHS_to_Zero<1U,Element>;
template class NumIntegral_SetRHS_to_Zero<2U,Element>;
template class NumIntegral_SetRHS_to_Zero<3U,Element>;

template class NumIntegral_SetRHS_to_Zero<1U,Face>;
template class NumIntegral_SetRHS_to_Zero<2U,Face>;
template class NumIntegral_SetRHS_to_Zero<3U,Face>;

} // csmp
