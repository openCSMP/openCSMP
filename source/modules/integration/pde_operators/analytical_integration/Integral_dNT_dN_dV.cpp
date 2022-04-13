#include "Integral_dNT_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim,class CELL>
Integral_dNT_dN_dV<dim,CELL>::Integral_dNT_dN_dV( const PropertyDatabase<dim>& pref,
                                                     const char* test )
  : MathOperatorRHS<dim>(pref,test),
    DN(2,3), DNT(3,2), UNITY(3,1)
{
    MathOperatorRHS<dim>::Name("Integral_dNT_dN_dV", test );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_dNT_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}



template<uint32_t dim,class CELL>
void Integral_dNT_dN_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );
    
    // setting up the unity vector
    UNITY.Resize(e.Nodes(),1);
    UNITY = 1.;

    // calculate the element contribution to RHS (lumping into vector format)
    DNT *= MathOperatorRHS<dim>::MTRL[0];
    DNT *= DN;
    DNT *= UNITY;
    
    // assigning element contribution & integrating the matrix
    double volume = e.Volume();
    
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    for ( auto i{0U}; i<e.Nodes(); i++ )
      MathOperatorRHS<dim>::RHS[i] = DNT(i,0) * volume;

} // end ComputeContribution



template class Integral_dNT_dN_dV<1U,Element<1U> >;
template class Integral_dNT_dN_dV<2U,Element<2U> >;
template class Integral_dNT_dN_dV<3U,Element<3U> >;

template class Integral_dNT_dN_dV<1U,Face<1U> >;
template class Integral_dNT_dN_dV<2U,Face<2U> >;
template class Integral_dNT_dN_dV<3U,Face<3U> >;

} // csmp
