#include "NumIntegral_SetRHS_to_One.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<size_t dim,class SIMPLEX>
NumIntegral_SetRHS_to_One<dim,SIMPLEX>::~NumIntegral_SetRHS_to_One() {}


template<size_t dim,class SIMPLEX>
NumIntegral_SetRHS_to_One<dim,SIMPLEX>::NumIntegral_SetRHS_to_One( const PropertyDatabase<dim>& pref, const char* test )
  : MathOperatorRHS<dim>(pref,test)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_SetRHS_to_One", test );
 
    if ( MathOperatorRHS<dim>::TestOperandType() != SCALAR ||
         MathOperatorRHS<dim>::TestOperandPlacement() != NODE )
      throw csmp::Exception( ERROR, "MathOperatorRHS->NumIntegral_SetRHS_to_One<dim>::(constructor):",
                             test, "Dependent-variable must be a scalar variable placed on the nodes." );
                      
 } // end constructor


/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  void NumIntegral_SetRHS_to_One<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Generats a RHS vector of ones that is added to the global solution matrix. <p>

<H4>Input Arguments:</H4><!----------------------------------------------->

The element for which the RHS is generated <p>

<H4>Output Arguments &amp; Return Value</H4><!---------------------------->

The result is returned into the MathOperatorRHS vector<fT> V.<p>

<H4>Implementation:</H4><!------------------------------------------------>


<H4>Application:</H4><!--------------------------------------------------->

To assign stress boundary conditions to the boundary of the model. <p>

<!------------------------------------------------------------------------>
tested: O.K.  */
template<size_t dim,class SIMPLEX>
void NumIntegral_SetRHS_to_One<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
   // create a RHS vector of zeros
   MathOperatorRHS<dim>::RHS.resize(e.Nodes());
   for ( size_t i=0; i<e.Nodes(); i++ ) MathOperatorRHS<dim>::RHS[i] = 1.0/static_cast<double64>(e.N(i)->Parents());

} // end ComputeContribution


template class NumIntegral_SetRHS_to_One<1U,Element<1U> >;
template class NumIntegral_SetRHS_to_One<2U,Element<2U> >;
template class NumIntegral_SetRHS_to_One<3U,Element<3U> >;

template class NumIntegral_SetRHS_to_One<1U,Face<1U> >;
template class NumIntegral_SetRHS_to_One<2U,Face<2U> >;
template class NumIntegral_SetRHS_to_One<3U,Face<3U> >;

} // csmp
