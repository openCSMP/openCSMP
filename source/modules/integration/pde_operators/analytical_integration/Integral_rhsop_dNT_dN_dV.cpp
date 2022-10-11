#include "Integral_rhsop_dNT_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL,typename var>
Integral_rhsop_dNT_dN_dV<dim,CELL,var>::Integral_rhsop_dNT_dN_dV( const PropertyDatabase<dim>& pref,
                                                                  const char*             oper,
                                                                  const char*             basic,
                                                                  const char*             test,
                                                                  double prefactor )
  : MathOperatorRHS<dim,CELL>(pref,oper,basic),
    DN(2,3), DNT(3,2),
    basic_(pref.Parameter(test)),
    prefactor_(prefactor)
{
    string name = "Integral_rhsop_dNT_dN_dV ";
    name += test;
    char* cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorRHS<dim,CELL>::Name( cname, oper, test );

    // testing the Operands
    if ( basic_.key.place != NODE || basic_.key.type != SCALAR )
    throw csmp::Exception( ERROR, "Integral_NT_rhsop_N_dV<dim>::(constructor)",
                           basic, "Basic variable must be a scalar property placed on the node." );
                   
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT and MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_rhsop_dNT_dN_dV::(constructor)", 
                           oper, "Operand must be placed on the element or group.");

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_rhsop_dNT_dN_dV::(constructor)", 
                           test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}



/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL,typename var>
void Integral_rhsop_dNT_dN_dV<dim,CELL,var>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

   // reading oper and basic
   MathOperatorRHS<dim,CELL>::GetOperands(e);
   e.NodePropertyVector( basic_.key, basic_var_ );
}



template<uint32_t dim, template<uint32_t> class CELL,typename var>
void Integral_rhsop_dNT_dN_dV<dim,CELL,var>::ComputeContribution( const CELL<dim>& e )
 {
    e.dN( DN );
    
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );

    // calculate the element contribution to RHS (lumping into vector format)
    DNT *= MathOperatorRHS<dim,CELL>::MTRL[0];
    DNT *= DN;

    // assigning element contribution & integrating the matrix
    fill(MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.);
    double volume = e.Volume();
    
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    for (auto i{0U}; i < e.Nodes(); ++i) {
      for (auto j{0U}; j < e.Nodes(); ++j) {
           MathOperatorRHS<dim,CELL>::RHS[i] += DNT(i,j) * basic_var_[j]() * volume * prefactor_;
        }
    }

} // end ComputeContribution




template class Integral_rhsop_dNT_dN_dV<1U,Element,ScalarVariable>;
template class Integral_rhsop_dNT_dN_dV<2U,Element,ScalarVariable>;
template class Integral_rhsop_dNT_dN_dV<3U,Element,ScalarVariable>;

template class Integral_rhsop_dNT_dN_dV<1U,Face,ScalarVariable>;
template class Integral_rhsop_dNT_dN_dV<2U,Face,ScalarVariable>;
template class Integral_rhsop_dNT_dN_dV<3U,Face,ScalarVariable>;

} // csmp
