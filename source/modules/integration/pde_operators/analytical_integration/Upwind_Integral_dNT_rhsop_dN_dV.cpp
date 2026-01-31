#include "Upwind_Integral_dNT_rhsop_dN_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
    basic and test variable exchanged so that the line is indicated by the basic variable (as it is
    for the lhs operators)
*/
template<uint32_t dim, template<uint32_t> class CELL>
Upwind_Integral_dNT_rhsop_dN_dV<dim,CELL>::Upwind_Integral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& pref,
                                                                            const char* oper,
                                                                            const char* basic,
                                                                            const char* test,
                                                                            const char* upwind,
                                                                            const char* trigger,
                                                                            const double prefactor )
  : MathOperatorRHS<dim,CELL>(pref,oper,basic),
    DN(3,3),
    DNT(3,3),
    basic_(pref.Parameter(test)),
    upwind_(pref.Parameter(upwind)),
    trigger_(pref.Parameter(trigger)),
    prefactor_(prefactor)
 {
    string name = "Upwind_Integral_dNT_rhsop_dN_dV ";
    name += test;
    name += ", ";
    name += upwind;
    name += ", ";
    name += trigger;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorRHS<dim,CELL>::Name(cname, oper, basic );

    if ( basic_.key.place != NODE || basic_.key.type != SCALAR )
      throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV<dim>::(constructor)",
                   basic, "Basic variable must be a scalar property placed on the node." );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT || 
         MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV<dim>::(constructor)", 
                   oper, "Operand must be a scalar property placed on the element." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
       throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
     
    if ( upwind_.key.place != NODE || upwind_.key.type != SCALAR )
       throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV<dim>::(constructor)", 
                   upwind, "Upwind variable must be a scalar property placed on the nodes." );
                   
    if ( trigger_.key.place != NODE || trigger_.key.type != SCALAR )
       throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV<dim>::(constructor)", 
                   trigger, "Trigger variable must be a scalar property placed on the nodes." );
 }







template<uint32_t dim, template<uint32_t> class CELL>
void Upwind_Integral_dNT_rhsop_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    // reading operand
    MathOperatorRHS<dim,CELL>::GetOperands(e);
     
    // basic, upwind and trigger
    e.NodePropertyVector(basic_.key, basic_var_);
    e.NodePropertyVector(upwind_.key, upwind_var_);
    e.NodePropertyVector(trigger_.key, trigger_var_);
}








template<uint32_t dim, template<uint32_t> class CELL>
void Upwind_Integral_dNT_rhsop_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );

    // consistent formulation
    if ( !MathOperatorRHS<dim,CELL>::LumpedFormulation() ) {
      e.dN( DN );
      // transpose the shape function derivative matrix
      DNT.Resize(e.Nodes(),dim); 
      DN.Transposed( DNT );
    
      // calculate the element contribution to RHS (E_OP is the Basic Operand)
      DNT *= MathOperatorRHS<dim,CELL>::MTRL[0];
      DNT *= DN;
      
      
      // calculate upwinding coefficients and multiply them with operand matrix   
      for (uint32_t i = 0; i < e.Nodes(); ++i) {
          for (uint32_t j = 0; j < e.Nodes(); ++j) {;
              if (i != j) {
                  const double decision = DNT(i, j)*(trigger_var_[j]() - trigger_var_[i]());
                  if      (decision > 0) DNT(i, j) *= upwind_var_[i]();
                  else if (decision < 0) DNT(i, j) *= upwind_var_[j]();
                  else                   DNT(i, j) *= 0.5*(upwind_var_[i]() + upwind_var_[j]());
              }
              else                     DNT(i, j) = static_cast<double>(0.0);
          }
      }
    
      for ( auto i{0U}; i < e.Nodes(); ++i ) {
          DNT(i, i) = -DNT.RowSum(i);
      }
      
      // the matrix is contracted into a vector by multiplying with the basis vector
      for (uint32_t i{0U}; i < e.Nodes(); ++i)
        for (uint32_t j{0U}; j < e.Nodes(); ++j)
          MathOperatorRHS<dim,CELL>::RHS[i] += DNT(i,j) * basic_var_[j]();
          
      // scaling with prefactor
      for (uint32_t i = 0; i < e.Nodes(); ++i)
        MathOperatorRHS<dim,CELL>::RHS[i] *= e.Volume() * prefactor_;
    }
    // lumped formulation  
    else {
      throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV::(ComputeContribution)",
                             "Lumped formulation not allowed for this operator"); 
    }
} // end ComputeContribution


template class Upwind_Integral_dNT_rhsop_dN_dV<1U,Element>;
template class Upwind_Integral_dNT_rhsop_dN_dV<2U,Element>;
template class Upwind_Integral_dNT_rhsop_dN_dV<3U,Element>;

template class Upwind_Integral_dNT_rhsop_dN_dV<1U,Face>;
template class Upwind_Integral_dNT_rhsop_dN_dV<2U,Face>;
template class Upwind_Integral_dNT_rhsop_dN_dV<3U,Face>;

} // csmp
