#include "Integral_var_NT_rhsop_N_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** basic and test variable exchanged so that the row in the matrix is indicated by the basic variable (as it is
    for the lhs operators)
*/
template<uint32_t dim, template<uint32_t> class CELL>
Integral_var_NT_rhsop_N_dV<dim,CELL>::Integral_var_NT_rhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                                  const char* oper,
                                                                  const char* basic,
                                                                  const char* test,
                                                                  const char* var,
                                                                  const double prefactor )
  : MathOperatorRHS<dim,CELL>(pref,oper,basic),
    basic_(pref.Parameter(test)),
    var_(pref.Parameter(var)),
    prefactor_(prefactor)
 {
    string name = "Integral_var_NT_rhsop_N_dV ";
    name += test;
    name += ", ";
    name += var;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorRHS<dim,CELL>::Name(cname, oper, basic );
    // testing the Operands 
    if ( basic_.key.place != NODE || basic_.key.type != SCALAR )
       throw csmp::Exception( ERROR, "Integral_var_NT_rhsop_N_dV<dim>::(constructor)",
                              basic, "Basic variable must be a scalar property placed on the node." );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT and
         MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != REGION )
       throw csmp::Exception( ERROR, "Integral_var_NT_rhsop_N_dV<dim>::(constructor)", 
                              oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
       throw csmp::Exception( ERROR, "Integral_var_NT_rhsop_N_dV<dim>::(constructor)", 
                              test, "Dependent variable must be a scalar property placed on the nodes." );
     
    if ( var_.key.place != NODE || var_.key.type != SCALAR )
      throw csmp::Exception( ERROR, "Integral_var_NT_rhsop_N_dV<dim>::(constructor)", 
                             var, "Additional variable must be a scalar property placed on the nodes." );
                   
    
 }





/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_var_NT_rhsop_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
   // this integral is only for analytically integrated finite elements
   assert( e.FE()->UsesLocalCoordinates() == false );
  
  // reading basic
  e.NodePropertyVector(basic_.key, basic_var_);
   // reading oper and var
  e.Read(MathOperatorRHS<dim,CELL>::MaterialOperandKey(), op_);
  e.NodePropertyVector(var_.key, vvar_);
}





/**
 
Computes the volume (area) integral over the basic function products 
multiplied with the Operand and stores the result in the test function
part of the right hand side vector.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_var_NT_rhsop_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
  if ( !MathOperatorRHS<dim,CELL>::LumpedFormulation() ) { // consistent formulation
    if (e.FE_Type() != LINEAR_TRIANGLE)
       throw csmp::Exception(ERROR, 
                      "Integral_var_NT_lhsop_N_dV::ComputeContribution",
                      "Consistent formulation only provided for elements of type LinearTriangle");
     ComputeIntegral(e);
  } else { // lumped formulation
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
     fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );
  
    DenseMatrix< DM_MIN> mat;
    e.IntegralNN(mat);
    
    for (uint32_t i = 0; i < e.Nodes(); ++i) {
      for (uint32_t k = 0; k < e.Nodes(); ++k) {
        MathOperatorRHS<dim,CELL>::RHS[i] += mat(i, k) * vvar_[k]();
      }
    }
    for (uint32_t i = 0; i < e.Nodes(); ++i) {
      MathOperatorRHS<dim,CELL>::RHS[i] *= op_() * prefactor_ * basic_var_[i]();
    }
  }
} // end ComputeContribution




template<uint32_t dim, template<uint32_t> class CELL>
void Integral_var_NT_rhsop_N_dV<dim,CELL>::ComputeIntegral( const CELL<dim>& e ) {
  assert(e.Nodes() == 3);
  assert(vvar_.size() == 3);

  DenseMatrix< DM_MIN> elMat;
  MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
  fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );

  elMat(0, 0) = 6.0*vvar_[0]() + 2.0*vvar_[1]() + 2.0*vvar_[2]();
  elMat(0, 1) = 2.0*vvar_[0]() + 2.0*vvar_[1]() + vvar_[2]();
  elMat(0, 2) = 2.0*vvar_[0]() + vvar_[1]() + 2.0*vvar_[2]();
  elMat(1, 0) = elMat(0, 1);
  elMat(1, 1) = 2.0*vvar_[0]() + 6.0*vvar_[1]() + 2.0*vvar_[2]();
  elMat(1, 2) = vvar_[0]() + 2.0*vvar_[1]() + 2.0*vvar_[2]();
  elMat(2, 0) = elMat(0, 2);
  elMat(2, 1) = elMat(1, 2);
  elMat(2, 2) = 2.0*vvar_[0]() + 2.0*vvar_[1]() + 6.0*vvar_[2]();
  
  for (uint32_t i = 0; i < e.Nodes(); i++) 
    for (uint32_t j = 0; j < e.Nodes(); j++) 
      MathOperatorRHS<dim,CELL>::RHS[i] += elMat(i,j) * basic_var_[j]();
  
  for (uint32_t i = 0; i < e.Nodes(); ++i) {
    MathOperatorRHS<dim,CELL>::RHS[i] *= e.Volume() * prefactor_ * op_() / 60.;
  }
}



template class Integral_var_NT_rhsop_N_dV<2U,Element>;
template class Integral_var_NT_rhsop_N_dV<3U,Element>;

template class Integral_var_NT_rhsop_N_dV<2U,Face>;
template class Integral_var_NT_rhsop_N_dV<3U,Face>;

} // csmp
