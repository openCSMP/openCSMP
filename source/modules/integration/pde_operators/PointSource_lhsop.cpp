// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PointSource_lhsop.h"
#include "ErrorHandler.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
PointSource_lhsop<dim,CELL>::PointSource_lhsop( const PropertyDatabase<dim>& pref,
                                                const char* oper, const char* basic, const char* test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    SRC_(8)
 {
    MathOperatorLHS<dim,CELL>::Name("PointSource_lhsop", oper, test );
 
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
     if ( MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR ||
          MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE )
       csmp_error.Note( ERROR, "MathOperatorLHS->PointSource_lhsop<dim>::(constructor)",
                          test, "must be a scalar variable placed on the node" ); 
                          
     if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR ||
          MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != NODE )
       csmp_error.Note( ERROR, "MathOperatorLHS->PointSource_lhsop<dim>::(constructor)",
                          oper, "must be a scalar variable placed on the node" ); 
 }




/// reads the values of basic operand from element nodes
template<uint32_t dim, template<uint32_t> class CELL>
void PointSource_lhsop<dim,CELL>::GetOperands( const CELL<dim>& e )
   { 
      e.NodePropertyVector( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), SRC_ );
       
      // taking into account that the point source contributes to several elements
      for ( auto i{0}; i<e.Nodes(); i++ )
        SRC_[i]() /= static_cast<double>(e.N(i)->Parents());
        
   } // end GetOperands





/**
 
Applies point (nodal) source terms = Neumann boundary conditions to 
any kind of FEM models. The contribution is divided by the number of 
parent elements of each node because it will be accumulated by each one
of of them. This is already done in the GetOperands() method. 

@param e The current element from which the nodal values are accumulated.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void PointSource_lhsop<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
   MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
   MathOperatorLHS<dim,CELL>::LHS.Zero();
   for ( auto i{0}; i<e.Nodes(); i++ )
     MathOperatorLHS<dim,CELL>::LHS(i,i) = SRC_[i]();
     
} // end ComputeContribution


template class PointSource_lhsop<1U>;
template class PointSource_lhsop<2U>;
template class PointSource_lhsop<3U>;

template class PointSource_lhsop<1U,Face>;
template class PointSource_lhsop<2U,Face>;
template class PointSource_lhsop<3U,Face>;

} // csmp
