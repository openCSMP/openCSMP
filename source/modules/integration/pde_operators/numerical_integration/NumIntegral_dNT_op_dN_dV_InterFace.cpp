#include "NumIntegral_dNT_op_dN_dV_InterFace.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
    Constructor gets the finite element types that  NumIntegral_dNT_op_dN_dV_InterFace  needs to turn the InterFace into the equivalent
    of a volumetric element from the FiniteElementManager. Since InterFaces can either be triangles or quadrilaterals, their extruded volumetric
    pendants are prism and hexahedral elements.
    
    @attention the variable 'thickness' needs to be defined on the InterFace, else the pde operator will not work.
*/
template<uint32_t dim>
NumIntegral_dNT_op_dN_dV_InterFace<dim>::NumIntegral_dNT_op_dN_dV_InterFace( const Model<dim>& model,
                                                                              const char*       oper,
                                                                              const char*       basic,
                                                                              const char*       test )
  : MathOperatorLHS<dim,InterFace>(model.Database(),oper,basic,test),
    thi_key_(model.Database().StorageKey("thickness")),
    pris_ptr_(model.Mesh().FiniteElements().E( ISOPARAMETRIC_LINEAR_PRISM ) ),
    hexa_ptr_(model.Mesh().FiniteElements().E( ISOPARAMETRIC_LINEAR_HEXAHEDRON ) )
{
    MathOperatorLHS<dim,InterFace>::Name("NumIntegral_dNT_op_dN_dV_InterFace", oper, basic, test );
    
    if ( MathOperatorLHS<dim,InterFace>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,InterFace>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_InterFace<dim>::(constructor)",
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,InterFace>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,InterFace>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_InterFace<dim>::(constructor)",
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}




/**
    Reads the material operand and the interface thickness from the current interface cell.
    
    @attention this is a first draft of the method that only works for piecewise constant scalar variables.
    
    @todo fix for vector or tensor material operands which may be needed if the permeability is a vector or tensor
    @todo use harmonic mean of the permeabilities of the parent elements on the inside and outside.
    @todo implement for integration point variables
*/
template<uint32_t dim>
void NumIntegral_dNT_op_dN_dV_InterFace<dim>::GetOperands( const InterFace<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    if ( MathOperatorLHS<dim,InterFace>::MaterialOperandPlacement() != ELEMENT or
         MathOperatorLHS<dim,InterFace>::MaterialOperandPlacement() != INTER_FACE )
    throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_InterFace<dim,CELL>::GetOperands",
                                  "material operand placements other than on the CELL are not handled yet. " );

    // only if the property is an InterFace property something is done here
     MathOperatorLHS<dim,InterFace>::MTRL[0].Resize(dim,dim);
     MathOperatorLHS<dim,InterFace>::MTRL[0].Zero();
     // identity matrix scaled with the element thickness
     MathOperatorLHS<dim,InterFace>::MTRL[0].AssignToDiagonal( (thickness_ = e.Read( thi_key_ )) );
  
     // the operand variables are taken as the mean of the variables on the adjacent higher dim elements
     // (note that higher-dim parents are guaranteed to exist on inside and outside of InterFace)
     if ( MathOperatorLHS<dim,InterFace>::MaterialOperandType() == SCALAR ) {
          MathOperatorLHS<dim,InterFace>::MTRL[0](0,0) = (e.InnerParent()->Read( MathOperatorLHS<dim,InterFace>::MaterialOperandKey() ) +
                                                     e.OuterParent()->Read( MathOperatorLHS<dim,InterFace>::MaterialOperandKey() ) ) / 2. ;
          if constexpr ( dim != 1U ) MathOperatorLHS<dim,InterFace>::MTRL[0](1,1) = MathOperatorLHS<dim,InterFace>::MTRL[0](0,0);
          if constexpr ( dim == 3U ) MathOperatorLHS<dim,InterFace>::MTRL[0](2,2) = MathOperatorLHS<dim,InterFace>::MTRL[0](0,0);
       }
     else throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_InterFace<dim,CELL>::GetOperands",
                                "'vector' or 'tensor' material operands are not handled yet. " );
                                
     /* TODO: later on, together with integration point variables
     if ( MathOperatorLHS<dim>::MaterialOperandType() == VECTOR ) {
          VectorVariable<dim>  vc;
          e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
          MathOperatorLHS<dim,CELL>::MTRL[0](0,0) = vc[0];
          if ( dim != 1U ) MathOperatorLHS<dim,CELL>::MTRL[0](1,1) = vc[1];
          if ( dim == 3U ) MathOperatorLHS<dim,CELL>::MTRL[0](2,2) = vc[2];
       }
     if ( MathOperatorLHS<dim>::MaterialOperandType() == TENSOR ) {
          TensorVariable<dim>  ts;
          e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), ts );
          MathOperatorLHS<dim>::MTRL[0] = ts;
       }
     */

 } // end GetOperands




/**
    Spatially separates opposing nodes in the direction of the unit normal of the InterFace.
    It follows that the extrusion direction is perpendicular to the InterFace and the extrusion amount is given by the InterFace "thickness.
    The method writes the new coordinates to the XY matrix of the InterFace cell, following the node numbering convention
    of the corresponding (resulting) prism or hexahedral element.
    
    @todo Before using this in computations test that the finite element swap works and that the node numbering of the contribution matches that of the prism or hexa.
 */
template<uint32_t dim>
void NumIntegral_dNT_op_dN_dV_InterFace<dim>::InitialiseCoordinateMatrix( const InterFace<dim>& iface,
                                                                          DenseMatrix<DM_MIN>& XYX ) const
 {
     // retrieving the Interface normal for the extrusion of the interface from the INSIDE to the outside
     Point<dim> unrml = iface.UnitNormal();
     // scaling unit normal by element thickness
     unrml[0U] *= thickness_;
     if constexpr ( dim == 2U ) unrml[1U] *= thickness_;
     if constexpr ( dim == 3U ) unrml[2U] *= thickness_;
     
     assert( iface.FE()->OrderOfShapeFunctions() == 1U );
     
     // node coordinate extrusion; procedure does not require a know of the InterFace element
     // -------------------------------------------------------------------------------------
     // the triangular face is extruded into prism element retaining the inside node positions
     // resizing the coordinate matrix that lives in the FE base class
     XYX.Resize( iface.Nodes(), dim );

     // assigning the inside nodes
     // --------------------------
     const auto n_side_nodes{ iface.FE()->Nodes() };
     for ( uint32_t i{0U}; i<n_side_nodes; i++ )
       XYX.AssignRow( i, iface.N(i,INSIDE)->Coordinate() );

     // assigning the outside nodes extruding their position by interface thickness along the unit normal
     // -------------------------------------------------------------------------------------------------
      uint32_t node_idx = 0;
      const auto total_nodes = iface.Nodes();

      // Create a view of indices from (total_nodes - 1) down to 0
      auto reversed_indices = views::iota(0u, total_nodes) | views::reverse;

      for (auto i : reversed_indices) {
          // Check if the node exists before accessing to prevent crashes
          if (auto* node_ptr = iface.N(static_cast<uint32_t>(i), OUTSIDE)) {
              Point<dim> extr_coord = node_ptr->Coordinate() + unrml;
              XYX.AssignRow(node_idx++, extr_coord);
          }
      }

 } // end InitialiseCoordinateMatrix







/** Computes integral over shape function derivatives squared. These are the shape functions of the volumetric elements that match
    the lower-dimensional InterFace.

    @attention The connections could also be established using line elements, but this would require a more complicated
    accumulation process.
    
    @todo in GetOperands implement the collection of material operands from integration points otherwise second part of method will not work
*/
template<uint32_t dim>
void NumIntegral_dNT_op_dN_dV_InterFace<dim>::ComputeContribution( const InterFace<dim>& iface )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( iface.UsesLocalCoordinates() == true );
    InitialiseCoordinateMatrix( iface, iface.FE()->XY );

    // initialize output matrix, InterFace already has the right number of nodes to match prism or hexa
    MathOperatorLHS<dim,InterFace>::LHS.Resize( iface.Nodes(), iface.Nodes() );
    MathOperatorLHS<dim,InterFace>::LHS.Zero();

    // thus far this only works for constant element or interface properties
    const bool piecewise_constant_material( this->MaterialOperandPlacement() == ELEMENT or
                                            this->MaterialOperandPlacement() == INTER_FACE or
                                            this->MaterialOperandPlacement() == REGION or
                                            this->MaterialOperandPlacement() == FACE );
    assert( piecewise_constant_material );
    
    // 0. Finding the volumetric finite element that corresponds to the Face
    // ---------------------------------------------------------------------
    FiniteElement* fptr = iface.FE();
    if ( isTriangular( iface.FE_Type() ) )
      // ugly way to get to the base class!
      FiniteElementPolicy<dim,InterFace>(iface).Assign( pris_ptr_ );
    else
      FiniteElementPolicy<dim,InterFace>(iface).Assign( hexa_ptr_ );
   
    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    if ( piecewise_constant_material )
      {
        // looping over the integration points of the prism or hexahedral element
        for ( uint32_t i{0U}; i<iface.IntegrationPoints(); i++ ) {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (DN is already in global coordinates)
             const double detJ = iface.dN_AtIntegrationPoint( DN_, i, SCALAR );
             // transposing B -> BT
             DN_.Transposed( DNT_ );
             // multiply  BT . MTRL
             DNT_ *= MathOperatorLHS<dim,InterFace>::MTRL[0];
             // multiplying BT . B
             DNT_ *= DN_;
             // multiplying with determinant and weights
             DNT_ *= iface.WeightAtIntegrationPoint(i) * detJ;
             // accumulating ME Gauss point integral contributions into element
             // contribution to global conductance matrix
             MathOperatorLHS<dim,InterFace>::LHS += DNT_;
          }
        return;
      }
      
    // NODE or ELEMENT_INTEGRATION_POINT material placements
    for ( uint32_t i{0U}; i<iface.FE()->IntegrationPoints(); i++ ) {
         const double detJ = iface.dN_AtIntegrationPoint( DN_, i, SCALAR );
         DN_.Transposed( DNT_ );
         DNT_ *= MathOperatorLHS<dim,InterFace>::MTRL[i];
         DNT_ *= DN_;
         DNT_ *= iface.WeightAtIntegrationPoint(i) * detJ;
         MathOperatorLHS<dim,InterFace>::LHS += DNT_;
      }

   // 3. reassigning the FE pointer
   FiniteElementPolicy<dim,InterFace>(iface).Assign( fptr );

} // end ComputeContribution


// FOR DEBUGGING
//cout <<"\nNumIntegral_dNT_op_dN_dV_InterFace: on InterFace "<< e.Idx() << endl;
//MathOperatorLHS<dim>::LHS.Out();

// for ( auto i{0U}; i<this->LHS.Rows(); i++ )
//   if ( this->LHS(i,i) < numeric_limits::epsilon() ) 
//     cout <<"\nNumIntegral_dNT_op_dN_dV_InterFace: zero element in diagonal of element matrix.";


template class NumIntegral_dNT_op_dN_dV_InterFace<1U>;
template class NumIntegral_dNT_op_dN_dV_InterFace<2U>;
template class NumIntegral_dNT_op_dN_dV_InterFace<3U>;

} // csmp











