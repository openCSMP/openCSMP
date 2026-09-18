// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature.h"
#include "Index.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
    Constructor gets the finite element types that  NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature  needs to turn the InterFace into the equivalent
    of a volumetric element from the FiniteElementManager. Since InterFaces can either be triangles or quadrilaterals, their extruded volumetric
    pendants are prism and hexahedral elements.
    
    @attention the variable 'thickness' needs to be defined on the InterFace, else the pde operator will not work.
*/
template<uint32_t dim>
NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature( const Model<dim>& model,
                                                                                                                                  const char*       oper,
                                                                                                                                  const char*       basic,
                                                                                                                                  const char*       test)
    : MathOperatorLHS<dim,InterFace>(model.Database(),oper,basic,test),
    thickness_key_(model.Database().StorageKey("thickness"))//Benoit add

{
    MathOperatorLHS<dim,InterFace>::Name("NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature", oper, basic, test );
    
    // element interpolation function operand i
    if ( MathOperatorLHS<dim,InterFace>::BasicOperandPlacement() != NODE ||
        MathOperatorLHS<dim,InterFace>::BasicOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::(constructor)",
                              basic, "Operand (basic) must be a scalar property placed on the nodes." );

    // element interpolation function operand j
    if ( MathOperatorLHS<dim,InterFace>::TestOperandPlacement() != NODE ||
        MathOperatorLHS<dim,InterFace>::TestOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::(constructor)",
                              test, "Operand (test) must be a scalar property placed on the nodes." );

    // oper = material operand
    if ( MathOperatorLHS<dim,InterFace>::MaterialOperandPlacement() != ELEMENT ||
        MathOperatorLHS<dim,InterFace>::MaterialOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim,CELL>::GetOperands",
                              "material operand placements other than on the CELL are not handled yet. " );
}

/**
   Computes difference between the variable values at the topologically collocated nodes, including that in the intervening
   lower-dimensional region, using it to compute  the transfer coefficient (=coupling coefficient) at the time level t + delta t.
   
   The hydraulic conductivity (Operand value) is read from the higher-dimensional neighbors of the lower-dimensional region.
*/
template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::GetOperands( const InterFace<dim>& iface )
{

    assert( iface.HasInterveningElement() );
    
    // (the area of the interface is taken into account later)
    const auto n_nodes{ iface.InterveningElement()->Nodes() };
    transfer_coefficients_.resize( n_nodes * 3U );//WHY IS THIS NOT ZEROING VECTOR??

    for ( auto i{0U}; i<transfer_coefficients_.size(); i++ )
    {
        transfer_coefficients_[i] = 0.;
    }

    // 1. First-order finite-difference approximation of the interface fluxes
    // ----------------------------------------------------------------------
    if ( order1_FD_approximation_ ) {

        double thickness = iface.InterveningElement()->Read( thickness_key_ );


        // computing node-related fractions of the area of the intervening finite element
        assert( iface.FE()->MidSideNodes() == 0U );
        const double area3 = iface.InterveningElement()->Volume() / static_cast<double>( iface.InterveningElement()->Nodes() );
        double conductivity = iface.InterveningElement()->Read( MathOperatorLHS<dim,InterFace>::MaterialOperandKey() );

        for ( auto i{0U}; i<n_nodes; i++ )
        {
            Point<dim> in_mid = iface.MatchingN(i,INSIDE)->Coordinate() - iface.MatchingN(i,MIDDLE)->Coordinate();
            Point<dim> out_mid = iface.MatchingN(i,OUTSIDE)->Coordinate() - iface.MatchingN(i,MIDDLE)->Coordinate();


            double distance_in_mid (0.), distance_out_mid (0.);
            if(iface.MatchingN(i,INSIDE) == iface.MatchingN(i,OUTSIDE))
            {
                distance_in_mid = thickness/2;
                distance_out_mid = thickness/2;
            }
            else
            {
                distance_in_mid  = in_mid.Length();
                distance_out_mid = out_mid.Length();
                // We need to pull apart the SB by thickness for this to work!!!!!! Alternative: use thickness/2 instead
            }

            if(distance_in_mid > 0. && (iface.MatchingN(i,INSIDE) != iface.MatchingN(i,OUTSIDE)))
            {
                transfer_coefficients_[i] += conductivity / (distance_in_mid ) * area3;
            }


            if( distance_out_mid > 0. /*&& (iface.MatchingN(i,INSIDE) != iface.MatchingN(i,OUTSIDE))*/)
            {
                transfer_coefficients_[i + n_nodes] += conductivity / (distance_out_mid ) * area3;
            }


            // transfer into the intervening element
            transfer_coefficients_[i + n_nodes + n_nodes] = transfer_coefficients_[i] + transfer_coefficients_[i+n_nodes];
        }

        return;
    }

} // end GetOperands


/** Computes integral over shape function derivatives squared. These are the shape functions of the volumetric elements that match
    the lower-dimensional InterFace.

    @attention The connections could also be established using line elements, but this would require a more complicated
    accumulation process.
    
    @todo in GetOperands implement the collection of material operands from integration points otherwise second part of method will not work
*/
template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::ComputeContribution( const InterFace<dim>& iface )
{
    //cerr<<endl<<"NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::ComputeContribution";
    // this integral is only for numerically integrated isoparametric finite elements
    assert( iface.UsesLocalCoordinates() == true );
    assert( iface.HasInterveningElement() );

    // initialize output matrix, for creating the cross couplings between the face nodes with the intervening element
    // (ordering interface-inside-nodes, interface-outside nodes, intervening-element nodes
    const uint32_t n_total_nodes{ iface.Nodes() + iface.InterveningElement()->Nodes() },
        n_elmt_nodes{ iface.InterveningElement()->Nodes() };
    MathOperatorLHS<dim,InterFace>::LHS.Resize( n_total_nodes, n_total_nodes );
    MathOperatorLHS<dim,InterFace>::LHS.Zero();

    // populating the coupling matrix with the transfer coefficients
    // matrix diagonal
    for ( auto i{0U}; i<n_total_nodes; i++ )
        MathOperatorLHS<dim,InterFace>::LHS(i,i) = transfer_coefficients_[i];

    // off-diagonal terms depend on node-node connections
    // interface-inside-nodes with intervening-element nodes
    for ( auto i{0U}; i<n_elmt_nodes; i++ ) {
        // upper diagonal
        MathOperatorLHS<dim,InterFace>::LHS(i,i+n_elmt_nodes*2U) = -transfer_coefficients_[i];
        // lower diagonal
        MathOperatorLHS<dim,InterFace>::LHS(i+n_elmt_nodes*2U,i) = -transfer_coefficients_[i];
    }
    // interface-outside-nodes with intervening-element nodes
    for ( auto i{0U}; i<n_elmt_nodes; i++ ) {
        // upper diagonal
        MathOperatorLHS<dim,InterFace>::LHS(i+n_elmt_nodes*2U,i+n_elmt_nodes) = -transfer_coefficients_[i+n_elmt_nodes];
        // lower diagonal
        MathOperatorLHS<dim,InterFace>::LHS(i+n_elmt_nodes,i+n_elmt_nodes*2U) = -transfer_coefficients_[i+n_elmt_nodes];
    }

    // cerr<<endl;
    // for ( auto i{0U}; i<this->LHS.Rows(); i++ ){
    //     cerr<<"LHS("<<i<<"):";
    //     for ( auto j{0U}; j<this->LHS.Cols(); j++ )
    //         cerr<<this->LHS(i,j) << "  ";
    //     cerr << endl;
    // }

} // end ComputeContribution



/** AssignToGlobal matrix function
 
     Ordering of nodes in the vectors is  inner face nodes, outer face nodes, and intervening element nodes
     
    @author SKM - modified from LKT's version that eliminates Dirichlet constraints from G
*/
template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::AssignToGlobal( const InterFace<dim>& iface,
                                                                                        SparseMatrix& G,
                                                                                        vector<double>& pivotVector,
                                                                                        const vector<size_t>& DOF_indexes )
{
    const uint32_t n_total_nodes{ iface.Nodes() + iface.InterveningElement()->Nodes() },
        n_elmt_nodes{ iface.InterveningElement()->Nodes() };


    // vectors for mapping local to global node indices for test and basic operands
    vector<size_t> IDT( n_total_nodes ); // ii
    vector<size_t> IDB( n_total_nodes ); // jj

    const uint32_t n_elmt_nodes2{ n_elmt_nodes * 2 };

    // interface nodes + intervening element
    for ( auto i{0U}; i < n_elmt_nodes; i++ ) {
        IDT[i] = iface.MatchingN(i,INSIDE)->Idx();
        IDB[i] = IDT[i];
        IDT[i+n_elmt_nodes] = iface.MatchingN(i,OUTSIDE)->Idx();
        IDB[i+n_elmt_nodes] = IDT[i+n_elmt_nodes];
        IDT[i+n_elmt_nodes2] = iface.MatchingN(i,MIDDLE)->Idx();
        IDB[i+n_elmt_nodes2] = IDT[i+n_elmt_nodes2];
    }

    // ONLY SCALARS are used in the pde operator
    assert( this->BasicOperandType() == SCALAR );
    assert( this->TestOperandType() == SCALAR );
    
    // applying offset to interpolation function operand
    for (auto i{0U}; i < IDT.size(); i++) {
        IDT[i] += this->TestOperandOffset();
        IDT[i] = DOF_indexes[ IDT[i] ];
    }

    // applying offset to weighting function operand
    for (auto i{0U}; i < IDB.size(); i++) {
        IDB[i] += this->BasicOperandOffset();
        IDB[i] = DOF_indexes[ IDB[i] ];
    }

    // get local value from the given (local) node
    vector<double> nodal_values(IDB.size(),1.);
    {
        // interface and intervening element
        for ( auto nIdx{0U}; nIdx < n_elmt_nodes; ++nIdx ) {
            nodal_values[nIdx] = iface.MatchingN(nIdx,INSIDE)->Read(this->TestOperandKey());
            nodal_values[nIdx+n_elmt_nodes] = iface.MatchingN(nIdx,OUTSIDE)->Read(this->TestOperandKey());
            nodal_values[nIdx+n_elmt_nodes2] = iface.MatchingN(nIdx,MIDDLE)->Read(this->TestOperandKey());
        }
    }

    if ( this->multiply_accumulate_ )
    {
        throw csmp::Exception( ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::AssignToGlobal(InterFace):",
                              "Multiply Accumulate not supported yet");
    }
    /*
        Luat's factor_ = 1, always in the applications throughout code
    */
    else if ( this->add_accumulate_ || this->add_accumulate_later_ )
    {
        for (auto i{0U}; i < this->LHS.Rows(); i++) {
            if ( IDT[i] != NULL_IDX ) {
                for (auto j{0U}; j < this->LHS.Cols(); j++) {
                    if ( IDB[j] == NULL_IDX ) {
                        pivotVector[ IDT[i] ] -= this->LHS(i,j) * nodal_values[j] * this->factor_; // notice sign iface.N(j / this->TestOperandOffset())->Read(this->TestOperand())
                    }
                    else {
                        G.Add( IDT[i],
                              IDB[j],
                              this->LHS(i, j) * this->factor_ );
                    }
                }
            }
        }
    }
    else if ( this->subtract_accumulate_ || this->subtract_accumulate_later_ )
    {
        for ( auto i{0U}; i < this->LHS.Rows(); i++ ) {
            if ( IDT[i] != NULL_IDX ) {
                for (auto j{0U}; j < this->LHS.Cols(); j++ ) {
                    if ( IDB[j] == NULL_IDX ) {
                        pivotVector[ IDT[i] ] += this->LHS(i,j) * nodal_values[j] * this->factor_;   // notice the sign LHS(i, j) * e.N(j)->Read(TestOperandKey());
                    }
                    else {
                        G.Add( IDT[i],
                              IDB[j],
                              -this->LHS(i, j) * this->factor_ );
                    }
                }
            }
        }
    }
    else
        throw csmp::Exception( ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>::AssignToGlobal(InterFace):",
                              "accumulation instructions could not be parsed.");

} // end AssignToGlobal (InterFace)


template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<1U>;
template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<2U>;
template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<3U>;

} // csmp











