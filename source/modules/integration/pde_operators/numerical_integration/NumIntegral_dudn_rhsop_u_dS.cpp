#include "NumIntegral_dudn_rhsop_u_dS.h"
#include "Index.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
    Constructor gets the finite element types that  NumIntegral_dudn_rhsop_u_dS  needs to turn the InterFace into the equivalent
    of a volumetric element from the FiniteElementManager. Since InterFaces can either be triangles or quadrilaterals, their extruded volumetric
    pendants are prism and hexahedral elements.
    
    @attention the variable 'thickness' needs to be defined on the InterFace, else the pde operator will not work.
*/
template<uint32_t dim>
NumIntegral_dudn_rhsop_u_dS<dim>::NumIntegral_dudn_rhsop_u_dS( const Model<dim>& model,
                                                               const char* interface_oper,
                                                               const char* elmt_oper,
                                                               const char* test )
:  MathOperatorRHS<dim,InterFace>(model.Database(),interface_oper,test),
   diffusivity_key_(model.Database().StorageKey("hydraulic diffusivity")),
   conductivity_key_(model.Database().StorageKey("conductivity")),
   delta_t_(numeric_limits<double>::quiet_NaN())
{
    MathOperatorRHS<dim,InterFace>::Name("NumIntegral_dudn_rhsop_u_dS", interface_oper, test );
    
    // element interpolation function operand j
    if ( MathOperatorRHS<dim,InterFace>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,InterFace>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dudn_rhsop_u_dS<dim>::(constructor)",
                      test, "Operand (test) must be a scalar property placed on the nodes." );

    // oper = material operand
    if ( MathOperatorRHS<dim,InterFace>::MaterialOperandPlacement() != INTER_FACE ||
         MathOperatorRHS<dim,InterFace>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_dudn_rhsop_u_dS<dim,InterFace>::GetOperands",
                                  "material operand placements other than on the INTER_FACE are not handled yet. " );
}



template<uint32_t dim>
NumIntegral_dudn_rhsop_u_dS<dim>*  NumIntegral_dudn_rhsop_u_dS<dim>::clone() const
 {
    return new NumIntegral_dudn_rhsop_u_dS<dim>(*this);
 }
      



/**
   Computes difference between the variable values at the topologically collocated nodes, including that in the intervening
   lower-dimensional region, using it to compute  the transfer coefficient (=coupling coefficient) at the time level t + delta t.
   
   The hydraulic conductivity (Operand value) is read from the higher-dimensional neighbors of the lower-dimensional region.
*/
template<uint32_t dim>
void NumIntegral_dudn_rhsop_u_dS<dim>::GetOperands( const InterFace<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.HasInterveningElement() );
    
    // retrieving conductivity from inner and outer higher-dimensional parent elements
    // (their values for the intervening element are ignored as they are accumulated separately)
    const double conductivity_inner = e.InnerParent()->Read( conductivity_key_ );
//    const double diffusivity_inner  = e.InnerParent()->Read( diffusivity_key_ );
    const double conductivity_outer = e.OuterParent()->Read( conductivity_key_ );
//    const double diffusivity_outer  = e.OuterParent()->Read( diffusivity_key_ );

    // (the area of the interface is taken into account later)
    const auto n_nodes{ e.InterveningElement()->Nodes() };
    transfer_coefficients_.resize( n_nodes * 3U );


    // 1. First-order finite-difference approximation of the interface fluxes
    // ----------------------------------------------------------------------
    if ( order1_FD_approximation_ ) {
         double     val_intervening = e.PropertyValueAtBaryCenter( MathOperatorRHS<dim,InterFace>::TestOperandKey() );
         Point<dim> unrml( e.InterveningElement()->UnitNormal() );
         
         // gradient of test function-variable is computed from the facet-normal difference
         // -------------------------------------------------------------------------------
         const double val_inner = e.InnerParent()->PropertyValueAtBaryCenter( MathOperatorRHS<dim,InterFace>::TestOperandKey() );
         // separation of barycenters of inner and intervening cells
         Point<dim> bctr_distance = e.InnerParent()->BaryCenter() - e.InterveningElement()->BaryCenter();
         // projecting this distance on the normal of the intervening element
         const double distance_to_inner_bctr = fabs( dotProduct( unrml, bctr_distance ) );
         // the same gradient is used for all nodes
         double grad_var0 = (val_inner - val_intervening) / distance_to_inner_bctr;
         for ( auto i{0U}; i<n_nodes; i++ )
           // source term due to pressure difference (sink if the gradient is positive)
           transfer_coefficients_[i] = -conductivity_inner * grad_var0;
 
cout <<"\ninterface "<< e.Idx() <<": inner dist-n: "<< distance_to_inner_bctr <<", grad p: "<< grad_var0;
 
         const double val_outer = e.OuterParent()->PropertyValueAtBaryCenter( MathOperatorRHS<dim,InterFace>::TestOperandKey() );
         // separation of barycenters of outer and intervening cells
         bctr_distance = e.OuterParent()->BaryCenter() - e.InterveningElement()->BaryCenter();
         // projecting this distance on the normal of the intervening element
         const double distance_to_outer_bctr = fabs( dotProduct( unrml, bctr_distance ) );
         // the same gradient is used for all nodes
         grad_var0 = (val_outer - val_intervening) / distance_to_outer_bctr;
         for ( auto i{0U}; i<n_nodes; i++ )
           // source term due to pressure difference (sink if gradient is positive)
           transfer_coefficients_[i+n_nodes] = -conductivity_outer * grad_var0;

         // transfer into the intervening element
         for ( auto i{0U}; i<n_nodes; i++ )
           // source term due to pressure difference (opposite sign than for inside and outside)
           transfer_coefficients_[i + n_nodes + n_nodes] = -(transfer_coefficients_[i] + transfer_coefficients_[i+n_nodes]);
 
cout <<", outer dist-n: "<< distance_to_outer_bctr <<", grad p: "<< grad_var0;
  
         return; // EXIT FROM THE GETOPERANDS
      }

 } // end GetOperands









    // and the values of the transported variable at the interface faces and the intervening element
/*
    vector<ScalarVariable> node_var_iface, node_var_elmt;
    e.NodePropertyVector( MathOperatorLHS<dim,InterFace>::TestOperandKey(), node_var_iface );
    e.InterveningElement()->NodePropertyVector( MathOperatorLHS<dim,InterFace>::TestOperandKey(), node_var_elmt );
*/
    // computing transfer coefficient values using the variable values at the interface nodes
    // and the intervening element (inside - middle = 3, outside - middle = 3, 6 in total.
    // --------------------------------------------------------------------------------------
/*
    // inside face and intervening element (nodes have the same indices)
    for ( auto i{0U}; i<n_nodes; i++ ) {
         // nodes of inner face come first
         double delta_var0 = e.N(i,INSIDE)->Read( MathOperatorRHS<dim,InterFace>::TestOperandKey() ) -
                             e.InterveningElement()->N(i)->Read( MathOperatorRHS<dim,InterFace>::TestOperandKey() );
         // if there is no gradient, the existing hydraulic conductivity is set as the transfer coefficient
         if ( fabs(delta_var0) < numeric_limits<double>::epsilon() )
           // @note sources (not sinks) are negative if added to LHS
           transfer_coefficients_[i] = -conductivity_inner;
         else {
              // using 1D analytical solution to compute pressure gradient normal to inner interface at t+dt
              double gradient = Grad_Var_AtX0( delta_var0, diffusivity_inner, delta_t_ );
              // computing the transfer coefficient that will give rise to this node-to-node flux at t+dt
              transfer_coefficients_[i] = -(delta_var0 * conductivity_inner) / gradient;
           }
      }
    // outside face and intervening element (nodes have different indices)
    for ( auto i{0U}; i<n_nodes; i++ ) {
         // nodes of outer face
         double delta_var0 = e.MatchingN(i,OUTSIDE)->Read( MathOperatorRHS<dim,InterFace>::TestOperandKey() ) -
                             e.InterveningElement()->N(i)->Read( MathOperatorRHS<dim,InterFace>::TestOperandKey() );
         if ( fabs(delta_var0) < numeric_limits<double>::epsilon() )
           transfer_coefficients_[i+n_nodes] = 0.; // conductivity_outer;
         else {
              // using 1D analytical solution to compute pressure gradient normal to outer interface at t+dt
              double gradient = Grad_Var_AtX0( delta_var0, diffusivity_outer, delta_t_ );
              // computing the transfer coefficient that will give rise to this node-to-node flux at t+dt
              transfer_coefficients_[i+n_nodes] = 0.; // (delta_var0 * conductivity_outer) / gradient;
           }
      }
  
    // TODO: flux from country rock into the fracture must still be applied as source/sink on the fracture element

*/




/** Computes integral over shape function derivatives squared. These are the shape functions of the volumetric elements that match
    the lower-dimensional InterFace.

    @attention The connections could also be established using line elements, but this would require a more complicated
    accumulation process.
    
    @todo in GetOperands implement the collection of material operands from integration points otherwise second part of method will not work
*/
template<uint32_t dim>
void NumIntegral_dudn_rhsop_u_dS<dim>::ComputeContribution( const InterFace<dim>& iface )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( iface.UsesLocalCoordinates() == true );
    assert( iface.HasInterveningElement() );

    // initialize output vector, for creating the cross couplings between the face nodes with the intervening element
    // (ordering interface-inside-nodes, interface-outside nodes, intervening-element nodes
    const uint32_t n_total_nodes{ iface.Nodes() + iface.InterveningElement()->Nodes() };
    MathOperatorRHS<dim,InterFace>::RHS.resize( n_total_nodes );
    fill( MathOperatorRHS<dim,InterFace>::RHS.begin(), MathOperatorRHS<dim,InterFace>::RHS.end(), 0. );
    
    // computing node-related fractions of the area of the intervening finite element
    assert( iface.FE()->MidSideNodes() == 0U );
    const double area3 = iface.InterveningElement()->Volume() / static_cast<double>( iface.InterveningElement()->Nodes() );
    
    // populating the coupling matrix with the area-weighted transfer coefficients
    // matrix diagonal
    for ( auto i{0U}; i<n_total_nodes; i++ )
      MathOperatorRHS<dim,InterFace>::RHS[i] = transfer_coefficients_[i] * area3;

// TESTING
/*
MathOperatorRHS<dim,InterFace>::RHS[0] = -1.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[1] = -1.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[2] = -1.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[3] = -1.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[4] = -1.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[5] = -1.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[6] = 0.; // 2.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[7] = 0.; // 2.0e-5 * area3;
MathOperatorRHS<dim,InterFace>::RHS[8] = 0.; // 2.0e-5 * area3;
*/

// MONITORING MIN/MAX VALUES OF TRANSFER TERMS
for ( const auto& tc : transfer_coefficients_ ) {
      transfer_min_ = min( transfer_min_, tc );
      transfer_max_ = max( transfer_max_, tc );
  }

//   cout <<"\ninterface "<< iface.Idx() <<": ";
//   for ( const auto& rit : MathOperatorRHS<dim,InterFace>::RHS ) cout << rit <<" ";
//   cout << endl;

} // end ComputeContribution





/** AssignToGlobal matrix function
 
     Ordering of nodes in the vectors is  inner face nodes, outer face nodes, and intervening element nodes
     
    @author SKM - modified from LKT's version that eliminates Dirichlet constraints from G
*/
template<uint32_t dim>
void NumIntegral_dudn_rhsop_u_dS<dim>::AssignToGlobal( const InterFace<dim>& iface,
                                                       vector<double>& rhs,
                                                       const vector<size_t>& DOF_indexes )
	{
    // OK: cerr <<"\nNumIntegral_dudn_rhsop_u_dS<dim>::AssignToGlobal: is getting called"<< endl;
    const uint32_t n_total_nodes{ iface.Nodes() + iface.InterveningElement()->Nodes() },
                   n_elmt_nodes{ iface.InterveningElement()->Nodes() };

		// vectors for mapping local to global node indices for test and basic operands
		vector<size_t> IDT( n_total_nodes ); // ii
    for ( uint32_t n{0U}; n<n_total_nodes; ++n )
      // var1_comp1 var1_comp2 ... var2_comp1 ordering
      for (uint32_t i{0U}; i<iface.Nodes(); ++i )
        IDT[n] = iface.N(n)->Idx();

    // mapping local to global node indices
    for ( auto i{0U}; i < n_elmt_nodes; i++ )
      IDT[i] = iface.N(i,INSIDE)->Idx();
    for ( auto i{0U}; i < n_elmt_nodes; i++ )
      IDT[i+n_elmt_nodes] = iface.N(i,OUTSIDE)->Idx();
    for ( auto i{0U}; i < n_elmt_nodes; i++ )
      IDT[i+n_elmt_nodes+n_elmt_nodes] = iface.InterveningElement()->N(i)->Idx();

    // taking offsets into account
    for (auto i{0U}; i < n_total_nodes; i++) {
         IDT[i] += this->TestOperandOffset();
         IDT[i]  = DOF_indexes[ IDT[i] ]; //-> to the global index
      }

    // perform assignment from local matrix to global matrix
    if (this->multiply_accumulate_) {
      for (auto i{0U}; i < n_total_nodes; i++)
        if (IDT[i] != NULL_IDX)
          rhs[IDT[i]] *= this->RHS[i] * this->factor_;
    }
    else if (this->add_accumulate_ || this->add_accumulate_later_){
        for (auto i{0U}; i < n_total_nodes; i++)
          if (IDT[i] != NULL_IDX)
            rhs[IDT[i]] += this->RHS[i] * this->factor_;
      }
    else if (this->subtract_accumulate_ || this->subtract_accumulate_later_) {
        for (auto i{0U}; i < n_total_nodes; i++)
          if (IDT[i] != NULL_IDX)
            rhs[IDT[i]] -= this->RHS[i] * this->factor_;
      }
    else throw csmp::Exception(ERROR,
                              "NumIntegral_dudn_rhsop_u_dS<dim>::AssignToGlobal(InterFace)",
                              "accumulation instructions could not be parsed.");

	} // end AssignToGlobal (InterFace)


// FOR DEBUGGING
//cout <<"\nNumIntegral_dudn_rhsop_u_dS: on InterFace "<< e.Idx() << endl;
//MathOperatorLHS<dim>::LHS.Out();

// for ( auto i{0U}; i<this->LHS.Rows(); i++ )
//   if ( this->LHS(i,i) < numeric_limits::epsilon() ) 
//     cout <<"\nNumIntegral_dudn_rhsop_u_dS: zero element in diagonal of element matrix.";



/** Analytic one-dimensional temperature profile in homogeneous material.
 
    @param val_farfield difference of variable value in medium before boundary value was changed (oC) and the new value assigned at X0
    @param diffusivity  diffusivity (m2/s)
    @param x distance from interface (m)
    @param t time (s)
*/
template<uint32_t dim>
double NumIntegral_dudn_rhsop_u_dS<dim>::Diffusion1D( double val_farfield, double diffusivity, double x, double t )
 {
    return val_farfield * erf( x / (2. * sqrt( diffusivity * t ) ) );
 }
 

/**
     Derivative of dependent variable at  position x0, as a function of time, t, and the initial difference between the variable value in the farfield
     and at the X0 boundary
*/
template<uint32_t dim>
double NumIntegral_dudn_rhsop_u_dS<dim>::Grad_Var_AtX0( double val_farfield, double diffusivity, double t )
 {
    // derivative pf previous expression
    constexpr double pi = 3.14159265358979323846;
    const double     sqrtPi{ sqrt(pi) }; // x{0.};
    // const double     x2{ x * x };
    // exp(0) = 1! double dTdx  = T0 * exp( -x2 ) / (kappa * t);
    double dvaldx  = val_farfield / (diffusivity * t);
    dvaldx /= sqrtPi * sqrt( diffusivity * t );
 
    // derivative of dependent variable at time t
    const double max_grad(1e20);
    assert( fabs(dvaldx) < max_grad );
    return dvaldx;
 }


/**
     Derivative *  conductivity product at  position x0, as a function of time and the initial difference between the variable value in the farfield
     and at the X0 boundary
*/
template<uint32_t dim>
double NumIntegral_dudn_rhsop_u_dS<dim>::Flux1DAtX0( double val_farfield, double diffusivity, double t, double conductivity )
 {
    // derivative pf previous expression
    constexpr double pi = 3.14159265358979323846;
    const double     sqrtPi{ sqrt(pi) }; // x{0.};
    // const double     x2{ x * x };
    // exp(0) = 1! double dTdx  = T0 * exp( -x2 ) / (kappa * t);
    double dvaldx  = val_farfield / (diffusivity * t);
    dvaldx /= sqrtPi * sqrt( diffusivity * t );
 
    // heat flux
    return dvaldx * conductivity;
 }



template class NumIntegral_dudn_rhsop_u_dS<1U>;
template class NumIntegral_dudn_rhsop_u_dS<2U>;
template class NumIntegral_dudn_rhsop_u_dS<3U>;

} // csmp











