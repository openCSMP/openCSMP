#include "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region.h"
#include "Index.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
    Constructor gets the finite element types that  NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region  needs to turn the InterFace into the equivalent
    of a volumetric element from the FiniteElementManager. Since InterFaces can either be triangles or quadrilaterals, their extruded volumetric
    pendants are prism and hexahedral elements.
    
    @attention the variable 'thickness' needs to be defined on the InterFace, else the pde operator will not work.
*/
template<uint32_t dim>
NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region( const Model<dim>& model,
                                                                                                           const char*       oper,
                                                                                                           const char*       basic,
                                                                                                           const char*       test,
                                                                                                           double dt )
  : MathOperatorLHS<dim,InterFace>(model.Database(),oper,basic,test),
    diffusivity_key_(model.Database().StorageKey("hydraulic diffusivity")),
    delta_t_(dt)
{
    MathOperatorLHS<dim,InterFace>::Name("NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region", oper, basic, test );
    
    // element interpolation function operand i
    if ( MathOperatorLHS<dim,InterFace>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,InterFace>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::(constructor)",
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    // element interpolation function operand j
    if ( MathOperatorLHS<dim,InterFace>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,InterFace>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::(constructor)",
                      test, "Operand (test) must be a scalar property placed on the nodes." );

    // oper = material operand
    if ( MathOperatorLHS<dim,InterFace>::MaterialOperandPlacement() != ELEMENT ||
         MathOperatorLHS<dim,InterFace>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim,CELL>::GetOperands",
                                  "material operand placements other than on the CELL are not handled yet. " );
}




/**
   Computes difference between the variable values at the topologically collocated nodes, including that in the intervening
   lower-dimensional region, using it to compute  the transfer coefficient (=coupling coefficient) at the time level t + delta t.
   
   The hydraulic conductivity (Operand value) is read from the higher-dimensional neighbors of the lower-dimensional region.
*/
template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::GetOperands( const InterFace<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.HasInterveningElement() );
    
    // retrieving conductivity from inner and outer higher-dimensional parent elements
    // (their values for the intervening element are ignored as they are accumulated separately)
    const double conductivity_inner = e.InnerParent()->Read( MathOperatorLHS<dim,InterFace>::MaterialOperandKey() );
    const double diffusivity_inner  = e.InnerParent()->Read( diffusivity_key_ );
    const double conductivity_outer = e.OuterParent()->Read( MathOperatorLHS<dim,InterFace>::MaterialOperandKey() );
    const double diffusivity_outer  = e.OuterParent()->Read( diffusivity_key_ );

    // (the area of the interface is taken into account later)
    const auto n_nodes{ e.InterveningElement()->Nodes() };
    transfer_coefficients_.resize( n_nodes * 3U );


    // 1. First-order finite-difference approximation of the interface fluxes
    // ----------------------------------------------------------------------
    if ( order1_FD_approximation_ ) {
         double     val_intervening = e.PropertyValueAtBaryCenter( MathOperatorLHS<dim,InterFace>::TestOperandKey() );
         Point<dim> unrml( e.InterveningElement()->UnitNormal() );
         
         // gradient of test function-variable is computed from the facet-normal difference
         // -------------------------------------------------------------------------------
         const double val_inner = e.InnerParent()->PropertyValueAtBaryCenter( MathOperatorLHS<dim,InterFace>::TestOperandKey() );
         // separation of barycenters of inner and intervening cells
         Point<dim> bctr_distance = e.InnerParent()->BaryCenter() - e.InterveningElement()->BaryCenter();
         // projecting this distance on the normal of the intervening element
         const double distance_to_inner_bctr = fabs( dotProduct( unrml, bctr_distance ) );
         // the same gradient is used for all nodes
         double grad_var0 = (val_inner - val_intervening) / distance_to_inner_bctr;
         for ( auto i{0U}; i<n_nodes; i++ )
           // source term due to pressure difference
           transfer_coefficients_[i] = -conductivity_inner * grad_var0;
           
         const double val_outer = e.OuterParent()->PropertyValueAtBaryCenter( MathOperatorLHS<dim,InterFace>::TestOperandKey() );
         // separation of barycenters of outer and intervening cells
         bctr_distance = e.OuterParent()->BaryCenter() - e.InterveningElement()->BaryCenter();
         // projecting this distance on the normal of the intervening element
         const double distance_to_outer_bctr = fabs( dotProduct( unrml, bctr_distance ) );
         // the same gradient is used for all nodes
         grad_var0 = (val_outer - val_intervening) / distance_to_outer_bctr;
         for ( auto i{0U}; i<n_nodes; i++ )
           // source term due to pressure difference
           transfer_coefficients_[i+n_nodes] = -conductivity_outer * grad_var0;

         // transfer into the intervening element
         for ( auto i{0U}; i<n_nodes; i++ )
           // source term due to pressure difference (negative because incoming)
           transfer_coefficients_[i + n_nodes + n_nodes] = -(transfer_coefficients_[i] + transfer_coefficients_[i+n_nodes]);

         return;
      }


    // and the values of the transported variable at the interface faces and the intervening element
/*
    vector<ScalarVariable> node_var_iface, node_var_elmt;
    e.NodePropertyVector( MathOperatorLHS<dim,InterFace>::TestOperandKey(), node_var_iface );
    e.InterveningElement()->NodePropertyVector( MathOperatorLHS<dim,InterFace>::TestOperandKey(), node_var_elmt );
*/
    // computing transfer coefficient values using the variable values at the interface nodes
    // and the intervening element (inside - middle = 3, outside - middle = 3, 6 in total.
    // --------------------------------------------------------------------------------------

    // inside face and intervening element (nodes have the same indices)
    for ( auto i{0U}; i<n_nodes; i++ ) {
         // nodes of inner face come first
         double delta_var0 = e.N(i,INSIDE)->Read( MathOperatorLHS<dim,InterFace>::TestOperandKey() ) -
                             e.InterveningElement()->N(i)->Read( MathOperatorLHS<dim,InterFace>::TestOperandKey() );
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
         double delta_var0 = e.MatchingN(i,OUTSIDE)->Read( MathOperatorLHS<dim,InterFace>::TestOperandKey() ) -
                             e.InterveningElement()->N(i)->Read( MathOperatorLHS<dim,InterFace>::TestOperandKey() );
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

 } // end GetOperands









/** Computes integral over shape function derivatives squared. These are the shape functions of the volumetric elements that match
    the lower-dimensional InterFace.

    @attention The connections could also be established using line elements, but this would require a more complicated
    accumulation process.
    
    @todo in GetOperands implement the collection of material operands from integration points otherwise second part of method will not work
*/
template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::ComputeContribution( const InterFace<dim>& iface )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( iface.UsesLocalCoordinates() == true );
    assert( iface.HasInterveningElement() );

    // initialize output matrix, for creating the cross couplings between the face nodes with the intervening element
    // (ordering interface-inside-nodes, interface-outside nodes, intervening-element nodes
    const uint32_t n_total_nodes{ iface.Nodes() + iface.InterveningElement()->Nodes() },
                   n_elmt_nodes{ iface.InterveningElement()->Nodes() };
    MathOperatorLHS<dim,InterFace>::LHS.Resize( n_total_nodes, n_total_nodes );
    MathOperatorLHS<dim,InterFace>::LHS.Zero();
    
    // computing node-related fractions of the area of the intervening finite element
    assert( iface.FE()->MidSideNodes() == 0U );
    const double area3 = iface.InterveningElement()->Volume() / static_cast<double>( iface.InterveningElement()->Nodes() );
    
    // weighting the transfer coefficients by element area corresponding to node
    for ( auto& tc : transfer_coefficients_ ) tc *= area3;

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
         MathOperatorLHS<dim,InterFace>::LHS(n_elmt_nodes*3U-1-i,i+n_elmt_nodes) = -transfer_coefficients_[i+n_elmt_nodes];
         // lower diagonal
         MathOperatorLHS<dim,InterFace>::LHS(i+n_elmt_nodes,n_elmt_nodes*3U-1-i) = -transfer_coefficients_[i+n_elmt_nodes];
      }

// MONITORING MIN/MAX VALUES OF TRANSFER TERMS
for ( const auto& tc : transfer_coefficients_ ) {
      transfer_min_ = min( transfer_min_, tc );
      transfer_max_ = max( transfer_max_, tc );
  }

//   MathOperatorLHS<dim,InterFace>::LHS.Out();
//   cout << endl;

} // end ComputeContribution





/** AssignToGlobal matrix function
 
     Ordering of nodes in the vectors is  inner face nodes, outer face nodes, and intervening element nodes
     
    @author SKM - modified from LKT's version that eliminates Dirichlet constraints from G
*/
template<uint32_t dim>
void NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::AssignToGlobal( const InterFace<dim>& iface,
                                                                             SparseMatrix& G,
                                                                             vector<double>& pivotVector,
                                                                             const vector<size_t>& DOF_indexes )
	{
    // OK: cerr <<"\nNumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::AssignToGlobal: is getting called"<< endl;
    const uint32_t n_total_nodes{ iface.Nodes() + iface.InterveningElement()->Nodes() },
                   n_elmt_nodes{ iface.InterveningElement()->Nodes() };

		// map local to global indexes for test and basic operands
		this->IDT.resize(n_total_nodes); // ii
		this->IDB.resize(n_total_nodes); // jj
    
    const uint32_t n_elmt_nodes2{ n_elmt_nodes * 2 };
    
    // interface nodes + intervening element
		for ( auto i{0U}; i < n_elmt_nodes; i++ ) {
         this->IDT[i] = iface.N(i,INSIDE)->Idx();
         this->IDB[i] = this->IDT[i];
         this->IDT[i+n_elmt_nodes] = iface.N(i,OUTSIDE)->Idx();
         this->IDB[i+n_elmt_nodes] = this->IDT[i+n_elmt_nodes];
         this->IDT[i+n_elmt_nodes2] = iface.InterveningElement()->N(i)->Idx();
         this->IDB[i+n_elmt_nodes2] = this->IDT[i+n_elmt_nodes2];
      }

    // ONLY SCALARS are used in the pde operator
    assert( this->BasicOperandType() == SCALAR );
    assert( this->TestOperandType() == SCALAR );
    
    // applying offset to interpolation function operand
		for (auto i{0U}; i < this->IDT.size(); i++) {
        this->IDT[i] += this->TestOperandOffset();
        this->IDT[i] = DOF_indexes[ this->IDT[i] ];
      }

    // applying offset to weighting function operand
		for (auto i{0U}; i < this->IDB.size(); i++) {
        this->IDB[i] += this->BasicOperandOffset();
        this->IDB[i] = DOF_indexes[ this->IDB[i] ];
      }

		// get local value from the given (local) node
		vector<double> nodal_values(this->IDB.size(),1.);
//  /* causes errratic values
      {
        // interface and intervening element
        for ( auto nIdx{0U}; nIdx < n_elmt_nodes; ++nIdx ) {
             nodal_values[nIdx] = iface.N(nIdx,INSIDE)->Read(this->TestOperandKey());
             nodal_values[nIdx+n_elmt_nodes] = iface.N(nIdx,OUTSIDE)->Read(this->TestOperandKey());
             nodal_values[nIdx+n_elmt_nodes2] = iface.InterveningElement()->N(nIdx)->Read(this->TestOperandKey());
          }
      }
//  */
		// perform assignment from local matrix to global matrix


		if ( this->multiply_accumulate_ )
      {
        throw csmp::Exception( ERROR,
                              "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::AssignToGlobal(InterFace):",
                              "Multiply Accumulate not supported yet");

        //if (IDT[i] != NULL_IDX) {
        //  for (size_t j{0U}; j<LHS.Cols(); j++)
        //    if (IDB[j] != NULL_IDX) {
        //      G.MultiplyEntryWith(IDT[i],
        //        IDB[j],
        //        LHS(i, j) * factor_);
        //    }
        //}
      }
    /*
        Luat's factor_ = 1, always in the applications throughout code
    */
		else if ( this->add_accumulate_ || this->add_accumulate_later_ )
		{
			for (auto i{0U}; i < this->LHS.Rows(); i++) {
				if ( this->IDT[i] != NULL_IDX ) {
					for (auto j{0U}; j < this->LHS.Cols(); j++) {
						if ( this->IDB[j] == NULL_IDX ) {
							pivotVector[ this->IDT[i] ] -= this->LHS(i,j) * nodal_values[j] * this->factor_;  // notice sign iface.N(j / this->TestOperandOffset())->Read(this->TestOperand())
						}
						else {
							G.Add( this->IDT[i],
								     this->IDB[j],
								     this->LHS(i, j) * this->factor_ );
						}
					}
				}
			}
		}
		else if ( this->subtract_accumulate_ || this->subtract_accumulate_later_ )
      {
        for ( auto i{0U}; i < this->LHS.Rows(); i++ ) {
          if ( this->IDT[i] != NULL_IDX ) {
            for (auto j{0U}; j < this->LHS.Cols(); j++ ) {
              if ( this->IDB[j] == NULL_IDX ) {
                pivotVector[ this->IDT[i] ] += this->LHS(i,j) * nodal_values[j] * this->factor_;   // notice the sign LHS(i, j) * e.N(j)->Read(TestOperandKey());
              }
              else {
                G.Add( this->IDT[i],
                       this->IDB[j],
                       -this->LHS(i, j) * this->factor_ );
              }
            }
          }
        }
		}
		else
			throw csmp::Exception( ERROR,
				                    "NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::AssignToGlobal(InterFace):",
				                    "accumulation instructions could not be parsed.");

	} // end AssignToGlobal (InterFace)


// FOR DEBUGGING
//cout <<"\nNumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region: on InterFace "<< e.Idx() << endl;
//MathOperatorLHS<dim>::LHS.Out();

// for ( auto i{0U}; i<this->LHS.Rows(); i++ )
//   if ( this->LHS(i,i) < numeric_limits::epsilon() ) 
//     cout <<"\nNumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region: zero element in diagonal of element matrix.";



/** Analytic one-dimensional temperature profile in homogeneous material.
 
    @param val_farfield difference of variable value in medium before boundary value was changed (oC) and the new value assigned at X0
    @param diffusivity  diffusivity (m2/s)
    @param x distance from interface (m)
    @param t time (s)
*/
template<uint32_t dim>
double NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::Diffusion1D( double val_farfield, double diffusivity, double x, double t )
 {
    return val_farfield * erf( x / (2. * sqrt( diffusivity * t ) ) );
 }
 

/**
     Derivative of dependent variable at  position x0, as a function of time, t, and the initial difference between the variable value in the farfield
     and at the X0 boundary
*/
template<uint32_t dim>
double NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::Grad_Var_AtX0( double val_farfield, double diffusivity, double t )
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
double NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>::Flux1DAtX0( double val_farfield, double diffusivity, double t, double conductivity )
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




template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<1U>;
template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<2U>;
template class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<3U>;

} // csmp











