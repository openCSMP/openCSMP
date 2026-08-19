#include "NumIntegral_BT_D_op_dV.h"
#include "PropertyDatabase.h"
#include "mechanics.h"
#include "Face.h"
#include "CSMP_mathUtilities.h"
#include "Element.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_BT_D_op_dV<dim,CELL>::NumIntegral_BT_D_op_dV( const PropertyDatabase<dim>& pref,
                                                          const char* oper,
                                                          const char* youngs,
                                                          const char* poissons,
                                                          const char* test,
                                                          bool        plane_strain )
  : MathOperatorRHS<dim,CELL>( pref, oper, test ),
    E_key_( pref.StorageKey( youngs ) ),
    nu_key_( pref.StorageKey( poissons ) ),
    B_( 3U, dim ),
    D_( (dim-1U)*3U, (dim-1U)*3U ),
    E_( 1U ),
    nu_( 1U ),
    plane_strain_( plane_strain )
 {
    MathOperatorRHS<dim,CELL>::Name( "NumIntegral_BT_D_op_dV", oper, test );

    if ( E_key_.place != nu_key_.place )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim>::(constructor)",
                             "Young's modulus and Poisson's ratio must have the same placement." );

    if ( (E_key_.place != ELEMENT && E_key_.place != ELEMENT_INTEGRATION_POINT) || E_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim>::(constructor)",
                             "Young's modulus must be a scalar placed on the element or integration point." );

    if ( (nu_key_.place != ELEMENT && nu_key_.place != ELEMENT_INTEGRATION_POINT) || nu_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim>::(constructor)",
                             "Poisson's ratio must be a scalar placed on the element or integration point." );

 } // end constructor



template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_op_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    assert( e.FE()->Isoparametric() == true );

    // Young's modulus and Poisson's ratio
    if ( E_key_.place == ELEMENT ) {
		     // reading Young's modulus (must be an element variables)
		     E_[0U] = e.Read( E_key_ );
		     // getting Poisson's ratio
	       nu_[0U] = e.Read( nu_key_ );
      }
    else if ( E_key_.place == ELEMENT_INTEGRATION_POINT ) { // Integration Point
         E_.resize(e.IntegrationPoints());
         nu_.resize(e.IntegrationPoints());
         for ( uint32_t i{0U}; i<e.IntegrationPoints(); ++i ) {
              E_[i]  = e.Read( i, E_key_ );
              nu_[i] = e.Read( i, nu_key_ );
           }
      }
    else if ( E_key_.place == NODE ) { // Integration Point
         E_.resize(e.IntegrationPoints());
         nu_.resize(e.IntegrationPoints());
         for ( uint32_t i{0U}; i<e.IntegrationPoints(); ++i ) {
              E_[i]  = e.PropertyValueAtIntegrationPoint( E_key_, i );
              nu_[i] = e.PropertyValueAtIntegrationPoint( nu_key_, i );
           }
      }
    else throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim,CELL>::GetOperands",
                               "Young's modulus and Poisson's ratio must be placed on ELEMENT, its integration points or NODE");

    // material Operand 'dilation' or other
    const auto placement = MathOperatorRHS<dim,CELL>::MaterialOperandPlacement();
    const auto op_key = MathOperatorRHS<dim,CELL>::MaterialOperandKey();
    
    const uint32_t n_ip = e.IntegrationPoints();
    initial_strains_.assign(n_ip, {0.0, 0.0, 0.0}); // Stack-friendly storage

    if (placement == ELEMENT) {
        // Read once for the whole element
        if (MathOperatorRHS<dim,CELL>::MaterialOperandType() == SCALAR) {
            e.Read(op_key, sc_);
            for (uint32_t i = 0; i < dim; ++i) initial_strains_[0][i] = sc_();
        } else if (MathOperatorRHS<dim,CELL>::MaterialOperandType() == VECTOR) {
            e.Read(op_key, vc_);
            for (uint32_t i = 0; i < dim; ++i) initial_strains_[0][i] = vc_[i];
        }
        // If it's a full tensor, you'd extract the diagonal components here
    } 
    else { // ELEMENT_INTEGRATION_POINT or NODE
        for (uint32_t ip = 0; ip < n_ip; ++ip) {
            // Assuming the base class or element has a method to get the vector at the IP
            // If it returns a scalar, assign it to all 'dim' components
            double val = ( op_key.place == ELEMENT_INTEGRATION_POINT ) ? e.Read(ip, op_key) :
                           e.PropertyValueAtIntegrationPoint( op_key, ip );
            for (uint32_t i = 0; i < dim; ++i) initial_strains_[ip][i] = val;
        }
    }
}


template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_op_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    const uint32_t num_nodes = e.Nodes();
    const uint32_t num_dofs  = dim * num_nodes;
    const uint32_t voigt_dim = (dim == 1U) ? 1U : (dim == 2U) ? 3U : 6U;

    auto& RHS = MathOperatorRHS<dim,CELL>::RHS;
    RHS.assign(num_dofs, 0.0);

    const bool is_elem_prop   = (E_key_.place == ELEMENT);
    const bool is_elem_strain = (this->MaterialOperandPlacement() == ELEMENT);

    // 1. Pre-compute stiffness matrix D if constant
    if (is_elem_prop) {
        if constexpr (dim == 2U) {
            plane_strain_ ? planeStrainMatrix(E_[0U], nu_[0U], D_) : planeStressMatrix(E_[0U], nu_[0U], D_);
        } else if constexpr (dim == 3U) {
            stiffnessMatrix(E_[0U], nu_[0U], D_);
        }
    }

    // Accumulation lambda using direct stack arrays (no heap allocs, no matrix transpositions)
    // ----------------------------------------------------------------------------------------
    auto accumulate_point = [&](uint32_t ip, double w_detJ) {
        if (!is_elem_prop) {
            if constexpr (dim == 2U) {
                plane_strain_ ? planeStrainMatrix(E_[ip], nu_[ip], D_) : planeStressMatrix(E_[ip], nu_[ip], D_);
            } else if constexpr (dim == 3U) {
                stiffnessMatrix(E_[ip], nu_[ip], D_);
            }
        }

        const uint32_t strain_idx = is_elem_strain ? 0U : ip;

        // Compute local stress vector: sigma_0 = D * eps_0
        // We only loop to 'dim' because shear initial strains are zero
        std::array<double, 6> sigma_0 = {0.0};
        for (uint32_t row = 0; row < voigt_dim; ++row) {
            for (uint32_t col = 0; col < dim; ++col) {
                sigma_0[row] += D_(row, col) * initial_strains_[strain_idx][col];
            }
        }

        // Compute B^T * sigma_0 directly without building BT
        // F_n = sum_j ( B_jn * sigma_{0,j} )
        for (uint32_t n = 0; n < num_dofs; ++n) {
            double force_n = 0.0;
            for (uint32_t j = 0; j < voigt_dim; ++j) {
                force_n += B_(j, n) * sigma_0[j];
            }
            RHS[n] += force_n * w_detJ;
        }
    }; // end lambda

    // 2. Simplex shortcut (Barycenter)
    if (e.FE()->IsSimplex() && e.Interpolation() == 1U) {
        const double detJ = e.dN_AtBaryCenter(B_, dim);
        if (detJ <= 0.0)
          throw csmp::Exception(FATAL_ERROR, "NumIntegral_BT_D_op_dV<dim,CELL>::ComputeContribution", "Jacobian failed.");
        
        accumulate_point(0U, e.Volume());
        return;
    }

    // 3. Standard Gauss quadrature
    const uint32_t n_ip = e.IntegrationPoints();
    for (uint32_t i = 0; i < n_ip; ++i) {
        const double detJ = e.dN_AtIntegrationPoint(B_, i, dim);
        if (detJ <= 0.0)
          throw csmp::Exception(FATAL_ERROR, "NumIntegral_BT_D_op_dV<dim,CELL>::ComputeContribution", "Jacobian failed.");
        
        accumulate_point(i, e.WeightAtIntegrationPoint(i) * detJ);
    }
}


template class NumIntegral_BT_D_op_dV<1U,Element>;
template class NumIntegral_BT_D_op_dV<2U,Element>;
template class NumIntegral_BT_D_op_dV<3U,Element>;

template class NumIntegral_BT_D_op_dV<1U,Face>;
template class NumIntegral_BT_D_op_dV<2U,Face>;
template class NumIntegral_BT_D_op_dV<3U,Face>;

} // namespace csmp



/* COMPLICATED EARLIER VERSION */

/*
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_op_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    assert( e.FE()->Isoparametric() == true );

    // Young's modulus and Poisson's ratio
    if ( E_key_.place == ELEMENT ) {
		     // reading Young's modulus (must be an element variables)
		     E_[0U] = e.Read( E_key_ );
		     // getting Poisson's ratio
	       nu_[0U] = e.Read( nu_key_ );
      }
    else if ( E_key_.place == ELEMENT_INTEGRATION_POINT ) { // Integration Point
         E_.resize(e.IntegrationPoints());
         nu_.resize(e.IntegrationPoints());
         for ( uint32_t i{0U}; i<e.IntegrationPoints(); ++i ) {
              E_[i]  = e.PropertyValueAtIntegrationPoint( E_key_, i );
              nu_[i] = e.PropertyValueAtIntegrationPoint( nu_key_, i );
           }
      }
    else throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim,CELL>::GetOperands",
                               "Young's modulus and Poisson's ratio must be placed on ELEMENT or its integration points");

    // input strain / dilatation
    switch ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() )
      {
        case ELEMENT:
          switch ( MathOperatorRHS<dim,CELL>::MaterialOperandType() )
            {
              case SCALAR:
                e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), sc_ );
                MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonalAndZeroOffDiagonal( dim, sc_() );
                break;
              case VECTOR:
                e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), vc_ );
                MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonal( vc_ );
                break;
              case TENSOR:
                e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), ts_ );
                MathOperatorRHS<dim,CELL>::MTRL[0] = ts_;
                break;
              default:
                throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_op_dV::GetOperands",
                                       "Unsupported MaterialOperandType for ELEMENT placement." );
            }
          break;

        case ELEMENT_INTEGRATION_POINT:
          {
            const auto n_ip{ e.IntegrationPoints() };
            MathOperatorRHS<dim,CELL>::MTRL.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              for ( uint32_t j{0U}; j < dim; ++j )
                MathOperatorRHS<dim,CELL>::MTRL[i](j,j) = e.Read( i, MathOperatorRHS<dim,CELL>::MaterialOperandKey() );
          }
          break;

        default: // NODE
          {
            const auto n_ip{ e.IntegrationPoints() };
            MathOperatorRHS<dim,CELL>::MTRL.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              MathOperatorRHS<dim,CELL>::PropertyAtIntegrationPoint( e,
                MathOperatorRHS<dim,CELL>::MaterialOperandKey(), i,
                MathOperatorRHS<dim,CELL>::MTRL[i] );
          }
          break;
      }

 } // end GetOperands





template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_op_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    const uint32_t num_nodes = e.Nodes();
    const uint32_t num_dofs  = dim * num_nodes;
    
    // Voigt notation size: 1D -> 1, 2D -> 3 (exx, eyy, gxy), 3D -> 6 (exx, eyy, ezz, gyz, gzx, gxy)
    constexpr uint32_t voigt_dim = (dim == 1U) ? 1U : (dim == 2U) ? 3U : 6U;

    // Initialize RHS vector
    auto& RHS = MathOperatorRHS<dim,CELL>::RHS;
    RHS.assign(num_dofs, 0.0);

    const bool is_elem_prop   = (E_key_.place == ELEMENT);
    const bool is_elem_strain = (this->MaterialOperandPlacement() == ELEMENT);

    // 1. Pre-compute element-constant material stiffness matrix D
    if (is_elem_prop) {
        if constexpr (dim == 2U) {
            if (plane_strain_) planeStrainMatrix(E_[0U], nu_[0U], D);
            else               planeStressMatrix(E_[0U], nu_[0U], D);
        } 
        else if constexpr (dim == 3U) {
            stiffnessMatrix(E_[0U], nu_[0U], D);
        }
    }

    // 2. Pre-fill strain vector if uniform over element
    STR.Resize(voigt_dim, 1U);
    STR.Zero();
    if (is_elem_strain) {
        for (uint32_t i{0U}; i < dim; ++i) {
            STR(i, 0U) = MathOperatorRHS<dim,CELL>::MTRL[0U](i, i);
        }
    }

    // Lambda to evaluate B^T * D * strain at a single evaluation point
    auto accumulate_point = [&](uint32_t ip, double weight_detJ) {
        // Rebuild D if properties vary spatially
        if (!is_elem_prop) {
            if constexpr (dim == 2U) {
                if (plane_strain_) planeStrainMatrix(E_[ip], nu_[ip], D);
                else               planeStressMatrix(E_[ip], nu_[ip], D);
            } 
            else if constexpr (dim == 3U) {
                stiffnessMatrix(E_[ip], nu_[ip], D);
            }
        }

        // Update strain vector if operand varies spatially
        if (!is_elem_strain) {
            for (uint32_t j{0U}; j < dim; ++j) {
                STR(j, 0U) = MathOperatorRHS<dim,CELL>::MTRL[ip](j, j);
            }
        }

        // Compute equivalent initial stress: {sigma_0} = [D] * {e_0}
        TEMP = D;
        TEMP *= STR.Data();

        // Map to nodal forces: {F_e} += [B]^T * {sigma_0} * w_detJ
        B.Transposed(BT);
        BT *= TEMP;

        for (uint32_t n{0U}; n < num_dofs; ++n) {
            RHS[n] += BT(n, 0U) * weight_detJ;
        }
    };

    // 3. Simplex shortcut (constant B matrix evaluated at barycenter)
    if (e.FE()->IsSimplex() && e.Interpolation() == 1U) {
        const double detJ = e.dN_AtBaryCenter(B, dim);
        if (detJ <= 0.0) {
            throw csmp::Exception(FATAL_ERROR, "NumIntegral_BT_D_op_dV::ComputeContribution",
                                  "Jacobian transformation failed at barycenter.");
        }
        accumulate_point(0U, e.Volume());
        return;
    }

    // 4. Standard Gauss quadrature loop
    const uint32_t n_ip = e.IntegrationPoints();
    for (uint32_t i{0U}; i < n_ip; ++i) {
        const double detJ = e.dN_AtIntegrationPoint(B, i, dim);
        if (detJ <= 0.0) {
            throw csmp::Exception(FATAL_ERROR, "NumIntegral_BT_D_op_dV::ComputeContribution",
                                  "Jacobian transformation failed at integration point.");
        }
        
        const double w_detJ = e.WeightAtIntegrationPoint(i) * detJ;
        accumulate_point(i, w_detJ);
    }
}
*/


/**
Reads the operand variable (strain/dilatation), Young's modulus and
Poisson's ratio. Only the x, y, z diagonal components of the strain
will be used — shear strains are not considered.

The righthand contribution:  {F} = ElementIntegral [BT][D]{e} dV

With the strain source term: {e} = strain * {1,1,0} (2D)

is accumulated into the righthand vector.
*/
/*
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_op_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    assert( e.FE()->Isoparametric() == true );

    // Young's modulus and Poisson's ratio
    if ( E_key_.place == ELEMENT ) {
		     // reading Young's modulus (must be an element variables)
		     E_[0U] = e.Read( E_key_ );
		     // getting Poisson's ratio
	       nu_[0U] = e.Read( nu_key_ );
      }
    else if ( E_key_.place == ELEMENT_INTEGRATION_POINT ) { // Integration Point
         E_.resize(e.IntegrationPoints());
         nu_.resize(e.IntegrationPoints());
         for ( uint32_t i{0U}; i<e.IntegrationPoints(); ++i ) {
              E_[i]  = e.PropertyValueAtIntegrationPoint( E_key_, i );
              nu_[i] = e.PropertyValueAtIntegrationPoint( nu_key_, i );
           }
      }
    else throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim,CELL>::GetOperands",
                               "Young's modulus and Poisson's ratio must be placed on ELEMENT or its integration points");

    // input strain / dilatation
    switch ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() )
      {
        case ELEMENT:
          switch ( MathOperatorRHS<dim,CELL>::MaterialOperandType() )
            {
              case SCALAR:
                e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), sc_ );
                MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonalAndZeroOffDiagonal( dim, sc_() );
                break;
              case VECTOR:
                e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), vc_ );
                MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonal( vc_ );
                break;
              case TENSOR:
                e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), ts_ );
                MathOperatorRHS<dim,CELL>::MTRL[0] = ts_;
                break;
              default:
                throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_op_dV::GetOperands",
                                       "Unsupported MaterialOperandType for ELEMENT placement." );
            }
          break;

        case ELEMENT_INTEGRATION_POINT:
          {
            const auto n_ip{ e.IntegrationPoints() };
            MathOperatorRHS<dim,CELL>::MTRL.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              for ( uint32_t j{0U}; j < dim; ++j )
                MathOperatorRHS<dim,CELL>::MTRL[i](j,j) = e.Read( i, MathOperatorRHS<dim,CELL>::MaterialOperandKey() );
          }
          break;

        default: // NODE
          {
            const auto n_ip{ e.IntegrationPoints() };
            MathOperatorRHS<dim,CELL>::MTRL.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              MathOperatorRHS<dim,CELL>::PropertyAtIntegrationPoint( e,
                MathOperatorRHS<dim,CELL>::MaterialOperandKey(), i,
                MathOperatorRHS<dim,CELL>::MTRL[i] );
          }
          break;
      }

 } // end GetOperands
*/

/**
Computes the spatial integral over the matrix product:

[B]^T [D] {e} detJ

as the element righthand contribution, accumulated into the base
class protected member vector {RHS}.

Application: linear elasticity computations.
*/
/*
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_op_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // detect simplex elements where B is constant
    const bool is_simplex( e.FE()->IsSimplex() && e.Interpolation() == 1 );

    // build D for element-constant properties
    if ( E_key_.place == ELEMENT )
      {
        if constexpr ( dim == 2U )
          {
            if ( plane_strain_ ) planeStrainMatrix( E_[0U], nu_[0U], D );
            else                 planeStressMatrix( E_[0U], nu_[0U], D );
          }
        else if constexpr ( dim == 3U ) stiffnessMatrix( E_[0U], nu_[0U], D );
      }

    // initialise output vector
    MathOperatorRHS<dim,CELL>::RHS.resize( dim * e.Nodes() );
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );

    // map isostatic diagonal components of strain matrix MTRL into STR column vector
    STR.Resize(dim,1);
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
      for ( uint32_t i{0U}; i < dim; ++i )
        STR(i,0U) = MathOperatorRHS<dim,CELL>::MTRL[0](i,i);

    // for simplex elements B is constant — evaluate dN once at the barycentre
    if ( is_simplex )
      {
        const double detJ = e.dN_AtBaryCenter( B, dim );

        if ( detJ <= 0.0 )
          {
            cerr <<"\n\tElement "<< e.Idx() <<": determinant of Jacobian at barycentre: "<< detJ << endl;
            throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_op_dV<dim>::ComputeContribution",
                                   "Jacobian transformation failed. Element nodes are perhaps not numbered correctly." );
          }

        B.Transposed( BT );

        TEMP  = D;
        TEMP *= STR.Data();

        // BT(dim * nodes x 3) * TEMP(3x1) -> BT(dim * nodes x 1)
        BT *= TEMP;
        BT *= e.Volume();

        for ( uint32_t n{0U}; n < BT.Rows(); ++n )
          MathOperatorRHS<dim,CELL>::RHS[n] += BT(n,0U);

        return;
      }

    // numerical integration over Gauss points
    for ( uint32_t i{0U}; i < e.IntegrationPoints(); ++i )
      {
        // rebuild D at each integration point if properties vary spatially
        if ( E_key_.place == ELEMENT_INTEGRATION_POINT || E_key_.place == NODE )
          {
            if constexpr ( dim == 2U )
              {
                if ( plane_strain_ ) planeStrainMatrix( E_[i], nu_[i], D );
                else                 planeStressMatrix( E_[i], nu_[i], D );
              }
            else if constexpr ( dim == 3U ) stiffnessMatrix( E_[i], nu_[i], D );
          }
      
        // dim argument assures that the matrix gets transformed for a solution variable that is a vector
        const double detJ = e.dN_AtIntegrationPoint( B, i, dim );

        if ( detJ <= 0.0 )
          {
            cerr <<"\n\tElement "<< e.Idx() <<": determinant of Jacobian at Gauss point "<< i
                 <<": "<< detJ << endl;
            throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_op_dV<dim>::ComputeContribution",
                                   "Jacobian transformation failed. Element nodes are perhaps not numbered correctly." );
          }

        B.Transposed( BT );

        // update STR for spatially varying operands
        if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
             MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE )
          {
            STR.Resize(dim,1);
            for ( uint32_t j{0U}; j < dim; ++j )
              STR(j,0U) = MathOperatorRHS<dim,CELL>::MTRL[i](j,j);
          }

        // D(3x3) * STR(3x1) -> TEMP(3x1)
        TEMP  = D;
        TEMP *= STR.Data();

        // BT(dim*nodes x 3) * TEMP(3x1) -> BT(dim*nodes x 1)
        BT *= TEMP;

        const double w = e.WeightAtIntegrationPoint(i) * detJ;
        for ( uint32_t n{0U}; n < BT.Rows(); ++n )
          {
            BT(n,0U) *= w;
            MathOperatorRHS<dim,CELL>::RHS[n] += BT(n,0U);
          }
      }

 } // end ComputeContribution

*/
