/*
 *  MohrCoulombFailure_Visitor.cpp
 *  Revised and corrected.
 */

#include "MohrCoulombFailure_Visitor.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
#include "Node.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "StressInvariants.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

// delete For2D for 3D and 1D
template<>
MohrCoulombFailure_Visitor<3U>
MohrCoulombFailure_Visitor<3U>::For2D( Model<3U>&, PlaneAssumption, PLACEMENT ) = delete;

template<>
MohrCoulombFailure_Visitor<1U>
MohrCoulombFailure_Visitor<1U>::For2D( Model<1U>&, PlaneAssumption, PLACEMENT ) = delete;

// provide For2D only for 2D
template<>
MohrCoulombFailure_Visitor<2U>
MohrCoulombFailure_Visitor<2U>::For2D( Model<2U>&       model,
                                       PlaneAssumption  assumption,
                                       PLACEMENT        target_placement )
{
    return MohrCoulombFailure_Visitor<2U>( model, assumption, target_placement );
}

// ----------------------------------------------------------------------------
//  Shared constructor — 3D and 1D (no PlaneAssumption)
// ----------------------------------------------------------------------------

template<uint32_t dim>
MohrCoulombFailure_Visitor<dim>::MohrCoulombFailure_Visitor( Model<dim>& model, PLACEMENT analysis_var_placement )
    : Visitor<dim>( MODEL, ELEMENT ),
      plane_assumption_( PlaneAssumption::PLANE_STRESS ),  // irrelevant for 3D/1D
      analysis_var_placement_( analysis_var_placement                    ),
      stress_key_      ( model.Database().StorageKey( "stress"           ) ),
      cohesion_key_    ( model.Database().StorageKey( "cohesion"         ) ),
      friction_key_    ( model.Database().StorageKey( "friction angle"   ) ),
      tensile_str_key_ ( model.Database().StorageKey( "tensile strength" ) ),
      biot_alpha_key_  ( model.Database().StorageKey( "Biot alpha"       ) ),
      nu_key_          ( csmp::Index()                                      ),
      failure_key_     ( model.Database().StorageKey( "failure"          ) ),
      failure01_key_   ( model.Database().StorageKey( "tensile failure"  ) ),
      pf_key_          ( model.Database().IsDefined( "fluid pressure" )
                         ? model.Database().StorageKey( "fluid pressure" )
                         : csmp::Index()                                   )
{
    const char* caller = "MohrCoulombFailure_Visitor (constructor)";

    if ( analysis_var_placement_ != ELEMENT &&
         analysis_var_placement_ != ELEMENT_INTEGRATION_POINT )
        throw csmp::Exception( ERROR, caller,
                               "target placement must be ELEMENT or "
                               "ELEMENT_INTEGRATION_POINT" );

    if ( stress_key_.place != analysis_var_placement_ )
        throw csmp::Exception( ERROR, caller,
                               "'stress' placement does not match "
                               "target placement" );
    if ( failure_key_.place != analysis_var_placement_ )
        throw csmp::Exception( ERROR, caller,
                               "'failure' placement does not match "
                               "target placement" );
    if ( failure01_key_.place != analysis_var_placement_ )
        throw csmp::Exception( ERROR, caller,
                               "'tensile failure' placement does not match "
                               "target placement" );

    if ( stress_key_.type != TENSOR )
        throw csmp::Exception( ERROR, caller,
                               "'stress' must be a tensor variable" );
    if ( cohesion_key_.type    != SCALAR ||
         friction_key_.type    != SCALAR ||
         tensile_str_key_.type != SCALAR ||
         biot_alpha_key_.type  != SCALAR ||
         failure_key_.type     != SCALAR ||
         failure01_key_.type   != SCALAR )
        throw csmp::Exception( ERROR, caller,
                               "one or more scalar variables have the "
                               "wrong type" );

    auto requireVariablePlacement = [&]( const csmp::Index& key,
                                         const char*        name )
    {
        if ( key.place != analysis_var_placement_ )
            throw csmp::Exception( ERROR, caller,
                                   name, parsePlacement(analysis_var_placement_) );
    };

    requireVariablePlacement( cohesion_key_,    "cohesion"         );
    requireVariablePlacement( friction_key_,    "friction angle"   );
    requireVariablePlacement( tensile_str_key_, "tensile strength" );
    requireVariablePlacement( biot_alpha_key_,  "Biot alpha"       );

    if ( pf_key_.place != UNDEFINED && pf_key_.place != NODE )
        throw csmp::Exception( ERROR, caller,
                               "'fluid pressure' must be placed at NODE" );

    if constexpr ( verbose_ )
        cout << "\nMohrCoulombFailure_Visitor: construction checks passed\n"
             << "  placement of analysis variables: "
             << ( analysis_var_placement_ == ELEMENT_INTEGRATION_POINT
                  ? "ELEMENT_INTEGRATION_POINT" : "ELEMENT" ) << "\n";
                  
    // setting failure criteria to zero
    model.InputPropertyValue("failure", makeScalar(ANY,0.) );
    model.InputPropertyValue("tensile failure", makeScalar(ANY,0.) );
}



// ----------------------------------------------------------------------------
//  Private constructor — 2D only, called by For2D factory
// ----------------------------------------------------------------------------

template<>
MohrCoulombFailure_Visitor<2U>::MohrCoulombFailure_Visitor( Model<2U>&      model,
                                                            PlaneAssumption assumption,
                                                            PLACEMENT       target_placement )
    : Visitor<2U>( MODEL, ELEMENT ),
      plane_assumption_( assumption       ),
      analysis_var_placement_( target_placement ),
      stress_key_      ( model.Database().StorageKey( "stress"           ) ),
      cohesion_key_    ( model.Database().StorageKey( "cohesion"         ) ),
      friction_key_    ( model.Database().StorageKey( "friction angle"   ) ),
      tensile_str_key_ ( model.Database().StorageKey( "tensile strength" ) ),
      biot_alpha_key_  ( model.Database().StorageKey( "Biot alpha"       ) ),
      nu_key_          ( model.Database().IsDefined( "Poissons ratio" )
                         ? model.Database().StorageKey( "Poissons ratio" )
                         : csmp::Index()                                   ),
      failure_key_     ( model.Database().StorageKey( "failure"          ) ),
      failure01_key_   ( model.Database().StorageKey( "tensile failure"  ) ),
      pf_key_          ( model.Database().IsDefined( "fluid pressure" )
                         ? model.Database().StorageKey( "fluid pressure" )
                         : csmp::Index()                                   )
{
    const char* caller = "MohrCoulombFailure_Visitor<2U> (constructor)";

    if ( analysis_var_placement_ != ELEMENT &&
         analysis_var_placement_ != ELEMENT_INTEGRATION_POINT )
        throw csmp::Exception( ERROR, caller,
                               "target placement must be ELEMENT or "
                               "ELEMENT_INTEGRATION_POINT" );

    if ( stress_key_.place != analysis_var_placement_ )
        throw csmp::Exception( ERROR, caller,
                               "'stress' placement does not match "
                               "target placement" );
    if ( failure_key_.place != analysis_var_placement_ )
        throw csmp::Exception( ERROR, caller,
                               "'failure' placement does not match "
                               "target placement" );
    if ( failure01_key_.place != analysis_var_placement_ )
        throw csmp::Exception( ERROR, caller,
                               "'tensile failure' placement does not match "
                               "target placement" );

    if ( stress_key_.type != TENSOR )
        throw csmp::Exception( ERROR, caller,
                               "'stress' must be a tensor variable" );
    if ( cohesion_key_.type    != SCALAR ||
         friction_key_.type    != SCALAR ||
         tensile_str_key_.type != SCALAR ||
         biot_alpha_key_.type  != SCALAR ||
         failure_key_.type     != SCALAR ||
         failure01_key_.type   != SCALAR )
        throw csmp::Exception( ERROR, caller,
                               "one or more scalar variables have the "
                               "wrong type" );

    auto requireVariablePlacement = [&]( const csmp::Index& key,
                                        const char*        name )
    {
        if ( key.place != analysis_var_placement_ )
            throw csmp::Exception( ERROR, caller,
                                   name, parsePlacement(analysis_var_placement_) );
    };

    requireVariablePlacement( cohesion_key_,    "cohesion"         );
    requireVariablePlacement( friction_key_,    "friction angle"   );
    requireVariablePlacement( tensile_str_key_, "tensile strength" );
    requireVariablePlacement( biot_alpha_key_,  "Biot alpha"       );

    if ( pf_key_.place != UNDEFINED && pf_key_.place != NODE )
        throw csmp::Exception( ERROR, caller,
                               "'fluid pressure' must be placed at NODE" );

    if constexpr ( verbose_ )
        cout << "\nMohrCoulombFailure_Visitor<2U>: construction checks passed\n"
             << "  placement of analysis variables: "
             << ( analysis_var_placement_ == ELEMENT_INTEGRATION_POINT
                  ? "ELEMENT_INTEGRATION_POINT" : "ELEMENT" ) << "\n"
             << "  plane assumption: "
             << ( plane_assumption_ == PlaneAssumption::PLANE_STRAIN
                  ? "PLANE_STRAIN" : "PLANE_STRESS" ) << "\n";

    // setting failure criteria to zero
    model.InputPropertyValue("failure", makeScalar(ANY,0.) );
    model.InputPropertyValue("failure01", makeScalar(ANY,0.) );
}





// ----------------------------------------------------------------------------
//  evaluateFailure — single-point computation
// ----------------------------------------------------------------------------

template<uint32_t dim>
void MohrCoulombFailure_Visitor<dim>::ComputeFailure(
    const StressInvariants& si,
    double                  cohesion,
    double                  friction_deg,
    double                  tensile_str,
    double                  biot_alpha,
    double                  pf,
    double&                 Fmc,
    double&                 F01 ) noexcept
{
    const double p_eff   = si.MeanStress() - biot_alpha * pf;
    const double K_theta = si.SmoothMohrCoulombYieldEnvelope( friction_deg );
    const double phi_rad = friction_deg * ( CSMP_PI / 180. );
    const double sin_phi = std::sin( phi_rad );
    const double cos_phi = std::cos( phi_rad );

    // Mohr-Coulomb shear failure criterion:
    //   F_mc = p_eff · sin φ - c · cos φ + q · K(θ)
    //   F_mc ≥ 0 means shear failure
    Fmc = p_eff * sin_phi - cohesion * cos_phi + si.DeviatoricStress() * K_theta;

    const double sigma3_eff = si.LeastPrincipalStress3() - biot_alpha * pf;

    // cap tensile strength at the Mohr-Coulomb apex c/tan(φ) [Pa]:
    // beyond the apex the envelope has no physical meaning.
    // guard against φ = 0 (frictionless) which gives an infinite apex.
    double max_tensile = tensile_str;
    if ( friction_deg > 0.001 ) {
        const double apex_limit = cohesion / std::tan( phi_rad );
        if ( tensile_str > apex_limit )
            max_tensile = apex_limit;
    }

    // tensile failure criterion (continuous measure, Pa):
    //   F01 > 0 means tensile failure
    //   F01 = -σ₃_eff - T_s
    //   example: σ₃_eff = -5 MPa, T_s = 2 MPa → F01 = 5 - 2 = 3 MPa (failing)
    F01 = -sigma3_eff - max_tensile;
}




template<uint32_t dim>
void MohrCoulombFailure_Visitor<dim>::EvaluateFailure(
                          const TensorVariable<dim>& stress,
                          double                     cohesion,
                          double                     friction_deg,
                          double                     tensile_str,
                          double                     biot_alpha,
                          double                     nu,
                          PlaneAssumption            condition,
                          double                     pf,
                          double&                    Fmc,
                          double&                    F01 ) const noexcept
{
    if constexpr ( dim == 1U ) {
        // Mohr-Coulomb failure is not defined for 1D stress states (false && always triggers)
        assert( false && "MohrCoulombFailure_Visitor: Mohr-Coulomb failure "
               "is not defined for 1D stress states" );
        Fmc = 0.;
        F01 = 0.;
        return;
    }
    else if constexpr ( dim == 2U ) {
        ComputeFailure( StressInvariants( stress, condition, nu ),
                        cohesion, friction_deg, tensile_str,
                        biot_alpha, pf, Fmc, F01 );
    }
    else if constexpr ( dim == 3U ) {
        ComputeFailure( StressInvariants( stress ),
                        cohesion, friction_deg, tensile_str,
                        biot_alpha, pf, Fmc, F01 );
    }
}



// ----------------------------------------------------------------------------
//  Visit
// ----------------------------------------------------------------------------

template<uint32_t dim>
void MohrCoulombFailure_Visitor<dim>::Visit( Element<dim>* e )
{
    // Mohr-Coulomb failure is only meaningful for volume elements —
    // face and interface elements do not carry a full stress tensor
    // TODO: this restriction is dangerous and must be removed after testing
    if constexpr (dim == 3 ) if ( !e->IsVolume() ) return;
    if constexpr (dim == 2 ) if ( !e->IsSurface() ) return;

    // --- integration point placement ---
    if ( stress_key_.place == ELEMENT_INTEGRATION_POINT ) {

        const uint32_t n_ip = e->IntegrationPoints();
        for ( uint32_t i{ 0U }; i < n_ip; ++i ) {

            TensorVariable<dim> stress;
            e->Read( i, stress_key_, stress );

            ScalarVariable cohesion_var, friction_var,
                           tensile_var,  biot_var;
            e->Read( i, cohesion_key_,    cohesion_var );
            e->Read( i, friction_key_,    friction_var );
            e->Read( i, tensile_str_key_, tensile_var  );
            e->Read( i, biot_alpha_key_,  biot_var     );

            // Poisson's ratio is only needed for 2D plane strain
            double nu = 0.;
            if constexpr ( dim == 2U ) {
                ScalarVariable nu_var;
                e->Read( i, nu_key_, nu_var );
                nu = nu_var();
            }

            double pf = 0.;
            if ( pf_key_.place != UNDEFINED )
                pf = e->PropertyValueAtIntegrationPoint( pf_key_, i );

            double Fmc = 0.;
            double F01 = 0.;
            EvaluateFailure( stress,
                             cohesion_var(),
                             friction_var(),
                             tensile_var(),
                             biot_var(),
                             nu,
                             plane_assumption_,
                             pf,
                             Fmc, F01 );

            e->Store( i, failure_key_,   makeScalar( PLAIN, Fmc ) );
            e->Store( i, failure01_key_, makeScalar( PLAIN, F01 ) );

            if ( verbose_ )
                cout << "  ip " << i
                     << "  F_mc=" << Fmc
                     << "  F01="  << F01 << "\n";
        }
        return;
    }

    // --- element placement ---
    if ( stress_key_.place == ELEMENT ) {

        TensorVariable<dim> stress;
        e->Read( stress_key_, stress );

        ScalarVariable cohesion_var, friction_var,
                       tensile_var,  biot_var;
        e->Read( cohesion_key_,    cohesion_var );
        e->Read( friction_key_,    friction_var );
        e->Read( tensile_str_key_, tensile_var  );
        e->Read( biot_alpha_key_,  biot_var     );

        // Poisson's ratio is only needed for 2D plane strain
        double nu = 0.;
        if constexpr ( dim == 2U ) {
            ScalarVariable nu_var;
            e->Read( nu_key_, nu_var );
            nu = nu_var();
        }

        double pf = 0.;
        if ( pf_key_.place != UNDEFINED ) {
            ScalarVariable sc;
            e->PropertyValueAtBaryCenter( pf_key_, sc );
            pf = sc();
        }

        double Fmc = 0.;
        double F01 = 0.;
        EvaluateFailure( stress,
                         cohesion_var(),
                         friction_var(),
                         tensile_var(),
                         biot_var(),
                         nu,
                         plane_assumption_,
                         pf,
                         Fmc, F01 );

        e->Store( failure_key_,   makeScalar( PLAIN, Fmc ) );
        e->Store( failure01_key_, makeScalar( PLAIN, F01 ) );

        if ( verbose_ )
            cout << "  element " << e->Idx()
                 << "  F_mc=" << Fmc
                 << "  F01="  << F01 << "\n";
    }

} // end Visit


template class MohrCoulombFailure_Visitor<1U>;
template class MohrCoulombFailure_Visitor<2U>;
template class MohrCoulombFailure_Visitor<3U>;

} // end csmp

