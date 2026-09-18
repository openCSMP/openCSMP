/*
 *  MohrCoulombFailure_Visitor_Test.cpp
 *  csmp_core
 */

#include "MohrCoulombFailure_Visitor_Test.h"

#include "vsetMakers.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "MohrCoulombFailure_Visitor.h"
#include "StressInvariants.h"
#include "TensorVariable.h"
#include "ScalarVariable.h"
#include "CSMP_physical_constants.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// ----------------------------------------------------------------------------
//  Construction
// ----------------------------------------------------------------------------

MohrCoulombFailure_Visitor_Test::MohrCoulombFailure_Visitor_Test( bool verbose )
    : verbose_( verbose )
{
}


// ----------------------------------------------------------------------------
//  run
// ----------------------------------------------------------------------------

void MohrCoulombFailure_Visitor_Test::run()
{
    _test( TestOverburdenStressAssignment() );
    _test( TestVisitorRuns()               );
    _test( TestTensileFailureAtTop()       );
    _test( TestShearFailureAtDepth()       );
    _test( TestStableInterior()            );
    _test( TestOutputPlacement()           );
}


// ----------------------------------------------------------------------------
//  Helpers
// ----------------------------------------------------------------------------

namespace {

/**
    @brief Builds the FracBox model with all required variables.

    Returns a fully initialised Model<3> with:
    - stress tensor at element integration points
    - friction angle, cohesion, tensile strength, Biot alpha at elements
    - fluid pressure at nodes (zero throughout)
    - failure and tensile failure at element integration points
*/
Model<3U> buildModel()
{
    VSet<3U>      vset;
    const string  var_file( "MohrCoulombFailure_Visitor_Test-variables.txt" );
    ModelTopology topology = create_FracBox( vset );
    vset.EstablishElementConnectivity3D();
    vset.RemoveData("element number");
    vset.RemoveData("node number");
    return Model<3U>( topology, vset, var_file.c_str(), true );
}


/**
    @brief Assigns overburden stress to all volume elements.

    Stress is total (not effective) — the visitor applies the Biot
    effective stress correction internally using the fluid pressure
    and Biot coefficient.

    Stress is computed from the element barycentre depth z:

        σ_zz = ρ·g·z          (vertical overburden, compression positive)
        σ_xx = σ_yy = K₀·σ_zz (horizontal stress with earth pressure K₀)
        σ_xy = σ_xz = σ_yz = 0

    The model top is taken as the maximum y-coordinate in the model.
    Depth is measured downward from the top.

    @param model    The model to assign stress to.
    @param rho      Bulk density [kg/m³].
    @param K0       Earth pressure coefficient at rest [-] is used here for the testing.
*/
void assignOverburdenStress( Model<3U>& model,
                             double     rho,
                             double     K0 )
{
    const csmp::Index stress_key = model.Database().StorageKey( "stress" );

    Region<3U>& domain = model.Region( "Model" );

    // find the model top (maximum y-coordinate)
    double y_top = -std::numeric_limits<double>::max();
    for ( auto nit  = domain.NodesBegin();
               nit != domain.NodesEnd(); ++nit )
        if ( (*nit)->y() > y_top ) y_top = (*nit)->y();

    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit ) {
        auto* e = *eit;
        if ( !e->IsVolume() ) continue;

        const Point<3U> bc    = e->BaryCenter();
        const double    depth = y_top - bc[1];   // depth below model top
        const double    szz   = rho * ACC_GRAVITY * depth;  // compression positive
        const double    sxx   = K0 * szz;
        const double    syy   = K0 * szz;

        TensorVariable<3U> stress;
        stress(0,0) = sxx;  stress(1,1) = syy;  stress(2,2) = szz;
        stress(0,1) = 0.;   stress(1,0) = 0.;
        stress(0,2) = 0.;   stress(2,0) = 0.;
        stress(1,2) = 0.;   stress(2,1) = 0.;

        // store at all integration points
        for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i )
            e->Store( i, stress_key, stress );
    }
}


/**
    @brief Assigns uniform material properties to all volume elements.

    @param model          The model.
    @param friction_deg   Friction angle [degrees].
    @param cohesion       Cohesion [Pa].
    @param tensile_str    Tensile strength [Pa].
    @param biot_alpha     Biot coefficient [-].
*/
void assignMaterialProperties( Model<3U>& model,
                                double     friction_deg,
                                double     cohesion,
                                double     tensile_str,
                                double     biot_alpha )
{
    Region<3U>& domain = model.Region( "Model" );

    domain.InputPropertyValue( "friction angle",   makeScalar( ANY, friction_deg ) );
    domain.InputPropertyValue( "cohesion",         makeScalar( ANY, cohesion      ) );
    domain.InputPropertyValue( "tensile strength", makeScalar( ANY, tensile_str   ) );
    domain.InputPropertyValue( "Biot alpha",       makeScalar( ANY, biot_alpha    ) );
    domain.InputPropertyValue( "fluid pressure",   makeScalar( ANY, 100325.       ) ); // minimum from database
}

} // anonymous namespace


// ----------------------------------------------------------------------------
//  TestOverburdenStressAssignment
// ----------------------------------------------------------------------------

bool MohrCoulombFailure_Visitor_Test::TestOverburdenStressAssignment()
{
    Model<3U> model = buildModel();

    const double rho = 2500.;  // kg/m³ — typical rock density
    const double K0  = 0.5;    // earth pressure coefficient at rest

    assignOverburdenStress( model, rho, K0 );
    assignMaterialProperties( model, 30., 1.0e6, 0.1e6, 1.0 );
    
   if ( verbose_ )
        cout << "TestOverburdenStressAssignment: "
             << "fluid pressure correctly set to atmospheric\n";

    const csmp::Index stress_key = model.Database().StorageKey( "stress" );
    Region<3U>& domain = model.Region( "Model" );

    // verify fluid pressure is at atmospheric baseline
    const csmp::Index pf_key = model.Database().StorageKey( "fluid pressure" );
    for ( auto nit  = domain.NodesBegin();
               nit != domain.NodesEnd(); ++nit ) {
        ScalarVariable pf;
        (*nit)->Read( pf_key, pf );
        _equal( pf(), 100325., 1.0 );
    }

     // verify that stress has been stored and is non-zero at depth
    uint32_t nonzero_count = 0U;
    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit ) {
        auto* e = *eit;
        if ( !e->IsVolume() ) continue;

        TensorVariable<3U> stress;
        e->Read( 0U, stress_key, stress );

        // σ_zz must be non-negative (compression positive)
        _test( stress(2,2) >= 0. );

        // σ_xx and σ_yy must equal K0 * σ_zz
        _test( std::abs( stress(0,0) - K0 * stress(2,2) ) < 1.0e-6 );
        _test( std::abs( stress(1,1) - K0 * stress(2,2) ) < 1.0e-6 );

        // all shear components must be zero
        _test( stress(0,1) == 0. );
        _test( stress(0,2) == 0. );
        _test( stress(1,2) == 0. );

        if ( stress(2,2) > 0. ) ++nonzero_count;
    }

    // at least some elements must have non-zero stress (those at depth)
    _test( nonzero_count > 0U );

    if ( verbose_ )
        cout << "TestOverburdenStressAssignment: "
             << nonzero_count << " elements with non-zero overburden stress\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestVisitorRuns
// ----------------------------------------------------------------------------

bool MohrCoulombFailure_Visitor_Test::TestVisitorRuns()
{
    Model<3U> model = buildModel();
    assignOverburdenStress( model, 2500., 0.5 );
    assignMaterialProperties( model, 30., 1.0e6, 0.1e6, 1.0 );

    // diagnostic: count elements before visitor runs
    uint32_t element_count = 0U;
    Region<3U>& domain = model.Region( "Model" );
    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit )
        if ( (*eit)->IsVolume() ) ++element_count;

    if ( verbose_ )
        cout << "TestVisitorRuns: volume element count = "
             << element_count << "\n";

    _test( element_count > 0U );

    bool threw = false;
    try {
        MohrCoulombFailure_Visitor<3U> visitor( model, ELEMENT_INTEGRATION_POINT );

        if ( verbose_ )
            cout << "TestVisitorRuns: visitor constructed, "
                 << "ApplicationLevel=" << visitor.ApplicationLevel()
                 << " ApplicationTarget=" << visitor.ApplicationTarget() << "\n";

        model.Accept( visitor );

        if ( verbose_ )
            cout << "TestVisitorRuns: Accept returned\n";
    }
    catch ( const csmp::Exception& ex ) {
        threw = true;
        cerr << "TestVisitorRuns: caught csmp::Exception: "
             << ex.what() << "\n";
    }
    catch ( const std::exception& ex ) {
        threw = true;
        cerr << "TestVisitorRuns: caught std::exception: "
             << ex.what() << "\n";
    }

    _test( threw == false );

    // verify at least one integration point was written
    const csmp::Index failure_key = model.Database().StorageKey( "failure" );
    uint32_t written = 0U;
    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit ) {
        if ( !(*eit)->IsVolume() ) continue;
        for ( uint32_t i{ 0U }; i < (*eit)->IntegrationPoints(); ++i ) {
            ScalarVariable fmc;
            (*eit)->Read( i, failure_key, fmc );
            if ( std::isfinite( fmc() ) ) ++written;
        }
    }

    if ( verbose_ )
        cout << "TestVisitorRuns: integration points written = "
             << written << "\n";

    _test( written > 0U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestTensileFailureAtTop
// ----------------------------------------------------------------------------

bool MohrCoulombFailure_Visitor_Test::TestTensileFailureAtTop()
{
    // At the model top, depth = 0, so σ_zz = 0 and σ_xx = σ_yy = 0.
    // The least principal stress σ₃ = 0, which equals the tensile strength
    // threshold of 0.1 MPa only if σ₃ < -T_s. With zero stress there is no
    // tensile failure. However, with a very small tensile strength (near zero)
    // and K0 < 1, the horizontal stress is less than the vertical, so elements
    // near the top with low confinement are closest to tensile failure.
    //
    // To guarantee tensile failure at the top, we use a high density so that
    // even shallow elements develop enough deviatoric stress, combined with
    // a very small tensile strength.

    Model<3U> model = buildModel();

    // high density and low K0 to maximise deviatoric stress
    const double rho         = 2700.;
    const double K0          = 0.3;   // low horizontal confinement
    const double tensile_str = 10.;   // very small tensile strength [Pa]

    assignOverburdenStress( model, rho, K0 );
    assignMaterialProperties( model, 30., 1.0e6, tensile_str, 1.0 );

    MohrCoulombFailure_Visitor<3U> visitor( model, ELEMENT_INTEGRATION_POINT );
    model.Accept( visitor );

    const csmp::Index failure01_key = model.Database().StorageKey( "tensile failure" );
    Region<3U>& domain = model.Region( "Model" );

    // find the model top
    double y_top = -std::numeric_limits<double>::max();
    for ( auto nit  = domain.NodesBegin();
               nit != domain.NodesEnd(); ++nit )
        if ( (*nit)->y() > y_top ) y_top = (*nit)->y();

    // collect elements in the top 10% of the model height
    double y_min = std::numeric_limits<double>::max();
    double y_max = -std::numeric_limits<double>::max();
    for ( auto nit  = domain.NodesBegin();
               nit != domain.NodesEnd(); ++nit ) {
        if ( (*nit)->y() < y_min ) y_min = (*nit)->y();
        if ( (*nit)->y() > y_max ) y_max = (*nit)->y();
    }
    const double model_height = y_max - y_min;
    const double top_threshold = y_max - 0.1 * model_height;

    uint32_t tensile_failure_count = 0U;
    uint32_t top_element_count     = 0U;

    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit ) {
        auto* e = *eit;
        if ( !e->IsVolume() ) continue;

        const Point<3U> bc = e->BaryCenter();
        if ( bc[1] < top_threshold ) continue;

        ++top_element_count;

        // check all integration points
        for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i ) {
            ScalarVariable f01;
            e->Read( i, failure01_key, f01 );
            // tensile failure if f01 is greater than zero
            if ( f01() > 0. ) {
                ++tensile_failure_count;
                break;
            }
        }
    }

    if ( verbose_ )
        cout << "TestTensileFailureAtTop: "
             << tensile_failure_count << " of "
             << top_element_count
             << " top elements show tensile failure\n";

    // at least some top elements must show tensile failure
    _test( top_element_count     > 0U );
    _test( tensile_failure_count > 0U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestShearFailureAtDepth
// ----------------------------------------------------------------------------

bool MohrCoulombFailure_Visitor_Test::TestShearFailureAtDepth()
{
    // At depth, σ_zz is large and σ_xx = σ_yy = K0·σ_zz with K0 < 1,
    // so the deviatoric stress q = σ_zz - σ_xx is large.
    // With low cohesion, shear failure (F_mc >= 0) should occur at depth.

    Model<3U> model = buildModel();

    const double rho      = 2700.;
    const double K0       = 0.3;
    const double cohesion = 100.;   // very low cohesion [Pa] to ensure failure

    assignOverburdenStress( model, rho, K0 );
    assignMaterialProperties( model, 30., cohesion, 10., 1.0 );

    MohrCoulombFailure_Visitor<3U> visitor( model, ELEMENT_INTEGRATION_POINT );
    model.Accept( visitor );

    const csmp::Index failure_key = model.Database().StorageKey( "failure" );
    Region<3U>& domain = model.Region( "Model" );

    // find model bottom 10%
    double y_min = std::numeric_limits<double>::max();
    double y_max = -std::numeric_limits<double>::max();
    for ( auto nit  = domain.NodesBegin();
               nit != domain.NodesEnd(); ++nit ) {
        if ( (*nit)->y() < y_min ) y_min = (*nit)->y();
        if ( (*nit)->y() > y_max ) y_max = (*nit)->y();
    }
    const double model_height    = y_max - y_min;
    const double bottom_threshold = y_min + 0.1 * model_height;

    uint32_t shear_failure_count  = 0U;
    uint32_t bottom_element_count = 0U;

    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit ) {
        auto* e = *eit;
        if ( !e->IsVolume() ) continue;

        const Point<3U> bc = e->BaryCenter();
        if ( bc[1] > bottom_threshold ) continue;

        ++bottom_element_count;

        // check all integration points
        for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i ) {
            ScalarVariable fmc;
            e->Read( i, failure_key, fmc );
            if ( fmc() >= 0. ) {
                ++shear_failure_count;
                break;
            }
        }
    }

    if ( verbose_ )
        cout << "TestShearFailureAtDepth: "
             << shear_failure_count << " of "
             << bottom_element_count
             << " bottom elements show shear failure\n";

    _test( bottom_element_count > 0U );
    _test( shear_failure_count  > 0U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestStableInterior
// ----------------------------------------------------------------------------

bool MohrCoulombFailure_Visitor_Test::TestStableInterior()
{
    // With high cohesion and moderate stress, interior elements should
    // remain stable (F_mc < 0) under the overburden loading.

    Model<3U> model = buildModel();

    const double rho      = 2500.;
    const double K0       = 0.8;      // near-isotropic horizontal stress
    const double cohesion = 100.0e6;  // very high cohesion [Pa]

    assignOverburdenStress( model, rho, K0 );
    assignMaterialProperties( model, 30., cohesion, 10.0e6, 1.0 );

    MohrCoulombFailure_Visitor<3U> visitor( model, ELEMENT_INTEGRATION_POINT );
    model.Accept( visitor );

    const csmp::Index failure_key   = model.Database().StorageKey( "failure"        );
    const csmp::Index failure01_key = model.Database().StorageKey( "tensile failure" );
    Region<3U>& domain = model.Region( "Model" );

    // find model middle 20%
    double y_min = std::numeric_limits<double>::max();
    double y_max = -std::numeric_limits<double>::max();
    for ( auto nit  = domain.NodesBegin();
               nit != domain.NodesEnd(); ++nit ) {
        if ( (*nit)->y() < y_min ) y_min = (*nit)->y();
        if ( (*nit)->y() > y_max ) y_max = (*nit)->y();
    }
    const double model_height = y_max - y_min;
    const double mid_low      = y_min + 0.4 * model_height;
    const double mid_high     = y_min + 0.6 * model_height;

    uint32_t stable_count  = 0U;
    uint32_t failed_count  = 0U;
    uint32_t middle_count  = 0U;

    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit ) {
        auto* e = *eit;
        if ( !e->IsVolume() ) continue;

        const Point<3U> bc = e->BaryCenter();
        if ( bc[1] < mid_low || bc[1] > mid_high ) continue;

        ++middle_count;

        bool element_failed = false;
        for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i ) {
            ScalarVariable fmc, f01;
            e->Read( i, failure_key,   fmc );
            e->Read( i, failure01_key, f01 );
            if ( fmc() >= 0. || f01() > 0.5 ) {
                element_failed = true;
                break;
            }
        }

        if ( element_failed ) ++failed_count;
        else                  ++stable_count;
    }

    if ( verbose_ )
        cout << "TestStableInterior: "
             << stable_count << " stable, "
             << failed_count << " failed, of "
             << middle_count << " middle elements\n";

    _test( middle_count > 0U  );
    _test( stable_count > 0U  );
    _test( failed_count == 0U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestOutputPlacement
// ----------------------------------------------------------------------------

bool MohrCoulombFailure_Visitor_Test::TestOutputPlacement()
{
    Model<3U> model = buildModel();

    // must assign stress and material properties before running visitor
    assignOverburdenStress( model, 2500., 0.5 );
    assignMaterialProperties( model, 30., 1.0e6, 0.1e6, 1.0 );

    MohrCoulombFailure_Visitor<3U> visitor( model, ELEMENT_INTEGRATION_POINT );
    model.Accept( visitor );

    const csmp::Index stress_key    = model.Database().StorageKey( "stress"          );
    const csmp::Index failure_key   = model.Database().StorageKey( "failure"         );
    const csmp::Index failure01_key = model.Database().StorageKey( "tensile failure" );

    // stress, failure, and tensile failure must all be at integration points
    _test( stress_key.place    == ELEMENT_INTEGRATION_POINT );
    _test( failure_key.place   == ELEMENT_INTEGRATION_POINT );
    _test( failure01_key.place == ELEMENT_INTEGRATION_POINT );

    Region<3U>& domain = model.Region( "Model" );

    uint32_t total_ip   = 0U;
    uint32_t written_ip = 0U;

    for ( auto eit  = domain.CellsBegin();
               eit != domain.CellsEnd(); ++eit ) {
        auto* e = *eit;
        if ( !e->IsVolume() ) continue;

        for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i ) {
            ++total_ip;

            ScalarVariable fmc, f01;
            e->Read( i, failure_key,   fmc );
            e->Read( i, failure01_key, f01 );

            // F01 > 0 means tensile failure, F01 <= 0 means stable
            // no constraint on exact value — just checking if it is finite
            _test( std::isfinite( f01() ) );

            // F_mc must be a finite number
            _test( std::isfinite( fmc() ) );

            ++written_ip;
        }
    }

    _test( total_ip   > 0U        );
    _test( written_ip == total_ip );

    if ( verbose_ )
        cout << "TestOutputPlacement: "
             << written_ip << " of " << total_ip
             << " integration points written\n";

    return true;
}

} // end csmp

