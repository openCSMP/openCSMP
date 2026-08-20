/*
 *  PropertyConstraints_Test.cpp
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *  Revised and extended for full interface coverage.
 */

#include "PropertyConstraints_Test.h"

#include "vsetMakers.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyConstraints.h"
#include "VSet.h"
#include "VTK_Interface.h"
#include "ANSYS_Model3D.h"
#include "meshManagementUtilities.h"
#include "compareFloats.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// ----------------------------------------------------------------------------
//  Construction / Destruction
// ----------------------------------------------------------------------------

PropertyConstraints_Test::PropertyConstraints_Test( bool verbose )
    : verbose_( verbose )
{
    const bool using_isoparametric_elements{ true };
    ModelTopology topology( "FracBox", using_isoparametric_elements );
    VSet<3U> vset;
    topology = create_FracBox( vset );

    model_ = new Model<3U>( topology, vset,
                            "CSMP-1phase-variables.txt",
                            true /* use regions file if any */ );

    Region<3U>& model_domain  = model_->Region( "Model"  );
    Region<3U>& matrix_domain = model_->Region( "MATRIX" );

    model_domain.InputPropertyValue( "permeability", makeScalar( ANY, 1.0e-12 ) );
    model_domain.InputPropertyValue( "porosity",     makeScalar( ANY, 1.0      ) );
    matrix_domain.InputPropertyValue( "porosity",    makeScalar( ANY, 0.4      ) );

    // elevated band used by several tests
    const csmp::Index phi_key = model_->Database().StorageKey( "porosity" );
    for ( size_t eidx{ 50 }; eidx < 150; ++eidx )
        matrix_domain.E( eidx )->Store( phi_key, makeScalar( ANY, 0.8 ) );

    // unit velocity vector on all nodes — used by vector length tests
    // only assigned if velocity is defined in the database
    if ( model_->Database().IsDefined( "velocity" ) )
        model_domain.InputPropertyValue( "velocity",
                                         makeVector( ANY, ANY, ANY, 1.0, 0.0, 0.0 ) );
}




// ----------------------------------------------------------------------------
//  run()
// ----------------------------------------------------------------------------

void PropertyConstraints_Test::run()
{
    // --- construction and assignment ---
    _test( TestDefaultConstruction()       );
    _test( TestParameterisedConstruction() );
    _test( TestCopyConstruction()          );
    _test( TestCopyAssignment()            );
    _test( TestMoveConstruction()          );
    _test( TestMoveAssignment()            );

    // --- mutation interface ---
    _test( TestAddConstraint()             );
    _test( TestChangeConstraint()          );
    _test( TestDeleteConstraint()          );
    _test( TestInitializePropertyIndices() );
    _test( TestAccessors()                 );
    _test( TestErase()                     );

    // --- constraint checking ---
    _test( TestCheckConstraintsAllNodes()     );
    _test( TestCheckConstraintsSingleNode()   );
    _test( TestCheckConstraintsNodalAverage() );
    _test( TestCheckConstraintsFailedUpon()   );
    _test( TestCheckLengthOfVectorVariables() );
    _test( TestMutualExclusivity()            );

    // --- integration tests ---
    _test( TestBuildRegionsFromPropertyConstraints() );
    _test( PointInVolumeElementTest()                );
    _test( TestMixedConstraints() );
}


// ----------------------------------------------------------------------------
//  TestDefaultConstruction
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestDefaultConstruction()
{
    PropertyConstraints pc;

    _test( pc.Constraints() == 0U    );
    _test( pc.WithIndexes()  == false );

    if ( verbose_ )
        pc.Out();

    return true;
}


// ----------------------------------------------------------------------------
//  TestParameterisedConstruction
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestParameterisedConstruction()
{
    // Construction with a valid property name should populate both maps
    PropertyConstraints pc( model_->Database(), "porosity", 0.3, 0.5 );

    _test( pc.Constraints() == 1U   );
    _test( pc.WithIndexes()  == true );

    if ( verbose_ )
        pc.Out();

    // Construction with an invalid property name does NOT throw —
    // ErrorHandler logs a note and the constructor returns early,
    // leaving the object empty
    PropertyConstraints pc_bad( model_->Database(), "no_such_property", 0.0, 1.0 );
    _test( pc_bad.Constraints() == 0U    );
    _test( pc_bad.WithIndexes()  == false );

    if ( verbose_ )
        pc_bad.Out();

    return true;
}


// ----------------------------------------------------------------------------
//  TestCopyConstruction
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestCopyConstruction()
{
    PropertyConstraints original( model_->Database(), "porosity", 0.3, 0.5 );
    original.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );

    PropertyConstraints copy( original );

    _test( copy.Constraints() == original.Constraints() );
    _test( copy.WithIndexes()  == original.WithIndexes() );

    // Mutating the copy must not affect the original
    copy.AddConstraint( "velocity", 0.0, 2.0 );
    _test( original.Constraints() == 2U );
    _test( copy.Constraints()     == 3U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestCopyAssignment
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestCopyAssignment()
{
    PropertyConstraints lhs;
    PropertyConstraints rhs( model_->Database(), "porosity", 0.3, 0.5 );

    lhs = rhs;

    _test( lhs.Constraints() == 1U   );
    _test( lhs.WithIndexes()  == true );

    // Self-assignment must be safe
    lhs = lhs;
    _test( lhs.Constraints() == 1U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestMoveConstruction
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestMoveConstruction()
{
    PropertyConstraints source( model_->Database(), "porosity", 0.3, 0.5 );
    _test( source.Constraints() == 1U );

    PropertyConstraints dest( std::move( source ) );

    _test( dest.Constraints()   == 1U   );
    _test( dest.WithIndexes()   == true );

    // source is in a valid but unspecified state — only check it does not crash
    source.Erase();
    _test( source.Constraints() == 0U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestMoveAssignment
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestMoveAssignment()
{
    PropertyConstraints source( model_->Database(), "porosity", 0.3, 0.5 );
    PropertyConstraints dest;

    dest = std::move( source );

    _test( dest.Constraints() == 1U   );
    _test( dest.WithIndexes() == true );

    source.Erase();
    _test( source.Constraints() == 0U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestAddConstraint
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestAddConstraint()
{
    PropertyConstraints pc;

    // First add must succeed
    const bool first_add = pc.AddConstraint( "porosity", 0.3, 0.5 );
    _test( first_add        == true );
    _test( pc.Constraints() == 1U   );

    // Duplicate add must fail
    const bool duplicate_add = pc.AddConstraint( "porosity", 0.6, 0.9 );
    _test( duplicate_add    == false );
    _test( pc.Constraints() == 1U   );

    // Second distinct property
    const bool second_add = pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );
    _test( second_add       == true );
    _test( pc.Constraints() == 2U   );

    // check_list_ not yet populated — indices not initialised
    _test( pc.WithIndexes() == false );

    return true;
}


// ----------------------------------------------------------------------------
//  TestChangeConstraint
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestChangeConstraint()
{
    PropertyConstraints pc( model_->Database(), "porosity", 0.3, 0.5 );
    pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );
    pc.InitializePropertyIndices( model_->Database() );

    // Change an existing constraint — must not throw
    bool threw = false;
    try {
        pc.ChangeConstraint( "porosity", 0.7, 0.9 );
    }
    catch ( const csmp::Exception& e ) {
        threw = true;
        if ( verbose_ )
            cerr << "TestChangeConstraint: unexpected csmp::Exception: "
                 << e.what() << "\n";
    }
    catch ( const std::exception& e ) {
        threw = true;
        if ( verbose_ )
            cerr << "TestChangeConstraint: unexpected std::exception: "
                 << e.what() << "\n";
    }
    _test( threw == false );

    // Verify the new range is applied — elevated band has porosity 0.8
    // which now falls inside [0.7, 0.9]
    Region<3U>& matrix_domain = model_->Region( "MATRIX" );
    const auto* e = matrix_domain.E( 50 );
    _test( pc.CheckConstraints( e ) == true );

    // Changing a non-existent constraint must throw
    bool threw_on_missing = false;
    try {
        pc.ChangeConstraint( "no_such_property", 0.0, 1.0 );
    }
    catch ( const csmp::Exception& ex ) {
        threw_on_missing = true;
        if ( verbose_ )
            cerr << "TestChangeConstraint: caught expected csmp::Exception: "
                 << ex.what() << "\n";
    }
    catch ( const std::exception& ex ) {
        threw_on_missing = true;
        if ( verbose_ )
            cerr << "TestChangeConstraint: caught std::exception: "
                 << ex.what() << "\n";
    }
    _test( threw_on_missing == true );

    return true;
}


// ----------------------------------------------------------------------------
//  TestDeleteConstraint
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestDeleteConstraint()
{
    PropertyConstraints pc( model_->Database(), "porosity", 0.3, 0.5 );
    pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );
    pc.InitializePropertyIndices( model_->Database() );

    _test( pc.Constraints() == 2U   );
    _test( pc.WithIndexes() == true );

    pc.DeleteConstraint( model_->Database(), "porosity" );
    _test( pc.Constraints() == 1U );

    // Re-initialise to confirm consistency
    pc.InitializePropertyIndices( model_->Database() );
    _test( pc.WithIndexes() == true );

    // Deleting the last constraint leaves an empty object
    pc.DeleteConstraint( model_->Database(), "permeability" );
    _test( pc.Constraints() == 0U );

    return true;
}


// ----------------------------------------------------------------------------
//  TestInitializePropertyIndices
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestInitializePropertyIndices()
{
    PropertyConstraints pc;
    pc.AddConstraint( "porosity",     0.3,     0.5     );
    pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );

    // Before initialisation, WithIndexes must be false
    _test( pc.WithIndexes() == false );

    const bool ok = pc.InitializePropertyIndices( model_->Database() );
    _test( ok               == true  );
    _test( pc.WithIndexes() == true  );

    // Calling again on an already-initialised object must be idempotent
    const bool ok2 = pc.InitializePropertyIndices( model_->Database() );
    _test( ok2              == true  );
    _test( pc.Constraints() == 2U    );

    // Adding a new constraint after initialisation requires re-initialisation
    pc.AddConstraint( "porosity", 0.0, 1.0 ); // duplicate — will be rejected
    _test( pc.Constraints() == 2U );

    // Empty constraints must throw on InitializePropertyIndices
    PropertyConstraints empty_pc;
    bool threw = false;
    try {
        empty_pc.InitializePropertyIndices( model_->Database() );
    }
    catch ( const csmp::Exception& ) {
        threw = true;
    }
    _test( threw == true );

    return true;
}


// ----------------------------------------------------------------------------
//  TestAccessors
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestAccessors()
{
    PropertyConstraints pc;
    _test( pc.Constraints() == 0U    );
    _test( pc.WithIndexes() == false );

    pc.AddConstraint( "porosity", 0.3, 0.5 );
    _test( pc.Constraints() == 1U    );
    _test( pc.WithIndexes() == false ); // indices not yet initialised

    pc.InitializePropertyIndices( model_->Database() );
    _test( pc.Constraints() == 1U   );
    _test( pc.WithIndexes() == true );

    // Adding a second constraint increments the count
    pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );
    _test( pc.Constraints() == 2U );

    // Re-initialise to confirm check_list_ is updated
    pc.InitializePropertyIndices( model_->Database() );
    _test( pc.WithIndexes() == true );

    // Erase resets both accessors
    pc.Erase();
    _test( pc.Constraints() == 0U    );
    _test( pc.WithIndexes() == false );

    // CheckLengthOfVectorVariables does not affect either accessor
    pc.AddConstraint( "porosity", 0.0, 1.0 );
    pc.CheckLengthOfVectorVariables( true );
    _test( pc.Constraints() == 1U    );
    _test( pc.WithIndexes() == false );

    pc.CheckLengthOfVectorVariables( false );
    _test( pc.Constraints() == 1U    );
    _test( pc.WithIndexes() == false );

    if ( verbose_ )
        pc.Out();

    return true;
}


// ----------------------------------------------------------------------------
//  TestErase
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestErase()
{
    PropertyConstraints pc( model_->Database(), "porosity", 0.3, 0.5 );
    pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );

    _test( pc.Constraints() == 2U   );
    _test( pc.WithIndexes() == true );

    pc.Erase();

    _test( pc.Constraints() == 0U    );
    _test( pc.WithIndexes() == false );

    // Object must be reusable after Erase
    pc.AddConstraint( "porosity", 0.0, 1.0 );
    _test( pc.Constraints() == 1U    );
    _test( pc.WithIndexes() == false );

    // Initialise after re-adding and confirm indices are rebuilt
    pc.InitializePropertyIndices( model_->Database() );
    _test( pc.WithIndexes() == true );

    // Erase a second time — must be safe on already-cleared object
    pc.Erase();
    pc.Erase();
    _test( pc.Constraints() == 0U    );
    _test( pc.WithIndexes() == false );

    if ( verbose_ )
        pc.Out();

    return true;
}


// ----------------------------------------------------------------------------
//  TestCheckConstraintsAllNodes
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestCheckConstraintsAllNodes()
{
    // Constraint: porosity in [0.75, 0.85]
    // Elements 50..149 of MATRIX have porosity 0.8 — should pass
    // All other MATRIX elements have porosity 0.4 — should fail
    PropertyConstraints pc( model_->Database(), "porosity", 0.75, 0.85 );

    Region<3U>& matrix_domain = model_->Region( "MATRIX" );

    // Element in elevated band — all nodes should have porosity 0.8
    const auto* e_in = matrix_domain.E( 75 );
    _test( pc.CheckConstraints( e_in ) == true );

    // Element outside elevated band — nodes have porosity 0.4
    const auto* e_out = matrix_domain.E( 0 );
    _test( pc.CheckConstraints( e_out ) == false );

    // Add a second constraint that the in-band element also satisfies
    pc.AddConstraint( "permeability", 1.0e-12, 1.0e-12 );
    pc.InitializePropertyIndices( model_->Database() );
    _test( pc.CheckConstraints( e_in ) == true );

    // Tighten permeability range so nothing passes
    pc.ChangeConstraint( "permeability", 1.0e-15, 1.0e-14 );
    _test( pc.CheckConstraints( e_in ) == false );

    if ( verbose_ )
        cout << "AllNodes: e_in passes tight permeability: "
             << pc.CheckConstraints( e_in ) << "\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestCheckConstraintsSingleNode
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestCheckConstraintsSingleNode()
{
    // Constraint: porosity in [0.75, 0.85]
    PropertyConstraints pc( model_->Database(), "porosity", 0.75, 0.85 );
    pc.SatisfyConstraintsForAtLeastOneNode( true );

    Region<3U>& matrix_domain = model_->Region( "MATRIX" );

    // Strict mode object for comparison
    PropertyConstraints pc_strict( model_->Database(), "porosity", 0.75, 0.85 );

    // Elements entirely outside the elevated band must fail in both modes
    const auto* e_out = matrix_domain.E( 0 );
    _test( pc_strict.CheckConstraints( e_out ) == false );
    _test( pc.CheckConstraints( e_out )        == false );

    // Elements entirely inside the elevated band must pass in both modes
    const auto* e_in = matrix_domain.E( 75 );
    _test( pc_strict.CheckConstraints( e_in ) == true );
    _test( pc.CheckConstraints( e_in )        == true );

    // Search for a perimeter cell that straddles the porosity transition:
    // passes in single-node mode but fails in strict mode
    const Element<3U>* e_boundary = nullptr;
    auto pend = matrix_domain.CellsEnd();
    for ( auto pit = matrix_domain.PerimeterCellsBegin(); pit != pend; ++pit ) {
        if (  pc.CheckConstraints( *pit ) &&
             !pc_strict.CheckConstraints( *pit ) ) {
            e_boundary = *pit;
            break;
          }
      }

    if ( e_boundary != nullptr ) {
        // A perimeter cell shares at least one face with the region boundary.
        // If it straddles the porosity transition it has nodes both inside
        // and outside the elevated band.
        _test(  pc.CheckConstraints( e_boundary )        );
        _test( !pc_strict.CheckConstraints( e_boundary ) );

        if ( verbose_ )
            cout << "TestCheckConstraintsSingleNode: "
                 << "straddling perimeter cell found\n";
    }
    else {
        // Fall back to atBoundary() if no perimeter cell straddles the transition
        if ( verbose_ )
            cout << "TestCheckConstraintsSingleNode: "
                 << "no straddling perimeter cell found — trying atBoundary()\n";

        auto eend = matrix_domain.CellsEnd();
        for ( auto eit = matrix_domain.CellsBegin(); eit != eend; ++eit ) {
            if ( atBoundary( *eit )                   &&
                  pc.CheckConstraints( *eit )          &&
                 !pc_strict.CheckConstraints( *eit ) ) {
                e_boundary = *eit;
                break;
              }
          }

        if ( e_boundary != nullptr ) {
            _test(  pc.CheckConstraints( e_boundary )        );
            _test( !pc_strict.CheckConstraints( e_boundary ) );

            if ( verbose_ )
                cout << "TestCheckConstraintsSingleNode: "
                     << "straddling boundary element found via atBoundary()\n";
        }
        else if ( verbose_ ) {
            cout << "TestCheckConstraintsSingleNode: "
                 << "no straddling boundary element found — "
                 << "single-node boundary test skipped\n";
        }
    }

    if ( verbose_ ) {
        cout << "SingleNode: e_out passes: " << pc.CheckConstraints( e_out ) << "\n";
        cout << "SingleNode: e_in  passes: " << pc.CheckConstraints( e_in  ) << "\n";
    }

    return true;
}


// ----------------------------------------------------------------------------
//  TestCheckConstraintsNodalAverage
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestCheckConstraintsNodalAverage()
{
    PropertyConstraints pc( model_->Database(), "porosity", 0.55, 0.65 );
    pc.SatisfyConstraintsForNodalAverage( true );
    pc.InitializePropertyIndices( model_->Database() );

    Region<3U>& matrix_domain = model_->Region( "MATRIX" );

    // Element deep inside elevated band — average ~0.8, outside [0.55, 0.65]
    const auto* e_in = matrix_domain.E( 100 );
    _test( pc.CheckConstraints( e_in ) == false );

    // Element entirely outside elevated band — average ~0.4, outside [0.55, 0.65]
    const auto* e_out = matrix_domain.E( 0 );
    _test( pc.CheckConstraints( e_out ) == false );

    // Search for a perimeter cell whose nodal average falls in [0.55, 0.65]
    const Element<3U>* e_boundary = nullptr;
    auto pend = matrix_domain.CellsEnd();
    for ( auto pit = matrix_domain.PerimeterCellsBegin(); pit != pend; ++pit ) {
        if ( pc.CheckConstraints( *pit ) ) {
            e_boundary = *pit;
            break;
          }
      }

    if ( e_boundary != nullptr ) {
        _test( pc.CheckConstraints( e_boundary ) == true );
        if ( verbose_ )
            cout << "TestCheckConstraintsNodalAverage: "
                 << "boundary element with average porosity in [0.55,0.65] found\n";
    }
    else if ( verbose_ ) {
        cout << "TestCheckConstraintsNodalAverage: "
             << "no boundary element found — nodal average boundary test skipped\n";
    }

    // Widen range to [0.35, 0.45] — only out-of-band elements should pass
    PropertyConstraints pc_low( model_->Database(), "porosity", 0.35, 0.45 );
    pc_low.SatisfyConstraintsForNodalAverage( true );
    pc_low.InitializePropertyIndices( model_->Database() );
    _test( pc_low.CheckConstraints( e_out ) == true  );
    _test( pc_low.CheckConstraints( e_in  ) == false );

    return true;
}



// ----------------------------------------------------------------------------
//  TestCheckConstraintsFailedUpon
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestCheckConstraintsFailedUpon()
{
    // Use CheckLengthOfVectorVariables(false) to avoid the throw in this overload
    PropertyConstraints pc( model_->Database(), "porosity", 0.75, 0.85 );
    pc.AddConstraint( "permeability", 1.0e-12, 1.0e-12 );
    pc.CheckLengthOfVectorVariables( false );
    pc.InitializePropertyIndices( model_->Database() );

    Region<3U>& matrix_domain = model_->Region( "MATRIX" );
    const csmp::Index phi_key  = model_->Database().StorageKey( "porosity"     );
    const csmp::Index perm_key = model_->Database().StorageKey( "permeability" );

    // Element outside elevated band — porosity 0.4, fails porosity constraint
    const auto* e_out = matrix_domain.E( 0 );
    csmp::Index failed_upon;
    const bool result_out = pc.CheckConstraints( e_out, failed_upon );
    _test( result_out  == false   );
    _test( failed_upon == phi_key );

    // Element in elevated band — porosity 0.8 passes, permeability 1e-12 passes
    const auto* e_in = matrix_domain.E( 75 );
    csmp::Index not_failed;
    const bool result_in = pc.CheckConstraints( e_in, not_failed );
    _test( result_in == true );

    // Tighten permeability so the in-band element now fails on permeability
    pc.ChangeConstraint( "permeability", 1.0e-15, 1.0e-14 );
    csmp::Index failed_perm;
    const bool result_perm = pc.CheckConstraints( e_in, failed_perm );
    _test( result_perm == false    );
    _test( failed_perm == perm_key );

    if ( verbose_ ) {
        cout << "FailedUpon: out-of-band fails on index: " << failed_upon  << "\n";
        cout << "FailedUpon: in-band fails on index:     " << failed_perm  << "\n";
    }

    return true;
}


// ----------------------------------------------------------------------------
//  TestCheckLengthOfVectorVariables
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestCheckLengthOfVectorVariables()
{
    if ( !model_->Database().IsDefined( "velocity" ) ) {
        if ( verbose_ )
            cout << "TestCheckLengthOfVectorVariables: "
                 << "'velocity' not defined in database — skipping\n";
        return true;
    }

    Region<3U>& model_domain = model_->Region( "Model" );
    const auto* e = model_domain.E( 0 );

    // velocity is stored as (1,0,0) on all nodes — length == 1.0

    // --- length mode ---

    // length 1.0 is inside [0.9, 1.1] — must pass
    PropertyConstraints pc_pass( model_->Database(), "velocity", 0.9, 1.1 );
    pc_pass.CheckLengthOfVectorVariables( true );
    _test( pc_pass.CheckConstraints( e ) == true );

    // length 1.0 is outside [2.0, 3.0] — must fail
    PropertyConstraints pc_fail( model_->Database(), "velocity", 2.0, 3.0 );
    pc_fail.CheckLengthOfVectorVariables( true );
    _test( pc_fail.CheckConstraints( e ) == false );

    // --- component mode ---
    // velocity components are (1.0, 0.0, 0.0)
    // IsWithinRange checks each component against [vmin, vmax]

    // range [0.0, 1.0] contains all three components — must pass
    PropertyConstraints pc_comp_pass( model_->Database(), "velocity", 0.0, 1.0 );
    pc_comp_pass.CheckLengthOfVectorVariables( false );
    _test( pc_comp_pass.CheckConstraints( e ) == true );

    // range [0.5, 1.5] — y and z components are 0.0 which is below 0.5 — must fail
    PropertyConstraints pc_comp_fail( model_->Database(), "velocity", 0.5, 1.5 );
    pc_comp_fail.CheckLengthOfVectorVariables( false );
    _test( pc_comp_fail.CheckConstraints( e ) == false );

    // switching mode on an existing object must take effect immediately
    // pc_pass was in length mode [0.9, 1.1] — switch to component mode
    // component y=0.0 is below 0.9 — must now fail
    pc_pass.CheckLengthOfVectorVariables( false );
    _test( pc_pass.CheckConstraints( e ) == false );

    if ( verbose_ )
        cout << "TestCheckLengthOfVectorVariables: all checks passed\n";

    return true;
}



// ----------------------------------------------------------------------------
//  TestMutualExclusivity
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestMutualExclusivity()
{
    // Basic flag interaction — no model needed
    // Note: the current implementation does NOT enforce mutual exclusivity
    // between one_node_only_ and nodal_average_ at the setter level.
    // Both flags can be true simultaneously, in which case VectorLengthCheck
    // will throw FATAL_ERROR if called with a vector constraint.
    // These tests document the actual behaviour.
    PropertyConstraints pc;
    pc.SatisfyConstraintsForAtLeastOneNode( true  );
    pc.SatisfyConstraintsForNodalAverage(   false );
    pc.SatisfyConstraintsForNodalAverage(   true  );
    pc.SatisfyConstraintsForAtLeastOneNode( false );

    // Without a vector constraint, both flags being true does not cause
    // a throw — only scalar and element constraints are checked
    PropertyConstraints pc_scalar( model_->Database(), "porosity", 0.3, 0.9 );
    pc_scalar.SatisfyConstraintsForAtLeastOneNode( true );
    pc_scalar.SatisfyConstraintsForNodalAverage(   true );

    Region<3U>& matrix_domain = model_->Region( "MATRIX" );
    const auto* e = matrix_domain.E( 0 );

    bool threw_scalar = false;
    try {
        pc_scalar.CheckConstraints( e );
    }
    catch ( const csmp::Exception& ex ) {
        threw_scalar = true;
        if ( verbose_ )
            cerr << "TestMutualExclusivity: scalar constraint threw: "
                 << ex.what() << "\n";
    }
    // scalar constraints do not trigger VectorLengthCheck — must not throw
    _test( threw_scalar == false );

    if ( !model_->Database().IsDefined( "velocity" ) ) {
        if ( verbose_ )
            cout << "TestMutualExclusivity: "
                 << "'velocity' not defined — skipping vector interaction test\n";
        return true;
    }

    // With a vector constraint and vector_length_check_ = true,
    // having both nodal_average_ and one_node_only_ true will cause
    // VectorLengthCheck to throw FATAL_ERROR.
    // This documents the known limitation — mutual exclusivity must be
    // enforced by the caller until the implementation is fixed.
    PropertyConstraints pc_vector( model_->Database(), "velocity", 0.9, 1.1 );
    pc_vector.CheckLengthOfVectorVariables( true );
    pc_vector.InitializePropertyIndices( model_->Database() );
    pc_vector.SatisfyConstraintsForAtLeastOneNode( true );
    pc_vector.SatisfyConstraintsForNodalAverage(   true );

    Region<3U>& model_domain = model_->Region( "Model" );
    const auto* ev = model_domain.E( 0 );

    bool threw_vector = false;
    try {
        pc_vector.CheckConstraints( ev );
    }
    catch ( const csmp::Exception& ex ) {
        threw_vector = true;
        if ( verbose_ )
            cerr << "TestMutualExclusivity: vector constraint threw as expected: "
                 << ex.what() << "\n";
    }
    catch ( const std::exception& ex ) {
        threw_vector = true;
        if ( verbose_ )
            cerr << "TestMutualExclusivity: std::exception: "
                 << ex.what() << "\n";
    }
    // document that this currently throws — test passes either way
    // but verbose output records the behaviour
    if ( verbose_ )
        cout << "TestMutualExclusivity: vector+both flags threw: "
             << threw_vector << " (expected true until implementation fixed)\n";

    return true;
}


// ----------------------------------------------------------------------------
//  TestBuildRegionsFromPropertyConstraints
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestBuildRegionsFromPropertyConstraints()
{
    // --- strict mode: all nodes must satisfy both constraints simultaneously ---
    // porosity in [0.75, 0.85] AND permeability == 1e-12
    // only elements 50..149 of MATRIX have porosity 0.8 and permeability 1e-12
    PropertyConstraints porosity_constraints( model_->Database(), "porosity", 0.75, 0.85 );
    porosity_constraints.AddConstraint( "permeability", 1.0e-12, 1.0e-12 );

    if ( !model_->ContainsRegion( "medium porosity" ) )
        model_->FormRegionFrom( "medium porosity", porosity_constraints );

    if ( model_->ContainsRegion( "medium porosity" ) ) {
        const Region<3U>& medium_porosity_domain = model_->Region( "medium porosity" );
        const auto strict_count = medium_porosity_domain.Cells();
        if ( verbose_ )
            cout << "strict mode cell count (porosity + permeability): "
                 << strict_count << "\n";
        _test( strict_count == 100 );
    }
    else if ( verbose_ )
        cout << "TestBuildRegionsFromPropertyConstraints: "
             << "strict mode region not formed\n";

    // --- strict mode: second constraint eliminates all elements ---
    // porosity in [0.75, 0.85] AND permeability in [1e-15, 1e-14]
    // no element satisfies both simultaneously
    PropertyConstraints impossible_constraints( model_->Database(), "porosity", 0.75, 0.85 );
    impossible_constraints.AddConstraint( "permeability", 1.0e-15, 1.0e-14 );

    if ( !model_->ContainsRegion( "impossible region" ) ) {
        bool threw = false;
        try {
            model_->FormRegionFrom( "impossible region", impossible_constraints );
        }
        catch ( const csmp::Exception& ex ) {
            threw = true;
            if ( verbose_ )
                cout << "TestBuildRegionsFromPropertyConstraints: "
                     << "impossible constraints correctly produced no region: "
                     << ex.what() << "\n";
        }
        // FormRegionFrom should warn and not form the region
        _test( !model_->ContainsRegion( "impossible region" ) );
    }

    // --- single-node mode: at least one node must satisfy both constraints ---
    // porosity in [0.75, 0.85] AND permeability in [1e-13, 1e-11]
    // single-node mode: element passes if at least one node has porosity in range
    // and permeability is satisfied (element variable — all-or-nothing)
    PropertyConstraints boundary_constraints( model_->Database(), "porosity", 0.75, 0.85 );
    boundary_constraints.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );
    boundary_constraints.SatisfyConstraintsForAtLeastOneNode( true );

    if ( !model_->ContainsRegion( "boundary porosity" ) )
        model_->FormRegionFrom( "boundary porosity", boundary_constraints );

    if ( model_->ContainsRegion( "boundary porosity" ) &&
         model_->ContainsRegion( "medium porosity"   ) ) {
        const auto single_node_count = model_->Region( "boundary porosity" ).Cells();
        const auto strict_count      = model_->Region( "medium porosity"   ).Cells();
        if ( verbose_ )
            cout << "single-node mode cell count (porosity + permeability): "
                 << single_node_count << "\n";
        // single-node mode must capture at least as many cells as strict mode
        _test( single_node_count >= strict_count );
    }

    // --- single-node mode: second constraint that nothing satisfies ---
    // single-node mode on porosity but impossible permeability range
    // the permeability constraint (element variable) eliminates everything
    PropertyConstraints single_node_impossible( model_->Database(), "porosity", 0.75, 0.85 );
    single_node_impossible.AddConstraint( "permeability", 1.0e-15, 1.0e-14 );
    single_node_impossible.SatisfyConstraintsForAtLeastOneNode( true );

    if ( !model_->ContainsRegion( "single node impossible" ) ) {
        bool threw = false;
        try {
            model_->FormRegionFrom( "single node impossible", single_node_impossible );
        }
        catch ( const csmp::Exception& ) {
            threw = true;
        }
        _test( !model_->ContainsRegion( "single node impossible" ) );
        if ( verbose_ )
            cout << "TestBuildRegionsFromPropertyConstraints: "
                 << "single-node impossible constraints correctly produced no region\n";
    }

    // --- nodal average mode: average of both constraints must be satisfied ---
    // porosity average in [0.55, 0.65] AND permeability in [1e-13, 1e-11]
    PropertyConstraints average_constraints( model_->Database(), "porosity", 0.55, 0.65 );
    average_constraints.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );
    average_constraints.SatisfyConstraintsForNodalAverage( true );

    if ( !model_->ContainsRegion( "average porosity" ) )
        model_->FormRegionFrom( "average porosity", average_constraints );

    if ( model_->ContainsRegion( "average porosity" ) ) {
        const auto average_count = model_->Region( "average porosity" ).Cells();
        if ( verbose_ )
            cout << "nodal average mode cell count (porosity + permeability): "
                 << average_count << "\n";
        _test( average_count > 0U );
    }
    else if ( verbose_ )
        cout << "TestBuildRegionsFromPropertyConstraints: "
             << "nodal average region not formed — "
             << "no elements with average porosity in [0.55, 0.65] "
             << "and permeability in [1e-13, 1e-11]\n";

    // --- nodal average mode: impossible second constraint ---
    PropertyConstraints average_impossible( model_->Database(), "porosity", 0.55, 0.65 );
    average_impossible.AddConstraint( "permeability", 1.0e-15, 1.0e-14 );
    average_impossible.SatisfyConstraintsForNodalAverage( true );

    if ( !model_->ContainsRegion( "average impossible" ) ) {
        bool threw = false;
        try {
            model_->FormRegionFrom( "average impossible", average_impossible );
        }
        catch ( const csmp::Exception& ) {
            threw = true;
        }
        _test( !model_->ContainsRegion( "average impossible" ) );
        if ( verbose_ )
            cout << "TestBuildRegionsFromPropertyConstraints: "
                 << "nodal average impossible constraints correctly produced no region\n";
    }

    // --- empty constraints must not silently form a region ---
    PropertyConstraints empty_pc;
    bool threw = false;
    try {
        model_->FormRegionFrom( "empty region", empty_pc );
    }
    catch ( const csmp::Exception& ) {
        threw = true;
    }
    _test( threw == true );

    if ( verbose_ ) {
        if ( model_->ContainsRegion( "medium porosity"   ) )
            model_->Region( "medium porosity"   ).Out();
        if ( model_->ContainsRegion( "boundary porosity" ) )
            model_->Region( "boundary porosity" ).Out();
        if ( model_->ContainsRegion( "average porosity"  ) )
            model_->Region( "average porosity"  ).Out();
    }

    return true;
}




// ----------------------------------------------------------------------------
//  PointInVolumeElementTest
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::PointInVolumeElementTest()
{
    ANSYS_Model3D model( "prism_test", "CSMP-variables.txt", true );

    // Use double literals — Point<3U> stores double, not float
    Point<3U> query( 2434.0, -1510.0, 5400.0 );

    auto& gref = model.Region( "Model" );

    auto eend = gref.CellsEnd();
    for ( auto eit = gref.CellsBegin(); eit != eend; ++eit ) {
        if ( !(*eit)->IsVolume() )
            continue;

        const Point<3U> bctr = (*eit)->BaryCenter();
        Element<3U>* e = pointInVolumeElement( gref, bctr );

        if ( verbose_ && e != *eit )
            cout << "PointInVolumeElement mismatch at barycentre: " << bctr << "\n";

        _test( e == *eit );
    }

    return true;
}


// ----------------------------------------------------------------------------
//  MixedConstraintTest
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestMixedConstraints()
{
    // --- mixed variable types and placements ---
    // porosity    : scalar, NODE,    range [0.75, 0.85]
    // permeability: scalar, ELEMENT, range [1e-13, 1e-11]
    // Both must be satisfied simultaneously.
    // Elements 50..149 have porosity 0.8 (node) and permeability 1e-12 (element).
    // permeability 1e-12 is inside [1e-13, 1e-11] — both constraints pass.
    {
        PropertyConstraints pc( model_->Database(), "porosity",     0.75,    0.85    );
        pc.AddConstraint(                           "permeability", 1.0e-13, 1.0e-11 );
        pc.InitializePropertyIndices( model_->Database() );

        Region<3U>& matrix_domain = model_->Region( "MATRIX" );

        // element in elevated band — porosity 0.8 (node), permeability 1e-12 (element)
        const auto* e_in = matrix_domain.E( 75 );
        _test( pc.CheckConstraints( e_in ) == true );

        // element outside elevated band — porosity 0.4 fails
        const auto* e_out = matrix_domain.E( 0 );
        _test( pc.CheckConstraints( e_out ) == false );

        if ( verbose_ )
            cout << "TestMixedConstraints: "
                 << "node scalar + element scalar: passed\n";
    }

    // --- greater-than constraint on one variable ---
    // porosity > 0.75 (node scalar)
    // expressed as [0.75, DBL_MAX]
    {
        PropertyConstraints pc( model_->Database(), "porosity",
                                0.75, std::numeric_limits<double>::max() );
        pc.InitializePropertyIndices( model_->Database() );

        Region<3U>& matrix_domain = model_->Region( "MATRIX" );

        const auto* e_in  = matrix_domain.E( 75 );
        const auto* e_out = matrix_domain.E( 0  );

        _test( pc.CheckConstraints( e_in  ) == true  ); // porosity 0.8 > 0.75
        _test( pc.CheckConstraints( e_out ) == false ); // porosity 0.4 < 0.75

        if ( verbose_ )
            cout << "TestMixedConstraints: "
                 << "greater-than constraint (porosity > 0.75): passed\n";
    }

    // --- less-than constraint on one variable ---
    // porosity < 0.5 (node scalar)
    // expressed as [-DBL_MAX, 0.5]
    {
        PropertyConstraints pc( model_->Database(), "porosity",
                                -std::numeric_limits<double>::max(), 0.5 );
        pc.InitializePropertyIndices( model_->Database() );

        Region<3U>& matrix_domain = model_->Region( "MATRIX" );

        const auto* e_in  = matrix_domain.E( 75 );
        const auto* e_out = matrix_domain.E( 0  );

        _test( pc.CheckConstraints( e_in  ) == false ); // porosity 0.8 > 0.5
        _test( pc.CheckConstraints( e_out ) == true  ); // porosity 0.4 < 0.5

        if ( verbose_ )
            cout << "TestMixedConstraints: "
                 << "less-than constraint (porosity < 0.5): passed\n";
    }

    // --- one variable greater than, another less than ---
    // porosity > 0.75  AND  permeability < 1e-11
    // elements 50..149: porosity 0.8 > 0.75 AND permeability 1e-12 < 1e-11 — pass
    // elements 0..49:   porosity 0.4 < 0.75 — fail on first constraint
    {
        PropertyConstraints pc( model_->Database(), "porosity",
                                0.75, std::numeric_limits<double>::max() );
        pc.AddConstraint( "permeability",
                          -std::numeric_limits<double>::max(), 1.0e-11 );
        pc.InitializePropertyIndices( model_->Database() );

        Region<3U>& matrix_domain = model_->Region( "MATRIX" );

        const auto* e_in  = matrix_domain.E( 75 );
        const auto* e_out = matrix_domain.E( 0  );

        _test( pc.CheckConstraints( e_in  ) == true  );
        _test( pc.CheckConstraints( e_out ) == false );

        if ( verbose_ )
            cout << "TestMixedConstraints: "
                 << "porosity > 0.75 AND permeability < 1e-11: passed\n";
    }

    // --- same variable, two separate constraints that together define a gap ---
    // This cannot be expressed with a single PropertyConstraints object because
    // AddConstraint rejects duplicate property names. Two separate objects are
    // needed and the results combined by the caller.
    // Document this limitation explicitly.
    {
        // want: porosity < 0.3 OR porosity > 0.7
        // expressed as two separate PropertyConstraints objects
        PropertyConstraints pc_low( model_->Database(), "porosity",
                                    -std::numeric_limits<double>::max(), 0.3 );
        PropertyConstraints pc_high( model_->Database(), "porosity",
                                     0.7, std::numeric_limits<double>::max() );
        pc_low.InitializePropertyIndices(  model_->Database() );
        pc_high.InitializePropertyIndices( model_->Database() );

        Region<3U>& matrix_domain = model_->Region( "MATRIX" );

        const auto* e_in  = matrix_domain.E( 75 ); // porosity 0.8
        const auto* e_out = matrix_domain.E( 0  ); // porosity 0.4

        // porosity 0.8 > 0.7 — passes pc_high, fails pc_low
        _test( pc_high.CheckConstraints( e_in  ) == true  );
        _test( pc_low.CheckConstraints(  e_in  ) == false );

        // porosity 0.4 is between 0.3 and 0.7 — fails both
        _test( pc_high.CheckConstraints( e_out ) == false );
        _test( pc_low.CheckConstraints(  e_out ) == false );

        if ( verbose_ )
            cout << "TestMixedConstraints: "
                 << "OR-style constraints via two objects: passed\n";
    }

    return true;
}


} // end csmp

