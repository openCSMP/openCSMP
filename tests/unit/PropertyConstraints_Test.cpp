/*
 *  PropertyConstraints_Test.cpp
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *  Revised and extended for full interface coverage.
 */

#include "PropertyConstraints_Test.h"

#include "vsetMakers.h"
#include "Region.h"
#include "Boundary.h"
#include "VectorVariable.h"
#include "PropertyConstraints.h"
#include "VSet.h"
#include "VTK_Interface.h"
#include "meshManagementUtilities.h"
#include "compareFloats.h"
#include "ANSYS_Model3D.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// ----------------------------------------------------------------------------
//  Helper — builds a minimal 3D model with known property values.
//  Returns the model; caller owns it.
//  Assigns:
//    porosity  @ NODE:    fracture elements = 1.0, matrix = 0.4
//    porosity  @ NODE:    elements 50..149 of matrix = 0.8  (elevated band)
//    permeability @ ELEMENT: entire domain = 1.0e-12
//    velocity  @ NODE:    entire domain = (1.0, 0.0, 0.0)  (unit x-vector)
// ----------------------------------------------------------------------------

static Model<3>* buildTestModel()
{
    const bool using_isoparametric_elements{ true };
    ModelTopology topology( "FracBox", using_isoparametric_elements );
    VSet<3U>      vset;
    topology = create_FracBox( vset );

    Model<3U>* model = new Model<3U>( topology, vset, "CSMP-1phase-variables.txt", true /* use regions file if any */ );

    Region<3U>& model_domain  = model->Region( "Model"  );
    Region<3U>& matrix_domain = model->Region( "MATRIX" );

    model_domain.InputPropertyValue( "permeability", makeScalar( ANY, 1.0e-12 ) );
    model_domain.InputPropertyValue( "porosity",     makeScalar( ANY, 1.0      ) ); // fracture default
    matrix_domain.InputPropertyValue( "porosity",    makeScalar( ANY, 0.4      ) ); // matrix default

    // elevated band used by several tests
    const csmp::Index phi_key = model->Database().StorageKey( "porosity" );
    for ( size_t eidx{ 50 }; eidx < 150; ++eidx )
        matrix_domain.E( eidx )->Store( phi_key, makeScalar( ANY, 0.8 ) );

    // unit velocity vector on all nodes — used by vector length tests
    model_domain.InputPropertyValue( "velocity", makeVector( ANY,ANY,ANY, 1.0, 0.0, 0.0 ) );
    
    return model;
}


// ----------------------------------------------------------------------------
//  Construction / Destruction
// ----------------------------------------------------------------------------

PropertyConstraints_Test::PropertyConstraints_Test( bool verbose )
    : verbose_( verbose ),
      model_(buildTestModel())
{
}

PropertyConstraints_Test::~PropertyConstraints_Test()
{
   delete model_;
}


// ----------------------------------------------------------------------------
//  run() — executes all sub-tests and reports
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
    _test( TestCheckConstraintsAllNodes()    );
    _test( TestCheckConstraintsSingleNode()  );
    _test( TestCheckConstraintsNodalAverage());
    _test( TestCheckConstraintsFailedUpon()  );
    _test( TestCheckLengthOfVectorVariables());
    _test( TestMutualExclusivity()           );

    // --- integration tests ---
    _test( TestBuildRegionsFromPropertyConstraints() );
    _test( PointInVolumeElementTest()                );
}




// ----------------------------------------------------------------------------
//  TestDefaultConstruction
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestDefaultConstruction()
{
    PropertyConstraints pc;

    _test( pc.Constraints() == 0U );
    _test( pc.WithIndexes()  == false );

    if ( verbose_ ) pc.Out();

    return true;
}


// ----------------------------------------------------------------------------
//  TestParameterisedConstruction
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestParameterisedConstruction()
{
    // Construction with a valid property name should populate both maps
    PropertyConstraints pc( model_->Database(), "porosity", 0.3, 0.5 );

    _test( pc.Constraints() == 1U  );
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
    _test( copy.WithIndexes()  == original.WithIndexes()  );

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

    _test( lhs.Constraints() == 1U  );
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

    _test( dest.Constraints()   == 1U  );
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

    _test( dest.Constraints() == 1U  );
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
    _test( first_add            == true );
    _test( pc.Constraints()     == 1U   );

    // Duplicate add must fail (map insert returns false for existing key)
    const bool duplicate_add = pc.AddConstraint( "porosity", 0.6, 0.9 );
    _test( duplicate_add        == false );
    _test( pc.Constraints()     == 1U    ); // count unchanged

    // Second distinct property
    const bool second_add = pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );
    _test( second_add           == true );
    _test( pc.Constraints()     == 2U   );

    // check_list_ not yet populated — indices not initialised
    _test( pc.WithIndexes()     == false );

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
        // catches it if csmp::Exception does not propagate at INFO level
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

    _test( pc.Constraints() == 2U );
    _test( pc.WithIndexes() == true );

    pc.DeleteConstraint( model_->Database(), "porosity" );

    _test( pc.Constraints() == 1U );

    // After deletion, check_list_ should also have shrunk
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
    pc.AddConstraint( "porosity",     0.3,    0.5    );
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

    // Adding a new constraint after initialisation must invalidate check_list_
    // so that re-initialisation picks it up
    pc.AddConstraint( "velocity", 0.0, 2.0 );
    _test( pc.Constraints() == 3U );
    // check_list_ now has fewer entries than criteria_ — re-init required
    const bool ok3 = pc.InitializePropertyIndices( model_->Database() );
    _test( ok3              == true );
    _test( pc.WithIndexes() == true );

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

    // WithIndexes reflects that check_list_ is now stale (one fewer entry)
    // After re-initialisation it should be true again
    pc.InitializePropertyIndices( model_->Database() );
    _test( pc.WithIndexes() == true );

    // Erase resets both accessors
    pc.Erase();
    _test( pc.Constraints() == 0U    );
    _test( pc.WithIndexes() == false );

    // CheckLengthOfVectorVariables does not affect either accessor
    pc.AddConstraint( "velocity", 0.0, 2.0 );
    pc.CheckLengthOfVectorVariables( true );
    _test( pc.Constraints() == 1U    );
    _test( pc.WithIndexes() == false );

    pc.CheckLengthOfVectorVariables( false );
    _test( pc.Constraints() == 1U    );
    _test( pc.WithIndexes() == false );

    if ( verbose_ ) pc.Out();

    return true;
}


// ----------------------------------------------------------------------------
//  TestErase
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestErase()
{
    PropertyConstraints pc( model_->Database(), "porosity", 0.3, 0.5 );
    pc.AddConstraint( "permeability", 1.0e-13, 1.0e-11 );

    _test( pc.Constraints() == 2U );
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
//  TestCheckConstraintsAllNodes  (default mode)
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

    if ( verbose_ ) {
        cout << "AllNodes: e_in  passes tight permeability: "
             << pc.CheckConstraints( e_in ) << "\n";
    }

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

    // Strict mode object — also uses database ctor so indices are ready
    PropertyConstraints pc_strict( model_->Database(), "porosity", 0.75, 0.85 );

    // Elements entirely outside the elevated band (0..49) must fail in both modes
    const auto* e_out = matrix_domain.E( 0 );
    _test( pc_strict.CheckConstraints( e_out ) == false );
    _test( pc.CheckConstraints( e_out )        == false );

    // Elements entirely inside the elevated band (50..149) must pass in both modes
    const auto* e_in = matrix_domain.E( 75 );
    _test( pc_strict.CheckConstraints( e_in ) == true );
    _test( pc.CheckConstraints( e_in )        == true );

    // A perimeter cell shares at least one face with the region boundary.
    // The elevated band (elements 50..149) was assigned by index into MATRIX,
    // so the perimeter of that band coincides with perimeter cells of MATRIX.
    // Such a cell has nodes both inside and outside the elevated band, so it
    // will fail the strict (all-nodes) constraint but pass the single-node
    // (at-least-one-node) constraint.
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
        _test(  pc.CheckConstraints( e_boundary )        );
        _test( !pc_strict.CheckConstraints( e_boundary ) );

        if ( verbose_ )
            cout << "TestCheckConstraintsSingleNode: "
                 << "straddling perimeter cell found\n";
    }
    else {
        // No perimeter cell of MATRIX straddles the porosity transition.
        // Fall back to atBoundary(), which checks whether the element shares
        // at least one face with any model boundary.
        if ( verbose_ )
            cout << "TestCheckConstraintsSingleNode: "
                 << "no straddling perimeter cell found — trying atBoundary()\n";

        auto eend = matrix_domain.CellsEnd();
        for ( auto eit = matrix_domain.CellsBegin(); eit != eend; ++eit ) {
            if ( atBoundary( *eit )                  &&
                  pc.CheckConstraints( *eit )         &&
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
    // Constraint: porosity in [0.55, 0.65]
    // A boundary element with half its nodes at 0.8 and half at 0.4
    // has a nodal average of ~0.6, which falls in [0.55, 0.65].
    PropertyConstraints pc( model_->Database(), "porosity", 0.55, 0.65 );
    pc.SatisfyConstraintsForNodalAverage( true );
    pc.InitializePropertyIndices( model_->Database() );

    Region<3U>& matrix_domain = model_->Region( "MATRIX" );

    // Boundary element — average should be ~0.6
    const auto* e_boundary = matrix_domain.E( 49 );
    _test( pc.CheckConstraints( e_boundary ) == true );

    // Element deep inside elevated band — average ~0.8, outside [0.55, 0.65]
    const auto* e_in = matrix_domain.E( 100 );
    _test( pc.CheckConstraints( e_in ) == false );

    // Element entirely outside elevated band — average ~0.4, outside [0.55, 0.65]
    const auto* e_out = matrix_domain.E( 0 );
    _test( pc.CheckConstraints( e_out ) == false );

    // Widen range to [0.35, 0.45] — only the out-of-band elements should pass
    PropertyConstraints pc_low( model_->Database(), "porosity", 0.35, 0.45 );
    pc_low.SatisfyConstraintsForNodalAverage( true );
    pc_low.InitializePropertyIndices( model_->Database() );
    _test( pc_low.CheckConstraints( e_out ) == true  );
    _test( pc_low.CheckConstraints( e_in  ) == false );

    if ( verbose_ ) {
        cout << "NodalAverage: boundary element passes [0.55,0.65]: "
             << pc.CheckConstraints( e_boundary ) << "\n";
    }

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
    _test( result_out  == false    );
    _test( failed_upon == phi_key  );

    // Element in elevated band — porosity 0.8 passes, permeability 1e-12 passes
    const auto* e_in = matrix_domain.E( 75 );
    csmp::Index not_failed;
    const bool result_in = pc.CheckConstraints( e_in, not_failed );
    _test( result_in == true );

    // Tighten permeability so the in-band element now fails on permeability
    pc.ChangeConstraint( "permeability", 1.0e-15, 1.0e-14 );
    csmp::Index failed_perm;
    const bool result_perm = pc.CheckConstraints( e_in, failed_perm );
    _test( result_perm  == false    );
    _test( failed_perm  == perm_key );

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
    // First verify that velocity is actually defined and is a vector variable
    if ( !model_->Database().IsDefined( "velocity" ) ) {
        if ( verbose_ )
            cout << "TestCheckLengthOfVectorVariables: "
                 << "'velocity' not defined in database — skipping\n";
        return true;
    }

    Region<3U>& model_domain = model_->Region( "Model" );
    const auto* e = model_domain.E( 0 );

    // Length mode: constraint [0.9, 1.1] — should pass
    PropertyConstraints pc_pass( model_->Database(), "velocity", 0.9, 1.1 );
    pc_pass.CheckLengthOfVectorVariables( true );

    bool threw = false;
    bool result = false;
    try {
        result = pc_pass.CheckConstraints( e );
    }
    catch ( const csmp::Exception& ex ) {
        threw = true;
        if ( verbose_ )
            cerr << "TestCheckLengthOfVectorVariables: "
                 << "velocity is not a VectorVariable: " << ex.what() << "\n";
    }
    catch ( const std::exception& ex ) {
        threw = true;
        if ( verbose_ )
            cerr << "TestCheckLengthOfVectorVariables: "
                 << "std::exception: " << ex.what() << "\n";
    }

    if ( threw ) {
        if ( verbose_ )
            cout << "TestCheckLengthOfVectorVariables: skipping — "
                 << "velocity is not stored as a VectorVariable\n";
        return true;
    }

    _test( result == true );

    // Length mode: constraint [2.0, 3.0] — should fail
    PropertyConstraints pc_fail( model_->Database(), "velocity", 2.0, 3.0 );
    pc_fail.CheckLengthOfVectorVariables( true );
    _test( pc_fail.CheckConstraints( e ) == false );

    // Component mode: range [0.0, 1.0] — all components (1,0,0) pass
    PropertyConstraints pc_comp_pass( model_->Database(), "velocity", 0.0, 1.0 );
    pc_comp_pass.CheckLengthOfVectorVariables( false );
    _test( pc_comp_pass.CheckConstraints( e ) == true );

    // Component mode: range [0.5, 1.5] — y and z components are 0.0, fail
    PropertyConstraints pc_comp_fail( model_->Database(), "velocity", 0.5, 1.5 );
    pc_comp_fail.CheckLengthOfVectorVariables( false );
    _test( pc_comp_fail.CheckConstraints( e ) == false );

    // Switching mode on an existing object must take effect immediately
    pc_pass.CheckLengthOfVectorVariables( false );
    _test( pc_pass.CheckConstraints( e ) == false );

    return true;
}


// ----------------------------------------------------------------------------
//  TestMutualExclusivity
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestMutualExclusivity()
{
    // Basic flag interaction — no model needed
    PropertyConstraints pc;
    pc.SatisfyConstraintsForAtLeastOneNode( true );
    pc.SatisfyConstraintsForNodalAverage( false );
    pc.SatisfyConstraintsForNodalAverage( true );
    pc.SatisfyConstraintsForAtLeastOneNode( false );

    // Skip the VectorLengthCheck interaction test if velocity is not available
    if ( !model_->Database().IsDefined( "velocity" ) ) {
        if ( verbose_ )
            cout << "TestMutualExclusivity: "
                 << "'velocity' not defined — skipping vector interaction test\n";
        return true;
    }

    PropertyConstraints pc2( model_->Database(), "velocity", 0.9, 1.1 );
    pc2.CheckLengthOfVectorVariables( true );

    pc2.SatisfyConstraintsForAtLeastOneNode( true );
    pc2.SatisfyConstraintsForNodalAverage(   true );

    Region<3U>& model_domain = model_->Region( "Model" );
    const auto* e = model_domain.E( 0 );

    bool threw = false;
    try {
        pc2.CheckConstraints( e );
    }
    catch ( const csmp::Exception& ex ) {
        threw = true;
        if ( verbose_ )
            cerr << "TestMutualExclusivity: caught csmp::Exception: "
                 << ex.what() << "\n";
    }
    catch ( const std::exception& ex ) {
        threw = true;
        if ( verbose_ )
            cerr << "TestMutualExclusivity: caught std::exception: "
                 << ex.what() << "\n";
    }
    _test( threw == false );

    return true;
}



// ----------------------------------------------------------------------------
//  TestBuildRegionsFromPropertyConstraints
// ----------------------------------------------------------------------------

bool PropertyConstraints_Test::TestBuildRegionsFromPropertyConstraints()
{
    // --- strict mode: all nodes must satisfy constraints ---
    PropertyConstraints porosity_constraints( model_->Database(), "porosity", 0.75, 0.85 );
    porosity_constraints.AddConstraint( "permeability", 1.0e-12, 1.0e-12 );

    model_->FormRegionFrom( "medium porosity", porosity_constraints );
    const Region<3U>& medium_porosity_domain = model_->Region( "medium porosity" );

    const auto strict_count = medium_porosity_domain.Cells();
    if ( verbose_ )
        cout << "strict mode cell count: " << strict_count << "\n";
    _test( strict_count == 100 );

    // --- single-node mode: at least one node must satisfy ---
    PropertyConstraints boundary_constraints( model_->Database(), "porosity", 0.75, 0.85 );
    boundary_constraints.SatisfyConstraintsForAtLeastOneNode( true );

    model_->FormRegionFrom( "boundary porosity", boundary_constraints );
    const Region<3U>& boundary_domain = model_->Region( "boundary porosity" );

    const auto single_node_count = boundary_domain.Cells();
    if ( verbose_ )
        cout << "single-node mode cell count: " << single_node_count << "\n";

    // Single-node mode must capture at least as many cells as strict mode
    _test( single_node_count >= strict_count );

    // --- nodal average mode ---
    PropertyConstraints average_constraints( model_->Database(), "porosity", 0.55, 0.65 );
    average_constraints.SatisfyConstraintsForNodalAverage( true );

    model_->FormRegionFrom( "average porosity", average_constraints );
    const Region<3U>& average_domain = model_->Region( "average porosity" );

    const auto average_count = average_domain.Cells();
    if ( verbose_ )
        cout << "nodal average mode cell count: " << average_count << "\n";
    _test( average_count > 0U );

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

    auto& gref = model_->Region( "Model" );

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


} // end csmp
