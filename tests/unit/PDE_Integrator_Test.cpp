#include "PDE_Integrator_Test.h"
#include "PDE_Integrator.h"
#include "Attorney.h"

#include "VSet.h"
#include "ModelTopology.h"
#include "vsetMakers.h"

#include "GaussJordan_Solver.h"
#include "Model.h"
#include "Region.h"

// TODO: distinguish this UNIT TEST from an INTEGRATION TEST
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_dNT_op_dV.h"
#include "VelocityAndVolumeFlux.h"
    // thermal equation
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"

using namespace std;

namespace csmp {

// Attorney design pattern gives access to protected / private member variables
// and methods
template<uint32_t dim>
class PDE_Integrator_Attorney : public Attorney<class PDE_Integrator<dim>> {
public:
    using Attorney<PDE_Integrator<dim>>::Attorney; // Inherit Attorney constructor
    using PDE_Integrator<dim>::Accumulate;
    using PDE_Integrator<dim>::EstablishMatrixSetup;
    using PDE_Integrator<dim>::lhs_operators_;
    using PDE_Integrator<dim>::rhs_operators_;
    using PDE_Integrator<dim>::G_;
    using PDE_Integrator<dim>::rh_;
};



/**

Computes contribution from fixed value input in constructor.

This method sets @a value_for_matrix_ on the diagonal and -@a value_for_matrix_
off the diagonal. If the TestOperand is a vector, then off-diagonal terms are
set separately for each component of the vector (no coupling).
*/
template<uint32_t dim,class CELL>
void LHS_FixedValueMatrix<dim,CELL>::ComputeContribution( const CELL& e ) {
    const size_t nodes = e.Nodes();
    size_t size = nodes;
    if( MathOperatorLHS<dim>::TestOperandType() == VECTOR ) {
        size *= 2;
    }
    MathOperatorLHS<dim>::LHS.Resize( size, size );
    MathOperatorLHS<dim>::LHS.Zero();
    if( MathOperatorLHS<dim>::TestOperandType() == SCALAR ) {
        MathOperatorLHS<dim>::LHS = -value_for_matrix_; // off-diagonal is negative
        MathOperatorLHS<dim>::LHS.AssignToDiagonal( value_for_matrix_ ); // diagonal is positive
        return;
    } else if( MathOperatorLHS<dim>::TestOperandType() == VECTOR ) {
        for( size_t i = 0; i < nodes; ++i ) {
            for( size_t j = 0; j < nodes; ++j ) {
                if( i == j ) {
                    this->LHS(i*2, i*2) = value_for_matrix_;
                    this->LHS(i*2+1, i*2+1) = value_for_matrix_;
                } else {
                    this->LHS(i*2, j*2) = -value_for_matrix_;
                    this->LHS(i*2+1, j*2+1) = -value_for_matrix_;
                }
            }
        }
        this->LHS.Out();
        return;
    }
    throw csmp::Exception( WARNING, "LHS_FixedValueMatrix<dim,CELL>::ComputeContribution", "TestOperandType not handled yet." );
}




/**

Computes contribution from value input in constructor.

*/
template<uint32_t dim,class CELL>
void RHS_FixedValueMatrix<dim,CELL>::ComputeContribution( const CELL& e ) {
    size_t size = e.Nodes();
    if( MathOperatorRHS<dim>::TestOperandType() == VECTOR ) {
        size *= 2;
    }
    this->RHS.resize( size );
    for ( auto i{0U}; i<size; i++ )
        this->RHS[i] = value_for_vector_;
}




/**

Constructor of PDE_Integrator_Test without input Model.

This creates a simple test model and stores it in @a model_.

*/
PDE_Integrator_Test::PDE_Integrator_Test()
: model_( nullptr ), delete_model_( true ) {
    // 0. creates test model with 4 elements and 6 Face objects for the box boundaries
    VSet<2U> vset;
    ModelTopology topo = test_CreateSimplestPolyElement2DModel( vset );
    const bool treat_regions_as_boundaries{false};
    model_ = new Model<2U>( topo, vset, "CSMP-1phase-variables.txt", treat_regions_as_boundaries );
}

/**

Constructor of PDE_Integrator_Test with input Model.

*/
PDE_Integrator_Test::PDE_Integrator_Test( Model<2U>& model )
: model_( &model ), delete_model_( false ) {
}

PDE_Integrator_Test::~PDE_Integrator_Test() {
    if( delete_model_ ) {
        delete model_;
    }
}



// PDE_Integrator test functions
namespace {

/**

Generic function to test the matrix accumulated in PDE_Integrator from
LHS_FixedValueMatrix
@param[in] attorney Attorney to the PDE_Integrator being tested
@param[in] mesh Mesh of the model used in this test
@param[in] dirich BOX_BOUNDARY codes on which Dirichlet conditions are set
@param[in] variable_type Type of Test variable (currently SCALAR or VECTOR only)

This function tests that
  1. diagonal terms are positive
  2. off-diagonal terms for nodes belonging to the same Element are negative
  3. all other terms are zero
 
*/
    template<uint32_t dim> void testMatrix(
        const PDE_Integrator_Attorney<dim>& attorney,
        const MeshManager<dim>& mesh,
        const set<BOX_BOUNDARY>& dirich,
        const VARIABLE_TYPE variable_type
    ) {
        const size_t nb_nodes = mesh.Nodes();

        // Get indices of Dirichlet Nodes, for which
        // corresponding line and column are missing in matrix
        set<size_t> dirich_nodes;
        for( auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n ) {
            if( dirich.find( n->AtBoundary() ) != dirich.end() ) {
                dirich_nodes.insert( n->Idx() );
            }
        }

        size_t nb_dof = nb_nodes;
        size_t var_size = 1;
        switch( variable_type ) {
            case VECTOR: {
                nb_dof *= 2;
                var_size = 2;
                break;
            }
            case SCALAR:
            default:
                break;
        }

        // Matrix should have positive diagonal terms
        size_t n = 0;
        vector<size_t> mapping( nb_dof );
        for( size_t i = 0; i < nb_nodes; i++ ) {
            if( dirich_nodes.find( i ) != dirich_nodes.end() ) {
                for( size_t v = 0; v < var_size; ++v ) {
                    mapping[i*var_size + v] = NULL_IDX;
                }
                continue;
            }
            for( size_t v = 0; v < var_size; ++v ) {
                mapping[i*var_size + v] = n;
                _test( attorney.G_.At( n, n ) > 0. );
                n++;
            }
        }

        // Matrix should have negative off-diagonal terms for Nodes belonging to the
        // same Element
        set<pair<size_t, size_t>> negative;
        for( auto el = mesh.ElementsBegin(); el != mesh.ElementsEnd(); ++el ) {
            for( auto ni = el->NodesBegin(); ni != el->NodesEnd(); ++ni ) {
                if( dirich_nodes.find( (*ni)->Idx() ) != dirich_nodes.end() ) {
                    continue;
                }
                const size_t mi = mapping[(*ni)->Idx()*var_size];
                assert( mi != NULL_IDX );
                for( auto nj = ni + 1; nj != el->NodesEnd(); ++nj ) {
                    if( dirich_nodes.find( (*nj)->Idx() ) != dirich_nodes.end() ) {
                        continue;
                    }
                    const size_t mj = mapping[(*nj)->Idx()*var_size];
                    assert( mj != NULL_IDX );
                    _test( attorney.G_.At( mi, mj ) < 0. );
                    _test( attorney.G_.At( mj, mi ) < 0. );
                    negative.insert( {mi, mj} );
                    negative.insert( {mj, mi} );
                    if( variable_type == SCALAR ) {
                        continue;
                    }
                    const size_t mip = mi + 1;
                    const size_t mjp = mj + 1;
                    _test( attorney.G_.At( mip, mjp ) < 0. );
                    _test( attorney.G_.At( mjp, mip ) < 0. );
                    negative.insert( {mip, mjp} );
                    negative.insert( {mjp, mip} );
                }
            }
        }

        // All other terms should be zero
        for( size_t i = 0; i < nb_nodes; ++i ) {
            if( dirich_nodes.find( i ) != dirich_nodes.end() ) {
                continue;
            }
            const size_t mi = mapping[i*var_size];
            assert( mi != NULL_IDX );
            for( size_t j = 0; j < nb_nodes; ++j ) {
                if( i == j || dirich_nodes.find( j ) != dirich_nodes.end() ) {
                    continue;
                }
                const size_t mj = mapping[j*var_size];
                assert( mj != NULL_IDX );
                if( negative.find( {mi, mj} ) != negative.end() ) {
                    continue;
                }
                _test( attorney.G_( mi, mj ) == 0. );
                if( variable_type == SCALAR ) {
                    continue;
                }
                const size_t mip = mi + 1;
                const size_t mjp = mj + 1;
                if( negative.find( {mip, mjp} ) != negative.end() ) {
                    continue;
                }
                _test( attorney.G_( mip, mjp ) == 0. );
            }
        }
    }



/**

Generic function to test the RHS vector accumulated in PDE_Integrator from
RHS_FixedValueMatrix
@param[in] attorney Attorney to the PDE_Integrator being tested
@param[in] mesh Mesh of the model used in this test
@param[in] dirich BOX_BOUNDARY codes on which Dirichlet conditions are set
@param[in] variable_type Type of Test variable (currently SCALAR or VECTOR only)

This function tests that RHS value for a given non-Dirichlet node is the value
input in the constructor times the number of parent elements for this node.
 
*/
    template<uint32_t dim> void testRHSVector(
        const PDE_Integrator_Attorney<dim>& attorney,
        const MeshManager<dim>& mesh,
        const double val,
        const set<BOX_BOUNDARY>& dirich,
        const VARIABLE_TYPE variable_type
    ) {
        size_t ni = 0;
        for( auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n ) {
            if( dirich.find( n->AtBoundary() ) != dirich.end() ) {
                continue;
            }
            _test( attorney.rh_[ni++] == val*n->Parents() );
            if( variable_type == SCALAR ) {
                continue;
            }
            _test( attorney.rh_[ni++] == val*n->Parents() );
        }
    }
} // end anonymous namespace





/**

Tests assembly of the system from LHS and RHS operators, in various
configurations

*/
void PDE_Integrator_Test::TestAssembly()
 {
    // 0. check model
    assert( model_ != nullptr );
    
    const bool debug{true};

    // --------------------------------------------------------------------------
    // 1. setting up PDE_Integrator for simple case with a single scalar variable
    // --------------------------------------------------------------------------
    TestAssemblySingleScalarNoDirichlet( debug );

    // --------------------------------------------------------------------------
    // 2. setting up PDE_Integrator for simple case with a single scalar variable
    //    and a Dirichlet boundary condition
    // --------------------------------------------------------------------------
    TestAssemblySingleScalarDirichlet( debug );

    // --------------------------------------------------------------------------
    // 3. setting up PDE_Integrator for case with a single vector solution variable
    // ----------------------------------------------------------------------------
    TestAssemblySingleVectorNoDirichlet( debug );

    // --------------------------------------------------------------------------
    // 4. setting up PDE_Integrator for case with a single vector solution variable
    //    and a Dirichlet boundary condition
    // ----------------------------------------------------------------------------
    TestAssemblySingleVectorDirichlet( debug );

    // 5. setting up PDE_Integrator for coupled system of 2 scalars
    // ----------------------------------------------------------------------------
    

    // 6. setting up PDE_Integrator for coupled system of 1 scalar + 1 vector variable
    // -------------------------------------------------------------------------------
    
    cout <<"\nPDE_Integrator_Test::TestAssembly: finished test."<< endl;
    
 } // end TestAssembly




/**

Testing a simple case with a single scalar variable

*/
void PDE_Integrator_Test::TestAssemblySingleScalarNoDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarNoDirichlet: starting test."<< endl;
    Reset();

    // No Dirichlet conditions
    set<BOX_BOUNDARY> dirich;

    GaussJordan_Solver         solver;
    PDE_Integrator<2U,Element> pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<2U> attorney( pde_integrator );

    // Create pde-operators
    const double val = 1.;
    LHS_FixedValueMatrix<2U>  lhs( model_->Database(), "permeability", "fluid pressure", "fluid pressure", val );
    RHS_FixedValueMatrix<2U>  rhs( model_->Database(), "fluid volume source", "fluid pressure", val );

    // assign them to pde integrator
    attorney.Add(&lhs);
    attorney.Add(&rhs);
    
    Region<2U>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 1 );
    _test( attorney.rhs_operators_.size() == 1 );

    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<2U>& mesh = model_->Mesh();
    const size_t nb_nodes = mesh.Nodes();
    _test( attorney.G_.Rows() == nb_nodes );
    _test( attorney.G_.Cols() == nb_nodes );

    model_region.RenumberNodes();
    attorney.Accumulate( model_region );

    testMatrix<2U>( attorney, mesh, dirich, SCALAR );
    testRHSVector<2U>( attorney, mesh, val, dirich, SCALAR );

    if( debug ) {
        attorney.Out();
        attorney.OutputGlobals();
    }
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarNoDirichlet: finished test."<< endl;
}




/**

Testing a simple case with a single scalar variable and a Dirichlet boundary
condition

*/
void PDE_Integrator_Test::TestAssemblySingleScalarDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarDirichlet: starting test."<< endl;
    Reset();

    // Dirichlet condition
    model_->InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.) );
    set<BOX_BOUNDARY> dirich;
    dirich.insert( RIGHT );
    dirich.insert( CNR2 );
    dirich.insert( CNR6 );
    dirich.insert( CNR3 );
    dirich.insert( CNR8 );

    const MeshManager<2U>& mesh = model_->Mesh();

    // Get indices of Dirichlet Nodes, for which
    // corresponding line and column are missing in matrix
    set<size_t> dirich_nodes;
    for( auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n ) {
        if( dirich.find( n->AtBoundary() ) != dirich.end() ) {
            dirich_nodes.insert( n->Idx() );
        }
    }

    GaussJordan_Solver         solver;
    PDE_Integrator<2U,Element> pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<2U> attorney( pde_integrator );

    // Create pde-operators
    const double val = 1.;
    LHS_FixedValueMatrix<2U>  lhs( model_->Database(), "permeability", "fluid pressure", "fluid pressure", val );
    RHS_FixedValueMatrix<2U>  rhs( model_->Database(), "fluid volume source", "fluid pressure", val );

    // assign them to pde integrator
    attorney.Add(&lhs);
    attorney.Add(&rhs);
    
    Region<2U>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 1 );
    _test( attorney.rhs_operators_.size() == 1 );

    // compare matrix with expected matrix,
    // taking into account Dirichlet conditions
    const size_t nb_nodes = mesh.Nodes() - dirich_nodes.size();
    _test( attorney.G_.Rows() == nb_nodes );
    _test( attorney.G_.Cols() == nb_nodes );

    model_region.RenumberNodes();
    attorney.Accumulate( model_region );

    testMatrix<2U>( attorney, mesh, dirich, SCALAR );
    testRHSVector<2U>( attorney, mesh, val, dirich, SCALAR );

    if( debug ) {
        attorney.Out();
        attorney.OutputGlobals();
    }
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarDirichlet: finished test."<< endl;
}




/**

Test case with a single vector solution variable, no boundary conditions and no
coupling between the variable components

*/
void PDE_Integrator_Test::TestAssemblySingleVectorNoDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleVectorNoDirichlet: starting test."<< endl;
    Reset();

    // No boundary condition
    set<BOX_BOUNDARY> dirich;

    GaussJordan_Solver         solver;
    PDE_Integrator<2U,Element> pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<2U> attorney( pde_integrator );

    // Create pde-operators
    const double val = 1.;
    LHS_FixedValueMatrix<2U>  lhs( model_->Database(), "Young's modulus", "displacement", "displacement", val );
    RHS_FixedValueMatrix<2U>  rhs( model_->Database(), "force", "displacement", val );

    // assign them to pde integrator
    attorney.Add(&lhs);
    attorney.Add(&rhs);
    
    Region<2U>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 1 );
    _test( attorney.rhs_operators_.size() == 1 );

    // compare matrix with expected matrix
    const MeshManager<2U>& mesh = model_->Mesh();
    const size_t nb_nodes2 = mesh.Nodes()*2;
    _test( attorney.G_.Rows() == nb_nodes2 );
    _test( attorney.G_.Cols() == nb_nodes2 );

    model_region.RenumberNodes();
    attorney.Accumulate( model_region );

    testMatrix<2U>( attorney, mesh, dirich, VECTOR );
    testRHSVector<2U>( attorney, mesh, val, dirich, VECTOR );

    if( debug ) {
        attorney.Out();
        attorney.OutputGlobals();
    }
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleVectorNoDirichlet: finished test."<< endl;
}



/**

Test case with a single vector solution variable, a boundary condition and no
coupling between the variable components

*/
void PDE_Integrator_Test::TestAssemblySingleVectorDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleVectorDirichlet: starting test."<< endl;
    Reset();

    // Dirichlet condition
    model_->InputBoundaryValue( RIGHT, "displacement", makeVector( DIRICH, DIRICH, 0.1, 0.5 ) );
    set<BOX_BOUNDARY> dirich;
    dirich.insert( RIGHT );
    dirich.insert( CNR2 );
    dirich.insert( CNR6 );
    dirich.insert( CNR3 );
    dirich.insert( CNR8 );

    const MeshManager<2U>& mesh = model_->Mesh();

    // Get indices of Dirichlet Nodes, for which
    // corresponding line and column are missing in matrix
    set<size_t> dirich_nodes;
    for( auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n ) {
        if( dirich.find( n->AtBoundary() ) != dirich.end() ) {
            dirich_nodes.insert( n->Idx() );
        }
    }

    GaussJordan_Solver         solver;
    PDE_Integrator<2U,Element> pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<2U> attorney( pde_integrator );

    // Create pde-operators
    const double val = 1.;
    LHS_FixedValueMatrix<2U>  lhs( model_->Database(), "Young's modulus", "displacement", "displacement", val );
    RHS_FixedValueMatrix<2U>  rhs( model_->Database(), "force", "displacement", val );

    // assign them to pde integrator
    attorney.Add(&lhs);
    attorney.Add(&rhs);
    
    Region<2U>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 1 );
    _test( attorney.rhs_operators_.size() == 1 );

    // compare matrix with expected matrix
    const size_t nb_nodes = mesh.Nodes() - dirich_nodes.size();
    const size_t nb_nodes2 = nb_nodes*2;
    _test( attorney.G_.Rows() == nb_nodes2 );
    _test( attorney.G_.Cols() == nb_nodes2 );

    model_region.RenumberNodes();
    attorney.Accumulate( model_region );

    testMatrix<2U>( attorney, mesh, dirich, VECTOR );
    testRHSVector<2U>( attorney, mesh, val, dirich, VECTOR );

    if( debug ) {
        attorney.Out();
        attorney.OutputGlobals();
    }
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleVectorDirichlet: finished test."<< endl;
}



    // Create pde-operators
//    NumIntegral_dNT_op_dN_dV<2U> pressureLHS(model_->Database(), "permeability", "fluid pressure", "fluid pressure");
//    NumIntegral_NT_op_N_dV<2U> sourceVolume(model_->Database(), "fluid volume source", "fluid pressure");
//    pde_integrator.Add(&pressureLHS);
//    pde_integrator.Add(&sourceVolume);





/**

Reset stored model between two tests

*/
void PDE_Integrator_Test::Reset()
  {
    // Initialise properties and remove boundary conditions
    model_->InputPropertyValue("fluid pressure", makeScalar(PLAIN, 1.0));
    model_->InputPropertyValue("permeability", makeScalar(PLAIN, 1.));
    model_->InputPropertyValue("fluid volume source", makeScalar(PLAIN, 1.));

    // Find or create additional required properties
    if( !model_->Database().IsDefined( "force" ) ) {
        model_->CreateProperty( "force", "G", "Pa", VECTOR, NODE );
    }
    model_->InputPropertyValue( "force", makeVector( ANY, ANY, 0.0, 0.0 ) );
    if( !model_->Database().IsDefined( "displacement" ) ) {
        model_->CreateProperty( "displacement", "displ", "m", VECTOR, NODE );
    }
    model_->InputPropertyValue( "displacement", makeVector( ANY, ANY, 0.1, 0.5 ) );
    if( !model_->Database().IsDefined( "Young's modulus" ) ) {
        model_->CreateProperty( "Young's modulus", "E", "Pa", SCALAR, NODE );
    }
    model_->InputPropertyValue( "Young's modulus", makeScalar(PLAIN, 4.5e+10) );
  }




  void PDE_Integrator_Test::run()
  {
    //=======================================
    // test single variable
    //=======================================
    TestAssembly();

    // test scalar variable
    //TestSingleVariable();
    //TestOutputSingleVariable();
    // test vector variable

    // test array variable

    // test flagged array variable

    // test tensor variable




    //=======================================
    // test multiple single variable
    //=======================================

    // test 2 scalar variables
    //TestTwoScalarVariables();

    // test scalar variable and vector variable

    // test vector variable and vector variable

  }

/*
  void PDE_Integrator_Test::TestSingleVariable() {
    // 0. prepare
    Region<2U>& region = model_->Region("Model");
    region.RenumberNodes();

    GaussJordan_Solver solver;

    // Create pde-operators
    NumIntegral_dNT_op_dN_dV<2U> pressureLHS(model_->Database(),
      "permeability",
      "fluid pressure",
      "fluid pressure");

    NumIntegral_NT_op_N_dV<2U> sourceVolume(model_->Database(),
      "fluid volume source", "fluid pressure");

    PointSource_rhsop<2U> sourcePoint(model_->Database(),
      "nodal fluid point source", "fluid pressure");
    
    NumIntegral_dNT_op_dV<2U> gravityTerm(model_->Database(),
      "gravity vector",
      "fluid pressure");

    //fluid_velocity = new VelocityAndVolumeFlux<2U, Element<2U>>(*model, "total mobility", "porosity", "fluid pressure", true, "velocity");

    NumIntegral_dNT_op_dN_dV<2U> temperatureLHS(model_->Database(),
      "permeability",
      "temperature",
      "temperature");

    NumIntegral_NT_op_N_dV<2U> sourceHeat(model_->Database(),
      "thermal volume source", "temperature");


    pde_reference_ = new PDE_Integrator<2U,Region>(solver);
    pde_test_      = new PDE_Integrator<2U,Region>(solver);

    pde_reference_->Add(&pressureLHS);
    pde_reference_->Add(&sourceVolume);
    pde_reference_->Add(&sourcePoint);
    pde_reference_->Add(&gravityTerm);
    pde_test_->Add(&pressureLHS);
    pde_test_->Add(&sourceVolume);
    pde_test_->Add(&sourcePoint);
    pde_test_->Add(&gravityTerm);

    const vector<size_t>& DOF_indexes = pde_test_->DOF_indexes_;

    // 1. test matrix establish and enumerate DOFs
    pde_reference_->EstablishMatrixSetup(region);
    pde_test_->EstablishMatrixSetup(region);

    size_t pressureDirchletDOFs(0);
    Index pressureKey = model_->Database().StorageKey("fluid pressure");
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) == DIRICH) {
        pressureDirchletDOFs += 1;
      }
    }

    _test((pde_reference_->GetRH()->size()) == (pde_test_->GetRH()->size()));

    pde_test_->EnumerateAndFixMatrixSize(region);
    _test(pde_reference_->GetRH()->size() == DOF_indexes.size());
    _test(pde_reference_->GetRH()->size() == pde_test_->GetRH()->size() + pressureDirchletDOFs);
    _test(pde_reference_->GetG()->Rows() == pde_test_->GetG()->Rows() + pressureDirchletDOFs);
    _test(pde_reference_->GetG()->Cols() == pde_test_->GetG()->Cols() + pressureDirchletDOFs);

    // 2. test accumulate
    pde_reference_->Accumulate(region);
    pde_test_->Accumulate(region);

    // 3. test AssignInitialConditions & LateAccumulate
    pde_reference_->TimeIncrement(0.1);
    pde_test_->TimeIncrement(0.1);
    pde_reference_->AssignInitialConditions(region);
    pde_test_->AssignInitialConditions(region);
	_test(pde_reference_->Transient() == true);
	_test(pde_test_->Transient() == true);
	pde_reference_->LateAccumulate(region);
	pde_test_->LateAccumulate(region);

    for (size_t i(0); i < DOF_indexes.size(); ++i) {
      if (DOF_indexes[i] != NULL_IDX) {
        // _test conductance matrix
        for (size_t j(0); j < DOF_indexes.size(); ++j) {
          if (DOF_indexes[j] != NULL_IDX) {
            _test(pde_test_->GetG()->At(DOF_indexes[i], DOF_indexes[j]) == pde_reference_->GetG()->At(i, j));
          }
        }
        // check load vector
        _test((*pde_test_->GetRH())[DOF_indexes[i]] == (*pde_reference_->GetRH())[i]);
      }
    }

    delete pde_reference_;
    delete pde_test;

  } // end method
*/



#if 0

void PDE_Integrator_Test::TestTwoScalarVariables() {
    Region<2U>& region = model_->Region("Model");
    region.RenumberNodes();
    
    GaussJordan_Solver solver;

    pde_reference_ = new PDE_Integrator_<2U, Region>(solver);
    pde_test = new PDE_Integrator<2U, Region>(solver);

    pde_reference_->Add(pressureLHS);
    pde_reference_->Add(sourceVolume);
    pde_reference_->Add(sourcePoint);
    pde_reference_->Add(gravityTerm);
    pde_reference_->Add(temperatureLHS);
    pde_reference_->Add(sourceHeat);

    pde_test_->Add(pressureLHS);
    pde_test_->Add(sourceVolume);
    pde_test_->Add(sourcePoint);
    pde_test_->Add(gravityTerm);
    pde_test_->Add(temperatureLHS);
    pde_test_->Add(sourceHeat);

    const std::vector<size_t>& DOF_indexes = pde_test_->GetDOFIndex();
    // 1. test matrix establish and enumerate DOFs
    pde_reference_->EstablishMatrixSetup(region);
    pde_test_->EstablishMatrixSetupTest(region);
    size_t dirchletDOFs(0);
    Index pressureKey = model_->Database().StorageKey("fluid pressure");
    Index temperatureKey = model_->Database().StorageKey("temperature");
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) == DIRICH) {
        dirchletDOFs += 1;
      }
      if ((*nIter)->Status(temperatureKey) == DIRICH) {
        dirchletDOFs += 1;
      }
    }

    _test((pde_reference_->GetRH()->size()) == (pde_test_->GetRH()->size()));
    pde_test_->EnumerateAndFixMatrixSize(region);
    _test(pde_reference_->GetRH()->size() == DOF_indexes.size());
    _test(pde_reference_->GetRH()->size() == pde_test_->GetRH()->size() + dirchletDOFs);
    _test(pde_reference_->GetG()->Rows() == pde_test_->GetG()->Rows() + dirchletDOFs);
    _test(pde_reference_->GetG()->Cols() == pde_test_->GetG()->Cols() + dirchletDOFs);

    // 2. test accumulate
    pde_reference_->Accumulate(region);
    pde_test_->AccumulateTest(region);

	// 3. test AssignInitialConditions & LateAccumulate
	pde_reference_->TimeIncrement(0.1);
	pde_test_->TimeIncrement(0.1);
	pde_reference_->AssignInitialConditions(region);
	pde_test_->AssignInitialConditionsTest(region);
	_test(pde_reference_->Transient() == true);
	_test(pde_test_->Transient() == true);
	pde_reference_->LateAccumulate(region);
	pde_test_->LateAccumulateTest(region);
    for (size_t i(0); i < DOF_indexes.size(); ++i) {
      if (DOF_indexes[i] != NULL_IDX) {
        // _test conductance matrix
        for (size_t j(0); j < DOF_indexes.size(); ++j) {
          if (DOF_indexes[j] != NULL_IDX) {
            _test(pde_test_->GetG()->At(DOF_indexes[i], DOF_indexes[j]) == pde_reference_->GetG()->At(i, j));
          }
        }
        // check load vector
        _test((*pde_test_->GetRH())[DOF_indexes[i]] == (*pde_reference_->GetRH())[i]);
      }
    }

    delete pde_reference_;
    delete pde_test;
  }


  void PDE_Integrator_Test::TestOutputSingleVariable() {
    Region<2U>& region = model_->Region("Model");
    region.RenumberNodes();
    Reset();

    Index pressureKey = model_->Database().StorageKey("fluid pressure");
    Index pressureValKey = model_->Database().StorageKey("fluid pressure previous");
    
    GaussJordan_Solver solver;

    pde_test = new PDE_Integrator_UoM_Mock<2U, Region>(solver);
    pde_test_->Add(pressureLHS);
    pde_test_->Add(sourceVolume);
    pde_test_->EstablishMatrixSetupTest(region);
    pde_test_->EnumerateAndFixMatrixSize(region);
    std::map<size_t, double> result;
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) != DIRICH) {
        result[(*nIter)->Idx()] = (*nIter)->Read(pressureKey);
      } 
    }

    std::vector<double>* valid_x = pde_test_->GetX();
    size_t idx(0);
    for (auto& it : result) {
      (*valid_x)[idx] = it.second;
      ++idx;
    }

    pde_test_->OutputResultsTest(region);
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      _test((*nIter)->Read(pressureKey) == (*nIter)->Read(pressureValKey));     
    }
    
    delete pde_test;
  }

#endif

} // end namespace csmp

