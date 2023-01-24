#include "PDE_Integrator_Test.h"
#include "PDE_Integrator.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "Attorney.h"

#include "VSet.h"
#include "ModelTopology.h"
#include "vsetMakers.h"

#include "LinearSolver.h"
#include "Model.h"
#include "Region.h"

using namespace std;

namespace csmp {

/**

Initialises Element integral contribution to global solution matrix with a fixed value as input in constructor.

This method sets @a value_for_matrix_ on the diagonal and -@a value_for_matrix_
off the diagonal. If the TestOperand is a vector or array then off-diagonal terms are
set separately for each component of the vector (no coupling).

@note transformNodeIndexVector( uint32_t dim, const csmp::Index&, vector<size_t>& ); in MathOperatorRHS changes matrices according to the degrees of freedom of the solution variables
*/
template<uint32_t dim,class CELL>
void LHS_FixedValueMatrix<dim,CELL>::ComputeContribution( const CELL& e )
  {
    assert( MathOperatorLHS<dim>::MaterialOperandType() == SCALAR );
    assert( MathOperatorLHS<dim>::MTRL[0].Rows() == dim );
    assert( MathOperatorLHS<dim>::MTRL[0].Cols() == dim );
    
    const uint32_t nodes{ e.Nodes() };
    // note: element matrix is not necessarily square, e.g., when coupling matrix blocks for scalar with vector degrees of freedom
    MathOperatorLHS<dim>::LHS.Resize( nodes * this->BasicOperandDataDepth(), nodes * this->TestOperandDataDepth() );
    MathOperatorLHS<dim>::LHS.Zero();
    
    // SCALAR variables: fixed values are directly written to element matrix
    if ( MathOperatorLHS<dim>::TestOperandType() == SCALAR ) {
         MathOperatorLHS<dim>::LHS = -value_label_for_matrix_entry_; // off-diagonal is negative
         MathOperatorLHS<dim>::LHS.AssignToDiagonal( value_label_for_matrix_entry_ ); // diagonal is positive
         // this->LHS.Out(0);
         return;
      }
      
    // VECTOR and other variable types: are expanded to achieve an ordering comp1, comp2, comp3... in LH matrix and RH vector
    // (the different components are tagged with integers, e.g., if scalar=3, the component becomes 3x; x={0..dim}
    else {
        for ( uint32_t i{0U}; i < nodes; ++i )
          {
             for ( uint32_t j{0U}; j < nodes; ++j )
               {
                  // create positive diagonal elements
                  if ( i == j ) {
                      for ( uint32_t k{0U}; k<this->BasicOperandDataDepth(); ++k )
                        for ( uint32_t l{0U}; l<this->TestOperandDataDepth(); ++l )
                          this->LHS( i * this->BasicOperandDataDepth() + k,
                                     j * this->TestOperandDataDepth() + l ) = value_label_for_matrix_entry_ + l + 10;
                    }
                  // create negative off-diagonal elements
                  else {
                       for ( uint32_t k{0U}; k<this->BasicOperandDataDepth(); ++k )
                         for ( uint32_t l{0U}; l<this->TestOperandDataDepth(); ++l )
                           this->LHS( i * this->BasicOperandDataDepth() + k,
                                      j * this->TestOperandDataDepth() + l ) = -(value_label_for_matrix_entry_ + l + 10);
                    }
               }
          }
        // this->LHS.Out(0);
        return;
      }
      
    throw csmp::Exception( ERROR, "LHS_FixedValueMatrix<dim,CELL>::ComputeContribution", "TestOperandType not handled yet." );
    
} // end ComputeContribution




/** Computes contribution from value input in constructor.
*/
template<uint32_t dim,class CELL>
void RHS_FixedValueMatrix<dim,CELL>::ComputeContribution( const CELL& e )
 {
    uint32_t nodes{ e.Nodes() };
    this->RHS.resize( nodes * this->TestOperandDataDepth() );

    for ( uint32_t i{0U}; i<this->RHS.size(); i++ )
        this->RHS[i] = value_label_for_vector_entry_;
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
template<uint32_t dim>
void PDE_Integrator_Test::TestMatrix(
        const PDE_Integrator_Attorney<dim>& attorney,
        const MeshManager<dim>& mesh,
        const set<BOX_BOUNDARY>& dirich,
        const VARIABLE_TYPE variable_type )
  {
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
                  nb_dof  *= dim;
                  var_size = dim;
                }
              break;
            case SCALAR:
              break;
            default:
              throw csmp::Exception( ERROR, "TestMatrix", "Tensor, Array and FlaggedArray variables not handled yet");
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
        
    } // end TestMatrix



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
template<uint32_t dim>
void PDE_Integrator_Test::TestRHSVector(  const PDE_Integrator_Attorney<dim>& attorney,
                                          const MeshManager<dim>& mesh,
                                          const double val,
                                          const set<BOX_BOUNDARY>& dirich,
                                          const VARIABLE_TYPE variable_type )
  {
        size_t ni = 0;
        for( auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n ) {
            if( dirich.find( n->AtBoundary() ) != dirich.end() ) {
                continue;
            }
            _test( attorney.rh_[ni++] == val * n->Parents() );
            if( variable_type == SCALAR ) {
                continue;
            }
            _test( attorney.rh_[ni++] == val * n->Parents() );
        }
        
  } // end TestVector





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

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<2U,Element>  pde_integrator(solver);

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

    TestMatrix<2U>( attorney, mesh, dirich, SCALAR );
    TestRHSVector<2U>( attorney, mesh, val, dirich, SCALAR );

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

    CSMP_DEFAULT_LINEAR_SOLVER solver;
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

    TestMatrix<2U>( attorney, mesh, dirich, SCALAR );
    TestRHSVector<2U>( attorney, mesh, val, dirich, SCALAR );

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

    CSMP_DEFAULT_LINEAR_SOLVER solver;
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

    TestMatrix<2U>( attorney, mesh, dirich, VECTOR );
    TestRHSVector<2U>( attorney, mesh, val, dirich, VECTOR );

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

    CSMP_DEFAULT_LINEAR_SOLVER solver;
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

    TestMatrix<2U>( attorney, mesh, dirich, VECTOR );
    TestRHSVector<2U>( attorney, mesh, val, dirich, VECTOR );

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
//    TestAssembly();

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
    TestAssemblyTwoScalarVariablesNoDirichlet( true /* debug */ );
    TestAssemblyScalarAndVectorVariableNoDirichlet( true /* debug */ );
    
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



/**

Testing a  case with two scalar variables that are coupled together.
Since the pde operators are stored in alphabetical order, 'concentration' (scalar 2, valued 2),
gets accumulated into the first blocked matrix, and 'fluid pressure' (scalar 1, valued 1)
gets accumlated into the second block.

SKM - printing matrix in integer format to better see pattern.

test_CreateSimplestPolyElement2DModel() - 2 triangles + 1 quadrilateral

*/
void PDE_Integrator_Test::TestAssemblyTwoScalarVariablesNoDirichlet( bool debug )
 {
    cout <<"\nPDE_Integrator_Test::TestAssemblyTwoScalarVariablesNoDirichlet: starting test."<< endl;
    Reset();

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<2U,Element>  pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<2U> attorney( pde_integrator );

    // Create pde-operators for a system with 2 coupled scalar solution variables (placed on the node)
    const double val1 = 1.;
    // scalar 1: 'fluid pressure'
    LHS_FixedValueMatrix<2U>  lhs1( model_->Database(), "permeability", "fluid pressure", "fluid pressure", val1 );
    RHS_FixedValueMatrix<2U>  rhs1( model_->Database(), "fluid volume source", "fluid pressure", val1 );
    // assign them to pde integrator
    attorney.Add(&lhs1);
    attorney.Add(&rhs1);

    // scalar 2: 'concentration'
    const double val2 = 2.;
    LHS_FixedValueMatrix<2U>  lhs2( model_->Database(), "diffusivity", "concentration", "concentration", val2 );
    RHS_FixedValueMatrix<2U>  rhs2( model_->Database(), "fluid volume source", "concentration", val2 );
    attorney.Add(&lhs2);
    attorney.Add(&rhs2);
    
    // cross-coupling terms (lower diagonal of matrix)
    // scalar 1: 'fluid pressure' with 'concentration'
    LHS_FixedValueMatrix<2U>  lhs12( model_->Database(), "element number", "concentration", "fluid pressure", val1 );
    attorney.Add(&lhs12);

    // cross-coupling terms (upper diagonal of matrix)
    // scalar 1: 'concentration' with 'fluid pressure'
    LHS_FixedValueMatrix<2U>  lhs21( model_->Database(), "element number", "fluid pressure", "concentration", val1 );
    attorney.Add(&lhs21);

    
    // TESTING
    Region<2U>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 4 );
    _test( attorney.rhs_operators_.size() == 2 );

    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<2U>& mesh = model_->Mesh();
    const size_t nb_nodes = mesh.Nodes() * 2U;
    _test( attorney.G_.Rows() == nb_nodes );
    _test( attorney.G_.Cols() == nb_nodes );

    model_region.RenumberNodes();
    attorney.Accumulate( model_region );

    set<BOX_BOUNDARY>  dirich; // which boundary flags do the solution variables have?
    TestMatrix<2U>( attorney, mesh, dirich, SCALAR );
    // Note: can only test first block of the matrix which is occupied by the 'concentration' values
    TestRHSVector<2U>( attorney, mesh, val2, dirich, SCALAR );

    if ( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 /* print zero decimal places */ );
     }
    cout <<"\nPDE_Integrator_Test::TestAssemblyTwoScalarVariablesNoDirichlet: finished test."<< endl;
    
} // end TestAssemblyTwoScalarVariablesNoDirichlet







/**

Testing a  case with a  scalar variable and a vector variable coupled together.

SKM - printing matrix in integer format to better see pattern.

test_CreateSimplestPolyElement2DModel() - 2 triangles + 1 quadrilateral

*/
void PDE_Integrator_Test::TestAssemblyScalarAndVectorVariableNoDirichlet( bool debug )
 {
    cout <<"\nPDE_Integrator_Test::TestAssemblyScalarAndVectorVariableNoDirichlet: starting test."<< endl;
    Reset();
    const uint32_t model_dimensions{2U};

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<model_dimensions,Element>  pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<model_dimensions> attorney( pde_integrator );

    // Create pde-operators for a system with 2 coupled scalar solution variables (placed on the node)
    const double val1 = 1.;
    // scalar 1: 'fluid pressure'
    LHS_FixedValueMatrix<model_dimensions>  lhs1( model_->Database(), "permeability", "nodal velocity", "nodal velocity", val1 );
    RHS_FixedValueMatrix<model_dimensions>  rhs1( model_->Database(), "element number", "nodal velocity", val1 );
    // assign them to pde integrator
    attorney.Add(&lhs1);
    attorney.Add(&rhs1);

    // scalar 2: 'concentration'
    const double val2 = 2.;
    LHS_FixedValueMatrix<model_dimensions>  lhs2( model_->Database(), "diffusivity", "concentration", "concentration", val2 );
    RHS_FixedValueMatrix<model_dimensions>  rhs2( model_->Database(), "fluid volume source", "concentration", val2 );
    attorney.Add(&lhs2);
    attorney.Add(&rhs2);
    
    // cross-coupling terms assembled into upper diagonal of matrix
    // ------------------------------------------------------------
    // coupling basic operand 'concentration' rows with test operand 'nodal velocity' columns
    LHS_FixedValueMatrix<model_dimensions>  lhs12( model_->Database(), "element number", "concentration", "nodal velocity", val1 );
    attorney.Add(&lhs12);

    // cross-coupling terms assembled into lower diagonal of matrix
    // ------------------------------------------------------------
    // coupling 'nodal velocity' with 'concentration'
    LHS_FixedValueMatrix<model_dimensions>  lhs21( model_->Database(), "element number", "nodal velocity", "concentration", val1 );
    attorney.Add(&lhs21);

    
    // TESTING
    const Region<model_dimensions>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 4 );
    _test( attorney.rhs_operators_.size() == 2 );

    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<model_dimensions>& mesh = model_->Mesh();
    //                                  nodes * ( 1-for-scalars + dim-for-vectors)
    const size_t n_degrees_of_freedom = mesh.Nodes() + mesh.Nodes() * model_dimensions;
    _test( attorney.G_.Rows() == n_degrees_of_freedom );
    _test( attorney.G_.Cols() == n_degrees_of_freedom );

    model_region.RenumberNodes();
    attorney.Accumulate( model_region );

/*
    set<BOX_BOUNDARY>  dirich; // which boundary flags do the solution variables have?
    TestMatrix<2U>( attorney, mesh, dirich, SCALAR );
    // Note: can only test first block of the matrix which is occupied by the 'concentration' values
    TestRHSVector<2U>( attorney, mesh, val2, dirich, SCALAR );
*/
    if ( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 /* print zero decimal places */ );
     }
    cout <<"\nPDE_Integrator_Test::TestAssemblyScalarAndVectorVariableNoDirichlet: finished test."<< endl;
    
} // end TestAssemblyScalarAndVectorVariableNoDirichlet







#if 0

void PDE_Integrator_Test::TestTwoScalarVariables() {
    Region<2U>& region = model_->Region("Model");
    region.RenumberNodes();
    
    CSMP_DEFAULT_LINEAR_SOLVER  solver;

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

    const vector<size_t>& DOF_indexes = pde_test_->GetDOFIndex();
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
    map<size_t, double> result;
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) != DIRICH) {
        result[(*nIter)->Idx()] = (*nIter)->Read(pressureKey);
      } 
    }

    vector<double>* valid_x = pde_test_->GetX();
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

