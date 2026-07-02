#include "PDE_Integrator_Test.h"
#include "PDE_Integrator.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "Attorney.h"
#include "compareFloats.h"

#include "NumIntegral_BT_D_B_dV.h"
#include "NumJacobianIntegral_BT_m_N_dV.h"
#include "NumJacobianIntegral_N_mT_B_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_PT_op_dV.h"

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
    const uint32_t nodes{ e.Nodes() };
    const uint32_t row_depth = this->TestOperandDataDepth(); // Rows (Test)
    const uint32_t col_depth = this->BasicOperandDataDepth(); // Cols (Basic)

    // note: element matrix is not necessarily square, e.g., when coupling matrix blocks for scalar with vector degrees of freedom
    MathOperatorLHS<dim>::LHS.Resize( nodes * row_depth, nodes * col_depth );
    MathOperatorLHS<dim>::LHS.Zero();
    
    // Universal loop handles all coupling permutations (Scalar, Vector, Mixed)
    for ( uint32_t i{0U}; i < nodes; ++i ) {
        for ( uint32_t j{0U}; j < nodes; ++j ) {
            for ( uint32_t k_row{0U}; k_row < row_depth; ++k_row ) {
                for ( uint32_t k_col{0U}; k_col < col_depth; ++k_col ) {
                    
                    if (k_row == k_col) {
                        const double val = value_label_for_matrix_entry_;
                        this->LHS( i * row_depth + k_row, j * col_depth + k_col ) = 
                            (i == j) ? val : -val;
                    }
                    // Else is handled by LHS.Zero()
                }
            }
        }
    }
    
} // end ComputeContribution




/** Computes contribution from value input in constructor.
*/
template<uint32_t dim,class CELL>
void RHS_FixedValueMatrix<dim,CELL>::ComputeContribution( const CELL& e )
{
    const uint32_t nodes{ e.Nodes() };
    const uint32_t dofs_per_node = this->TestOperandDataDepth();
    
    this->RHS.resize( nodes * dofs_per_node );

    for ( uint32_t i{0U}; i < nodes; ++i )
    {
        for ( uint32_t k{0U}; k < dofs_per_node; ++k )
        {
            // Differentiate the value by the component index (k)
            // Example: If base is 1.0, DOF 0 gets 1.0, DOF 1 gets 2.0
            this->RHS[i * dofs_per_node + k] = value_label_for_vector_entry_ + static_cast<double>(k);
        }
    }
}



/**

Constructor of PDE_Integrator_Test without input Model.

This creates a simple test model and stores it in @a model_.

*/
PDE_Integrator_Test::PDE_Integrator_Test()
: model_( nullptr ), delete_model_( true ) {
    // 0. creates test model with 4 elements and 6 Face objects for the box boundaries
    VSet<2U> vset;
    ModelTopology topo = create_SimplePolyElement2DModel( vset );
    const bool treat_regions_as_boundaries{false};
    model_ = new Model<2U>( topo, vset, "CSMP-1phase-variables.txt", treat_regions_as_boundaries );
}



/** Constructor of PDE_Integrator_Test with input Model.
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
  
  @attention This test only works for single solution variables
 
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator_Test::TestMatrix(
    const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
    const MeshManager<dim>&       mesh,
    const set<BOX_BOUNDARY>& dirich,
    const VARIABLE_TYPE           variable_type )
{
   // systems of equations
   if ( attorney.test_operands_.size() > 1 ) {
        ErrorHandler::Instance().Note( WARNING, "PDE_Integrator_Test::TestMatrix", "test works only for global matrices with a single solution variable");
        return;
     }

    const size_t nb_nodes = mesh.Nodes();

    // 1. Identify Dirichlet Nodes
    set<size_t> dirich_nodes;
    for (auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n) {
        if (dirich.find(n->AtBoundary()) != dirich.end()) {
            dirich_nodes.insert(n->Idx());
        }
    }

    // 2. Determine variable size dynamically
    uint32_t var_size{1U};
    switch (variable_type) {
        case SCALAR:        var_size = 1U; break;
        case VECTOR:        var_size = dim; break;
        case TENSOR:        var_size = dim * dim; break;
        case ARRAY:
        case FLAGGEDARRAY:  var_size = 0; /* Handle specifically if needed */ break;
        default:
            throw csmp::Exception(ERROR, "TestMatrix", "Variable type not supported");
    }

    const size_t nb_dof = nb_nodes * var_size;

    // 3. Build Active Equation Mapping Matrix
    size_t active_eq_counter = 0;
    vector<size_t> mapping(nb_dof, NULL_IDX);
    
    for (size_t i = 0; i < nb_nodes; ++i) {
        if (dirich_nodes.find(i) != dirich_nodes.end()) continue;
        
        for (size_t v = 0; v < var_size; ++v) {
            mapping[i * var_size + v] = active_eq_counter++;
        }
    }

    // 4. Build Node Adjacency Map (Which nodes share an element)
    set<pair<size_t, size_t>> adjacent_nodes;
    for (auto el = mesh.ElementsBegin(); el != mesh.ElementsEnd(); ++el) {
        for (auto ni = el->NodesBegin(); ni != el->NodesEnd(); ++ni) {
            for (auto nj = el->NodesBegin(); nj != el->NodesEnd(); ++nj) {
                if ((*ni)->Idx() != (*nj)->Idx()) {
                    adjacent_nodes.insert({(*ni)->Idx(), (*nj)->Idx()});
                }
            }
        }
    }

    // 5. Comprehensive Matrix Structural Verification Loop
    for (size_t i = 0; i < nb_nodes; ++i) {
        if (dirich_nodes.find(i) != dirich_nodes.end()) continue;

        for (size_t j = 0; j < nb_nodes; ++j) {
            if (dirich_nodes.find(j) != dirich_nodes.end()) continue;

            const bool share_element = (adjacent_nodes.find({i, j}) != adjacent_nodes.end());

            // Check every component permutation between node i and node j
            for (size_t vi = 0; vi < var_size; ++vi) {
                const size_t mi = mapping[i * var_size + vi];

                for (size_t vj = 0; vj < var_size; ++vj) {
                    const size_t mj = mapping[j * var_size + vj];

                    if (i == j) {
                        // CASE A: Same Node
                        if (vi == vj) {
                            // Main Diagonal Terms: must be strictly positive
                            _test(attorney.G_.At(mi, mj) > 0.);
                        } else {
                            // Intra-node component cross-coupling: must be zero (uncoupled)
                            _test(attorney.G_.At(mi, mj) == 0.);
                        }
                    } 
                    else if (share_element) {
                        // CASE B: Different Nodes Sharing an Element
                        if (vi == vj) {
                            // Corresponding component interaction: must be negative
                            _test(attorney.G_.At(mi, mj) < 0.);
                        } else {
                            // Inter-node cross-component coupling: must be zero (uncoupled)
                            _test(attorney.G_.At(mi, mj) == 0.);
                        }
                    } 
                    else {
                        // CASE C: Unconnected Nodes
                        // Any interaction between independent components must be zero
                        _test(attorney.G_.At(mi, mj) == 0.);
                    }
                }
            }
        }
    }
} // end TestMatrix





template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator_Test::TestMultiVariableMatrix(
                            const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
                            const MeshManager<dim>&                                    mesh,
                            const vector<size_t>&                                      DOF_indexes )
{
    // -----------------------------------------------------------------------
    // Build per-variable Dirichlet node sets
    // A node is Dirichlet for variable V if its boundary flag matches
    // the Dirichlet boundary for V.
    // -----------------------------------------------------------------------
    // For this test:
    //   "fluid pressure"  Dirichlet on RIGHT
    //   "concentration"   Dirichlet on LEFT
    // We determine this per-variable from the actual node Status flags.
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // Build node adjacency map (nodes sharing at least one element)
    // -----------------------------------------------------------------------
    set<pair<size_t,size_t>> adjacent_nodes;
    for ( auto el=mesh.ElementsBegin(); el!=mesh.ElementsEnd(); ++el )
        for ( auto ni=el->NodesBegin(); ni!=el->NodesEnd(); ++ni )
            for ( auto nj=el->NodesBegin(); nj!=el->NodesEnd(); ++nj )
                if ( (*ni)->Idx() != (*nj)->Idx() )
                    adjacent_nodes.insert({ (*ni)->Idx(), (*nj)->Idx() });

    // -----------------------------------------------------------------------
    // Helper: get component count for a variable type
    // -----------------------------------------------------------------------
    auto comp_count = [&]( const csmp::Index& key ) -> uint32_t {
        switch ( key.type ) {
            case SCALAR:     return 1U;
            case VECTOR:     return dim;
            case TENSOR:     return dim * dim;
            case ARRAY:
            case FLAGGEDARRAY: return key.dataDepth;
            default:         return 1U;
        }
    };

    // -----------------------------------------------------------------------
    // Helper: flat DOF position for node idx, variable key, component vi
    // -----------------------------------------------------------------------
    auto flat_pos = [&]( size_t node_idx,
                          const csmp::Index& key,
                          size_t offset,
                          uint32_t vi ) -> size_t {
        if ( key.type == SCALAR )
            return node_idx + offset;
        return node_idx * comp_count(key) + vi + offset;
    };

    // -----------------------------------------------------------------------
    // Helper: is this (node, variable, component) Dirichlet?
    // -----------------------------------------------------------------------
    auto is_dirich = [&]( const auto* node,
                           const csmp::Index& key,
                           uint32_t vi ) -> bool {
        if ( key.type == SCALAR )
            return node->Status(key) == DIRICH;
        return node->Status(key, vi) == DIRICH;
    };

    // -----------------------------------------------------------------------
    // Check whether any cross-coupling LHS operators exist
    // -----------------------------------------------------------------------
    bool has_cross_coupling = false;
    for ( const auto& [name, lhs_op] : attorney.lhs_operators_ )
    {
        const auto& basic_key = lhs_op->BasicOperandKey();
        const auto& test_key  = lhs_op->TestOperandKey();
        // Two keys are the same variable if they have the same index
        if ( basic_key.index != test_key.index )
        {
            has_cross_coupling = true;
            break;
        }
    }

    // -----------------------------------------------------------------------
    // Main verification loop
    // -----------------------------------------------------------------------
    for ( auto ni=mesh.NodesBegin(); ni!=mesh.NodesEnd(); ++ni )
    {
        const size_t i = (*ni).Idx();

        for ( auto nj=mesh.NodesBegin(); nj!=mesh.NodesEnd(); ++nj )
        {
            const size_t j = (*nj).Idx();
            const bool share_element =
                ( i == j ) ||
                ( adjacent_nodes.find({i,j}) != adjacent_nodes.end() );

            for ( const auto& [param_i, offset_i] : attorney.test_operands_ )
            {
                const auto& key_i = param_i.key;
                const uint32_t nc_i = comp_count(key_i);

                for ( uint32_t vi=0; vi<nc_i; ++vi )
                {
                    // Skip Dirichlet row DOFs
                    if ( is_dirich( &(*ni), key_i, vi) ) continue;

                    const size_t pos_i = flat_pos(i, key_i, offset_i, vi);
                    const size_t mi    = DOF_indexes[pos_i];
                    if ( mi == NULL_IDX ) continue;

                    for ( const auto& [param_j, offset_j] : attorney.test_operands_ )
                    {
                        const auto& key_j = param_j.key;
                        const uint32_t nc_j = comp_count(key_j);

                        for ( uint32_t vj=0; vj<nc_j; ++vj )
                        {
                            // Skip Dirichlet column DOFs
                            if ( is_dirich( &(*nj), key_j, vj) ) continue;

                            const size_t pos_j = flat_pos(j, key_j, offset_j, vj);
                            const size_t mj    = DOF_indexes[pos_j];
                            if ( mj == NULL_IDX ) continue;

                            // Same variable and same component?
                            const bool same_var =
                                ( key_i == key_j ) && ( vi == vj );

                            const double G_ij = attorney.G_.At(mi, mj);

                            if ( !share_element )
                            {
                                // --- CASE C: Unconnected nodes ---
                                // No element connects i and j — must be zero
                                // (entry should not even be in sparsity pattern)
                                _test( G_ij == 0.0 );
                            }
                            else if ( i == j )
                            {
                                // --- CASE A: Same node ---
                                if ( same_var )
                                {
                                    // Diagonal: must be strictly positive
                                    _test( G_ij > 0.0 );
                                }
                                else if ( has_cross_coupling )
                                {
                                    // Cross-variable same-node: non-zero only if cross-coupling exists
                                    _test( G_ij != 0.0 );
                                }
                                // else: no cross-coupling — legitimately zero, no assertion
                            }
                            else
                            {
                                // --- CASE B: Different nodes sharing an element ---
                                if ( same_var )
                                {
                                    // Off-diagonal same-variable: must be strictly negative
                                    _test( G_ij < 0.0 );
                                }
                                else if ( has_cross_coupling )
                                {
                                    // Cross-variable off-diagonal: positive only if cross-coupling exists
                                    _test( G_ij > 0.0 );
                                }
                                // else: no cross-coupling — legitimately zero, no assertion
                            }
                        }
                    }
                }
            }
        }
    }

} // end TestMultiVariableMatrix




/**
Generic function to test the (single solution variable) RHS vector accumulated in PDE_Integrator from
RHS_FixedValueMatrix
@param[in] attorney Attorney to the PDE_Integrator being tested
@param[in] mesh Mesh of the model used in this test
@param[in] val Base value input in the constructor
@param[in] dirich BOX_BOUNDARY codes on which Dirichlet conditions are set

This function tests that RHS value for a given non-Dirichlet node is the value
input in the constructor (offset by the component index for vectors) times the 
number of parent elements for this node.

  @attention This test only works for single solution variables and a prescribed material op value of 'val'

*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator_Test::TestRHSVector(
    const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
    const MeshManager<dim>&                                    mesh,
    const double                                               val,
    const set<BOX_BOUNDARY>&                                   dirich )
{
    // Single-variable test only
    if ( attorney.test_operands_.size() > 1 )
        throw csmp::Exception( FATAL_ERROR,
            "PDE_Integrator_Test::TestRHSVector",
            "test works only for a single solution variable" );

    // -----------------------------------------------------------------------
    // Get the single test operand key and offset
    // -----------------------------------------------------------------------
    const auto& [param, offset] = *attorney.test_operands_.begin();
    const csmp::Index& key = param.key;

    const size_t var_size = [&]() -> size_t {
        switch ( key.type ) {
            case SCALAR:     return 1U;
            case VECTOR:     return dim;
            case TENSOR:     return dim * dim;
            case ARRAY:
            case FLAGGEDARRAY: return key.dataDepth;
            default:         return 1U;
        }
    }();

    constexpr double tol = 1.0e-10;

    // -----------------------------------------------------------------------
    // For each active (non-Dirichlet) DOF, compute expected RHS value.
    //
    // rh_[i] = source_term[i] + pivot_contribution[i]
    //
    // source_term[i]:
    //   RHS_FixedValueMatrix contributes val per parent element per component.
    //   => source_term = val * node->Parents()
    //
    // pivot_contribution[i]:
    //   For each Dirichlet column j sharing an element with free row i:
    //     rh_[i] -= K_ij * x_j_prescribed
    //   For LHS_FixedValueMatrix:
    //     K_ij = +val (diagonal, same node)  or  -val (off-diagonal)
    //   Prescribed value x_j = val (from InputBoundaryValue in the test)
    //
    // For the standard single-variable test with all Dirichlet values = val:
    //   Each Dirichlet neighbour contributes: -(-val)*val = +val^2
    //   The Dirichlet self-entry (if node is its own Dirichlet):
    //     not applicable since Dirichlet nodes are skipped
    // -----------------------------------------------------------------------

    // Build Dirichlet node set
    set<size_t> dirich_node_ids;
    for ( auto n=mesh.NodesBegin(); n!=mesh.NodesEnd(); ++n )
        if ( dirich.find(n->AtBoundary()) != dirich.end() )
            dirich_node_ids.insert(n->Idx());

    // Build per-node count of Dirichlet neighbours (nodes sharing an element)
    map<size_t, size_t> dirich_neighbour_count;
    for ( auto el=mesh.ElementsBegin(); el!=mesh.ElementsEnd(); ++el )
    {
        for ( auto ni=el->NodesBegin(); ni!=el->NodesEnd(); ++ni )
        {
            if ( dirich_node_ids.count((*ni)->Idx()) ) continue; // skip Dirichlet rows

            for ( auto nj=el->NodesBegin(); nj!=el->NodesEnd(); ++nj )
            {
                if ( !dirich_node_ids.count((*nj)->Idx()) ) continue; // only Dirichlet cols
                if ( (*ni)->Idx() == (*nj)->Idx() ) continue;         // skip self

                dirich_neighbour_count[(*ni)->Idx()]++;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Validate each active DOF
    // -----------------------------------------------------------------------
    for ( auto n=mesh.NodesBegin(); n!=mesh.NodesEnd(); ++n )
    {
        const size_t node_idx = n->Idx();

        // Skip Dirichlet nodes
        if ( dirich_node_ids.count(node_idx) ) continue;

        for ( size_t v=0; v<var_size; ++v )
        {
            // Compute reduced DOF index
            const size_t pos = ( key.type == SCALAR )
                             ? node_idx + offset
                             : node_idx * var_size + v + offset;

            const size_t reduced_idx = attorney.DOF_indexes_[pos];
            if ( reduced_idx == NULL_IDX ) continue;

            // Source term: val per parent element
            const double source = val * static_cast<double>( n->Parents() );

            // Pivot contribution: each Dirichlet neighbour contributes
            // -K_ij * x_j = -(-val) * val = +val^2
            // (K_ij = -val for off-diagonal LHS_FixedValueMatrix,
            //  x_j = val prescribed)
            const size_t n_dirich_neighbours =
                dirich_neighbour_count.count(node_idx)
                ? dirich_neighbour_count.at(node_idx)
                : 0U;

            const double pivot = static_cast<double>(n_dirich_neighbours)
                               * val * val;

            const double expected = source + pivot;
            const double actual   = attorney.rh_[reduced_idx];
            const double diff     = std::abs( actual - expected );

            if ( diff >= tol )
            {
                std::cout << "RHS mismatch at node " << node_idx
                          << " component " << v
                          << " reduced_idx " << reduced_idx
                          << ": actual="   << actual
                          << " expected="  << expected
                          << " diff="      << diff << "\n";
                _test( diff < tol );
            }
        }
    }

} // end TestRHSVector(val)




// same as above but for values read by the pde_operators
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator_Test::TestRHSVector(
    const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
    const MeshManager<dim>&                                    mesh )
{
    if ( attorney.test_operands_.size() > 1 )
        throw csmp::Exception( FATAL_ERROR,
            "PDE_Integrator_Test::TestRHSVector",
            "test works only for a single solution variable" );

    const auto& [param, offset] = *attorney.test_operands_.begin();
    const csmp::Index& key = param.key;

    const size_t var_size = [&]() -> size_t {
        switch ( key.type ) {
            case SCALAR:       return 1U;
            case VECTOR:       return static_cast<size_t>(dim);
            case TENSOR:       return static_cast<size_t>(dim * dim);
            case ARRAY:
            case FLAGGEDARRAY: return key.dataDepth;
            default:           return 1U;
        }
    }();

    double val = 1.0;
    for ( const auto& [name, lhs_op] : attorney.lhs_operators_ )
        { val = lhs_op->MultiplyBy(); break; }

    constexpr double tol = 1.0e-10;

    // Per-component Dirichlet status — uses node Status flags directly
    // This is independent of renumbering
    auto is_dirich = [&]( const auto* node, size_t comp ) -> bool {
        if ( key.type == SCALAR )
            return node->Status(key) == DIRICH;
        return node->Status(key, static_cast<uint32_t>(comp)) == DIRICH;
    };

    // Read prescribed value for component comp from node
    auto prescribed = [&]( const auto* node, size_t comp ) -> double {
        if ( key.type == SCALAR )
            return node->Read(key);
        VectorVariable<dim> vc;
        node->Read(key, vc);
        return vc.Component(static_cast<uint32_t>(comp));
    };

    // -----------------------------------------------------------------------
    // For each free DOF, compute expected RHS:
    //
    //   expected = val * n_parents                    (source term)
    //            + val * sum(x_prescribed[comp])      (pivot from Dirichlet neighbours)
    //
    // Dirichlet neighbours are identified by their Status flags — not by
    // membership in a dirich_nodes set — so renumbering does not affect this.
    // -----------------------------------------------------------------------
    for ( auto n=mesh.NodesBegin(); n!=mesh.NodesEnd(); ++n )
    {
        const size_t node_idx = n->Idx();

        for ( size_t v=0; v<var_size; ++v )
        {
            if ( is_dirich( &(*n), v) ) continue;

            const size_t pos = ( key.type == SCALAR )
                             ? node_idx + offset
                             : node_idx * var_size + v + offset;

            if ( pos >= attorney.DOF_indexes_.size() ) continue;
            const size_t ridx = attorney.DOF_indexes_[pos];
            if ( ridx == NULL_IDX ) continue;

            // Source term
            double expected = ( val + static_cast<double>(v) ) * static_cast<double>( n->Parents() );
                
            // Pivot: iterate all parent elements
            for ( auto el=mesh.ElementsBegin(); el!=mesh.ElementsEnd(); ++el )
            {
                // Does element contain this node?
                bool has_node = false;
                for ( auto ni=el->NodesBegin(); ni!=el->NodesEnd(); ++ni )
                    if ( (*ni)->Idx() == node_idx )
                        { has_node = true; break; }
                if ( !has_node ) continue;

                // Sum Dirichlet contributions for component v
                for ( auto nj=el->NodesBegin(); nj!=el->NodesEnd(); ++nj )
                {
                    if ( (*nj)->Idx() == node_idx ) continue;
                    if ( !is_dirich(*nj, v) ) continue;

                    // K_ij = -val (off-diagonal LHS_FixedValueMatrix)
                    // pivot contribution = -K_ij * x_j = +val * x_j
                    expected += val * prescribed(*nj, v);
                }
            }

            const double actual = attorney.rh_[ridx];
            const double diff   = std::abs(actual - expected);

            if ( diff >= tol )
            {
                std::cout << "RHS mismatch at node " << node_idx
                          << " comp " << v
                          << " ridx " << ridx
                          << ": actual=" << actual
                          << " expected=" << expected
                          << " diff=" << diff << "\n";
                _test( diff < tol );
            }
        }
    }

} // end TestRHSVector






template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator_Test::TestMultiVariableRHSVector(
    const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
    const map<Parameter, size_t>&                             test_operands,
    const vector<size_t>&                                     DOF_indexes,
    const ModelSubDomain<dim, CELLTYPE>&                      model_region )
{
    // -----------------------------------------------------------------------
    // The RHS must satisfy:  A * x_free + A_FD * x_dirich = b_source
    // => b = b_source - A_FD * x_dirich
    //
    // For LHS_FixedValueMatrix and RHS_FixedValueMatrix, we verify two
    // simple structural properties that hold regardless of the specific
    // operator values:
    //
    // 1. SIGN: For val > 0, each free DOF's RHS entry must be positive
    //    (source term dominates for nodes far from Dirichlet boundaries,
    //     and Dirichlet values are positive in these tests)
    //
    // 2. CONSISTENCY: A * rh_ should be computable without NaN/Inf
    //    (verifies the assembled system is numerically sound)
    //
    // 3. SYMMETRY of pivot: rh_[i] is modified by Dirichlet neighbours
    //    consistently with the assembled matrix off-diagonal entries:
    //    rh_[i] += -G(i,j_dirich) * x_j  for each Dirichlet col j
    // -----------------------------------------------------------------------

    const size_t n_free = attorney.rh_.size();

    // --- Check 1: No NaN or Inf in RHS ---
    bool no_nan = true;
    for ( size_t i=0; i<n_free; ++i )
        if ( !std::isfinite(attorney.rh_[i]) )
        {
            std::cout << "Non-finite RHS at idx " << i
                      << ": " << attorney.rh_[i] << "\n";
            no_nan = false;
        }
    _test( no_nan );

    // --- Check 2: Matrix-vector product is computable ---
    std::vector<double> Ax;
    attorney.G_.MultiplyWith( attorney.rh_, Ax );
    _test( Ax.size() == n_free );
    bool Ax_finite = true;
    for ( size_t i=0; i<n_free; ++i )
        if ( !std::isfinite(Ax[i]) )
            { Ax_finite = false; break; }
    _test( Ax_finite );

    // --- Check 3: Pivot consistency ---
    // For each free DOF i, verify:
    //   rh_[i] = source_i + sum_j( -G(i,j) * x_j )  for Dirichlet j
    // where source_i > 0 for RHS_FixedValueMatrix with val > 0
    //
    // We verify this indirectly: the assembled matrix G already encodes
    // the connectivity. For each free row i and each column j that IS
    // in the sparsity pattern but corresponds to a Dirichlet DOF,
    // the pivot was applied. We check the sign is consistent.
    //
    // Specifically: if G(i,j) < 0 (off-diagonal) and x_j > 0 (prescribed),
    // then the pivot contribution +(-G(i,j))*x_j > 0 increases rh_[i].
    // The total rh_[i] must therefore be >= source_i > 0.
    bool rhs_positive = true;
    for ( size_t i=0; i<n_free; ++i )
        if ( attorney.rh_[i] <= 0.0 )
        {
            std::cout << "Non-positive RHS at idx " << i
                      << ": " << attorney.rh_[i] << "\n";
            rhs_positive = false;
        }
    _test( rhs_positive );

} // end TestMultiVariableRHSVector




/**
    Is the matrix symmetric and positive definite (SPD)?
    
    Symmetry ($A_{ij} = A_{ji}$): This is straightforward to test in a Compressed Row format by searching for the transposed column index within the corresponding row.
    
    Positive Definite ($x^T A x > 0$): Rigorously proving a matrix is Positive Definite requires computing a Cholesky decomposition or finding all its eigenvalues,
    which is too heavy for a quick unit test. Instead, the standard industry practice for FEM unit tests is a two-part heuristic:
    
    Strictly Positive Diagonals: Every diagonal entry $A_{ii}$ must exist and be $> 0$.
    Stochastic Energy Test: Generate a pseudo-random vector $x$, compute $y = Ax$, and check that the scalar "energy" $x^T y$ is strictly greater than zero.
    
    Doing this a few times gives extremely high confidence that the matrix is PD.
 */
void PDE_Integrator_Test::TestMatrixIsSPD( const CompressedRowMatrix& G )
 {
    const size_t rows = G.Rows();
    _test(rows > 0);

    // Fetch references to the internal storage using friends
    const auto& ia = G.ia;
    const auto& ja = G.ja;
    const auto& a  = G.a;

    // ==========================================
    // 1. SYMMETRY CHECK: A_ij == A_ji
    // ==========================================
    for (size_t i = 0; i < rows; ++i) {
        const int32_t row_start = ia[i];
        const int32_t row_end   = ia[i+1];

        bool diagonal_found = false;

        for (int32_t k = row_start; k < row_end; ++k) {
            const int32_t j = ja[k];
            const double v_ij = a[k];

            if (static_cast<size_t>(j) == i) {
                diagonal_found = true;
                _test(v_ij > 0.0); // Diagonal must be positive
                continue;
            }

            // Find transposed element A_ji via public arrays
            bool transposed_found = false;
            double v_ji = 0.0;
            
            const int32_t t_row_start = ia[j];
            const int32_t t_row_end   = ia[j+1];
            
            for (int32_t m = t_row_start; m < t_row_end; ++m) {
                if (ja[m] == static_cast<int32_t>(i)) {
                    transposed_found = true;
                    v_ji = a[m];
                    break;
                }
            }

            _test(transposed_found);
            // const double epsilon{ 1.0e-12 };
            // _test( approximatelyEqual(v_ij, v_ji, epsilon) );
            _test( approximatelyEqual(v_ij, v_ji) );
        }
        _test(diagonal_found);
    }

    // ==========================================
    // 2. POSITIVE DEFINITE CHECK: x^T * A * x > 0
    // ==========================================
    {
    mt19937 rng(42);
    uniform_real_distribution<double> dist(-1.0, 1.0);

    for (int iter = 0; iter < 3; ++iter) {
        vector<double> x(rows);
        vector<double> y(rows, 0.0);

        for (size_t i = 0; i < rows; ++i) {
            x[i] = dist(rng);
            if (x[i] == 0.0) x[i] = 0.1; 
        }

        // SpMV using public getters
        for (size_t i = 0; i < rows; ++i) {
            const int32_t row_start = ia[i];
            const int32_t row_end   = ia[i+1];
            for (int32_t k = row_start; k < row_end; ++k) {
                y[i] += a[k] * x[ja[k]];
            }
        }

        double energy = 0.0;
        for (size_t i = 0; i < rows; ++i) {
            energy += x[i] * y[i];
        }

        _test(energy > 0.0);
     }
   }

    // ==========================================
    // 2. POSITIVE DEFINITE CHECK: x^T * A * x > 0
    // ==========================================
    {
    // Set up a C++ random number generator for vector x
    mt19937 rng(42); // Fixed seed for deterministic, repeatable tests
    uniform_real_distribution<double> dist(-1.0, 1.0);

    // Run the energy test for a few iterations (e.g., 3 different random vectors)
    for (int iter = 0; iter < 3; ++iter) {
        vector<double> x(rows);
        vector<double> y(rows, 0.0);

        // Populate x with random non-zero values
        for (size_t i = 0; i < rows; ++i) {
            x[i] = dist(rng);
            // Ensure no zeros to guarantee a valid test vector
            if (x[i] == 0.0) x[i] = 0.1; 
        }

        // Perform Sparse Matrix-Vector Multiplication: y = A * x
        for (size_t i = 0; i < rows; ++i) {
            const int32_t row_start = G.ia[i];
            const int32_t row_end   = G.ia[i+1];
            
            for (int32_t k = row_start; k < row_end; ++k) {
                y[i] += G.a[k] * x[G.ja[k]];
            }
        }

        // Calculate the dot product: energy = x^T * y
        double energy = 0.0;
        for (size_t i = 0; i < rows; ++i) {
            energy += x[i] * y[i];
        }

        // If the matrix is Positive Definite, the energy must be strictly > 0
        _test(energy > 0.0);
    }
  }
} // end SPD test




/**
  Generic function to test the pivotVector_ containing accumulated
  matrix-rhs-vector products from eliminated Dirichlet columns.

  @param[in] attorney Attorney to the PDE_Integrator being tested
  @param[in] mesh Mesh of the model used in this test
  @param[in] mat_val Base value input in the matrix constructor (value_label_for_matrix_entry_)
  @param[in] dirich_val Presumed uniform boundary value applied to Dirichlet DOFs
  @param[in] dirich BOX_BOUNDARY codes on which Dirichlet conditions are set
  @param[in] variable_type Type of Test variable (SCALAR, VECTOR, TENSOR)

  @attention This test must be applied after the application of boundary conditions and matrix assembly

  @note Tracking the column data pulled from Dirichlet nodes and moved to the right-hand side (RHS) is tested here.
  Without coupling matrices, an active node $i$ only accumulates pivot entries from an adjacent Dirichlet node $j$  only if they share an element,
  and only within the same component index $k$. Since they are different nodes ($i \neq j$),
  the interaction always uses the off-diagonal element value, which is negative: $-(v_{\text{label}} + k + 10)$.
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator_Test::TestPivotVector(
    const PDE_Integrator_Attorney<dim, CELLTYPE, MATRIXTYPE>& attorney,
    const MeshManager<dim>&       mesh,
    const double                  mat_val,
    const double                  dirich_val,
    const set<BOX_BOUNDARY>& dirich,
    const VARIABLE_TYPE           variable_type )
{
    const size_t nb_nodes = mesh.Nodes();

    // 1. Identify Dirichlet Nodes
    set<size_t> dirich_nodes;
    for (auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n) {
        if (dirich.find(n->AtBoundary()) != dirich.end()) {
            dirich_nodes.insert(n->Idx());
        }
    }

    // 2. Determine variable size dynamically
    uint32_t var_size{1U};
    switch (variable_type) {
        case SCALAR:        var_size = 1U; break;
        case VECTOR:        var_size = dim; break;
        case TENSOR:        var_size = dim * dim; break;
        default:
            throw csmp::Exception(ERROR, "TestPivotVector", "Variable type not supported");
    }

    // 3. Count shared elements between pairs of nodes (N_ij topology mapping)
    // Using a map keyed by a sorted pair of node IDs
    map<pair<size_t, size_t>, size_t> shared_element_count;
    for (auto el = mesh.ElementsBegin(); el != mesh.ElementsEnd(); ++el) {
        for (auto ni = el->NodesBegin(); ni != el->NodesEnd(); ++ni) {
            for (auto nj = el->NodesBegin(); nj != el->NodesEnd(); ++nj) {
                if ((*ni)->Idx() != (*nj)->Idx()) {
                    shared_element_count[{(*ni)->Idx(), (*nj)->Idx()}]++;
                }
            }
        }
    }

    // 4. Validate pivot vector entries for active equations
    size_t active_eq_idx = 0;

    for (size_t i = 0; i < nb_nodes; ++i) {
        // Dirichlet nodes are completely eliminated; skip their equations
        if (dirich_nodes.find(i) != dirich_nodes.end()) {
            continue;
        }

        // Check every degree of freedom for the active node
        for (size_t k = 0; k < var_size; ++k) {
            double expected_pivot_contribution = 0.0;
            const double local_off_diagonal_value = -(mat_val + static_cast<double>(k) + 10.0);

            // Sum up contributions from all adjacent Dirichlet neighbors
            for (const size_t j : dirich_nodes) {
                auto it = shared_element_count.find({i, j});
                if (it != shared_element_count.end()) {
                    const size_t num_shared_elements = it->second;
                    
                    // Product: (Matrix Entry) * (Dirichlet Boundary Value)
                    const double global_matrix_entry = static_cast<double>(num_shared_elements) * local_off_diagonal_value;
                    expected_pivot_contribution += global_matrix_entry * dirich_val;
                }
            }

            // NOTE: If your pipeline stores the pivot vector pre-subtracted 
            // for the RHS modification (b = b - pivot), invert the sign below:
            // expected_pivot_contribution = -expected_pivot_contribution;

            // Safe validation execution outside macro expansion evaluation
            _test(attorney.pivotVector_[active_eq_idx] == expected_pivot_contribution);
            
            // Increment to match the global active system layout tracking
            ++active_eq_idx;
        }
    }
} // end TestPivotVector



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
    TestAssemblySingleVectorNoDirichlet( debug ); // fails

    // --------------------------------------------------------------------------
    // 4. setting up PDE_Integrator for case with a single vector solution variable
    //    and a Dirichlet boundary condition
    // ----------------------------------------------------------------------------
    TestAssemblySingleVectorDirichlet( debug ); // fails

    // 5. setting up PDE_Integrator for coupled system of 2 scalars
    // ----------------------------------------------------------------------------
    TestAssemblyTwoScalarVariablesNoDirichlet( true /* debug */ );
    TestAssemblyTwoScalarVariablesDirichlet( true /* debug */ );

    // 6. setting up PDE_Integrator for coupled system of 1 scalar + 1 vector variable
    // -------------------------------------------------------------------------------
    TestAssemblyScalarAndVectorVariableNoDirichlet( true /* debug */ );
    TestAssemblyScalarAndVectorVariableDirichlet( true /* debug */ );
    
    cout <<"\nPDE_Integrator_Test::TestAssembly: finished test."<< endl;
    
 } // end TestAssembly




/**

Testing a simple case with a single scalar variable

*/
void PDE_Integrator_Test::TestAssemblySingleScalarNoDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarNoDirichlet: starting test."<< endl;
    ResetModelProperties();

    // No Dirichlet conditions
    set<BOX_BOUNDARY> dirich;

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
//    PDE_Integrator<2U,Element,SparseMatrix>  pde_integrator(solver); // OK
    PDE_Integrator<2U>  pde_integrator(solver); // CompressedRowMatrix

    // Enable access to private PDE_Integrator members and methods
    // TODO: do not attempt to access pde_integrator after this, just use attorney
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
    //              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    _test( attorney.lhs_operators_.size() == 1 );
    _test( attorney.rhs_operators_.size() == 1 );
    if ( debug ) {
        cout <<"\nNode neighbors:\n";
        for ( const auto& node : model_region.NodeVector() )
          printNeighbors(node);
        attorney.OutputGlobals();
      }

    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<2U>& mesh = model_->Mesh();
    const size_t nb_nodes = mesh.Nodes();
    _test( attorney.G_.Rows() == nb_nodes );
    _test( attorney.G_.Cols() == nb_nodes );
    _test( attorney.G_.VerifySparsityPattern() );

    model_region.RenumberNodes();
    attorney.Accumulate( model_region );

    TestMatrix<2U>( attorney, mesh, dirich, SCALAR );
    TestRHSVector<2U>( attorney, mesh, val, dirich );

    if( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 );
    }
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarNoDirichlet: finished test."<< endl;
}




/**

Testing a simple case with a single scalar variable and a Dirichlet boundary
condition

*/
void PDE_Integrator_Test::TestAssemblySingleScalarDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarDirichlet: starting test."<< endl;
    ResetModelProperties();

    // Dirichlet condition
    model_->InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.) );
    set<BOX_BOUNDARY> dirich;
    dirich.insert( RIGHT );
    dirich.insert( CNR2 );
    dirich.insert( CNR3 );

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

    // the order here is very important because the pivot vector depends on the material matrices that get accumulated
    model_region.RenumberNodes();
    attorney.Accumulate( model_region );
    attorney.AssignEssentialConditions( model_region );

    TestMatrix<2U>( attorney, mesh, dirich, SCALAR );
    TestRHSVector<2U>( attorney, mesh, val, dirich );

    TestMatrixIsSPD( attorney.G_ );

    if( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 );
    }
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleScalarDirichlet: finished test."<< endl;
}




/**

Test case with a single vector solution variable, no boundary conditions and no
coupling between the variable components

*/
void PDE_Integrator_Test::TestAssemblySingleVectorNoDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleVectorNoDirichlet: starting test."<< endl;
    ResetModelProperties();

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
    model_region.RenumberNodes();
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 1 );
    _test( attorney.rhs_operators_.size() == 1 );
    _test( attorney.G_.VerifySparsityPattern() == 1 );
    
    // sparsity pattern
    if ( debug ) {
        cout <<"\nNode neighbors:\n";
        for ( const auto& node : model_region.NodeVector() )
          printNeighbors(node);
        attorney.OutputGlobals();
      }

    // compare matrix with expected matrix
    const MeshManager<2U>& mesh = model_->Mesh();
    const size_t nb_nodes2 = mesh.Nodes()*2;
    _test( attorney.G_.Rows() == nb_nodes2 );
    _test( attorney.G_.Cols() == nb_nodes2 );

    attorney.Accumulate( model_region );

    TestMatrix<2U>( attorney, mesh, dirich, VECTOR );
    TestRHSVector<2U>( attorney, mesh );

    TestMatrixIsSPD( attorney.G_ );

    // accumulated pattern
    if( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 );
    }
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleVectorNoDirichlet: finished test."<< endl;
}



/**

Test case with a single vector solution variable, a boundary condition and no
coupling between the variable components

*/
void PDE_Integrator_Test::TestAssemblySingleVectorDirichlet( bool debug ) {
    cout <<"\nPDE_Integrator_Test::TestAssemblySingleVectorDirichlet: starting test."<< endl;
    ResetModelProperties();

    // 1. Dirichlet condition assignment
    // ---------------------------------
    model_->InputBoundaryValue( RIGHT, "displacement", makeVector( DIRICH, DIRICH, 0.1, 0.5 ) );
    set<BOX_BOUNDARY> dirich;
    // BOX_BOUNDARY flags constituting the right boundary in 3D models
    dirich.insert( RIGHT );
    dirich.insert( CNR2 ); // within model TINY
    dirich.insert( CNR3 ); // within model TINY

    // Generate node indices 0..n-1
    const Region<2>& model_domain = model_->Region("Model");

    // Get indices of Dirichlet Nodes, for which
    // corresponding line and column are missing in matrix
    set<size_t> dirich_nodes1;
    const MeshManager<2U>& mesh = model_->Mesh();
    for( auto n = mesh.NodesBegin(); n != mesh.NodesEnd(); ++n ) {
        if( dirich.find( n->AtBoundary() ) != dirich.end() ) {
            dirich_nodes1.insert( n->Idx() );
        }
    }
    // searching nodes with Dirichlet variable flags in model
    set<size_t> dirich_nodes2;
    for( const auto& n : model_domain.NodeVector() ) {
        if( dirich.find( n->AtBoundary() ) != dirich.end() ) {
            dirich_nodes2.insert( n->Idx() );
        }
    }
    // are these the same indices found in model domain
    _test( dirich_nodes1 == dirich_nodes2 );



    // 2. PDE_Integrator and FE-assembly
    // ---------------------------------
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
    _test( attorney.EstablishMatrixSetup( model_region ) ); // calls ReduceSystemSizeEliminatingEssentialConditions()
    _test( attorney.lhs_operators_.size() == 1 );
    _test( attorney.rhs_operators_.size() == 1 );
    _test( attorney.G_.VerifySparsityPattern() == 1 );

    // compare matrix with expected matrix
    const size_t remaining_nodes = mesh.Nodes() - dirich_nodes1.size(); // non-boundary nodes
    _test( attorney.G_.Rows() == remaining_nodes * 2 );
    _test( attorney.G_.Cols() == remaining_nodes * 2 );

    attorney.Accumulate( model_region );
    attorney.AssignEssentialConditions( model_region );

    // matrix after accumulation
    if ( debug ) {
        attorney.Out();
        attorney.OutputGlobals();
      }

    TestMatrix<2U>( attorney, mesh, dirich, VECTOR );
    TestRHSVector<2U>( attorney, mesh );

    TestMatrixIsSPD( attorney.G_ );
    
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
void PDE_Integrator_Test::ResetModelProperties()
  {
    // Initialise properties and remove boundary conditions
    model_->InputPropertyValue("fluid pressure", makeScalar(PLAIN, 1.));
    model_->InputPropertyValue("concentration", makeScalar(PLAIN, 1.));
    model_->InputPropertyValue("permeability", makeScalar(PLAIN, 1.));
    model_->InputPropertyValue("fluid volume source", makeScalar(PLAIN, 0.));

    // Find or create additional required properties
    if( !model_->Database().IsDefined( "diffusivity" ) ) {
        model_->CreateProperty( "diffusivity", "D", "m2/s", SCALAR, ELEMENT );
    }
    model_->InputPropertyValue( "diffusivity", makeScalar( ANY, 2. ) );
    
    if( !model_->Database().IsDefined( "Biot alpha" ) ) {
        model_->CreateProperty( "Biot alpha", "alpha", "none", SCALAR, ELEMENT );
    }
    model_->InputPropertyValue( "Biot alpha", makeScalar( ANY, 1.0 ) );

    if( !model_->Database().IsDefined( "Biot modulus" ) ) { // actually the inverse of M
        model_->CreateProperty( "Biot modulus", "M", "none", SCALAR, ELEMENT );
    }
    //                           1/M = (alpha - phi) / beta_r + phi/beta_f = (1-0.2) / 1.0e-10 + 0.2/1.0e-9 =~ 2e-10
    model_->InputPropertyValue( "Biot modulus", makeScalar( ANY, 2.0e-5 ) );

    if( !model_->Database().IsDefined( "Young's modulus" ) ) {
        model_->CreateProperty( "Young's modulus", "E", "Pa", SCALAR, ELEMENT );
    }
    model_->InputPropertyValue( "Young's modulus", makeScalar( ANY, 1.0e+4 ) );

    if( !model_->Database().IsDefined( "Poisson's ratio" ) ) {
        model_->CreateProperty( "Poisson's ratio", "mu", "none", SCALAR, ELEMENT );
    }
    model_->InputPropertyValue( "Poisson's ratio", makeScalar( ANY, 0.25 ) );

    if( !model_->Database().IsDefined( "force" ) ) {
        model_->CreateProperty( "force", "G", "Pa", VECTOR, NODE );
    }
    model_->InputPropertyValue( "force", makeVector( ANY, ANY, 0.0, 0.0 ) );

    if( !model_->Database().IsDefined( "displacement" ) ) {
        model_->CreateProperty( "displacement", "displ", "m", VECTOR, NODE );
    }
    model_->InputPropertyValue( "displacement", makeVector( ANY, ANY, 0.1, 0.5 ) );

    if( !model_->Database().IsDefined( "gravity term" ) ) {
        model_->CreateProperty( "gravity term", "grho", "Pa", VECTOR, ELEMENT );
    }
    // model_->InputPropertyValue( "gravity term", makeVector( ANY, ANY, 0., -9.8 * 2700. ) );
    model_->InputPropertyValue( "gravity term", makeVector( ANY, ANY, 0., 0. ) );

  } // end ResetModelProperties




void PDE_Integrator_Test::run()
  {
    //==========================================================
    // test for a single and two solution variables of same type
    //==========================================================

    // test scalar variable with and without BCs
    // test vector variable and or without BCs
    // test scalar variable and scalar variable with and without BCs
    // test scalar variable and vector variable with and without BCs
    TestAssembly();

    // TODO: testing of accumulation with split boundaries

    // test array variable
    // test flagged array variable


    //=================================================================
    // test two different solution variables (non-square elmt matrices)
    // with cross-coupling terms
    //=================================================================
    TestBlockStructuredCoupledFEM_AssemblyWithDirichletElimination( true /* debug */ );

  } // end run




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
    ResetModelProperties();

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
    LHS_FixedValueMatrix<2U>  lhs2( model_->Database(), "diffusivity", "concentration", "concentration", val1 );
    RHS_FixedValueMatrix<2U>  rhs2( model_->Database(), "fluid volume source", "concentration", val1 );
    attorney.Add(&lhs2);
    attorney.Add(&rhs2);
    
    /// @attention: not testing cross-coupling terms, but they ares tested in the displacement-fluid pressure poroelasticity case 
    
    // TESTING
    Region<2U>& model_region = model_->Region( "Model" );
    model_region.RenumberNodes();
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 2 );
    _test( attorney.rhs_operators_.size() == 2 );

    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<2U>& mesh = model_->Mesh();
    const size_t nb_nodes = mesh.Nodes() * 2U;
    _test( attorney.G_.Rows() == nb_nodes );
    _test( attorney.G_.Cols() == nb_nodes );
    _test( attorney.G_.VerifySparsityPattern() == 1 );

    attorney.Accumulate( model_region );

    TestMultiVariableMatrix<2U>( attorney, model_->Mesh(), attorney.DOF_indexes_ );

    TestMultiVariableRHSVector( attorney,
                                attorney.test_operands_,
                                attorney.DOF_indexes_,
                                model_region );

    if ( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 /* print zero decimal places */ );
     }
    cout <<"\nPDE_Integrator_Test::TestAssemblyTwoScalarVariablesNoDirichlet: finished test."<< endl;
    
} // end TestAssemblyTwoScalarVariablesNoDirichlet





// As method above but with Dirichlet constraints for both variables
void PDE_Integrator_Test::TestAssemblyTwoScalarVariablesDirichlet( bool debug )
 {
    cout <<"\nPDE_Integrator_Test::TestAssemblyTwoScalarVariablesDirichlet: starting test."<< endl;
    ResetModelProperties();
    
    // boundary conditions for both scalar variable
    model_->InputBoundaryValue( RIGHT, "fluid pressure", makeScalar( DIRICH, 2e+6 ) );
    model_->InputBoundaryValue( LEFT, "concentration", makeScalar( DIRICH, 2e-5 ) );

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<2U,Element>  pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<2U> attorney( pde_integrator );

    // Create pde-operators for a system with 2 coupled scalar solution variables (placed on the node)
    const double val = 2.0; // use a single test value because testing functionality only supports this
   
    // scalar 1: 'fluid pressure'
    LHS_FixedValueMatrix<2U>  lhs1( model_->Database(), "permeability", "fluid pressure", "fluid pressure", val );
    RHS_FixedValueMatrix<2U>  rhs1( model_->Database(), "fluid volume source", "fluid pressure", val );
    // assign them to pde integrator
    attorney.Add(&lhs1);
    attorney.Add(&rhs1);

    // scalar 2: 'concentration'
    LHS_FixedValueMatrix<2U>  lhs2( model_->Database(), "diffusivity", "concentration", "concentration", val );
    RHS_FixedValueMatrix<2U>  rhs2( model_->Database(), "fluid volume source", "concentration", val );
    attorney.Add(&lhs2);
    attorney.Add(&rhs2);
    
    // TESTING
    Region<2U>& model_region = model_->Region( "Model" );
    model_region.RenumberNodes();
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 2 );
    _test( attorney.rhs_operators_.size() == 2 );

    // calculating DOFs
    set<BOX_BOUNDARY>  dirich; // which boundary flags do the solution variables have?
    size_t             Dirich_nodes{0};
    for ( const auto& node : model_region.NodeVector() ) {
         if ( isLEFT(node->AtBoundary()) || isRIGHT(node->AtBoundary()) ) {
              dirich.insert( node->AtBoundary() );
              Dirich_nodes++;
           }
      }

    // compare matrix with expected matrix
    const size_t nb_nodes = model_region.Nodes() * 2U - Dirich_nodes;
    _test( attorney.G_.Rows() == nb_nodes );
    _test( attorney.G_.Cols() == nb_nodes );
    _test( attorney.G_.VerifySparsityPattern() == 1 );

    attorney.Accumulate( model_region );
    attorney.AssignEssentialConditions( model_region );

    if ( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 /* print zero decimal places */ );
     }

    TestMultiVariableMatrix<2U>( attorney, model_->Mesh(), attorney.DOF_indexes_ );

    TestMultiVariableRHSVector( attorney,
                                attorney.test_operands_,
                                attorney.DOF_indexes_,
                                model_region );

    cout <<"\nPDE_Integrator_Test::TestAssemblyTwoScalarVariablesDirichlet: finished test."<< endl;
    
} // end TestAssemblyTwoScalarVariablesDirichlet





/**

Testing a  case with a  scalar variable and a vector variable coupled together.

SKM - printing matrix in integer format to better see pattern.

test_CreateSimplestPolyElement2DModel() - 2 triangles + 1 quadrilateral

*/
void PDE_Integrator_Test::TestAssemblyScalarAndVectorVariableNoDirichlet( bool debug )
 {
    cout <<"\nPDE_Integrator_Test::TestAssemblyScalarAndVectorVariableNoDirichlet: starting test."<< endl;
    ResetModelProperties();
    const uint32_t model_dimensions{2U};

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<model_dimensions,Element>  pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<model_dimensions> attorney( pde_integrator );

    // Create pde-operators for a system with 1 vector coupled with 1 scalar solution variable
    const double val1 = 1.;

    // vector 1: 'nodal velocity'
    LHS_FixedValueMatrix<model_dimensions>  lhs1( model_->Database(), "permeability", "nodal velocity", "nodal velocity", val1 );
    RHS_FixedValueMatrix<model_dimensions>  rhs1( model_->Database(), "fluid volume source", "nodal velocity", val1 );
    // assign them to pde integrator
    attorney.Add(&lhs1);
    attorney.Add(&rhs1);

    // scalar 1: 'concentration'
    LHS_FixedValueMatrix<model_dimensions>  lhs2( model_->Database(), "diffusivity", "concentration", "concentration", val1 );
    RHS_FixedValueMatrix<model_dimensions>  rhs2( model_->Database(), "fluid volume source", "concentration", val1 );
    attorney.Add(&lhs2);
    attorney.Add(&rhs2);
    
    // TESTING
    const Region<model_dimensions>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 2 );
    _test( attorney.rhs_operators_.size() == 2 );

    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<model_dimensions>& mesh = model_->Mesh();
    //                                  nodes * ( 1-for-scalars + dim-for-vectors)
    const size_t n_degrees_of_freedom = mesh.Nodes() + mesh.Nodes() * model_dimensions;
    _test( attorney.G_.Rows() == n_degrees_of_freedom );
    _test( attorney.G_.Cols() == n_degrees_of_freedom );

    attorney.Accumulate( model_region );

    _test( attorney.G_.VerifySparsityPattern() == 1 );

    // allows for cross-coupling terms etc.
    TestMultiVariableMatrix<2U>( attorney, model_->Mesh(), attorney.DOF_indexes_ );
    TestMultiVariableRHSVector( attorney, attorney.test_operands_, attorney.DOF_indexes_, model_region );

    if ( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 0 /* print zero decimal places */ );
     }
    cout <<"\nPDE_Integrator_Test::TestAssemblyScalarAndVectorVariableNoDirichlet: finished test."<< endl;
    
} // end TestAssemblyScalarAndVectorVariableNoDirichlet





void PDE_Integrator_Test::TestAssemblyScalarAndVectorVariableDirichlet( bool debug )
 {
    cout <<"\nPDE_Integrator_Test::TestAssemblyScalarAndVectorVariableDirichlet: starting test."<< endl;
    ResetModelProperties();
    const uint32_t model_dimensions{2U};
    
        // boundary conditions for both scalar variable
    model_->InputBoundaryValue( RIGHT, "nodal velocity", makeVector( DIRICH,PLAIN, 2e+6, 0. ) );
    model_->InputBoundaryValue( LEFT, "concentration", makeScalar( DIRICH, 2e-5 ) );

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<model_dimensions,Element>  pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<model_dimensions> attorney( pde_integrator );

    // Create pde-operators for a system with 1 vector coupled with 1 scalar solution variable
    const double val1 = 3.;

    // vector 1: 'nodal velocity'
    LHS_FixedValueMatrix<model_dimensions>  lhs1( model_->Database(), "permeability", "nodal velocity", "nodal velocity", val1 );
    RHS_FixedValueMatrix<model_dimensions>  rhs1( model_->Database(), "fluid volume source", "nodal velocity", val1 );
    // assign them to pde integrator
    attorney.Add(&lhs1);
    attorney.Add(&rhs1);

    // scalar 1: 'concentration'
    LHS_FixedValueMatrix<model_dimensions>  lhs2( model_->Database(), "diffusivity", "concentration", "concentration", val1 );
    RHS_FixedValueMatrix<model_dimensions>  rhs2( model_->Database(), "fluid volume source", "concentration", val1 );
    attorney.Add(&lhs2);
    attorney.Add(&rhs2);
    
    // TODO: to extend this test for cross-coupling one needs to implement suitable Jacobian PDE operators
    // cross-coupling terms assembled into upper diagonal of matrix
    // ------------------------------------------------------------
    // coupling basic operand 'concentration' rows with test operand 'nodal velocity' columns
//    LHS_FixedValueMatrix<model_dimensions>  lhs12( model_->Database(), "diffusivity", "concentration", "nodal velocity", val1 );
//    attorney.Add(&lhs12);

    // cross-coupling terms assembled into lower diagonal of matrix
    // ------------------------------------------------------------
    // coupling 'nodal velocity' with 'concentration'
//    LHS_FixedValueMatrix<model_dimensions>  lhs21( model_->Database(), "diffusivity", "nodal velocity", "concentration", val1 );
//    attorney.Add(&lhs21);

    
    // TESTING
    const Region<model_dimensions>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 2 );
    _test( attorney.rhs_operators_.size() == 2 );
    _test( attorney.G_.VerifySparsityPattern() == 1 );

    // counting the DOFs eliminated due to Dirichlet boundary conditions
    const csmp::Index v_key = model_->Database().StorageKey("nodal velocity");
    const csmp::Index c_key = model_->Database().StorageKey("concentration");
    size_t Dirich_constraints{0};
    for ( const auto* node : model_region.NodeVector() ) {
      for ( uint32_t i{0}; i<model_dimensions; ++i )
        if ( node->Status(v_key,i) == DIRICH ) Dirich_constraints++;
      if ( node->Status(c_key) == DIRICH ) Dirich_constraints++;
    }

    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<model_dimensions>& mesh = model_->Mesh();
    //                                  nodes * ( 1-for-scalars + dim-for-vectors)
    const size_t n_degrees_of_freedom = mesh.Nodes() + mesh.Nodes() * model_dimensions - Dirich_constraints;
    _test( attorney.G_.Rows() == n_degrees_of_freedom );
    _test( attorney.G_.Cols() == n_degrees_of_freedom );

    attorney.Accumulate( model_region );
    attorney.AssignEssentialConditions( model_region );

    if ( debug ) {
        attorney.Out();
        attorney.OutputGlobals( 2 /* print zero decimal places */ );
     }

    set<BOX_BOUNDARY>  dirich = { LEFT, RIGHT };

    TestMultiVariableMatrix<2U>( attorney, model_->Mesh(), attorney.DOF_indexes_ );
    TestMultiVariableRHSVector( attorney, attorney.test_operands_, attorney.DOF_indexes_, model_region );

    cout <<"\nPDE_Integrator_Test::TestAssemblyScalarAndVectorVariableDirichlet: finished test."<< endl;
    
} // end TestAssemblyScalarAndVectorVariableDirichlet




/**
We are assembling a Biot-type poromechanics system coupling:

Displacement , u
u (vector, 2 DOFs per node:  ux,uy,ux,ux = interleaved)
Pore pressure  p (scalar, 1 DOF per node)
For the TINY mesh: 6 nodes

Governing System

The coupled block system is:

|K   Q|   {u}        |fu|
|         |           =
|QT  S|   { p}      |fp|

K	Mechanical stiffness matrix 12 x 12
Q	Coupling matrix — volumetric strain drives pressure 12 x 6
S	Compressibility / storage matrix  6 x 6
fu  Mechanical force vector
fp  Fluid source vector

Element Contributions

Stiffness Matrix 
K
K (Mechanical Block)

For a linear triangle (plane strain), the element stiffness is:

Ke = t Ae BT D B

where

B is the strain-displacement matrix, 

D is the elastic constitutive matrix, 

Ae   is the element area, and
t is thickness.

For the quadrilateral (Elmt2), standard 2×2
2×2 Gauss quadrature applies:

Ke=t∫ΩeBTDBdΩ

Coupling Matrix
Q (Biot Coupling)
The coupling matrix links volumetric strain to pore pressure:

Qe=α∫ΩeBTmNpdΩ

where:

α is the Biot coefficient

m=[1,1,0]T
extracts volumetric strain

Np  are the scalar pressure shape functions
B operates on displacement shape functions

For a constant-strain triangle with area Ae:

Qe=αtAeBTmNp

where Nˉp=1/3[1,1,1] (centroid evaluation for linear triangle).

Storage Matrix S

Se=1/M ∫ΩeNpTNp dΩ

where M is the Biot modulus. For a linear triangle:

Se= t Ae/12M [211...

The reduced system to solve is:

•	AFFxF=bF−AFDxD
 
The rows and columns corresponding to  D are fully removed from the assembled matrix.

Assembly Algorithm

# Global sizes (before elimination)
n_disp = 2 * n_nodes   # = 12
n_pres = n_nodes       # = 6
n_total = n_disp + n_pres  # = 18

# Allocate global blocks
K  = zeros(n_disp, n_disp)
Q  = zeros(n_disp, n_pres)
S  = zeros(n_pres, n_pres)
f_u = zeros(n_disp)
f_p = zeros(n_pres)

for element in [Elmt0, Elmt1, Elmt2]:
    nodes = element.connectivity          # local node indices
    
    # Local displacement DOF indices (interleaved)
    u_dofs = interleave(nodes)            # [2*n, 2*n+1 for n in nodes]
    p_dofs = nodes                        # pressure DOFs = node indices
    
    Ke, Qe, Se = compute_element_matrices(element)
    fe_u, fe_p = compute_element_forces(element)
    
    # Scatter into global blocks
    K[u_dofs, u_dofs]  += Ke
    Q[u_dofs, p_dofs]  += Qe
    S[p_dofs, p_dofs]  += Se
    f_u[u_dofs]        += fe_u
    f_p[p_dofs]        += fe_p

# --- Assemble full block system ---
A = block([[K,   Q  ],
           [Q.T, S  ]])
b = concatenate([f_u, f_p])

# --- Apply Dirichlet BCs by full elimination ---
free_u = [i for i in range(n_disp) if i not in dirichlet_u_dofs]
free_p = [i for i in range(n_pres) if i not in dirichlet_p_dofs]

# Map free DOFs into global block indices
free_global = free_u + [n_disp + i for i in free_p]

# Modify RHS for non-zero Dirichlet values
b[free_global] -= A[free_global, :][:, dirichlet_global] @ x_dirichlet

# Extract reduced system
A_red = A[free_global, :][:, free_global]
b_red = b[free_global]

# Solve
x_free = solve(A_red, b_red)

Interleaved vs. Block Displacement Storage

The RHS is interleaved for the displacement part, but the global system is block-structured (
u
u block then 
p
p block). These are compatible — the interleaving only affects the internal ordering within the 
K
K block.

Symmetry of the Coupled System

Note that QTQ T  appears in the lower-left block. For static (undrained) problems the full system is symmetric. For transient Biot (with a time derivative on
S), the system remains symmetric after time discretisation.

Line Element (Elmt — not shown in mesh diagram)

The 1D line element contributes only to K and f (e.g., a structural reinforcement or interface element)
it carries no pressure DOFs unless explicitly modelled as a fracture/conduit.

 */
void PDE_Integrator_Test::TestBlockStructuredCoupledFEM_AssemblyWithDirichletElimination( bool )
 {
    cout <<"\nPDE_Integrator_Test::TestBlockStructuredCoupledFEM_AssemblyWithDirichletElimination: starting test."<< endl;
    ResetModelProperties(); // includes displacement, Biot alpha etc.
    constexpr uint32_t DIM{2U};
    
    // boundary conditions for displacement and pressure
    model_->InputBoundaryValue( LEFT, "displacement", makeVector( DIRICH,PLAIN, 0., 0. ) ); // displ. left 0.,0.
    Region<DIM>& model_domain = model_->Region("Model");

    const csmp::Index displ_key = model_->Database().StorageKey("displacement");
    // ascertain the same node numbers as in the Claude test dataset
    for ( size_t i{0}; i<model_domain.Nodes(); ++i ) model_domain.N(i)->Idx(i);
    // LEFT
    model_domain.N(0)->Store( displ_key, makeVector( DIRICH,DIRICH, 0., 0. ) ); // displ. left 0.,0.
    model_domain.N(3)->Store( displ_key, makeVector( DIRICH,DIRICH, 0., 0. ) );
    // RIGHT
    model_domain.N(2)->Store( displ_key, makeVector( DIRICH,PLAIN, 0.01, 0. ) );
    model_domain.N(5)->Store( displ_key, makeVector( DIRICH, DIRICH, 0.01, 0. ) );
    // BOTTOM
    model_domain.N(4)->Store( displ_key, makeVector( PLAIN, DIRICH, 0., 0. ) );

    const csmp::Index p_key = model_->Database().StorageKey("fluid pressure");

    model_domain.N(0)->Store( p_key, makeScalar( DIRICH, 0. ) );
    model_domain.N(1)->Store( p_key, makeScalar( DIRICH, 0. ) );
    model_domain.N(2)->Store( p_key, makeScalar( DIRICH, 500. ) );
    model_domain.N(3)->Store( p_key, makeScalar( DIRICH, 0. ) );
    model_domain.N(5)->Store( p_key, makeScalar( DIRICH, 500. ) );

    CSMP_DEFAULT_LINEAR_SOLVER   solver;
    PDE_Integrator<DIM,Element>  pde_integrator(solver);

    // Enable access to private PDE_Integrator members and methods
    PDE_Integrator_Attorney<DIM> attorney( pde_integrator );

    // Create pde-operators for a system with 1 vector coupled with 1 scalar solution variable
    
    // Lefthandside LHS
    // ----------------
    
    // primary variables (plain strain assumption)
    
    // stiffness, K
    NumIntegral_BT_D_B_dV<DIM> lhs_K( model_->Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
    attorney.Add(&lhs_K);
    
    // storage matrix S = 1/M NumInt NT N dV            // 1/M = (alpha - phi) / Ks + phi/Kf ~= compressibility
    NumIntegral_NT_lhsop_N_dV<DIM>  lhs_S( model_->Database(), "Biot modulus", "fluid pressure", "fluid pressure" );
    attorney.Add(&lhs_S);

    // cross-coupling terms

    // Biot (cross-) coupling matrix, Q (Jacobian in upper diagonal)
    NumJacobianIntegral_BT_m_N_dV<DIM> lhs_Q( model_->Database(), "Biot alpha", "fluid pressure", "displacement" );
    lhs_Q.SubtractAccumulate(); // subtract the pore pressure from compressive stress to get effective stress
    attorney.Add(&lhs_Q);
    
    // Biot (cross-) coupling matrix, QT (=Q transposed Jacobian in lower diagonal)
    NumJacobianIntegral_N_mT_B_dV<DIM> lhs_QT( model_->Database(), "Biot alpha", "displacement", "fluid pressure" );
    lhs_QT.SubtractAccumulate(); // subtract the pore pressure from compressive stress to get effective stress
    attorney.Add(&lhs_QT);

    // Righthandside RHS
    // -----------------
    // for each test operand, a RHS operator with this test operand is needed
    
    // body forces
    NumIntegral_PT_op_dV<DIM>  Fu( model_->Database(), "gravity term", "displacement" );
    attorney.Add(&Fu);
    
    NumIntegral_NT_op_N_dV<DIM>  Fp( model_->Database(), "fluid volume source", "fluid pressure" );
    attorney.Add(&Fp);


    // TESTING
    const Region<DIM>& model_region = model_->Region( "Model" );
    _test( attorney.EstablishMatrixSetup( model_region ) );
    _test( attorney.lhs_operators_.size() == 4 );
    _test( attorney.rhs_operators_.size() == 2 );

    // counting the DOFs eliminated due to Dirichlet boundary conditions
    size_t Dirich_constraints{0};
    for ( const auto* node : model_region.NodeVector() ) {
      for ( uint32_t i{0}; i<DIM; ++i )
        if ( node->Status(displ_key,i) == DIRICH ) Dirich_constraints++;
      if ( node->Status(p_key) == DIRICH ) Dirich_constraints++;
    }
    if ( verbose_ ) {
       cout <<"\nNodal boundary conditions applied to solution variables in model domain: "<< endl;
       // displacement
       cout <<"\n"<<"displacement:"<< endl;
       for ( size_t i{0}; i<model_region.Nodes(); ++i ) {
             cout <<"\t"<< model_region.N(i)->Idx() <<": "<< parseBoundary( model_region.N(i)->AtBoundary() ) <<": ";
             for ( uint32_t j{0}; j<DIM; ++j )
               cout << parseStatus( model_region.N(i)->Status(displ_key,j) ) <<" ";
             cout << endl;
          }
       // fluid pressure
       cout <<"\n"<<"fluid pressure:"<< endl;
       for ( size_t i{0}; i<model_region.Nodes(); ++i ) {
             cout <<"\t"<< model_region.N(i)->Idx() <<": "<< parseBoundary( model_region.N(i)->AtBoundary() );
             cout <<": "<< parseStatus( model_region.N(i)->Status(p_key) ) <<" ";
             cout << endl;
          }
    }
    
    // compare matrix with expected matrix
    // No Dirichlet conditions: Matrix should be nb nodes x nb nodes
    const MeshManager<DIM>& mesh = model_->Mesh();
    //                                  nodes * ( 1-for-scalars + dim-for-vectors)
    const size_t n_degrees_of_freedom = mesh.Nodes() + mesh.Nodes() * DIM - Dirich_constraints;
    _test( attorney.G_.Rows() == n_degrees_of_freedom );
    _test( attorney.G_.Cols() == n_degrees_of_freedom );

    // ascertain the same node numbers as in the Claude test dataset
    for ( size_t i{0}; i<model_domain.Nodes(); ++i ) model_domain.N(i)->Idx(i);
    
    attorney.Accumulate( model_region );
    attorney.AssignEssentialConditions( model_region );

    _test( attorney.G_.VerifySparsityPattern() == 1 );

    set<BOX_BOUNDARY>  dirich = { LEFT, RIGHT };

                                                        
    // TESTING
    // test actual matrix and vector entries against independently computed values
    CompressedRowMatrix  A_ref;
    vector<double>       b_ref;
    GenerateTinyPoromechanicsTestSystem( A_ref, b_ref );
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if ( verbose_ ) {
        // A reference solution
        cout <<"\n"<<"Claude test matrix:"<< endl;
        PrintTinyTestSystem( A_ref, b_ref );
        // A accumulated by PDE_Integrator
        printRangeOfVariable( *model_, "Young's modulus", true );
        printRangeOfVariable( *model_, "Poisson's ratio", true );
        printRangeOfVariable( *model_, "Biot alpha", true );
        printRangeOfVariable( *model_, "Biot modulus", true );
        printRangeOfVariable( *model_, "displacement", true );
        printRangeOfVariable( *model_, "fluid pressure", true );
        cout <<"\n"<<"CSMP matrix:"<< endl;
        PrintTinyTestSystem( attorney.G_, attorney.rh_ );
        attorney.OutputGlobals( 1 /* print zero decimal places */ );
     }
     
    // testing
    // -----------------------------------------------------------------------
    // Basic size checks
    // -----------------------------------------------------------------------
    _test( A_ref.Rows() == 5U );
    _test( A_ref.Cols() == 5U );
    _test( b_ref.size() == 5U );
    _test( A_ref.VerifySparsityPattern() );

    // -----------------------------------------------------------------------
    // Symmetry check
    // -----------------------------------------------------------------------
    constexpr double sym_tol = 1.0e-8;
    for ( size_t i=0; i<5U; ++i )
        for ( size_t j=0; j<5U; ++j )
        {
            std::ostringstream msg;
            msg << "Symmetry A(" << i << "," << j << ") vs A(" << j << "," << i << ")";
            _test( std::abs( A_ref(i,j) - A_ref(j,i) ) < sym_tol );
        }

    // -----------------------------------------------------------------------
    // Positive diagonal check (coercive system)
    // -----------------------------------------------------------------------
    for ( size_t i=0; i<4U; ++i )   // mechanical DOFs only — storage diagonal ~0
    {
        std::ostringstream msg;
        msg << "Positive diagonal A(" << i << "," << i << ") = " << A_ref(i,i);
        _test( A_ref(i,i) > 0.0 );
    }

    // -----------------------------------------------------------------------
    // Reference values
    //
    // Free DOF layout:
    //   [0] u_x(1)   [1] u_y(1)   [2] u_y(2)   [3] u_x(4)   [4] p(4)
    //
    // Material: E=10000, nu=0.25, alpha=1, plane strain
    // Mesh:     TINY — 2 triangles + 1 quad + 1 line element
    // BCs:      Left u=0,p=0 | Bottom u_y=0 | Right u_x=0.01,p=500
    // -----------------------------------------------------------------------
    constexpr double mech_tol    = 1.0e-3;   // mechanical stiffness entries
    constexpr double coupling_tol= 1.0e-4;   // coupling entries (small values)
    constexpr double rhs_tol     = 1.0e-3;   // RHS entries

    // --- Mechanical block K (rows/cols 0-3) ---

    // Diagonal
    _equal( A_ref(0,0), 18990.1876, mech_tol );   // u_x(1): tri0+tri1+quad+line
    _equal( A_ref(1,1), 18990.1876, mech_tol );   // u_y(1): tri0+tri1+quad+line
    _equal( A_ref(2,2),  5333.3333, mech_tol );   // u_y(2): quad only
    _equal( A_ref(3,3), 13333.3333, mech_tol );   // u_x(4): tri1+quad

    // Off-diagonal mechanical
    _equal( A_ref(0,1),   828.4271, mech_tol );   // u_x(1) - u_y(1)
    _equal( A_ref(1,0),   828.4271, mech_tol );   // symmetric
    _equal( A_ref(0,3), -1333.3333, mech_tol );   // u_x(1) - u_x(4)
    _equal( A_ref(3,0), -1333.3333, mech_tol );   // symmetric
    _equal( A_ref(1,2),   666.6667, mech_tol );   // u_y(1) - u_y(2)
    _equal( A_ref(2,1),   666.6667, mech_tol );   // symmetric
    _equal( A_ref(1,3),  2000.0000, mech_tol );   // u_y(1) - u_x(4)
    _equal( A_ref(3,1),  2000.0000, mech_tol );   // symmetric
    _equal( A_ref(2,3), -2000.0000, mech_tol );   // u_y(2) - u_x(4)
    _equal( A_ref(3,2), -2000.0000, mech_tol );   // symmetric
    _equal( A_ref(0,2),   0.0000,   mech_tol );   // u_x(1) - u_y(2): no connection
    _equal( A_ref(2,0),   0.0000,   mech_tol );   // symmetric

    // --- Coupling block Q (col 4 = p(4), rows 0-3) ---
    // Q = -alpha * integral( B^T m Np dV )
    // sigma' = sigma - alpha*p*I  =>  compression drives positive pressure

    _equal( A_ref(0,4),  0.0833, coupling_tol );   // Q(u_x(1), p(4))
    _equal( A_ref(4,0),  0.0833, coupling_tol );   // Q^T symmetric
    _equal( A_ref(1,4), -0.3333, coupling_tol );   // Q(u_y(1), p(4))
    _equal( A_ref(4,1), -0.3333, coupling_tol );   // Q^T symmetric
    _equal( A_ref(2,4), -0.0833, coupling_tol );   // Q(u_y(2), p(4))
    _equal( A_ref(4,2), -0.0833, coupling_tol );   // Q^T symmetric
    _equal( A_ref(3,4),  0.0000, coupling_tol );   // Q(u_x(4), p(4)): cancels
    _equal( A_ref(4,3),  0.0000, coupling_tol );   // Q^T symmetric

    // --- Storage block S (row/col 4 = p(4)) ---
    // S = (1/M) * integral( Np^T Np dV )
    // Very small relative to mechanical stiffness (~1e-5)
    _equal( A_ref(4,4),  0.0000, 1.0e-4 );   // ~3e-6, prints as 0 at 4dp

    // -----------------------------------------------------------------------
    // RHS vector b_red
    //
    // Driven by: -A_FD * x_D
    // Non-zero Dirichlet: u_x(2)=0.01, u_x(5)=0.01, p(2)=500, p(5)=500
    // -----------------------------------------------------------------------
    _equal( b_ref[0],  -65.0000, rhs_tol );   // u_x(1)
    _equal( b_ref[1],   63.3333, rhs_tol );   // u_y(1)
    _equal( b_ref[2],  146.6667, rhs_tol );   // u_y(2)
    _equal( b_ref[3],  -65.0000, rhs_tol );   // u_x(4)
    _equal( b_ref[4],    0.0017, 1.0e-4  );   // p(4): small storage contribution

    // -----------------------------------------------------------------------
    // Cross-check: matrix-vector product A*x should be computable
    // (just verifies the CRM is functional, not the solution)
    // -----------------------------------------------------------------------
    std::vector<double> x_test(5, 1.0);
    std::vector<double> result;
    A_ref.MultiplyWith( x_test, result );
    _test( result.size() == 5U );

    _succeed();
     
    cout <<"\nPDE_Integrator_Test::TestAssemblyScalarAndVectorVariableDirichlet: finished test."<< endl;
    
} // end TestBlockStructuredCoupledFEM_AssemblyWithDirichletElimination


// ==================================================================================================
// GENERATE REDUCED SYSTEM FOR TESTING (Claude Sonnet 4.6)
// ==================================================================================================


/**
 * @brief Prints all intermediate assembly quantities for diagnosis.
 *        Insert immediately after the three assemble calls and
 *        before building A_full.
 */
static void diagnosePrintAssembly( const vector<double>& K_gl,
                                   const vector<double>& Q_gl,
                                   const vector<double>& S_gl,
                                   size_t N_DISP, size_t N_PRES )
{
    cout << fixed << setprecision(6);

    // --- K diagonal ---
    cout << "\nK diagonal (displacement block):\n";
    for ( size_t i=0; i<N_DISP; ++i )
        cout << "  K(" << i << "," << i << ") = "
                  << K_gl[i*N_DISP+i] << "\n";

    // --- Q column sums (should be non-zero where coupling exists) ---
    cout << "\nQ matrix (disp x pres):\n";
    for ( size_t r=0; r<N_DISP; ++r ) {
        cout << "  row " << r << ": ";
        for ( size_t c=0; c<N_PRES; ++c )
            cout << setw(12) << Q_gl[r*N_PRES+c];
        cout << "\n";
    }

    // --- S diagonal ---
    cout << "\nS diagonal (pressure block):\n";
    for ( size_t i=0; i<N_PRES; ++i )
        cout << "  S(" << i << "," << i << ") = "
                  << S_gl[i*N_PRES+i] << "\n";
}


/**
 * @brief Generates the reduced (Dirichlet-eliminated) coupled poromechanics
 *        system for the TINY test mesh.
 *
 * Mesh layout:
 *
 *   0 (0,1) ------- 1 (1,1) ------- 2 (2,1)
 *   |  \                      |                      |
 *   |       \ [Elmt0]  |                      |
 *   |               \         |   [Elmt2]     |
 *   |   [Elmt1]        |                      |
 *   |                     \  |                      |
 *   3 (0,0) ------- 4 (1,0) ------- 5 (2,0)
 *
 * DOF layout (18 total before elimination):
 *   Displacement (interleaved): u_x(i) = 2*i,  u_y(i) = 2*i+1   [0..11]
 *   Pressure:                   p(i)   = 12+i                     [12..17]
 *
 * Dirichlet constraints:
 *   DOF  0  -> u_x node 0 = 0.0
 *   DOF  1  -> u_y node 0 = 0.0
 *   DOF  4  -> u_x node 2 = 0.01
 *   DOF  6  -> u_x node 3 = 0.0
 *   DOF  7  -> u_y node 3 = 0.0
 *   DOF  9  -> u_y node 4 = 0.0
 *   DOF 10  -> u_x node 5 = 0.01
 *   DOF 11  -> u_y node 5 = 0.0
 *   DOF 12  -> p   node 0 = 0.0
 *   DOF 13  -> p   node 1 = 0.0
 *   DOF 14  -> p   node 2 = 500.0
 *   DOF 15  -> p   node 3 = 0.0
 *   DOF 17  -> p   node 5 = 500.0
 *
 * Free DOFs (reduced system, 5 unknowns):
 *   Reduced 0 -> Global  2  (u_x node 1)
 *   Reduced 1 -> Global  3  (u_y node 1)
 *   Reduced 2 -> Global  5  (u_y node 2)
 *   Reduced 3 -> Global  8  (u_x node 4)
 *   Reduced 4 -> Global 16  (p   node 4)
 *
 * @param[out] A_red   Reduced 5x5 system matrix (CompressedRowMatrix)
 * @param[out] b_red   Reduced RHS vector (length 5)
 */
void PDE_Integrator_Test::GenerateTinyPoromechanicsTestSystem( CompressedRowMatrix& A_red,
                                                               vector<double>& b_red )
{
    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------
    constexpr size_t N_NODES = 6;
    constexpr size_t N_DISP  = 2 * N_NODES;   // 12
    constexpr size_t N_PRES  = N_NODES;        // 6
    constexpr size_t N_TOTAL = N_DISP + N_PRES; // 18

    // -----------------------------------------------------------------------
    // Material properties
    // -----------------------------------------------------------------------
    constexpr double E     = 10000.0;
    constexpr double nu    = 0.25;
    constexpr double alpha = 1.0;
    constexpr double M_bio = 50000.0;

    // Plane strain D matrix (3x3, stored row-major)
    // D = E*(1-nu)/((1+nu)*(1-2*nu)) * [[1, nu/(1-nu), 0],
    //                                    [nu/(1-nu), 1, 0],
    //                                    [0, 0, (1-2*nu)/(2*(1-nu))]]
    const double fac    = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0*nu));
    const double D_diag = fac;                          // 12000
    const double D_off  = fac * nu / (1.0 - nu);       //  4000
    const double D_shear= fac * (1.0-2.0*nu)/(2.0*(1.0-nu)); // 4000

    // D[row][col]
    const array<array<double,3>,3> D = {{
        {{ D_diag, D_off,   0.0     }},
        {{ D_off,  D_diag,  0.0     }},
        {{ 0.0,    0.0,     D_shear }}
    }};

    // Volumetric coupling vector m = [1, 1, 0]^T
    const std::array<double,3> m_vec = { -1.0, -1.0, 0.0 };  // compression is positive

    // -----------------------------------------------------------------------
    // Node coordinates [node][x,y]
    // -----------------------------------------------------------------------
    const array<array<double,2>, N_NODES> coords = {{
        {{ 0.0, 1.0 }},  // node 0
        {{ 1.0, 1.0 }},  // node 1
        {{ 2.0, 1.0 }},  // node 2
        {{ 0.0, 0.0 }},  // node 3
        {{ 1.0, 0.0 }},  // node 4
        {{ 2.0, 0.0 }}   // node 5
    }};

    // -----------------------------------------------------------------------
    // Global matrices (dense, row-major, before elimination)
    // -----------------------------------------------------------------------
    // Using flat arrays: index = row*cols + col
    vector<double> K_gl(N_DISP  * N_DISP,  0.0);
    vector<double> Q_gl(N_DISP  * N_PRES,  0.0);
    vector<double> S_gl(N_PRES  * N_PRES,  0.0);
    vector<double> f_u (N_DISP,            0.0);
    vector<double> f_p (N_PRES,            0.0);

    // Accessors
    auto K = [&](size_t r, size_t c) -> double& { return K_gl[r*N_DISP + c]; };
    auto Q = [&](size_t r, size_t c) -> double& { return Q_gl[r*N_PRES + c]; };
    auto S = [&](size_t r, size_t c) -> double& { return S_gl[r*N_PRES + c]; };

    // -----------------------------------------------------------------------
    // Utility: interleaved displacement DOF indices for a node list
    // -----------------------------------------------------------------------
    auto u_dofs = []( const vector<size_t>& nodes ) -> vector<size_t>
    {
        vector<size_t> dofs;
        dofs.reserve( 2 * nodes.size() );
        for ( size_t n : nodes ) { dofs.push_back(2*n); dofs.push_back(2*n+1); }
        return dofs;
    };

    // -----------------------------------------------------------------------
    // Matrix multiply helpers (small dense, no external library needed)
    // -----------------------------------------------------------------------

    // C(m,p) = A(m,n) * B(n,p)
    auto matmul = []( const vector<double>& A, size_t m, size_t n,
                      const vector<double>& B, size_t p )
        -> vector<double>
    {
        vector<double> C(m*p, 0.0);
        for ( size_t i=0; i<m; ++i )
            for ( size_t k=0; k<n; ++k )
                for ( size_t j=0; j<p; ++j )
                    C[i*p+j] += A[i*n+k] * B[k*p+j];
        return C;
    };

    // Transpose: A(m,n) -> AT(n,m)
    auto transpose = []( const vector<double>& A, size_t m, size_t n )
        -> vector<double>
    {
        vector<double> AT(n*m);
        for ( size_t i=0; i<m; ++i )
            for ( size_t j=0; j<n; ++j )
                AT[j*m+i] = A[i*n+j];
        return AT;
    };


    // -----------------------------------------------------------------------
    // Triangle assembly
    // Nodes: node_ids[0,1,2] with coordinates from coords[]
    // -----------------------------------------------------------------------
    auto assemble_triangle = [&]( const std::vector<size_t>& nids )
    {
        constexpr size_t nn = 3;
        constexpr size_t ud = 2 * nn;
        constexpr size_t pd = nn;

        // 3-point Gauss rule matching CSMP
        const std::array<double,3> rr = { 0.5, 0.5, 0.0 };
        const std::array<double,3> ss = { 0.0, 0.5, 0.5 };
        const std::array<double,3> ww = { 1.0/6.0, 1.0/6.0, 1.0/6.0 };

        // XY matrix: XY(j,i) = coordinate i of node j
        // matches CSMP's CoordinateMatrix() layout
        const std::array<std::array<double,2>,3> XY = {{
            {{ coords[nids[0]][0], coords[nids[0]][1] }},
            {{ coords[nids[1]][0], coords[nids[1]][1] }},
            {{ coords[nids[2]][0], coords[nids[2]][1] }}
        }};

        std::vector<double> Ke(ud*ud, 0.0);
        std::vector<double> Qe(ud*pd, 0.0);
        std::vector<double> Se(pd*pd, 0.0);

        for ( size_t gp=0; gp<3; ++gp )
        {
            const double r = rr[gp];
            const double s = ss[gp];
            const double w = ww[gp];

            // Shape function derivatives in reference coords
            // N0=1-r-s, N1=r, N2=s
            const std::array<double,3> DNR = { -1.0,  1.0,  0.0 };
            const std::array<double,3> DNS = { -1.0,  0.0,  1.0 };

            // Shape functions
            const std::array<double,3> Np = { 1.0-r-s, r, s };

            // Jacobian: JAC(0,i) = sum_j DNR[j]*XY(j,i)
            //           JAC(1,i) = sum_j DNS[j]*XY(j,i)
            double JAC[2][2] = {{0,0},{0,0}};
            for ( size_t j=0; j<nn; ++j ) {
                JAC[0][0] += DNR[j] * XY[j][0];
                JAC[0][1] += DNR[j] * XY[j][1];
                JAC[1][0] += DNS[j] * XY[j][0];
                JAC[1][1] += DNS[j] * XY[j][1];
            }

            const double detJ = JAC[0][0]*JAC[1][1] - JAC[0][1]*JAC[1][0];
            if ( detJ <= 0.0 )
                throw std::runtime_error("assemble_triangle: non-positive detJ.");

            // Jacobian inverse
            const double JINV[2][2] = {{
                 JAC[1][1]/detJ, -JAC[0][1]/detJ },
               { -JAC[1][0]/detJ,  JAC[0][0]/detJ }
            };

            // Local DN matrix (2x3): B_local(0,j)=DNR[j], B_local(1,j)=DNS[j]
            // Global DN = JINV * B_local  (2x3)
            double dN_dx[2][3];
            for ( size_t j=0; j<nn; ++j ) {
                dN_dx[0][j] = JINV[0][0]*DNR[j] + JINV[0][1]*DNS[j];
                dN_dx[1][j] = JINV[1][0]*DNR[j] + JINV[1][1]*DNS[j];
            }

            // B matrix (3x6) from dN_To2DOF
            std::vector<double> B_mat(3*ud, 0.0);
            auto B = [&](size_t r, size_t c) -> double& { return B_mat[r*ud+c]; };
            for ( size_t i=0; i<nn; ++i ) {
                B(0, 2*i  ) = dN_dx[0][i];   // epsilon_xx
                B(1, 2*i+1) = dN_dx[1][i];   // epsilon_yy
                B(2, 2*i  ) = dN_dx[1][i];   // gamma_xy
                B(2, 2*i+1) = dN_dx[0][i];
            }

            std::vector<double> D_flat(9);
            for (size_t rp=0; rp<3; ++rp)
                for (size_t c=0; c<3; ++c)
                    D_flat[rp*3+c] = D[rp][c];

            auto BT_mat = transpose(B_mat, 3, ud);
            auto DB     = matmul(D_flat, 3, 3, B_mat, ud);
            auto Ke_gp  = matmul(BT_mat, ud, 3, DB, ud);
            for ( size_t k=0; k<ud*ud; ++k ) Ke[k] += w * detJ * Ke_gp[k];

            // BTm = BT * m_vec
            for ( size_t rp=0; rp<ud; ++rp ) {
                double BTm_r = B_mat[0*ud+rp]*m_vec[0]
                             + B_mat[1*ud+rp]*m_vec[1]
                             + B_mat[2*ud+rp]*m_vec[2];
                for ( size_t c=0; c<pd; ++c )
                    Qe[rp*pd+c] += w * detJ * alpha * BTm_r * Np[c];
            }

            for ( size_t rp=0; rp<pd; ++rp )
                for ( size_t c=0; c<pd; ++c )
                    Se[rp*pd+c] += w * detJ * (1.0/M_bio) * Np[rp] * Np[c];
        }

        // Scatter
        const auto uid = u_dofs(nids);
        for ( size_t i=0; i<ud; ++i )
            for ( size_t j=0; j<ud; ++j )
                K(uid[i], uid[j]) += Ke[i*ud+j];
        for ( size_t i=0; i<ud; ++i )
            for ( size_t j=0; j<pd; ++j )
                Q(uid[i], nids[j]) += Qe[i*pd+j];
        for ( size_t i=0; i<pd; ++i )
            for ( size_t j=0; j<pd; ++j )
                S(nids[i], nids[j]) += Se[i*pd+j];
    };


    // -----------------------------------------------------------------------
    // Quad assembly (2x2 Gauss quadrature)
    // -----------------------------------------------------------------------
    auto assemble_quad = [&]( const vector<size_t>& node_ids )
    {
        const int nn = 4;
        const int ud = 2 * nn;   // 8
        const int pd = nn;       // 4

        // Gauss points and weights for 2x2 rule
        const double gp = 1.0 / sqrt(3.0);
        const array<array<double,2>,4> gauss_pts = {{
            {{-gp, -gp}}, {{ gp, -gp}},
            {{-gp,  gp}}, {{ gp,  gp}}
        }};
        constexpr double w = 1.0;   // all weights = 1 for 2x2 rule

        vector<double> Ke(ud*ud, 0.0);
        vector<double> Qe(ud*pd, 0.0);
        vector<double> Se(pd*pd, 0.0);

        for ( const auto& gpt : gauss_pts )
        {
            const double xi  = gpt[0];
            const double eta = gpt[1];

            // Standard bilinear quad shape functions and derivatives
            // for CCW node ordering: (-1,-1), (+1,-1), (+1,+1), (-1,+1)
            // Matches CSMP convention: {1, 4, 5, 2}

            // Shape functions
            const std::array<double,4> N = {{
                0.25*(1.0-xi)*(1.0-eta),   // N0: node at (-1,-1)
                0.25*(1.0+xi)*(1.0-eta),   // N1: node at (+1,-1)
                0.25*(1.0+xi)*(1.0+eta),   // N2: node at (+1,+1)
                0.25*(1.0-xi)*(1.0+eta)    // N3: node at (-1,+1)
            }};

            // Derivatives: dN_dxi[i][0] = dNi/dxi,  dN_dxi[i][1] = dNi/deta
            const std::array<std::array<double,2>,4> dN_dxi = {{
                {{ -0.25*(1.0-eta),  -0.25*(1.0-xi) }},  // dN0/dxi, dN0/deta
                {{  0.25*(1.0-eta),  -0.25*(1.0+xi) }},  // dN1/dxi, dN1/deta
                {{  0.25*(1.0+eta),   0.25*(1.0+xi) }},  // dN2/dxi, dN2/deta
                {{ -0.25*(1.0+eta),   0.25*(1.0-xi) }}   // dN3/dxi, dN3/deta
            }};

            // Jacobian J (2x2): J = dN_dxi^T * x_coords
            double J[2][2] = {{0,0},{0,0}};
            for ( size_t i=0; i<nn; ++i ) {
                J[0][0] += dN_dxi[i][0] * coords[node_ids[i]][0];
                J[0][1] += dN_dxi[i][0] * coords[node_ids[i]][1];
                J[1][0] += dN_dxi[i][1] * coords[node_ids[i]][0];
                J[1][1] += dN_dxi[i][1] * coords[node_ids[i]][1];
            }

            const double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
            if ( detJ <= 0.0 )
                throw runtime_error(
                    "assemble_quad: non-positive Jacobian determinant.");

            // J inverse
            const double Jinv[2][2] = {{
                 J[1][1]/detJ, -J[0][1]/detJ },
               { -J[1][0]/detJ, J[0][0]/detJ }
            };

            // Global shape function derivatives dN_dx (4x2)
            array<array<double,2>,4> dN_dx;
            for ( size_t i=0; i<nn; ++i ) {
                dN_dx[i][0] = Jinv[0][0]*dN_dxi[i][0] + Jinv[0][1]*dN_dxi[i][1];  // row 0 of Jinv
                dN_dx[i][1] = Jinv[1][0]*dN_dxi[i][0] + Jinv[1][1]*dN_dxi[i][1];  // row 1 of Jinv
            }
            // B matrix (3 x 8)
            vector<double> B_mat(3*ud, 0.0);
            auto B = [&](size_t r, size_t c) -> double& { return B_mat[r*ud+c]; };
            for ( size_t i=0; i<nn; ++i ) {
                B(0, 2*i  ) = dN_dx[i][0];   // epsilon_xx
                B(1, 2*i+1) = dN_dx[i][1];   // epsilon_yy
                B(2, 2*i  ) = dN_dx[i][1];   // gamma_xy (symmetric)
                B(2, 2*i+1) = dN_dx[i][0];
            }

            // D matrix flat
            vector<double> D_flat(9);
            for (size_t r=0; r<3; ++r)
                for (size_t c=0; c<3; ++c)
                    D_flat[r*3+c] = D[r][c];

            auto BT_mat = transpose(B_mat, 3, ud);

            // Ke += w * detJ * BT * D * B
            auto DB  = matmul(D_flat, 3, 3, B_mat, ud);
            auto Ke_gp = matmul(BT_mat, ud, 3, DB, ud);
            for ( size_t k=0; k<ud*ud; ++k ) Ke[k] += w * detJ * Ke_gp[k];

            // BTm = BT * m_vec  (8x1)
            vector<double> BTm(ud, 0.0);
            for ( size_t r=0; r<ud; ++r )
                for ( size_t k=0; k<3; ++k )
                    BTm[r] += BT_mat[r*3+k] * m_vec[k];

            // Np row vector (1x4)
            // Qe += w * detJ * alpha * outer(BTm, N)
            for ( size_t r=0; r<ud; ++r )
                for ( size_t c=0; c<pd; ++c )
                    Qe[r*pd+c] += w * detJ * alpha * BTm[r] * N[c];

            // Se += w * detJ * (1/M) * Np^T * Np
            for ( size_t r=0; r<pd; ++r )
                for ( size_t c=0; c<pd; ++c )
                    Se[r*pd+c] += w * detJ * (1.0/M_bio) * N[r] * N[c];

        } // end for gp..

        // Scatter into global matrices
        auto u_idx = u_dofs(node_ids);
        const auto& p_idx = node_ids;

        for ( size_t i=0; i<ud; ++i )
            for ( size_t j=0; j<ud; ++j )
                K(u_idx[i], u_idx[j]) += Ke[i*ud+j];

        for ( size_t i=0; i<ud; ++i )
            for ( size_t j=0; j<pd; ++j )
                Q(u_idx[i], p_idx[j]) += Qe[i*pd+j];

        for ( size_t i=0; i<pd; ++i )
            for ( size_t j=0; j<pd; ++j )
                S(p_idx[i], p_idx[j]) += Se[i*pd+j];
    };


    // -----------------------------------------------------------------------
    // Line element assembly
    // Nodes {1, 3}: (1,1) -> (0,0)
    // -----------------------------------------------------------------------
    auto assemble_line = [&]( const std::vector<size_t>& nids )
    {
        // Replicates CSMP IsoparametricLinearLineElement::dN_AtIntegrationPoint
        // for dim=2, single Gauss point at xi=0, weight=2.
        //
        // Global derivatives are computed as:
        //   cosa = dxdr / detJ,  sina = dydr / detJ
        //   DN(0,i) = DNR[i] * cosa * detJinv   (x-derivative)
        //   DN(1,i) = DNR[i] * sina * detJinv   (y-derivative)
        //
        // Then dN_To2DOF expands DN (2x2) -> B (3x4)

        const int ud = 4;

        const double x0 = coords[nids[0]][0], y0 = coords[nids[0]][1];
        const double x1 = coords[nids[1]][0], y1 = coords[nids[1]][1];

        // Reference derivatives at xi=0
        const double dN0_dxi = -0.5;
        const double dN1_dxi = +0.5;

        // Jacobian components
        const double dxdr  = dN0_dxi*x0 + dN1_dxi*x1;   // (x1-x0)/2
        const double dydr  = dN0_dxi*y0 + dN1_dxi*y1;   // (y1-y0)/2
        const double detJ  = std::hypot(dxdr, dydr);      // L/2
        const double detJinv = 1.0 / detJ;

        if ( detJ <= 0.0 )
            throw std::runtime_error("assemble_line_csmp: zero length element.");

        const double cosa = detJinv * dxdr;   // direction cosine x
        const double sina = detJinv * dydr;   // direction cosine y

        // Global derivatives (matches DN matrix from dN_AtIntegrationPoint)
        const double DN00 = dN0_dxi * cosa * detJinv;   // dN0/dx = +0.5
        const double DN01 = dN1_dxi * cosa * detJinv;   // dN1/dx = -0.5
        const double DN10 = dN0_dxi * sina * detJinv;   // dN0/dy = +0.5
        const double DN11 = dN1_dxi * sina * detJinv;   // dN1/dy = -0.5

        // B matrix (3x4) from dN_To2DOF
        std::vector<double> B_mat(3*ud, 0.0);
        auto B = [&](size_t r, size_t c) -> double& { return B_mat[r*ud+c]; };

        B(0,0) = DN00;  B(0,2) = DN01;          // epsilon_xx
        B(1,1) = DN10;  B(1,3) = DN11;          // epsilon_yy
        B(2,0) = DN10;  B(2,1) = DN00;          // gamma_xy
        B(2,2) = DN11;  B(2,3) = DN01;

        std::vector<double> D_flat(9);
        for (size_t r=0; r<3; ++r)
            for (size_t c=0; c<3; ++c)
                D_flat[r*3+c] = D[r][c];

        auto BT_mat = transpose(B_mat, 3, ud);
        auto DB     = matmul(D_flat, 3, 3, B_mat, ud);
        auto Ke     = matmul(BT_mat, ud, 3, DB, ud);

        // w=2 (single Gauss point), scale by detJ
        const double scale = 2.0 * detJ;
        for ( auto& v : Ke ) v *= scale;

        // Scatter into global K only
        const auto uid = u_dofs(nids);
        for ( size_t i=0; i<ud; ++i )
            for ( size_t j=0; j<ud; ++j )
                K(uid[i], uid[j]) += Ke[i*ud+j];
    };
   
    auto assemble_line_coupling = [&]( const std::vector<size_t>& nids )
    {
        // Line element Biot coupling using CSMP's dN_AtIntegrationPoint convention
        // Single Gauss point at xi=0, weight=2, detJ=L/2

        const int ud = 4;   // 2 nodes x 2 DOFs
        const int pd = 2;   // 2 pressure DOFs (one per node)

        const double x0 = coords[nids[0]][0], y0 = coords[nids[0]][1];
        const double x1 = coords[nids[1]][0], y1 = coords[nids[1]][1];

        const double dN0_dxi = -0.5,  dN1_dxi = +0.5;
        const double dxdr  = dN0_dxi*x0 + dN1_dxi*x1;
        const double dydr  = dN0_dxi*y0 + dN1_dxi*y1;
        const double detJ  = std::hypot(dxdr, dydr);
        const double detJinv = 1.0 / detJ;
        const double cosa  = detJinv * dxdr;
        const double sina  = detJinv * dydr;

        // Global derivatives from dN_AtIntegrationPoint (dim=2)
        const double DN00 = dN0_dxi * cosa * detJinv;   // dN0/dx
        const double DN01 = dN1_dxi * cosa * detJinv;   // dN1/dx
        const double DN10 = dN0_dxi * sina * detJinv;   // dN0/dy
        const double DN11 = dN1_dxi * sina * detJinv;   // dN1/dy

        // B matrix (3x4) from dN_To2DOF
        // Row 0: [DN00,  0,    DN01,  0   ]
        // Row 1: [ 0,   DN10,   0,   DN11 ]
        // Row 2: [DN10, DN00, DN11, DN01  ]
        const std::array<double,3*4> B_mat = {
            DN00,  0.0,  DN01,  0.0,
            0.0,  DN10,  0.0,  DN11,
            DN10, DN00, DN11, DN01
        };

        // Np at xi=0: both nodes = 0.5
        const std::array<double,2> Np = { 0.5, 0.5 };

        // scale = w * detJ * alpha = 2 * detJ * 1
        const double scale = 2.0 * detJ;

        // Qe (4x2): BTm outer product Np
        // Using compression-positive m = [1,1,0]^T to match CSMP
        std::vector<double> Qe(ud*pd, 0.0);
        for ( size_t r=0; r<ud; ++r )
        {
            // Use m_vec consistently — same as triangle and quad assemblers
            const double BTm_r = scale * ( B_mat[0*4+r] * m_vec[0]
                                         + B_mat[1*4+r] * m_vec[1]
                                         + B_mat[2*4+r] * m_vec[2] );
            for ( size_t c=0; c<pd; ++c )
                Qe[r*pd+c] = BTm_r * Np[c];
        }

        // Scatter into global Q
        const auto uid = u_dofs(nids);
        for ( size_t i=0; i<ud; ++i )
            for ( size_t j=0; j<pd; ++j )
                Q(uid[i], nids[j]) += Qe[i*pd+j];
    };

    // -----------------------------------------------------------------------
    // Run element assembly
    // -----------------------------------------------------------------------
    // node numbering following CSMP conventions
    assemble_triangle({ 0, 3, 1 });    // CCW and starting with the lowest node number
    assemble_triangle({ 1, 3, 4 });
    assemble_quad    ({ 1, 4, 5, 2 });
    assemble_line({ 1, 3 });
    assemble_line_coupling({ 1, 3 });  // pressure coupling that was originally ignored
    

    // -----------------------------------------------------------------------
    // Build full 18x18 block system A_full
    // -----------------------------------------------------------------------
    vector<double> A_full(N_TOTAL * N_TOTAL, 0.0);
    auto A = [&](size_t r, size_t c) -> double& { return A_full[r*N_TOTAL+c]; };

    // K block  [0..11, 0..11]
    for ( size_t r=0; r<N_DISP; ++r )
        for ( size_t c=0; c<N_DISP; ++c )
            A(r, c) = K_gl[r*N_DISP+c];

    // Q block  [0..11, 12..17]
    for ( size_t r=0; r<N_DISP; ++r )
        for ( size_t c=0; c<N_PRES; ++c )
            A(r, N_DISP+c) = Q_gl[r*N_PRES+c];

    // Q^T block [12..17, 0..11]
    for ( size_t r=0; r<N_PRES; ++r )
        for ( size_t c=0; c<N_DISP; ++c )
            A(N_DISP+r, c) = Q_gl[c*N_PRES+r];

    // S block  [12..17, 12..17]
    for ( size_t r=0; r<N_PRES; ++r )
        for ( size_t c=0; c<N_PRES; ++c )
            A(N_DISP+r, N_DISP+c) = S_gl[r*N_PRES+c];

    // see function definition above
    if ( verbose_ )
      diagnosePrintAssembly(  K_gl, Q_gl, S_gl, N_DISP, N_PRES );


    // -----------------------------------------------------------------------
    // Dirichlet BCs
    // -----------------------------------------------------------------------
    const std::map<size_t,double> dirichlet = {
        {  0,   0.0  },   // u_x node 0  (Left)
        {  1,   0.0  },   // u_y node 0  (Left)
        {  4,   0.01 },   // u_x node 2  (Right)
        {  6,   0.0  },   // u_x node 3  (Left)
        {  7,   0.0  },   // u_y node 3  (Left)
        {  9,   0.0  },   // u_y node 4  (Bottom)
        { 10,   0.01 },   // u_x node 5  (Right)
        { 11,   0.0  },   // u_y node 5  (Bottom+Right)
        { 12,   0.0  },   // p   node 0  (Left)
        { 13,   0.0  },   // p   node 1  (Top — no conflict, interior top)
        { 14, 500.0  },   // p   node 2  (Right dominates over Top)
        { 15,   0.0  },   // p   node 3  (Left)
        { 17, 500.0  }    // p   node 5  (Right)
    };

    // Build prescribed value vector x_d
    vector<double> x_d(N_TOTAL, 0.0);
    for ( const auto& [idx, val] : dirichlet )
        x_d[idx] = val;

    // Identify free DOFs
    set<size_t> fixed_set;
    for ( const auto& [idx, val] : dirichlet )
        fixed_set.insert(idx);

    vector<size_t> free_dofs;
    free_dofs.reserve(N_TOTAL - fixed_set.size());
    for ( size_t i=0; i<N_TOTAL; ++i )
        if ( !fixed_set.contains(i) )
            free_dofs.push_back(i);

    // free_dofs should be: {2, 3, 5, 8, 16}
    const size_t n_free = free_dofs.size();

    // -----------------------------------------------------------------------
    // Modify RHS: b_mod = b_full - A_full * x_d
    // (b_full is zero here — no body forces or fluid sources)
    // -----------------------------------------------------------------------
    vector<double> b_mod(N_TOTAL, 0.0);
    for ( size_t r=0; r<N_TOTAL; ++r )
        for ( size_t c=0; c<N_TOTAL; ++c )
            b_mod[r] -= A(r,c) * x_d[c];

    // -----------------------------------------------------------------------
    // Extract reduced system into output arguments
    // -----------------------------------------------------------------------
    b_red.assign(n_free, 0.0);
    for ( size_t i=0; i<n_free; ++i )
        b_red[i] = b_mod[free_dofs[i]];

    // Populate CompressedRowMatrix
    // First pass: collect all non-zero (i,j) entries
    vector<double> A_dense(n_free * n_free, 0.0);
    
    for ( size_t i=0; i<n_free; ++i )
    {
        b_red[i] = b_mod[free_dofs[i]];
        for ( size_t j=0; j<n_free; ++j )
            A_dense[i*n_free+j] = A_full[free_dofs[i]*N_TOTAL + free_dofs[j]];
    }
    
    A_red = makeCompressedRowMatrix(A_dense, n_free);
            
} // end GenerateTinyPoromechanicsTestSystem



/**
 * @brief Prints the reduced matrix and RHS to stdout for visual inspection.
 *        Call after GenerateTinyPoromechanicsTestSystem().
 */
void PDE_Integrator_Test::PrintTinyTestSystem( const CompressedRowMatrix& A_red,
                                               const vector<double>& b_red )
{
    // Free DOF labels — useful when reading output
    constexpr array<const char*, 5> labels = {{
        "u_x(1)", "u_y(1)", "u_y(2)", "u_x(4)", "p(4)"
    }};

    cout << scientific << setprecision(6);

    // The pressure diagonal must be strictly positive
    if ( A_red(4,4) <= 0.0 ) {
        cout << scientific << setprecision(6);
        cout <<"Error: Pressure diagonal S(4,4) must be positive: " << A_red(4,4);
      }

    cout << "\nFree DOFs:\n";
    for ( size_t i=0; i<5; ++i )
        cout << "  [" << i << "] global "
                  << array<int,5>{2,3,5,8,16}[i]
                  << "  " << labels[i] << "\n";

    cout << "\nReduced 5x5 matrix A_red:\n";
    for ( size_t i=0; i<5; ++i ) {
        for ( size_t j=0; j<5; ++j )
            cout << setw(14) << scientific << fixed << setprecision(4) << A_red(i,j);
        cout << "\n";
    }

    cout << "\nReduced RHS b_red:\n";
    for ( size_t i=0; i<5; ++i )
        cout << "  [" << i << "] " << scientific << fixed << setprecision(4) << b_red[i]
                  << "  (" << labels[i] << ")\n";
}





} // end namespace csmp


