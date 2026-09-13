#include "PeriodicBoundaryCondition_Test.h"

#include "Model.h"
#include "PDE_Integrator.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "EigenSolver.h"
#include "Exception.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

namespace csmp {

namespace {

/**
    Rebuilds the periodic pairing for one pair of opposite box boundaries at
    runtime (see PeriodicBoundaryCondition_Example.cpp for the full
    explanation of why the imported "master_node_id" values can't be trusted
    as-is): matches boundary1 nodes to boundary2 nodes by their shared
    coordinate, designates boundary1 as the self-referencing anchor of each
    pair and boundary2 as the slave, then flags both PERIODIC.
*/
void SetPeriodicBoundaryConditions( Model<2U>& model, BOX_BOUNDARY boundary1, BOX_BOUNDARY boundary2, double pressure_jump )
{
    auto& region = model.Region( "Model" );
    region.RenumberNodes();

    const Index mkey = model.Database().StorageKey( "master_node_id" );

    vector<Node<2U>*> nodes1, nodes2;
    for ( auto it = region.PerimeterNodesBegin(); it != region.NodesEnd(); ++it ) {
        if      ( (*it)->AtBoundary() == boundary1 ) nodes1.push_back( *it );
        else if ( (*it)->AtBoundary() == boundary2 ) nodes2.push_back( *it );
    }

    if ( nodes1.size() != nodes2.size() )
        throw csmp::Exception( ERROR, "PeriodicBoundaryCondition_Test::SetPeriodicBoundaryConditions",
            "Periodic node counts differ; cannot pair periodic nodes." );

    auto by_coord = [boundary1]( Node<2U>* a, Node<2U>* b ) {
        return (boundary1 == TOP || boundary1 == BOTTOM) ?
            a->Coordinate()[0] < b->Coordinate()[0] : a->Coordinate()[1] < b->Coordinate()[1];
    };
    sort( nodes1.begin(), nodes1.end(), by_coord );
    sort( nodes2.begin(), nodes2.end(), by_coord );

    for ( size_t i{0U}; i < nodes1.size(); ++i ) {
        nodes1[i]->Store( mkey, makeScalar( ANY, static_cast<double>(nodes1[i]->Idx()) ) );
        nodes2[i]->Store( mkey, makeScalar( ANY, static_cast<double>(nodes1[i]->Idx()) ) );
    }

    model.InputBoundaryValue( boundary1, "fluid pressure", makeScalar(PERIODIC,pressure_jump) );
    model.InputBoundaryValue( boundary2, "fluid pressure", makeScalar(PERIODIC,pressure_jump) );
}

/// Collects (shared coordinate, solved fluid pressure) for every node on 'boundary', sorted by coordinate.
vector<pair<double,double>> BoundaryPressures( Model<2U>& model, BOX_BOUNDARY boundary )
{
    auto& region = model.Region( "Model" );
    const Index pkey = model.Database().StorageKey( "fluid pressure" );
    const bool use_x = (boundary == TOP || boundary == BOTTOM);

    vector<pair<double,double>> values;
    for ( auto it = region.PerimeterNodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->AtBoundary() == boundary ) {
            const double coord = use_x ? (*it)->Coordinate()[0] : (*it)->Coordinate()[1];
            values.emplace_back( coord, (*it)->Read(pkey) );
        }
    sort( values.begin(), values.end() );
    return values;
}

} // anonymous namespace



/**
    Periodic top/bottom only; left/right remain ordinary Dirichlet edges.
    Every top/bottom node pair at the same x-coordinate must end up with
    exactly the same solved fluid pressure.
*/
void PeriodicBoundaryCondition_Test::TestSinglePeriodicDirection()
{
    const string model_name = "periodic_fracture_network"; // std::string: avoids the protected Model(const char*) overload
    Model<2U> model( model_name );
    model.InputPropertyValue( "nodal fluid volume source", makeScalar(ANY,0.0) );

    SetPeriodicBoundaryConditions( model, BOTTOM, TOP, 0.0 );
    model.InputBoundaryValue( LEFT,  "fluid pressure", makeScalar(DIRICH,1.0) );
    model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,0.0) );

    EigenSolver  solver;
    PDE_Integrator<2U,Element>  darcy_equation( solver );
    NumIntegral_dNT_lhsop_dN_dV<2U> darcyTerm( model.Database(), "conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<2U>   source( model.Database(), "nodal fluid volume source", "fluid pressure" );
    darcy_equation.Add( &darcyTerm );
    darcy_equation.Add( &source );
    model.Apply( darcy_equation );

    const auto bottom = BoundaryPressures( model, BOTTOM );
    const auto top    = BoundaryPressures( model, TOP );
    _test( bottom.size() == top.size() );

    double max_diff{0.0};
    for ( size_t i{0U}; i < bottom.size() && i < top.size(); ++i ) {
        _test( isfinite(bottom[i].second) && isfinite(top[i].second) );
        max_diff = max( max_diff, fabs(bottom[i].second - top[i].second) );
    }
    // pure periodicity (no jump): top and bottom must match exactly (to solver tolerance)
    _equal( max_diff, 0.0, 1.0e-9 );
}



/**
    Both pairs of model edges periodic: a pressure jump is imposed across
    LEFT/RIGHT to drive net flow, while TOP/BOTTOM remains purely periodic.
    The four corners are Dirichlet-pinned, since a fully periodic domain has
    no boundary left to anchor the absolute pressure level otherwise.

    Every top/bottom pair must match exactly (no jump). Every left/right pair
    must differ by the *same* constant amount everywhere — a non-uniform
    difference across pairs would indicate the periodic node pairing itself
    is wrong (e.g. mismatched nodes, or a mutual self-reference that leaves
    both sides eliminated instead of tied to one master DOF).
*/
void PeriodicBoundaryCondition_Test::TestDoublePeriodicDirections()
{
    const double pressure_jump{-1.0};

    const string model_name = "periodic_fracture_network"; // std::string: avoids the protected Model(const char*) overload
    Model<2U> model( model_name );
    model.InputPropertyValue( "nodal fluid volume source", makeScalar(ANY,0.0) );

    SetPeriodicBoundaryConditions( model, BOTTOM, TOP,   0.0 );
    SetPeriodicBoundaryConditions( model, LEFT,   RIGHT, pressure_jump );

    model.InputBoundaryValue( CNR1, "fluid pressure", makeScalar(DIRICH,1.0) );
    model.InputBoundaryValue( CNR2, "fluid pressure", makeScalar(DIRICH,0.0) );
    model.InputBoundaryValue( CNR3, "fluid pressure", makeScalar(DIRICH,0.0) );
    model.InputBoundaryValue( CNR4, "fluid pressure", makeScalar(DIRICH,1.0) );

    EigenSolver  solver;
    PDE_Integrator<2U,Element>  darcy_equation( solver );
    NumIntegral_dNT_lhsop_dN_dV<2U> darcyTerm( model.Database(), "conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<2U>   source( model.Database(), "nodal fluid volume source", "fluid pressure" );
    darcy_equation.Add( &darcyTerm );
    darcy_equation.Add( &source );
    model.Apply( darcy_equation );

    // TOP/BOTTOM: no jump, must match exactly
    const auto bottom = BoundaryPressures( model, BOTTOM );
    const auto top    = BoundaryPressures( model, TOP );
    _test( bottom.size() == top.size() );
    double max_diff{0.0};
    for ( size_t i{0U}; i < bottom.size() && i < top.size(); ++i ) {
        _test( isfinite(bottom[i].second) && isfinite(top[i].second) );
        max_diff = max( max_diff, fabs(bottom[i].second - top[i].second) );
    }
    _equal( max_diff, 0.0, 1.0e-9 );

    // LEFT/RIGHT: jump-driven, but the offset must be the *same* for every pair
    const auto left  = BoundaryPressures( model, LEFT );
    const auto right = BoundaryPressures( model, RIGHT );
    _test( left.size() == right.size() );
    _test( !left.empty() );
    if ( !left.empty() ) {
        const double reference_offset = left.front().second - right.front().second;
        double max_offset_deviation{0.0};
        for ( size_t i{0U}; i < left.size() && i < right.size(); ++i ) {
            _test( isfinite(left[i].second) && isfinite(right[i].second) );
            const double offset = left[i].second - right[i].second;
            max_offset_deviation = max( max_offset_deviation, fabs(offset - reference_offset) );
        }
        _equal( max_offset_deviation, 0.0, 1.0e-9 );
    }
}



void PeriodicBoundaryCondition_Test::run()
{
    TestSinglePeriodicDirection();
    TestDoublePeriodicDirections();
}

} // csmp
