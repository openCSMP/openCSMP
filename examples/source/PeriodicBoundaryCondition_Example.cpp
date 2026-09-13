#include "PeriodicBoundaryCondition_Example.h"

#include "Model.h"
#include "PDE_Integrator.h"

// fluid pressure algorithm and velocity computation
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "EigenSolver.h"

// interfaces
#include "VTU_Interface.h"

// utilities
#include "CSMP_highLevelUtilities.h"

#include <algorithm>
#include <vector>

using namespace std;

namespace csmp {

void PeriodicBoundaryCondition_Example::Specifications()
{
   SetTitle( "Steady-state Darcy flow with a periodic boundary condition" );
   SetDifficulty( 3 );
   SetCategory( "Simulation of Physical Processes" );
   AddAuthor( "SKM" );
   AddDescription( "2D steady-state fluid pressure distribution through a fracture network" );
   AddDescription( "choice at runtime between a single periodic direction (periodic top/bottom, "
                   "Dirichlet left/right) and doubly-periodic (both pairs of edges periodic, "
                   "with a pressure jump driving flow and the four corners Dirichlet-pinned)" );
   AddDescription( "demonstrates PDE_Integrator's periodic ('master_node_id'-linked) DOF handling" );
   AddDescription( "source in: PeriodicBoundaryCondition_Example.cpp" );
   AddRequirement( "file set: 'periodic_fracture_network'" );
}



/**
    The "master_node_id" property imported with a CSMP native model is only ever
    meaningful for nodes that a previous step actually flagged PERIODIC — and even
    then its values are node *indices*, which are only valid for the exact node
    numbering that was in effect when they were computed. Copying such a dataset
    between builds (or even just across a mesh re-numbering) invalidates them,
    since node numbering is not guaranteed to be portable.

    This helper rebuilds the periodic pairing for one pair of opposite box
    boundaries at runtime instead of trusting whatever is stored on disk: it
    matches boundary1 nodes to boundary2 nodes by their shared coordinate (x
    for TOP/BOTTOM, y for LEFT/RIGHT), then designates boundary1 as the anchor
    of each pair (self-referencing its own Idx()) and boundary2 as the slave
    that points at it. Both boundaries are then flagged PERIODIC, with
    'pressure_jump' as the nominal prescribed value (0.0 for pure periodicity,
    nonzero to drive a net pressure difference across the pair while keeping
    fluctuations periodic).

    @attention PDE_Integrator::EliminateEssentialConditions() only allocates a
    real degree of freedom for a periodic pair when exactly one side satisfies
    position == master_position (its "master_node_id" resolves to its own
    Idx()). Writing a *mutual* reference — both sides pointing at each other —
    leaves neither side its own master, so neither ever gets a real DOF: both
    stay eliminated (as if Dirichlet) and the periodic coupling has no effect.
*/
static void SetPeriodicBoundaryConditions( Model<2U>& model, BOX_BOUNDARY boundary1, BOX_BOUNDARY boundary2, double pressure_jump )
{
    auto& region = model.Region( "Model" );
    region.RenumberNodes(); // fixes Idx() now, so it matches what Apply() will reuse

    const Index mkey = model.Database().StorageKey( "master_node_id" );

    vector<Node<2U>*> nodes1, nodes2;
    for ( auto it = region.PerimeterNodesBegin(); it != region.NodesEnd(); ++it ) {
        if      ( (*it)->AtBoundary() == boundary1 ) nodes1.push_back( *it );
        else if ( (*it)->AtBoundary() == boundary2 ) nodes2.push_back( *it );
    }

    if ( nodes1.size() != nodes2.size() )
        throw csmp::Exception( ERROR, "PeriodicBoundaryCondition_Example::SetPeriodicBoundaryConditions",
            "Periodic node counts differ; cannot pair periodic nodes." );

    // TOP/BOTTOM pairs share an x-coordinate; LEFT/RIGHT pairs share a y-coordinate.
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



/// Periodic top/bottom only; left/right remain ordinary Dirichlet edges.
static void SetSinglePeriodicBoundaryCondition( Model<2U>& model )
{
    SetPeriodicBoundaryConditions( model, BOTTOM, TOP, 0.0 ); // pure periodicity in y

    model.InputBoundaryValue( LEFT,  "fluid pressure", makeScalar(DIRICH,1.0) );
    model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,0.0) );
}

/**
    Both pairs of model edges periodic: a pressure jump of -1.0 is imposed
    across LEFT/RIGHT to drive net flow through the otherwise-periodic
    domain, while TOP/BOTTOM remains purely periodic (no jump). Since a fully
    periodic domain has no boundary left to anchor the absolute pressure
    level, the four model corners are pinned with Dirichlet values instead.
*/
static void SetDoublePeriodicBoundaryCondition( Model<2U>& model )
{
    SetPeriodicBoundaryConditions( model, BOTTOM, TOP,   0.0 ); // pure periodicity in y
    SetPeriodicBoundaryConditions( model, LEFT,   RIGHT, -1.0 ); // pressure jump drives flow in x

    model.InputBoundaryValue( CNR1, "fluid pressure", makeScalar(DIRICH,1.0) );
    model.InputBoundaryValue( CNR2, "fluid pressure", makeScalar(DIRICH,0.0) );
    model.InputBoundaryValue( CNR3, "fluid pressure", makeScalar(DIRICH,0.0) );
    model.InputBoundaryValue( CNR4, "fluid pressure", makeScalar(DIRICH,1.0) );
}



/**
    Steady-state 2D Darcy flow through a heterogeneous fracture network, with
    a choice — asked interactively — between two periodic boundary condition
    setups:

      1) single periodic direction: TOP/BOTTOM periodic, LEFT/RIGHT Dirichlet;
      2) double periodic: both pairs of edges periodic, with a pressure jump
         across LEFT/RIGHT driving net flow and the four corners
         Dirichlet-pinned to anchor the absolute pressure level (a fully
         periodic domain has no boundary left to do this otherwise).

    Illustrates:
      - how to flag a boundary PERIODIC via Model::InputBoundaryValue(), tying
        pairs of nodes to the same degree of freedom instead of eliminating
        them like a Dirichlet condition;
      - PDE_Integrator's need for a valid "master_node_id" node-pairing (see
        SetPeriodicBoundaryConditions() above for why this can't simply be
        reused from an imported dataset without rebuilding it for the
        current run).

    Use ParaView to visualize the VTU output; fluid pressure should match
    exactly between every pair of periodic nodes at the same shared
    coordinate (allowing for the imposed pressure jump on LEFT/RIGHT, in the
    double-periodic case).
*/
void PeriodicBoundaryCondition_Example::Run()
{
    // 1.  Load the fracture-network model
    // ------------------------------------
    const string model_name = "periodic_fracture_network";
    const string file_name  = GetExampleFileName( __FILE__ );
    CreateWorkingDirectoryAndCopyInputModelFiles( file_name, model_name, "" );
    Model<2U>  model( model_name );

    printModelDimensions( model, true );


    // 2.  Ensure the righthand-side source term has a well-defined value
    //     everywhere (it is otherwise left unset in the imported dataset)
    // ---------------------------------------------------------------------
    model.InputPropertyValue( "nodal fluid volume source", makeScalar(ANY,0.0) );


    // 3.  Let the user choose the boundary condition setup
    // -------------------------------------------------------
    int choice{0};
    while ( choice != 1 && choice != 2 ) {
        cout << "\nChoose the boundary condition setup:\n";
        cout << "  1) single periodic direction (periodic top/bottom, Dirichlet left/right)\n";
        cout << "  2) double periodic (both pairs of edges periodic, pressure-jump-driven, corners pinned)\n";
        cout << "Enter choice [1 or 2]: ";
        cin >> choice;
        if ( choice != 1 && choice != 2 )
            cout << "Invalid choice, please try again.\n";
    }

    if ( choice == 1 ) SetSinglePeriodicBoundaryCondition( model );
    else                SetDoublePeriodicBoundaryCondition( model );


    // 4.  Building the steady-state FE algorithm "darcy_equation"
    // ---------------------------------------------------------------
    EigenSolver  solver;
    PDE_Integrator<2U,Element>  darcy_equation( solver );

    NumIntegral_dNT_lhsop_dN_dV<2U> darcyTerm( model.Database(), "conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<2U>   source( model.Database(), "nodal fluid volume source", "fluid pressure" );
    VelocityAndVolumeFlux<2U>       flow( model, "tensor conductivity", "porosity", "fluid pressure", true );

    darcy_equation.Add( &darcyTerm );
    darcy_equation.Add( &source );
    darcy_equation.AddPostProcess( &flow );


    // 5.  Solve
    // ---------
    model.Apply( darcy_equation );


    // 6.  Report and output
    // ----------------------
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "velocity" );
    printRangeOfVariable( model, "volume flux" );

    VTU_Interface<2U>  vtu_output{ model };
    vtu_output.OutputDataToVTU( "fluid-pressure",
        list<string>{
            "master_node_id",
            "permeability",
            "porosity",
            "fluid pressure",
            "velocity",
            "pore velocity",
            "volume flux",
        }, "Model", 0 );

    cout << "\nPeriodicBoundaryCondition_Example: done. Fluid pressure should match exactly\n";
    cout << "between every periodic node pair at the same shared coordinate (allowing for\n";
    cout << "the imposed pressure jump on the LEFT/RIGHT pair).\n";

    filesystem::current_path( "../../example_inputs/" );

} // end Run

} // csmp
