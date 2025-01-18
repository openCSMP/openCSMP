#include "NimbleRegion_Test.h"

#include "PropertyDatabase.h"
#include "NimbleRegion.h"
#include "Region.h"

// FE algorithm
#include "PDE_Integrator.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "SAMG_Exception.h"
#else
#include "LinearSolver.h"
#endif

#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "compute2PhaseMobilityAtBaryCenter.h"

// FV algorithms
#include "TwoPhaseExplicitNodeCenteredFVTransport.h"
#include "ExplicitStencilProcessor.h"

// relative permeability calculations
#include "BrooksCorey.h"

// interfaces
#include "InputDataManager.h"
#include "VTU_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"
#include "VTK_Interface.h"
#include "Quadrilaterator.h"

using namespace std;

namespace csmp {

NimbleRegion_Test::NimbleRegion_Test( const string& model, const string& variables_file )
 : model2D_(nullptr)
 {
    Quadrilaterator    quadrilaterator; // simple FE mesher
    VSet<2U>           mesh_container;  // container to store the input mesh
    const string       file_name("tutorial1_input");

    quadrilaterator.QuadrilateralsFromRegularGrid(mesh_container, model.c_str(), 100., 100. );

    model2D_ = new Model<2U>( mesh_container, variables_file.c_str() );

    model2D_->InputPropertyValue("permeability", makeScalar(PLAIN, 1e-12));
    model2D_->InputPropertyValue("porosity", makeScalar(PLAIN, 0.25));
    model2D_->InputPropertyValue("fluid volume source", makeScalar(PLAIN, 0.0));
    model2D_->InputPropertyValue("fluid pressure", makeScalar(PLAIN,0.0));
    model2D_->InputPropertyValue("saturation oil", makeScalar(PLAIN, 0.0));
    model2D_->InputPropertyValue("saturation water", makeScalar(PLAIN, 1.0));
    model2D_->InputPropertyValue("viscosity oil", makeScalar(PLAIN, 0.001));
    model2D_->InputPropertyValue("viscosity water", makeScalar(PLAIN, 0.001));
    model2D_->InputPropertyValue("density oil", makeScalar(PLAIN, 1000.));
    model2D_->InputPropertyValue("density water", makeScalar(PLAIN, 1000.));
    model2D_->InputPropertyValue("residual saturation wetting phase", makeScalar(PLAIN, 0.));
    model2D_->InputPropertyValue("residual saturation non-wetting phase", makeScalar(PLAIN, 0.));
    model2D_->InputPropertyValue("brooks corey parameter", makeScalar(PLAIN, 3.));
    model2D_->InputPropertyValue("entry pressure", makeScalar(PLAIN, 1.e-20));
    model2D_->InputPropertyValue("total system compressibility", makeScalar(PLAIN, 1.0e-9));
    model2D_->InputPropertyValue("nodal fluid volume source", makeScalar(PLAIN, 0.));
    model2D_->InputPropertyValue("velocity", makeVector(PLAIN, PLAIN, 0.0, 0.));
    model2D_->InputPropertyValue("pore velocity", makeVector(PLAIN, PLAIN, 0.0, 0.0));
    model2D_->InputPropertyValue("nodal velocity", makeVector(PLAIN, PLAIN, 0.0, 0.0));
    model2D_->InputPropertyValue("nodal pore velocity", makeVector(PLAIN, PLAIN, 0.0, 0.0));
    model2D_->InputPropertyValue("nodal volume flux", makeScalar(PLAIN, 0.0));
   
    model2D_->InputBoundaryValue(LEFT, "fluid pressure", makeScalar(DIRICH, 0.));
    model2D_->InputBoundaryValue(RIGHT, "fluid pressure", makeScalar(DIRICH, 0.));
    model2D_->InputBoundaryValue(LEFT, "saturation oil", makeScalar(DIRICH, 0.));
    model2D_->InputBoundaryValue(LEFT, "saturation oil", makeScalar(DIRICH, 0.));

    // creating a 10x10m region in the center of the 2D model, which is initialised with elevated non-wetting phase saturation
    model2D_->FormRectangularRegion("central", Point<2U>(45., 45.), Point<2U>(55., 55.));
    Region<2U>& central = model2D_->Region("central");
    central.InputPropertyValue("fluid pressure", makeScalar(DIRICH, 2.e5));
    central.InputPropertyValue("saturation oil", makeScalar(DIRICH, 0.9));
    central.InputPropertyValue("saturation water", makeScalar(DIRICH, 0.1));
   
   // flagging the nimble region nodes and elements as 1, which is a value that gets incremented as the nimble egion evolves
    double nimble_val(1.);
    model2D_->InputPropertyValue("nimble nodes", makeScalar(PLAIN, 0.));
    model2D_->InputPropertyValue("nimble elements", makeScalar(PLAIN, 0.));
    central.InputPropertyValue("nimble nodes", makeScalar(PLAIN, nimble_val));
    central.InputPropertyValue("nimble elements", makeScalar(PLAIN, nimble_val));

} // end constructor





void NimbleRegion_Test::run()
	{
		// -------------------------------------------------------------------------
		// 1. creating a flexible NimbleRegion that moves with the saturation front
		// -------------------------------------------------------------------------
    Region<2U>&  model_domain(model2D_->Region("Model"));
    Region<2U>&  sub_domain(model2D_->Region("central"));
    cout <<"\nrun: creating NimbleRegion object from region '"<< sub_domain.Name() <<"' with "<< sub_domain.Nodes() <<" nodes and "<< sub_domain.Cells() <<" elements.\n";
		NimbleRegion<2U>  plume_region( model2D_->Database(), sub_domain.NodesBegin(), sub_domain.NodesEnd() ); // will have a halo of one element extra
    plume_region.Out();
    // tested: O.K.
    // TODO: wrap into automatic test
    
    // first visual check
    VTU_Interface<2U>  vtu(*model2D_);
    // augmenting the extra perimeter nodes by 2, to make sure they are visible
    const csmp::Index  nn_key(model2D_->Database().StorageKey("nimble nodes"));
    const csmp::Index  ne_key(model2D_->Database().StorageKey("nimble elements"));
    
    for ( auto nit=plume_region.PerimeterNodesBegin(); nit!=plume_region.NodesEnd(); ++nit ) {
         double nval = (*nit)->Read( nn_key ) + 2.;
         (*nit)->Store( nn_key, makeScalar(PLAIN,nval) );
      }
    vtu.OutputDataToVTU("nimble-nodes", "nimble nodes", "Model", 0 );


    // -------------------------------------------------------------------------
    // 2. recreating same NimbleRegion
    // -------------------------------------------------------------------------
    plume_region.Update( sub_domain.NodesBegin(), sub_domain.NodesEnd() );
    // tested: O.K.
    // TODO: wrap into automatic test


    // -------------------------------------------------------------------------
    // 3. building a new Nimble region just from every 5th node
    // -------------------------------------------------------------------------
    vector<Node<2U>*>  every25thNode;
    every25thNode.reserve( model_domain.Nodes() / 5 );
    size_t node(0U);
    while ( node < model_domain.Nodes() )
      {
         every25thNode.push_back( model_domain.N(node) );
         node += 25U;
      }
    plume_region.Update( every25thNode.begin(), every25thNode.end() );
    plume_region.Out();
    // tested: O.K.
    // TODO: wrap into automatic test

    // creating a footprint of the revised region
    model2D_->InputPropertyValue("nimble elements", makeScalar(PLAIN, 0.));
    for ( auto it=plume_region.CellsBegin(); it!=plume_region.CellsEnd(); ++it ) {
         (*it)->Store( ne_key, makeScalar(PLAIN,1.) );
      }
    model2D_->InputPropertyValue("nimble nodes", makeScalar(PLAIN, 0.));
    // OK - interior nodes are correctly detected
    for ( auto nit=plume_region.NodesBegin(); nit!=plume_region.PerimeterNodesBegin(); ++nit ) {
          (*nit)->Store( nn_key, makeScalar(PLAIN,1.) );
      }
    for ( auto nit=plume_region.PerimeterNodesBegin(); nit!=plume_region.NodesEnd(); ++nit ) {
         double nval = (*nit)->Read( nn_key ) + 2.;
         (*nit)->Store( nn_key, makeScalar(PLAIN,nval) );
      }
    vtu.OutputDataToVTU("nimble-nodes", "nimble nodes", "Model", 1 );
    vtu.OutputDataToVTU("nimble-elements", "nimble elements", "Model", 1 );


    // The idea is that the NimbleRegion receives a vector of nodes that it consists of and rebuilds
    // itself such that it contains all the elements that contain at least one of these nodes
    cout <<"\nrun: more testing to be implemented here.\n";

} // end run




} // csmp
