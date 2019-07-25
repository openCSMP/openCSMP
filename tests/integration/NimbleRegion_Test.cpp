#include "NimbleRegion_Test.h"

// FE algorithm
#include "PDE_Integrator_UoM.h"
#include "SAMG_Solver.h"
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

void NimbleRegion_Test::run()
	{
		// --------------------------------------------
		// 1. Configure simulation model from bitmap file
		// --------------------------------------------
		Quadrilaterator    quadrilaterator; // simple FE mesher
		VSet<2U>           mesh_container;  // container to store the input mesh
		const string       file_name("tutorial1_input");
		double64           x(100.); double64 y(100.);
		quadrilaterator.QuadrilateralsFromRegularGrid(mesh_container, file_name.c_str(), x, y);
		Model<2U>  model(mesh_container, "NimbleRegion_Test-variables.txt", true);
		const PropertyDatabase<2>& p_ref(model.Database());  // constand reference to the property database

		model.InputPropertyValue("permeability", makeScalar(PLAIN, 1e-12));
		model.InputPropertyValue("porosity", makeScalar(PLAIN, 0.25));
		model.InputPropertyValue("fluid volume source", makeScalar(PLAIN, 0.0));
		model.InputPropertyValue("fluid pressure", makeScalar(PLAIN,0.0));
		model.InputPropertyValue("saturation oil", makeScalar(PLAIN, 0.0));
		model.InputPropertyValue("saturation water", makeScalar(PLAIN, 1.0));
		model.InputPropertyValue("viscosity oil", makeScalar(PLAIN, 0.001));
		model.InputPropertyValue("viscosity water", makeScalar(PLAIN, 0.001));
		model.InputPropertyValue("density oil", makeScalar(PLAIN, 1000.));
		model.InputPropertyValue("density water", makeScalar(PLAIN, 1000.));
		model.InputPropertyValue("residual saturation wetting phase", makeScalar(PLAIN, 0.));
		model.InputPropertyValue("residual saturation non-wetting phase", makeScalar(PLAIN, 0.));
		model.InputPropertyValue("brooks corey parameter", makeScalar(PLAIN, 3.));
		model.InputPropertyValue("entry pressure", makeScalar(PLAIN, 1.e-20));
		model.InputPropertyValue("total system compressibility", makeScalar(PLAIN, 1.0e-9));
		model.InputPropertyValue("nodal fluid volume source", makeScalar(PLAIN, 0.));
		model.InputPropertyValue("velocity", makeVector(PLAIN, PLAIN, 0.0, 0.));
		model.InputPropertyValue("pore velocity", makeVector(PLAIN, PLAIN, 0.0, 0.0));
		model.InputPropertyValue("nodal velocity", makeVector(PLAIN, PLAIN, 0.0, 0.0));
		model.InputPropertyValue("nodal pore velocity", makeVector(PLAIN, PLAIN, 0.0, 0.0));
		model.InputPropertyValue("nodal volume flux", makeScalar(PLAIN, 0.0));

    // creating a 10x10m region in the center of the 2D model, which is initialised with elevated non-wetting phase saturation
		model.FormRectangularRegion("central", Point<2U>(45, 45), Point<2U>(55, 55));
		Region<2U>& central = model.Region("central");
		central.InputPropertyValue("fluid pressure", makeScalar(DIRICH, 2.e5));
		central.InputPropertyValue("saturation oil", makeScalar(DIRICH, 0.9));
    central.InputPropertyValue("saturation water", makeScalar(DIRICH, 0.1));
		model.InputBoundaryValue(LEFT, "fluid pressure", makeScalar(DIRICH, 0.));
		model.InputBoundaryValue(RIGHT, "fluid pressure", makeScalar(DIRICH, 0.));
		model.InputBoundaryValue(LEFT, "saturation oil", makeScalar(DIRICH, 0.));
		model.InputBoundaryValue(LEFT, "saturation oil", makeScalar(DIRICH, 0.));
   // flagging the nimble region nodes and elements as 1, which is a value that gets incremented as the nimble egion evolves
    double64 nimble_val(1.);
    model.InputPropertyValue("nimble nodes", makeScalar(PLAIN, 0.));
    model.InputPropertyValue("nimble elements", makeScalar(PLAIN, 0.));
    central.InputPropertyValue("nimble nodes", makeScalar(PLAIN, nimble_val));
    central.InputPropertyValue("nimble elements", makeScalar(PLAIN, nimble_val));
    nimble_val += 1.;

		// ---------------------------------------------------------------------
		// 2. Use Brooks-Corey model to compute the relative permeabilities
		// ---------------------------------------------------------------------
		BrooksCorey<2U> relperm_model(p_ref, "brooks corey parameter", "entry pressure");
		Region<2U>& model_domain = model.Region("Model");
		ComputeTotalMobility( model_domain, p_ref, relperm_model );

		// output the range of the result variable
		printRangeOfVariable(model, "total mobility");


		// ------------------------------------------------------------------------------------------
		// 3. Setting up FE algorithm to solve pressure diffusion equation c dp/dt = div(K grad p) + S
		//     p = fluid pressure
		//     c = compressibility (fluid and rock)
		//     K = k/mu = hydraulic conductivity (from above)
		//     S = volumetric source term
		//
		//     We solve the discretised equation full implict as
		//
		//     ([c]/dt + [K]){p}t+dt = {c}/dt{p}t + {S}t+dt
		//
		//     Note: [] denotes a matrix, {} a vector
		//
		//     This results in the linear system [A] * {x} = {b}
		//     where [A] is the discretisation of div(K grad p) and c dp/dt
		//     {b} contains the known pressure at time t and the unknown source at
		//     time t+dt; {x} is the unknown pressure at time t+dt that we are solving for
		//
		// ------------------------------------------------------------------------------------------
		// create the CSMP FE Algorithm with SAMG solver to invert linear system
		//GaussJordan_Solver*  solver_test = new GaussJordan_Solver();
		//GaussJordan_Solver*  solver_sub_region = new GaussJordan_Solver();
		SAMG_Settings settings;
		SAMG_Solver   solver_test;
    // settings nrc=11 invokes Intel's PARDISO, set-levelx=1 (see SAMG manual, p. 104)
    /*
        call samg_set_levelx(1)            ! intentional one-level method
        call samg_set_clsolver_finest(111) ! select pardiso
        iswtch=41                          ! cleanup & save decomposition
        call samg(..,ia,ja,a,f1,u,..)      ! use pardiso with RHS f1
        iswtch=11                          ! re-use decomposition
        call samg(..,ia,ja,a,f2,u,..)      ! use pardiso with RHS f2
        call samg(..,ia,ja,a,f3,u,..) .... ! use pardiso with RHS f3
     */
    settings.ExplicitSecondary( true ); // to get access to SAMG's hidden parameters
    settings.Set_levelx( 1 ); // single-level solver
//    settings.Set_clsolver_finest (111 ); NOT IMPLEMENTED YET
		SAMG_Solver   solver_sub_region(&settings);
		PDE_Integrator_UoM<2U, Region>        pde_model(solver_test);
		PDE_Integrator_UoM<2U, NimbleRegion>  pde_plume(solver_sub_region);

		// LHS conductance matrix                            operand         basis function    test function
		NumIntegral_dNT_op_dN_dV<2U, Element<2U> >  stiffness_matrix(p_ref, "total mobility", "fluid pressure", "fluid pressure");
		// LHS mass matrix
		NumIntegral_NT_lhsop_N_dV<2U, Element<2U> > mass_matrix_lhs(p_ref, "total system compressibility", "fluid pressure", "fluid pressure");

		// RHS mass vector
		NumIntegral_NT_op_N_dV<2U, Element<2U> >    mass_matrix_rhs(p_ref, "total system compressibility", "fluid pressure");

		// RHS mass vector for integrating source term
		NumIntegral_NT_op_N_dV<2U, Element<2U> >    source_term(p_ref, "fluid volume source", "fluid pressure");

		// mass matrices for dp/dt term must be divided by time increment
		mass_matrix_lhs.MultiplyWithTimeIncrement(true);
		mass_matrix_rhs.MultiplyWithTimeIncrement(true);

		// use lumped formulation for all mass matrices (i.e., diagonalise matrices)
		mass_matrix_lhs.LumpedFormulation(true);
		mass_matrix_rhs.LumpedFormulation(true);
		source_term.LumpedFormulation(true);

		// evalute source term last
		source_term.AddAccumulateLater();
  
		// define a post-processing step that computes the velocity in each finite element by solving Darcy's law
//		TwoPhaseVelocityAndVolumeFlux<2U>  velo( model, relperm_model ); // "total mobility", "porosity", "fluid pressure", false); // true = extrapolate element velocities to nodes
  
		pde_model.Add(&stiffness_matrix);
		pde_model.Add(&source_term);
		pde_model.Add(&mass_matrix_lhs);
		pde_model.Add(&mass_matrix_rhs);
//		pde_model.AddPostProcess(&velo);

		pde_plume.Add(&stiffness_matrix);
		pde_plume.Add(&source_term);
		pde_plume.Add(&mass_matrix_lhs);
		pde_plume.Add(&mass_matrix_rhs);
//		pde_plume.AddPostProcess(&velo);
  
  
		// ------------------------------
		// 4. Setting up transport scheme
		// ------------------------------
		TwoPhaseExplicitNodeCenteredFVTransport<2U, ExplicitStencilProcessor> transport("Model", model, "porosity",
                                                                                    "diffusivity",
                                                                                    "saturation water",
                                                                                    "saturation oil",
                                                                                    "velocity",
                                                                                    "nodal fluid volume source",
                                                                                    false);
		// define some constant variables
		double64 model_time(0);
		const double64    day(86400.0);
		const double64    max_time(9.0*day);     // run for 10 days
		size_t            frequency(3);
		double64          time_increment(frequency * day);      // timestep 3 day
		const double64 sub_time_increment(time_increment/ frequency);


		// -------------------------------------------------------------------------
		// 5. creating a flexible NimbleRegion that moves with the saturation front
		// -------------------------------------------------------------------------
    Region<2U>&       sub_domain(model.Region("central"));
    cout <<"\nrun: creating NimbleRegion object from region '"<< sub_domain.Name() <<"' with "<< sub_domain.Nodes() <<" nodes and "<< sub_domain.Elements() <<" elements.\n";
		NimbleRegion<2U>  plume_region( sub_domain.NodesBegin(), sub_domain.NodesEnd() ); // will have a halo of one element extra
    plume_region.Out(); // tested: O.K.
    
    // first check
    VTU_Interface<2U>  vtu(model);
    // augmenting the extra perimeter nodes by 2, to make sure they are visible
    const csmp::Index  nn_key(model.Database().StorageKey("nimble nodes"));
    for ( auto nit=plume_region.PerimeterNodesBegin(); nit!=plume_region.NodesEnd(); ++nit ) {
         double64 nval = (*nit)->Read( nn_key ) + 2.;
         (*nit)->Store( nn_key, makeScalar(PLAIN,nval) );
      }

    vtu.OutputDataToVTU("nimble-nodes", "nimble nodes", "Model", 0 );

    
    // The idea is that the NimbleRegion receives a vector of nodes that it consists of and rebuilds
    // itself such that it contains all the elements that contain at least one of these nodes
		std::vector<Node<2U>*>  plume_nodes;
		const double64          increment(10.);

		while (model_time < max_time)
      {
        // compute 2phase properties
        ComputeTotalMobility( model_domain, p_ref, relperm_model);

        // compute pressure
        pde_model.IntegrateOver( model_domain );
        compute_vt_AtBaryCenter( model, "Model" );
        printRangeOfVariable( model, "fluid pressure" );
        printRangeOfVariable( model, "total velocity" );
        vtu.OutputDataToVTU("fluid_pressure", "fluid pressure", "Model", 1000000000);
      
        // for loop to compute saturation and pressure
        for (size_t i(0); i < frequency; ++i)
          {
              model_time += sub_time_increment;
              Point<2U>bottom_left(40. - i * increment, 40. - i * increment);
              Point<2U>top_right(60. + i * increment, 60. + i * increment);

              for ( auto nIter = model_domain.NodesBegin(); nIter != model_domain.NodesEnd(); ++nIter )
                if ((*nIter)->Coordinate().IsBetween(bottom_left, top_right))
                  plume_nodes.push_back(*nIter);

              // NIMBLE REGION IS REBUILD TO FOLLOW SATURATION FRONT
//              plume_region.Initialise( plume_nodes.begin(), plume_nodes.begin() );

              // compute2PhaseMobilityAtBaryCenter( model_domain, relperm_model );
              //ComputeTotalMobility(plume_region, p_ref, relperm_model);
              pde_plume.IntegrateOver(plume_region);
              compute_vt_AtBaryCenter( model, "central" );

              // compute advection of phases
              transport.TransportPhase(relperm_model, sub_time_increment);
              vtu.OutputDataToVTU("volume_flux", "nodal volume flux", "Model", model_time);
              vtu.OutputDataToVTU("fluid_pressure", "fluid pressure", "Model", model_time);
              vtu.OutputDataToVTU("saturation oil", "saturation oil", "Model", model_time);

              // augmenting the perimeter nodes by one, to make sure they are visible
              for ( auto nit=plume_region.PerimeterNodesBegin(); nit!=plume_region.NodesEnd(); ++nit )
                (*nit)->Store( nn_key, makeScalar(PLAIN,(*nit)->Read( nn_key ) + 1. ) );

              vtu.OutputDataToVTU("nimble-nodes", "nimble nodes", "Model", model_time );
              vtu.OutputDataToVTU("nimble-elements", "nimble elements", "Model", model_time );
          }

        // runtime info
        cout << "\n\nrun: RUNTIME (DAYS): " << model_time / day << endl << endl;
     }

} // end run





void NimbleRegion_Test::ConstrainPlumeBoundary( NimbleRegion<2U>& region, const PropertyDatabase<2U>& p_ref)
{
		const Index pressureKey = p_ref.StorageKey("fluid pressure");
		size_t i(0);
		for (auto node = region.PerimeterNodesBegin() ; node != region.NodesEnd(); ++ node) {
			(*node)->Status(pressureKey, DIRICH);
			std::cout << ++i << std::endl;
		}
}




void NimbleRegion_Test::ReleasePlumeBoundary( NimbleRegion<2U>& region, const PropertyDatabase<2U>& p_ref) {
		const Index pressureKey = p_ref.StorageKey("fluid pressure");
		for (auto node = region.PerimeterNodesBegin(); node != region.NodesEnd(); ++node) {
			(*node)->Status(pressureKey, PLAIN);
		}
	}


// AUXILIARY FUNCTIONS

void  NimbleRegion_Test::ComputeTotalMobility( NimbleRegion<2U>& region, const PropertyDatabase<2U>& p_ref, TwoPhaseModel<2U>& relperm )
	{
		// keys to properties
		const Index mobt_key =  p_ref.StorageKey("total mobility");
		const Index satw_key =  p_ref.StorageKey("saturation water");
		const Index satn_key =  p_ref.StorageKey("saturation oil");

		const double64  one(1.);
		double64        sw;
		ScalarVariable  mob_t;

		// 1. Computing the saturation of water = 1 - So
		//    loop over the FE nodes
		vector<Node<2U>* >::const_iterator nit;
		for (nit = region.NodesBegin(); nit != region.NodesEnd(); nit++)
		{
			// read in So, compute Sw and store back to nodes along with the flag of So
			sw = one - (*nit)->Read(satn_key);
			(*nit)->Store(satw_key, makeScalar((*nit)->Status(satn_key), sw));
		}


		// 2. Computing the multiphase flow properties
		//    loop over finite elements
		vector<Element<2U>* >::const_iterator eit;
		for (eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++)
		{
			// 1. setting up the relative permeability model
			// ---------------------------------------------
			relperm.Initialize(*(*eit));
			relperm.InitializeForBaryCenter(*(*eit));
			relperm.EffectiveSaturation();

			// 2. total mobility
			// -----------------
			mob_t = relperm.TotalMobility();
			(*eit)->Store(mobt_key, mob_t);
		}


	} // end Compute2PhaseFlowProperties





void  NimbleRegion_Test::ComputeTotalMobility(Region<2U>& region, const PropertyDatabase<2U>& p_ref, TwoPhaseModel<2U>& relperm)
	{
		// keys to properties
		const Index mobt_key = p_ref.StorageKey("total mobility");
		const Index satw_key = p_ref.StorageKey("saturation water");
		const Index satn_key = p_ref.StorageKey("saturation oil");

		const double64  one(1.);
		ScalarVariable  mob_t;

		// 1. Computing the saturation of water = 1 - So
		//    loop over the FE nodes
		for ( vector<Node<2U>* >::const_iterator nit = region.NodesBegin(); nit != region.NodesEnd(); nit++)
		{
			// read in So, compute Sw and store back to nodes along with the flag of So
			double64 sw = one - (*nit)->Read(satn_key);
			(*nit)->Store(satw_key, makeScalar((*nit)->Status(satn_key), sw));
		}


		// 2. Computing the multiphase flow properties
		//    loop over finite elements
		for ( vector<Element<2U>* >::const_iterator eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++)
		{
			// 1. setting up the relative permeability model
			// ---------------------------------------------
			relperm.Initialize(*(*eit));
			relperm.InitializeForBaryCenter(*(*eit));
			relperm.EffectiveSaturation();

			// 2. total mobility
			// -----------------
			mob_t = relperm.TotalMobility();
			(*eit)->Store(mobt_key, mob_t);
		}

	} // end compute2PhaseFlowProperties






void  compute_vt_AtBaryCenter( Model<2U>& model, const char* flow_domain )
 {
    ErrorHandler&       csmp_error( ErrorHandler::Instance() );
    ScalarVariable      flux;
    VectorVariable<2U>  velo, gvt; // total velocity and gravity term
    double64            vmin, vmax;
    DenseMatrix<DM_MIN> DERIV;
   
    // variable keys
    const csmp::Index lt_key(model.Database().StorageKey("total mobility"));
    const csmp::Index pf1_key(model.Database().StorageKey("fluid pressure"));
//    const csmp::Index gt_key(model.Database().StorageKey("nimble nodes"));
    const csmp::Index vt_key(model.Database().StorageKey("total velocity"));
    const csmp::Index vf1_key(model.Database().StorageKey("volume flux"));
   
    // establishing permissible range of output variable
    model.Database().RangeOf( "total velocity", vmin, vmax );
   
    Region<2U>&  gref(model.Region(flow_domain));
 
    for ( vector<Element<2U>*>::iterator
          it=gref.ElementsBegin(); it!=gref.ElementsEnd(); it++ )
      {
         // computing the total velocity: vt = -k t (lt grad p - (lw rho_w + lo rho_o) g)
         // (taking into account the thickness attribute of 2D elements)
         const double64 total_mobility = (*it)->Read( lt_key );
         velo = 0.;
         (*it)->dN_AtBaryCenter( DERIV, 1U );
         for ( size_t i=0U; i<(*it)->Nodes(); i++ ) {
              double64 pf = (*it)->N(i)->Read( pf1_key );
              velo(0)  += pf  * -DERIV(0,i) * total_mobility;
              velo(1)  += pf  * -DERIV(1,i) * total_mobility;
              //velo(2)  += pf  * -DERIV(2,i) * total_mobility;
           }
         // taking into account the gravitational flow component for calculation of velocity of non-wetting phase
         // (reading following vector variable from element)
         //       gvelo = t * k * g * (lambda_w * rho_w + lambda_o * rho_o) grad z
         //       in case of an inclined element it is sin alpha (=vertical / vector length) * gvelo
 //        (*it)->Read( gt_key, gvt );
        
         // adding viscous + gravitional velocity components and applying thickness multiplier
 //        velo += gvt;
        
         if ( !velo.IsWithinRange(vmin,vmax) ) {
               (*it)->Out();
               velo.Out();
               csmp_error.notice( WARNING, "compute_vt_AtBaryCenterIncludingGravity",
                                 "computed 'total velocity' is out of range specified in variables file.");
           }
        
         // storing the computed velocity (checking the variable status flags)
         (*it)->Store( vt_key, velo );
        
        // volume flux
        flux = velo.Length();
        (*it)->Store( vf1_key, flux );
      }

    printRangeOfVariable( model, "total velocity" );
   
} // end Compute_vt_AtBaryCenterIncludingGravity





} // csmp
