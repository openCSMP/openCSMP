#include "NimbleRegion_Test.h"

// the CSMP model
#include "ANSYS_Model2D.h"
#include "ModelTopology.h"

// FE algorithm
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"

// FV algorithms
#include "TwoPhaseImplicitNodeCenteredFVTransport.h"
#include "TwoPhaseExplicitNodeCenteredFVTransport.h"
#include "ExplicitStencilProcessor.h"

// relative permeability calculations
#include "BrooksCorey.h"

// monitoring individual regions
#include "RegionMonitor.h"

// interfaces
#include "InputDataManager.h"
#include "MatlabInterface.h"
#include "VTU_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"
#include "GaussJordan_Solver.h"
#include "PDE_Integrator_UoM_Mock.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
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
		model.InputBoundaryValue(LEFT, "fluid pressure", makeScalar(DIRICH, 0.));
		model.InputBoundaryValue(RIGHT, "fluid pressure", makeScalar(DIRICH, 0.));
		model.InputBoundaryValue(LEFT, "saturation oil", makeScalar(DIRICH, 0.));
		model.InputBoundaryValue(LEFT, "saturation oil", makeScalar(DIRICH, 0.));

  

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

		// LHS stiffness matrix                            operand         basis function    test function
		NumIntegral_dNT_op_dN_dV<2U, Element<2U> >  stiffness_matrix(p_ref, "total mobility", "fluid pressure", "fluid pressure");
		//NumIntegral_dNT_dN_dV<2U, Element<2U> >  stiffness_matrix(p_ref, "fluid pressure", "fluid pressure");
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
		VelocityAndVolumeFlux<2U>  velo(model, "total mobility", "porosity", "fluid pressure", false); // true = extrapolate element velocities to nodes
  
		pde_model.Add(&stiffness_matrix);
		pde_model.Add(&source_term);
		pde_model.Add(&mass_matrix_lhs);
		pde_model.Add(&mass_matrix_rhs);
		pde_model.AddPostProcess(&velo);

		pde_plume.Add(&stiffness_matrix);
		pde_plume.Add(&source_term);
		pde_plume.Add(&mass_matrix_lhs);
		pde_plume.Add(&mass_matrix_rhs);
		pde_plume.AddPostProcess(&velo);
  
  
		// ------------------------------
		// 4. Setting up transport scheme
		// ------------------------------
		TwoPhaseExplicitNodeCenteredFVTransport<2U, ExplicitStencilProcessor>* transport = new TwoPhaseExplicitNodeCenteredFVTransport<2U, ExplicitStencilProcessor>("Model", model,
			"porosity",
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
		NimbleRegion<2U>  plume_region( sub_domain.NodesBegin(), sub_domain.PerimeterNodesBegin() );
    
		std::vector<Element<2U>*> plume_elements;
		std::vector<Node<2U>*> plume_nodes;
		VTU_Interface<2U> vtu(model);
		double64 inc(10.);

		while (model_time < max_time)
      {
        // compute 2phase properties
        ComputeTotalMobility( model_domain, p_ref, relperm_model);

        // compute pressure
        pde_model.IntegrateOver( model_domain );
        vtu.OutputDataToVTU("fluid_pressure", "fluid pressure", "Model", 1000000000);
      
        // for loop to compute saturation and pressure
        for (size_t i(0); i < frequency; ++i)
          {
              model_time += sub_time_increment;
              Point<2U>bottom_left(40. - i*inc, 40. - i*inc);
              Point<2U>top_right(60. + i*inc, 60. + i*inc);

              for ( auto nIter = model_domain.NodesBegin(); nIter != model_domain.NodesEnd(); ++nIter ) {
                if ((*nIter)->Coordinate().IsBetween(bottom_left, top_right)) {
                  plume_nodes.push_back(*nIter);
                }
          }

          plume_region.Rebuild( plume_nodes.begin(), plume_nodes.begin() );

          ComputeTotalMobility(plume_region, p_ref, relperm_model);
          
          // UPDATING THE SHAPE OF THE NIMBLE REGION
          ConstrainPlumeBoundary(plume_region, p_ref);
          pde_plume.IntegrateOver(plume_region);
          ReleasePlumeBoundary(plume_region, p_ref);
          plume_region.Erase();
          
          // compute advection of phases
          transport->TransportPhase(relperm_model, sub_time_increment);
          vtu.OutputDataToVTU("volume_flux", "nodal volume flux", "Model", model_time);
          vtu.OutputDataToVTU("fluid_pressure", "fluid pressure", "Model", model_time);
          vtu.OutputDataToVTU("saturation oil", "saturation oil", "Model", model_time);

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


} // csmp
