#include "Hani_Example.h"
#include "VTK_Interface.h"
#include "string.h"
#include "ConstantFactor.h"
// relative permeability calculations
#include "BrooksCorey.h"
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"
#include "VTU_Interface.h"
#include "InputDataManager.h"
#include "Integral_NT_op_N_dV.h"
#include "Integral_dNT_op_dN_dV.h"
#include "Integral_dNT_op_dN_dV_Analytical.h"
#include "PDE_Integrator_CRM.h"


using namespace std;

namespace csmp {
	void Hani_Example::Specifications()
	{
		SetTitle("Hani_Example");
		SetDifficulty(0);
		SetCategory("0000 (Hani) Test Category");
		AddAuthor("H. Akbari");
		AddDescription("Still needs NOTHING");
	}

	void Hani_Example::Run()
	{
		// return;
		cout << " This is Hani Run \n";
		const string  model_path("D:\\myCSMP_Data\\hex2_10_single_frac_model\\");
		bool isoparametric(false);
		bool CRM(true);
		bool LumpedRHS(true);

		#ifdef CSMP_WITH_SAMG_SOLVER
			SAMG_Solver  linear_solver;
		#else
			cout << " ***NOT*** in SAMG\n ";
			CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
		#endif

		ANSYS_Model3D model((model_path + "hex2_10").c_str(), 
							(model_path + "hex2_10-variables.txt").c_str(),
                            false,
                            true,
                            true,
                            true,
                            true,
                            isoparametric);
		InputDataManager<3U>  model_configuration;
		model_configuration.Configure_ANSYS_ModelFromFile(model, (model_path + "hex2_10").c_str());		
		
		const double64 pressure_left(5.0e6); // create pressure drop of 50 bars/500 m
		const double64 pressure_right(1.0e6); // create pressure drop of 50 bars/500 m
		model.InputBoundaryValue(LEFT, "fluid pressure", makeScalar(DIRICH, pressure_left));
		model.InputBoundaryValue(RIGHT, "fluid pressure", makeScalar(DIRICH, pressure_right));

		const PropertyDatabase<3>&   p_ref = model.Database();		
		const double64 fluid_viscosity(1.0);
		ConstantFactor<3U, divides>  conductivity(p_ref, "conductivity", "permeability", fluid_viscosity);
		model.Apply(conductivity);
		// Use the Brooks-Corey model to compute the relative permeabilities
		/*
		BrooksCorey<3U> relperm_model(p_ref, "brooks corey parameter", "entry pressure");
		computeTotalMobility(model, relperm_model);
		// output the range of the result variable
		printRangeOfVariable(model, "total mobility");
		*/
		// binary VTK output, creates much smaller files than VTK
		VTU_Interface<3U>  vtu(model);


		if (isoparametric) {
			SteadyStateDiffusor<3U, Region> fluid_pressure(model,
						"conductivity", "fluid pressure", "fluid volume source", LumpedRHS);
			// calculate initial pressure distribution
			fluid_pressure.ComputeSteadyState(model);
			cout << "\n\n\n Steady state calculated.\n";
			vtu.OutputDataToVTU((model_path + "fluid_pressure_ISO").c_str(), 
				                                 "fluid pressure", "Model", 0);
		}
		else {
			// LHS stiffness matrix       operand         basis function    test function
			Integral_dNT_op_dN_dV_Analytical<3U, Element<3U>>
				stiffness_matrix(p_ref, "conductivity", "fluid pressure", "fluid pressure");
			// Integral_dNT_op_dN_dV<3U, Element<3U>>  
			//stiffness_matrix(p_ref, "conductivity", "fluid pressure", "fluid pressure");
			// RHS mass vector for integrating source term
			Integral_NT_op_N_dV<3U, Element<3U>>
				source_term(p_ref, "fluid volume source", "fluid pressure");
			source_term.LumpedFormulation(LumpedRHS);
			if (CRM) {
				PDE_Integrator_CRM<3U, Region>  fluid_pressure(linear_solver);
				dynamic_cast<SAMG_Settings*>(linear_solver.GetSolverSettings())->Set_isym(1);
				fluid_pressure.Add(&stiffness_matrix);
				fluid_pressure.Add(&source_term);
				// pressure diffusion is computed as the Model applies the FE algorithm
				model.Apply(fluid_pressure);
				vtu.OutputDataToVTU((model_path + "_fluid_pressure_CRM_NON").c_str(),
					"fluid pressure", "Model", 0);
			}
			else {
				PDE_Integrator<3U, Region>  fluid_pressure(linear_solver);
				fluid_pressure.Add(&stiffness_matrix);
				fluid_pressure.Add(&source_term);
				// pressure diffusion is computed as the Model applies the FE algorithm
				model.Apply(fluid_pressure);
				vtu.OutputDataToVTU((model_path + "_fluid_pressure_SM_NON").c_str(),
					"fluid pressure", "Model", 0);
			}
		}

		// operation to compute velocity
		// VelocityAndVolumeFlux<3U, Element<3U>>  
		   //velo(model,"total mobility", "porosity", "fluid pressure", true);

		// add velocity calculation as post process
		// fluid_pressure.AddPostProcess(&velo);
		


		cerr << " Test is finished here currently. Hit return!\n ";
		getchar();
		exit(EXIT_SUCCESS);
		/*
		std::string Path("E:\\CSMP\\csmp-api-library\\examples\\data\\TwoPhase_H\\");
		cout << " Generate Model and Variables \n";
		SYS_Model2D model((Path + "box2d_fault").c_str(), 
		    (Path + "Hani_Test_variables.txt").c_str());
		K_Interface<2U>  vtk;
		k.OutputDataToVTK(model, (Path+"Pressure").c_str(), "fluid pressure", 0);
		*/
		/*
		cout << "\n TEST LinearHexahedron \n";
		// Set 3 coordinates of 8 Nodes 
		double pXY[8][3] = { -1,-1,1, 1,-1,1,  1,-1,-1, -1,-1,-1, -1,1,1, 1,1,1, 1,1,-1, -1,1,-1 };
		// An object 
		DenseMatrix<DM12> XY;
		XY.Resize(8, 3);
		for (auto i = 0; i < 8; ++i)
			for (auto j = 0; j < 3; ++j)
				XY(i, j) = pXY[i][j];
		LinearCuboid_Test E(true, XY);
		E.run();
		*/
		// Integrals of shape functions and derivatives
		//E.TestElementIntegrals(XY);
	}

	void  Hani_Example::computeTotalMobility(ANSYS_Model3D& mdl, TwoPhaseModel<3U>& relperm)
	{
		// keys to properties
		static Index  mobt_key(mdl.Database().StorageKey("total mobility"));
		static Index  satw_key(mdl.Database().StorageKey("saturation water"));
		static Index  sato_key(mdl.Database().StorageKey("saturation oil"));

		const double64  one(1.);
		double64        sw;
		ScalarVariable  mob_t;

		// 1. Computing the saturation of water = 1 - So
		//    loop over the FE nodes
		static const Region<3U>& mref = mdl.Region("Model");
		vector<Node<3U>* >::const_iterator nit;
		for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++)
		{
			// read in So, compute Sw and store back to nodes along with the flag of So
			sw = one - (*nit)->Read(sato_key);
			(*nit)->Store(satw_key, makeScalar((*nit)->Status(sato_key), sw));
		}
		// 2. Computing the multiphase flow properties
		//    loop over finite elements
		vector<Element<3U>* >::const_iterator eit;
		for (eit = mref.ElementsBegin(); eit != mref.ElementsEnd(); eit++)
		{
			// 1. setting up the relative permeability model
			// ---------------------------------------------
			relperm.Initialize(*(*eit));
			relperm.InitializeForBaryCenter(*(*eit));
			relperm.EffectiveSaturation();

			// 2. total mobility
			// -----------------
			mob_t = relperm.TotalMobility();
			if (mob_t() > 1.0e-8 || mob_t() < 0 || isnan(mob_t())) { 
				cout << (*eit)->FE_Type() << "\t" << mob_t() << "\n"; 
			}
			(*eit)->Store(mobt_key, mob_t);
		}


	} // end compute2PhaseFlowProperties
	
}
