#include "Tutorial4_Example_Revisited.h"

// CSMP model
#include "ANSYS_Model2D.h"
// FE algorithm
#include "PDE_Integrator_CRM.h"
// PDE operators building the FE algorithm
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_dNi_dV_sc.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "PointSource_rhsop.h"
#include "LinearSolver.h"
#include "NumIntegral_NT_op_N_dS.h"
// interfaces
#include "VTU_Interface.h"
// utility functions
#include "CSMP_highLevelUtilities.h"
// visitors
#include "StabilizationParameterVisitor.h"

using namespace std;

namespace csmp {
	// This constructor gets address of input files and model
	Tutorial4_Example_Revisited::Tutorial4_Example_Revisited(string inpf, string outf){
		_inpdata.set_ = true;
		_inpdata._model_name = inpf;
		_inpdata._output_name = outf;
	}

	void Tutorial4_Example_Revisited::Specifications()
	{
		SetTitle("Tutorial 4 Revisited: Stokes lubrication equation");
		SetDifficulty(3);
		SetCategory("Simulation of Physical Processes by CRM integrator");
		AddAuthor("Sebastian Geiger. Edited by Hani Akbari");
		AddDescription("A CSMP main file that uses an solves the coupled Stokes equation in a 2D geometry representing");
		AddDescription("pores and grains in a carbonate rock. The solution of the Stokes equation provides the pressure");
		AddDescription("and velocity fields inside the connected pore scale. Note that special SAMG settings are needed");
		AddDescription("for solving the Stokes equation using a standard FE method. Furthermore, a stabilisation parameter");
		AddDescription("is introduced to allow the discretisation of velocity and pressure in the Stokes equation using");
		AddDescription("the same FE basis functions. CAUTION: We set source term on left, no slip in y direction on top and bottom");
		AddDescription("and no slip around grains.");
		AddRequirement("source code in: 'Tutorial4_Example.cpp'");
		AddRequirement("pores Binary files, stokes_variables.txt, no need of configuration");
	} // Initialize()

	void Tutorial4_Example_Revisited::Run()
	{
		// ----------------------------------------------------------------------------------
		// 0.0 Set variables used throughout the simulation
		// ----------------------------------------------------------------------------------
		clock_t start(clock()); // record the CPU time
		string input_file, outputName("");
		if (!_inpdata.set_) {
			cout << "\nTutorial4_Example: Please enter the name of input mesh. Use double backslash ( default: pores ): ";
			cin >> input_file;
		}
		else {
			input_file = _inpdata._model_name;
			outputName = _inpdata._output_name;
		}
		// ------------------------------
		// 1.0 Build the CSMP Model
		// ------------------------------
		ANSYS_Model2D model((input_file + "pores").c_str(), (input_file + "stokes_variables.txt").c_str());
		const PropertyDatabase<2U>&  p_ref = model.Database();
		Region<2U>&                  r_ref = model.Region("PORES"); // reference to PORE region where calculations are performed
		scaleRegion(model, 20000.0); // scale the model size by 1/value
		printModelDimensions(model, true);
		// ------------------------------------
		// 2.0 Boundary and initial conditions
		// ------------------------------------

		// input material and fluid properties
		const double64 viscosity(0.001); // in Pa sec
		model.InputPropertyValue("viscosity", makeScalar(PLAIN, viscosity));
		model.InputPropertyValue("porosity", makeScalar(PLAIN, 1.0));
		model.InputPropertyValue("zero", makeScalar(PLAIN, 0.0)); // dummy variable needed to close coupled FE algorithm

	    // assigning initial conditions
		model.InputPropertyValue("fluid pressure", makeScalar(PLAIN, 0.0));
		model.InputPropertyValue("nodal velocity x", makeScalar(PLAIN, 0.0));
		model.InputPropertyValue("nodal velocity y", makeScalar(PLAIN, 0.0));
		model.InputPropertyValue("Neumann traction x", makeScalar(PLAIN, 0.0));
		model.InputPropertyValue("Neumann traction y", makeScalar(PLAIN, 0.0));

		// Dirichlet (no-slip) conditions for velocity at boundaries of mineral grains
		model.Region("GRAINS").ChangePropertyStatus("nodal velocity x", DIRICH);
		model.Region("GRAINS").ChangePropertyStatus("nodal velocity y", DIRICH);

		// Neumann condition for traction term. This uses a work around to translate the flux
		// into a nodal point source because the FE algorithm integration of surface fluxes
		// does not work yet and should be replaced once the PDE_Integrator can be used
		// to compute flux boundary conditions
		// model.InputBoundaryValue(LEFT, "Neumann traction x", makeScalar(DIRICH,  200.0));
		// translate into nodal source term. We only do that for LEFT X. Nothing for others
		// assignFluxToPointSource(model, "Neumann traction x"); 

		// visitor for calculation of the stabilization parameter, which is needed to compute
		// Stokes flow using FEs which have the same basis functions for velocity and pressure
		StabilizationParameterVisitor<2U> stabpar_visitor(model, viscosity, 2.0, "stabilization parameter");
		model.Accept(stabpar_visitor);
		printRangeOfVariable(model, "stabilization parameter");

		// ------------------------------------------------------------
		// 3.0 Coupled FE algorithm for solving Stokes equation
		// ------------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
		SAMG_Settings settings;
		settings.Set_nxtyp(2);
		settings.Set_ncycle(100);		
		//settings.Set_napproach(2);
		// FE algorithm with specialised SAMG settings
		SAMG_Solver  samg_solver(&settings);
		PDE_Integrator_CRM<2U, Region>  stokes_flow(samg_solver);
#else
		CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
		PDE_Integrator<2U, Region>  stokes_flow(linear_solver);
#endif
		// Stokes lubrication equation

		// LHS operators                                        operand       basic function      test function
		NumIntegral_dNT_op_dN_dV<2U>  viscosity_matr_x(p_ref, "viscosity", "nodal velocity x", "nodal velocity x");
		NumIntegral_NT_dNi_dV_sc<2U>  gradient_x(p_ref, "zero", "fluid pressure", "nodal velocity x");
		gradient_x.SpatialDerivative(X_DIRECTION);
		//    gradient_x.Transposed();
		gradient_x.MultiplyBy(-1.0);

		NumIntegral_dNT_op_dN_dV<2U>  viscosity_matr_y(p_ref, "viscosity", "nodal velocity y", "nodal velocity y");
		NumIntegral_NT_dNi_dV_sc<2U>  gradient_y(p_ref, "zero", "fluid pressure", "nodal velocity y");
		gradient_y.SpatialDerivative(Y_DIRECTION);
		//    gradient_y.Transposed();
		gradient_y.MultiplyBy(-1.0);

		// LHS operators                                    operand  basic function      test function
		NumIntegral_NT_dNi_dV_sc<2U>  divergence_x(p_ref, "zero", "nodal velocity x", "fluid pressure");
		divergence_x.SpatialDerivative(X_DIRECTION);
		divergence_x.Transposed();
		//    divergence_x.MultiplyBy(-1.0);

		NumIntegral_NT_dNi_dV_sc<2U>  divergence_y(p_ref, "zero", "nodal velocity y", "fluid pressure");
		divergence_y.SpatialDerivative(Y_DIRECTION);
		divergence_y.Transposed();
		//    divergence_y.MultiplyBy(-1.0);

		NumIntegral_dNT_op_dN_dV<2U>  stab_matrix(p_ref, "stabilization parameter", "fluid pressure", "fluid pressure");

		// RHS operators                              operand               test function
		
		

		// RHS operators                           operand   test function
		NumIntegral_NT_op_N_dV<2U>  dummy(p_ref, "zero", "fluid pressure");

		// add each PDE Operator to the FE algorithm

		// Navier-Stokes equation
		stokes_flow.Add(&viscosity_matr_x);
		stokes_flow.Add(&gradient_x);
		stokes_flow.Add(&viscosity_matr_y);
		stokes_flow.Add(&gradient_y);
		
		

		// Continuity equation
		stokes_flow.Add(&divergence_x);
		stokes_flow.Add(&divergence_y);
		stokes_flow.Add(&stab_matrix);
		stokes_flow.Add(&dummy);

		// Face Source
		/*********************/
		/*
		model.InputBoundaryValue(LEFT, "nodal velocity x", makeScalar(DIRICH, 0.0));
		model.InputBoundaryValue(RIGHT, "nodal velocity x", makeScalar(DIRICH, 0.0));
		//model.Boundary("TOP").InputPropertyValue("face source", makeScalar(NEUMANN, -1.0e+1));
		model.Boundary("BOTTOM").InputPropertyValue("face source", makeScalar(NEUMANN, 1.0e+1));
		NumIntegral_NT_op_N_dS<2U, Face>  face_source(p_ref, "face source", "nodal velocity y");
		PointSource_rhsop<2U>  extpressure_x(p_ref, "Neumann traction x", "nodal velocity x");
		stokes_flow.Add(&extpressure_x);
		face_source.LumpedFormulation();
		stokes_flow.AddBoundaryIntegrals(&face_source);
		*/

		// Dirichlet (no-slip) conditions for velocity at Top and Bottom
		
		model.InputBoundaryValue(TOP, "nodal velocity y", makeScalar(DIRICH, 0.0));
		model.InputBoundaryValue(BOTTOM, "nodal velocity y", makeScalar(DIRICH, 0.0));
		//model.Boundary("LEFT").InputPropertyValue("face source", makeScalar(NEUMANN, 1.0e+1));
		model.Boundary("RIGHT").InputPropertyValue("face source", makeScalar(NEUMANN, -1.0e+1));
		NumIntegral_NT_op_N_dS<2U, Face>  face_source(p_ref, "face source", "nodal velocity x");
		PointSource_rhsop<2U>  extpressure_y(p_ref, "Neumann traction y", "nodal velocity y");
		stokes_flow.Add(&extpressure_y);
		face_source.LumpedFormulation();
		stokes_flow.AddBoundaryIntegrals(&face_source);
		

		

		stokes_flow.IntegrateOver(model,r_ref);
		// solve the Stokes equation and release memory
		//stokes_flow.IntegrateOver(r_ref);
		stokes_flow.Reset();

		// output resulting variable ranges
		printRangeOfVariable(model, "nodal velocity x");
		printRangeOfVariable(model, "nodal velocity y");
		printRangeOfVariable(model, "fluid pressure");

		// construct nodal velocity vector
		constructVelocityVector(model);

		// record timing of the simulation
		clock_t end(clock());
		cerr << "\n\nmain: CPU time was " << static_cast<double64>((end - start) / CLOCKS_PER_SEC) << " seconds " << endl << endl;
		// output results to VTU
		cout << "Writing VTU ..." << endl;
		VTU_Interface<2U>   vtu(model);
		std::list<std::string> outputProps;
		outputProps.push_back("fluid pressure");
		outputProps.push_back("nodal velocity");

		// output pressure and (nodal) velocity to the same VTU file
		vtu.OutputDataToVTU(outputName + "PoreScaleVariables", outputProps, "PORES", static_cast<int>(0));
		// terminate
		cerr << "\n\nmain: That's it, run completed successfully..." << endl;
	} // Run()


	/**
	 No change here by Hani
	*/

	/**
	Quick hack that loops over a 2D boundary, reads in a nodal source term and scales it by 0.5 of the length
	of the edges of the connected FEs wich lie at the model boundary (here it is the LEFT boundary)
	*/
	void Tutorial4_Example_Revisited::assignFluxToPointSource(Model<2U>& mdl, const char* flux)
	{

		// backup the fluxes that are orginally assigned in a temporary variable
		std::string         temp_flux("flux");

		if (!mdl.Database().IsDefined(temp_flux.c_str()))
			mdl.CreateProperty(temp_flux.c_str(), "X", SCALAR, NODE);
		mdl.CopyReplace(flux, temp_flux.c_str());

		// define variables needed to compute the length of the FE edges that lie at the boundary
		// and read the temporary flux and store the final flux
		csmp::Index       tf_key(mdl.Database().StorageKey(temp_flux.c_str())),
			f_key(mdl.Database().StorageKey(flux));
		ScalarVariable    tf, f;
		double64          y[2], area, length(0.0);
		size_t            j;

		cout << "\nassignFluxToPointSource: Translating '" << flux << "' into a nodal point source" << endl;

		// final flux is zero initionally
		mdl.InputPropertyValue(flux, makeScalar(PLAIN, 0.0)); // set to zero initially

															  // loop over all finite elements and identify elements that lie at the same model boundary of interest (here LEFT)
		const Region<2U>&   mref = mdl.Region("Model");
		vector<Element<2U>* >::const_iterator eit;
		for (eit = mref.ElementsBegin(); eit != mref.ElementsEnd(); eit++) {
			if ((*eit)->AtBoundary() == LEFT or
				(*eit)->AtBoundary() == CNR1 or
				(*eit)->AtBoundary() == CNR4) {
				j = 0;
				// first loop to calculate length of the FE edge that lies at the boundary
				for (size_t i = 0; i<(*eit)->Nodes(); i++) {
					if ((*eit)->N(i)->AtBoundary() == LEFT or
						(*eit)->N(i)->AtBoundary() == CNR1 or
						(*eit)->N(i)->AtBoundary() == CNR4) {
						y[j] = (*eit)->N(i)->y();
						j++;
					}
				}
				// provide a warning if the edge has less than 2 nodes (as it is a 2D models, 2D elements should have 2 nodes
				// at the model boundary)
				if (j != 2) {
					cerr << "\nassignFluxToPointSource: ERROR: Counted less than two boundary nodes for element " << endl;
					(*eit)->Out();
					area = 0.0;
				}
				// calculate length (area)
				else {
					area = y[0] - y[1];
					if (area < 0.0) area *= -1.0;
					length += area;
					area /= static_cast<double64>(j); // 2 nodes per triangle or quadrilateral
				}
				// second loop to calculate and scale nodal flux
				for (size_t i = 0; i<(*eit)->Nodes(); i++) {
					if ((*eit)->N(i)->AtBoundary() == LEFT or
						(*eit)->N(i)->AtBoundary() == CNR1 or
						(*eit)->N(i)->AtBoundary() == CNR4) {
						// read existing flux at node i
						tf = (*eit)->N(i)->Read(tf_key);
						// read existing source at node i
						f = (*eit)->N(i)->Read(f_key);
						// add current heat flux to this value
						f += tf() * area;
						// store it back to node
						(*eit)->N(i)->Store(f_key, f);
					}
				}
			}
		}

		printRangeOfVariable(mdl, flux);
		cout << "\nassignFluxToPointSource: Successfully assigned '" << flux << "' to the nodes... " << endl;

	}


	/// simple combine x and y component of computed velocities into a vector and saves it to the grod
	void Tutorial4_Example_Revisited::constructVelocityVector(Model<2U>& mdl)
	{
		VectorVariable<2U> v;
		ScalarVariable     p;
		csmp::Index vx_key(mdl.Database().StorageKey("nodal velocity x")),
			vy_key(mdl.Database().StorageKey("nodal velocity y")),
			v_key(mdl.Database().StorageKey("nodal velocity"));

		const Region<2U>& mref = mdl.Region("Model");
		vector<Node<2U>* >::const_iterator nit;
		for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++)
		{
			v(0) = (*nit)->Read(vx_key);
			v(1) = (*nit)->Read(vy_key);
			(*nit)->Store(v_key, v);
		}
	}


	/// scale the size of the CSMP model (NB divides by the provided factor!)
	void Tutorial4_Example_Revisited::scaleRegion(Model<2U>& mdl, double64 scale_factor)
	{
		double64       x_, y_;
		const double64 factor(scale_factor);

		static const Region<2U>& mref = mdl.Region("Model");
		vector<Node<2U>* >::const_iterator nit;
		for (nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++)
		{
			x_ = (*nit)->x();
			(*nit)->x(x_ / factor);
			y_ = (*nit)->y();
			(*nit)->y(y_ / factor);
		}
	}

} // csmp
