#include "Transient_Test.h"
#include "PDE_Integrator_UoM_Mock.h"

#include "GaussJordan_Solver.h"

// the CSMP model
#include "Model.h"
#include "ModelTime.h"

// a simple FE mesh generator
#include "Triangulator.h"

// the FE algorithm
#include "PDE_Integrator.h"

// PDE operators building the FE algorithm
#include "Integral_NT_op_N_dV.h"
#include "Integral_NT_lhsop_N_dV.h"
#include "Integral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

// output interfaces
#include "VTK_Interface.h"
#include "MatlabInterface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"
#include "ConstantFactor.h"



#include "PDE_Integrator_UoM.h"

// PDE operators building the FE algorithm
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_dNT_dN_dV.h"

// FE grid generation
#include "Quadrilaterator.h"

using namespace std;

namespace csmp {

	void Transient_Test::run()
	{
		double& model_time(ModelTime::Instance().modelTime);
		model_time = 0.;
		Quadrilaterator    quadrilaterator; // simple FE mesher
		VSet<2U>           mesh_container;  // container to store the input mesh
		string             file_name("tutorial1_input");
		double           x(10); double y(10);
		cout << "\nmain: Enter the pixel-based input geometry for the quadrilaterator: " << endl;
		cout << "\nmain: The x- and y-dimensions of your model (in m): " << endl;

		// read in file and generate mesh
		quadrilaterator.QuadrilateralsFromRegularGrid(mesh_container, file_name.c_str(), x, y);


		//VSet<2U>   vset=readTextPixelData();
		Model<2U>  model(mesh_container, "tutorial1_variables.txt", true);
		const PropertyDatabase<2>& p_ref(model.Database());  // constand reference to the property database
		Region<2U>& region = model.Region("Model");
		// give the model dimensions
		printModelDimensions(model, true);

		// -----------------------------------------------------------------------
		// 2.0 Now we apply boundary and initial conditions (this can also be done,
		//     more conveniently, in a configuration file for more realistic runs
		// -----------------------------------------------------------------------

		// assigning material properties
		model.InputPropertyValue("porosity", makeScalar(PLAIN, 0.1));     // always as a fraction
		model.InputPropertyValue("permeability", makeScalar(PLAIN, 1.0e-15)); // always in m2 (comment out if heterogeneous k-field is used in input file)
		model.InputPropertyValue("compressibility", makeScalar(PLAIN, 5.0e-10)); // for fluid and rock, in Pa-1

		// assigning initial conditions
		model.InputPropertyValue("fluid pressure", makeScalar(PLAIN, 1.0e+07));  // always in Pascal
		model.InputPropertyValue("fluid volume source", makeScalar(PLAIN, 0.0));      // no sources/sinks (units m3 m-2 s-1)

		// assigning boundary conditions for fluid pressure at the LEFT and RIGHT model boundaries
		// such that a pressure wave travels from left to right through the model
		model.InputBoundaryValue(LEFT, "fluid pressure", makeScalar(DIRICH, 3.0e+07));
		model.InputBoundaryValue(RIGHT, "fluid pressure", makeScalar(DIRICH, 1.0e+07));


		// ----------------------------------------------------------------------------------------------
		// 3.0 Now we use an Interrelation (ConstantFactor, inherited from base class Interrelation
		//     to compute the hydraulic conductivity K = k/mu (k = permeability, mu = viscosity) at each
		//     finite element
		// ----------------------------------------------------------------------------------------------
		//ConstantFactor<2U, divides>  conductivity(p_ref, "conductivity", "permeability", 0.001); // viscosity 1 cp = 0.001 Pa s

		// the Model applies the object "conductivity", which is instantiated from class ConstantFactor
		// the result variable "conductivity" is computed automatically and its range is checked
		//model.Apply(conductivity);
		const Index conductKey = p_ref.StorageKey("conductivity");
		for ( auto eIter = region.ElementsBegin(); eIter != region.ElementsEnd(); eIter++ ) {
			(*eIter)->Store(conductKey, makeScalar(DIRICH, 1.));
		}
		// output the range of the result variable
		printRangeOfVariable(model, "conductivity");

		// ------------------------------------------------------------------------------------------
		// 4.0 Setting up an FE algorithm to solve the diffusion equation c dp/dt = div(K grad p) + S
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

    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;

		//PDE_Integrator<2U,Region>  fluid_pressure(samg_solver);
		PDE_Integrator_UoM_Mock<2U, Region>  pde_validate(linear_solver);
		PDE_Integrator_UoM_Mock<2U, Region>  pde_test(linear_solver);

		NumIntegral_dNT_dN_dV<2U, Element<2U> >  stiffness_matrix(p_ref, "fluid pressure", "fluid pressure");
		// LHS mass matrix
		NumIntegral_NT_lhsop_N_dV<2U, Element<2U> > mass_matrix_lhs(p_ref, "compressibility", "fluid pressure", "fluid pressure");

		// RHS mass vector
		NumIntegral_NT_op_N_dV<2U, Element<2U> >    mass_matrix_rhs(p_ref, "compressibility", "fluid pressure");

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
		VelocityAndVolumeFlux<2U, Element<2U> >     velo(model, "conductivity", "porosity", "fluid pressure", true); // true = extrapolate element velocities to nodes

		// now add each FE operation (i.e., PDE Operator) to the FE algorithm
		pde_validate.Add(&stiffness_matrix);
		pde_validate.Add(&source_term);
		pde_validate.Add(&mass_matrix_lhs);
		pde_validate.Add(&mass_matrix_rhs);
		pde_validate.AddPostProcess(&velo);


		pde_test.Add(&stiffness_matrix);
		pde_test.Add(&source_term);
		pde_test.Add(&mass_matrix_lhs);
		pde_test.Add(&mass_matrix_rhs);
		pde_test.AddPostProcess(&velo);

		// -----------------------
		// 5.0 Time Loop Variables
		// -----------------------
		// write output in VTK format and for Matlab
		VTK_Interface<2U>  vtk_output;
		vtk_output.OutputDataToVTK(model, "fluid_pressure", "fluid pressure", 0);

		// define some constant variables
		const double     hour(3600.0);
		double           time_increment(2.0 * hour); // timestep 2 hours

		// set the time increment for the FE algorithm
		pde_validate.TimeIncrement(1.0 / time_increment);
		pde_test.TimeIncrement(1.0 / time_increment);
		vtk_output.OutputDataToVTK(model, "fluid_pressure", "fluid pressure", 1);
		const std::vector<size_t>& index = pde_test.GetDOFIndex();
		const vector<double> & rh_validate = *pde_validate.GetRH();
		const vector<double> & rh_test = *pde_test.GetRH();


		// ====================== TESTING PROCESS ======================
		//
		// 1. test matrix establish and enumerate DOFs
		pde_validate.EstablishMatrixSetup(region);
		pde_test.EstablishMatrixSetupTest(region);

		size_t pressureDirchletDOFs(0);
		Index pressureKey = model.Database().StorageKey("fluid pressure");
		for ( auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter ) {
			if ((*nIter)->Status(pressureKey) == DIRICH) {
				pressureDirchletDOFs += 1;
			}
		}

		_test(pde_validate.GetRH()->size() == pde_test.GetRH()->size());

		pde_test.EnumerateAndFixMatrixSize(region);
		_test(pde_validate.GetRH()->size() == index.size());
		_test(pde_validate.GetRH()->size() == pde_test.GetRH()->size() + pressureDirchletDOFs);
		_test(pde_validate.GetG()->Rows() == pde_test.GetG()->Rows() + pressureDirchletDOFs);
		_test(pde_validate.GetG()->Cols() == pde_test.GetG()->Cols() + pressureDirchletDOFs);

		// 2. test accumulate
		pde_validate.Accumulate(region);
		pde_test.AccumulateTest(region);

		// 3. test AssignInitialConditions & LateAccumulate
		pde_validate.AssignInitialConditions(region);
		pde_test.AssignInitialConditionsTest(region);
		_test(pde_validate.Transient() == true);
		_test(pde_test.Transient() == true);
		pde_validate.LateAccumulate(region);
		pde_test.LateAccumulateTest(region);

		for (size_t i(0); i < index.size(); ++i) {
			if (index[i] != NULL_IDX) {
				// _test conductance matrix
				for (size_t j(0); j < index.size(); ++j) {
					if (index[j] != NULL_IDX) {
						_test(pde_test.GetG()->At(index[i], index[j]) == pde_validate.GetG()->At(i, j));
					}
				}
				// check load vector
				_test(rh_test[index[i]] == rh_validate[i]);
			}
		}

		pde_validate.AssignEssentialConditions(region);
		pde_test.AssignEssentialConditionsTest(region);
		pde_test.GetG()->OutForMatlab("pde_test");
		pde_validate.GetG()->OutForMatlab("pde_validate");
		outVector(rh_test, "rh_test");
		outVector(rh_validate, "rh_validate");

	}


	void Transient_Test::outVector(const vector<double>& vector, std::string file) {
		ofstream  ofs(file);
		long         prec;
		const long   digits(3);

		if (digits != 0) {
			ofs.setf(ios::scientific);
			prec = ofs.precision(digits);
		}


		for (auto& ditc : vector) {
			ofs << ditc << "\n";
		}

		if (digits != 0) {
			ofs.unsetf(ios::scientific);
			ofs.precision(prec);
		}

	}

	/*
	void Transient_Test::oldIntegrate(PDE_Integrator_UoM_Mock<2U, Region>& pde, Region<2U>& domain ) {
		// 1. configure algorithm
		pde.EstablishMatrixSetup(domain);

		// 2. Accumulation: Note that the conditions that pertain to the group must be input !                                 
		pde.Accumulate(domain);

		// 3. If the computation is transient initial conditions must be input into the righthand vector
		if (pde.Transient() == true) pde.AssignInitialConditions(domain);

		// 4. If the computation is transient initial conditions must be input into the righthand vector
		if (pde.Transient() == true) pde.LateAccumulate(domain);

		// 5. assign conditions like Dirichlet or Neumann boundary conditions etc.
		pde.AssignEssentialConditions(domain);

		// 6. diagnostics
		pde.GetG()->OutForMatlab("pde_validate");
		Transient_Test::outVector(*pde.GetRH(), "rh_validate");	
	}

	void Transient_Test::newIntegrate(PDE_Integrator_UoM_Mock<2U, Region>& pde, Region<2U>& domain) {

		// 1. configure algorithm
		// => this is importance, since in dynamic changing of DIRICHLET BCs i.e: coupling and decoupling process
		pde.EstablishMatrixSetupTest(domain);
		pde.EnumerateAndFixMatrixSize(domain);

		// 2. Accumulation: Note that the conditions that pertain to the group must be input !                                 
		pde.AccumulateTest(domain);

		// 3. If the computation is transient initial conditions must be input into the righthand vector
		if (pde.Transient() == true) pde.AssignInitialConditionsTest(domain);

		// 4. If the computation is transient initial conditions must be input into the righthand vector
		if (pde.Transient() == true) pde.LateAccumulateTest(domain);

		// 5. assign conditions like Dirichlet or Neumann boundary conditions etc.
		pde.AssignEssentialConditionsTest(domain);

		// 6. diagnostics
		pde.GetG()->OutForMatlab("pde_test");
		Transient_Test::outVector(*pde.GetRH(), "rh_test");
	}

*/

} // csmp
