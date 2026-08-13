#include "PDE_Integrator_Transient_Test.h"

// the CSMP model
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "ModelTime.h"

// a simple FE mesh generator
#include "Triangulator.h"

// the FE algorithm
#include "PDE_Integrator.h"

// PDE operators building the FE algorithm
// for comparison between analytic and numeric integration
#include "Integral_NT_op_N_dV.h"
#include "Integral_NT_lhsop_N_dV.h"
#include "Integral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

// output interfaces
#include "VTK_Interface.h"
#include "MatlabInterface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"
#include "ConstantFactor.h"

// PDE operators building the FE algorithm
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_dNT_dN_dV.h"

// FE grid generation
#include "Quadrilaterator.h"

using namespace std;

namespace csmp {

void PDE_Integrator_Transient_Test::run()
{
#ifdef CSMP_WITH_PETSC_SOLVER
      PetscInitializeNoArguments();
#endif

    // -----------------------------------------------------------------------
    // 1. Mesh and model setup
    // -----------------------------------------------------------------------
    double& model_time(ModelTime::Instance().modelTime);
    model_time = 0.0;

    Quadrilaterator quadrilaterator;
    VSet<2U>        mesh_container;
    const string    file_name("tutorial1_input");
    const double    x(10), y(10);

    quadrilaterator.QuadrilateralsFromRegularGrid(
        mesh_container, file_name.c_str(), x, y);

    Model<2U> model(mesh_container, "Tutorial1-variables.txt");
    const PropertyDatabase<2>& p_ref(model.Database());
    Region<2U>& region = model.Region("Model");

    // -----------------------------------------------------------------------
    // 2. Material properties and boundary conditions
    // -----------------------------------------------------------------------
    model.InputPropertyValue("porosity",        makeScalar(PLAIN, 0.1));
    model.InputPropertyValue("permeability",    makeScalar(PLAIN, 1.0e-15));
    model.InputPropertyValue("compressibility", makeScalar(PLAIN, 5.0e-10));
    model.InputPropertyValue("fluid pressure",  makeScalar(PLAIN, 1.0e+07));
    model.InputPropertyValue("fluid volume source", makeScalar(PLAIN, 0.0));

    model.InputBoundaryValue(LEFT,  "fluid pressure",
                             makeScalar(DIRICH, 3.0e+07));
    model.InputBoundaryValue(RIGHT, "fluid pressure",
                             makeScalar(DIRICH, 1.0e+07));

    // -----------------------------------------------------------------------
    // 3. Hydraulic conductivity
    // -----------------------------------------------------------------------
    ConstantFactor<2U, divides> conductivity(
        p_ref, "conductivity", "permeability", 0.001);
    model.Apply(conductivity);

    // -----------------------------------------------------------------------
    // 4. PDE operators
    // -----------------------------------------------------------------------
    NumIntegral_dNT_dN_dV<2U> stiffness_matrix(
        p_ref, "fluid pressure", "fluid pressure");

    NumIntegral_NT_lhsop_N_dV<2U> mass_matrix_lhs(
        p_ref, "compressibility", "fluid pressure", "fluid pressure");

    NumIntegral_NT_rhsop_N_dV<2U> mass_matrix_rhs(
        p_ref, "compressibility", "fluid pressure");

    NumIntegral_NT_rhsop_N_dV<2U> source_term(
        p_ref, "fluid volume source", "fluid pressure");

    mass_matrix_lhs.MultiplyWithTimeIncrement(true);
    mass_matrix_rhs.MultiplyWithTimeIncrement(true);
    mass_matrix_lhs.LumpedFormulation(true);
    mass_matrix_rhs.LumpedFormulation(true);
    source_term.LumpedFormulation(true);
    source_term.AddAccumulateLater();

    VelocityAndVolumeFlux<2U> velo( model, "conductivity", "porosity", "fluid pressure", true );

    // -----------------------------------------------------------------------
    // 5. Single PDE_Integrator — accessed via attorney
    // -----------------------------------------------------------------------
    // ([C] + dt[K]){p}t+dt = [C]{p}t + dt {Q}t+dt
    PDE_Integrator<2U,Element>  pde;
  #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver  samg_solver;
    pde.SetSolver( samg_solver );
  #else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    pde.SetSolver( linear_solver );
  #endif
    PDE_Integrator_Attorney1<2U> attorney(pde);

    attorney.Add(&stiffness_matrix);
    attorney.Add(&source_term);
    attorney.Add(&mass_matrix_lhs);
    attorney.Add(&mass_matrix_rhs);
    attorney.AddPostProcess(&velo);

    const double hour(3600.0);
    const double dt(2.0 * hour);
    attorney.TimeIncrement(1.0 / dt);

    // -----------------------------------------------------------------------
    // 6. Count Dirichlet DOFs
    // -----------------------------------------------------------------------
    const Index pressureKey = p_ref.StorageKey("fluid pressure");
    size_t n_dirich = 0;
    for ( auto nIter=region.NodesBegin(); nIter!=region.NodesEnd(); ++nIter )
        if ( (*nIter)->Status(pressureKey) == DIRICH )
            ++n_dirich;

    const size_t n_nodes = region.Nodes();
    const size_t n_free  = n_nodes - n_dirich;

    // -----------------------------------------------------------------------
    // 7. EstablishMatrixSetup
    // -----------------------------------------------------------------------
    _test( attorney.EstablishMatrixSetup(region) );

    // --- Structural tests ---
    _test( attorney.G_.Rows()    == n_free );
    _test( attorney.G_.Cols()    == n_free );
    _test( attorney.rh_.size()   == n_free );
    _test( attorney.G_.VerifySparsityPattern() );
    _test( attorney.DOF_indexes_.size() == n_nodes );

    // -----------------------------------------------------------------------
    // 8. Accumulate stiffness and mass matrices
    // -----------------------------------------------------------------------
    attorney.Accumulate(region);

    // --- Numerical tests on assembled matrix ---
    TestSymmetry( attorney, n_free );

    bool diag_positive = true;
    for ( size_t i=0; i<n_free; ++i )
        if ( attorney.G_.At(i,i) <= 0.0 )
            { diag_positive = false; break; }
    _test( diag_positive );

    TestLumpedMassMatrix( attorney, region, model, dt );

    // -----------------------------------------------------------------------
    // 9. AssignInitialConditions — populates RHS with p^n / dt
    // -----------------------------------------------------------------------
    attorney.AssignInitialConditions(region);
    _test( attorney.Transient() == true );

    TestInitialConditionRHS( attorney, region, model, dt );

    // RHS must be non-zero after initial conditions (p_initial = 1e7 Pa)
    bool rhs_nonzero = false;
    for ( size_t i=0; i<n_free; ++i )
        if ( std::abs(attorney.rh_[i]) > 0.0 )
            { rhs_nonzero = true; break; }
    _test( rhs_nonzero );

    // -----------------------------------------------------------------------
    // 10. LateAccumulate — adds source term
    // -----------------------------------------------------------------------
    attorney.LateAccumulate(region);

    // -----------------------------------------------------------------------
    // 11. AssignEssentialConditions — applies Dirichlet BCs to RHS
    // -----------------------------------------------------------------------
    attorney.AssignEssentialConditions(region);

    TestDirichletRHSModification( attorney, region, model );

    // RHS entries must be finite
    bool rhs_finite = true;
    for ( size_t i=0; i<n_free; ++i )
        if ( !std::isfinite(attorney.rh_[i]) )
            { rhs_finite = false; break; }
    _test( rhs_finite );

    // --- Steady state solution test ---
    TestSteadyState( attorney, region, model );

    // -----------------------------------------------------------------------
    // 12. VTK output (non-critical, just verify no crash)
    // -----------------------------------------------------------------------
    VTK_Interface<2U> vtk_output;
    vtk_output.OutputDataToVTK(model, "fluid_pressure", "fluid pressure", 0);

    cout << "\nPDE_Integrator_Transient_Test: all tests passed.\n";

#ifdef CSMP_WITH_PETSC_SOLVER
    PetscFinalize();
#endif

} // end run



  
  

// Verify steady state: solve with very large dt
// Solution should be linear: p(x) = 3e7 - 2e7 * x/L
void PDE_Integrator_Transient_Test::TestSteadyState( PDE_Integrator_Attorney1<2U>& attorney,
                                                     Region<2U>&              region,
                                                     Model<2U>&               model )
{
    // Use huge time increment -> mass term negligible
    constexpr double dt_large = 1.0e+20;
    attorney.TimeIncrement(1.0 / dt_large);
    attorney.EstablishMatrixSetup(region);
    attorney.Accumulate(region);
    attorney.AssignInitialConditions(region);
    attorney.LateAccumulate(region);
    attorney.AssignEssentialConditions(region);

    // 7. Solve linear algebraic system of equations
    attorney.Solve();

    // 8. Write results from the solution vector back to Model
    attorney.OutputResults( region );
                           
    // 9. Calculation of result-dependent properties
    attorney.PostProcess( region );

    // Check solution is linear between boundary values
    const Index pKey = model.Database().StorageKey("fluid pressure");
    const double p_left  = 3.0e+07;
    const double p_right = 1.0e+07;
    const double L       = 10.0;   // model length

    for ( auto n=region.NodesBegin(); n!=region.NodesEnd(); ++n )
    {
        if ( (*n)->Status(pKey) == DIRICH ) continue;
        const double x        = (*n)->x();
        const double p_exact  = p_left + (p_right - p_left) * x / L;
        const double p_solved = (*n)->Read(pKey);
        _equal( p_solved, p_exact, 1.0e+02 );   // 100 Pa tolerance
    }
}



void PDE_Integrator_Transient_Test::TestLumpedMassMatrix( const PDE_Integrator_Attorney1<2U>& attorney,
                                                          const Region<2U>&                   region,
                                                          const Model<2U>&                    model,
                                                          double                              dt )
{
    const Index cKey = model.Database().StorageKey("compressibility");
    const Index pKey = model.Database().StorageKey("fluid pressure");
    const double inv_dt = 1.0 / dt;

    for ( auto n=region.NodesBegin(); n!=region.NodesEnd(); ++n )
    {
        if ( (*n)->Status(pKey) == DIRICH ) continue;

        const size_t ridx = attorney.DOF_indexes_[(*n)->Idx()];
        if ( ridx == NULL_IDX ) continue;

        // Expected lumped mass diagonal = c * sum(vol/n_nodes) * inv_dt
        double expected_mass = 0.0;
        for ( auto el=region.CellsBegin(); el!=region.CellsEnd(); ++el )
        {
            bool has_node = false;
            for ( auto ni=(*el)->NodesBegin(); ni!=(*el)->NodesEnd(); ++ni )
                if ( (*ni)->Idx() == (*n)->Idx() )
                    { has_node = true; break; }
            if ( !has_node ) continue;

            const double c   = (*el)->Read(cKey);
            const double vol = (*el)->Volume();
            expected_mass += c * vol / static_cast<double>((*el)->Nodes());
        }
        expected_mass *= inv_dt;

        // The diagonal of G includes both stiffness and mass contributions.
        // Isolate mass by comparing G with and without mass matrix.
        // Here we just check the diagonal is at least as large as mass alone.
        _test( attorney.G_.At(ridx, ridx) >= expected_mass );
    }
}



void PDE_Integrator_Transient_Test::TestSymmetry( const PDE_Integrator_Attorney1<2U>& attorney,
                                                  size_t                               n_free )
{
    constexpr double tol = 1.0e-10;
    bool symmetric = true;
    for ( size_t i=0; i<n_free; ++i )
        for ( size_t j=0; j<n_free; ++j )
        {
            const double diff = std::abs(
                attorney.G_.At(i,j) - attorney.G_.At(j,i));
            if ( diff > tol )
            {
                std::cout << "Symmetry violation at ("
                          << i << "," << j << "): "
                          << attorney.G_.At(i,j) << " vs "
                          << attorney.G_.At(j,i) << "\n";
                symmetric = false;
            }
        }
    _test( symmetric );
}



void PDE_Integrator_Transient_Test::TestInitialConditionRHS( const PDE_Integrator_Attorney1<2U>& attorney,
                               const Region<2U>&                   region,
                               const Model<2U>&                    model,
                               double                              dt )
{
    const Index pKey = model.Database().StorageKey("fluid pressure");
    const Index cKey = model.Database().StorageKey("compressibility");
    const double p_initial = 1.0e+07;
    const double inv_dt    = 1.0 / dt;

    for ( auto n=region.NodesBegin(); n!=region.NodesEnd(); ++n )
    {
        if ( (*n)->Status(pKey) == DIRICH ) continue;

        const size_t ridx = attorney.DOF_indexes_[(*n)->Idx()];
        if ( ridx == NULL_IDX ) continue;

        // Lumped mass * p_initial * inv_dt
        double expected = 0.0;
        for ( auto el=region.CellsBegin(); el!=region.CellsEnd(); ++el )
        {
            bool has_node = false;
            for ( auto ni=(*el)->NodesBegin(); ni!=(*el)->NodesEnd(); ++ni )
                if ( (*ni)->Idx() == (*n)->Idx() )
                    { has_node = true; break; }
            if ( !has_node ) continue;

            const double c   = (*el)->Read(cKey);
            const double vol = (*el)->Volume();
            expected += c * vol / static_cast<double>((*el)->Nodes());
        }
        expected *= inv_dt * p_initial;

        _equal( attorney.rh_[ridx], expected, expected * 1.0e-6 );
    }
}



void PDE_Integrator_Transient_Test::TestDirichletRHSModification( const PDE_Integrator_Attorney1<2U>& attorney,
                                    const Region<2U>&                   region,
                                    const Model<2U>&                    model )
{
    const Index pKey = model.Database().StorageKey("fluid pressure");

    // For each free node adjacent to a Dirichlet node,
    // the RHS must have been modified (non-zero contribution from BC)
    bool found_modification = false;
    for ( auto n=region.NodesBegin(); n!=region.NodesEnd(); ++n )
    {
        if ( (*n)->Status(pKey) == DIRICH ) continue;
        const size_t ridx = attorney.DOF_indexes_[(*n)->Idx()];
        if ( ridx == NULL_IDX ) continue;

        // Check if this node has a Dirichlet neighbour
        for ( auto el=region.CellsBegin(); el!=region.CellsEnd(); ++el )
        {
            bool has_free = false, has_dirich = false;
            for ( auto ni=(*el)->NodesBegin(); ni!=(*el)->NodesEnd(); ++ni )
            {
                if ( (*ni)->Idx() == (*n)->Idx() ) has_free = true;
                if ( (*ni)->Status(pKey) == DIRICH ) has_dirich = true;
            }
            if ( has_free && has_dirich )
            {
                // This node's RHS should have a non-trivial value
                // from both IC and Dirichlet modification
                _test( std::abs(attorney.rh_[ridx]) > 0.0 );
                found_modification = true;
                break;
            }
        }
    }
    _test( found_modification );   // at least one such node must exist
}



void PDE_Integrator_Transient_Test::OutVector(const vector<double>& vector, std::string file) {
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

	} // end OutVector
  


} // csmp
