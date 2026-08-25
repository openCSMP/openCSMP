#include "ANSYS_Model3D.h"
#include "Model.h"
#include "Boundary.h"
#include "VTU_Interface.h"
#include "PDE_Integrator.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_PT_op_dS.h"
#include "StressesAndStrains.h"
#include "ExtractTensorVariableComponent.h"
#include "BoreHole_stability_VerticalWell3D_VVCase.h"
#include "Test.h"
#include "VSet.h"
#include "VSetConverter.h"
#include "ANSYS_Interface.h"
#include "MohrCoulombFailure_Visitor.h"
#include "VTK_Interface.h"
#include"ScalarVariable.h"
#include "PropertyHandle.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "SAMG_Exception.h"
#else
#include "LinearSolver.h"
#endif

using namespace std;

namespace csmp {

BoreHole_stability_VerticalWell3D_VVCase::BoreHole_stability_VerticalWell3D_VVCase(const char* prefix)
  {
    this->setName("BoreHole_stability_VerticalWell3D_VVCase");
    prefix_=prefix;
  }



    // =========================================================================
    //  Borehole Stability — Vertical Well, 3D Linear Elastic Analysis
    // =========================================================================
    //
    //  This test simulates the stress state around a vertical borehole in a
    //  3D linear elastic rock mass under in-situ stress conditions.
    //
    //  The model geometry consists of:
    //    - Region "ROCK"          : the rock matrix (volume elements)
    //    - Region "BOUNDARY_WELL" : the borehole wall surface (face elements
    //                               shared between ROCK and the borehole void)
    //
    //  Boundary conditions:
    //    - Far-field stresses applied as Neumann tractions on the six outer
    //      model boundaries (RIGHT, LEFT, TOP, BOTTOM, FRONT, BACK)
    //    - Mud pressure applied as a Neumann traction on the borehole wall
    //      (BOUNDARY_WELL), representing the support provided by drilling fluid
    //    - Displacement constrained on the BOTTOM boundary to prevent rigid
    //      body motion (roller condition: no vertical displacement)
    //    - One bottom node fully fixed to prevent horizontal rigid body motion
    //
    //  In-situ stress state (compression positive, geotechnical convention):
    //    SHmax = 90.0 MPa  (maximum horizontal stress, x-direction)
    //    SHmin = 51.5 MPa  (minimum horizontal stress, z-direction)
    //    SV    = 88.2 MPa  (vertical overburden stress, y-direction)
    //    Pp    = 31.5 MPa  (pore pressure)
    //    Pmud  = 31.5 MPa  (mud pressure at the borehole wall)
    //
    //  Rock properties (carbonate):
    //    E     = 24.0 GPa  (Young's modulus)
    //    nu    = 0.29      (Poisson's ratio)
    //    rho*g = 2.4e4 N/m³ (bulk unit weight, y-direction downward)
    //
    //  Failure assessment (Mohr-Coulomb):
    //    cohesion        = 3.0 MPa
    //    friction angle  = 31 degrees
    //    tensile strength = 0.3 MPa
    //    Biot alpha      = 1.0 (fully coupled pore pressure)
    // =========================================================================

void BoreHole_stability_VerticalWell3D_VVCase::run()
{
    // --- in-situ stress and pressure parameters [Pa] ---
    const double SHmax( 90.0e6  );
    const double SHmin( 51.5e6  );
    const double Pmud ( 31.5e6  );
    const double SV   ( 88.2e6  );
    const double Pp   ( 31.5e6  );

    // --- zero initialisers ---
    const ScalarVariable      zeroScalar( PLAIN, 0. );
    const VectorVariable<3U>  zeroVector( PLAIN, 0. );

    // --- displacement constraint for BOTTOM boundary ---
    // roller condition: nodes can move horizontally but not vertically (y)
    const VectorVariable<3U> DisplacementVectorBOTTOM( PLAIN, DIRICH, PLAIN,
                                                        0., 0., 0. );

    // --- bulk unit weight vector [N/m³], acting downward in y-direction ---
    const VectorVariable<3U> Rho_g( PLAIN, PLAIN, PLAIN, 0., -2.4e4, 0. );

    // =========================================================================
    //  Model setup
    // =========================================================================

    string input_file_name( prefix_ );

    ANSYS_Model3D model( input_file_name.data(),
                         this->getName().c_str(),
                         "BoreHole_stability_VerticalWell3D_VVCase-variables.txt" );

    VTU_Interface<3U> vtu( model );
    vtu.OmitZeroInFileName( true );

    // --- initialise all output variables to zero before applying BCs ---
    model.InputPropertyValue( "mean stress",    zeroScalar );
    model.InputPropertyValue( "displacement",   zeroVector );
    model.InputPropertyValue( "force",          zeroVector );
    model.InputPropertyValue( "Neumann stress", zeroVector );
    model.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, Pp ) );

    // output the unloaded (zero stress) state for reference
    list<string> outputProps;
    vtu.OutputDataToVTU( input_file_name + "UnLoaded",
                         outputProps, "Model", static_cast<int>(0) );

    // =========================================================================
    //  Material properties
    // =========================================================================

    // rock matrix (carbonate)
    model.InputPropertyValue( "Young's modulus",  makeScalar( PLAIN, 24.0e9 ) );
    model.InputPropertyValue( "Poisson's ratio",  makeScalar( PLAIN, 0.29   ) );
    model.InputPropertyValue( "gravity force",    Rho_g                        );

    // failure criterion parameters — must be set before constructing the visitor
    model.InputPropertyValue( "cohesion",          makeScalar( PLAIN, 3.0e6  ) );
    model.InputPropertyValue( "friction angle",    makeScalar( PLAIN, 31.    ) );
    model.InputPropertyValue( "tensile strength",  makeScalar( PLAIN, 3.0e5  ) );
    model.InputPropertyValue( "Biot alpha",        makeScalar( PLAIN, 1.0    ) );

    // =========================================================================
    //  Boundary conditions
    // =========================================================================

    // retrieve boundary objects
    Boundary<3U>& rightBoundary  ( model.Boundary( "RIGHT"  ) );
    Boundary<3U>& leftBoundary   ( model.Boundary( "LEFT"   ) );
    Boundary<3U>& topBoundary    ( model.Boundary( "TOP"    ) );
    Boundary<3U>& bottomBoundary ( model.Boundary( "BOTTOM" ) );
    Boundary<3U>& frontBoundary  ( model.Boundary( "FRONT"  ) );
    Boundary<3U>& backBoundary   ( model.Boundary( "BACK"   ) );

    // --- far-field traction boundary conditions ---
    // Neumann tractions representing the in-situ stress state.
    // The NEUMANN flag is mandatory — AccumulateBoundaryIntegrals only
    // accumulates face integrals where the material operand is flagged NEUMANN.
    // Sign convention: compression positive (geotechnical).
    // The traction vector on each face points inward (into the model).
    const VectorVariable<3U> StressRight ( NEUMANN,NEUMANN,NEUMANN, -SHmax,  0.,     0.    );
    const VectorVariable<3U> StressLeft  ( NEUMANN,NEUMANN,NEUMANN,  SHmax,  0.,     0.    );
    const VectorVariable<3U> StressTop   ( NEUMANN,NEUMANN,NEUMANN,  0.,    -SV,     0.    );
    const VectorVariable<3U> Stressfront ( NEUMANN,NEUMANN,NEUMANN,  0.,     0.,    -SHmin );
    const VectorVariable<3U> Stressback  ( NEUMANN,NEUMANN,NEUMANN,  0.,     0.,     SHmin );

    rightBoundary.InputPropertyValue(  "Neumann stress", StressRight );
    leftBoundary.InputPropertyValue(   "Neumann stress", StressLeft  );
    topBoundary.InputPropertyValue(    "Neumann stress", StressTop   );
    frontBoundary.InputPropertyValue(  "Neumann stress", Stressfront );
    backBoundary.InputPropertyValue(   "Neumann stress", Stressback  );

    // --- mud pressure on the borehole wall ---
    // The borehole wall is represented by "BOUNDARY_WELL", a lower-dimensional
    // region of face elements shared between the rock matrix ("ROCK") and the
    // borehole void. Fluid pressure is not a primary variable — its effect on
    // the borehole wall is applied as a Neumann traction via the surface
    // integral operator NumIntegral_PT_op_dS, which integrates the pressure
    // over the borehole wall faces and adds the resulting force to the RHS.
    constexpr bool check_topo_attributes_of_nodes = false;
    pair<string,bool> wellBoundary = model.CreateExternalBoundaryFrom( "BOUNDARY_WELL", check_topo_attributes_of_nodes );
    // The NEUMANN flag is required for AccumulateBoundaryIntegrals to pick
    // up this condition.
    auto& wellboreBoundary( model.Boundary( wellBoundary.first ) );
    const ScalarVariable MudPressure( NEUMANN, Pmud );
    model.CreateProperty("mud pressure", "Pfm", "Pa", SCALAR, FACE );
    model.InputPropertyValue( "mud pressure", makeScalar(ANY,0.) );
    wellboreBoundary.InputPropertyValue( "fluid pressure", MudPressure );

    // --- displacement constraints ---
    // BOTTOM boundary: roller condition — no vertical (y) displacement.
    // Nodes are free to move horizontally, preventing stress locking.
    bottomBoundary.InputPropertyValue( "displacement", DisplacementVectorBOTTOM );

    // Fix one bottom node in all three directions to prevent horizontal
    // rigid body translation. Without this, the system matrix is singular.
    const csmp::Index displ_key = model.Database().StorageKey( "displacement" );
    bottomBoundary.N(0)->Store( displ_key,
                                makeVector( DIRICH,DIRICH,DIRICH, 0.,0.,0. ) );

    // =========================================================================
    //  Diagnostic output of input state
    // =========================================================================

    printRangeOfVariable( model, "Neumann stress" );

    const list<string> inputProps = {
        "Young's modulus", "Poisson's ratio", "displacement",
        "gravity force",   "Neumann stress",  "mean stress",
        "fluid pressure"
    };
    vtu.OutputDataToVTU( "BoreHole_stability_VerticalWell3D_VVCase_input",
                         inputProps, "ROCK", 0 );

    // =========================================================================
    //  Finite element solver setup
    // =========================================================================

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_napproach( 2 );
    SAMG_Solver solver( &settings );
    PDE_Integrator<3U,Element> deformation( solver );
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    PDE_Integrator<3U,Element> deformation( solver );
#endif

    // =========================================================================
    //  FE operators
    // =========================================================================

    // Stiffness matrix: ∫ Bᵀ D B dV
    // Assembles the elastic stiffness matrix from Young's modulus and
    // Poisson's ratio using the strain-displacement matrix B and the
    // constitutive matrix D.
    NumIntegral_BT_D_B_dV<3U> stiffness( model.Database(),
                                          "Young's modulus",
                                          "Poisson's ratio",
                                          "displacement",
                                          "displacement" );

    // Body force RHS: ∫ Nᵀ f dV
    // Integrates the gravitational body force vector over the element volume.
    // "gravity force" = rho * g stored as a VectorVariable at element level.
    // Note: PT_op (pre-computed nodal force) is not used here since "force"
    // is zero — only the numerically integrated body force is needed.
    NumIntegral_PT_op_dV<3U> bodyforce( model.Database(),
                                         "gravity force",
                                         "displacement" );

    // Far-field traction RHS: ∫ Nᵀ t dS  (outer boundaries)
    // Integrates the applied traction vector over the boundary faces.
    // Only faces where "Neumann stress" is flagged NEUMANN are accumulated.
    NumIntegral_PT_op_dS<3U> AppliedStress( model.Database(),
                                             "Neumann stress",
                                             "displacement" );

    // Wellbore pressure RHS: ∫ Nᵀ p dS  (borehole wall)
    // Integrates the mud pressure over the borehole wall faces.
    // Fluid pressure acts as an isotropic normal traction on the borehole
    // surface. Only faces where "fluid pressure" is flagged NEUMANN are
    // accumulated — hence the NEUMANN flag set on BOUNDARY_WELL above.
    NumIntegral_PT_op_dS<3U> WellBorePressure( model.Database(),
                                                "mud pressure",
                                                "displacement" );

    // register operators with the integrator
    deformation.Add( &stiffness );
    deformation.Add( &bodyforce );
    deformation.AddBoundaryIntegral( &AppliedStress );
    // deformation.AddBoundaryIntegral( &WellBorePressure );
    deformation.AddBoundaryIntegral( &WellBorePressure );

    // =========================================================================
    //  Post-processing: stress and strain invariants
    // =========================================================================

    // StressesAndStrains computes principal stresses (sigma1, sigma2, sigma3),
    // principal strains, mean stress, and von Mises stress from the computed
    // displacement field. Results are stored at element integration points.
    constexpr bool output_principal_vectors    { true  };
    constexpr bool extrapolate_results_to_nodes{ false };
    StressesAndStrains<3U> postpro( model,
                                    "Young's modulus",
                                    "Poisson's ratio",
                                    "displacement",
                                    output_principal_vectors,
                                    extrapolate_results_to_nodes );
    deformation.AddPostProcess( &postpro );

    // =========================================================================
    //  Solve
    // =========================================================================

    // IntegrateOver assembles and solves the global system over the entire
    // model, restricting element-level computations to the "Model" region.
    // Boundary integrals are accumulated from all registered boundaries.
    deformation.IntegrateOver( model, model.Region("ROCK") );

    // =========================================================================
    //  Failure assessment (Mohr-Coulomb)
    // =========================================================================

    // MohrCoulombFailure_Visitor evaluates shear and tensile failure criteria
    // at element integration points using the computed stress field and the
    // Biot effective stress convention:
    //   F_mc = p_eff * sin(phi) - c * cos(phi) + q * K(theta)
    // where p_eff = p - alpha * Pp is the Biot effective mean stress.
    // Results are stored in "failure" (F_mc value) and "tensile failure" (F01).
    MohrCoulombFailure_Visitor<3U> failureTest( model, ELEMENT );

    // restrict failure assessment to the rock matrix — the borehole void
    // ("WELL") does not contain meaningful stress values
    model.Region( "ROCK" ).Accept( failureTest );

    // =========================================================================
    //  Output
    // =========================================================================

    const set<string> output_variables {
        "strain1",       // minimum principal strain
        "strain2",       // intermediate principal strain
        "strain3",       // maximum principal strain
        "sigma1",        // maximum principal stress [Pa]
        "sigma2",        // intermediate principal stress [Pa]
        "sigma3",        // minimum principal stress [Pa]
        "mean stress",   // isotropic stress p = (sigma1+sigma2+sigma3)/3 [Pa]
        "displacement",  // displacement vector [m]
        "failure",       // Mohr-Coulomb criterion F_mc (>=0 means failure)
        "stress",        // full Cauchy stress tensor [Pa]
        "strain",        // full strain tensor [-]
        "failure01",     // tensile failure measure [Pa]
        "stress node",   // stress extrapolated to nodes [Pa]
        "strain node"    // strain extrapolated to nodes [-]
    };

    vtu.OutputDataToVTU( "BoreholeStability", output_variables, "ROCK", 0 );

} // end run

} // csmp
