#include "Tractions_Example.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"

#include "ANSYS_Model3D.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"

#include "InputDataManager.h"
#include "ComputationalSettings.h"

#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"  // conductance matrix (LHS)
#include "NumIntegral_NT_rhsop_N_dV.h" // fluid volume source
#include "NumIntegral_NT_op_N_dS.h"    // boundary integral

#include "VelocityAndVolumeFlux.h"     // post-processing of Darcy velocity

#include "CSMP_EigenGlue.h" // conversion between CSMP and Eigen vectors and matrices
#include "compareFloats.h"

#include "VTK_Interface.h"


using namespace std;

namespace csmp {

void Tractions_Example::Specifications()
  {
     SetTitle( "Tractions_Example" );
     SetDifficulty( 1 );
     SetCategory( "Simulation of Physical Processes" );
     AddAuthor( "SKM" );
     AddDescription( "Computes 'fluid pressure' and 'pore velocity' (shear stress) induced surface tractions and forces acting on a sphere in a user-specified fluid velocity field" );
     AddDescription( "source in: Tractions_Example.cpp" );
     AddRequirement( "user specified flow field" );
     AddRequirement( "Model in one_sphere_0.45_tetra" );
  }


void Tractions_Example::Run()
 {
  // Testing the components of the traction calculation (will be executed in debug mode only);
  assert( Test_NormalStressCalculation() == true );
  assert( Test_ViscousShearStressCalculation() == true );
  assert( Test_TractionsCalculation() == true );
  assert( Test_ComputationOfDragForces() == true );
  assert( Test_NoSlipTangentialTractionComputation() == true );
  // illustration / unit test
  DragForceOnInclinedPlane();
 
 
  // 0. Building computational model from ANSYS input files
  // ------------------------------------------------------
  constexpr uint32_t dim{3};
  const string       model_name("one_sphere_0.45_tetra");
  ANSYS_Model3D      model( model_name.c_str(), "Tractions_Example-variables.txt");
  printModelDimensions(model);


  // 1. Assigning properties from file
  // ---------------------------------
  /* NOTE:
     - the sphere in the model is a hole. Thus, region "GRAIN" is an outside surface
       which gets converted into a BOUNDARY object, replacing surface Elements with Faces that
       have access to their higher-dimensional neighbors in the fluid domain.
  */
  const double k1D{1.0e-12};
  model.InputPropertyValue("permeability", makeScalar(ANY,k1D) );
  model.InputPropertyValue("porosity", makeScalar(ANY,1.0) );
  const double dynamic_viscosity{1.6e-3}; // of water 1.6 centi-Poise (in Pa.s)
  model.InputPropertyValue("conductivity", makeScalar(ANY,k1D/dynamic_viscosity) );
  model.InputPropertyValue("fluid volume source", makeScalar(ANY,0.) );
  // induce left->right fluid pressure gradient
  model.Boundary("LEFT").InputPropertyValue("fluid pressure", makeScalar(DIRICH,2. * 100325.) ); // patm
  model.Boundary("RIGHT").InputPropertyValue("fluid pressure", makeScalar(DIRICH, 100325.) ); // patm


  // 2. Fluid pressure computation algorithm and velocity post-processing
  // --------------------------------------------------------------------
  // (Darcy velocities for simplicity as opposed to more correct Stokes-flow velocities)
  PDE_Integrator<dim,Element>  fluid_pressure;
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Solver  samg_solver;
  fluid_pressure.SetSolver( samg_solver );
#else
  CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
  fluid_pressure.SetSolver( linear_solver );
#endif

  NumIntegral_dNT_lhsop_dN_dV<dim> conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
  NumIntegral_NT_rhsop_N_dV<dim> source( model.Database(), "fluid volume source",  "fluid pressure" );

  VelocityAndVolumeFlux<dim>  velocity( model,  "conductivity", "porosity", "fluid pressure", false );

  fluid_pressure.Add( &conductance );
  fluid_pressure.Add( &source );
  fluid_pressure.AddPostProcess( &velocity );
  
  Region<dim>& flow_domain = model.Region("PORES");
  // computation and post processing of flow velocities
  fluid_pressure.IntegrateOver( flow_domain );
  // check
  printRangeOfVariable( model, "fluid pressure" );
  printRangeOfVariable( model, "velocity" );
  printRangeOfVariable( model, "pore velocity" );
  printRangeOfVariable( model, "volume flux" );

  VTK_Interface<dim>  vtk_output;
  vtk_output.OutputDataToVTK( model, "PORES", "fluid-pressure", "fluid pressure", 0, true );
  vtk_output.OutputDataToVTK( model, "PORES", "velocity", "pore velocity", 0, true );


  // 3. Calculation of forces acting on the sphere due to pressure-gradient and shear forces
  // ---------------------------------------------------------------------------------------
  // converting the surface region 'GRAIN' into a boundary
  const bool check_topo_attributes_of_nodes{false};
  pair<string,bool>  boundary = model.CreateExternalBoundaryFrom( "GRAIN", check_topo_attributes_of_nodes );
  assert( boundary.second == true );
  cout <<"\nName of newly created boundary: "<< boundary.first << endl;
  
  // (result should be equivalent to pressure-gradient * cross-sectional area of spherical grain, i.e. 0.45 * grad_p. Y and Z components of force vector must be zero)
  pressureForceActingOnObject( model, boundary.first );


  // 4. Calculating normal and shear tractions acting on boundary + the torque that is arising from these
  // ----------------------------------------------------------------------------------------------------
// Normal forces only:  ForceAndTorque<3> tractions = computeNormalAndShearTractionsOnBoundary( model, "GRAIN_BOUNDARY", dynamic_viscosity );
   const Boundary<dim>& grain_boundary = model.Boundary( boundary.first ); // "GRAIN_BOUNDARY"
   Eigen::Matrix<double,dim,1> obj_centroid = eigen_glue::ToEigen(grain_boundary.Centroid());
   Eigen::Matrix<double,dim,1> total_force  = Eigen::Matrix<double,dim,1>::Zero();
   Eigen::Matrix<double,dim,1> total_torque = Eigen::Matrix<double,dim,1>::Zero();
    
   const csmp::Index   p_key = model.Database().StorageKey("fluid pressure");
   const csmp::Index   u_key = model.Database().StorageKey("pore velocity");
   
   const double fluid_viscosity{1.6e-3}; // Pa.s = water
    
   for ( const auto& face : grain_boundary.CellVector() )
     {
        // computation of tractions:
        // --- Viscous stress: mu * (∇u + (∇u)^T);   ∇u is a dim * dim, since u is a vector
        // --- Full traction tensor, tau = -pI + tau_visc
        auto tractions = integrateFluidTractionOnFaceUsingCellBarycentre<dim>( face, u_key, p_key, fluid_viscosity, obj_centroid );
        
        // --- Tractions
        total_force += tractions.force;
        
        // --- Torque, M (lever: bctr - subdomain_centroid)
        total_torque += tractions.torque;
     }
     
  cout <<"\n"<<" Summary actions on '"<< grain_boundary.Name() <<"'"<< endl;
  cout <<"\t"<<"tractions:"<< endl;
  cout << total_force;
  cout <<"\n\t"<<"torque:"<< endl;
  cout << total_torque << endl;
  
  // TODO: display shear and normal forces separately
  // TODO: test against analytic solution (making this a formal "unit test"

  cout <<"\nRun: end of program"<< endl;
 
} // end Run



// ======================================================================================================
//
//     DEVELOPMENT OF FINITE ELEMENT FORMULATION FOR SEEPAGE FORCES
//
// ======================================================================================================


// COMPONENT TESTS AND CALCULATION EXAMPLES

/**
    Using Eigen matrices and vectors, tests the computation of pressure related normal stress acting on a surface.
    With Z pointing upward, the NE unit vector in the x–y plane is  eNE = {\sqrt2}(1,1,0).  Since the  surface normal is tilted by 45 degrees from the vertical toward NE,
    the unit normal is
                   un = cos(theta) z } +  sin(theta) eNE = (0.5,0.5,1/sqrt(2))

    With a velocity field  u(x) = alpha ( n . x) n, and alpha = 0.1 per second
    the velocity gradient is the outer product, nabla u = alpha n cross n.

 */
bool Tractions_Example::Test_NormalStressCalculation()
 {
    constexpr uint32_t dim{3u};
    const double       pf{1.0e7};
    Point<dim>         nrml(0.,1.,0.);
    Eigen::Matrix<double,dim,dim> pressure_stress = Eigen::Matrix<double,dim,dim>::Identity() * pf;
    // pressure_traction = pressure_stress * nrml;
    // cout <<"\n"<<"face-normal stress at ip "<< ip <<":\n"<< pressure_stress * nrml;
    auto normal_stress = pressure_stress * eigen_glue::ToEigen(nrml);
    cout <<"\n"<<"face-normal stress:\n";
    cout << normal_stress;
    if ( !approximatelyEqual( normal_stress(1), pf) ) return false;
    return true;
 }


/**
   Test case of a surface inclined by 45° toward the NE and a flow velocity that increases along the surface normal at 0.1 s⁻¹ (i.e., 10 cm/s per meter along the normal).
 */
bool Tractions_Example::Test_ViscousShearStressCalculation()
 {
    const double alpha = 0.1;    // 10 cm/s per m along the normal
    const double s2 = sqrt(0.5); // 1/sqrt(2)

    // Surface normal tilted 45° from vertical toward NE
    Eigen::Vector3d n(0.5, 0.5, s2); // (1/2, 1/2, 1/sqrt(2))
    // (already unit length: sqrt(0.25+0.25+0.5)=1)

    // ∇u = gradU = alpha * (n ⊗ n)
    Eigen::Matrix3d gradU = alpha * (n * n.transpose());

    cout << "n = " << n.transpose() << "\n";
    cout << "gradU =\n" << gradU << "\n";

    // viscous stress tau = mu * (gradU + gradU^T) = mu * (∇u + (∇u)^T)
    const double mu = 1.0;  // choose any viscosity for testing
    Eigen::Matrix3d tau = mu * (gradU + gradU.transpose());
    cout << "tau (mu=1) =\n" << tau << "\n";
    
    // full Cauchy stress tensor: sigma = -pI + tau
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double pf{1.0e7};
    Eigen::Matrix3d sigma = -pf * I + tau;

    // Output
    cout << "Surface normal n = " << n.transpose() << "\n\n";
    cout << "Velocity gradient gradU =\n" << gradU << "\n\n";
    cout << "Viscous stress tau =\n" << tau << "\n\n";
    cout << "Full Cauchy stress sigma =\n" << sigma << "\n";
 
    return true;
 }

/**
    Shear and normal tractions acting on the 45-degree inclined surface
 */
bool Tractions_Example::Test_TractionsCalculation()
 {
    // Parameters
    const double alpha = 0.1;     // velocity gradient rate (1/s)
    const double mu    = 1.0;     // viscosity (Pa·s)
    const double p     = 1.0;     // pressure (Pa) - test value

    // Surface normal tilted 45° from vertical toward NE (unit length)
    Eigen::Vector3d n(0.5, 0.5, sqrt(0.5));   // (1/2, 1/2, 1/sqrt(2))
    n.normalize(); // ensure unit length (idempotent here)

    // Velocity gradient: gradU = alpha * (n ⊗ n), note: tensor product symbol (\otimes in LaTeX)
    Eigen::Matrix3d gradU = alpha * (n * n.transpose());

    // Symmetric viscous stress: tau = mu * (gradU + gradU^T)
    Eigen::Matrix3d tau = mu * (gradU + gradU.transpose());

    // Full Cauchy stress: sigma = -p * I + tau
    Eigen::Matrix3d sigma = -p * Eigen::Matrix3d::Identity() + tau;

    // Traction vector on the surface (force per unit area)
    Eigen::Vector3d traction = sigma * n;
    //              ^^^^^^^^^^^^^^^^^^^^^
    // Decompose into normal and tangential parts
    double t_normal_mag = traction.dot(n);        // signed normal traction (positive in n-direction)
    Eigen::Vector3d t_normal = t_normal_mag * n;  // vector normal component
    Eigen::Vector3d t_tangential = traction - t_normal;

    // Print results
    cout << fixed << setprecision(6);
    cout << "Surface normal n           : " << n.transpose() << "\n\n";
    cout << "Velocity gradient gradU    :\n" << gradU << "\n\n";
    cout << "Viscous stress tau         :\n" << tau << "\n\n";
    cout << "Full Cauchy stress sigma   :\n" << sigma << "\n\n";
    cout << "Traction t = sigma * n     : " << traction.transpose() << "\n";
    cout << "  - normal component t_n   : " << t_normal.transpose()
         << "  (magnitude = " << t_normal.norm() << ")\n";
    cout << "  - tangential component t_t: " << t_tangential.transpose()
         << "  (magnitude = " << t_tangential.norm() << ")\n";

    return true;
 }



/**
     Introduction of a no-slip boundary condition on a horizontal surface, U=0 (z=0).
     while far from the surface, the fluid velocity approaches some value U.
     That implies a velocity gradient normal to the surface, which creates the viscous shear traction.
     
     Fluid above (z>0) flows in x direction with velocity increasing linearly with z.
     Let Kdelta be a flow boundary-layer thickness, then gradU(0,2) = U/Kdelta (for constant flow velo all other entries in gradU are zero).
     The wall traction, t,  that arises from this is,
     
         t = sigma * n = tau * n =  [muUKdelta, 0, 0 ].
     
     Its calculation involves the "rate-of-strain tensor"
      
         D = 1/2(gradU + gradU^T) * 1/2 * [viscous stress]
         
     The viscous stress matrix, tau is symmetric, and valued 2 * mu * D, ie.,
     
            [ 0 0 muU/Kdelta ]
      tau = [ 0 0     0      ]
            [ muU/Kdelta 0 0 ]
            
     The traction tangential to the wall is called the drag force, which, in this example,
     is acting parallel to the x direction.
 */
bool Tractions_Example::Test_ComputationOfDragForces()
 {
    double U = 0.1;      // m/s far-field velocity
    double delta = 0.01; // m, boundary layer thickness
    double mu = 1.0e-3;  // Pa·s, viscosity of water

    // Velocity gradient
    Eigen::Matrix3d gradU = Eigen::Matrix3d::Zero();
    gradU(0,2) = U/delta; // ∂u_x/∂z

    // Rate-of-strain tensor
    Eigen::Matrix3d D = 0.5 * (gradU + gradU.transpose());

    // Viscous stress
    Eigen::Matrix3d tau = 2.0 * mu * D;

    // Wall normal
    Eigen::Vector3d n(0,0,1);

    // Traction on the wall
    Eigen::Vector3d t = tau * n;

    cout << "gradU:\n" << gradU << "\n\n";
    cout << "D:\n" << D << "\n\n";
    cout << "tau:\n" << tau << "\n\n";
    // drag force per unit area acting on the surface
    cout << "Traction vector:\n" << t.transpose() << endl;

    // test: drag-force vector should only have a non-zero entry in the x direction
    if ( !approximatelyEqual(t(0),0.01) || !approximatelyEqual(t(1),0.) || !approximatelyEqual(t(2),0.) ) return false;

    return true;
 }





/**
    Eigen example that implements a no-slip wall on an inclined surface (tilted 45° toward NE),
    chooses a sensible tangential flow direction along the surface,
    imposes a linear velocity profile normal to the surface (wall at n·x=0 is stationary).

    Example computes the velocity gradient, the viscous stress, the full Cauchy stress (with pressure),
    and the traction on the inclined surface. It also decomposes the traction into normal and tangential components
    and prints the results.

    Key formula used (with alpha = U/δ):
      •	 velocity (u) field: \mathbf{u}(\mathbf{x}) = \alpha ( \mathbf{n}\cdot\mathbf{x} ) \; \mathbf{t}
        where \mathbf{n} is the unit surface normal and \mathbf{t} is a unit tangent direction.
        
      •	velocity gradient: \nabla\mathbf{u} = \alpha \; \mathbf{t}\otimes\mathbf{n}
      
      •	viscous stress: \tau = \mu(\nabla \mathbf{u} + \nabla\mathbf{u}^\top)
      
      •	full Cauchy stress: \sigma = -pI + \tau
      
      •	traction: \mathbf{T} = \sigma \cdot \mathbf{n}
 */
void Tractions_Example::DragForceOnInclinedPlane()
 {
    using Eigen::Vector3d;
    using Eigen::Matrix3d;
    cout << fixed << setprecision(8);

    // Physical parameters
    const double U     = 0.1;      // m/s far-field tangential velocity scale
    const double delta = 0.01;     // m, characteristic thickness where u rises from 0 to U
    const double mu    = 1.0e-3;   // Pa·s (viscosity)
    const double p     = 1.0;      // Pa (example pressure)
    const double alpha = U / delta; // velocity gradient scale (1/s)

    // Surface normal: 45° tilt toward NE (as in earlier messages)
    Vector3d n(0.5, 0.5, sqrt(0.5)); // (1/2, 1/2, 1/sqrt(2))
    n.normalize(); // defensive normalization

    // Choose a tangential flow direction: project global x-axis onto the surface tangent plane
    Vector3d ex(1.0, 0.0, 0.0);             // flow direction
    Vector3d t_hat = ex - (ex.dot(n)) * n;  // remove normal velocity component
    double t_norm = t_hat.norm();
    if (t_norm < 1e-12) {
        // fallback: use global y if ex is parallel to n (very unlikely here)
        Vector3d ey(0.0, 1.0, 0.0);
        t_hat = ey - (ey.dot(n)) * n;
        t_norm = t_hat.norm();
        if (t_norm < 1e-12) {
            cerr << "Could not find a tangent direction (degenerate normal)." << endl;
            return;
        }
    }
    t_hat /= t_norm; // unit tangent vector

    // Velocity gradient for u(x) = alpha * (n·x) * t_hat  -> gradU = alpha * (t_hat ⊗ n)
    Matrix3d gradU = alpha * (t_hat * n.transpose()); // outer product t_hat n^T

    // Symmetric part (rate-of-strain) and viscous stress
    Matrix3d D   = 0.5 * (gradU + gradU.transpose());
    Matrix3d tau = mu * (gradU + gradU.transpose()); // == 2*mu*D

    // Full Cauchy stress sigma = -p I + tau
    Matrix3d sigma = -p * Matrix3d::Identity() + tau;

    // Traction on the surface: T = sigma * n
    Vector3d traction = sigma * n;

    // Decompose traction into normal and tangential components
    double t_normal_scalar = traction.dot(n);           // scalar normal traction (signed)
    Vector3d t_normal_vec = t_normal_scalar * n;        // vector normal comp
    Vector3d t_tangential = traction - t_normal_vec;    // vector tangential comp

    // For diagnostics also compute tangential component magnitude along chosen tangent
    double t_tangential_along_t = t_tangential.dot(t_hat);

    // Output
    cout << "Parameters:\n";
    cout << " U      = " << U << " m/s\n";
    cout << " delta  = " << delta << " m\n";
    cout << " alpha  = U/delta = " << alpha << " 1/s\n";
    cout << " mu     = " << mu << " Pa*s\n";
    cout << " p      = " << p << " Pa\n\n";

    cout << "Surface normal n       = " << n.transpose() << "\n";
    cout << "Surface tangent t_hat  = " << t_hat.transpose() << " (unit)\n\n";

    cout << "Velocity gradient gradU =\n" << gradU << "\n\n";
    cout << "Rate-of-strain D =\n" << D << "\n\n";
    cout << "Viscous stress tau =\n" << tau << "\n\n";
    cout << "Full Cauchy stress sigma =\n" << sigma << "\n\n";

    cout << "Traction T = sigma * n = " << traction.transpose() << "\n";
    cout << "  - normal component (scalar)    : " << t_normal_scalar << "\n";
    cout << "  - normal component (vector)    : " << t_normal_vec.transpose() << "\n";
    cout << "  - tangential component (vector): " << t_tangential.transpose()
              << "  (magnitude = " << t_tangential.norm() << ")\n";
    cout << "  - tangential component along chosen t_hat: " << t_tangential_along_t << "\n";

    // If you want force on a patch of area A at this traction:
    double A = 0.5; // example area (m^2)
    Vector3d force_on_patch = traction * A;
    Vector3d moment_about_origin = ( /* use centroid position r (example) */ Vector3d(1.0,0.0,0.0) ).cross(force_on_patch);
    cout << "\nExample patch area A = " << A << " m^2\n";
    cout << "  -> Force on patch = " << force_on_patch.transpose() << " N\n";
    cout << "  -> Moment about origin (for centroid at [1,0,0]) = " << moment_about_origin.transpose() << " N·m\n";
    
} // end DragForceOnInclinedPlane







/**
  Computes fluid induced shear stress acting on surface with no-slip boundary conditions

  1.	User provides the flow velocity U at a point (with position) near the surface.
	2.	We know the surface normal n
	3.	We estimate the tangential gradient based on the shortest distance of point from surface

\f[ \frac{\partial \mathbf{u}t}{\partial n} \approx \frac{\mathbf{u}\text{sample} -
    \mathbf{u}\text{wall}}{d} = \frac{\mathbf{u}\text{sample}}{d} \quad \text{(assuming no-slip wall: } \mathbf{u}_\text{wall}=0\text{)} \f]
    
  4.	Construct a rank-one velocity gradient:

\nabla \mathbf{u} \approx \frac{\mathbf{u}_\text{tangent} \otimes \mathbf{n}}{d}
	•	u_tangent = u_sample - (u_sample·n) n ensures only tangential components contribute.

  @attention Normal component is ignored because it does not produce tangential traction.
  
  This creates its own rank-one approximation of gradU = u_tangent ⊗ n / d i.e., the dominant shear near the wall, which produces tangential traction.
	•	Works for any actual velocity vector measured near the surface.
	•	Normal velocity component does not contribute to tangential shear, only to normal stress.
	•	For the FE calculation, the traction must be integrated over the Face objects to get total force and torque.
  
  Obviously, in this planar surface demo, the torque computation is not included. To compute the torque just use the centroid of the object, x_ref in:
  
      torque += (x - x_ref).cross(traction * dS);
*/
TractionResult Tractions_Example::FlowInducedTangentialTraction( const Eigen::Vector3d& u_sample, // actual velocity at point
                                                                 const Eigen::Vector3d& n,        // surface normal
                                                                 double d,                        // distance from surface along n
                                                                 double mu,                       // dynamic viscosity
                                                                 double p )                       // pressure at point
{
    TractionResult res;

    Eigen::Vector3d n_unit = n.normalized();

    // tangential velocity component
    Eigen::Vector3d u_tangent = u_sample - (u_sample.dot(n_unit)) * n_unit;

    // approximate velocity gradient (rank-one)
    Eigen::Matrix3d gradU = (1.0/d) * (u_tangent * n_unit.transpose());

    // rate-of-strain tensor
    Eigen::Matrix3d D = 0.5 * (gradU + gradU.transpose());

    // viscous stress
    Eigen::Matrix3d tau = 2.0 * mu * D;

    // full Cauchy stress
    Eigen::Matrix3d sigma = -p * Eigen::Matrix3d::Identity() + tau;

    // traction vector
    Eigen::Vector3d traction = sigma * n_unit;
    Eigen::Vector3d t_normal = (traction.dot(n_unit)) * n_unit;
    Eigen::Vector3d t_tangential = traction - t_normal;

    // fill struct
    res.velocity = u_sample;
    res.gradU = gradU;
    res.viscousStress = tau;
    res.cauchyStress = sigma;
    res.traction = traction;
    res.tractionNormal = t_normal;
    res.tractionTangential = t_tangential;

    return res;
}

/**
    Tests FlowInducedTangentialTraction()  for an example of flow past a plane inclined 45o to the NE
 */
bool Tractions_Example::Test_NoSlipTangentialTractionComputation()
 {
    using namespace Eigen;
    std::cout << std::fixed << std::setprecision(6);

    // Example: arbitrary surface and sample velocity
    Vector3d n(0.5, 0.5, std::sqrt(0.5));
    Vector3d u_sample(0.05, 0.02, 0.0); // measured near-surface flow
    double d = 0.005; // distance from wall (m)
    double mu = 1e-3;
    double p = 1.0;

    auto result = FlowInducedTangentialTraction(u_sample, n, d, mu, p);
    //TODO: write test

    std::cout << "Sample velocity: " << result.velocity.transpose() << "\n";
    std::cout << "Velocity gradient gradU:\n" << result.gradU << "\n";
    std::cout << "Viscous stress tau:\n" << result.viscousStress << "\n";
    std::cout << "Cauchy stress sigma:\n" << result.cauchyStress << "\n";
    std::cout << "Traction vector: " << result.traction.transpose() << "\n";
    std::cout << "  Normal component: " << result.tractionNormal.transpose() << "\n";
    std::cout << "  Tangential component: " << result.tractionTangential.transpose() 
              << " (magnitude=" << result.tractionTangential.norm() << ")\n";
              
    return true;
}



// ===================================================================================================
//
//          INTEGRATION INTO CSMP++
//
// ===================================================================================================


/** Integrates viscous + pressure tractions over a face using
    the barycentre of the parent cell for velocity gradients.

    @param face       The boundary face (with outward-pointing normal)
    @param vel_key    Index for velocity (VectorVariable, NODE placement)
    @param press_key  Index for pressure (ScalarVariable, NODE or ELEMENT)
    @param viscosity  Dynamic viscosity μ
    @param obj_centroid  Reference point for torque calculation
*/
template<uint32_t dim>
ForceAndTorque<dim> integrateFluidTractionOnFaceUsingCellBarycentre(
                                            const Face<dim>* face,
                                            const csmp::Index& vel_key,
                                            const csmp::Index& press_key,
                                            double viscosity,
                                            const Eigen::Matrix<double,dim,1>& obj_centroid )
{
    auto& csmp_error( ErrorHandler::Instance() );

    ForceAndTorque<dim> result;

    // parent element
    // TODO: InnerParent() works only for one-sphere example; change to OuterParent() for volume-meshed particle
    const auto* cell = face->InnerParent();
    if (!cell) {
        csmp_error.Note(ERROR, "IntegrateFluidTractionOnFaceUsingCellBarycentre",
                                                "Face has no inner parent cell.");
        return ForceAndTorque<dim>();
    }
    if ( press_key.place != NODE ) {
        csmp_error.Note(ERROR,
            "IntegrateFluidTractionOnFaceUsingCellBarycentre",
            "Unsupported placement for pressure variable.");
        return ForceAndTorque<dim>();
    }

    // velocity gradient at (volumetric) parent element barycentre
    double detJ{0.};
    Eigen::Matrix<double,dim,dim> gradU = cell->PropertyGradientAtBaryCenter( vel_key, detJ );

    // NOTE: we must still loop over the integration points of the (surface) Face to complete the FE integral
    const uint32_t nip = face->IntegrationPoints();
    for (uint32_t ip = 0; ip < nip; ++ip)
      {
        // mapped Face area
        detJ = face->det_J_AtIntegrationPoint(ip);

        // quadrature weight
        double iweight = face->WeightAtIntegrationPoint(ip);

        // outward unit normal
        Eigen::Matrix<double,dim,1> normal = eigen_glue::ToEigen(face->UnitNormal());

        // barycentre coordinates
        Eigen::Matrix<double,dim,1> bctr = eigen_glue::ToEigen(face->BaryCenter());

        // fluid pressure at Face barycentre (interpolated)
        double pressure = face->PropertyValueAtIntegrationPoint( press_key, ip );

        // stress tensor: -p I + μ (∇u + ∇u^T)
        Eigen::Matrix<double,dim,dim> sigma =
              (-pressure) * Eigen::Matrix<double,dim,dim>::Identity()
            + viscosity * (gradU + gradU.transpose());

        // traction t = σ ⋅ n
        Eigen::Matrix<double,dim,1> traction = sigma * normal;

        // area contribution
        double dA = iweight * detJ;

        // accumulate force
        result.force += traction * dA;

        // accumulate torque
        Eigen::Matrix<double,dim,1> r = bctr - obj_centroid;
        result.torque += r.cross(traction) * dA;
    }

    return result;
    
} // end IntegrateFluidTractionOnFaceUsingCellBarycentre

// explicit function template instantiation
template ForceAndTorque<3> integrateFluidTractionOnFaceUsingCellBarycentre<3>( const Face<3>*, const Index&, const Index&, double, const Eigen::Matrix<double,3,1>& );



/**
     Piecewise (face-by-face) surface integration results in a single (summary) force acting on the object:
     
      \f[ \boldsymbol{\sigma}_f = -p(\mathbf{x}) \, \mathbf{I} + \mu \left( \nabla \mathbf{u} + \nabla \mathbf{u}^\top \right)   \f]
       
         where
         \f[ \mathbf{F} = -\int_S \boldsymbol{\sigma}_f \cdot \mathbf{n} \, dS \f]
         
     so that
     \f[ \mathbf{t}(\mathbf{x}) = \boldsymbol{\sigma}_f \cdot \mathbf{n} = -p \, \mathbf{n} + \mu (\nabla \mathbf{u} + \nabla \mathbf{u}^\top) \cdot \mathbf{n} \f]
*/
array<double,3u> pressureForceActingOnObject( const Model<3u>& model, const string& model_subdomain )
 {
    ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( !model.ContainsBoundary(model_subdomain) ) {
         csmp_error.Note( ERROR, "pressureForceActingOnObject", "computation cannot be carried out without further pre-processing");
         return array<double,3>();
      }
      
    // Computing Face(surface)-normal pressures and summing them into a force vector
    const csmp::Index  p_key = model.Database().StorageKey("fluid pressure");
    array<double,3>    sum_of_forces{ {0., 0., 0.} };
    double             obj_surface_area{0.};
    const Boundary<3>& surface = model.Boundary(model_subdomain);
    
    for ( const auto& fit : surface.CellVector() ) {
          // face area
          auto face_area    = fit->Area();
          obj_surface_area += face_area;
          // face unit normal
          Point<3> unrml    = fit->UnitNormal();
          // interpolating nodal forces to cell barycentres
          auto pf = fit->PropertyValueAtBaryCenter( p_key );
          // pressure force integrated over face (pf acts opposite unit normal direction)
          // (computing it at the barycentre rather than integrating it integration point by integration point)
          unrml *= (face_area * pf); // no minus pf, because no nrml projection is performed
          for ( uint32_t i{0u}; i<3u; ++i ) sum_of_forces[i] += unrml[i];
      }
      
    cout <<"\n"<<"pressureForceActingOnObject: integrated pressure forces over surface area of '";
    cout << model_subdomain <<"'' ("<< obj_surface_area <<")"<< endl;
    cout <<"\t"<< sum_of_forces[0] <<","<< sum_of_forces[1] <<","<< sum_of_forces[2] << endl;
    
    return sum_of_forces;
 
 } // end pressureForceActingOnObject








/**
  Computation of combined tractions acting on the surface of the object as well as the resulting torque acting on it.
  
     \f[ \mathbf{t}_q = -p_q \, \mathbf{n}_q + \mu \left( \nabla \mathbf{u}_q + \nabla \mathbf{u}_q^\top \right) \cdot \mathbf{n}_q \f]
 
   The two integrals (pressure-term and viscous-term) are physically distinct and additive.
   
      \f[ \mathbf{F}=\int_S \mathbf{t}\,dS = -\int_S p\,\mathbf{n}\,dS \;+\; \int_S (\boldsymbol{\tau}_{\text{visc}}\cdot\mathbf{n})\,dS \f],
      \f[ \boldsymbol{\Tau}=\int_S (\mathbf{x}-\mathbf{x}_0)\times\mathbf{t}\,dS.\f]

   Computing both and summing them is the correct way to get the full fluid force.
   The torque is just the moment of the same traction distribution.
   For the computation of the torque, the barycentre of the region must be computed.
   
    @note all CSMP variables are converted to Eigen Vector and Matrix types to make use of their special functionality.
 */
template<uint32_t dim>
ForceAndTorque<dim> computeNormalAndShearTractionsOnBoundary( const Model<dim>& model, const string& region_boundary, double mu )
  {
    ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( !model.ContainsBoundary(region_boundary) ) {
         csmp_error.Note( ERROR, "computeNormalAndShearTractionsOnBoundary",
                                 "computation cannot be carried out without further pre-processing");
         return ForceAndTorque<dim>();
      }
    constexpr bool verbose{true};

    const Boundary<dim>& boundary = model.Boundary(region_boundary);

    Eigen::Matrix<double,dim,1> subdomain_centroid = eigen_glue::ToEigen(boundary.Centroid());
    Eigen::Matrix<double,dim,1> total_force  = Eigen::Matrix<double,dim,1>::Zero();
    Eigen::Matrix<double,dim,1> total_torque = Eigen::Matrix<double,dim,1>::Zero();
    
    const csmp::Index   p_key = model.Database().StorageKey("fluid pressure");
    const csmp::Index   u_key = model.Database().StorageKey("pore velocity");
    const csmp::Index   U_key = model.Database().StorageKey("pore velocity face");
    
    VectorVariable<dim> face_velo;
    
    for ( const auto& face : boundary.CellVector() ) {
          auto bctr = eigen_glue::ToEigen(face->BaryCenter());
          auto nrml = eigen_glue::ToEigen(face->UnitNormal());
          for ( uint32_t ip{0u}; ip<face->IntegrationPoints(); ++ip )
            {
                // --- Pressure forces: -pI
                double pf_ip  = face->PropertyValueAtIntegrationPoint( p_key, ip );
                Eigen::Matrix<double,dim,dim> pressure_stress = Eigen::Matrix<double,dim,dim>::Identity() * pf_ip;
                // pressure_traction = pressure_stress * nrml;
                if ( verbose ) cout <<"\n"<<"face-normal stress at ip "<< ip <<":\n"<< pressure_stress * nrml;
                
                // mapping 'pore velocity' from the Face's inside parent to Face
                face->InnerParent()->Read( u_key, face_velo );
                face->Store( U_key, face_velo );

                // --- Viscous stress: mu * (∇u + (∇u)^T);   ∇u is a dim * dim, since u is a vector
                double detJ; // determinant of Jacobian at ip
                auto grad_u = face->PropertyGradientAtIntegrationPoint( U_key, ip, detJ );  // ∇u
                auto shear_stress = mu * (grad_u + grad_u.transpose()); // symmetric 2nd-order shear-stress tensor
                // shear_traction = shear_stress * nrml;
                if ( verbose ) cout <<"\n"<<"shear stress acting on face at ip "<< ip <<":\n"<< shear_stress * nrml;

                // --- Full traction tensor, tau = -pI + tau_visc
                auto traction_tensor = pressure_stress + shear_stress;
                
                // --- Projected traction, t = sigma . n
                Eigen::Matrix<double,dim,1> tn = traction_tensor * nrml;
                
                // --- Torque, M
                auto r = bctr - subdomain_centroid;
                total_torque += r.cross(shear_stress * nrml);
                
                // Surface integration
                total_force += tn * face->WeightAtIntegrationPoint(ip) * detJ;
            }
      }

    if ( verbose ) {
         cout <<"\n"<<"computeNormalAndShearTractionsOnBoundary: force acting on '"<< region_boundary <<"':\n";
         cout <<"\t"<< total_force;
         cout <<"\n"<<"torque acting on '"<< region_boundary <<"':\n";
         cout <<"\t"<< total_torque;
      }

    return { total_force, total_torque };
}

template ForceAndTorque<3> computeNormalAndShearTractionsOnBoundary<3>( const Model<3>&, const string&, double );
//template ForceAndTorque<2> computeNormalAndShearTractionsOnBoundary<2>( const Model<2>&, const string&, double ); //cross product only in 3D





/**
     For computing the viscous shear stress acting on the boundary, we need the full viscous stress tensor, which depends on the velocity gradient (not just pressure).
     For a Newtonian fluid, the stress tensor is:

     \f[ \sigma = -p \mathbf{I} + 2\mu \mathbf{D} \quad \text{with} \quad \mathbf{D} = \frac{1}{2}\left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right) \f]

where:
	•	\mathbf{u} = velocity vector,
	•	\nabla \mathbf{u} = velocity gradient tensor,
	•	\mu = dynamic viscosity,
	•	\mathbf{D} = rate-of-deformation tensor (symmetric part of velocity gradient),
	•	\mathbf{I} = identity tensor,
	•	p = pressure.

     * @brief Compute the viscous stress tensor at a given integration point.
     *
     * @param cell_ptr pointer to the face we are operating on
     * @param velo_key Index of the velocity field variable (VectorVariable)
     * @param mu    Dynamic viscosity
     * @param ip    Integration point index
     * @return csmp::Tensor<dim, dim> Full viscous stress tensor
 */
template<uint32_t dim>
csmp::TensorVariable<dim> computeViscousStressAtIntegrationPoint( const Face<dim>* cell_ptr, const csmp::Index& velo_key, double mu, uint32_t ip )
 {
    // 1. Compute velocity gradient ∇u at integration point
    double detJ;
    Eigen::Matrix<double,dim,Eigen::Dynamic> grad_u = cell_ptr->PropertyGradientAtIntegrationPoint( velo_key, ip, detJ );

    // 2. Compute symmetric part: D = 0.5 * (grad_u + grad_u^T)
    csmp::TensorVariable<dim> D;
    for (uint32_t i = 0; i < dim; ++i) {
        for (uint32_t j = 0; j < dim; ++j) {
            D(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
        }
    }

    // 3. Build viscous stress: σ = 2μD  (pressure handled separately)
    csmp::TensorVariable<dim> viscous_stress;
    for (uint32_t i = 0; i < dim; ++i) {
        for (uint32_t j = 0; j < dim; ++j) {
            viscous_stress(i, j) = 2.0 * mu * D(i, j);
        }
    }

    return viscous_stress;
}

template TensorVariable<3> computeViscousStressAtIntegrationPoint( const Face<3>*, const csmp::Index&, double, uint32_t );

} // csmp
