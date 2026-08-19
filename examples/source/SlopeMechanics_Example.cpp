//
//  SlopeMechanics_Example
//  CSMP_API_library
//
//  Created by Stephan Matthai on 1/27/17.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "SlopeMechanics_Example.h"

#include "compareFloats.h"
#include "ANSYS_Model2D.h"
#include "MeshDiagnostics.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "PDE_Integrator.h"
#include "Face.h"
#include "NumIntegral_NT_op_N_dS.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "LinearSolver.h"
#endif

#include "PropertyHandle.h"
#include "InputDataManager.h"
#include "VTU_Interface.h"

// include any header files that you need here...
#include "SteadyStateDiffusor.h"
#include "CSMP_highLevelUtilities.h"
#include "ComputationalSettings.h"
#include "IAPWS_H2OPropertiesVisitor.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_PT_op_dS.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_BT_D_op_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "StressesAndStrains2.h"
#include "ExtractTensorVariableComponent.h"

#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

void SlopeMechanics_Example::Specifications()
  {
     SetTitle( "Pore pressure and effective stress distribution in a hill slope" );
     SetDifficulty( 2 );
     SetCategory( "Simulation of Physical Processes" );
     AddAuthor( "Stephan Matthai" );
     AddDescription( "source in: SlopeMechanics_Example.cpp" );
     AddDescription( "Empty example for the user to experiment with" );
     AddRequirement( "Input file suite: Jura-slope1, CSMP_field_scale_mechanics_variables.txt" );
  }





/**
     Put CSMP code that you would like to test here and run it as part of the
     example suite.
 
     Systems Modelling and Design model that calculates the effective stress in a dam,
     using the gravitational loading and the plane stress assumption.
 
     input models:
      
     Exercise: extend this example:
                
     TODO: define interrelation between water-content and volume of clays
     TODO: add computation of seepage force
     TODO: consider potential stress changes due to insolation of slope etc.
     TODO: using the stress invariants apply failure criteria: tensile, shear...; assess slope stability iteratively
     TODO: run same model using QFEM = quadratic finite element approximation
     TODO: work with an anisotropic permeability
     TODO: change placement of the strength related variables to the element integration points
     
*/
void SlopeMechanics_Example::Run()
{
  /*
  // ----------------------------------------------------------
  // 1. building and configering model from ANSYS - csp dataset
  // ----------------------------------------------------------
   string input_file("Jura-slope1");
   ANSYS_Model2D   model( input_file.c_str(), "CSMP_field_scale_mechanics_variables.txt" );
   */

   // ------------------------------------------------------------
   // 1. Load CSMP native format model and checking it
   // ------------------------------------------------------------
   string model_name;
   cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'Jura-slope1':"<<endl;
   cin.ignore();
   getline(cin, model_name);
   if (model_name.length() == 0) model_name = "Jura-slope1";

   //find the name of current example source file
   string file_name = GetExampleFileName(__FILE__);
   string variable_file = "CSMP_field_scale_mechanics_variables.txt";
   string config_file = model_name;
   //create of directory with current example name, go into this directory, and copy input files into it.
   CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file, config_file);
   //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
   Model<DIM>  model(model_name, variable_file);


   Region<DIM>& model_domain(model.Region("Model"));

   printModelDimensions( model, true );
   
   MeshDiagnostics<DIM> mesh_quality;
   mesh_quality.FixFiniteElementNeighborOrientationOfSurfaceMeshes( model );
   mesh_quality.ScrutinizeMesh( model, 10., 1.0e-7 );

   InputDataManager<DIM>  model_configuration;
   ComputationalSettings settings;

   model_configuration.ConfigureFromFile( model, config_file.c_str(),
                                          false,           ///< region name from parameter range
                                          true,            ///< default property values
                                          true,            ///< regional property values
                                          true,            ///< boundary conditions for box-shaped model
                                          true,           ///< essential conditions for regions
                                          true,           ///< boundary conditions for arbitrary-shaped model
                                          settings );

  
  // ----------------------------------------------------------
  // 2. computing initial steady state temperature distribution
  // ----------------------------------------------------------
    {
      SteadyStateDiffusor<DIM,Element>  temperature( model, "thermal conductivity", "temperature", "energy source");
      printRangeOfVariable( model, "thermal conductivity" );
      printRangeOfVariable( model, "energy source");
      model.InputBoundaryValue( TOP, "temperature", makeScalar(DIRICH,15.));
      //model.InputBoundaryValue( BOTTOM, "temperature", makeScalar(DIRICH,80.));
      temperature.ComputeSteadyState( model.Region("Model") );
      printRangeOfVariable( model, "temperature");
    }
  
  // Output to VTU (XML encoded VTK files in text format without compression)
  VTU_Interface<DIM> vtu_output( model );
  
  vtu_output.OutputDataToVTU( string(model.Name()) +"_temperature", "temperature", string("Model"), static_cast<long>(0) );

  
  // -------------------------------------------------------------------------
  // 3. computing initial hydrostatic pressure (no loading via grain skeleton)
  // -------------------------------------------------------------------------
  //  Visitor computes thermodynamic properties of water; density is interpolated to
  //     element barycentre for later vertical integration
  // --------------------------------------------------------------------------------------
  IAPWS_H2OPropertiesVisitor<DIM>  properties_visitor( model, "fluid pressure" ,"fluid density", "fluid viscosity");
  // computing fluid density and viscosity
  model.Accept( properties_visitor );
  printRangeOfVariable( model, "fluid density");
  printRangeOfVariable( model, "fluid viscosity");

  // hydraulic conductivity K = k/mu
  const csmp::Index perm_key(model.Database().StorageKey("permeability"));
  const csmp::Index visc_key(model.Database().StorageKey("fluid viscosity"));
  const csmp::Index cond_key(model.Database().StorageKey("conductivity"));
  ScalarVariable  visc;
  for ( auto& it : model_domain.CellVector() ) {
       it->PropertyValueAtBaryCenter( visc_key, visc );
       it->Store( cond_key, makeScalar( it->Status(cond_key), it->Read(perm_key) / visc() ) );
    }
  printRangeOfVariable( model, "conductivity" );
// testing
vtu_output.OutputDataToVTU( string(model.Name()) + "_initial_conductivity", "conductivity", string("Model"), static_cast<long>(0) );
  
  // Setting up the FE algorithm to compute the initial hydrostatic fluid pressure and velocities
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings  samg_settings;
  SAMG_Solver    solver(&samg_settings);
  // minimizing screen output
  samg_settings.Set_iout1( 0 );
  samg_settings.Set_iout2( 0 );
  samg_settings.Set_idmp( -1 );
  samg_settings.SetSolverInstance(2);
#else
  CSMP_DEFAULT_LINEAR_SOLVER  solver;
#endif

  PDE_Integrator<DIM,Element>  hydrostatic_pressure( solver );

  NumIntegral_dNT_lhsop_dN_dV<DIM>  hydrostatic_conductance( model.Database(),
                                                            "conductivity",
                                                            "fluid pressure", "fluid pressure" );

  NumIntegral_dNT_op_dV<DIM>  fluid_gravity( model.Database(), "fluid body force", "fluid pressure" );
                                                                 
  hydrostatic_pressure.Add( &hydrostatic_conductance );
  hydrostatic_pressure.Add( &fluid_gravity );

  // total dissolved solids (ppm = g/tonne; normal seawater=12000 g/t are just added to the fluid density
  double total_dissolved_solids(0.); // g->kg (157500-ppm = 157kg salt)
  total_dissolved_solids /= 1000.; // gets kg/m3

  // Apply Dirichlet (essential) boundary conditions to all sides of the model
  // (irregular and left boundaries where already handled in config file)
  constexpr double p_atm = 100325.;
  /* auto bottom_pf = */ computeAndFreezeBoundaryPressure<DIM>( model, "RIGHT", p_atm, csmp::ACC_GRAVITY );
  // creates issue with pressure constraint on the left
  // model.Boundary("BOTTOM").InputPropertyValue("fluid pressure", makeScalar(DIRICH,bottom_pf) );

  cout <<"\nrun: initial guess of fluid pressure range.";
  printRangeOfVariable( model, "fluid pressure" );

  cout << "\n\n\nmain: Iterating fluid pressure to find correct fluid density and viscosity... " << endl;
  PropertyHandle<DIM>  fluid_bf( model, "fluid body force", VECTOR, ELEMENT );
  const csmp::Index    rhof_key = model.Database().StorageKey("fluid density");
  const csmp::Index    K_key    = model.Database().StorageKey("conductivity");
  const csmp::Index    pf_key(model.Database().StorageKey("fluid pressure"));
  
  for ( uint32_t i=0; i<3U; i++ )
    {
       // computing the element vector property "fluid body force": rho_f(node)_interpolated * g
       for ( auto& it : model_domain.CellVector() ) {
            assert( it->IsLine() == false );
            VectorVariable<DIM> vc( ANY, ANY, 0., -1. * csmp::ACC_GRAVITY * it->Read(K_key) * (it->PropertyValueAtBaryCenter(rhof_key) + total_dissolved_solids) );
            it->Store( fluid_bf.Key(), vc );
         }
       printRangeOfVariable( model, "fluid body force" );
       cout <<"\n\t"<<"Fluid pressure iteration "<< i+1U <<":"<< endl;
       model.Apply( hydrostatic_pressure );
       // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       hydrostatic_pressure.Reset(false);
       auto pmin = printRangeOfVariable( model, "fluid pressure", false );
       // bracketing negative fluid pressures that might arise near the atmospheric pressure boundary due to vertical integration inaccuracy
       // (circumventing this problem would require working with relative fluid density)
       if ( pmin < p_atm ) {
          for ( auto& nit : model_domain.NodeVector() ) {
               double pf = nit->Read( pf_key );
               if ( pf < p_atm ) nit->Store( pf_key, makeScalar(nit->Status(pf_key),p_atm) );
            }
       }
       // updating hydraulic conductivity
       model.Accept( properties_visitor );
       for ( auto& it : model_domain.CellVector() ) {
            it->PropertyValueAtBaryCenter( visc_key, visc );
            it->Store( cond_key, makeScalar( it->Status(cond_key), it->Read(perm_key) / visc() ) );
         }
       printRangeOfVariable( model, "fluid viscosity" );
       printRangeOfVariable( model, "fluid density" );
    }

  vtu_output.OutputDataToVTU( string(model.Name()) + "_iterated-fluid-pressure", "fluid pressure", string("Model"), static_cast<long>(0) );

  // -----------------------------------------------------------
  // 4. computing the Darcy velocity
  // -----------------------------------------------------------
  assignNodeCoordinatesTo( model, 'Y', "elevation" );
  vtu_output.OutputDataToVTU( string(model.Name()) + "_elevation", "elevation", string("Model"), static_cast<long>(0) );
  
  // 4.1 computing the hydraulic head, h = pf/rho*g + z   from the fluid pressure distribution
  const csmp::Index E_key(model.Database().StorageKey("elevation"));
  const csmp::Index h_key(model.Database().StorageKey("hydraulic head"));
  const csmp::Index vD_key(model.Database().StorageKey("Darcy velocity"));

  for ( auto& nit : model_domain.NodeVector() )
    {
        double head = nit->Read(pf_key) / (nit->Read(rhof_key) * csmp::ACC_GRAVITY) + nit->Read(E_key);
       nit->Store( h_key, makeScalar( nit->Status(h_key), head ) );
    }
  printRangeOfVariable( model, "hydraulic head" );
  vtu_output.OutputDataToVTU( string(model.Name()) + "_hydraulic-head", "hydraulic head", string("Model"), static_cast<long>(0) );

  // 4.2 computing hydraulic head, its gradients and the ensuing Darcy velocity
  DenseMatrix<DM_MIN> DERIV;
  VectorVariable<DIM> velo;
  // 4.3 computing ensuing fluid flow
  for ( auto& it : model_domain.CellVector() )
    {
       const double K = it->Read( K_key );
       velo = 0.;
       it->dN_AtBaryCenter( DERIV, 1U );
       for ( uint32_t i=0U; i<it->Nodes(); i++ ) {
            double h = it->N(i)->Read( h_key );
            velo(0)  += h  * -DERIV(0,i) * K;
            velo(1)  += h  * -DERIV(1,i) * K;
            velo.Flag(0) = it->Status(vD_key,0);
            velo.Flag(1) = it->Status(vD_key,1);
         }
      
       it->Store( vD_key, velo );
    }
  printRangeOfVariable( model, "Darcy velocity" );
  vtu_output.OutputDataToVTU( string(model.Name()) + "final-", "Darcy velocity", string("Model"), static_cast<long>(0) );
  model.ExtrapolateCellToNodeProperty( "Darcy velocity", "Darcy velocity node" );
  printRangeOfVariable( model, "Darcy velocity node" );
  vtu_output.OutputDataToVTU( string(model.Name()) + "_Darcy-velocity", "Darcy velocity", string("Model"), static_cast<long>(0) );



  // -----------------------------------------------------------
  // 5. initial rock effective stress due to gravitional loading
  // -----------------------------------------------------------
  // mechanical properties are placed on the element integration points
  // computes the increase in pore-pressure due to the gravitational loading
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  samg_mechanics_settings;
    samg_mechanics_settings.Set_napproach(2); // sorts rhs vector [x1, y1, x2, y2, ..., xn, yn]
    samg_mechanics_settings.SetSolverInstance(2);
                               // which is needed for deformation simulations
    SAMG_Solver                 solver2(&samg_mechanics_settings);
#else
    CSMP_DEFAULT_LINEAR_SOLVER  solver2;
#endif
    PDE_Integrator<DIM,Element>  deformation(solver2);

    // this will also include boundary stresses translated into nodal forces
    const bool principal_vectors(true), plane_strain(false);
    PT_op<DIM>                  bforces( model.Database(), "force", "displacement" );
    NumIntegral_BT_D_B_dV<DIM>  stiffness( model.Database(),
                                          "Youngs modulus", "Poissons ratio", "displacement", "displacement", plane_strain);
    NumIntegral_PT_op_dV<DIM>   bodyforce( model.Database(), "gravity force", "displacement");
    NumIntegral_BT_D_op_dV<DIM> volstrain( model.Database(),
                                          "fluid volume source", "Youngs modulus", "Poissons ratio", "displacement");
    NumIntegral_BT_op_dV<DIM>   porepressure( model.Database(), "Biot term", "displacement");

    // computes the 'Biot term'
    porePressureBiotAlphaProduct( model, "Model" );
    cout <<"\nmain: taking into account (in the displacement equation), the effects of pore pressure.";
    printRangeOfVariable( model, "Biot alpha" );
    printRangeOfVariable( model, "Biot term" );
  
    // computing the gravity force
    gravityForce( model, csmp::ACC_GRAVITY );
    printRangeOfVariable( model, "gravity force" );
  
    deformation.Add( &porepressure );
    deformation.Add( &stiffness );
    deformation.Add( &bforces );  // force vector must always be there so that Dirichlet conditions are accumulated
    deformation.Add( &bodyforce );
    deformation.Add( &volstrain );
    deformation.Add( &porepressure );

    //StressesAndStrains<DIM>  postpro( model, "Youngs modulus", "Poissons ratio", "displacement", principal_vectors );
    StressesAndStrains<DIM>  postpro( model, "Youngs modulus", "Poissons ratio", "displacement", plane_strain, principal_vectors );
    deformation.AddPostProcess( &postpro );
 
    // -------------------------------------------
    // displacement, strain and stress computation
    // -------------------------------------------
    printRangeOfVariable( model, "Youngs modulus" );
    printRangeOfVariable( model, "Poissons ratio" );
    printRangeOfVariable( model, "gravity force" );
    deformation.IntegrateOver( model, model_domain );
    // among other things the stresses and strains operator computes the (FE-based) (elastic) dilatation
    printRangeOfVariable( model, "displacement" );
  
  
   // ---------------------------------------------------------------------------------------
   // extracting the vertical stress component
   // ---------------------------------------------------------------------------------------
    if ( DIM == 2 ) {
         ExtractTensorVariableComponent<DIM>  ystress(   model.Database(), "stress", "stress-y",  1,1 );
         ExtractTensorVariableComponent<DIM>  stress_xy( model.Database(), "stress", "stress-xy", 0,1 );
         model.Apply( ystress );
         model.Apply( stress_xy );
      }
  
   // get a FEM-based estimate of the effects of a loading induced divergence of velocity field
    strainInducedChangeInPorePressure<DIM>( model );
    printRangeOfVariable( model, "fluid pressure" );

    std::list<string> output_props = {"displacement","fluid pressure","strain","stress","mean stress","stress-y","stress-xy"};

    vtu_output.OutputDataToVTU( (string(model.Name()) + "_mechanics-results").c_str(), output_props, string("Model"), static_cast<long>(0) );
  
    cout <<"\nSlopeMechanics_Example: That's it!\n";

   filesystem::current_path("../../example_inputs/");
  
} // end Run




// OTHER FUNCTION RELATIONSHIPS



/** 
    Computation of gravity-equilibrated fixed pressure at vertical boundaries of 2D models,
    by top down piecewise integration. The new pressure values are flagged Dirichlet.
*/
/**
    @brief Computes gravity-equilibrated fluid pressure at a named vertical
    boundary of a 2D model and flags all boundary nodes as Dirichlet.

    Integration proceeds top-down using the trapezoidal rule, approximating
    fluid density as the piecewise average between adjacent boundary nodes.

    @tparam  dim         Spatial dimension of the model (typically 2).
    @param  model       The CSMP model containing the boundary and database.
    @param  edge_boundary Name of the boundary on which pressure is set.
    @param  pressure_at_top Atmospheric (or other) pressure at the topmost node [Pa].
    @param  acc_gravity Gravitational acceleration, positive magnitude [m/s^2].
    @return the maximum pressure found by the top-down integration

    @pre  "fluid density"  must be defined and initialised on boundary nodes.
    @pre  "fluid pressure" must be defined on boundary nodes.
    @pre  acc_gravity > 0.
    
    @post All nodes on edge_boundary have fluid pressure stored with DIRICH flag.
*/
template<uint32_t dim>
double computeAndFreezeBoundaryPressure( Model<dim>&  model,
                                         const char*  edge_boundary,
                                         double       pressure_at_top,
                                         double       acc_gravity )
{
    assert( acc_gravity > 0. );

    // --- validate inputs ---
    if ( !model.Database().IsDefined( "fluid density" ) )
        throw csmp::Exception( ERROR, "computeAndFreezeBoundaryPressure",
                               "fluid density is not defined in the database" );
    if ( !model.Database().IsDefined( "fluid pressure" ) )
        throw csmp::Exception( ERROR, "computeAndFreezeBoundaryPressure",
                               "fluid pressure is not defined in the database" );

    const Boundary<dim>& bref( model.Boundary( edge_boundary ) );

    const csmp::Index rho_key( model.Database().StorageKey( "fluid density"  ) );
    const csmp::Index pf_key ( model.Database().StorageKey( "fluid pressure" ) );

    // --- order boundary nodes top-down by Y coordinate ---
    // multimap handles the case where two nodes share the same Y value
    multimap<double, Node<dim>*, greater<double>> nodes_map;
    for ( auto nit = bref.NodesBegin(); nit != bref.NodesEnd(); ++nit )
        nodes_map.insert( make_pair( (*nit)->y(), *nit ) );

    // --- top-down trapezoidal integration ---
    set<double> pressure_range;
    ScalarVariable sc;

    auto nit = nodes_map.begin();

    // store atmospheric pressure at the topmost node
    (*nit).second->Store( pf_key, makeScalar( DIRICH, pressure_at_top ) );
    pressure_range.insert( pressure_at_top );

    (*nit).second->Read( rho_key, sc );
    double density1  = sc();
    double pressure1 = pressure_at_top;
    double elevation = (*nit).first;
    ++nit;

    while ( nit != nodes_map.end() ) {
        (*nit).second->Read( rho_key, sc );
        const double density2  = sc();
        const double rho_avg   = ( density1 + density2 ) / 2.;
        const double dz        = elevation - (*nit).first; // positive going downward
        const double pressure2 = pressure1 + acc_gravity * rho_avg * dz;

        (*nit).second->Store( pf_key, makeScalar( DIRICH, pressure2 ) );
        pressure_range.insert( pressure2 );

        pressure1 = pressure2;
        density1  = density2;
        elevation = (*nit).first;
        ++nit;
      }

    cout << "\ncomputeAndFreezeBoundaryPressure: pressure range at boundary '"
         << edge_boundary << "': "
         << *pressure_range.begin()  << " to "
         << *pressure_range.rbegin() << " Pa.\n";
         
   return *pressure_range.rbegin();

} // computeAndFreezeBoundaryPressure

template double computeAndFreezeBoundaryPressure<DIM>( Model<DIM>&, const char*, double, double );





/**
    @brief Computes the instantaneous change in pore pressure induced by
    volumetric strain at element integration points.

    Rather than deforming the mesh or using the finite volume framework,
    this function evaluates the volumetric strain directly at the integration
    points of each element — where it is exact for isoparametric elements —
    and distributes the resulting pore pressure change to the element nodes
    via shape function weighting.

    The poromechanical coupling is:

    @f[
        \Delta p_f = -\frac{\epsilon_v}{c_f}
    @f]

    where @f$ \epsilon_v = \nabla \cdot \mathbf{u} @f$ is the volumetric
    strain (trace of the strain tensor) and @f$ c_f @f$ is the fluid
    compressibility. Compression (@f$ \epsilon_v < 0 @f$) increases pore
    pressure (@f$ \Delta p_f > 0 @f$).

    The nodal pressure update is computed as the integration-point-weighted
    average over all elements sharing the node:

    @f[
        \Delta p_f^{node} = -\frac{1}{c_f} \frac{\sum_e \sum_i w_i |J_i| \epsilon_v^{e,i} N_j^{e,i}}
                                                  {\sum_e \sum_i w_i |J_i| N_j^{e,i}}
    @f]

    where the outer sum is over all elements @f$ e @f$ sharing the node,
    the inner sum is over integration points @f$ i @f$, @f$ w_i @f$ are
    the quadrature weights, @f$ |J_i| @f$ is the Jacobian determinant,
    and @f$ N_j^{e,i} @f$ is the shape function value at integration
    point @f$ i @f$ for the local node corresponding to the global node.

    @note The mesh is not deformed. Displacements are read from the
    displacement variable but node coordinates are not updated.

    @note Boundary nodes retain their Dirichlet status. If the boundary
    pressure should respond to the strain field, remove the Dirichlet
    guard before calling this function.

    @param model  The model on which to compute the pressure update.

    @throws csmp::Exception if fluid compressibility is not strictly
            positive at any node, or if the strain tensor variable is
            not defined at element integration points.
*/
template<uint32_t dim>
void strainInducedChangeInPorePressure( Model<dim>& model )
{
    const csmp::Index pf_key  ( model.Database().StorageKey( "fluid pressure"       ) );
    const csmp::Index bf_key  ( model.Database().StorageKey( "fluid compressibility") );
    const csmp::Index eps_key ( model.Database().StorageKey( "strain"               ) );

    // confirm strain is stored at integration points — this is where it is exact
    if ( model.Database().Placement("strain") != ELEMENT_INTEGRATION_POINT )
        throw csmp::Exception( ERROR,
            "strainInducedChangeInPorePressure",
            "strain must be stored at element integration points" );

    Region<dim>& model_domain( model.Region( "Model" ) );

    // --- accumulate weighted dilatation and weight at each node ---
    // Two passes: first accumulate, then normalise and update pressure.
    // This avoids double-counting when a node is shared between elements.

    // node index -> accumulated weighted dilatation
    std::map<size_t, double> weighted_dilatation;
    // node index -> accumulated weight
    std::map<size_t, double> accumulated_weight;

    for ( auto eit  = model_domain.CellsBegin();
               eit != model_domain.CellsEnd(); ++eit ) {
        const Element<dim>* e = *eit;

        if ( !e->UsesLocalCoordinates() ) continue;

        const uint32_t n_integration_points = e->IntegrationPoints();
        const uint32_t n_nodes              = e->Nodes();

        // shape function values at each integration point — n_nodes x n_ip
        // dN_AtIntegrationPoint fills the gradient matrix and returns detJ;
        // we need N (not dN) here, so we use N_AtIntegrationPoint
        for ( uint32_t i{ 0U }; i < n_integration_points; ++i ) {

            // shape function values at integration point i
            std::vector<double> N( n_nodes );
            e->N_AtIntegrationPoint( i, N );

            const double weight = e->WeightAtIntegrationPoint( i );

            // Jacobian determinant at integration point i
            DenseMatrix<DM_MIN> gradientMatrix;
            const double detJ = e->dN_AtIntegrationPoint( gradientMatrix, i, SCALAR );

            const double quadrature_weight = weight * detJ;

            // read strain tensor at this integration point
            TensorVariable<dim> strain;
            e->Read( i, eps_key, strain );

            // volumetric strain = trace of strain tensor
            double volumetric_strain = 0.;
            for ( uint32_t d{ 0U }; d < dim; ++d )
                volumetric_strain += strain( d, d );

            // distribute to nodes via shape function weighting
            for ( uint32_t j{ 0U }; j < n_nodes; ++j ) {
                const size_t global_node_idx = e->N(j)->Idx();
                const double nodal_weight    = N[j] * quadrature_weight;

                weighted_dilatation[ global_node_idx ] += volumetric_strain * nodal_weight;
                accumulated_weight [ global_node_idx ] += nodal_weight;
            }
        }
    }

    // --- second pass: normalise and update nodal pore pressure ---
    for ( auto& nit : model_domain.NodeVector() ) {
        const size_t node_idx = nit->Idx();

        // skip nodes with no accumulated weight (e.g. isolated nodes)
        const auto weight_it = accumulated_weight.find( node_idx );
        if ( weight_it == accumulated_weight.end() ) continue;
        if ( weight_it->second <= 0. )               continue;

        // read fluid compressibility at this node
        ScalarVariable bf_var;
        nit->Read( bf_key, bf_var );
        const double fluid_compressibility = bf_var();

        if ( fluid_compressibility <= 0. )
            throw csmp::Exception( ERROR,
                "strainInducedChangeInPorePressure",
                "fluid compressibility must be strictly positive" );

        // normalised volumetric strain at this node
        const double nodal_dilatation = weighted_dilatation.at( node_idx )
                                      / weight_it->second;

        // compression (negative dilatation) increases pore pressure
        const double delta_pf = -nodal_dilatation / fluid_compressibility;

        // read current pressure and update, preserving Dirichlet status
        ScalarVariable pf_var;
        nit->Read( pf_key, pf_var );

        // do not overwrite Dirichlet boundary conditions
        if ( pf_var.Flag() == DIRICH ) continue;

        nit->Store( pf_key, makeScalar( ANY, pf_var() + delta_pf ) );
    }

} // end strainInducedChangeInPorePressure

// explicit instantiations
template void strainInducedChangeInPorePressure( Model<2U>& );
template void strainInducedChangeInPorePressure( Model<3U>& );




void permeabilityPorosityCorrelation( Model<DIM>& model )
 {
    const csmp::Index phi_key(model.Database().StorageKey("porosity"));
    const csmp::Index k_key(model.Database().StorageKey("permeability"));
 
    Region<DIM>& model_domain(model.Region("Model"));
   
    for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it )
      {
         double porosity = (*it)->Read( phi_key );
         // TODO: exercise: introduce credible relationship here
         double permeability = porosity * 1.0e-12;
         (*it)->Store( k_key, makeScalar((*it)->Status(k_key), permeability) );
      }
 
 } // end permeabilityPorosityCorrelation





/**
    Computes the product of the Biot coefficient alpha with the pore pressure 
    for assignment to the pre-stress term in the governing equation
    
    The result is returned into the scalar element variable 'Biot term'
    
    @attention the fluid pressure effect is averaged on the element barycenter
*/
void porePressureBiotAlphaProduct( Model<DIM>& model, const char* target_region )
 {
    const csmp::Index pf_key    = model.Database().StorageKey("fluid pressure");
    const csmp::Index alpha_key = model.Database().StorageKey("Biot alpha");
    assert( alpha_key.place == ELEMENT );
    assert( alpha_key.type  == SCALAR );
    // result
    const csmp::Index res_key = model.Database().StorageKey("Biot term");
    assert( res_key.place == ELEMENT );
    assert( res_key.type  == SCALAR );

    csmp::Region<DIM>& ref = model.Region(target_region);
    ScalarVariable sc;
    
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); it++ )
      {
         const double alpha((*it)->Read( alpha_key ));
         if ( alpha < 0. or alpha > 1. ) {
              cerr <<"\n\tBiot coefficient alpha (1 - K_dry/K_grain): "<< alpha;
              throw csmp::Exception( ERROR, "porePressureBiotAlphaProduct:", "Biot coefficient alpha is out of range." );
           }
         (*it)->PropertyValueAtBaryCenter( pf_key, sc );
         sc *= alpha;
         (*it)->Store( res_key, makeScalar( (*it)->Status(res_key), sc() ) );
      }
      
 } // end porePressureBiotAlphaProduct





/**
    Computing the specific gravity force that acts on the rock skeleton from the "dry rock density"
*/
void gravityForce( Model<DIM>& model, double acc_gravity )
 {
    csmp::Region<DIM>& ref = model.Region("Model");

    const csmp::Index drd_key = model.Database().StorageKey("dry rock density");
    const csmp::Index gf_key  = model.Database().StorageKey("gravity force");
   
    VectorVariable<DIM> gforce(ANY,0. );
    model.InputPropertyValue( "gravity force", gforce );
   
    for ( auto it=ref.CellsBegin(); it!=ref.CellsEnd(); it++ )
      {
         (*it)->Read( gf_key, gforce );
         gforce(1) += (*it)->Read( drd_key ) * acc_gravity;
         (*it)->Store( gf_key, gforce );
      }
   
 } // end gravity force



} // csmp
