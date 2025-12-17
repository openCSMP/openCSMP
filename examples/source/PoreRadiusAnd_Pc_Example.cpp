#include "PoreRadiusAnd_Pc_Example.h"

// the CSMP model
#include "Index.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "geometricCalculations.h"

// the FE algorithm
#include "PDE_Integrator.h"
#include "PropertyHandle.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
#include "Node.h"

// PDE operators building the FE algorithm
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "NumIntegral_SetRHS_to_One.h"

// SAMG settings for solving the FE algorithm
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "LinearSolver.h"
#endif

// FV scheme
#include "NodeCenteredFiniteVolumeTransport.h"

// interfaces
#include "VTK_Interface.h"
#include "MatlabInterface.h"
#include "ExtractVectorVariableLength.h"

// utility functions
#include "CSMP_definitions.h"
#include "CSMP_highLevelUtilities.h"
#include "CSMP_mathUtilities.h"

//visitors
#include "VelocityAndVolumeFlux.h"

//#define DIM 2
#define DIM 3

using namespace std;

namespace csmp {

void PoreRadiusAnd_Pc_Example::Specifications()
{
  SetTitle( "Pore-radius and capillary pressure computation" );
  SetDifficulty( 1 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "SKM" );
  AddDescription( "as presented in Akanjii (2009, TIPM)" );
  AddDescription( "source in: PoreRadiusAnd_Pc_Example.cpp" );
  AddRequirement( "input file: 'tubes_quadra', 3D-model in ANSYS format");
}


/**
    Computes the pore-radii distribution, equivelent permeability, and local capillary entry pressure
    of the supplied sample.

    @attention current code only supports box-shaped models
    
    TODO: use Laplace equation surface-tension simulation approach in the multiphase flow formulation
    TODO: add possibility to enter model name and scale model with command-line arguments
    for this model
*/
void PoreRadiusAnd_Pc_Example::Run()
{
  // -----------------------------------------------------------------------
  // 1. CSMP Model from ANSYS input file set and property values
  // -----------------------------------------------------------------------
//  const string   input_file("constant_width_quadratic_projected");
//  const string   input_file("wedge2dtest_quadratic");
//  const string   input_file("constant_width_quadratic_projected");
//  const string   input_file("bcc_0366_coarse"); OK
//  const string   input_file("bcc_0366_coarse_quadra"); OK
//  const string   input_file("fcc_0366_fine_scaled");
  string input_file{"tubes_quadra_linearized"}; // LFEM model
  double scale_factor{1.};
	cout <<"\nmain: Enter name of pore-scale model and scaling factor: ";
  cin >> input_file >> scale_factor;      // TODO: use the model 'tubes'
  const double image_resolution(1.5e-6);  // microplugs (micrometer resolution)
  
#if DIM == 2
//  ANSYS_Model2D  model( input_file.c_str(), "pore_scale-variables.txt" );
  ANSYS_Model2D  model( input_file.c_str(), "pore_scale_quadratic-variables.txt" );
#elif DIM == 3
  ANSYS_Model3D model( input_file.c_str(), "pore_scale-variables.txt" );
// scale the model size by 1/value
//	if ( fabs(scaling_factor-1.) > numeric_limits<double>::epsilon() ) scaleModel( model, 1./scaling_factor );
#elif
  No model realisation exists for 1D
#endif
  
  const PropertyDatabase<DIM>&  p_ref = model.Database();
  printModelDimensions( model, true );

  const double dynamic_viscosity(1.0e-3);
  model.InputPropertyValue( "viscosity", makeScalar(PLAIN,dynamic_viscosity) );
  model.InputPropertyValue( "parabolic function", makeScalar(ANY,0.) );
  model.InputPropertyValue( "one", makeScalar(PLAIN,1.) );
  model.InputPropertyValue( "conductivity", makeScalar(PLAIN,0.) );
  
  const BOX_BOUNDARY  inflow_boundary(LEFT), outflow_boundary(RIGHT);
 
  VTK_Interface<DIM>  vtk_output;
  
  // -----------------------------------------------------------------------
  // 2. Parabolic profiles in all void spaces
  // -------------------------------------------------------------------------
  const string void_space("PORES");
  const string pore_edges("PORE_EDGE");
  const string solid_skeleton("GRAINS");
  Region<DIM>& pore_space(model.Region(void_space));
  
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings settings;
  SAMG_Solver   solver(&settings);
#else
  CSMP_DEFAULT_LINEAR_SOLVER  solver;
#endif
  PDE_Integrator<DIM,Element>  parabolic_profile(solver);

  NumIntegral_dNT_dN_dV<DIM>   laplacian( p_ref, "parabolic function", "parabolic function" );
  NumIntegral_NT_op_N_dV<DIM>  rhs( p_ref, "one", "parabolic function" );

  parabolic_profile.Add( &laplacian );
  parabolic_profile.Add( &rhs );           // End the numIntegration
  
  // 2.1 no-slip Dirichlet pore-wall conditions: parabolic function at walls = zero
  // -----------------------------------------------------------------------
  pore_space.InputPropertyValue( "parabolic function", makeScalar(DIRICH,0.), PERIMETER );
  
  // 2.2 boundary conditions: Dirichlet flags are released since where the boundary is open
  // -----------------------------------------------------------------------
  // inlet and outlet
  model.InputBoundaryValue( LEFT, "parabolic function",  makeScalar(PLAIN,1.) );
  model.InputBoundaryValue( RIGHT, "parabolic function",  makeScalar(PLAIN,1.) );
  // side boundaries
  model.InputBoundaryValue( TOP, "parabolic function",    makeScalar(PLAIN,1.) );
  model.InputBoundaryValue( BOTTOM, "parabolic function", makeScalar(PLAIN,1.) );
  model.InputBoundaryValue( FRONT, "parabolic function",  makeScalar(PLAIN,1.) );
  model.InputBoundaryValue( BACK, "parabolic function",   makeScalar(PLAIN,1.) );
  // except for the intersection of the pore walls with the in/outlet boundaries
  Region<DIM>& pore_edge(model.Region(pore_edges));
  pore_edge.InputPropertyValue( "parabolic function", makeScalar(DIRICH,0.) );
  // TODO: the nodes on the edges of the pore space are not freed up; fix!

  // 2.3 computation of parabolic function
  // -----------------------------------------------------------------------
  // TODO: create a criterion for the minimum refinement of the channel
  parabolic_profile.IntegrateOver( pore_space );
  parabolic_profile.Reset();
  printRangeOfVariable( model, void_space.c_str(), "parabolic function" );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "parabolic_function", "parabolic function", 1L, true );
  
  // 2.4 computing the pore diameter
  // -----------------------------------------------------------------------
  const PLACEMENT pgrad_funct_placement = model.Database().Placement("magnitude grad parabolic function");
  PropertyHandle<DIM>  grad( model, "grad parabolic function", VECTOR, pgrad_funct_placement );
  PropertyHandle<DIM>  magn( model, "magnitude grad parabolic function", SCALAR, pgrad_funct_placement );
  string  channel_width("channel width");

  // (fracture) cavity radius = 1/2-aperture (=derivative of "parabolic function")
  model.CopyGradientOfProperty_A_To_B( "parabolic function", "grad parabolic function" );
  ExtractVectorVariableLength<DIM>  grad_magnitude( model.Database(), "grad parabolic function", "magnitude grad parabolic function" );
  pore_space.Apply( grad_magnitude );
  printRangeOfVariable( model, void_space.c_str(), "magnitude grad parabolic function" );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "grad_parabolic_function", "grad parabolic function", 1L, true );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "gmagn_parabolic_function", "magnitude grad parabolic function", 1L, true );

  model.InputPropertyValue( channel_width.c_str(), makeScalar(PLAIN,0.) );

  // element to node vs. integration-point to node extrapolation of properties
  if ( pgrad_funct_placement == ELEMENT )
    pore_space.ExtrapolateCellToNodeProperty( "magnitude grad parabolic function", channel_width.c_str() );
  else
    pore_space.ExtrapolateIntegrationPointToNodeProperty( "magnitude grad parabolic function", channel_width.c_str() );

  printRangeOfVariable( model, void_space.c_str(), channel_width.c_str() );
  PropertyHandle<DIM>  cwidth( model, channel_width.c_str() );
  cwidth.OutputCondition(PLAIN);
  cwidth *= 2.;
  printRangeOfVariable( model, void_space.c_str(), channel_width.c_str() );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "channel_width", channel_width.c_str(), 1L, true );


  // 2.5 nodal conductivity values (parabolic function / viscosity) values
  //     for accumulation over integration points
  // -----------------------------------------------------------------------
  // (storing nodal values requires much less space than element values)
  const csmp::Index pfun_key = model.Database().StorageKey("parabolic function");
  const csmp::Index visc_key = model.Database().StorageKey("viscosity");    // node
  const csmp::Index cond_key = model.Database().StorageKey("conductivity"); // node

  for ( auto nit=pore_space.NodesBegin(); nit!=pore_space.NodesEnd(); ++nit )
    {
       double conductivity = (*nit)->Read( pfun_key ) / (*nit)->Read( visc_key );
       if ( conductivity < numeric_limits<double>::epsilon() ) conductivity = 1.0e-20;
       (*nit)->Store( cond_key, makeScalar(PLAIN,conductivity) );
    }
  printRangeOfVariable( model, void_space.c_str(), "conductivity" );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "conductivity", "conductivity", 1L, true );

  
  // -------------------------------------------------
  // 3. Calculating dynamic pore-pressure distribution
  // -------------------------------------------------
  model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) ); // for examination with Paraview
  model.InputBoundaryValue( inflow_boundary, "fluid pressure", makeScalar(DIRICH,1.1e5) );
  model.InputBoundaryValue( outflow_boundary, "fluid pressure", makeScalar(DIRICH,1.0e5) );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "fluid_pressure", "fluid pressure", 0L, true );

  PDE_Integrator<DIM,Element>      fluid_pressure(solver);
  NumIntegral_dNT_op_dN_dV<DIM>    conductance(  p_ref, "conductivity", "fluid pressure", "fluid pressure" );
  NumIntegral_SetRHS_to_Zero<DIM>  rhs0(  p_ref, "fluid pressure" );
  fluid_pressure.Add( &conductance );
  fluid_pressure.Add( &rhs0 );
  fluid_pressure.IntegrateOver( pore_space );
  fluid_pressure.Reset();
  printRangeOfVariable( model, void_space.c_str(), "fluid pressure" );
  
  // 3.1 post-processing of velocities
  // -----------------------------------------------------------------------
  const csmp::Index pf_key = model.Database().StorageKey("fluid pressure");
  const csmp::Index vt_key = model.Database().StorageKey("velocity");
  const csmp::Index qv_key = model.Database().StorageKey("volume flux");
  DenseMatrix<36>      DERIV;
  VectorVariable<DIM>  velo;
  ScalarVariable       conductivity;

  if ( pgrad_funct_placement == ELEMENT ) {
      for ( auto eit=pore_space.CellsBegin(); eit!=pore_space.CellsEnd(); ++eit )
        {
            (*eit)->PropertyValueAtBaryCenter( cond_key, conductivity );
            velo = 0.;
            (*eit)->dN_AtBaryCenter( DERIV, 1U );
            for ( uint32_t i=0U; i<(*eit)->Nodes(); i++ ) {
                 double pf = (*eit)->N(i)->Read( pf_key );
                 velo(0)  += pf  * -DERIV(0,i) * conductivity();
                 velo(1)  += pf  * -DERIV(1,i) * conductivity();
                 if ( DIM == 3 ) velo(2) += pf * -DERIV(2,i) * conductivity();
              }
            (*eit)->Store( vt_key, velo );
            (*eit)->Store( qv_key, makeScalar(PLAIN,velo.Length()) );
        }
    }
  // if the conductivity is placed at the element integration points
  else {
      for ( auto eit=pore_space.CellsBegin(); eit!=pore_space.CellsEnd(); ++eit )
        for ( uint32_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
          {
              (*eit)->PropertyValueAtIntegrationPoint( cond_key, i, conductivity );
              velo = 0.;
              (*eit)->dN_AtIntegrationPoint( DERIV, i, 1U );
              for ( uint32_t j=0U; j<(*eit)->Nodes(); j++ ) {
                   double pf = (*eit)->N(i)->Read( pf_key );
                   velo(0)  += pf  * -DERIV(0,j) * conductivity();
                   velo(1)  += pf  * -DERIV(1,j) * conductivity();
                   if ( DIM == 3 ) velo(2) += pf * -DERIV(2,j) * conductivity();
                }
              // integration point properties
              (*eit)->Store( i, vt_key, velo );
              (*eit)->Store( i, qv_key, makeScalar(PLAIN,velo.Length()) );
          }
    }
  
  printRangeOfVariable( model, void_space.c_str(), "velocity" );
  
  // velocity at the element integration points
/*
  const csmp::Index vp_key = model.Database().StorageKey("point velocity"); // vt at element integration points
  for ( auto eit=pore_space.ElementsBegin(); eit!=pore_space.ElementsEnd(); ++eit )
  for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
    {
        (*eit)->PropertyValueAtIntegrationPoint( cond_key, i, conductivity );
        velo = 0.;
        (*eit)->dN_AtIntegrationPoint( DERIV, i, 1U );
        for ( size_t i=0U; i<(*eit)->Nodes(); i++ ) {
             double pf = (*eit)->N(i)->Read( pf_key );
             velo(0)  += pf  * -DERIV(0,i) * conductivity();
             velo(1)  += pf  * -DERIV(1,i) * conductivity();
             if ( DIM == 3 ) velo(2) += pf * -DERIV(2,i) * conductivity();
          }
        // integration point properties
        (*eit)->Store( i, vp_key, velo );
    }
  printRangeOfVariable( model, void_space.c_str(), "point velocity" );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "ipoint_velocity", "point velocity", 1L, true );
*/
  
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "fluid_pressure", "fluid pressure", 1L, true );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "velocity", "velocity", 1L, true );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "volume_flux", "volume flux", 1L, true );


  // ----------------------------------------------------------------------------
  // 4. Computing the equivalent permeability in the flow direction
  // ----------------------------------------------------------------------------
  double flux_through_model;
  double k_equiv = equivalentPermeability( model, "PORES",
                                             inflow_boundary, outflow_boundary,
                                             flux_through_model );
  
  cout <<"\nmain: The equivalent permeability of the sample is: "<< k_equiv <<" m2.\n";


  // ----------------------------------------------------------------------------
  // 5. Non-wetting phase pressure calculation from pore-radius and wetting angle
  // ----------------------------------------------------------------------------
  /* 
     (non-wetting phase pressure at the pore wall is calculated from IFT, channel width, and wetting angle;
      this value is set to Dirich and then interpolated, using Laplace equation, across
      the width of the channel)
  */
  const double IFT(0.04); // interfacial tension of 40 mN/m converted into N/m (mN/m = dynes/cm)  = energy per unit area
  const double alpha(5.);  // wetting angle (between wetting phase and solid) of 5 degrees (out of 360)

  const csmp::Index nwp_key = model.Database().StorageKey("non-wetting phase entry pressure");
  const csmp::Index cwd_key = model.Database().StorageKey(channel_width.c_str());

  model.InputPropertyValue("non-wetting phase entry pressure", makeScalar(PLAIN,0.) );
  
  // computing and assigning pressure values as Dirichlet boundary conditions at the pore walls
  for ( auto nit=pore_space.PerimeterNodesBegin(); nit!=pore_space.NodesEnd(); ++nit )
    {
       // pnw = (IFT * cos(alpha)) / channel width (=2R, but capped by image resolution)
       double nw_p = (IFT * std::cos(degreesToRadians(alpha))) / std::max((*nit)->Read(cwd_key), image_resolution);
       (*nit)->Store( nwp_key, makeScalar(DIRICH,nw_p) );
    }
  printRangeOfVariable( model, void_space.c_str(), "non-wetting phase entry pressure" );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "nw_phase_pressure", "non-wetting phase entry pressure", 0L, true );

  PDE_Integrator<DIM,Element>      non_wetting_pressure(solver);
  NumIntegral_dNT_dN_dV<DIM>       laplacian2(  p_ref, "non-wetting phase entry pressure", "non-wetting phase entry pressure" );
  NumIntegral_SetRHS_to_Zero<DIM>  rhs2(  p_ref, "non-wetting phase entry pressure" );
  non_wetting_pressure.Add( &laplacian2 );
  non_wetting_pressure.Add( &rhs2 );
  non_wetting_pressure.IntegrateOver( pore_space );
  non_wetting_pressure.Reset();
  printRangeOfVariable( model, void_space.c_str(), "non-wetting phase entry pressure" );
  vtk_output.OutputDataToVTK( model, void_space.c_str(), "nw_phase_pressure", "non-wetting phase entry pressure", 1L, true );

   // terminate
  cout <<"\nmain: That's it..."<< endl;
  
} // end Run()




/**
    scale model with the inverse of the supplied value
*/
void scaleModel( Model<3>& mdl, double scale_factor )
  {
    assert( scale_factor > 0. );
    Region<3>& mref = mdl.Region("Model");
    for ( vector<Node<3>* >::const_iterator
          nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++)
      {
        (*nit)->x( (*nit)->x() / scale_factor);
        (*nit)->y( (*nit)->y() / scale_factor);
        (*nit)->z( (*nit)->z() / scale_factor);
      }
  }





/**
    Computes flux through the boundary of the model using finite-element integration.
    
    A single unit normal is used for the flux calculation. 
    It is taken from the Box class.
*/
template<uint32_t dim>
double boxBoundaryFluxFEM( const Region<dim>& flow_domain, const Index& velo_key, BOX_BOUNDARY target_boundary )
 {
    double             boundary_flux(0.);
    VectorVariable<dim>  velo;
    vector<size_t>       face_nids;
    vector<double>     nrml;
    Box().UnitNormalTo( target_boundary, dim, nrml );
    const Point<dim>     unrml(nrml);
    const bool           debug(false);
    // for integration point property extrapolation to face
    vector<double>     IVAR, NVAR, N;

    // for all boundary faces of the region model
    for ( size_t eidx=flow_domain.InteriorCells(); eidx<flow_domain.Cells(); ++eidx )
      {
         // retrieve the velocity that must already be scaled by thickness if lower dimensional elements are used
         if ( velo_key.place == ELEMENT ) flow_domain.E(eidx)->Read( velo_key, velo );
         // if the velocity is an integration point property
         else {
              assert( velo_key.place == ELEMENT_INTEGRATION_POINT );
              // extrapolating 'velocity' to 'nodal velocity'
              IVAR.resize( flow_domain.E(eidx)->IntegrationPoints() * dim );
              NVAR.resize( flow_domain.E(eidx)->Nodes() * dim );
              size_t  values(0U);
              for ( uint32_t i=0U; i<flow_domain.E(eidx)->IntegrationPoints(); ++i ) {
                    flow_domain.E(eidx)->Read( i, velo_key, velo );
                    for ( uint32_t j=0U; j<dim; ++j ) IVAR[++values] = velo[j];
                }
              flow_domain.E(eidx)->ExtrapolateIntegrationPointVariableToNodes( dim, IVAR, NVAR );
              // averaging nodal velocity
              velo = 0;
              size_t count(0U);
              for ( uint32_t i=0U; i<flow_domain.E(eidx)->Nodes(); ++i )
                for ( uint32_t j=0U; j<dim; ++j ) velo(j) += NVAR[++count];
              velo /= static_cast<double>(flow_domain.E(eidx)->Nodes());
           }

         // loop over the boundary faces of the elements
         for ( uint32_t j=0U; j<flow_domain.PerimeterFaces(eidx); ++j )
           {  // accessing the one or multiple faces of the element that lie on the model boundary
              const uint32_t face = flow_domain.PerimeterFace(eidx,j);
              // SKM TEST: lower dimensional faces
              if ( debug and face_nids.size() == 2 ) cerr <<".";

              for ( auto k : flow_domain.E(eidx)->FE()->NodesOfFace(face) )
                   if ( flow_domain.E(eidx)->N(k)->AtBoundary() == target_boundary ) {
                       // getting the unit normal and area of target face
                       double  Aface = faceArea( (*flow_domain.E(eidx)), face );
                       // computing the face flux Aj v . n
                       boundary_flux += Aface * velo.DotProduct( unrml );

              // SKM TEST: lower dimensional faces
              if ( debug and face_nids.size() == 2 ) cerr <<" A="<< Aface <<": "<< velo.DotProduct( unrml );

                       break;
                    }
           }
      }

    return boundary_flux;

 } // end BoxBoundaryFluxFEM

template double boxBoundaryFluxFEM<2>( const Region<2>&, const Index&, BOX_BOUNDARY );
template double boxBoundaryFluxFEM<3>( const Region<3>&, const Index&, BOX_BOUNDARY );





/**
    Loop over the faces of the boundary, computing normal fluxes 
    from the parent element velocties and integrating them over the model boundary.
    
    @return flux_integral - return the integral scalar flux value.
    
    @test gives big discrepancy between inflow and outflow.
*/
double boundaryFluxFEM( const Model<3U>& model, const char* boundary )
 {
    assert( model.ContainsBoundary(boundary) == true );
    const Boundary<3U>&  boundary_domain(model.Boundary(boundary));
   
    const Index   v_key(model.Database().StorageKey("velocity"));
   
    // looping over the faces of the boundary
    double              flux_integral(0.);
    VectorVariable<3U>  velo;
    for ( auto it=boundary_domain.CellsBegin(); it!=boundary_domain.CellsEnd(); ++it )
      {
         Point<3U>  nrml = (*it)->UnitNormal(); // outward pointing
         // get the velocity from inside parent element
         (*it)->InnerParent()->Read( v_key, velo );
         double face_flux = (nrml[0]*velo[0] +  nrml[1]*velo[1] + nrml[2]*velo[2]) * (*it)->Area();
         // accumulating the face fluxes into the integral
         flux_integral += face_flux;
       }

    return flux_integral;
 }






/**
    Assuming a box-shaped model, function calculates the equivalent permeability parallel to 
    its side boundaries.
*/
template<uint32_t dim>
double equivalentPermeability( Model<dim>& model, const char* region,
                               BOX_BOUNDARY inflow_boundary, BOX_BOUNDARY outflow_boundary,
                               double& flux_through_model )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const Region<dim>&  flow_domain(model.Region(region));
   
    // 1. calculating the cross-sectional area of the model
    // ---------------------------------------------------
    // (2 options: left-right and top bottom)
    Point<dim>  xyz_min, xyz_max;
    model.MinMaxCoordinates( xyz_min, xyz_max );
    // assumption is that flow boundaries are opposite to one another
    double area = xyz_max[0] - xyz_min[0];
    if ( inflow_boundary == LEFT or inflow_boundary == RIGHT )
      area = xyz_max[1] - xyz_min[1];
   
    const double xsect_area(area);
   
   
    // 2. calculating the far-field fluid pressure gradient
    // -----------------------------------------------------
    // O.K.
    const Index pf_key = model.Database().StorageKey("fluid pressure");
    Point<dim> in_coord, out_coord;
    in_coord = out_coord = 0.;
    double inflow_p(0.), outflow_p(0.);
    size_t   in_nodes(0U), out_nodes(0U);
    // for the entire flow domain
    for ( auto nit=flow_domain.PerimeterNodesBegin(); nit!=flow_domain.NodesEnd(); ++nit )
      {
         if ( (*nit)->AtBoundary() == inflow_boundary ) {
              inflow_p += (*nit)->Read( pf_key );
              in_coord += (*nit)->Coordinate();
              in_nodes++;
           }
         else if ( (*nit)->AtBoundary() == outflow_boundary ) {
              outflow_p += (*nit)->Read( pf_key );
              out_coord += (*nit)->Coordinate();
              out_nodes++;
           }
      }
    inflow_p  /= static_cast<double>(in_nodes);
    outflow_p /= static_cast<double>(out_nodes);
    in_coord  /= static_cast<double>(in_nodes);
    out_coord /= static_cast<double>(out_nodes);
   
    double farfield_pf_gradient = (inflow_p - outflow_p);
    if      ( inflow_boundary == LEFT )   farfield_pf_gradient /= fabs( out_coord[0] - in_coord[0] );
    else if ( inflow_boundary == BOTTOM ) farfield_pf_gradient /= fabs( out_coord[1] - in_coord[1] );
    else if ( inflow_boundary == BACK )   farfield_pf_gradient /= fabs( out_coord[2] - in_coord[2] );


	  // 3. fluid throughput through model as determined by velocity projections on model boundaries
    // -------------------------------------------------------------------------------------------
    const Index   vt_key = model.Database().StorageKey("velocity");
    const double  inflow  = fabs(boxBoundaryFluxFEM( flow_domain, vt_key, inflow_boundary ));
    const double  outflow = fabs(boxBoundaryFluxFEM( flow_domain, vt_key, outflow_boundary ));

    // computing the maximum of the in- and outflow
    flux_through_model = (fabs(inflow) + fabs(outflow)) / 2.;
    // making sure that the difference between in- and outflux from the model is less than 1%
    const double flux_discrepancy_tolerance(0.03); // 3 percent.
    if ( fabs(fabs(inflow)-fabs(outflow)) / flux_through_model > flux_discrepancy_tolerance ) {
         cerr <<"\ninflux versus outflux: "<< inflow <<" vs. "<< outflow <<" m3/s.";
         csmp_error.Note( ERROR, "EquivalentPermeability", "flux analysis is inaccurate; check integrity of boundary mesh." );
         // using the maximum flux estimate in this case
         flux_through_model = std::max( fabs(inflow), fabs(outflow) );
      }
   
	  cout <<"\nequivalentPermeability: total flux through model: "<< flux_through_model;
	  cout <<" m3/s, flux per m2: "<< flux_through_model / xsect_area <<" m3/s."<< endl;
    cout <<"\n\testimated error of calculation (%): <= "<< 100. * (fabs(fabs(inflow)-fabs(outflow)) / flux_through_model) << endl << endl;
   
    // recovering fluid viscosity
    double  mumin, fluid_viscosity;
    flow_domain.MinMaxOf( "viscosity", mumin, fluid_viscosity );
    // asserting that there is only a single viscosity value
    assert( fabs(fluid_viscosity - mumin) <= numeric_limits<double>::epsilon() );
   
    // equivalent permeability
	  return (flux_through_model * fluid_viscosity) / (xsect_area * farfield_pf_gradient);

 } // end equivalentPermeability (using box-model boundaries)

template double equivalentPermeability<2>( Model<2>&, const char*, BOX_BOUNDARY, BOX_BOUNDARY, double&  );
template double equivalentPermeability<3>( Model<3>&, const char*, BOX_BOUNDARY, BOX_BOUNDARY, double&  );



} // end csmp



