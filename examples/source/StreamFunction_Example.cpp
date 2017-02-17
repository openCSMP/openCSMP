#include "StreamFunction_Example.h"

#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "variableOperations.h"

// Mesh Import & Property Assignment
#include "TRIANGLE_Interface.h"
#include "VSet.h"
#include "VSetConverter.h"
//#include "BinaryFileInterface.h"
#include "InputDataManager.h"

// PDE Operators
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "NumIntegral_op_NT_dN_orthogonal_dV.h"

// Interrelations
#include "ConstantFactor.h"

// Other
#include "PropertyConstraints.h"
#include "VTK_Interface.h"
#include "Standard_IO_Handler.h"
#include "PropertyHandle.h"
#include "ModelTime.h"

// Analysis
#include "StatisticalAnalyzer.h"
#include "RegionMonitor.h"

using namespace std;

namespace csmp{

void StreamFunction_Example::Specifications()
{
  SetTitle( "Streamfunction from velocity field. Function can be contoured to get streamlines" );
  SetDifficulty( 3 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "post-processing of the velocity field using the streamfunction that can be contoured" );
  AddDescription( "steady-state analysis of flow and fluid pressure in heterogeneous medium" );
  AddDescription( "source in StreamFunction_Example.cpp" );
  AddRequirement( "file set: 'example21.1'");
}



/** *****************************************************************************************

  Analysis of steady-state fluid pressure and velocity distributions in fractured medium
  using quadratic FEM.
 
  Post-processing of the streamfunction that can be contoured to get streamlines.
 
  Example also allows to vary the permeability of model subregions (addressed by 
  their name) and calculates its effective permeability.

  Try 'example21.1'  input file set.

***************************************************************************************** */
void StreamFunction_Example::Run()
{
   TRIANGLE_Interface  mesh_interface;
   VSet<2U>            mesh_container;

   char   file_name[200];
   cout <<"\nmain: Enter name of 'Triangle' input file set: ";
   cin >> file_name;
   string model_name(file_name);
   model_name += "simulation.txt";

   Standard_IO_Handler  stdio( model_name.c_str() );

   mesh_interface.ReadTriangle2DMesh( file_name, mesh_container );

   VSetConverter<2U>  converter;
   converter.ConvertLinearToQuadraticTriangles( mesh_container );

 // ------------------------------------------------------------------------------------
 // 0. building Region object
 // ------------------------------------------------------------------------------------
   const bool isoparametric(true);
   Model<2U>  model( mesh_container, "example21.txt", isoparametric );
   mesh_container.Erase();


 // ------------------------------------------------------------------------------------
 // 1. Material properties, groups etc., initial & essential conditions &
 //    interrelations
 // ------------------------------------------------------------------------------------
   InputDataManager<2U>  model_configuration;
   model_configuration.ConfigureFromFile( model, file_name );

   printModelDimensions( model );


 // ------------------------------------------------------------------------------------
 // 2. Calculating hydraulic conductivity from permeability using Interrelation subclass
 // ------------------------------------------------------------------------------------
   const double64 fluid_viscosity(1.0e-03);
   ConstantFactor<2U,divides>  conductivity( model.Database(),
                                                     "conductivity", "permeability",
                                                      fluid_viscosity );
   model.Apply( conductivity );
   printRangeOfVariable( model, "conductivity" );


 // ------------------------------------------------------------------------------------
 // 3. Steady-state fluid pressure computation [K]{p} = {Q}
 // ------------------------------------------------------------------------------------
   #ifdef CSMP_WITH_SAMG_SOLVER
   PDE_Integrator<2U,Region>  fluid_pressure(new SAMG_Solver());
   #else
   PDE_Integrator<2U,Region>  fluid_pressure(new CSMP_DEFAULT_LINEAR_SOLVER());
   #endif

   NumIntegral_dNT_op_dN_dV<2U,Element<2U> >  conductance( model.Database(), "conductivity",   "fluid pressure", "fluid pressure" );
   NumIntegral_NT_op_N_dV<2U,Element<2U> >    source( model.Database(), "fluid volume source", "fluid pressure" );
   VelocityAndVolumeFlux<2U,Element<2U> >     velo( model, "conductivity", "porosity", "fluid pressure" );

   fluid_pressure.Add( &conductance );
   fluid_pressure.Add( &source );
   fluid_pressure.AddPostProcess( &velo );
   model.Apply( fluid_pressure );


 // ------------------------------------------------------------------------------------
 // 3.1 Analysis of results
 // ------------------------------------------------------------------------------------
   StatisticalAnalyzer<2U>                                     flux_histogram( model );
   vector<pair<double64,double64> >                            bins;
   map<string,pair<vector<pair<double64,double64> >,size_t> >  results;

   flux_histogram.DefineBins( "example21_histogram.bins", bins );
   flux_histogram.RegionPropertyHistograms( "volume flux", bins, results );
   flux_histogram.OutputRegionPropertyAbundancePolygonsMaple( "volume flux", bins, results, true ); // log10 of data

   printRangeOfVariable( model, stdio, "fluid pressure", true );
   printRangeOfVariable( model, stdio, "velocity", true );
   printRangeOfVariable( model, stdio, "volume flux", true );

   list<string> integral_properties;  integral_properties.push_back("volume flux");
   list<string> range_properties;     range_properties.push_back("fluid pressure");

   RegionMonitor<2U>  monitor( model, integral_properties, range_properties );
   bool normalize_integrated_values = false;
 //	    stdio.RecordLogicalChoice("Would you like to normalize region integral by region volume");
   double64& model_time( ModelTime::Instance().modelTime );  model_time = 0.;
   monitor.DivideIntegralPropertiesByRegionVolumes( normalize_integrated_values );
   monitor.ScalarPropertyIntegrals( model, model_time );
   monitor.ScalarPropertyRanges( model, model_time );
   monitor.Out( "CSMP-monitored-regions" );


// TESTING VTK OUTPUT OF QUADRATIC TRIANGLE ELEMENTS
const csmp::Index pf_key = model.Database().StorageKey("fluid pressure");
Region<2U>  rref = model.Region("Model");
rref.E(5)->CoordinateMatrix();
DenseMatrix<DM_MIN> DATA;
DATA.Resize(1,6);
// nodal fluid pressure
for ( size_t j=0; j<rref.E(5)->Nodes(); j++ ) DATA(0,j) = rref.E(5)->N(j)->Read( pf_key );
rref.E(5)->FE()->OutputNodeDataToVTK( "test_e", "fluid_pressure", DATA );


 // ------------------------------------------------------------------------------------
 // 4. effective permeability
 // ------------------------------------------------------------------------------------
   cout <<"\nmain: Enter direction of flow (parallel-to-X-axis=1, Y-axis=2): ";
   size_t  flow_direction;
   cin >> flow_direction;
   Point<2U>  xyz_min, xyz_max;
   model.Region("Model").MinMaxCoordinates( xyz_min, xyz_max );
   const double64  xsec_X = xyz_max[0] - xyz_min[0];
   const double64  xsec_Y = xyz_max[1] - xyz_min[1];
   double64  pmin, pmax, farfield_pf_gradient,
             // non signalling initialization to NAN value
             k_effective(std::numeric_limits<double64>::quiet_NaN());

   double64  model_throughput = integrateDomainBoundaryFlux( model );

   model.MinMaxOf( "fluid pressure", pmin, pmax );

   if ( flow_direction == 1U ) {
        farfield_pf_gradient = (pmax - pmin) / xsec_X;
        // from Darcy's law
        //        effective permeability = (q * mu) / grad_P
        k_effective  = model_throughput / farfield_pf_gradient;
        k_effective /= xsec_Y;
     }
   else if ( flow_direction == 2U ) {
        farfield_pf_gradient = (pmax - pmin) / xsec_Y;
        // from Darcy's law
        //        effective permeability = (q * mu) / grad_P
        k_effective  = model_throughput / farfield_pf_gradient;
        k_effective /= xsec_X;
     }
   k_effective *= fluid_viscosity;
   cout <<"\nmain: The (flux-integrated) effective permeability of the model is: "<< k_effective << endl;
   cout <<"\nmain: The permeability range of the model was:"<< endl;
   printRangeOfVariable( model, stdio, "permeability", true );

   // output results for visualization using VTK
   VTK_Interface<2U>  vtk_output;
   vtk_output.OutputDataToVTK( model, "conductivity", "conductivity", 0 );
   vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
   vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       0 );
   vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux", 0 );


 // ------------------------------------------------------------------------------------
 // 4. Compute the streamfunction which corresponds to solution
 // ------------------------------------------------------------------------------------
   if ( flow_direction == 1 ) computeStreamFunction( model, BOTTOM, TOP, model_throughput, "stream function" );
   else                       computeStreamFunction( model, LEFT, RIGHT, model_throughput, "stream function" );
   double64 max_stream = printRangeOfVariable( model, stdio, "stream function", true );
   cout <<"\nmain: Total throughput versus max indicated by streamfunction solution: ";
   cout << model_throughput <<" vs. "<< max_stream <<" (should be the same)."<< endl;

   // dependent variable ranges in groups
   map<string,Region<2> >::const_iterator  grit(model.RegionsBegin());
   while ( grit!=model.RegionsEnd() ) {
        printRangeOfVariable( model, stdio, (*grit).first.c_str(), "fluid pressure" );
        printRangeOfVariable( model, stdio, (*grit).first.c_str(), "velocity" );
        printRangeOfVariable( model, stdio, (*grit).first.c_str(), "volume flux" );
        grit++;
     }

   // output results for visualization using VTK
   vtk_output.OutputDataToVTK( model, "stream-function", "stream function",  0 );

 //  BinaryFileInterface<double64,2U>  bin_output;
 //  bin_output.WriteConnectivityFile( model, "model_connectivity" );
 //  bin_output.WriteDataTo( model, "fluid-pressure", "fluid pressure", 0 );


 // ------------------------------------------------------------------------------------
 // 4. Sensitivity test for a single group
 // ------------------------------------------------------------------------------------
   string group_name;

   if ( stdio.RecordLogicalChoice("Would you like to vary permeability in a single region") ) {
        cout <<"\nmain: Enter name of the desired region: ";
        cin >> group_name;

        printRangeOfVariable( model, stdio, group_name.c_str(), "fluid pressure" );
        printRangeOfVariable( model, stdio, group_name.c_str(), "velocity" );
        printRangeOfVariable( model, stdio, group_name.c_str(), "volume flux" );

        analyze_sensitivity( model, group_name.c_str(), stdio, conductivity, fluid_pressure );

        if ( stdio.RecordLogicalChoice("Output last result from sensitivity analysis ?") ) {
             vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure",    0 );
             vtk_output.OutputDataToVTK( model, "velocity",       "velocity",          0 );
             vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux",       0 );
          }
     }

   cout <<"\nmain: That's it..."<< endl;


} // Run()






// auxiliary methods

// this method assumes that permeability is a scalar
void  StreamFunction_Example::analyze_sensitivity( Model<2U>& sg, const char* group, Standard_IO_Handler& io,
                                                   Interrelation<2U>& itr, PDE_Integrator<2U,Region>& algo )
 {
    assert( sg.Database().Type("permeability") == SCALAR );
    for ( ;; ) {
          double64 average_k = sg.Region(group).Average( "permeability" );
          cout <<"\nanalyze_sensitivity: Current average permeability, k = "<< average_k;
          cout <<"; enter new k value (-1. to break loop): ";
          ScalarVariable  perm;
          cin >> perm();
          if ( perm() < 0. ) break;
          sg.Region(group).InputPropertyValue( "permeability", perm, COMPLETE );
          sg.Apply( itr );
          sg.Apply( algo );
          printRangeOfVariable( sg, io, group, "fluid pressure" );
          printRangeOfVariable( sg, io, group, "velocity" );
          printRangeOfVariable( sg, io, group, "volume flux" );
      }
 }

/**

bool integrateDomainBoundaryFlux( Region<3>& sg,
                                  double64& inflow, double64& outflow )

@section description Description

For all boundary elements of the current mesh,
from the fluid pressures at the nodes, velocities are found at the
integration points. These velocities are extrapolated to the nodes
points of each face and then integrated
by the faces over the faces.

Method returns total throughput through the model.

@return The results of the computation are returned into second and third
method arguments. They represent the total inflow and outflow of the
model, respectively.

@section application Application

This method is used to calculate the effective hydraulic conductivity
of the model domain.

tested: */
template<size_t dim>
double64 StreamFunction_Example::integrateDomainBoundaryFlux( Model<dim>& sg )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  const Region<dim>&  gref(sg.Region("Model"));
  PropertyHandle<dim>         face_normal( sg, "face normal", VECTOR, FACE );

  csmp::Index  pf_key=sg.Database().StorageKey("fluid pressure"),
               K_key= sg.Database().StorageKey("conductivity"),
               ve_key=sg.Database().StorageKey("nodal velocity"),
               fn_key=sg.Database().StorageKey("face normal");

  if ( pf_key.place != NODE || pf_key.type != SCALAR ) {
       throw csmp::Exception( ERROR, "integrateDomainBoundaryFlux",
                                  "fluid pressure variable must be a scalar placed on the nodes");
       return false;
    }
  if ( ve_key.place != NODE || ve_key.type != VECTOR ) {
       throw csmp::Exception( ERROR, "integrateDomainBoundaryFlux",
                                  "velocity variable must be a vector placed on the nodes");
       return false;
    }
  if ( fn_key.place != FACE || fn_key.type != VECTOR ) {
       throw csmp::Exception( ERROR, "integrateDomainBoundaryFlux",
                                  "face normal must be a vector placed on the face");
       return false;
    }

  // computing velocities at integration points and extrapolating them to the nodes
  // ------------------------------------------------------------------------------
  vector<ScalarVariable>  PF(10);
  ScalarVariable          K; // hydraulic conductivity
  vector<double64>        IPVF(4), NVF(10);
  DenseMatrix<DM_MIN>     DERIV;
  // storage for the extrapolated nodal velocities for each individual element
  // elements nodes velocity components
  map<size_t,vector<vector<double64> > >  NVELO;
  vector<vector<double64> >               dummy;
  pair<typename map<size_t,vector<vector<double64> > >::iterator,bool>  mit;

  for ( typename vector<Element<dim>*>::const_iterator
        eit=gref.ElementsBegin(); eit!=gref.ElementsEnd(); eit++ )
    if ( (*eit)->AtBoundary() != NOT )
      {
         // collecting nodal fluid pressures and elemental hydraulic conductivities from
         // each element
         (*eit)->NodePropertyVector( pf_key, PF );
         (*eit)->Read( K_key, K );

         DERIV.Resize(dim,(*eit)->Nodes());
         IPVF.resize( (*eit)->IntegrationPoints()*dim );
         fill( IPVF.begin(), IPVF.end(), 0. );

         for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
           {
              // getting shape function derivative matrix at the integration point
              (*(*eit)).dN_AtIntegrationPoint( DERIV, i, 1 );

              // getting the material properties, K and phi into the equation
              for ( size_t n=0; n<(*eit)->Nodes(); n++ )
                for ( size_t j=0; j<dim; j++ )
                  // -DERIV because fluid flows down pressure
                  IPVF[ i*dim + j ] += PF[n]() * -DERIV(j,n) * K();
           }

       // linear extrapolation of integration point velocities to the nodes
       NVF.resize( (*eit)->Nodes()*dim );
       (*eit)->ExtrapolateIntegrationPointVariableToNodes( dim, IPVF, NVF );

       // storing velocities in three dimensional map array for lookup by the faces (correct)
       mit = NVELO.insert( make_pair( (*eit)->Idx(), dummy ) );
       assert( mit.second );
       (*mit.first).second.resize((*eit)->Nodes());

       for ( size_t i=0U; i<(*eit)->Nodes(); i++ ) {
            (*mit.first).second[i].resize(dim);
            for ( size_t j=0; j<dim; j++ )
              (*mit.first).second[i][j] = NVF[ i*dim + j ];
         }

    } // end for all elements at model boundary


  // integrating fluxes over the faces and adding ensuing contributions to flux balances
  // -----------------------------------------------------------------------------------
  VectorVariable<dim>  normal, velo;
  double64             velocity, flux(0U);
  bool                 first_found(false), second_found(false);
  double64             sum_influx(0.), sum_outflux(0.);

  cout <<"\nintegrateDomainBoundaryFlux: Computing fluxes across faces..."<< endl;
  clock_t ticks = clock();

  for ( typename std::map<std::string,csmp::Boundary<dim> >::const_iterator
        bit=sg.BoundariesBegin(); bit!=sg.BoundariesEnd(); bit++ )
    for ( typename vector<Face<dim>*>::const_iterator
          fit=(*bit).second.ElementsBegin(); fit!=(*bit).second.ElementsEnd(); fit++ )
      /// Roman, 2014 ( Face&InterFace ): Should Face contain AtBoundary flag?
      //if ( (*fit)->AtBoundary() != NOT )
        {
           // computing and storing normal to face
           (*fit)->UnitNormal( normal );

           // finding the neighbor element, getting the velocity from it, and projecting it onto the facet normal
           if ( (*fit)->Parent(INSIDE) != NULL ) {
                first_found = true;
                (*fit)->Parent(INSIDE)->Read( ve_key, velo );
                // projecting flux onto facet normal
                velocity = dotProduct(normal,velo);
             }
           else { // the other neighbor pointer must not be zero
                assert( (*fit)->Parent(OUTSIDE) != NULL );
                second_found = true;
                (*fit)->Parent(OUTSIDE)->Read( ve_key, velo );
                velocity = dotProduct(normal,velo);
             }

           flux += (*fit)->Volume() * velocity;

           // recording the integrated fluxes (NOTE !: in and outflux distinction does not work this way)
           // -------------------------------------------------------------------------------------------
           // an outflux is given if the face has no neighbor in the direction where the flux
           // is going, else we are dealing with an influx
           // if there is no second element
           if ( first_found ) {
                if      ( flux > static_cast<double64>(0.) ) sum_outflux += flux;
                else if ( flux < static_cast<double64>(0.) ) sum_influx  += fabs(flux);
             }
           // second element, everything is just opposite
           else if ( second_found ) {
                if      ( flux > static_cast<double64>(0.) ) sum_outflux += flux;
                else if ( flux < static_cast<double64>(0.) ) sum_influx  += fabs(flux);
             }

      } // end for if faces

  // total throughput must be 1/2 of fabs() of what goes in and what goes out
  cout <<"\nintegrateDomainBoundaryFlux: Total cross-sectional throughput: "<< (sum_influx+sum_outflux)/2. << endl;

  ticks = clock() - ticks;
  cout <<"\nintegrateDomainBoundaryFlux: CPU ticks used for flux integration over faces: "<< ticks << endl;

   return (sum_influx + sum_outflux) / 2.;

} // end integrateDomainBoundaryFlux





// @test tested: O.K.
void StreamFunction_Example::computeStreamFunction( Model<2U>& sg,
                            BOX_BOUNDARY boundary0, BOX_BOUNDARY boundary1,
                            double64 total_flux, const char* stream_func_var )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();


    csmp::Index          res_key = sg.Database().StorageKey("resistivity");
    csmp::Index          con_key = sg.Database().StorageKey("conductivity");
    ScalarVariable       sc;
    Region<2>&  gref(sg.Region("Model"));

    for ( vector<Element<2U>*>::iterator
          eit=gref.ElementsBegin(); eit!=gref.ElementsEnd(); eit++ ) {
         (*eit)->Read( con_key, sc );
         sc = 1. / sc();
         (*eit)->Store( res_key, sc );
      }

    // assigning boundary conditions at one boundary which should be parallel to flow direction
    // from this boundary away the flux will be integrated by this method
    sg.InputBoundaryValue( boundary0, stream_func_var, makeScalar(DIRICH, 0.) );
    sg.InputBoundaryValue( boundary1, stream_func_var, makeScalar(DIRICH, total_flux) );

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver  samg_solver;
    PDE_Integrator<2U,Region>  stream_function(samg_solver);
#else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    PDE_Integrator<2U,Region>  stream_function(linear_solver);
#endif

    NumIntegral_dNT_op_dN_dV<2U,Element<2U> >  conductance( sg.Database(),
                                                        "resistivity", stream_func_var, stream_func_var );

// RENAME THIS OPERATOR into NumIntegral_NT_op_dN_orthogonal_dV
    NumIntegral_op_NT_dN_orthogonal_dV<2U,Element<2U> >  rhs( sg.Database(), "fluid pressure", stream_func_var );

    stream_function.Add( &conductance );
    stream_function.Add( &rhs );

    cout <<"\n\ncomputeStreamFunction: Computing stream function..."<< endl;
    sg.Apply( stream_function );

 } // end computeStreamFunction

} // csmp
