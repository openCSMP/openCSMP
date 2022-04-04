#include "IncompressibleSinglePhaseFlowFEM_VVCase.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"

using namespace std;


namespace csmp
{

IncompressibleSinglePhaseFlowFEM_VVCase::IncompressibleSinglePhaseFlowFEM_VVCase( const char* prefix, bool verbose )
  : verbose_(verbose)
{
  this->setName("IncompressibleSinglePhaseFlowFEM_VVCase");
  prefix_=prefix;
}

/* Finite Element Incompressible Single Phase Porous Media Flow Test Case
  =================================
  Mesh:       Linear Triangles
  Model:      40x2 m 2D
  Test:       Steady State Pressure (Incompressible Single Phase Porous Media Flow)
  BC:         1) Dirichlet (const Pressure) on RIGHT and LEFT boundaries. TOP and BOTTOM natural.
              2) Dirichlet (const Pressure) on bottom-left and top-right corners.
              3) Dirichlet (const Pressure) on bottom-left corner, NEUMANN (const rate)
  Criterion:  Comparison with analytical solution: pressure, flow velocity.

  =================================
  */

void IncompressibleSinglePhaseFlowFEM_VVCase::run()
{
  // some constants
  enum{DIM=2};
  const double tensileForceInX( 50000. );
  const ScalarVariable zeroScalar( PLAIN, 0. );
  const ScalarVariable zeroScalarDirichlet( DIRICH, 0. );
  const VectorVariable<DIM> zeroVector( PLAIN, 0. );
  const VectorVariable<DIM> zeroVectorDirichlet( DIRICH, 0. );
  string input_file_name(prefix_);
  cout <<"\nRunning Test Case Simulation - "<<this->getName()<<endl;


  // ----------------------------------
  // GEOMETRY SECTION
  // establishing model & output facility
  cout <<"Building Model..."<<endl;
  ANSYS_Model2D model( input_file_name.data(),this->getName().c_str(),(this->getName()+".txt").c_str() );
  cout <<"Finished reading mesh..."<<endl;

  // print model dimensions
  Point<DIM> min, max;
  model.MinMaxCoordinates( min, max );
  const double length( max[0]-min[0] ), height( max[1]-min[1] );
  printModelDimensions( model, true );
  double domain_volume = model.Region("Model").Volume();
  cout <<"\nThe model has a volume of: "<< domain_volume <<" m^3."<< endl;
  // END GEOMETRY SECTION
  // -----------------------------------

  // -----------------------------------
  // INPUT SECTION
  // model configuration
  InputDataManager<2U>  model_configuration;
  model_configuration.ConfigureFromFile( model, this->getName().c_str(),false, true, true, true, false );
  /*
  model.InputPropertyValue( "mobility", makeScalar( PLAIN, 1.0e-5 ) );
  model.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0e-5 ) );
  model.InputPropertyValue( "fluid volume source", makeScalar( PLAIN, 0.0e-5 ) );
  const ScalarVariable fluid_pressure_left( DIRICH, 101325.0 );
  Boundary<DIM>& rightBoundary( model.Boundary( "RIGHT" ) );
  Boundary<DIM>& leftBoundary( model.Boundary( "LEFT" ) );

  rightBoundary.InputPropertyValue( "fluid pressure", zeroScalarDirichlet );
  leftBoundary.InputPropertyValue( "fluid pressure",fluid_pressure_left);
  double min1,max1;
  rightBoundary.MinMaxOf("fluid pressure",min1,max1);
  std::cout<<"min max at right boundary"<<min1<<" "<<max1<<endl;
  leftBoundary.MinMaxOf("fluid pressure",min1,max1);
  std::cout<<"min max at left boundary"<<min1<<" "<<max1<<endl;
  */

  //output of initial conditions
  VTU_Interface<DIM> vtu( model );
  vtu.OmitZeroInFileName(true);
  list<string> outputProps;
  outputProps.push_back( "fluid pressure" );
  outputProps.push_back( "velocity" );
  vtu.OutputDataToVTU( "InitialConditions", outputProps, "Model", static_cast<int>(0) );

  // END INPUT SECTION
  // -----------------------------------


  // -----------------------------------
  // CORE SECTION
  // setting up & solving linear pressure diffusion
  #ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings settings;
  settings.Set_napproach(2);
  SAMG_Solver solver(&settings);
  #else
  CSMP_DEFAULT_LINEAR_SOLVER solver;
  #endif
  PDE_Integrator<DIM,Region> pressure_diffusion(solver);

  /*some old code
  NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  steady_conductance( p_ref, "conductivity", "hydrostatic pressure", "hydrostatic pressure" );
  NumIntegral_NT_op_dNi_dV<DIM,Element<DIM> >  gravity( p_ref, "element fluid density", "conductivity", "hydrostatic pressure" );
  NumIntegral_NT_op_N_dV<DIM,Element<DIM> >    fluid_src( p_ref,  "fluid volume source", "hydrostatic pressure" );
  */

  NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> > stiffness( model.Database(), "mobility", "fluid pressure",  "fluid pressure");
  printRangeOfVariable(model,"fluid pressure");
  printRangeOfVariable(model,"fluid volume source");
  printRangeOfVariable(model,"mobility");
  NumIntegral_NT_op_N_dV<DIM,Element<DIM> > fluid_src( model.Database(),  "fluid volume source", "fluid pressure" );
  //Assemble the integrator
  pressure_diffusion.Add( &stiffness );
  pressure_diffusion.Add( &fluid_src );
  pressure_diffusion.IntegrateOver( model.Region("Model") );
  singlePhaseVelocity(model,"Model","velocity","fluid pressure","mobility");
  // END CORE SECTION
  //------------------------------------

  //------------------------------------
  // OUTPUT SECTION
  //model.InputPropertyValue( "mean stress", zeroScalar );
  vtu.OutputDataToVTU( "FinalOutput", outputProps, "Model", static_cast<int>(0) );
  // END OUTPUT SECTION
  //------------------------------------

  //------------------------------------
  // VERIFICATION SECTION
  csmp::Index  p_key(model.Database().StorageKey("fluid pressure"));
  auto nodes_end=model.Region("Model").NodesEnd();
  auto nodes_begin=model.Region("Model").NodesBegin();
  double press;
  for ( auto npit= nodes_begin; npit!=nodes_end;npit++)
  {
    press=(*npit)->Read(p_key);
    if(verbose_)
      cout<<"Press "<<press<<" analytical "<<1000000.0-898675.0*((*npit)->x())<<" difference "<<press-(1000000.0-898675.0*((*npit)->x()))<<endl;
    _equal( press, 1000000.0-898675.0*((*npit)->x()), 10.e-8 );
  }

  // END VERIFICATION SECTION
  //------------------------------------

  return;
}

} // csmp
