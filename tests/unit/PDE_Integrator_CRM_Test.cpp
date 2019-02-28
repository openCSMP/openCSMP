#include "PDE_Integrator_CRM_Test.h"

using namespace std;

namespace csmp {

PDE_Integrator_CRM_Test::PDE_Integrator_CRM_Test()
{
  const bool       isoparametric( true );
  ANSYS_Interface  mesh_interface( isoparametric );  // true = isoparametric elements
  VSet<2U>         mesh_container;
  ModelTopology    mesh_topology( isoparametric );   // true = isoparametric elements

                                                     // Building Region object from ANSYS data files
  cout << "Reading mesh..." << endl;
  string mesh_name( "pde_integrator_test" );
  const bool binary_file( true );
  const bool irregular_mesh( false );
  mesh_interface.Read_ANSYS_Mesh( mesh_name.c_str(), mesh_container, mesh_topology, binary_file, irregular_mesh );
  cout << "Finished reading mesh..." << endl;
  cout << "Building Model..." << endl;
  sg_ = new Model<2U>( mesh_topology, mesh_container, "CSMP-2phase-variables.txt" );

  // Set values on nodes and elements
  sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 1. ) );
  sg_->InputPropertyValue( "diffusivity", makeScalar( PLAIN, 0. ) );
  sg_->InputPropertyValue( "permeability", makeScalar( PLAIN, 1. ) );

  // Set boundary conditions
  sg_->InputBoundaryValue( BOTTOM, "fluid pressure", makeScalar( DIRICH, 0.0 ) );
  sg_->InputBoundaryValue( RIGHT, "fluid pressure", makeScalar( DIRICH, 1.0 ) );
  sg_->InputBoundaryValue( TOP, "fluid pressure", makeScalar( DIRICH, 0.0 ) );
  sg_->InputBoundaryValue( LEFT, "fluid pressure", makeScalar( DIRICH, 0.0 ) );

  // Create pde-operators
  source_ = new Integral_NT_op_N_dV<2U, Element<2U> >( sg_->Database(),
                                                       "diffusivity", "fluid pressure" );

  stiff_ = new Integral_dNT_op_dN_dV<2U, Element<2U> >( sg_->Database(),
                                                        "permeability",
                                                        "fluid pressure",
                                                        "fluid pressure" );

  // Create algorithm object and add pde-operators
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings  settings;
  settings.Set_eps( 0. );
  alg_ = new PDE_Integrator_CRM<2U, Region>( new SAMG_Solver( &settings ) );
#else
  alg_ = new PDE_Integrator_CRM<2U, Region>( new CSMP_DEFAULT_LINEAR_SOLVER() );
#endif
  alg_->Add( source_ );
  alg_->Add( stiff_ );

}

PDE_Integrator_CRM_Test::~PDE_Integrator_CRM_Test()
{
  if ( sg_ != 0 )
    delete sg_;
  if ( source_ != 0 )
    delete source_;
  if ( stiff_ != 0 )
    delete stiff_;
  if ( alg_ != 0 )
    delete alg_;
}



void PDE_Integrator_CRM_Test::run()
{

  sameSolverTest();
  //constructorTest();
  //copyTest();
  //AddOperatorLHSTest();
  //AddOperatorRHSTest();
  //AddPostProcessTest();
  ////AddBoundaryIntegralsTest();
  //TimeIncrementTest();
  //TransientTest();
  IntegrateOverTest();


  //SetSolverTest();
  //GetSolverTest();


}

/*
void PDE_Integrator_Test::exchangeSolverTest() {
try {
alg_->SetSolver(new SAMG_Solver());
sg_->Apply(*alg_);
alg_->SetSolver(new SAMG_Solver());
sg_->Apply(*alg_);

} catch(...) {
_fail("Error while calling solver in exchangeSolverTest");
}
_succeed();
}
*/

void PDE_Integrator_CRM_Test::sameSolverTest()
{
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Solver* samg = new SAMG_Solver();
#else
  CSMP_DEFAULT_LINEAR_SOLVER* samg = new CSMP_DEFAULT_LINEAR_SOLVER();
#endif
  alg_->SetSolver( samg );
  alg_->SetSolver( samg );
  _test( samg != NULL );
  delete samg;
}


void PDE_Integrator_CRM_Test::constructorTest() {

  // default constructor
  {

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    settings.Set_eps( 0. );
    SAMG_Solver solver( &settings );
#else 
    CSMP_DEFAULT_LINEAR_SOLVER solver();
#endif

    PDE_Integrator_CRM<2U, Region> pde;
    //_test(typeid(pde.GetSolver()) == typeid(solver));

    //_test(4 == pde.dim2_());
    //_test(0 == pde.dof_per_node_());
    //_test(false == pde.get_setup_established_());
    //_test(false == pde.get_retain_matrix_());
    //_test(true == pde.get_newed_Solver_object());
    //_test(0. == pde.get_time_increment_());
    //_test(true == pde.GetVerbose());
  }


  {

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    settings.Set_eps( 0. );
    SAMG_Solver solver( &settings );
#else 
    CSMP_DEFAULT_LINEAR_SOLVER solver();
#endif

    PDE_Integrator_CRM<2U, Region> pde( solver );
    //_test(4 == pde.dim2_());
    //_test(0 == pde.dof_per_node_());
    //_test(false == pde.get_setup_established_());
    //_test(false == pde.get_retain_matrix_());
    //_test(false == pde.get_newed_Solver_object());    // using reference
    //_test(0. == pde.get_time_increment_());
    //_test(true == pde.GetVerbose());
  }


  {

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    settings.Set_eps( 0. );
    SAMG_Solver solver( &settings );
#else 
    CSMP_DEFAULT_LINEAR_SOLVER solver();
#endif

    PDE_Integrator_CRM<2U, Region> pde( &solver );
    //_test(4 == pde.dim2_());
    //_test(0 == pde.dof_per_node_());
    //_test(false == pde.get_setup_established_());
    //_test(false == pde.get_retain_matrix_());
    //_test(false == pde.get_newed_Solver_object());
    //_test(0. == pde.get_time_increment_());
    //_test(true == pde.GetVerbose());
  }


}



void PDE_Integrator_CRM_Test::copyTest() {

#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings  settings;
  settings.Set_eps( 0. );
  SAMG_Solver solver( &settings );
#else 
  CSMP_DEFAULT_LINEAR_SOLVER solver();
#endif

  PDE_Integrator_CRM<2U, Region> pde( solver );
  //PDE_Integrator_CRM<2U, Region> copy(pde);

  // check deep copy 
  //_test(copy.GetSolver() != pde.GetSolver());
  //_test(&(copy.get_G_()) != &(pde.get_G_()));
  //_test(&(copy.get_rh_()) != &(pde.get_rh_()));
  //_test(&(copy.get_x_()) != &(pde.get_x_()));
  //_test(&(copy.get_lhs_operators_()) != &(pde.get_lhs_operators_()));
  //_test(&(copy.get_rhs_operators_()) != &(pde.get_rhs_operators_()));
  //_test(&(copy.get_basic_operands_()) != &(pde.get_basic_operands_()));
  //_test(&(copy.get_test_operands_()) != &(pde.get_test_operands_()));
  //_test(&(copy.get_postpro_operators_()) != &(pde.get_postpro_operators_()));
  //_test(copy.dim2_() == pde.dim2_());
  //_test(copy.dof_per_node_() == pde.dof_per_node_());
  //_test(copy.get_setup_established_() == pde.get_setup_established_());
  //_test(copy.get_retain_matrix_() == pde.get_retain_matrix_());

  //_test(copy.get_newed_Solver_object() != pde.get_newed_Solver_object());   // false true
  //_test(copy.get_time_increment_() == pde.get_time_increment_());
  //_test(copy.get_retain_matrix_() == pde.get_retain_matrix_());


  // TODO: need to assert  && beware
  /*
  _test(&(copy.get_RH_()) != &(pde.get_RH_()));
  _test(); // scale_factor
  _      ; // target

  */

}


void PDE_Integrator_CRM_Test::AddOperatorLHSTest() {
  PDE_Integrator_CRM<2U, Region> pde;

  auto& name = stiff_->Name();
  pde.Add( stiff_ );

  //_test(pde.get_lhs_operators_().size() == 1);
  //_test(pde.get_lhs_operators_()[name] == stiff_);

  // need adding more LHS operand


}


void PDE_Integrator_CRM_Test::AddOperatorRHSTest() {
  PDE_Integrator_CRM<2U, Region> pde;

  auto& name = source_->Name();
  pde.Add( source_ );

  //_test(pde.get_rhs_operators_().size() == 1);
  //_test(pde.get_rhs_operators_()[name] == source_);

  // need adding more RHS operand

}

void PDE_Integrator_CRM_Test::AddPostProcessTest() {

  PDE_Integrator_CRM<2U, Region> pde;
  auto& name = fluid_velocity->Name();
  pde.AddPostProcess( fluid_velocity );

  //_test(pde.get_postpro_operators_().size() == 1);
  //_test(pde.get_postpro_operators_()[name] == fluid_velocity);

  // need adding more RHS operand

}


void PDE_Integrator_CRM_Test::TimeIncrementTest() {
  PDE_Integrator_CRM<2U, Region> pde;

  double64 dt( 100 );
  pde.TimeIncrement( dt );
  //_test(pde.get_time_increment_() == dt);

  dt = 0.;
  pde.TimeIncrement( dt );
  //_test(pde.get_time_increment_() == dt);

  dt = -100.;
  pde.TimeIncrement( dt );
  //_test(pde.get_time_increment_() == dt);

}


void PDE_Integrator_CRM_Test::TransientTest() {
  PDE_Integrator_CRM<2U, Region> pde;

  double64 dt( 100 );
  pde.TimeIncrement( dt );
  _test( pde.Transient() == true );

  dt = 0.;
  pde.TimeIncrement( dt );
  _test( pde.Transient() == false );

  dt = -100.;
  pde.TimeIncrement( dt );
  _test( pde.Transient() == false );

}

void PDE_Integrator_CRM_Test::SetSolverTest() {

#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings  settings;
  settings.Set_eps( 0. );
  SAMG_Solver solver( &settings );
#else 
  CSMP_DEFAULT_LINEAR_SOLVER solver();
#endif

  PDE_Integrator_CRM<2U, Region>* pde = new PDE_Integrator_CRM<2U, Region>( solver );
  PDE_Integrator_CRM<2U, Region>* another = new PDE_Integrator_CRM<2U, Region>( solver );

  _test( pde->GetSolver() == another->GetSolver() );
  //_test(pde->get_newed_Solver_object() == false);

  delete pde;
  delete another;

}


void PDE_Integrator_CRM_Test::GetSolverTest() {

#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings  settings;
  settings.Set_eps( 0. );
  SAMG_Solver solver( &settings );
#else 
  CSMP_DEFAULT_LINEAR_SOLVER solver();
#endif

  PDE_Integrator_CRM<2U, Region>* pde = new PDE_Integrator_CRM<2U, Region>( solver );
  _test( pde->GetSolver() == &solver );

  delete pde;

}

void PDE_Integrator_CRM_Test::EstablishMatrixSetupTest() {
  PDE_Integrator_CRM<2U, Region>* pde = new PDE_Integrator_CRM<2U, Region>();
  pde->Add( stiff_ );
  pde->Add( source_ );





  delete pde;
}

void PDE_Integrator_CRM_Test::IntegrateOverTest() {
  {
    PDE_Integrator_CRM<2U, Region>* pde = new PDE_Integrator_CRM<2U, Region>();
    pde->Add( stiff_ );
    pde->Add( source_ );


    (sg_->Region( "Model" )).E( 0 )->FE()->UsesLocalCoordinates();

    pde->IntegrateOver( sg_->Region( "Model" ), true );



    delete pde;
  }

  {
    PDE_Integrator_CRM<2U, Region>* pde = new PDE_Integrator_CRM<2U, Region>();
    pde->Add( stiff_ );
    pde->Add( source_ );
    pde->IntegrateOver( *sg_, sg_->Region( "Model" ), true );
    delete pde;
  }


}



////////////////////////////////////////////////////


} // end namespace csmp
