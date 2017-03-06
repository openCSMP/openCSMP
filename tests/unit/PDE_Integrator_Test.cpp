#include "PDE_Integrator_Test.h"

using namespace std;

namespace csmp {

PDE_Integrator_Test::PDE_Integrator_Test() {

    const bool       isoparametric(true);
    ANSYS_Interface  mesh_interface(isoparametric);  // true = isoparametric elements
    VSet<2U>         mesh_container;
    ModelTopology    mesh_topology(isoparametric);   // true = isoparametric elements

    // Building Region object from ANSYS data files
    cout <<"Reading mesh..."<<endl;
    string mesh_name("pde_integrator_test");
    const bool binary_file( true );
    const bool irregular_mesh( false );
    mesh_interface.Read_ANSYS_Mesh( mesh_name.c_str(), mesh_container, mesh_topology, binary_file, irregular_mesh );
    cout <<"Finished reading mesh..."<<endl;
    cout <<"Building Model..."<<endl;
    sg_= new Model<2U> ( mesh_topology, mesh_container, "CSMP-2phase-variables.txt");

    // Set values on nodes and elements
    sg_->InputPropertyValue("fluid pressure", makeScalar(PLAIN,1.));
    sg_->InputPropertyValue("diffusivity", makeScalar(PLAIN,0.));
    sg_->InputPropertyValue("permeability", makeScalar(PLAIN,1.));

    // Set boundary conditions
    sg_->InputBoundaryValue( BOTTOM, "fluid pressure", makeScalar(DIRICH, 0.0) );
    sg_->InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH, 1.0) );
    sg_->InputBoundaryValue( TOP, "fluid pressure", makeScalar(DIRICH, 0.0) );
    sg_->InputBoundaryValue( LEFT, "fluid pressure", makeScalar(DIRICH, 0.0) );

    // Create pde-operators
    source_ = new Integral_NT_op_N_dV<2U,Element<2U> > (sg_->Database(),
                                                     "diffusivity", "fluid pressure");
                                                                               
    stiff_ = new Integral_dNT_op_dN_dV<2U,Element<2U> > (sg_->Database(),
                                                     "permeability",
                                                     "fluid pressure",
                                                     "fluid pressure");
    
    // Create algorithm object and add pde-operators
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    settings.Set_eps(0.);
    alg_=new PDE_Integrator<2U,Region>(new SAMG_Solver(&settings));
    #else
    alg_=new PDE_Integrator<2U,Region>(new CSMP_DEFAULT_LINEAR_SOLVER() );
    #endif
    alg_->Add(source_);
    alg_->Add(stiff_);
}



PDE_Integrator_Test::~PDE_Integrator_Test()
{
  if( sg_ != 0 )
    delete sg_;
  if( source_ != 0 )
      delete source_;
  if( stiff_ != 0 )
      delete stiff_;
  if( alg_ != 0 )
    delete alg_;
}



void PDE_Integrator_Test::run()
{
//    exchangeSolverTest();
    sameSolverTest();
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

void PDE_Integrator_Test::sameSolverTest()
{
  #ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Solver* samg = new SAMG_Solver();
  #else
  CSMP_DEFAULT_LINEAR_SOLVER* samg = new CSMP_DEFAULT_LINEAR_SOLVER();
  #endif
  alg_->SetSolver(samg);
  alg_->SetSolver(samg);
  _test(samg != NULL);
  delete samg;
}

} // end namespace csmp
