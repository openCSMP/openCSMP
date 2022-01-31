#include "PDE_Integrator_Test.h"

using namespace std;

namespace csmp {

PDE_Integrator_Test::PDE_Integrator_Test()
 {
    const bool       isoparametric(true);
    ANSYS_Interface  mesh_interface(isoparametric);  // true = isoparametric elements
    VSet<2U>         mesh_container;
    ModelTopology    mesh_topology(isoparametric);   // true = isoparametric elements

    // Building Region object from ANSYS data files
    cout <<"Reading mesh..."<<endl;
    string mesh_name("pde_integrator_test");
    const bool binary_file( true );
    mesh_interface.Read_ANSYS_Mesh( mesh_name.c_str(), mesh_container, mesh_topology, binary_file, true );
    cout <<"Finished reading mesh..."<<endl;
    cout <<"Building Model..."<<endl;
    sg_= new Model<2U>( mesh_topology, mesh_container, "CSMP-2phase-variables.txt");

    // Set values on nodes and elements
    //sg_->InputPropertyValue("fluid pressure", makeScalar(PLAIN,1.));
    sg_->InputPropertyValue("diffusivity", makeScalar(PLAIN,1.0e-9));
    sg_->InputPropertyValue("permeability", makeScalar(PLAIN,1.0e-12));

    // Set boundary conditions
    sg_->InputBoundaryValue( BOTTOM, "fluid pressure", makeScalar(DIRICH, 0.0) );
    sg_->InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH, 1.0) );
    sg_->InputBoundaryValue( TOP, "fluid pressure", makeScalar(DIRICH, 0.0) );
    sg_->InputBoundaryValue( LEFT, "fluid pressure", makeScalar(DIRICH, 0.0) );

}



PDE_Integrator_Test::~PDE_Integrator_Test()
{
   delete sg_;
   delete source_;
   delete stiff_;
   delete alg_;
}



void PDE_Integrator_Test::run()
{
    SolveMatrixEquationWithSAMG();
}




void PDE_Integrator_Test::SolveMatrixEquationWithSAMG()
 {
    NumIntegral_NT_op_N_dV<2U>    diff(sg_->Database(), "diffusivity", "fluid pressure");
    NumIntegral_dNT_op_dN_dV<2U>  conductance(sg_->Database(), "permeability", "fluid pressure", "fluid pressure");
    
    // Create algorithm object and add pde-operators
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    settings.SetSolverInstance(1);
    settings.Set_eps(0.);
    SAMG_Solver samg(&settings);
    PDE_Integrator<2U,Region> integrator(samg);

    integrator.Add( &diff );
    integrator.Add( &conductance );

    integrator.IntegrateOver(sg_->Region("Model"));
    #endif
 }


} // end namespace csmp
