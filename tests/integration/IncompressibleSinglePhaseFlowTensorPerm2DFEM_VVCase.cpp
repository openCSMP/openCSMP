#include "IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase.h"
#include "ModelComparator.h"
using namespace std;


namespace csmp
{

IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase::IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase(const char* prefix)
{
    this->setName("IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase");
    prefix_=prefix;
}

/* Finite Element Incompressible Single Phase Porous Media Flow Test Case
  =================================
  Element type(s) : Linear Triangles
  Target Mesh set : simple_square_domain_unstructured_triangles
  Test            : Steady State Pressure (Incompressible Single Phase Porous Media Flow)
                    1) Dirichlet (const Pressure) on bottom-left and top-right corners.
                    Criterion:  Comparison with analytical solution: pressure.

  =================================
  */

void IncompressibleSinglePhaseFlowTensorPerm2DFEM_VVCase::run()
{
    // some constants
    enum{DIM=2};
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
    printModelDimensions( model, true );
    double64 domain_volume = model.Region("Model").Volume();
    cout <<"\nThe model has a volume of: "<< domain_volume <<" m^3."<< endl;
    // END GEOMETRY SECTION
    // -----------------------------------

    // -----------------------------------
    // INPUT SECTION
    // model configuration
    InputDataManager<DIM>  model_configuration;
    model_configuration.ConfigureFromFile( model, this->getName().c_str(),false, true, true, true, false );

    //output of initial conditions
    VTU_Interface<DIM> vtu( model );
    vtu.OmitZeroInFileName(true);
    list<string> outputProps;
    outputProps.push_back( "fluid pressure" );
    //outputProps.push_back( "velocity" );
    vtu.OutputDataToVTU( "InitialConditions", outputProps, "Model", static_cast<int>(0) );

    // END INPUT SECTION
    // -----------------------------------

    // -----------------------------------
    // CORE SECTION
    // setting up & solving linear pressure diffusion
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_napproach(2);
    settings.Set_eps(1.0e-14);
    settings.Set_rel_eps(1.0e-12);
    PDE_Integrator<DIM,Region> pressure_diffusion( new SAMG_Solver(&settings) );
    #else
    PDE_Integrator<DIM,Region> pressure_diffusion( new CSMP_DEFAULT_LINEAR_SOLVER() );
    #endif

    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> > stiffness( model.Database(), "mobility", "fluid pressure",  "fluid pressure");
    printRangeOfVariable(model,"fluid pressure");
    printRangeOfVariable(model,"fluid volume source");
    printRangeOfVariable(model,"mobility");
    NumIntegral_NT_op_N_dV<DIM,Element<DIM> > fluid_src( model.Database(),  "fluid volume source", "fluid pressure" );

    //Assemble the integrator
    pressure_diffusion.Add( &stiffness );
    pressure_diffusion.Add( &fluid_src );
    pressure_diffusion.IntegrateOver( model.Region("Model") );
    // END CORE SECTION
    //------------------------------------

    //------------------------------------
    // OUTPUT SECTION
    vtu.OutputDataToVTU(this->getName().c_str(), outputProps, "Model", static_cast<int>(0) );

    model.OutputToBinaryFile(this->getName().c_str());// csmp binary format.
    // END OUTPUT SECTION
    //------------------------------------

    //------------------------------------
    // VERIFICATION SECTION
    csmp::Index  p_key(model.Database().StorageKey("fluid pressure"));
    vector<Node<DIM>*>::iterator nodes_end=model.Region("Model").NodesEnd();
    vector<Node<DIM>*>::iterator nodes_begin=model.Region("Model").NodesBegin();
    double64 press;
    for (vector<Node<DIM>*>::iterator npit= nodes_begin; npit!=nodes_end;npit++)
    {
      press=(*npit)->Read(p_key);
      //if(verbose_)
        cout<<"Press "<<press<<" analytical "<<1000000.0-898675.0*((*npit)->x())<<" difference "<<press-(1000000.0-898675.0*((*npit)->x()))<<endl;
      _equal( press, 1000000.0-898675.0*((*npit)->x()), 10.e-9 );
    }

    // END VERIFICATION SECTION
    //------------------------------------


    return;
}

} // csmp
