#include "IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase.h"
#include "ModelComparator.h"

using namespace std;


namespace csmp
{

IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase::IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase(const char* prefix)
{
    this->setName("IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase");
    prefix_=prefix;
}

/* Finite Element Incompressible Single Phase Porous Media Flow Test Case
 * This variant of the test is in 3D and uses a tensor permeability
  =================================
  Element type(s) : Linear Tets
  Target Mesh set : UnitCubeFine
  Test            : Steady State Pressure (Incompressible Single Phase Porous Media Flow)
                    1) Dirichlet (const Pressure) on bottom-left and top-right corners.
                    Criterion:  Comparison with analytical solution: pressure.

  =================================
  */

void IncompressibleSinglePhaseFlowTensorPerm3DFEM_VVCase::run()
{
    // some constants
    enum{DIM=3};
    string input_file_name(prefix_);
    cout <<"\nRunning Test Case Simulation - "<<this->getName()<<endl;


    // ----------------------------------
    // GEOMETRY SECTION
    // establishing model & output facility
    cout <<"Building Model..."<<endl;
    ANSYS_Model3D model( input_file_name.data(),this->getName().c_str(),(this->getName()+".txt").c_str() );
    cout <<"Finished reading mesh..."<<endl;

    // print model dimensions
    Point<DIM> min, max;
    model.MinMaxCoordinates( min, max );
    printModelDimensions( model, true );
    double domain_volume = model.Region("Model").Volume();
    cout <<"\nThe model has a volume of: "<< domain_volume <<" m^3."<< endl;
    // END GEOMETRY SECTION
    // -----------------------------------

    // -----------------------------------
    // INPUT SECTION
    // model configuration

    InputDataManager<DIM>  model_configuration;
    model_configuration.ConfigureFromFile( model, this->getName().c_str(),false, true, true, true, false );
//    vector<Node<DIM>*>::iterator nodes_end=model.Region("Model").NodesEnd();
//    vector<Node<DIM>*>::iterator nodes_begin=model.Region("Model").NodesBegin();
//    csmp::Index  p_key(model.Database().StorageKey("fluid pressure"));
//    for (vector<Node<DIM>*>::iterator npit= nodes_begin; npit!=nodes_end;npit++)
//    {
//        Node<DIM> n;
//        if ((*npit)->x()<=0.01 && (*npit)->y()<=0.01 && (*npit)->z()<=0.01)
//            (*npit)->Store(p_key,makeScalar(DIRICH,101325.));

//        if ((*npit)->x()>=0.99 && (*npit)->y()>=0.99 && (*npit)->z()>=0.99)
//            (*npit)->Store(p_key,makeScalar(DIRICH,-101325.));

//    }


    //setup of output of initial conditions
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
    SAMG_Solver   solver( &settings );
    //samgsettings.Set_napproach(2);
    settings.Set_eps(1.0e-14);
    settings.Set_rel_eps(1.0e-12);
    #else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    #endif

    PDE_Integrator<DIM,Element> pressure_diffusion( solver );

    NumIntegral_dNT_op_dN_dV<DIM> stiffness( model.Database(), "mobility", "fluid pressure",  "fluid pressure");
    printRangeOfVariable(model,"fluid pressure");
    printRangeOfVariable(model,"fluid volume source");
    printRangeOfVariable(model,"mobility");
    NumIntegral_NT_op_N_dV<DIM> fluid_src( model.Database(),  "fluid volume source", "fluid pressure" );

    //Assemble the integrator
    pressure_diffusion.Add( &stiffness );
    pressure_diffusion.Add( &fluid_src );
    pressure_diffusion.IntegrateOver( model.Region("Model") );
    // END CORE SECTION
    //------------------------------------

    //------------------------------------
    // OUTPUT SECTION
    vtu.OutputDataToVTU(this->getName().c_str(), outputProps, "Model", static_cast<int>(0) );

    model.OutputToBinaryFile(this->getName().c_str());//csmp binary results.

    // END OUTPUT SECTION
    //------------------------------------

    //------------------------------------
    // VERIFICATION SECTION
    const csmp::Index  p_key(model.Database().StorageKey("fluid pressure"));
    auto nodes_end=model.Region("Model").NodesEnd();
    auto nodes_begin=model.Region("Model").NodesBegin();
    double press;
    for (auto npit= nodes_begin; npit!=nodes_end;npit++)
    {
      press=(*npit)->Read(p_key);
      //if(verbose_)
        cout<<"Press "<<press<<" analytical "<<1000000.0-898675.0*((*npit)->x())<<" difference "<<press-(1000000.0-898675.0*((*npit)->x()))<<endl;
      _equal( press, 1000000.0-898675.0*((*npit)->x()), 10.e-7 );
    }

    // END VERIFICATION SECTION
    //------------------------------------

    //------------------------------------
    // VERIFICATION SECTION (This commented section is left for reference to user the ModelComparator)

    // Can also compare models without output to disk (and simply load the comparison results from disk)
    // but this allows for flexibility when looking at why tests fail.

    //    ModelComparator<DIM> comparitor;
    //    double shouldBeZero( comparitor.CompareVSets( (this->getName()+".vset").c_str(), (this->getName()+"_Comparison.vset").c_str(),
    //                                                    "fluid pressure", "fluid pressure", (this->getName()+".txt").c_str(),
    //                                                    (this->getName()+".txt").c_str(), true ) );

    //    cout<<"L2 norm of difference for comparison: "<<shouldBeZero<<endl;

    //    _equal( shouldBeZero, 0., 1.0E-10 );

    // END VERIFICATION SECTION
    //------------------------------------


    return;
}

} // csmp
