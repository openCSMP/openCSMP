#include"SneddonCrack_VVCase.h"
#include "SplitBoundaryInterface_Test.h"
#include "string"
using namespace std;


namespace csmp
{



SneddonCrack_VVCase::SneddonCrack_VVCase(const char* prefix)
{
    this->setName("SneddonCrack_VVCase");
    prefix_=prefix;
}

/**

    Author:  E. Pezzulli  &  M. Nejati
    Linear Elasticity: Plane stress/strain problem with central crack under uniform tension; testing the aperture profile of the crack
    =================================
    Domain:     "Infinite" 2D rectangle shaped model length LxL meters with central crack of half length a where L/a > 20
    Mesh:       Linear/Quadratic Triangles tested OK; Quadrilaterials have to be tested
    Material:   Linear ELastic Isotropic Material: E=Young's Modulus, nu= Poisson's Ratio;
    Loading:    No stress on boundaries.
    BC:         Zero x and y displacements on LEFT BOTTOM Corner (CRN1), and zero y displacement on CRN2 (BOTTOM RIGHT)
    ConfigFiles:   (.asc, .dat)
    =================================

    The solution for the aperture w(x) of a horizontal fracture of half length a with origin x=0 representing middle of crack:

    Plane Strain:
    (aperture)    w = 2*Uy 4*(1 - v^2)*P/E * (a^2 - x^2 )^(1/2)
    E = youngs Modulus, v = Poisson Ration, a = crack half length, x = distance along crack with origin at midpoint


    Marji -  2006 - Crack Tip Elements - doi:10.1016/j.ijsolstr.2005.04.042
    Sneddon - 1946 - The opening of griffith cracks (factor of two wrong!)


    The following implementation is testing the creation of the Fracture Class, its initialisation and creation of InterFace Objects, constructed from
    opposite face pairs of the master and slave Boundaries of a 2D fracture. The calculation and calibration of a Unit Normal and the creation of midPoint
    Elements within the fracture, which are however not used/tested to work.

    @attention The following implemetation relies on LinearElasticityIsotropicDeformation2D_VVcase to work

    Assumptions of Test:
    1) Relies on Boundary Naming convention when created from lower dimensional region


*/




void SneddonCrack_VVCase::run()
{

  ///==========================================================================================
  /// Elasticity with Central Crack VV_Case
  ///==========================================================================================
  // Defining constants
    enum{dim=2};
    bool plane_strain = true; //no analy sol implemented for false
    const double ym (1.0e9), pr(0.3), P0(1.0e7);
                 //s_right(10.0e6), s_top(10.0e6), s_bottom(10.0e6), s_left(10.0e6);
    bool quarterpoint = false;

    // Model configuration:
    //Input file directory locations
    string input_dir         = "";
    string output_file       = "../Output/";
    output_file              += prefix_;

    //Variable inputs
    string vars_file         = input_dir + "SneddonConvergenceTest-variables.txt";
    string regions_file      = input_dir + "SneddonConvergenceTest";
    string mesh_file         = input_dir + prefix_;
    //string config_file            = input_dir +...;

    bool irregular(false), binary(true), regions(true);
    ANSYS_Model2D model( mesh_file.c_str(), regions_file.c_str(), vars_file.c_str(), irregular, binary, regions);    // Constractor for empty variables

    std::string sb_name = *(model.CreateSplitBoundaryFrom( "FRACTURE" ).first.begin()) ; //relies on split boundary naming convention
    std::string sb_reg  = *(model.InsertLowerDimensionalRegionsIntoSplitBoundaries(0).begin() );
    //_test( model.Region("Model").Contains( *(model.Region("FRACTURE_SPLIT_BOUNDARY_REGION").ElementsBegin()) ) );
    /// -------------------------------
    /// Setting up Fracture Configuration
    /// -------------------------------
    Fracture<dim> myFracture (model, sb_name, DC_TIP);            //Creates interface objects in Fracture, and initialises & configures Lubrication region
    //myFracture.TestSplitNodeAssignment(true);

    /// --------------------------------
    /// Elasto-Lubrication Parameter Configuration
    /// --------------------------------
    model.CreateProperty( "Neumann stress" , "tau", "SI", VECTOR, FACE);

    //Inputing Parameter Values

    model.InputPropertyValue("Young's modulus", ScalarVariable( PLAIN, ym ) );
    model.InputPropertyValue("Poisson's ratio", ScalarVariable( PLAIN, pr ) );
    model.InputPropertyValue("fluid pressure",  ScalarVariable( PLAIN, 0.0));     //setting all in model to Plain since decoupled
    model.InputPropertyValue("aperture",        ScalarVariable( PLAIN, 0.0 ) );
    model.InputPropertyValue("displacement",    VectorVariable<dim>(PLAIN,PLAIN,0.0,0.0)); //initialising displacement
    model.InputPropertyValue("Neumann stress",  VectorVariable<dim>(PLAIN,PLAIN,0.0,0.0)); //initialising displacement


    model.Region(sb_reg).InputPropertyValue("fluid pressure", ScalarVariable( PLAIN, P0),  COMPLETE);     //uniform fluid pressure accross fracture
    model.Region(sb_reg).InputPropertyValue("displacement",   VectorVariable<dim>(DIRICH,DIRICH, 0.0, 0.0),  COMPLETE);     //zero displacement on fracture


    /// --------------------------------
    /// Elasticity Boundary Conditions
    /// --------------------------------
    //Elasticity Dirichlet
    Node<dim> *cornerNodeLeftBottom, *cornerNodeRightBottom;
    for (vector<Node<dim>*>::const_iterator node (model.Region("Model").NodesBegin()); node != model.Region("Model").NodesEnd(); ++node )
    {
        if ( (*node)->AtBoundary() == CNR1 ) cornerNodeLeftBottom = (*node);
        if ( (*node)->AtBoundary() == CNR2 ) cornerNodeRightBottom = (*node);
    }

    cornerNodeLeftBottom->Store(model.Database().StorageKey("displacement"),VectorVariable<dim> (DIRICH, DIRICH , 0., 0.));   // when using boundary faces
    cornerNodeRightBottom->Store(model.Database().StorageKey("displacement"), VectorVariable<dim>(PLAIN, DIRICH, 0.0,0.0));   // for stopping rotation of body

    ///==================================================================
    /// Elasticity Equations
    ///==================================================================
    ///
    ///     K u = - Cp
    ///
    ///  K = int BT D B dV
    ///  C = int NT n N dS
    ///
    // setting up integrator
    //SAMG_Settings settings;
    //settings.Set_ncycle(10000);
    //settings.Set_napproach(2); // this is important because it sorts rhs vector [x1, y1, x2, y2, ..., xn, yn]
                               // which is needed for deformation simulations
    //SAMG_Solver samgSolver(&settings);
    EigenSolver eigen;

    PDE_Integrator<dim,Element> deformation (eigen);

    // Adding stiffness to the LHS list:
    NumIntegral_BT_D_B_dV<dim> stiffnessMatrix( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
    stiffnessMatrix.PlaneStress(plane_strain);
    deformation.Add( &stiffnessMatrix );

    // Adding surface tractions to the RHS (BoundaryIntegrals) list:
    //TODO: TEMPLATISE
    NumIntegral_PT_op_dS<dim> pressureTractions( model.Database(), "fluid pressure", "displacement" );
    pressureTractions.SubtractAccumulate();                       // because on RHS of equation, hence pressure traction term is  -int NT*n*N dx p^n  where n is outward pointing normal (pointing away from the element)
    //deformation.AddSplitBoundaryIntegral( &pressureTractions );


    InterFaceFractureVisitor<dim> aperture_visit(model.Database(), "displacement", "aperture", MIDDLE);



    printRangeOfVariable(model, "fluid pressure");
    printRangeOfVariable(model, "displacement");
    printRangeOfVariable(model, "Neumann stress");
    printRangeOfVariable(model, "aperture");


    ///==================================================================
    /// Displacement Solution & Verification
    ///==================================================================
    //Solve for displacements and move then nodes:
    deformation.IntegrateOver( model, model.Region("Model"));
    //Moving Nodes
    model.Accept(aperture_visit);
    model.MoveNodeCoordinatesBy("displacement");


    ///=========================================
    /// Outputing to VTU
    ///=========================================
    //Output variables:
    VTU_Interface<dim> vtu_elastic( model );
    vtu_elastic.OmitZeroInFileName(true);
    vtu_elastic.DeleteConnectivity();

    list<string> outputProps_elastic;
    outputProps_elastic.push_back("displacement");
    outputProps_elastic.push_back("fluid pressure");
    vtu_elastic.OutputDataToVTU( output_file, outputProps_elastic, "Model", static_cast<int>(1) );


    //Verify Aperture
    ///ANALYTICAL SOLUTION
    double a = 0.5 * myFracture.FractureLength(MIDDLE), max_percent_error(0.0), avg_percent_error(0.0);           // fracture half length // Maximum values obtained
    auto     F = [ P0, ym, pr, a] (double x) {return (4.0*P0*(1.0-pr*pr)/ym)*sqrt(a*a - x*x) ;} ;             //defining a lambda function to give analytical solution
    size_t ends = 0;
    std::vector<std::vector<double>> data(5);
    std::map<Point<dim>,Node<dim>* > NodeMapAfter = myFracture.NodeMap(MIDDLE);
    for (typename std::map<Point<dim>,Node<dim>*>::iterator it = NodeMapAfter.begin(); it != NodeMapAfter.end(); it++){
        assert( myFracture.DistanceFromTip(NodeMapAfter, it->second, RIGHT) <= 2.0*a);
        double x        =  a  -  myFracture.DistanceFromTip(NodeMapAfter, it->second, RIGHT) ;
        double w_num    = it->second->Read(model.Database().StorageKey("aperture"));
        double p_num    = it->second->Read(model.Database().StorageKey("fluid pressure"));
        double Sol_anal = F(x);
        double percent_error(0);

        std::cout << "\nDist from Origin: " << x << " Analytical w: " << Sol_anal << std::endl;
        std::cout << "Dist from Origin: " << x << " Numerical  w: " << w_num << std::endl;


        //Storing solutions in vector object
        data[0].push_back(x);
        data[1].push_back(Sol_anal);
        data[2].push_back(w_num);
        data[3].push_back(p_num);

        if (Sol_anal == 0.0){
            ends++;
            percent_error = 0.0;
            data[4].push_back(percent_error);
          }
        else{
            percent_error  = std::fabs( w_num  - Sol_anal ) / Sol_anal *100.0;
            data[4].push_back(percent_error);
            _equal(w_num, Sol_anal, 0.3 );
            avg_percent_error += percent_error;
            if (max_percent_error < percent_error )
              max_percent_error   = percent_error ;
          }
        std::cout << "Percent Error e = " << percent_error << std::endl;


      }
    assert(ends == 2);
    avg_percent_error = avg_percent_error / ( static_cast<double>(NodeMapAfter.size()) - static_cast<double>(ends) );
    //_equal(1.0, avg_percent_error < 10.0 ? 1.0:0.0, 0);
    //_equal( 1.0, max_percent_error < 15.0 ? 1.0 : 0.0 , 0);
    std::cout << "Average Percent error found for single node is: " << avg_percent_error << std::endl;
    std::cout << "Maximum Percent error found for single node is: " << max_percent_error << std::endl;
    _test( avg_percent_error < min_avg_error );


 /*   // Model dimensions:
    Point<dim> minPost, maxPost;
    model.MinMaxCoordinates( minPost, maxPost );
    const double l( maxPost[0]-minPost[0] ), h( maxPost[1]-minPost[1] );
    printModelDimensions( model, true );


    /// -------------------------------
    /// Setting up Fracture Configuration
    /// -------------------------------
    const std::string   split_fracture = "Split Fracture";
    model.ConvertBoundaryToSplitBoundary("FRACTURE", split_fracture);
    Fracture<dim> myFracture (model, model.SplitBoundary(split_fracture), quarterpoint, true );            //Creates interface objects in Fracture, and initialises & configures Lubrication region

    const std::string HydraulicFracture = myFracture.SplitBoundaryName();
    const std::string Lubrication       = myFracture.MidRegionName();

    /// --------------------------------
    /// Elasto-Lubrication Parameter Configuration
    /// --------------------------------
    // Elasticity Variables:
    model.CreateProperty( "displacement", "SI", VECTOR, NODE );
    model.CreateProperty( "Young's modulus", "SI", SCALAR, ELEMENT);
    model.CreateProperty( "Poisson's ratio", "SI", SCALAR, ELEMENT);
    // Fluid Flow variables
    model.CreateProperty( "fracture fluid pressure", "SI", SCALAR, NODE);                                       //makes use of updating PropertyDatabase Indexing
    model.CreateProperty("aperture", "SI", SCALAR, NODE, 1, 0.0, 1.0e+3) ;

    //Inputing Parameter Values
    model.InputPropertyValue("Young's modulus", ScalarVariable( PLAIN, ym ) );
    model.InputPropertyValue("Poisson's ratio", ScalarVariable( PLAIN, pr ) );


    /// --------------------------------
    /// Elasto-Lubrication Boundary Conditions
    /// --------------------------------
    //Elasticity Dirichlet
    Node<dim> *cornerNodeLeftBottom, *cornerNodeRightBottom;
    for (vector<Node<dim>*>::const_iterator node (model.Region("Model").NodesBegin()); node != model.Region("Model").NodesEnd(); ++node )
    {
        if ( (*node)->AtBoundary() == CNR1 ) cornerNodeLeftBottom = (*node);
        if ( (*node)->AtBoundary() == CNR2 ) cornerNodeRightBottom = (*node);
    }

    cornerNodeLeftBottom->Store(model.Database().StorageKey("displacement"),VectorVariable<dim> (DIRICH, DIRICH , 0., 0.));
    cornerNodeRightBottom->Store(model.Database().StorageKey("displacement"), VectorVariable<dim>(PLAIN, DIRICH, 0.,0.));


    /// --------------------------------
    /// Elasto-Lubrication Initial Conditions
    /// --------------------------------
    model.Region("Model").InputPropertyValue("fracture fluid pressure", ScalarVariable( PLAIN, P0),  COMPLETE);     //uniform fluid pressure accross fracture


    ///==================================================================
    /// Elasticity Equations
    ///==================================================================
    ///
    ///     K u = - Cp
    ///
    ///  K = int BT D B dV
    ///  C = int NT n N dS
    ///
    // setting up integrator
    SAMG_Settings settings;
    settings.Set_ncycle(10000); // Depending on the type of problem, it may take many iterations to converge for anisotropic cases.
    settings.Set_nxtyp(0);
    SAMG_Solver samgSolver(&settings);
    PDE_IntegratorExperimental<dim,Region> deformation (samgSolver);

    // Adding stiffness to the LHS list:
    NumIntegral_BT_C_B_dV<dim,Element<dim> > stiffnessMatrix( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
    stiffnessMatrix.OutOfPlane(oopc);
    deformation.Add( &stiffnessMatrix );

    NumIntegral_PT_op_dS<dim,InterFace<dim> > pressureTractions(model.Database(), "fracture fluid pressure", "displacement");
    pressureTractions.SubtractAccumulate();                       // because on RHS of equation, hence pressure traction term is  -int NT*n*N dx p^n  where n is outward pointing normal (pointing away from the element)
    deformation.AddBoundaryIntegrals( &pressureTractions );

    // Adding post processes (stress and strain calculation) to the integrator list:
    StressAndStrainOutput<dim> stressAndStrain( model, "Young's modulus", "Poisson's ratio", "displacement");
    stressAndStrain.OutOfPlane(oopc);
    deformation.AddPostProcess( &stressAndStrain );


    ///==================================================================
    /// Displacement Solution & Verification
    ///==================================================================
    //Solve for displacements and move then nodes:
    deformation.IntegrateOver( model, model.Region("Model"));
    //Moving Nodes
    myFracture.ImplicitlyCalculateAndStoreApertures("aperture", "displacement");
    model.MoveNodeCoordinatesBy("displacement");
    myFracture.AverageCoordinatesToMidPoint();



    ///ANALYTICAL SOLUTION
    double a = 0.5 * myFracture.FractureLength(), max_percent_error(0.0), avg_percent_error(0.0);           // fracture half length // Maximum values obtained
    auto     F = [ P0, ym, pr, a] (double x) {return (4.0*P0*(1.0-pr*pr)/ym)*sqrt(a*a - x*x) ;} ;             //defining a lambda function to give analytical solution
    size_t ends = 0;
    std::vector<std::vector<double>> data(5);
    std::map<Point<dim>,Node<dim>* > NodeMapAfter = myFracture.NodeMap();
    for (typename std::map<Point<dim>,Node<dim>*>::iterator it = NodeMapAfter.begin(); it != NodeMapAfter.end(); it++){
        assert( myFracture.DistanceFromTip(it->second, RIGHT) <= 2.0*a);
        double x        =  a  -  myFracture.DistanceFromTip(it->second, RIGHT) ;
        double w_num    = it->second->Read(model.Database().StorageKey("aperture"));
        double p_num    = it->second->Read(model.Database().StorageKey("fracture fluid pressure"));
        double Sol_anal = F(x);
        double percent_error(0);

        std::cout << "\nDist from Origin: " << x << " Analytical w: " << Sol_anal << std::endl;
        std::cout << "Dist from Origin: " << x << " Numerical  w: " << w_num << std::endl;


        //Storing solutions in vector object
        data[0].push_back(x);
        data[1].push_back(Sol_anal);
        data[2].push_back(w_num);
        data[3].push_back(p_num);

        if (w_num == 0.0){
            ends++;
            percent_error = 0.0;
            data[4].push_back(percent_error);
          }
        else{
            percent_error  = std::fabs( w_num  - Sol_anal ) / Sol_anal * 100.0;
            data[4].push_back(percent_error);
            _equal(percent_error, 0.0, 10.0 );
            avg_percent_error += percent_error;
            if (max_percent_error < percent_error )
              max_percent_error   = percent_error ;
          }
        std::cout << "Percent Error e = " << percent_error << std::endl;


      }
    assert(ends == 2);
    avg_percent_error = avg_percent_error / ( static_cast<double>(NodeMapAfter.size()) - static_cast<double>(ends) );
    //_equal(1.0, avg_percent_error < 10.0 ? 1.0:0.0, 0);
    //_equal( 1.0, max_percent_error < 15.0 ? 1.0 : 0.0 , 0);
    std::cout << "Average Percent error found for single node is: " << avg_percent_error << std::endl;
    std::cout << "Maximum Percent error found for single node is: " << max_percent_error << std::endl;






    ///=========================================
    /// Outputing to Matlab
    ///=========================================
    //Output variables:
    std::string meshname = prefix_;
    char buffer[20];
    meshname.copy(buffer, 13, 46);
    meshname = buffer;
    if (quarterpoint == true)
      meshname += "QPT";
    std::string results = "results" + meshname;
    if ( false ){
      WriteSolutionToFile(data, "Sneddon_" + meshname + ".dat");
      }


    ///=========================================
    /// Outputing to VTU
    ///=========================================
    VTU_Interface<dim> vtu_elastic( model );
    vtu_elastic.OmitZeroInFileName(true);
    stressAndStrain.SetComponentNamesForStressAndStrainArrays(vtu_elastic);
    vtu_elastic.DeleteConnectivity();

    list<string> outputProps_elastic;
    outputProps_elastic.push_back("displacement");
    outputProps_elastic.push_back("stress at element");
    outputProps_elastic.push_back("strain at element");
    outputProps_elastic.push_back("fracture fluid pressure");
    vtu_elastic.OutputDataToVTU( results, outputProps_elastic, "Model", static_cast<int>(1) );

*/
    return;
}














template<typename T>
void SneddonCrack_VVCase::WriteSolutionToFile(std::vector<std::vector<T>> data, std::string name) {
    std::ofstream SolutionFile (name);
    cout << "Saving file: " << name << " --- must be a .dat file" << endl;
    assert(data.size() >= 2);
    assert(data[0].size() == data[1].size());
    for (size_t ix = 0; ix < data[0].size(); ix++){
        for (size_t iy = 0; iy<data.size(); iy++){
            SolutionFile << data[iy][ix] << "      " ;
          }
        SolutionFile << std::endl;
    }

}



} // csmp
