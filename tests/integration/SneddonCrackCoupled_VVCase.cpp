#include"SneddonCrackCoupled_VVCase.h"
#include "SplitBoundaryInterface_Test.h"

#include "GlobalVerbose.h"

#include "CSMP_highLevelUtilities.h"
#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
//#include "NumIntegral_BT_C_B_dV.h"
//#include "StressAndStrainOutput.h"
#include "ExtractTensorVariableComponent.h"
#include "NumIntegral_PT_op_dS.h"
#include "ModelSubDomain.h"
#include "Fracture.h"
#include <functional>
#include <cmath>

//Interrellations for conductivity operator
#include "InterFaceFractureVisitor.h"

//Numerical Integrals Lubricaiton Equation
#include "NumIntegral_dNT_mixed_op_dN_dV.h"           // to interpolate the aperture
#include "NumIntegral_dNT_lhsop_dN_dV.h"                 // to interpolate the aperture cubed  (try with quadratic base functions...)
#include "NumIntegral_NT_lhsop_N_dV.h"                // Mass Matrix LHS
#include "NumIntegral_SetRHS_to_Zero.h"               // Zero right hand side
#include "NumIntegral_NT_rhsop_N_dV.h"                // Lumped Mass Matrix RHS or Source term with AccumulateLater()
#include "NumIntegral_PT_op_dV.h"
#include "PT_op.h"                                    //b force num int

//Coupled Case
#include "NumIntegral_BT_D_B_dV.h"                    //stiffness
#include "NumIntegral_PT_op_dS.h"                     //rhs traction terms (split)
//#include "NumIntegral_PT_n_N_dS.h"                    //lhs pressure traction term (splitboundary)
#include "NumIntegral_dNT_lhsop_dN_dV.h"                 //conductance (fluid pressure)
#include "PointSource_rhsop.h"                        //fluid rhs point source (on all nodes of element)
//#include "NumIntegral_NT_n_P_dS.h"                    // Aperture change LHS
//#include "NumIntegral_NT_n_rhs_P_dS.h"                // Aperture change RHS
//#include "NumIntegral_PT_n_N_dS.h"

//#include "RecoveryBasedOnDisplacement.h"

using namespace std;


namespace csmp
{



SneddonCrackCoupled_VVCase::SneddonCrackCoupled_VVCase(const char* prefix)
{
    this->setName("SneddonCrackCoupled_VVCase");
    prefix_=prefix;
}

/**

    Author:  E. Pezzulli
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




void SneddonCrackCoupled_VVCase::run()
{

  ///==========================================================================================
  /// Elasticity with Central Crack VV_Case
  ///==========================================================================================
  // Defining constants
    enum{dim=2};
    bool plane_strain = true;
    const double ym (1.0e9), pr(0.3), P0(1.0e7), mu(1.0e-3), /*s_unif(10000000.0), */ Q(0.01);
    //             s_right(s_unif), s_top(s_unif), s_bottom(10.0e6), s_left(10.0e6);


    // Model configuration:
    //Input file directory locations
    string input_dir         = "";
    string output_file       = "../Output/Workshop/";
    output_file              += prefix_;

    //Variable inputs
    string vars_file         = input_dir + "SneddonConvergenceTest-variables.txt";
    string regions_file      = input_dir + "SneddonConvergenceTest";
    //Possible meshes
    string  coarse_mesh      = "InternalBoundary_Test";
    string  coarse_mesh_quad = "InternalBoundary_Test_quadratic";
    string mesh_file         = input_dir + prefix_;

    bool binary(true), regions(true);
    ANSYS_Model2D model( mesh_file.c_str(), regions_file.c_str(), vars_file.c_str(),
                         binary, regions);    // Constractor for empty variables

    // SKM_FIX (and major improvement)
    const bool retain_elmts_as_intervening_elmts{ true };
    string sb_name = *(model.CreateSplitBoundaryFrom( "FRACTURE", retain_elmts_as_intervening_elmts ).first.begin() ) ; //relies on SB naming convention
    // string sb_reg  = *(model.InsertLowerDimensionalRegionsIntoSplitBoundaries(0).begin() );
    string sb_reg = "FRACTURE";
    // Note: this implies that the original connectivity of the lower-dim element mesh is not touched

    //testing region model has the new unique region
    //_test( model.Region("Model").Contains( *(model.Region("FRACTURE_SPLIT_BOUNDARY_REGION").ElementsBegin()) ) );

    /// -------------------------------
    /// Setting up Fracture Configuration
    /// -------------------------------
    Fracture<dim> myFracture (model, sb_name );            //Creates interface objects in Fracture, and initialises & configures Lubrication region
    //myFracture.TestSplitNodeAssignment(true);

    /// --------------------------------
    /// Elasto-Lubrication Parameter Configuration
    /// --------------------------------
    model.CreateProperty( "Neumann stress", "tau",   "SI",  VECTOR, FACE);
    model.CreateProperty( "viscosity",    "mu",    "SI",  SCALAR, NODE);
    model.CreateProperty( "conductance operator", "K", "SI", SCALAR, NODE);
    model.CreateProperty( "point source", "Q"  ,  "SI",  SCALAR, NODE);
    model.CreateProperty( "displacement old", "u_0" , "SI", VECTOR, NODE);

    //Initialising Parameter Values
    model.InputPropertyValue("Neumann stress",      VectorVariable<dim>(PLAIN, PLAIN, 0.0, 0.0));
    model.InputPropertyValue("displacement",        VectorVariable<dim>(PLAIN, PLAIN, 0.0, 0.0));
    model.InputPropertyValue("displacement old",    VectorVariable<dim>(PLAIN, PLAIN, 0.0, 0.0));
    model.InputPropertyValue("fluid pressure",      ScalarVariable( DIRICH, 0.0 ) );      //Initilising fluid pressure zero everywhere - so its ignored
    model.InputPropertyValue("aperture",            ScalarVariable(PLAIN, 0.01));
    model.InputPropertyValue("conductance operator",ScalarVariable( PLAIN, 0.0 ) ); //default initialisation
    model.InputPropertyValue("point source",        ScalarVariable(PLAIN,0.0));

    //Inputting material property values
    model.InputPropertyValue("Young's modulus", ScalarVariable( PLAIN, ym ) );
    model.InputPropertyValue("Poisson's ratio", ScalarVariable( PLAIN, pr ) );
    model.InputPropertyValue("viscosity",       ScalarVariable( PLAIN, mu ) );


    // Inputting fluid pressure values and Flags
    model.Region(sb_reg).InputPropertyValue("fluid pressure", ScalarVariable( DIRICH, P0),  COMPLETE);     //uniform fluid pressure accross fracture
    model.Region(sb_reg).InputPropertyValue("displacement", VectorVariable<dim>(DIRICH,DIRICH,0.0,0.0));
    /// --------------------------------
    /// Elasticity Boundary Conditions
    /// --------------------------------
    //Elasticity Dirichlet
    Node<dim> *cornerNodeLeftBottom{nullptr}, *cornerNodeRightBottom{nullptr};
    for (vector<Node<dim>*>::const_iterator node (model.Region("Model").NodesBegin()); node != model.Region("Model").NodesEnd(); ++node )
    {
        if ( (*node)->AtBoundary() == CNR1 ) cornerNodeLeftBottom = (*node);
        if ( (*node)->AtBoundary() == CNR2 ) cornerNodeRightBottom = (*node);
    }

    cornerNodeLeftBottom->Store(model.Database().StorageKey("displacement"),VectorVariable<dim> (DIRICH, DIRICH , 0., 0.));   // when using boundary faces
    cornerNodeRightBottom->Store(model.Database().StorageKey("displacement"), VectorVariable<dim>(PLAIN, DIRICH, 0.0,0.0));   // for stopping rotation of body


    /// --------------------------------
    /// Fluid Flow Boundary Conditions
    /// --------------------------------
    auto frac_elms = myFracture.ElementMap(MIDDLE);
    auto e_it = frac_elms.begin();
    advance( e_it, frac_elms.size()/2);
    double parents = static_cast<double>( e_it->second->N(0)->Parents()); //make sures we only apply Q and not Q*parents
    double nodes_per_element  = static_cast<double>(e_it->second->Nodes() );
    e_it->second->N(0)->Store(model.Database().StorageKey("point source"), ScalarVariable( NEUMANN, Q/parents*nodes_per_element  ));
    cout << "\nFlux " << Q/parents << " applied to Nodes with Coord" << e_it->second->N(0)->Coordinate() << std::endl;



    ///==================================================================
    /// Elasticity Equations
    ///==================================================================
    ///
    ///     K u + Cp = 0
    ///
    ///  K = int BT D B dV
    ///  C = int NT n N dS
    ///
    // setting up integrator
    //SAMG_Settings settings;
    //myFracture.SetSolverSettings(settings); //setting HF specific Solver settings
    //settings.Set_ncycle(5000);
    //settings.Set_napproach(2); // this is important because it sorts rhs vector [x1, y1, x2, y2, ..., xn, yn]
                               // which is needed for deformation simulations
    //SAMG_Solver samgSolver(&settings);
    CSMP_DEFAULT_LINEAR_SOLVER  eigen;
    PDE_Integrator<dim,Element> Coupled_HM (eigen);

    // Adding stiffness to the LHS list:
    NumIntegral_BT_D_B_dV<dim> stiffnessMatrix( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
    stiffnessMatrix.PlaneStress(plane_strain);
    Coupled_HM.Add( &stiffnessMatrix );

    //Adding coupled traction term LHS
    //NumIntegral_PT_n_N_dS<dim> Coupled_pressure_tractions(model.Database(), "fluid pressure", "displacement");
    //Coupled_HM.AddSplitBoundaryIntegral( &Coupled_pressure_tractions);

    //Adding surface tractions to the RHS (BoundaryIntegrals) list: (like zeroing right hand side)
    NumIntegral_SetRHS_to_Zero<dim> zero_load_rhs( model.Database(), "displacement");
    Coupled_HM.Add( &zero_load_rhs);

    //NumIntegral_PT_op_dS_Experimental<dim >  nodalTractions( model.Database(), "Neumann stress", "displacement" );
    //Coupled_HM.AddBoundaryIntegral( &nodalTractions );



    ///-------------------------------------------------------------------------------------------
    /// Lubrication Equation
    ///-------------------------------------------------------------------------------------------
    /// Steady State with Implicit Boundary Condition
    ///
    ///       dt* H p + C u  = Q + C u^(n-1)
    ///
    /// -----------------------------------------
    ///   H   =   int dNT * (w^3/(12*mu)) * dN dx
    ///   C   =   int NT * n * P dx
    ///
    //Interelations
    InterFaceFractureVisitor<dim> aperture_and_conductivity(model.Database(), "displacement", "aperture", MIDDLE, "conductivity", "viscosity", 1e-4);
    model.Accept(aperture_and_conductivity);

    //Conductance operator LHS
    NumIntegral_dNT_lhsop_dN_dV<dim> Conductance(model.Database(), "conductance operator", "fluid pressure", "fluid pressure");
    Coupled_HM.Add(&Conductance);

    //NumIntegral_NT_n_P_dS<dim> New_aperture(model.Database(), "displacement", "fluid pressure");
    //Coupled_HM.AddSplitBoundaryIntegral( &New_aperture);

    //NumIntegral_NT_n_rhs_P_dS<dim> Old_aperture(model.Database(), "displacement old", "fluid pressure");
    //Coupled_HM.AddSplitBoundaryIntegral( &Old_aperture );

    //Point source
    PointSource_rhsop<dim> point_source( model.Database(), "point source", "fluid pressure");
    Coupled_HM.Add( &point_source );


    /// checking Initialisation of properties
    printRangeOfVariable(model, "fluid pressure");
    printRangeOfVariable(model, "displacement");
    printRangeOfVariable(model, "displacement old");
    printRangeOfVariable(model, "Neumann stress");
    printRangeOfVariable(model, "aperture");
    printRangeOfVariable(model, "conductance operator");
    printRangeOfVariable(model, "point source");

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

    ///==================================================================
    /// Displacement Solution & Verification
    ///==================================================================
    //setting time steps
    double dt = 1.0;
    double t  = 0.0;
    double T  = 1.0;
    while ( t < T){
      t += dt;

      //Updating old displacements for rhs volume change
      model.CopyReplace("displacement", "displacement old");

      /// SOLVE ------------
      Coupled_HM.IntegrateOver( model, model.Region("Model"), false);

      //Updating conductivity
      model.Accept(aperture_and_conductivity);

      //Updating Volume of aperture
      double Volume = myFracture.Volume();

      //Outputing
      myFracture.SolOut(MIDDLE, "aperture", "fluid pressure");

      std::cout << "Fracture Volume and Flux = " << Volume << " " << Q*t << std::endl;
      vtu_elastic.OutputDataToVTU( output_file, outputProps_elastic, "Model", static_cast<int>(t) );
    }

    //Visualising node movement
    model.MoveNodeCoordinatesBy("displacement");
    // Moving mid region
    vtu_elastic.OutputDataToVTU( output_file, outputProps_elastic, "Model", static_cast<int>(t) );

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
        std::cout << "Dist from Origin: " << x << " Numerical  p: " << p_num << std::endl;

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
            percent_error  = std::fabs( w_num  - Sol_anal ) / Sol_anal * 100.0;
            data[4].push_back(percent_error);
            _equal(w_num,Sol_anal, 0.3 );
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
    _test( avg_percent_error < min_avg_error);
    _test(max_percent_error < 20.0); //too high


    std::cout << "Finished test SneddonCrackCoupled" << std::endl;

    return;


}



} // csmp
