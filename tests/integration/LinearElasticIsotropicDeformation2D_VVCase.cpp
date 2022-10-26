#include "LinearElasticIsotropicDeformation2D_VVCase.h"
#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
//#include "NumIntegral_BT_C_B_dV.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "StressesAndStrains.h"
//#include "StressAndStrainOutput.h"
#include "ExtractTensorVariableComponent.h"
#include "NumIntegral_PT_op_dS.h"
#include "ModelSubDomain.h"

using namespace std;


namespace csmp
  {

  LinearElasticIsotropicDeformation2D_VVCase::LinearElasticIsotropicDeformation2D_VVCase(const char* prefix)
    {

      this->setName("LinearElasticIsotropicDeformation2D_VVCase");

      prefix_=prefix; //THIS IS CHANGED LATER ANYWAY
    }

    /**

    Author:  M. Nejati -- Adapted by E.Pezzulli for Melbourne repository
    Linear Elastic: Plane stress/Strain problem under uniform traction
    =================================
    Domain:     Omega: w x h Rectangle (Provided meshes is for a 10x10 square plane)
    Mesh:       Square_10x10_LinearTrianglesCoarse                  (Quadrilaterials dont work bcause of incorrect Ansys mesh)
                Square_10x10_QuadraticTrianglesCoarse
    Var File:   LinearElasticIsotropicDeformation2D_VVCase-variables.txt
    Reg File:   LinearElasticIsotropicDeformation2D_VVCase-regions.txt
    Material:   Linear ELastic Isotropic Material: E=Young's Modulus, nu= Poisson's Ratio;
    Loading:    Uniaxial/Biaxial Tension/Compression on RIGHT and TOP Boundaries
    BC:         Zero x and y displacements on LEFT and BOTTOM Edges, rexpectively

    =================================

    Assume a rectangular plane of width and height of w, and h, subjected to uniform stresses S1 and S2 in x and y directions, respectively.
    If Ux and Uy are the solutions of displacements in x and y directions, respectively:
    BC: The displacement boundary conditions are Ux=0 on LEFT and Uy=0 on BOTTOM; The Neumann conditions are: Sx=S1 on RIGHT and Sy=S2 on TOP.
    Solution of stress and displacement at point (x,y) over the domain Omega:

    Plane Stress:
    Sx(x,y)=S1
    Sy(x,y)=S2
    Ux(x,y)=1/E * (S1-nuS2)x
    Uy(x,y)=1/E * (S2-nuS1)y

    Plain Strain:
    Sx(x,y)=S1
    Sy(x,y)=S2
    Ux(x,y)=1/E * (S1-nuS2-nu^2(S1+S2))x
    Uy(x,y)=1/E * (S2-nuS1-nu^2(S1+S2))y

    The following implemetation, together with the other elastic tests by MN, uses the capability to assign faces over the boundaries
    of the domain. In this case, all the 2D elements building the regions TOP, BOTTOM, LEFT and RIGHT are deleted, and equivalent
    faces are created instead. Using this functionality, one can apply boundary integrals only on the boundaries with known faces.
    An alternative method is to keep the elements, and apply the boundary interals as rhs terms. One can easily switch between these
    two methodologies below. However, the use of boundary faces has the following advantages:
    1) The lower dimension elements are deleted, and therefore they can not have any contribution to the LHS of the system. In the case of
    the presence of boundary elements, one has to make sure that they do not contribute to the stiffness matrix. This can be done, for example,
    by making Young's modulus zero on these boundary elements.
    2) When using RHS operators to add a boundary integral, an iteration happens over all the elements, including the domain higher dimension
    elements, to accumulate the boundary integrals. Although the contribution of domain elements can be forced to be zero (by making the operand
    zero for all domain elements), this iteration is very expensive compared to the one that just iterates over the boundary faces.
    */

  void LinearElasticIsotropicDeformation2D_VVCase::run(){


      //important - It is unclear if the setup of the first test affect the second.
      //Todo: checks tests are indeed independent. This may be due to how memory id allocated on the heap...
      //      For now, harder test (quadratic) is tested first.
      //Note: Inconsistent behaviour was observer when running the same executable. This was thought to be due a problem
      //with the calibraton of VData pfverts. After fixing that, consistency returned. But one day... it may come back...

      //Test 1 a) b)
      plane_strain_ = true;
      SetupLinearTriangles();
      Solve();

      /*
      SetupQuadraticTriangles();
      Solve();

      SetupLinearQuadrilaterals();
      Solve();

      //this doesnt work - bug in model setup
      SetupQuadraticQuadrilaterals();
      //Solve();

      //Test 2 a) b)
      plane_strain_ = false;
      SetupQuadraticTriangles();
      Solve();

      SetupLinearTriangles();
      Solve();
*/
      //All tests ran
  }

  void LinearElasticIsotropicDeformation2D_VVCase::Solve(){

      // Defining constants
      enum{DIM=2};
      const double s1(1.0e+6), s2(10.0e+6), ym (1.0e9), pr(0.3);

      // Model configuration:
      //Input file directory locations
      string input_dir    = "";
      string output_file  = "../Output/";
      output_file         += prefix_;
      //Variable inputs
      string vars_file     = input_dir + "LinearElasticIsotropicDeformation2D_VVCase-variables.txt";
      string regions_file  = input_dir + "LinearElasticIsotropicDeformation2D_VVCase";
      //mesh file
      string mesh_file     = input_dir + prefix_;

      // Model configuration:
      ANSYS_Model2D model( mesh_file.c_str(), regions_file.c_str(), vars_file.c_str(), false, true, true);                       // when using boundary faces

      //Testing correct construction of model
      Region<DIM>& MyModelRegion = model.Region("Model");

      /// ---------------------------------------------------------------------------------
      ///BEGIN TESTING OF MODEL CONSTRUCTION ------------------------------------------------
      /// ---------------------------------------------------------------------------------

      //Testing total amount of boundary elements
      if (isTriangle_){
        _test( MyModelRegion.PerimeterCells() == 40 ) ;        //Same for both linear and quadratic
        if (isLinear_ == true)
            _test( MyModelRegion.PerimeterNodes() == 40);         //Nodes on perimiter of model for Linear elm
        else
            _test( MyModelRegion.PerimeterNodes() == 80);         //Nodes on perim for Quadratic elmts
      } else {
        _test( MyModelRegion.PerimeterCells() == 48-4 ) ;
        if (isLinear_ == true)
            _test( MyModelRegion.PerimeterNodes() == 48);         //Nodes on perimiter of model for Linear elm
        else
            _test( MyModelRegion.PerimeterNodes() == 96);
      }


      //Testing each boundary face to have nodes on perimeter
      for (auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); it++){
          Boundary<DIM>& b_ref = it->second;
          double nx(0.0);
          double ny(0.0);
          if (b_ref.Name() == "LEFT")
              nx = -1.0;
          else if (b_ref.Name() == "BOTTOM")
              ny = -1.0;
          else if (b_ref.Name() == "RIGHT")
              nx = 1.0;
          else if (b_ref.Name() == "TOP")
              ny = 1.0;
          else
              throw csmp::Exception(ERROR, "LinearElasticityIsotropicDeformation2D_VVCase",
                                    "Model is not boxed shaped. Are you using Square_10x10_*****Coarse.*** mesh file? ");

          VectorVariable<DIM> nrml;
          for (auto b_elm_it = b_ref.CellsBegin(); b_elm_it != b_ref.CellsEnd(); b_elm_it++){
              //for each face object
              (*b_elm_it)->UnitNormal(nrml);        //obtain unit normal
              _equal( nrml[0],  nx, 0.001 );        //test x component
              _equal( nrml[1],  ny, 0.001 );        //test y component

              //testing order of elements
              if (isLinear_ == true)
                  _test( (*b_elm_it)->Nodes() == 2 );           //Testing faces are indeed linear bar elements
              else
                  _test( (*b_elm_it)->Nodes() == 3 );           // Testing faces are indeed quadratic bar elements

              //testing each node to be on its perimeter
              for ( size_t node = 0; node < (*b_elm_it)->Nodes() ; node++){
                  _test(MyModelRegion.IsPerimeterNode( (*b_elm_it)->N(node)) );
              }
          } //end of face iteration
      } // end of boundary iteration

      std::cout << "Model Tests finished - Now testing elasticity solver with traction terms" << std::endl;
      //END OF TESTING OF MODEL CONSTRUCTION ------------------------------------------------


      /// ---------------------------------------------------------------------------------
      /// BEGIN TESTING OF LINEAR ELASTICITY ------------------------------------------------
      /// ---------------------------------------------------------------------------------

      //checking unit normal direction
      Boundary<DIM> &b_ref = model.Boundary("TOP");
      for (auto it = b_ref.CellsBegin(); it != b_ref.CellsEnd(); it++){
          VectorVariable<DIM> nrml;
          (*it)->UnitNormal(nrml);
          std::cout << "\nUNIT NORMAL " << nrml[0] << " " << nrml[1] << std::endl;
          for ( size_t n = 0; n < (*it)->Nodes() ; n++)
              std::cout << "Node: " << n << " is " << model.Region("Model").IsPerimeterNode( (*it)->N(n) ) << std::endl;
       }

      // Model dimensions:
      Point<DIM> minPost (0,0), maxPost (0,0);
      model.MinMaxCoordinates( minPost, maxPost );
      const double l( maxPost[0]-minPost[0] ), h( maxPost[1]-minPost[1] );
      printModelDimensions( model, true );


      // Model Variables:
      model.InputPropertyValue("Young's modulus", ScalarVariable( PLAIN, ym ) );   // when using boundary faces
      model.InputPropertyValue("Poisson's ratio", ScalarVariable( PLAIN, pr ) );   // when using boundary faces
      model.InputPropertyValue("load", VectorVariable<DIM>(PLAIN, 0.0));

      // Model Boundary Conditions:
      Node<DIM> *cornerNodeLeftBottom, *cornerNodeRightTop;
      for (vector<Node<DIM>*>::const_iterator node (model.Region("Model").NodesBegin()); node != model.Region("Model").NodesEnd(); ++node )
      {
          if ( (*node)->AtBoundary() == CNR1 ) cornerNodeLeftBottom = (*node);
          if ( (*node)->AtBoundary() == CNR3 ) cornerNodeRightTop = (*node);        //needed for visualisation
      }

      model.Boundary("BOTTOM").InputPropertyValue("displacement", VectorVariable<DIM> (PLAIN , DIRICH, 0., 0.) );               // when using boundary faces
      model.Boundary("LEFT").InputPropertyValue(  "displacement", VectorVariable<DIM> (DIRICH, PLAIN , 0., 0.) );               // when using boundary faces
      cornerNodeLeftBottom->Store(model.Database().StorageKey("displacement"),VectorVariable<DIM> (DIRICH, DIRICH , 0., 0.));   // when using boundary faces
      model.Boundary("RIGHT").InputPropertyValue("Neumann stress", VectorVariable<DIM> (NEUMANN, NEUMANN, s1, 0.) );            // when using boundary faces
      model.Boundary("TOP").InputPropertyValue(  "Neumann stress", VectorVariable<DIM> (NEUMANN, NEUMANN, 0., s2) );            // when using boundary faces

      // setting up integrator
      //SAMG_Settings settings;
      //settings.Set_napproach(2); // needed apparently
      //settings.Set_ncycle(1000); // Depending on the type of problem, it may take many iterations to converge for anisotropic cases.
      //SAMG_Solver samgSolver(&settings);
      EigenSolver eigen;
      PDE_Integrator<DIM,Element> deformation (eigen);


      // Adding stiffness to the LHS list:
      NumIntegral_BT_D_B_dV<DIM > stiffnessMatrix( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
      stiffnessMatrix.PlaneStress(plane_strain_);  //E.P Note: Plane stress(true) makes the model Plane strain! Misleading. Should be relabelled PlaneStrain()
      deformation.Add( &stiffnessMatrix );

      // Adding surface tractions to the RHS (BoundaryIntegrals) list:
      NumIntegral_PT_op_dS<DIM> nodalTractions( model.Database(), "Neumann stress", "displacement" );
      deformation.AddBoundaryIntegral( &nodalTractions );                                                 // when using boundary faces

      StressesAndStrains<DIM>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement", plane_strain_, true );
      deformation.AddPostProcess( &postpro );

      //ExternalLoad<DIM> load_postpro(model.Database(), stiffnessMatrix, "load");
      //deformation.AddPostProcess( &load_postpro);

      // Solve for displacements and move then nodes:
      deformation.IntegrateOver( model, model.Region("Model"));                                             // when using boundary faces
      model.MoveNodeCoordinatesBy("displacement");

      // we extract stresses for later output
      ExtractTensorVariableComponent<DIM>  xstress( model.Database(), "stress", "stress-x",  0, 0 );
      model.Apply(xstress);
      ExtractTensorVariableComponent<DIM>  ystress( model.Database(), "stress", "stress-y",  1, 1 );
      model.Apply(ystress);
      ExtractTensorVariableComponent<DIM>  xystress( model.Database(), "stress", "stress-xy",  0, 1 );
      model.Apply(xystress);


      //Setting vtu file:
      VTU_Interface<DIM> vtu( model );
      vtu.OmitZeroInFileName(true);
      vtu.DeleteConnectivity();


      //Output variables:
      list<string> outputProps;
      outputProps.push_back( "displacement");
      outputProps.push_back( "stress" );
      outputProps.push_back( "stress-x" );
      outputProps.push_back( "stress-y" );
      outputProps.push_back( "stress-xy" );
      outputProps.push_back( "load");

      vtu.OutputDataToVTU( output_file, outputProps, "Model", static_cast<int>(1) );
      printRangeOfVariable(model, "stress-x");
      printRangeOfVariable(model, "stress-y");


      // Accuracy checks:
      double avg_stress_x = model.Region("Model").Average("stress-x");
      double avg_stress_y = model.Region("Model").Average("stress-y");
      cout << "\nAnalytical stress-x: " << -s1 << "\tAnalytical stress-y: " << -s2 << endl;
      cout << "\nNumerical  stress-x: " << avg_stress_x << "\tNumerical stress-y: " << avg_stress_y << endl;
      cout << "\nThe error of stress-x is: " << fabs(-s1-avg_stress_x)/fabs(avg_stress_x)*100 << "%" << endl;
      cout << "\nThe error of stress-y is: " << fabs(-s2-avg_stress_y)/fabs(avg_stress_y)*100 << "%" << endl;



      VectorVariable<DIM> dispN, dispA;
      cornerNodeRightTop->Read(model.Database().StorageKey("displacement"),dispN);
      if (!plane_strain_)
          dispA = VectorVariable<DIM> (PLAIN, PLAIN, l/ym*(s1-pr*s2),h/ym*(s2-pr*s1));
      else if (plane_strain_)
          dispA = VectorVariable<DIM> (PLAIN, PLAIN, l/ym*(s1-pr*s2-(pr*pr*(s1+s2))),h/ym*(s2-pr*s1-(pr*pr*(s1+s2))));
     double er = (dispA-dispN).Length()/dispA.Length()*100.;
     cout << "\nMax analytical displacement: " << dispA << endl;
     cout << "\nMax numerical displacement: " << dispN << endl;
     cout << "\nThe error of maximum displacement is: " << er << "%" << endl;

     _equal( dispA.Length(), dispN.Length(), dispA.Length()/1000. );
     _equal( -avg_stress_x/s1, 1.0, 0.05 );
     _equal( -avg_stress_y/s2, 1.0, 0.05 );

      return;
    } // End of Solve





  void LinearElasticIsotropicDeformation2D_VVCase::SetupLinearTriangles()
    {
      isLinear_   = true;
      isTriangle_ = true;
      prefix_     = "Square_10x10_LinearTrianglesCoarse";
    }

  void LinearElasticIsotropicDeformation2D_VVCase::SetupLinearQuadrilaterals(){
    isLinear_     = true;
    isTriangle_   = false;
    prefix_       = "Square_10x10_LinearQuadrilateralsCoarse";

  }

  void LinearElasticIsotropicDeformation2D_VVCase::SetupQuadraticTriangles(){
      isLinear_     = false;
      isTriangle_   = true;
      prefix_       = "Square_10x10_QuadraticTrianglesCoarse";
  }

  void LinearElasticIsotropicDeformation2D_VVCase::SetupQuadraticQuadrilaterals(){
    isLinear_     = false;
    isTriangle_   = false;
    prefix_       = "Square_10x10_QuadraticQuadrilateralsCoarse";
  }

  } // csmp
