#include "LinearElasticityA_VVCase.h"
#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "LUdcmp_Solver.h"
#include "PDE_Integrator.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "StressesAndStrains2.h"
#include "ExtractTensorVariableComponent.h"
#include "ANSYS_Model.h"

using namespace std;


namespace csmp
  {

  LinearElasticityA_VVCase::LinearElasticityA_VVCase(const char* prefix)
    {
      this->setName("LinearElasticityA_VVCase");
      prefix_=prefix;
    }

    /* Linear Elasticity Test Case A
    =================================
    Mesh:       Linear Triangles
    Model:      Beam 40x2 m 2D
    Test:       Tensile stress
    BC:         const displ / const force
    Criterion:  analytical solution: strain, stress, poisson deformation
    =================================
    */

  void LinearElasticityA_VVCase::run()
    {
      // some constants
      enum{DIM=2};
      const double tensileForceInX( 50000. );
      const double relativeToleranceStrain(2.0), relativeToleranceStress(.05), relativeTolerancePoisson(0.01);
      const ScalarVariable zeroScalar( PLAIN, 0. );
      const VectorVariable<DIM> zeroVector( PLAIN, 0. );
      const VectorVariable<DIM> zeroVectorDirichlet( DIRICH, 0. );
      
      // establishing model & output facility
      ANSYS_Model<2> model( prefix_, prefix_, "LinEl.txt" );
      VTU_Interface<DIM> vtu( model );
      vtu.OmitZeroInFileName(true);

      // model dimension
      Point<DIM> min, max;
      model.MinMaxCoordinates( min, max );
      const double length( max[0]-min[0] ), height( max[1]-min[1] );
      printModelDimensions( model, true );
      double64 volumePrior = model.Region("Model").Volume();
      cout <<"\nThe model has a volume of: "<< volumePrior <<" m^3 unloaded."<< endl;

      // output properties & initial output
      list<string> outputProps;
      outputProps.push_back( "mean stress" );
      outputProps.push_back( "stress-x" );
      outputProps.push_back( "force" );
      outputProps.push_back( "displacement" );
      model.InputPropertyValue( "mean stress", zeroScalar );
      

      // model configuration
      Boundary<DIM>& rightBoundary( model.Boundary( "RIGHT" ) );
      const VectorVariable<DIM> forceRight( DIRICH, DIRICH, tensileForceInX/rightBoundary.Nodes(), 0. );
      model.InputPropertyValue( "Poisson's ratio", makeScalar( PLAIN, 0.23 ) );
      model.InputPropertyValue( "Young's modulus", makeScalar( PLAIN, 1.0e+9 ) );
      model.InputPropertyValue( "displacement", zeroVector );
      model.InputPropertyValue( "force", zeroVector );
      model.Boundary( "LEFT" ).InputPropertyValue( "displacement", zeroVectorDirichlet );
      rightBoundary.InputPropertyValue( "force", forceRight );
      vtu.OutputDataToVTU( "BeamUnloaded", outputProps, "Model", static_cast<int>(0) );

      // setting up & solving linear elasticity fea problem
      #ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Settings settings;
      settings.Set_napproach(2); 
      PDE_Integrator<DIM,Region> deformation( new SAMG_Solver(&settings) );
      //PDE_Integrator<DIM,Region> deformation( new LUdcmp_Solver() );
      #else
      PDE_Integrator<DIM,Region> deformation( new CSMP_DEFAULT_LINEAR_SOLVER() );
      #endif

      PT_op<DIM,Element<DIM> > bforces( model.Database(), "force", "displacement" );
      NumIntegral_BT_D_B_dV<DIM,Element<DIM> > stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
      deformation.Add( &stiffness );
      deformation.Add( &bforces );
      StressesAndStrains<DIM>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement", true, true );
      deformation.AddPostProcess( &postpro );
      deformation.IntegrateOver( model.Region("Model") );

      // we extract stresses for later output
      ExtractTensorVariableComponent<DIM>  xstress( model.Database(), "stress", "stress-x",  0, 0 );
      model.Apply(xstress);

      // apply resulting displacement
      model.MoveNodeCoordinatesBy("displacement");

      // reset connectivity (displacement!) and output
      vtu.DeleteConnectivity();
      vtu.OutputDataToVTU( "BeamLoaded", outputProps, "Model", static_cast<int>(0) );

      // model dimension posteriori
      Point<DIM> minPost, maxPost;
      model.MinMaxCoordinates( minPost, maxPost );
      const double lengthPost( maxPost[0]-minPost[0] ), heightPost( maxPost[1]-minPost[1] );
      printModelDimensions( model, true );
      double64 volumePost = model.Region("Model").Volume();

      // analytic & io
      const double E( model.Region("Model").Average("Young's modulus") );
      const double strainX( tensileForceInX/E );
      const double sigmaX( tensileForceInX/heightPost );
      cout << endl << endl;
      cout <<"\nThe force in x is: "<< tensileForceInX << endl;
      cout <<"\nThe stiffness in x is: "<< E << endl;
      cout <<"\nThe cross sectional area in x is: "<< heightPost << endl;
      cout <<"\nThe model should have a strain in x of: "<< strainX << endl;
      cout <<"\nThe model should have a stress in x of: "<< sigmaX << " Pa" << endl;
      cout <<"\nThe model has a volume of: "<< volumePost <<" m^3 under load."<< endl;

      // tests & io
      const double toleranceA( strainX*relativeToleranceStrain );
      const double toleranceB( sigmaX*relativeToleranceStress );
      const double toleranceC( volumePrior*relativeTolerancePoisson );
      const double strainXnum( (lengthPost-length)/length );
      const double sigmaXnum( model.Region("Model").Average("stress-x") );
      cout <<"\nThe model has a strain in x of: "<< strainXnum << endl;
      cout <<"\nThe model has a stress in x of: "<< sigmaXnum << " Pa" << endl;
      cout << endl << endl;

      _equal( strainXnum, strainX, toleranceA );
      _equal( sigmaXnum, sigmaX, toleranceB );
      _equal( volumePost, volumePrior, toleranceC );
   }

  } // csmp
