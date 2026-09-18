#include "LinearElasticityB_VVCase.h"
#include "ANSYS_Model2D.h"
#include "Model.h"
#include "Boundary.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "StressesAndStrains2.h"
#include "ExtractTensorVariableComponent.h"

using namespace std;


namespace csmp
  {

  /* LinearElasticity Test Case B
  ================================
  Mesh:       2D Linear Triangles
  Model:      Beam 40x2 m
  Test:       Deflection
  BC:         const displ / const force
  Criterion:  analytical solution: deflection
  ================================
  */

  void LinearElasticityB_VVCase::run()
    {
      // some constants
      enum{DIM=2};
      const double forceInY( 5000. );
      const double relativeToleranceDeflection(1.0E3), relativeTolerancePoisson(0.01);
      const ScalarVariable zeroScalar( PLAIN, 0. );
      const VectorVariable<DIM> zeroVector( PLAIN, 0. );
      const VectorVariable<DIM> zeroVectorDirichlet( DIRICH, 0. );

      // establishing model & output facility
      string modelName( "BeamLin" );
      ANSYS_Model2D model( modelName.data(), "LinEl.txt" );
      VTU_Interface<DIM> vtu( model );
      vtu.OmitZeroInFileName(true);

      // model dimension
      Point<DIM> min, max;
      model.MinMaxCoordinates( min, max );
      const double length( max[0]-min[0] ), height( max[1]-min[1] );
      printModelDimensions( model, true );
      double volumePrior = model.Region("Model").Volume();
      cout <<"\nThe model has a volume of: "<< volumePrior <<" m^3 unloaded."<< endl;

      // output properties & initial output
      list<string> outputProps;
      outputProps.push_back( "mean stress" );
      outputProps.push_back( "stress-x" );
      outputProps.push_back( "force" );
      outputProps.push_back( "displacement" );
      model.InputPropertyValue( "mean stress", zeroScalar );
      vtu.OutputDataToVTU( "BeamUnloaded", outputProps, "Model", static_cast<int>(0) );

      // model configuration
      const VectorVariable<DIM> forceRight( DIRICH, DIRICH, 0., -forceInY );
      model.InputPropertyValue( "Poisson's ratio", makeScalar( PLAIN, 0.23 ) );
      model.InputPropertyValue( "Young's modulus", makeScalar( PLAIN, 1.0e+9 ) );
      model.InputPropertyValue( "displacement", zeroVector );
      model.InputPropertyValue( "force", zeroVector );
      model.Boundary( "LEFT" ).InputPropertyValue( "displacement", zeroVectorDirichlet );
      Index forceKey( model.Database().StorageKey("force") );
      Node<DIM>* upperRightCorner( *(model.Region("Model").NodesBegin()) );
      for( vector<Node<DIM>*>::const_iterator node( model.Region("Model").NodesBegin() ); node != model.Region("Model").NodesEnd(); ++node )
        if( (*node)->x() > upperRightCorner->x() || (*node)->y() > upperRightCorner->y() )
          upperRightCorner = (*node);
      upperRightCorner->Store( forceKey, forceRight );

      // setting up & solving linear elasticity fea problem
      #ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Settings settings;
      SAMG_Solver   solver(&settings);
      settings.Set_napproach(2);
      PDE_Integrator<DIM,Element> deformation( solver  );
      #else
      CSMP_DEFAULT_LINEAR_SOLVER solver;
      PDE_Integrator<DIM,Element> deformation( solver );
      #endif

      PT_op<DIM> bforces( model.Database(), "force", "displacement" );
      NumIntegral_BT_D_B_dV<DIM> stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
      deformation.Add( &stiffness );
      deformation.Add( &bforces );
      StressesAndStrains<DIM>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement", true, true );
      deformation.AddPostProcess( &postpro );
      model.Apply( deformation );

      // we extract stresses for later output
      ExtractTensorVariableComponent<DIM>  xstress( model.Database(), "stress", "stress-x",  0, 0 );
      model.Apply(xstress);

      // apply resulting displacement
      model.MoveNodeCoordinatesBy("displacement");

      // reset connectivity (displacement!) and output
      vtu.DeleteConnectivity();
      vtu.OutputDataToVTU( "BeamLoaded", outputProps, "Model", static_cast<int>(0) );

      // model dimension posteriori
      double volumePost = model.Region("Model").Volume();
      double minDisplacement(0.), maxDisplacement(0.);
      model.MinMaxOf( "displacement", minDisplacement, maxDisplacement );

      // analytic & io
      const double I( height*height*height*1./12. );
      const double E( model.Region( "Model" ).Average( "Young's modulus" ) );
      const double Wmax( forceInY * length * length * length / ( E*I*3 ) );
      cout << endl << endl;
      cout <<"\nThe force in y is: "<< forceInY << endl;
      cout <<"\nThe stiffness in x is: "<< E << endl;
      cout <<"\nThe moment of inertia in x is: "<< I << endl;
      cout <<"\nThe model should have a maximum deflection of in x of: "<< Wmax << " m" <<endl;
      cout <<"\nThe model has a volume of: "<< volumePost <<" m^3 under load."<< endl;

      // tests & io
      const double toleranceA( Wmax*relativeToleranceDeflection );
      const double toleranceB( volumePrior*relativeTolerancePoisson );
      cout <<"\nThe model has a maximum deflection in x of: "<< maxDisplacement << " m" << endl;
      cout << endl << endl;

      _equal( Wmax, maxDisplacement, toleranceA );
      _equal( volumePost, volumePrior, toleranceB );
   }

  } // csmp
