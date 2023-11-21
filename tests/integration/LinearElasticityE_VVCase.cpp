#include "LinearElasticityE_VVCase.h"
#include "Model1D.h"
#include "VTU_Interface.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "SAMG_Exception.h"
#else
#include "LinearSolver.h"
#endif

#include "PDE_Integrator.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "StressesAndStrains2.h"
#include "LinearSolver.h"
#include "ExtractTensorVariableComponent.h"

using namespace std;


namespace csmp
  {

  /* Linear Elasticity Test Case E
  =================================
  Mesh:       Linear Bars
  Model:      Beam 3 m 1D
  Test:       Tensile strain
  BC:         const displ / const displ
  Criterion:  analytical solution: strain
  =================================
  */

  void LinearElasticityE_VVCase::run()
    {
      // some constants
      enum{DIM=1};
      const size_t ELEMENTS(3);
      const double LENGTH(1.);
      const double displacementInX( 0.099 );
      const double relativeToleranceStrain(0.01);
      const ScalarVariable zeroScalar( PLAIN, 0. );
      const VectorVariable<DIM> zeroVector( PLAIN, 0. );
      const VectorVariable<DIM> zeroVectorDirichlet( DIRICH, 0. );

      // establishing model & output facility
      string modelName( "BeamPatch" );
      Model1D<1U> model( modelName.data(), "LinEl.txt", LENGTH, ELEMENTS );

      // model configuration
      const VectorVariable<DIM> displacementRight( DIRICH, displacementInX );
      model.InputPropertyValue( "Poisson's ratio", makeScalar( PLAIN, 0.23 ) );
      model.InputPropertyValue( "Young's modulus", makeScalar( PLAIN, 1.0e+9 ) );
      model.InputPropertyValue( "displacement", zeroVector );
      model.InputPropertyValue( "force", zeroVector );

      // finding boundaries & assigning conditions
      model.InputBoundaryValue( LEFT, "displacement", zeroVectorDirichlet );
      model.InputBoundaryValue( RIGHT, "displacement", displacementRight );

      // setting up & solving linear elasticity fea problem
#ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Settings settings;
      SAMG_Solver   solver(&settings);
      settings.Set_napproach(2);
#else
      CSMP_DEFAULT_LINEAR_SOLVER solver;
#endif
      PDE_Integrator<DIM,Element> deformation( solver );

      PT_op<DIM> bforces( model.Database(), "force", "displacement" );
      NumIntegral_BT_D_B_dV<DIM> stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
      deformation.Add( &stiffness );
      deformation.Add( &bforces );
      deformation.IntegrateOver( model.Region("Model"), true );

      /* this should look like

      | 1  -1   0   0 |          | u1 |   | R1 |
      |-1   2  -1   0 |  EA /    | u2 |   | R2 |
      | 0  -1   2  -1 |  L       | u3 |   | R3 |
      | 0   0  -1   1 |          | u4 |   | R4 |

      after assembly and like

      |  1      0       0       0 |        | u1 |   |  0    |
      | -1EA/L  2EA/L  -1EA/L   0 |        | u2 |   |  0    |
      |  0     -1EA/L   2EA/L   0 |        | u3 |   |  EA/L |
      |  0      0       0       1 |        | u4 |   |  U4   |

      after incorporating BCs

      */
      deformation.OutputGlobals();

      // apply resulting displacement
      model.MoveNodeCoordinatesBy("displacement");

      // model dimension posteriori
      Point<DIM> minPost, maxPost;
      model.MinMaxCoordinates( minPost, maxPost );
      const double lengthPost( maxPost[0]-minPost[0] );

      // analytic & io
      const double E( model.Region("Model").Average("Young's modulus") );
      const double strainX( displacementInX/LENGTH );
      const double sigmaX( E*strainX );
      const double tensileForceInX( sigmaX );
      cout << endl << endl;
      cout <<"\nThe force in x is: "<< tensileForceInX << endl;
      cout <<"\nThe stiffness in x is: "<< E << endl;
      cout <<"\nThe model should have a strain in x of: "<< strainX << endl;
      cout <<"\nThe model should have a stress in x of: "<< sigmaX << " Pa" << endl;
      cout <<"\nThe model has a length of: "<< lengthPost <<" m under load."<< endl;

      // tests & io
      const double toleranceA( strainX*relativeToleranceStrain );
      const double strainXnum( (lengthPost-LENGTH)/LENGTH );
      cout <<"\nThe model has a strain in x of: "<< strainXnum << endl;
      cout << endl << endl;

      _equal( strainXnum, strainX, toleranceA );
  }

} // csmp
