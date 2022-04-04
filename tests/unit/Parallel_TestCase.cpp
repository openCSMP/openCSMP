#include "Parallel_TestCase.h"
#include "ANSYS_Model3D.h"
#include "VTU_Interface.h"
#include "InterFace.h"
#include "SteadyStateDiffusor.h"
#include "TwoPhaseImplicitNodeCenteredFVTransport.h"
#include "BrooksCorey.h"

using namespace std;




namespace csmp
  {


  void Parallel_TestCase::run()
    {
      // working model
      ANSYS_Model3D model( "Overlap2D", "CSMP-transport-variables.txt", true );
      VTU_Interface<3> vtu( model );
      vtu.OmitZeroInFileName(true);

      // region refs
      /*
      Region<3>& modelRegion( model.Region( "Model" ) );
      Region<3>& leftRegion( model.Region( "LEFT" ) );
      Region<3>& rightRegion( model.Region( "RIGHT" ) );
      Region<3>& overlapRightRegion( model.Region( "RIGHT_OVERLAP" ) );
      Region<3>& overlapLeftRegion( model.Region( "LEFT_OVERLAP" ) );
      */
      Region<3>& boundary1( model.Region( "BOUNDARY_LEFT" ) );
      Region<3>& boundary2( model.Region( "BOUNDARY_RIGHT" ) );

      // create merged overlaps
      set<string> leftSide; leftSide.insert( "LEFT" ); leftSide.insert( "LEFT_OVERLAP" );
      set<string> rightSide; rightSide.insert( "RIGHT" ); rightSide.insert( "RIGHT_OVERLAP" );
      model.MergeRegions( leftSide, "LEFT_ENSEMBLE" );
      model.MergeRegions( rightSide, "RIGHT_ENSEMBLE" );
      /*
      Region<3>& leftEnsemble( model.Region( "LEFT_ENSEMBLE" ) );
      Region<3>& rightEnsemble( model.Region( "RIGHT_ENSEMBLE" ) );
      */

      // setting material and fluid properties
      model.InputPropertyValue( "permeability", makeScalar( PLAIN, 1.0E-13 ) );
      model.InputPropertyValue( "residual saturation oil", makeScalar( PLAIN, 0.1 ) );
      model.InputPropertyValue( "residual saturation water", makeScalar( PLAIN, 0.1 ) );
      model.InputPropertyValue( "viscosity oil", makeScalar( PLAIN, 1.E-3 ) );
      model.InputPropertyValue( "viscosity water", makeScalar( PLAIN, 1.E-3 ) );
      model.InputPropertyValue( "density oil", makeScalar( PLAIN, 1.E+3 ) );
      model.InputPropertyValue( "density water", makeScalar( PLAIN, 1.E+3 ) );
      model.InputPropertyValue( "porosity", makeScalar( PLAIN, 0.1 ) );
      model.InputPropertyValue( "entry pressure", makeScalar( PLAIN, 1.e-4 ) );
      model.InputPropertyValue( "brooks corey parameter", makeScalar( PLAIN, 2. ) );
      model.InputPropertyValue( "total system compressibility", makeScalar( PLAIN, 1.0E-9 ) );


      // initial conditions
      model.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
      model.InputPropertyValue( "saturation oil", makeScalar( PLAIN, 0.9 ) );
      model.InputPropertyValue( "saturation water", makeScalar( PLAIN, 0.1 ) );
      model.InputPropertyValue( "nodal fluid volume source", makeScalar( PLAIN, 0. ) );
      model.InputPropertyValue( "fluid volume source", makeScalar( PLAIN, 0.0 ) );

      // bc's
      boundary1.InputPropertyValue( "fluid pressure", makeScalar( DIRICH, 100.0E+5 ) );
      boundary1.InputPropertyValue( "saturation oil", makeScalar( DIRICH, 0.0 ) );
      boundary1.InputPropertyValue( "saturation water", makeScalar( DIRICH, 1.0 ) );
      boundary2.InputPropertyValue( "fluid pressure", makeScalar( DIRICH, 10.0E+5 ) );

      // pressure solver
      SteadyStateDiffusor<3,csmp::Region>  SSPS( model, "total mobility", "fluid pressure",
                                                                       "fluid volume source" );

      // transport algorithm
      // TwoPhaseImplicitNodeCenteredFVTransport<3,StencilProcessor>  advector( "Model", model );
      TwoPhaseImplicitNodeCenteredFVTransport<3,StencilProcessor>  advectorLeft( "LEFT_ENSEMBLE", model );
      TwoPhaseImplicitNodeCenteredFVTransport<3,StencilProcessor>  advectorRight( "RIGHT_ENSEMBLE", model );

      // saturation functions
      BrooksCorey<3> saturationFunctions( model.Database(), "permeability", "viscosity oil",
                                                             "viscosity water", "density oil", "density water",
                                                             "brooks corey parameter", "entry pressure", "saturation water",
                                                             "residual saturation water", "residual saturation oil" );

      // initial state and steady state pressure
      UpdateFlowProps( model, saturationFunctions );
      SSPS.ComputeSteadyState( model.Region("Model") );
      ComputeTotalVelocity( model );
      if ( verbose_ ) {
          vtu.OutputDataToVTU( "InitialPressureDistribution", "fluid pressure", "Model", static_cast<int>(0) );
          vtu.OutputDataToVTU( "InitialSaturationDistribution", "saturation oil", "Model", static_cast<int>(0) );
        }
      // extrapolating to visualize streamlines
      model.ExtrapolateCellToNodeProperty( "velocity", "nodal velocity" );
      if ( verbose_ ) vtu.OutputDataToVTU( "InitialVelocityField", "nodal velocity", "Model", static_cast<int>(0) );

      // simulation settings
      double model_time( 0. );
      const double TIME_INCREMENT( 5000. ), SIMULATION_DURATION( 9000000. );

      // transient loop
      //while( model_time < SIMULATION_DURATION )
      //  {
      //    // advecting and updating of saturations, recalculation of mobility,
      //    // solving for steady state pressure and computing velocity
      //    // incrementing model time, updating render window
      //    advector.TransportPhase( saturationFunctions, TIME_INCREMENT );
      //    UpdateFlowProps( model, saturationFunctions );
      //    model.Apply( SSPS );
      //    ComputeTotalVelocity( model );
      //    model_time += TIME_INCREMENT;
      //    vtu.OutputDataToVTU( "TransientSaturation", "saturation oil", "Model", (long)model_time );
      //    vtu.OutputDataToVTU( "TransientPressure", "fluid pressure", "Model", (long)model_time );
      //  }

      while( model_time < SIMULATION_DURATION )
        {
          // advecting and updating of saturations, recalculation of mobility,
          // solving for steady state pressure and computing velocity
          // incrementing model time, updating render window
          advectorLeft.TransportPhase( saturationFunctions, TIME_INCREMENT );
          advectorRight.TransportPhase( saturationFunctions, TIME_INCREMENT );
          //ImposeSaturation( model, leftEnsemble, rightEnsemble );
          UpdateFlowProps( model, saturationFunctions );
          model_time += TIME_INCREMENT;
          if ( verbose_ ) {
              vtu.OutputDataToVTU( "TransientSaturation", "saturation oil", "Model", static_cast<size_t>(model_time) );
              vtu.OutputDataToVTU( "TransientPressure", "fluid pressure", "Model", static_cast<size_t>(model_time));
            }
        }

    }



  void Parallel_TestCase::UpdateFlowProps( Model<3>& model, TwoPhaseModel<3>& saturationFunctions )
    {
    ScalarVariable  sc, mob_t;
    Index saturationOilKey( model.Database().StorageKey( "saturation oil" ) );
    Index saturationWaterKey( model.Database().StorageKey( "saturation water" ) );
    Index totalMobilityKey( model.Database().StorageKey( "total mobility" ) );

    csmp::Region<3>& mref( model.Region( "Model" ) );

    // updating the saturation of water
    const vector<Node<3>*>::const_iterator nodesEnd( mref.NodesEnd() );
    for ( vector<Node<3>*>::const_iterator it = mref.NodesBegin(); it != nodesEnd; ++it )
      {
      sc = 1. - (*it)->Read( saturationOilKey );
      (*it)->Store(  saturationWaterKey, sc );
      }

    // computing the total mobility
    const vector<Element<3>*>::const_iterator elementsEnd( mref.CellsEnd() );
    for ( vector<Element<3>*>::const_iterator it = mref.CellsBegin(); it != elementsEnd; ++it )
      {
      // setting up the relative permeability model
      saturationFunctions.Initialize( *(*it) );
      saturationFunctions.InitializeForBaryCenter( *(*it) );
      saturationFunctions.EffectiveSaturation();

      // total mobility
      mob_t = saturationFunctions.TotalMobility();
      (*it)->Store(  totalMobilityKey, mob_t );
      }
    } // updateSaturationsAndComputeTotalMobility



  void Parallel_TestCase::ComputeTotalVelocity( csmp::Model<3U>& model ) const
    {
      VectorVariable<3U>   velo;
      DenseMatrix<DM_MIN>  DERIV(3U,3U);
      Index velocityKey( model.Database().StorageKey( "velocity" ) );
      Index fluidPresssureKey( model.Database().StorageKey( "fluid pressure" ) );
      Index totalMobilityKey( model.Database().StorageKey( "total mobility" ) );

      csmp::Region<3>&  mref( model.Region( "Model" ) );

      const vector<Element<3U>*>::const_iterator elementsEnd = mref.CellsEnd();
      for ( vector<Element<3U>*>::const_iterator it = mref.CellsBegin(); it != elementsEnd; ++it )
        {
          // computing the total velocity: vt = -k (lt grad p )
          const double mob_t = (*it)->Read( totalMobilityKey );
          velo = 0.;
          (*it)->dN_AtBaryCenter( DERIV, 1U );
          for ( auto i = 0; i < (*it)->Nodes(); ++i )
            {
              double pf = (*it)->N(i)->Read( fluidPresssureKey );
              velo(0)  += pf  * -DERIV(0,i) * mob_t;
              velo(1)  += pf  * -DERIV(1,i) * mob_t;
              velo(2)  += pf  * -DERIV(2,i) * mob_t;
            }

          // storing the computed velocity
          (*it)->Store( velocityKey, velo );
        } // elements region
        
    } // end ComputTotalVelocity



  } // csmp
