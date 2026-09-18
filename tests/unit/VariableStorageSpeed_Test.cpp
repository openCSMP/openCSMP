#include "VariableStorageSpeed_Test.h"

#include "ANSYS_Model3D.h"
#include "Region.h"

#include "Timer.hpp"

using namespace std;

namespace csmp{

  void VariableStorageSpeed_Test::run()
    {
      Timer timer;
      double ansys_build_time, native_build_time;

      // COMPARING ANSYS model generation with reading model from CSMP native binary file
      {
        timer.Start();
        ANSYS_Model3D ansys_model( "HeuristicModel1coarse", "HeuristicModel1coarse", "VariablesTutorial.txt", true );
        cout << "\n\n\nTime taking to build model from ANSYS: " << (ansys_build_time=timer.StopClock()) << "\n\n\n";
        // saving model to disk
        ansys_model.OutputToBinaryFile("HeuristicModel1coarse");
      }
      // reading from binary file
      timer.Start();
      Model<3U> m0( string("HeuristicModel1coarse") );
      cout << "\n\n\nTime taking to build model from CSMP binary file set: " << (native_build_time=timer.StopClock()) << "\n\n\n";

      ScalarVariable sv( PLAIN, 1. );
      m0.InputPropertyValue( "nodal scalar", sv );
      m0.InputPropertyValue( "element scalar", sv );
      
      _test( native_build_time < ansys_build_time );
      cout <<"\nVariableStorageSpeed_TestCase::run: construction from CSMP native files is "<< ansys_build_time/native_build_time;
      cout <<" faster than via ANSYS_InterFace."<< endl;
      
      
      // Case 0: (base cases) reading vector of scalar variables (equivalent in length to Element vector
      // -----------------------------------------------------------------------------------------------
      {
        vector<pair<int8_t,double>>  elmt_pairs( m0.Mesh().Elements(), make_pair(PLAIN,3.) );
        timer.Start();
        pair<int8_t,double> val_flag_pair(PLAIN,7.);
        for ( auto& eit : elmt_pairs ) eit = val_flag_pair;
        cout << "\n\n"<<"Writing value-flag pairs: base case, range loop, writing vector: " << timer.StopClock() << endl;
      }
      {
        vector<ScalarVariable>  elmt_scalars( m0.Mesh().Elements(), makeScalar(PLAIN,4.) );
        timer.Start();
        for ( auto& eit : elmt_scalars ) eit = sv;
        cout << "\n\n"<<"Writing scalar variable: base case, range loop, writing vector of scalars: " << timer.StopClock() << endl;
      }
      // results
      // - shows that ScalarVariable is as efficient as in-built pair
      // - version with int8_t flags is twice as fast as generic 'int' version
      

      // Case 1: STANDARD ITERATORS, element scalars and including aquiring Index keys and access to Region
      // --------------------------------------------------------------------------------------------------
      // ( recommended way of writing variables in CSMP )
      auto sv1 = sv + 1.;
      timer.Start();
      const Index eKey0( m0.Database().StorageKey("element scalar") );
      Region<3U>& subdomain = m0.Region("Model");
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (range loop with index and region access): " << timer.StopClock() << endl;


      // Case 2: ITERATION, element scalars writing only
      // --------------------------------------------------------------------------------------------------
      sv1 = sv + 2.;
      timer.Start();
      const auto  elementsEnd = subdomain.CellsEnd();
      for( auto it=subdomain.CellsBegin(); it != elementsEnd; ++it ) (*it)->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (iteration and write only): " << timer.StopClock() << endl;


      // Case 3: RANGES, element scalars writing only
      // --------------------------------------------------------------------------------------------------
      sv1 = 5.;
      // RANGES
      timer.Start();
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (range loop and write): " << timer.StopClock() << endl;
      // Case 3b: RANGES, element scalars writing with templatized Index
      sv1 = 4.;
      const INDEX<SCALAR,ELEMENT> eKeyT( m0.Database().StorageKey("element scalar") );
      timer.Start();
      for( auto& it : subdomain.CellVector() ) it->Store( eKeyT, sv1 );
      cout << "\n\n"<<"Writing scalar variable using INDEX (range loop and write): " << timer.StopClock() << endl;


      // Case 4: RANGES, reading (accumulating) scalar variable
      // --------------------------------------------------------------------------------------------------
      const Region<3U>& subdomainC = m0.Region("Model");
      double sum{0.};
      // RANGES
      timer.Start(); // const-ness only has minor speed-up associated with it (addition needed, else code gets eliminated)
      for( const auto& it : subdomainC.CellVector() ) sum += it->Read( eKey0 );
      cout << "\n\n"<<"Reading value part of scalar variable (range loop and read): " << timer.StopClock() << endl;
      cout << sum << endl;


      // Case 5: RANGES, reading (accumulating) scalar variable
      // --------------------------------------------------------------------------------------------------
      sum = 0.;
      // RANGES
      timer.Start(); // const-ness only has minor speed-up associated with it (addition needed, else code gets eliminated)
      for( const auto& it : subdomainC.CellVector() ) { it->Read( eKey0, sv ); sum += sv(); }
      cout << "\n\n"<<"Reading scalar variable and summing values (range loop and read): " << timer.StopClock() << endl;
      cout << sum << endl;
      
      // results
      // - read returning a double is not signifcantly faster that Read() returning a ScalarVariable into its argument


      // Case 6: (base cases) writing vector of TENSOR variables (equivalent in length to Element vector
      // -----------------------------------------------------------------------------------------------
      TensorVariable<3U> ts(ANY,ANY,ANY,1.,0.,0.,0.,2.,0.,0.,0.,3.), ts2(ANY,ANY,ANY,4.,0.,0.,0.,8.,0.,0.,0.,12.);
      {
        vector<TensorVariable<3U>>  elmt_tensors( m0.Mesh().Nodes(), ts );
        timer.Start();
        // reading tensor, multiplying with tensor, storing tensor
        for ( auto& eit : elmt_tensors ) {
              eit *= ts2;
          }
        cout << "\n\n"<<"Reading, multiplying and storing TENSOR: base case, range loop: " << (native_build_time=timer.StopClock()) << endl;
      }

      // Case 7: STANDARD ITERATORS, nodal TENSORs and including aquiring Index keys and access to Region
      // --------------------------------------------------------------------------------------------------
      // ( recommended way of writing variables in CSMP )
      m0.InputPropertyValue( "nodal tensor", ts );
      timer.Start();
      const Index nKey0( m0.Database().StorageKey("nodal tensor") );
      const auto nodesEnd = subdomain.NodesEnd();
      for( auto nit=subdomain.NodesBegin(); nit != nodesEnd; ++nit ) {
           (*nit)->Read( nKey0, ts );
           ts *= ts2;
           (*nit)->Store( nKey0, ts );
        }
      cout << "\n\n"<<"Reading, multiplying and writing TENSOR variable (loop with index and region access): " << timer.StopClock() << endl;

      ScalarReadWriteWithSmallDataset();
      CompareIndexWithINDEX( m0 );
      TensorReadWithINDEXvsIndex( m0 );
      TestReadingArrayVariableVersusVector( m0 );
//      TestReadingArrayNew_vector_vs_reused_vector( m0 );
      
    } // end run
    
 
 
/*
   Uses minimalistic variable set to repeat tests from run() whether the size of the vector stored on the Element matters
*/
void VariableStorageSpeed_Test::ScalarReadWriteWithSmallDataset()
 {
      Timer timer;
      ANSYS_Model3D m0( "HeuristicModel1coarse", "HeuristicModel1coarse", "Minimum-variables.txt", true );

      ScalarVariable k( PLAIN, 1.0e-12 );
      m0.InputPropertyValue( "permeability", k );
      
      // Case 0: (base cases) reading vector of scalar variables (equivalent in length to Element vector
      // -----------------------------------------------------------------------------------------------
      // (no indexing overhead, multiplied by 5 because operation is repeated below)
      {
        vector<pair<int8_t,double>>  elmt_pairs( m0.Mesh().Elements(), make_pair(PLAIN,3.) );
        timer.Start();
        pair<int8_t,double> val_flag_pair(PLAIN,7.);
        for ( auto& eit : elmt_pairs ) eit = val_flag_pair;
        cout << "\n\n"<<"Writing value-flag pairs: base case, range loop, writing vector: " << timer.StopClock() << endl;
      }
      // (no indexing overhead, multiplied by 5 because operation is repeated below)
      {
        vector<ScalarVariable>  elmt_scalars( m0.Mesh().Elements(), makeScalar(PLAIN,4.) );
        timer.Start();
        for ( auto& eit : elmt_scalars ) eit = k;
        cout << "\n\n"<<"Writing scalar variable: base case, range loop, writing vector of scalars: " << timer.StopClock() << endl;
      }
      // results
      // - shows that ScalarVariable is nearly as efficient as in-built pair
      // - version with int8_t flags is twice as fast as generic 'int' version
      

      // Case 1: STANDARD ITERATORS, element scalars and including aquiring Index keys and access to Region
      // --------------------------------------------------------------------------------------------------
      // ( recommended way of writing variables in CSMP )
      auto sv1 = k + 1.;
      timer.Start();
      const Index eKey0( m0.Database().StorageKey("permeability") );
      Region<3U>& subdomain = m0.Region("Model");
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (range loop with index and region access): " << timer.StopClock() << endl;


      // Case 2: ITERATION, element scalars writing only
      // --------------------------------------------------------------------------------------------------
      sv1 = k + 2.;
      timer.Start();
      const auto  elementsEnd = subdomain.CellsEnd();
      for( auto it=subdomain.CellsBegin(); it != elementsEnd; ++it ) (*it)->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (iteration and write only): " << timer.StopClock() << endl;


      // Case 3: RANGES, element scalars writing only
      // --------------------------------------------------------------------------------------------------
      sv1 = 5.;
      // RANGES
      timer.Start();
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (range loop and write only): " << timer.StopClock() << endl;
      // Case 3b: RANGES, element scalars writing with templatized Index
      sv1 = 4.;
      const INDEX<SCALAR,ELEMENT> eKeyT( m0.Database().StorageKey("permeability") );
      timer.Start();
      for( auto& it : subdomain.CellVector() ) it->Store( eKeyT, sv1 );
      cout << "\n\n"<<"Writing scalar variable using INDEX (range loop and write only): " << timer.StopClock() << endl;


      // Case 4: RANGES, reading (accumulating) scalar variable
      // --------------------------------------------------------------------------------------------------
      const Region<3U>& subdomainC = m0.Region("Model");
      double sum{0.};
      // RANGES
      timer.Start(); // const-ness only has minor speed-up associated with it (addition needed, else code gets eliminated)
      for( const auto& it : subdomainC.CellVector() ) sum += it->Read( eKey0 );
      cout << "\n\n"<<"Reading value part of scalar variable (range loop and read only): " << timer.StopClock() << endl;
      cout << sum << endl;
      
 } // end ScalarReadWriteWithSmallDataset

 
 
 
 
/*
     Uses integration point variables as best case scenario for the application of templatized index variable.
*/
void VariableStorageSpeed_Test::CompareIndexWithINDEX( Model<3U>& model )
 {
     Timer timer;
     const Region<3U>& subdomainC = model.Region("Model");
       
      const csmp::Index                                   eKeyI( model.Database().StorageKey("element ip vector") );
      const csmp::INDEX<VECTOR,ELEMENT_INTEGRATION_POINT> eKeyT( model.Database().StorageKey("element ip vector") );
      
      model.InputPropertyValue( "element ip vector", makeVector(ANY,ANY,ANY,1.,2.,3.) );
      //printRangeOfVariable( model, "element ip vector", true );
 
      // Case 1: RANGES, reading (accumulating) vector variable from element integration points
      // --------------------------------------------------------------------------------------------------
      double             sum{0.};
      VectorVariable<3U> vc;
      // Index version
      timer.Start(); // const-ness only has minor speed-up associated with it (addition needed, else code gets eliminated)
      for ( const auto& it : subdomainC.CellVector() )
        for ( uint32_t i{0u}; i<it->IntegrationPoints(); ++i ) {
              it->Read( i, eKeyI, vc );
              sum += vc[0]; // + vc[1] + vc[3]; // making sure the variable gets used
          }
      cout << "\n\n"<<"Reading vector variable with Index: " << timer.StopClock() << endl;
      cout << sum << endl;

      // INDEX version
      sum = 0.;
      timer.Start(); // const-ness only has minor speed-up associated with it (addition needed, else code gets eliminated)
      for ( const auto& it : subdomainC.CellVector() )
        for ( uint32_t i{0u}; i<it->IntegrationPoints(); ++i ) {
              it->Read( i, eKeyT, vc );
              sum += vc[0]; //  + vc[1] + vc[3]; // making sure the variable gets used
          }
      cout << "\n\n"<<"Reading vector variable with INDEX: " << timer.StopClock() << endl;
      cout << sum << endl;
}




/*
    Learning: as expected, read/write with INDEX is faster than Index, but the difference is small.
*/
void VariableStorageSpeed_Test::TensorReadWithINDEXvsIndex( Model<3U>& model )
 {
     Timer timer;
     const Region<3U>& subdomain = model.Region("Model");
       
      const csmp::Index                                   eKeyI( model.Database().StorageKey("element ip vector") );
      const csmp::INDEX<VECTOR,ELEMENT_INTEGRATION_POINT> eKeyT( model.Database().StorageKey("element ip vector") );
      
      model.InputPropertyValue( "element ip vector", makeVector(ANY,ANY,ANY,1.,2.,3.) );
      printRangeOfVariable( model, "element ip vector", true );
 
      TensorVariable<3U> ts(ANY,ANY,ANY,1.,0.,0.,0.,2.,0.,0.,0.,3.), ts2(ANY,ANY,ANY,4.,0.,0.,0.,8.,0.,0.,0.,12.);
      model.InputPropertyValue( "nodal tensor", ts );

      // Case 1: iterators: reading tensor variable using Index
      // ------------------------------------------------------
      timer.Start();
      const Index nKeyI( model.Database().StorageKey("nodal tensor") );
      const auto nodesEnd = subdomain.NodesEnd();
      for( auto nit=subdomain.NodesBegin(); nit != nodesEnd; ++nit ) {
           (*nit)->Read( nKeyI, ts );
           ts(1,1) = 2;
           (*nit)->Store( nKeyI, ts );
        }
      cout << "\n\n"<<"Reading, adding to-, and writing TENSOR variable using Index: " << timer.StopClock() << endl;

      // Case 2: iterators: reading tensor variable using INDEX
      // ------------------------------------------------------
      timer.Start();
      const INDEX<TENSOR,NODE> nKeyT( model.Database().StorageKey("nodal tensor") );
      for( auto nit=subdomain.NodesBegin(); nit != nodesEnd; ++nit ) {
//           ts = (*nit)->Read( nKeyT );
           (*nit)->Read( nKeyT, ts );
           ts(1,1) = 2;
           (*nit)->Store( nKeyI, ts );
        }
      cout << "\n\n"<<"Reading, adding to-, and writing TENSOR variable using INDEX: " << timer.StopClock() << endl;

} // end TensorReadWithINDEXvsIndex



/*
    Learning: as expected, read/write with INDEX is faster than Index, but the difference is small.
*/
void VariableStorageSpeed_Test::TestReadingArrayVariableVersusVector( Model<3U>& model )
  {
     Timer timer;
     const Region<3U>& subdomain = model.Region("Model");
       
      ArrayVariable array(30), /* 30 is fixed in variables file */ ary;
      for ( unsigned int i{0}; i<30; ++i ) array(i) = static_cast<double>(i);
      ary = array;
      model.InputPropertyValue( "nodal array", array );
      printRangeOfVariable( model, "nodal array", true );

      // Case 1: iterators: reading array variable using Index
      // ------------------------------------------------------
      timer.Start();
      const Index nKeyI( model.Database().StorageKey("nodal array") );
      const auto nodesEnd = subdomain.NodesEnd();
      for( auto nit=subdomain.NodesBegin(); nit != nodesEnd; ++nit ) {
           (*nit)->Read( nKeyI, ary );
           (*nit)->Store( nKeyI, ary );
        }
      cout << "\n\n"<<"Reading, adding to-, and writing Array variable using Index: " << timer.StopClock() << endl;

      // Case 2: iterators: reading array variable using INDEX
      // ------------------------------------------------------
      timer.Start();
      const INDEX<ARRAY,NODE> nKeyT( model.Database().StorageKey("nodal array") );
      for( auto nit=subdomain.NodesBegin(); nit != nodesEnd; ++nit ) {
           (*nit)->Read( nKeyT, ary );
           (*nit)->Store( nKeyI, ary );
        }
      cout << "\n\n"<<"Reading, adding to-, and writing Array variable using INDEX: " << timer.StopClock() << endl;
      printRangeOfVariable( model, "nodal array", true );

      // Case 3: iterators: reading array variable vaue by value
      // -------------------------------------------------------
      // fastest!
      timer.Start();
      for( auto nit=subdomain.NodesBegin(); nit != nodesEnd; ++nit ) {
           for ( uint32_t i{0U}; i<nKeyI.dataDepth; ++i ) {
                double entry = (*nit)->ReadArrayEntry( nKeyI, i );
                (*nit)->StoreArrayEntry( nKeyI, entry, i );
             }
        }
      cout << "\n\n"<<"Reading & writing Array variable element by element: " << timer.StopClock() << endl;
      printRangeOfVariable( model, "nodal array", true );

} // end TestReadingArrayVariableVersusVector



} // end csmp
