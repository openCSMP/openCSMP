#include "VariableStorageSpeed_TestCase.h"

#include "ANSYS_Model3D.h"
#include "Region.h"

#include "Timer.hpp"

using namespace std;

namespace csmp{

  void VariableStorageSpeed_TestCase::run()
    {
      Timer timer;
      double ansys_build_time, native_build_time;

      // COMPARING ANSYS model generation with reading model from CSMP native binary file
      {
        timer.Start();
        ANSYS_Model3D ansys_model( "Clair", "Clair", "VariablesTutorial.txt", true ); // 799206 elmts, 145178 nodes
        cout << "\n\n\nTime taking to build model from ANSYS: " << (ansys_build_time=timer.StopClock()) << "\n\n\n";
        // saving model to disk
        ansys_model.OutputToBinaryFile("Clair");
      }
      // reading from binary file
      timer.Start();
      Model<3U> m0( string("Clair") );
      cout << "\n\n\nTime taking to build model from CSMP binary file set: " << (native_build_time=timer.StopClock()) << "\n\n\n";

      ScalarVariable sv( PLAIN, 1. );
      m0.InputPropertyValue( "nodal scalar", sv );
      m0.InputPropertyValue( "element scalar", sv );
      
      _test( native_build_time < ansys_build_time );
      cout <<"\nVariableStorageSpeed_TestCase::run: construction from CSMP native files is "<< ansys_build_time/native_build_time;
      cout <<" faster than via ANSYS_InterFace."<< endl;
      
      
      // Case 0: (base cases) reading vector of scalar variables (equivalent in length to Element vector
      // -----------------------------------------------------------------------------------------------
      // (no indexing overhead, multiplied by 5 because operation is repeated below)
      {
        vector<pair<double,int8_t>>  elmt_pairs( m0.Mesh().Elements() * 5, make_pair(4.,PLAIN) );
        timer.Start();
        pair<double,int8_t> val_flag_pair(2.,3);
        for ( auto& eit : elmt_pairs ) eit = val_flag_pair;
        cout << "\n\n"<<"Writing value-flag pairs: base case, range loop, writing vector: " << (native_build_time=timer.StopClock()) << endl;
      }
      // (no indexing overhead, multiplied by 5 because operation is repeated below)
      {
        vector<ScalarVariable>  elmt_scalars( m0.Mesh().Elements() * 5, makeScalar(PLAIN,4.) );
        timer.Start();
        for ( auto& eit : elmt_scalars ) eit = sv;
        cout << "\n\n"<<"Writing scalar variable: base case, range loop, writing vector of scalars: " << (native_build_time=timer.StopClock()) << endl;
      }
      // results
      // - shows that ScalarVariable is nearly as efficient as in-built pair
      // - version with int8_t flags is twice as fast as generic 'int' version
      

      // Case 1: STANDARD ITERATORS, element scalars and including aquiring Index keys and access to Region
      // --------------------------------------------------------------------------------------------------
      // ( recommended way of writing variables in CSMP
      auto sv1 = sv + 1.;
      timer.Start();
      const Index nKey0( m0.Database().StorageKey("nodal scalar") );
      const Index eKey0( m0.Database().StorageKey("element scalar") );
      Region<3U>& subdomain = m0.Region("Model");
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (range loop with index and region access): " << timer.StopClock() << endl;


      // Case 2: ITERATION, element scalars writing only
      // --------------------------------------------------------------------------------------------------
      sv1 = sv + 2.;
      timer.Start();
      const auto  elementsEnd = subdomain.CellsEnd();
      for( auto it=subdomain.CellsBegin(); it != elementsEnd; ++it ) (*it)->Store( eKey0, sv1 );
      for( auto it=subdomain.CellsBegin(); it != elementsEnd; ++it ) (*it)->Store( eKey0, sv1 );
      for( auto it=subdomain.CellsBegin(); it != elementsEnd; ++it ) (*it)->Store( eKey0, sv1 );
      for( auto it=subdomain.CellsBegin(); it != elementsEnd; ++it ) (*it)->Store( eKey0, sv1 );
      for( auto it=subdomain.CellsBegin(); it != elementsEnd; ++it ) (*it)->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (iteration and write only): " << timer.StopClock() << endl;


      // Case 3: RANGES, element scalars writing only
      // --------------------------------------------------------------------------------------------------
      sv1 = 5.;
      // RANGES
      timer.Start();
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      for( auto& it : subdomain.CellVector() ) it->Store( eKey0, sv1 );
      cout << "\n\n"<<"Writing scalar variable (range loop and write only): " << timer.StopClock() << endl;


      // Case 4: RANGES, reading (accumulating) scalar variable
      // --------------------------------------------------------------------------------------------------
      const Region<3U>& subdomainC = m0.Region("Model");
      double sum{0.};
      // RANGES
      timer.Start(); // const-ness only has minor speed-up associated with it (addition needed, else code gets eliminated)
      for( const auto& it : subdomainC.CellVector() ) sum += it->Read( eKey0 );
      for( const auto& it : subdomainC.CellVector() ) sum += it->Read( eKey0 );
      for( const auto& it : subdomainC.CellVector() ) sum += it->Read( eKey0 );
      for( const auto& it : subdomainC.CellVector() ) sum += it->Read( eKey0 );
      for( const auto& it : subdomainC.CellVector() ) sum += it->Read( eKey0 );
      cout << "\n\n"<<"Reading value part of scalar variable (range loop and read only): " << timer.StopClock() << endl;
      cout << sum << endl;


      // Case 5: RANGES, reading (accumulating) scalar variable
      // --------------------------------------------------------------------------------------------------
      sum = 0.;
      // RANGES
      timer.Start(); // const-ness only has minor speed-up associated with it (addition needed, else code gets eliminated)
      for( const auto& it : subdomainC.CellVector() ) { it->Read( eKey0, sv ); sum += sv(); }
      for( const auto& it : subdomainC.CellVector() ) { it->Read( eKey0, sv ); sum += sv(); }
      for( const auto& it : subdomainC.CellVector() ) { it->Read( eKey0, sv ); sum += sv(); }
      for( const auto& it : subdomainC.CellVector() ) { it->Read( eKey0, sv ); sum += sv(); }
      for( const auto& it : subdomainC.CellVector() ) { it->Read( eKey0, sv ); sum += sv(); }
      cout << "\n\n"<<"Reading scalar variable and summing values (range loop and read only): " << timer.StopClock() << endl;
      cout << sum << endl;
      
      // results
      // - read returning a double is not signifcantly faster that Read() returning a ScalarVariable into its argument

      // TODO: extra tests
      // TODO: - compared static dispatching with dynamic dispatching
      // TODO: - compare access of variables stored at FE integration points

      cin.get();
    }


  }
