#include "VariableStorageSpeed_TestCase.h"

#include "ANSYS_Model3D.h"

#include "Timer.hpp"

using namespace std;

namespace csmp{

  void VariableStorageSpeed_TestCase::run()
    {
      Timer timer;

      timer.Start();
      ANSYS_Model3D m0( "Clair", "Clair", "CSMP-variables-speedTest.txt", true );
      cout << "\n\n\nTime model building: " << timer.StopClocks() << "\n\n\n";

      Index nKey0( m0.Database().StorageKey("nodal variable") );
      Index eKey0( m0.Database().StorageKey("element variable") );
      ScalarVariable sv1( PLAIN, 1. );
      m0.InputPropertyValue( "nodal variable", sv1 );
      m0.InputPropertyValue( "element variable", sv1 );

      timer.Start();

      const vector<Element<3>*>::const_iterator elementsEnd( m0.Region("Model").ElementsEnd() );
      for( vector<Element<3>*>::const_iterator it( m0.Region("Model").ElementsBegin() ); it != elementsEnd; ++it )
        (*it)->Store( eKey0, sv1 );
      for( vector<Element<3>*>::const_iterator it( m0.Region("Model").ElementsBegin() ); it != elementsEnd; ++it )
        (*it)->Store( eKey0, sv1 );
      for( vector<Element<3>*>::const_iterator it( m0.Region("Model").ElementsBegin() ); it != elementsEnd; ++it )
        (*it)->Store( eKey0, sv1 );
      for( vector<Element<3>*>::const_iterator it( m0.Region("Model").ElementsBegin() ); it != elementsEnd; ++it )
        (*it)->Store( eKey0, sv1 );
      for( vector<Element<3>*>::const_iterator it( m0.Region("Model").ElementsBegin() ); it != elementsEnd; ++it )
        (*it)->Store( eKey0, sv1 );
      for( vector<Element<3>*>::const_iterator it( m0.Region("Model").ElementsBegin() ); it != elementsEnd; ++it )
        (*it)->Store( eKey0, sv1 );

      cout << "\n\n\nTime variable ops: " << timer.StopClocks() << "\n\n\n";

      cin.get();
    }


  }