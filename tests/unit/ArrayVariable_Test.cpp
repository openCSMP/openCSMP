#include "ArrayVariable_Test.h"

#include "ArrayVariable.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {


  void ArrayVariable_Test::run()
    {
      ArrayVariable av0( 100, 10. );

      _test( av0[21] == 10. );
      av0(21) = 21.;
      _test( av0[21] == 21. );
      _test( av0.Flag() == ANY );
      av0.Flag(DIRICH);
      _test( av0.Flag() == DIRICH );
      _test( av0.Size() == 100 );
      av0.Resize( 120 );
      _test( av0.Size() == 120 );
      _test( av0[21] == 21. );
      av0.Resize( 80 );
      _test( av0.Size() == 80 );
      _test( av0[21] == 21. );
      ArrayVariable av1( 100, 10., DIRICH );
      _test( !(av1==av0) );


       ArrayVariable av2( 81, 11. );
       _test( av2.IsWithinRange( 10., 12. ) );
       _test( !av2.IsWithinRange( 11.2, 12. ) );
       _test( !av2.IsWithinRange( 10.2, 10.8 ) );

       PropertyDatabase<3> pdb("CSMP-variables-vsTestLocked.txt");
       ArrayVariable av3( "element array 2", pdb, 1.3, DIRICH );
       _test( av3.Size() == 22 );
       _test( av3.Flag() == DIRICH );
       ArrayVariable av4( 22, 1.3, DIRICH );
       _test( av3 == av4 );


       // bin IO
       ArrayVariable avBinO( 500, 999., ANY );
       _test( av3 == av4 );
       fstream fp ("ArrayVariableBinaryIO", ios::out | ios::binary);
       if ( !fp.is_open() ) {
         cerr <<"\nArrayVariable_Test: Bindary file could not be created.";
         _test(false);
         return;
         }
       avBinO.Out(fp);
	   fp.close();
	   
	   fp.open("ArrayVariableBinaryIO", ios::in | ios::binary);
	   if (!fp.is_open()) {
         cerr <<"\nArrayVariable_Test: Bindary file could not be opened.";
         _test(false);
         return;
         }
       ArrayVariable avBinI;
       avBinI.In(fp);
	   fp.close();
       _test( avBinI == avBinO );

    }

  } // csmp
