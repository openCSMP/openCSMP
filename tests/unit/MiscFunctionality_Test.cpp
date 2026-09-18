#include "MiscFunctionality_Test.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;

namespace csmp {

	MiscFunctionality_Test::MiscFunctionality_Test()
    {
    }


MiscFunctionality_Test::~MiscFunctionality_Test()
  {
  }


void MiscFunctionality_Test::run()
  {
      // testing angle between edges
      double angle0 = angleBetweenEdges( make_pair( Point<2>(0.,0.),Point<2>(1.,0.) ), make_pair( Point<2>(1.,0.),Point<2>(2.,0.) ) );
      _equal( angle0, 0., 1.0e-3 );

      double angle45 = angleBetweenEdges( make_pair( Point<2>(0.,0.),Point<2>(1.,0.) ), make_pair( Point<2>(0.,0.),Point<2>(1.,1.) ) );
      _equal( angle45, 45., 1.0e-3 );

      double angle90 = angleBetweenEdges( make_pair( Point<2>(0.,0.),Point<2>(0.,1.) ), make_pair( Point<2>(0.,0.),Point<2>(1.,0.) ) );
      _equal( angle90, 90., 1.0e-3 );

      double angleX = angleBetweenEdges( make_pair( Point<2>(1.,1.),Point<2>(2.,0.9) ), make_pair( Point<2>(-1.,0.),Point<2>(1.,0.5) ) );
      _equal( angleX, 19.746, 1.0e-3 );
      
      // combinations
      vector<int64_t>          sequence{1,2,3,4,5};
      const int64_t            samples{2};
      deque<vector<int64_t> >  combinations;
      size_t n_combinations = createUniqueCombinations( sequence, samples, combinations );
      
      cout <<"\ncombination vector:\n";
      for ( auto i : combinations ) {
            
            for ( auto j : i )
               cout << j <<" ";
            cout << endl;
        }
          
      
  } // end run




} // end csmp
