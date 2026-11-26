// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "createUniqueCombinations_Test.h"
#include "CSMP_highLevelUtilities.h"
#include <algorithm>
#include <deque>
#include <map>
#include <set>
#include <vector>

using namespace std;

namespace csmp {

	createUniqueCombinations_Test::createUniqueCombinations_Test()
    {
    }


createUniqueCombinations_Test::~createUniqueCombinations_Test()
  {
  }


void createUniqueCombinations_Test::run()
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
      // ------------
      // 2 out of 5 numbers: 5!/(2! 3!) = 10 combinations
      {
        vector<int64_t>          sequence{1,2,3,4,5};
        const uint32_t           samples{2};
        deque<vector<int64_t> >  combinations;
        const size_t n_combinations = createUniqueCombinations( sequence, samples, combinations );

        _test( n_combinations == 10U );
        _test( n_combinations == combinations.size() );

        set<vector<int64_t> >  distinct_combinations;
        map<int64_t,size_t>    occurrences;
        bool right_size{true}, ascending_members_of_sequence{true};
        for ( const auto& combination : combinations ) {
             if ( combination.size() != samples ) right_size = false;
             for ( size_t i=0; i<combination.size(); ++i ) {
                  // members come from the sequence and are listed in ascending order, hence without repetition
                  if ( find( sequence.begin(), sequence.end(), combination[i] ) == sequence.end() ) ascending_members_of_sequence = false;
                  if ( i > 0 && combination[i] <= combination[i-1] ) ascending_members_of_sequence = false;
                  occurrences[ combination[i] ]++;
               }
             distinct_combinations.insert( combination );
          }
        _test( right_size );
        _test( ascending_members_of_sequence );
        // no combination occurs twice (combinations, not permutations: {1,2} and {2,1} are the same)
        _test( distinct_combinations.size() == n_combinations );
        // every number is combined with each of the 4 others exactly once
        _test( occurrences.size() == 5U );
        for ( const auto& number : occurrences ) _test( number.second == 4U );
        _test( distinct_combinations.count( vector<int64_t>{1,2} ) == 1U );
        _test( distinct_combinations.count( vector<int64_t>{4,5} ) == 1U );
        _test( distinct_combinations.count( vector<int64_t>{2,1} ) == 0U );

        if ( verbose_ ) {
             cout <<"\ncombination vector ("<< n_combinations <<" combinations):\n";
             for ( const auto& i : combinations ) {
                  for ( const auto& j : i )
                    cout << j <<" ";
                  cout << endl;
               }
          }
      }
      // 3 out of 5 numbers: also 10 combinations
      {
        vector<int64_t>          sequence{1,2,3,4,5};
        deque<vector<int64_t> >  combinations;
        _test( createUniqueCombinations( sequence, 3U, combinations ) == 10U );
        _test( combinations.size() == 10U );
        for ( const auto& combination : combinations ) _test( combination.size() == 3U );
      }
      // all numbers, or one number at a time
      {
        vector<size_t>          sequence{7,8,9};
        deque<vector<size_t> >  all, single;
        _test( createUniqueCombinations( sequence, 3U, all ) == 1U );
        _test( all.size() == 1U && all.front() == sequence );
        _test( createUniqueCombinations( sequence, 1U, single ) == 3U );
        _test( single.size() == 3U );
        for ( const auto& combination : single ) _test( combination.size() == 1U );
      }
      // unsorted input with repeated numbers is sorted and made unique first: {3,1,3,2} -> {1,2,3}
      {
        vector<int>          sequence{3,1,3,2};
        deque<vector<int> >  combinations;
        _test( createUniqueCombinations( sequence, 2U, combinations ) == 3U );
        _test( sequence == (vector<int>{1,2,3}) );
        const set<vector<int> > distinct_combinations( combinations.begin(), combinations.end() );
        _test( distinct_combinations == (set<vector<int> >{ {1,2}, {1,3}, {2,3} }) );
      }
      // results are appended to combinations that are already stored in the container
      {
        vector<int>          sequence{1,2,3};
        deque<vector<int> >  combinations;
        _test( createUniqueCombinations( sequence, 2U, combinations ) == 3U );
        _test( createUniqueCombinations( sequence, 2U, combinations ) == 6U );
      }
      // nothing to combine: empty input, or more samples requested than there are distinct numbers
      // (in the latter cases a message is written to the error stream)
      {
        vector<int>          empty_sequence, short_sequence{1,2}, repeated_number{7,7,7};
        deque<vector<int> >  combinations;
        _test( createUniqueCombinations( empty_sequence, 2U, combinations ) == 0U );
        _test( createUniqueCombinations( short_sequence, 3U, combinations ) == 0U );
        _test( createUniqueCombinations( repeated_number, 2U, combinations ) == 0U );
        _test( combinations.empty() );
      }

  } // end run




} // end csmp
