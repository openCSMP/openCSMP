//
//  CSMP_VariableBenchmarking_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 1/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include <chrono>
#include <random>
#include "CSMP_mathUtilities.h"
#include "CSMP_VariableBenchmarking_Test.h"

using namespace std;

namespace csmp {

// LOOSER - cannot use move constructor
inline const ScalarVariable&  makeScalarConstRef( const VARIABLE_FLAG flag, const double64 val )
  {
     return ScalarVariable(flag,val);
  }

// implicit invocation of move constructor
inline ScalarVariable  makeScalar2( const VARIABLE_FLAG flag, const double64 val )
  {
     return ScalarVariable(flag,val);
  }

// not inlined version
ScalarVariable  makeScalarNotInlined( const VARIABLE_FLAG flag, const double64 val )
  {
     return ScalarVariable(flag,val);
  }




void VariableBenchmarking_Test::run()
  {
     // create test array with a million variables
     vector<ScalarVariable>  scalars(1000000);
     vector<double>          random_values(1000000);
     vector<VARIABLE_FLAG>   random_flags(1000000);
    
     // initialising the integer->flag sequence
     using number_engine     = default_random_engine;
     using int_distribution  = uniform_int_distribution<>;
     number_engine     intgen;
     int_distribution  uniform1{0,8}; // eight flags in VARIABLE_FLAG
     auto random_int = bind(uniform1,intgen);
     for ( auto it=random_flags.begin(); it!=random_flags.end(); ++it )
       (*it) = static_cast<VARIABLE_FLAG>(random_int());

     // initialising the value sequence
     using value_distribution  = uniform_real_distribution<double>;
     number_engine       doublegen;
     value_distribution  uniform2{0.,1.}; // generate floats between 0 and 1
     auto random_double = bind(uniform2,doublegen);
     for ( auto it=random_values.begin(); it!=random_values.end(); ++it )
       (*it) = random_double();

// test flags: OK
//cerr <<"\nVariableBenchmarking_Test: run: random sequence of variable flags:\n";
//for ( auto it=random_flags.begin(); it!=random_flags.end(); ++it )
//  cerr << parseStatus(*it) <<" ";
//cerr <<"\n";
//out(random_values);
    
     size_t running_index;
	   cerr <<"\nVariableBenchmark_Test::run: testing assignment of variable values and flags"<< endl;
 		 auto t0 = chrono::high_resolution_clock::now();
     // assigning values to scalars, vectors, tensors
     // CURRENT SCHEME - in CSMP
     running_index = 0U;
     for ( auto it=scalars.begin(); it!=scalars.end(); ++it )
       scalars[running_index] = makeScalar( random_flags[running_index], random_values[running_index] );
    
 		 auto t1 = chrono::high_resolution_clock::now();
	   cerr <<"\n\tCPU clock ticks used for current CSMP scheme " << chrono::duration_cast<chrono::nanoseconds>(t1-t0).count() << " nanoseconds." << endl;
    
     // CONST REFERENCE (rvalue) MOVE SCHEME
     running_index = 0U;
 		 auto t2 = chrono::high_resolution_clock::now();
     for ( auto it=scalars.begin(); it!=scalars.end(); ++it )
       scalars[running_index] = makeScalarConstRef( random_flags[running_index], random_values[running_index] );
    
 		 auto t3 = chrono::high_resolution_clock::now();
	   cerr <<"\n\tCPU clock ticks used for const value reference scheme " << chrono::duration_cast<chrono::nanoseconds>(t3-t2).count() << " nanoseconds." << endl;
   
    // there is only a point for the new scheme if it is faster !
     _test( chrono::duration_cast<chrono::nanoseconds>(t1-t0).count() > chrono::duration_cast<chrono::nanoseconds>(t3-t2).count() );
   

     // NOT INLINED SCHEME
     running_index = 0U;
 		 auto t4 = chrono::high_resolution_clock::now();
     for ( auto it=scalars.begin(); it!=scalars.end(); ++it )
       scalars[running_index] = makeScalarNotInlined( random_flags[running_index], random_values[running_index] );
    
 		 auto t5 = chrono::high_resolution_clock::now();
	   cerr <<"\n\tCPU clock ticks used for not inlined scheme " << chrono::duration_cast<chrono::nanoseconds>(t5-t4).count() << " nanoseconds.\n";
     cout << endl;
   
    // there is only a point for the new scheme if it is faster !
     _test( chrono::duration_cast<chrono::nanoseconds>(t5-t4).count() < chrono::duration_cast<chrono::nanoseconds>(t3-t2).count() );

  } // end run

} // end csmp
