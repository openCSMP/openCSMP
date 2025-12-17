//
//  JaggedArray3D_Comparison_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 25/11/2024.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <random>

#include "JaggedArray3D_Comparison_Test.h"
#include "JaggedArray3D.h"
#include "DynamicArray3D.h"

using namespace std;

namespace csmp {

void JaggedArray3D_Comparison_Test::run()
 {
    // testing the basic functionality of JaggedArray
    // Initialize the jagged array
    JaggedArray3D<double> jaggedArray;

    // Insert elements into the jagged array
    jaggedArray( 0, 0, 0 ) =  10.0;    // Assign 10.0 to position (0, 0, 0)
    jaggedArray( 1, 2, 3 ) =  25.5;    // Assign 25.5 to position (1, 2, 3)
    jaggedArray( 5, 0, 2 ) =  -7.3;    // Assign -7.3 to position (5, 0, 2)
    jaggedArray( 2, 3, 4 ) = 100.0;    // Assign 100.0 to position (2, 3, 4)
    jaggedArray.set( 2, 3, 5, 200.0 ); // Assign 200.0 to position (2, 3, 5)
    
    _test( jaggedArray.slices()   == 6 );
    _test( jaggedArray.rows()     == 4 );
    _test( jaggedArray.max_cols() == 6 );
    auto dimensions = jaggedArray.max_kij();
    _test( dimensions == make_tuple(6,4,6) );

    // Access and modify an existing value
    auto key = make_tuple(1, 2, 3);  // Create a key for accessing (1, 2, 3)
    if ( jaggedArray.contains(key) ) {
        std::cout << "Value at {1, 2, 3}: " << jaggedArray(1,2,3) << std::endl;
        _test( jaggedArray.contains(key) );
        // modify the value
        jaggedArray(1,2,3) += 10.0;
        std::cout << "Updated value at {1, 2, 3}: " << jaggedArray(1,2,3) << std::endl;
        _test( jaggedArray(1,2,3) == 35.5 );
    } else {
        std::cout << "Key {1, 2, 3} not found!" << std::endl;
    }

    // Attempt to access a non-existent key
    key = {10, 10, 10};
    if ( isnan( jaggedArray(10,10,10) ) ) {
        std::cout << "Key {10, 10, 10} does not exist." << std::endl;
        _test( !jaggedArray.contains( make_tuple(10,10,10) ) );
    }

    // Iterate over all elements in the jagged array
    std::cout << "\nAll elements in the jagged array:\n";
    for (const auto& [Key, value] : jaggedArray) {
        auto [k, i, j] = Key;  // Decompose the tuple into k, i, j
        std::cout << "Value at {" << k << ", " << i << ", " << j << "} = " << value << std::endl;
    }

    // Remove an element
    jaggedArray.erase(2, 3, 4);
    std::cout << "\nAfter removing {2, 3, 4}:\n";
    for (const auto& [Key, value] : jaggedArray) {
        auto [x, y, z] = Key;
        std::cout << "Value at {" << x << ", " << y << ", " << z << "} = " << value << std::endl;
    }
    _test( !jaggedArray.contains( make_tuple(2,3,4) ) );
    
    // comparing its performance with other containers
    ComparePerformanceWithVectorOfVectors();
    
 } // end run
     

// Generate a 3D vector with the same dimensions and values as the Jagged3DArray
static vector<vector<vector<double>>> generateVector( int slices, int rows, int cols, double startValue )
 {
    vector<vector<vector<double>>> vec(slices, vector<vector<double>>(rows, vector<double>(cols)));
    double value = startValue;
    for ( int k = 0; k < slices; ++k ) {
        for ( int i = 0; i < rows; ++i ) {
            for (int j = 0; j < cols; ++j ) {
                vec[k][i][j] = value++;
            }
        }
    }
    return vec;
}



/**
Test Parameters:
	•	The test creates both Jagged3DArray and vector<vector<vector<double>>> structures with dimensions 6x4x3 and populates them with sequential values starting from 1.0.
	2.	Random Access:
	•	A loop performs 1,000,000 random accesses to elements in both structures.
	•	Random indices are generated using mt19937 and uniform_int_distribution.
	3.	Timing:
	•	The chrono::high_resolution_clock is used to measure the duration of the random access for both structures.
	4.	Volatile Keyword:
	•	The volatile keyword is used to prevent the compiler from optimizing away the loop by assuming the results are unused.
	5.	Output:
	•	The program prints the time taken (in microseconds) for both the vector and Jagged3DArray.
*/
void JaggedArray3D_Comparison_Test::ComparePerformanceWithVectorOfVectors()
 {
    const int32_t slices = 6, rows = 4, cols = 3;
    const double  startValue = 1.0;

    // Create the data structures
    DynamicArray3D<double> dyn( slices, rows, cols, startValue );
    auto                   vec = generateVector( slices, rows, cols, startValue );
    JaggedArray3D<double>  jaggedArray( vec );
    if ( verbose_ )
      jaggedArray.Out();

    // Random access test
    const int32_t numTests = 1000000;
    mt19937 rng(random_device{}());
    uniform_int_distribution<int32_t> sliceDist(0, slices - 1);
    uniform_int_distribution<int32_t> rowDist(0, rows - 1);
    uniform_int_distribution<int32_t> colDist(0, cols - 1);

    double value{0.};
    // Measure time for DynamicArray3D
    auto startDyn = chrono::high_resolution_clock::now();
    for ( int32_t i{0u}; i < numTests; ++i) {
          int32_t r = sliceDist(rng);
          int32_t rw = rowDist(rng);
          int32_t c = colDist(rng);
          value += dyn(r,rw,c); // volatile to prevent optimization
      }
    auto endDyn = chrono::high_resolution_clock::now();
    auto dynDuration = chrono::duration_cast<chrono::microseconds>(endDyn - startDyn);


    // Measure time for vector of vectors of vectors
    auto startVec = chrono::high_resolution_clock::now();
    for ( int32_t i{0u}; i < numTests; ++i) {
          int32_t r = sliceDist(rng);
          int32_t rw = rowDist(rng);
          int32_t c = colDist(rng);
          value += vec[r][rw][c]; // volatile to prevent optimization
      }
    auto endVec = chrono::high_resolution_clock::now();
    auto vecDuration = chrono::duration_cast<chrono::microseconds>(endVec - startVec);


    // Measure time for Jagged3DArray
    auto startJagged = chrono::high_resolution_clock::now();
    for ( int32_t i{0u}; i < numTests; ++i) {
          int32_t r = sliceDist(rng);
          int32_t rw = rowDist(rng);
          int32_t c = colDist(rng);
          value += jaggedArray(r, rw, c); // volatile to prevent optimization
      }
    auto endJagged = chrono::high_resolution_clock::now();
    auto jaggedDuration = chrono::duration_cast<chrono::microseconds>(endJagged - startJagged);

    // Print results
    cout << "\n"<<"DynamicArray3D access time: " << dynDuration.count() << " microseconds\n";
    cout << "Vector access time:         " << vecDuration.count() << " microseconds\n";
    cout << "Jagged3DArray access time:  " << jaggedDuration.count() << " microseconds\n";
    
    // just using the value so that the optimiser can not elimnate its computation
    cout <<"\n"<< value << endl;

} // end ComparePerformanceWithVectorOfVectors


} // end csmp
