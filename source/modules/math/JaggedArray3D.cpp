//
//  JaggedArray example
//  Open CSMP++
//
//  Created by Stephan Matthai on 24/11/2024.
//

#include <cassert>
#include "JaggedArray3D.h"

using namespace std;

namespace csmp {


/** Construction from a nested set  of initialiser lists
 
    Example:  Initialize jagged array with different column sizes for each row and rock
    
    @code
    Jagged3DArray jaggedArray = {
                                    {   // depth 0
                                        {1.1, 1.2},    // Row 0
                                        {2.1}           // Row 1
                                    },
                                    {   // depth 1
                                        {3.1, 3.2, 3.3}, // Row 0
                                        {4.1}            // Row 1
                                    },
                                    {   // depth 2
                                        {}               // Empty row
                                    }
                                };
    c @endcode
*/
template<class T>
JaggedArray3D<T>::JaggedArray3D( initializer_list<initializer_list<initializer_list<T>>> init )
{
  uint32_t depth = 0; // first dimension
  for ( const auto& slice : init ) {
        uint32_t row = 0; // second dimension
        for ( const auto& row_data : slice ) {
            uint32_t col = 0; // third dimension
            for ( double value : row_data ) {
                  // Assign the value to the map
                  data_[{depth, row, col}] = value;
                  ++col;
              }
            ++row;
        }
        ++depth;
    }
}


    /// construct from a monster of vector of vector of vectors
template<class T>
JaggedArray3D<T>::JaggedArray3D( const vector<std::vector<vector<T>>>& vec )
 {
    assert( !vec.empty() );
    for ( uint32_t k{0u}; k<vec.size(); ++k )
      for ( uint32_t i{0u}; i<vec[k].size(); ++i )
        for ( uint32_t j{0u}; j<vec[k][i].size(); ++j )
          set( i, j, k, vec[k][i][j] );
 }
 


   // DIAGNOSTICS
template<class T>
uint32_t JaggedArray3D<T>::slices() const {
     uint32_t slice{ 0u };
     for ( const auto& index : data_ )
       slice = max( slice, get<0>(index.first) );
     return slice + 1u;
  }
  
  
template<class T>
uint32_t JaggedArray3D<T>::rows() const {
     uint32_t row{ 0u };
     for ( const auto& index : data_ )
       row = max( row, get<1>(index.first) );
     return row + 1u;
  }
  
  
template<class T>
uint32_t JaggedArray3D<T>::max_cols() const {
     uint32_t col{ 0u };
     for ( const auto& index : data_ )
       col = max( col, get<2>(index.first) );
     return col + 1u;
  }


template<class T>
std::tuple<uint32_t, uint32_t, uint32_t> JaggedArray3D<T>::max_kij() const {
     uint32_t slice{ 0u };
     uint32_t row{ 0u };
     uint32_t col{ 0u };
     for ( const auto& index : data_ ) {
         slice = max( slice, get<0>(index.first) );
         row   = max( row, get<1>(index.first) );
         col   = max( col, get<2>(index.first) );
       }
     return make_tuple( slice+1u, row+1u, col+1u );
  }




static void example_JaggedArray()
 {
    // Initialize the jagged array
    JaggedArray3D<double> jaggedArray;

    // Insert elements into the jagged array
    jaggedArray( 0, 0, 0 ) =  10.0;    // Assign 10.0 to position (0, 0, 0)
    jaggedArray( 1, 2, 3 ) =  25.5;    // Assign 25.5 to position (1, 2, 3)
    jaggedArray( 5, 0, 2 ) =  -7.3;    // Assign -7.3 to position (5, 0, 2)
    jaggedArray( 2, 3, 4 ) = 100.0;    // Assign 100.0 to position (2, 3, 4)
    jaggedArray.set( 2, 3, 5, 200.0 ); // Assign 200.0 to position (2, 3, 5)

    // Access and modify an existing value
    auto key = make_tuple(1, 2, 3);  // Create a key for accessing (1, 2, 3)
    if ( jaggedArray.contains(key) ) {
        std::cout << "Value at {1, 2, 3}: " << jaggedArray(1,2,3) << std::endl;

        // Modify the value
        jaggedArray(1,2,3) += 10.0;
        std::cout << "Updated value at {1, 2, 3}: " << jaggedArray(1,2,3) << std::endl;
    } else {
        std::cout << "Key {1, 2, 3} not found!" << std::endl;
    }

    // Attempt to access a non-existent key
    key = {10, 10, 10};
    if ( isnan( jaggedArray(10,10,10) ) ) {
        std::cout << "Key {10, 10, 10} does not exist." << std::endl;
    }

    // Iterate over all elements in the jagged array
    std::cout << "\nAll elements in the jagged array:\n";
    for (const auto& [Key, value] : jaggedArray) {
        auto [x, y, z] = Key;  // Decompose the tuple into x, y, z
        std::cout << "Value at {" << x << ", " << y << ", " << z << "} = " << value << std::endl;
    }

    // Remove an element
    jaggedArray.erase(2, 3, 4);
    std::cout << "\nAfter removing {2, 3, 4}:\n";
    for (const auto& [Key, value] : jaggedArray) {
        auto [x, y, z] = Key;
        std::cout << "Value at {" << x << ", " << y << ", " << z << "} = " << value << std::endl;
    }

} // end jaggedArray example


template class JaggedArray3D<uint32_t>;
template class JaggedArray3D<float>;
template class JaggedArray3D<double>;

} // end csmp
