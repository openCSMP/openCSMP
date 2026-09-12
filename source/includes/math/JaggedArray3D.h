//
//  JaggedArray3D.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 25/11/2024.
//

#ifndef CSMP_JAGGED_ARRAY_3D_H
#define CSMP_JAGGED_ARRAY_3D_H

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <functional>
#include <cstdint>

namespace csmp {

/** custom hash for std::tuple<int, int, int>
 
    The STL uses just uses the function
    @code
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        std::size_t h1 = std::hash<int>{}(std::get<0>(t)); // Hash of the first element
        std::size_t h2 = std::hash<int>{}(std::get<1>(t)); // Hash of the second element
        std::size_t h3 = std::hash<int>{}(std::get<2>(t)); // Hash of the third element
        // Combine the hashes using XOR and bit-shifting
        return h1 ^ (h2 << 1) ^ (h3 << 2);
      }
    @endcode
*/
struct TupleHash {
    std::size_t operator()( const std::tuple<uint32_t, uint32_t, uint32_t>& t ) const {
        std::size_t h1 = std::hash<uint32_t>{}(std::get<0>(t));
        std::size_t h2 = std::hash<uint32_t>{}(std::get<1>(t));
        std::size_t h3 = std::hash<uint32_t>{}(std::get<2>(t));
        // combine the hashes using XOR and bit-shifting
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};


/**
    Three-dimensional jagged array (where the number of columns can vary),  implemented via an unordered_map plus hashing function 'TupleHash'
    Tested and compared in 'JaggedArray3D_Comparison_Test'.
    For dense regularly shaped matrix, JaggedArray3D is slower than DynamicArray3D and a vector of vecto of vectors.
    Good uses are for sparse 3D matrices the entries of which are in a flux (changing rapidly during a computation).
*/
template<class T>
class JaggedArray3D {
  public:
    using key = std::tuple<uint32_t, uint32_t, uint32_t>;
    using umap = std::unordered_map<key,T,TupleHash>;
    // ref
    using reference       = typename umap::reference;
    using const_reference = typename umap::const_reference;
    // iter
    using iterator       = typename umap::iterator;
    using const_iterator = typename umap::const_iterator;

    /// default constructor constructs an empty array
    JaggedArray3D() = default;
    
    /// construction using initialiser list (which must be known at compile time)
    explicit JaggedArray3D( std::initializer_list<std::initializer_list<std::initializer_list<T>>> );
    
    /// construct from a monster vector of vector of vectors
    explicit JaggedArray3D( const std::vector<std::vector<std::vector<T>>>& );
    
    // ACCESSORS
    
    /// iterators
    iterator begin() { return data_.begin(); }
    iterator end() { return data_.end(); }
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
  
    /// Insert or modify elements
    void set( uint32_t k, uint32_t i, uint32_t j, T value ) {
        data_[{k, i, j}] = value;
    }

    /// read elements using operator() return NaN if they don't exist
    T operator()( uint32_t k, uint32_t i, uint32_t j ) const {
        key index = {k, i, j};
        auto it = data_.find(index);
        if (it != data_.end()) {
            return it->second;  // Element exists, return its value
        } else {
            return std::numeric_limits<T>::quiet_NaN();  // Element does not exist
        }
    }
    
    /// read/modify elements (slower)
    T& operator()( uint32_t k, uint32_t i, uint32_t j ) {
        const key index = {k, i, j};
        return data_[ index ];
    }
    
    /// quickest way to test for element existance with a predefined tuple (use std::make_tuple)
    bool contains( const std::tuple<uint32_t, uint32_t, uint32_t>& index ) {
       if ( data_.find(index) != data_.end() ) return true;
       return false;
    }

   /// erase an element at a specific location
   void erase( uint32_t k, uint32_t i, uint32_t j ) {
        key index = {k, i, j};
        // key index = std::make_tuple(k, i, j ); // slower
        // check if the key exists in the map
        auto it = data_.find(index);
        if (it != data_.end()) {
            // key found, erase the element
            data_.erase(it);
        } else {
            // Key not found, throw an exception or silently ignore
            throw std::out_of_range("JaggedArray3D::erase:: array element does not exist at the specified indices.");
        }
    }

   // DIAGNOSTICS
   uint32_t slices() const;
   
   uint32_t rows() const;
   
   /// max because the number of columns can vary in a jagged array
   uint32_t max_cols() const;
   
   /// returns the  numbers of slices, rows, and maximum number of columns
   std::tuple<uint32_t, uint32_t, uint32_t> max_kij() const;

    /// print all elements to std out
    void Out() const {
        for (const auto& [index, value] : data_) {
            auto [k, i, j] = index;
            std::cout << "Value at (" << k << ", " << i << ", " << j << ") = " << value << std::endl;
        }
    }

  private:
    umap data_;
};

} // end csmp

#endif /* CSMP_JAGGED_ARRAY_3D_H */
