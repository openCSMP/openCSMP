// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DynamicArray3D.h - based on STL vector
//  Open CSMP++
//
//  Extended from PixelChemist (2013)'s DynamicArray2D that was shared via C++ StackExchange June 23,2013
//  extension by SKM
//

#ifndef DYNAMIC_ARRAY_3D_H
#define DYNAMIC_ARRAY_3D_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <iostream>
#include <cassert>

/**

@brief DynamicArray3D - three dimensional array based on vector

dataType arrayName[depth][row size][column size] = {
    {
        {depth0row0col0, depth0row0col1, ...},
        {depth0row1col0, depth0row1col1, ...},
    }
 
 usage with operator(depth,row,column)
      @code
        for (int i = 0; i < depth; i++) {
        printf("Depth %d:\n", i + 1); 
        for (int j = 0; j < rows; j++) { 
            for (int k = 0; k < cols; k++) { 
                printf("%d ", arr[i][j][k]); 
            } 
            printf("\n"); 
        } 
        printf("\n");
        }
   @endcode
       
   directly initializing 3D array at the time of declaration
   @code
        int arr1[2][3][4] = { { { 1, 2, 3, 4 },
                                { 5, 6, 7, 8 },
                                { 9, 10, 11, 12 } },
                              { { 13, 14, 15, 16 },
                                { 17, 18, 19, 20 },
                                { 21, 22, 23, 24 } } };
   @endcode
   
   Printing arr1
                  Depth 1:
                  1 2 3 4
                  5 6 7 8
                  9 10 11 12

                  Depth 2:
                  13 14 15 16
                  17 18 19 20
                  21 22 23 24
*/
template<class T>
class DynamicArray3D
  {
  public:
    // misc types
    using data_type  = std::vector<T>;
    using value_type = typename std::vector<T>::value_type;
    using size_type  = typename std::vector<T>::size_type;
    // ref
    using reference       = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;
    // iter
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    // reverse iter
    using reverse_iterator       = typename std::vector<T>::reverse_iterator;
    using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

    // empty construction
    DynamicArray3D() = default;

    // default-insert depth * row * column values
    DynamicArray3D( size_type depth, size_type rows, size_type cols )
      : m_depth(depth), m_rows(rows), m_cols(cols), m_data(depth*rows*cols)
    { /* no values to take no speed hit from initialisation */ }

    // initialized with matrix values
    DynamicArray3D( size_type depth, size_type rows, size_type cols, const_reference val )
      : m_depth(depth), m_rows(rows), m_cols(cols), m_data(depth*rows*cols, val)
    {}

    /**
       matrix constructor taking initialiser list
       supporting initialisation by:
              DynamicArray3D<double> (2 x 3 x 2)
      @code
            array[i0] => [j0] => [k0, k1]
                         [j1] => [k0, k1]
                         [j2] => [k0, k1]
            array[1]  => [0] => [0, 1]
                         [1] => [0, 1]
                         [2] => [0, 1]
            
           // (2 x 3 x 4)
           DynamicArray3D arr1 = { { { 1, 2, 3, 4 },
                                        { 5, 6, 7, 8 },
                                        { 9, 10, 11, 12 } },
                                      { { 13, 14, 15, 16 },
                                        { 17, 18, 19, 20 },
                                        { 21, 22, 23, 24 } } };
       @endcode
    */
    DynamicArray3D( std::initializer_list<std::initializer_list<data_type>> vals )
      : m_depth(vals.size()), // outer list
        m_rows((vals.begin())->size()), // inner list
        m_cols((*std::next(vals.begin(),1)).size()), // within inner list
        m_data(m_rows*m_cols*m_depth,0)
    {
       // assigning the values
       size_type dij = 0;
       for ( const auto& dep : vals )
         for ( const auto& row : dep )
           for ( const auto& col : row )
             // initialising to Nan since this method does not work yet
             m_data[dij++] = col; // std::numeric_limits<double>::quiet_NaN();
    }


    // 1d-iterators
    iterator begin() { return m_data.begin(); }
    iterator end() { return m_data.end(); }
    const_iterator begin() const { return m_data.begin(); }
    const_iterator end() const { return m_data.end(); }
    const_iterator cbegin() const { return m_data.cbegin(); }
    const_iterator cend() const { return m_data.cend(); }
    reverse_iterator rbegin() { return m_data.rbegin(); }
    reverse_iterator rend() { return m_data.rend(); }
    const_reverse_iterator rbegin() const { return m_data.rbegin(); }
    const_reverse_iterator rend() const { return m_data.rend(); }
    const_reverse_iterator crbegin() const { return m_data.crbegin(); }
    const_reverse_iterator crend() const { return m_data.crend(); }

    // element access (row major indexation)
    reference operator() ( size_type const depth,
                           size_type const row,
                           size_type const column )
      {
        // Multi-indexed array with dimensions (N,M,L) has one-dimensional (flat) access via
        //    array[i][j][k] = array(i,j,k) = array[M*L*i + L*j + k]
        return m_data[ m_rows * m_cols * depth + m_cols * row + column ];
      }
    // const element access (row major indexation)
    const_reference operator() ( size_type const depth,
                                 size_type const row,
                                 size_type const column ) const {
        return m_data[ m_rows * m_cols * depth + m_cols * row + column ];
      }
    reference at( size_type const depth, size_type const row, size_type const column ) {
        return m_data[ m_rows * m_cols * depth + m_cols * row + column ];
      }
    const_reference at( size_type const depth, size_type const row, size_type const column ) const {
        return m_data[ m_rows * m_cols * depth + m_cols * row + column ];
      }

    // resizing while retaining the existing elemments
    void resize( size_type new_depth, size_type new_rows, size_type new_cols ) {
        // new matrix new_rows times new_cols
        DynamicArray3D tmp( new_depth, new_rows, new_cols );
        // select smaller row and col size
        auto md = std::min(m_depth, new_depth);
        auto mr = std::min(m_rows, new_rows);
        auto mc = std::min(m_cols, new_cols);
        for ( size_type d(0U); d < md; ++d )
          for ( size_type i(0U); i < mr; ++i ) {
              // iterators to beginning of rows
              auto row     = begin() + i * m_cols + d * m_rows * m_cols;
              auto tmp_row = tmp.begin() + i * new_cols + d * new_rows * new_cols;
              // move mc elements from array to tmp array
              std::move( row, row + mc, tmp_row ); // move row by row
            }
        // move assignment to this (original array)
        *this = std::move(tmp);
      }

    // size and capacity
    size_type size() const { return m_data.size(); }
    size_type capacity() const { return m_data.capacity(); }
    bool      empty() const { return m_data.empty(); }
    
    // - maximum entries that could be help in data vector theoretically
    size_type max_size() const { return m_data.max_size(); }
    
    // dimensionality
    size_type depth() const { return m_depth; }
    size_type rows() const { return m_rows; }
    size_type cols() const { return m_cols; }
    
    // data swapping
    void swap( DynamicArray3D& rhs ) {
      using std::swap;
      m_data.swap(rhs.m_data);
      swap(m_depth, rhs.m_depth);
      swap(m_rows, rhs.m_rows);
      swap(m_cols, rhs.m_cols);
    }
    
    // print array to screen
    void out() const;
  
  private: // DATA MEMBERS
    size_type m_depth{ 0u };
    size_type m_rows{ 0u };
    size_type m_cols{ 0u };
    data_type m_data{};
  };
  

template<class T>
inline void swap(DynamicArray3D<T>& lhs, DynamicArray3D<T>& rhs)
  {
    lhs.swap(rhs);
  }

template<class T>
inline bool operator==( const DynamicArray3D<T>& a, const DynamicArray3D<T>& b )
  {
    if (a.rows() != b.rows() || a.cols() != b.cols() || a.depth() != b.depth() )
    {
      return false;
    }
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }
  
template<class T>
inline bool operator!=( DynamicArray3D<T> const &a, DynamicArray3D<T> const &b )
  {
    return !(a == b);
  }

// printing array to screen
template<class T>
void DynamicArray3D<T>::out() const
 {
    std::cout <<"\n"<<"DynamicArray3D: ("<< depth() <<" x "<< rows() <<" x " << cols() <<"):\n";
    for ( size_type i{0}; i<depth(); ++i ) {
        std::cout <<"\t"<<"i=depth="<< i << std::endl;
        for ( size_type j{0}; j<rows(); ++j ) {
            std::cout <<"\t";
            for ( size_type k{0}; k<cols(); ++k )
              std::cout <<" "<< (*this)(i,j,k);
            std::cout << std::endl;
          }
      }
    std::cout << std::endl;
 }

#endif /* DYNAMIC_ARRAY_3D_H */
