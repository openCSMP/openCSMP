//
//  DynamicArray2D.h - based on STL vector
//  Open CSMP++
//
//  Created by PixelChemist (2013), shared via C++ StackExchange June 23,2013
//

#ifndef DYNAMIC_ARRAY_2D_H
#define DYNAMIC_ARRAY_2D_H

#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>

template<class T>
class DynamicArray2D
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
    DynamicArray2D() = default;

    // default-insert rows*cols values
    DynamicArray2D(size_type rows, size_type cols)
      : m_rows(rows), m_cols(cols), m_data(rows*cols)
    {}

    // copy initialized matrix rows*cols
    DynamicArray2D(size_type rows, size_type cols, const_reference val)
      : m_rows(rows), m_cols(cols), m_data(rows*cols, val)
    {}

    /**
       matrix constructor taking vector-of-vectors-style initialiser list
       supporting initialisation by:
              DynamicArray2D<double> mat4{ {1., 0.}, {0., 2.} };
    */
    DynamicArray2D( std::initializer_list<data_type> vals )
      : m_rows(vals.size()), m_cols((*std::next(vals.begin(),0)).size()),
        m_data(m_rows*m_cols,0.)
    {
       assert( m_rows*m_cols == vals.size() );
       // assigning the values
       size_type ij = 0;
       for ( const auto& row : vals )
         for ( const auto& col : row )
           m_data[ij++] = col;
    }


    /// construction of square array via vector-style initialiser list
    /*
    DynamicArray2D( std::initializer_list<T> vals )
      : m_rows(std::sqrt(vals.size())), m_cols(std::sqrt(vals.size()), m_data(vals)
    {
       assert( m_rows*m_cols() == vals.size() );
    }
    */
    
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
    reference operator() (size_type const row,
                          size_type const column)
    {
      return m_data[m_cols*row + column];
    }
    const_reference operator() (size_type const row,
                                size_type const column) const
    {
      return m_data[m_cols*row + column];
    }
    reference at( size_type const row, size_type const column )
    {
      return m_data.at(m_cols*row + column);
    }
    const_reference at( size_type const row, size_type const column ) const
    {
      return m_data.at(m_cols*row + column);
    }

    // resizing
    void resize(size_type new_rows, size_type new_cols)
    {
      // new matrix new_rows times new_cols
      DynamicArray2D tmp(new_rows, new_cols);
      // select smaller row and col size
      auto mc = std::min(m_cols, new_cols);
      auto mr = std::min(m_rows, new_rows);
      for (size_type i(0U); i < mr; ++i)
      {
        // iterators to begin of rows
        auto row = begin() + i*m_cols;
        auto tmp_row = tmp.begin() + i*new_cols;
        // move mc elements to tmp
        std::move(row, row + mc, tmp_row);
      }
      // move assignment to this
      *this = std::move(tmp);
    }

    // size and capacity
    size_type size() const { return m_data.size(); }
    size_type max_size() const { return m_data.max_size(); }
    bool empty() const { return m_data.empty(); }
    // dimensionality
    size_type rows() const { return m_rows; }
    size_type cols() const { return m_cols; }
    // data swapping
    void swap(DynamicArray2D& rhs)
    {
      using std::swap;
      m_data.swap(rhs.m_data);
      swap(m_rows, rhs.m_rows);
      swap(m_cols, rhs.m_cols);
    }
  private:
    // content
    size_type m_rows{ 0u };
    size_type m_cols{ 0u };
    data_type m_data{};
  };
  template<class T>
  void swap(DynamicArray2D<T>& lhs, DynamicArray2D<T>& rhs)
  {
    lhs.swap(rhs);
  }
  template<class T>
  bool operator== (DynamicArray2D<T> const &a, DynamicArray2D<T> const &b)
  {
    if (a.rows() != b.rows() || a.cols() != b.cols())
    {
      return false;
    }
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }
  template<class T>
  bool operator!= (DynamicArray2D<T> const &a, DynamicArray2D<T> const &b)
  {
    return !(a == b);
  }

#endif /* DYNAMIC_ARRAY_2D_H */
