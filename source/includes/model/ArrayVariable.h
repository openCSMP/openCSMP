// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ARRAY_VARIABLE_H
#define ARRAY_VARIABLE_H

#include "CSMP_definitions.h"
#include "ScalarVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;

/**

  @brief Storage for an arbitrary number of variables with a single flag; array size
  is not predefined or dependent on model dimensions, but can be determined at runtime.

  @author P. Lang
  @date 2012

  In order to support arbitrary sized variables to be discretized in CSMP.

  @section usage Usage
  If used to read/store arrays from/on entities featuring a VariableStorage, the size of
  the ArrayVariable has to equal that specified in the variables file, since in the current
  version (9/4/2012) runtime changes regarding the PropertyDatabase and ArrayVariables are 
  not supported.
  
  To instantiate an ArrayVariable of arbitrary size:

  @code
  ArrayVariable av;                  // 0 entries, flag set to Plain
  ArrayVariable av(4, 10., DIRICH ); // 4 entries with value 10., flag set to Dirichlet
  @endcode

  When using ArrayVariable in Store/Read operations, we need to first create an instance of
  appropriate size, i.e. the same size as provided in the PropertyDatabase. A ctor which
  does exactly this is:

  @code
  ArrayVariable av( "element array", model.Database() ); // size as specified in PropertyDatabase, values = 0., flag = plain
  @endcode

  Again, the alternative would be to manually initiate the correct size.
  With an ArrayVariable of appropriate size, you may use functions like

  @code
  model.InputPropertyValue( "element array", av );
  @endcode

  Or something low level in the likes of

  @code
  Index eaKey( model.Database().StorageKey("element array") );
  for( vector<Element<3>*>::const_iterator it( model.Region("Model").CellsBegin() ); it != model.Region("Model").CellsEnd(); ++it )
    (*it)->Read( eaKey, av );
  @endcode

  ArrayVariable operations are possible as defined through the interface

  */
  class ArrayVariable
    {
    public:
      static constexpr VARIABLE_TYPE VariableType = ARRAY;

      typedef std::vector<double> ArrayContainer;

      ArrayVariable() noexcept;

      template<uint32_t dim>
      ArrayVariable( const char* arrayPropertyName,
                     const PropertyDatabase<dim>&,
                     double defaultValue=std::numeric_limits<double>::quiet_NaN(),
                     VARIABLE_FLAG flag = ANY );
      
      /// // use initializer list to initialize vector and flag: ArrayVariable ary( {1, 2, 3, 4}, ANY );
      ArrayVariable( std::initializer_list<double> initList, VARIABLE_FLAG flag )
        : data_(initList), flag_(flag) {}
        
      ArrayVariable( ArrayVariable&& other ) noexcept;

      ArrayVariable( uint32_t arraySize,
                     double defaultValue = 0.,
                     VARIABLE_FLAG flag = ANY ) noexcept;
      
      ArrayVariable( const Index& arrayKey,
                     double defaultValue = 0.,
                     VARIABLE_FLAG flag = ANY ) noexcept;
      
      ArrayVariable( const ArrayVariable& ) noexcept;

      /// assigments
      ArrayVariable& operator=( const ArrayVariable& ) noexcept;
      ArrayVariable& operator=( ArrayVariable&& other ) noexcept;
      void CopyValuesOnly( FlaggedArrayVariable& ) noexcept;
      ArrayVariable& operator=( double ) noexcept;
      ArrayVariable& operator=( const ScalarVariable& ) noexcept;

      /// standard operations
      ArrayVariable operator+( double ) const noexcept;
      ArrayVariable operator-( double ) const noexcept;
      ArrayVariable operator*( double ) const noexcept;
      ArrayVariable operator/( double ) const noexcept;
 
      ArrayVariable& operator+=( double ) noexcept;
      ArrayVariable& operator-=( double ) noexcept;
      ArrayVariable& operator*=( double ) noexcept;
      ArrayVariable& operator/=( double ) noexcept;

      ArrayVariable& operator+=( const ScalarVariable& ) noexcept;
      ArrayVariable& operator-=( const ScalarVariable& ) noexcept;
      ArrayVariable& operator*=( const ScalarVariable& ) noexcept;
      ArrayVariable& operator/=( const ScalarVariable& ) noexcept;

      /// value by value operations
      ArrayVariable& operator+=( const ArrayVariable& ) noexcept;
      ArrayVariable& operator-=( const ArrayVariable& ) noexcept;
      ArrayVariable& operator*=( const ArrayVariable& ) noexcept;
      ArrayVariable& operator/=( const ArrayVariable& ) noexcept;

      /// element-by-element operations
      ArrayVariable  operator-( const ArrayVariable& ) const noexcept;
      ArrayVariable  operator*( const ArrayVariable& ) const noexcept;
      ArrayVariable  operator+( const ArrayVariable& ) const noexcept;
      ArrayVariable  operator/( const ArrayVariable& ) const noexcept;
      ArrayVariable  Pow( double exponent ) const noexcept;
      ArrayVariable& PowInPlace( double exponent ) noexcept;

      /// comparison of flags and values
      bool           operator==( const ArrayVariable& ) const noexcept;
      bool           operator!=( const ArrayVariable& ) const noexcept;

      /// comparison of length (to allow ordering in containers)
      bool           operator<( const ArrayVariable& ) const noexcept;

      /// accessors to the data
      double&        operator()( uint32_t ) noexcept;
      double         operator[]( uint32_t ) const noexcept;
      
      /// functions needed to make CSMP variables interoperable; do not delete
      void           Component( uint32_t, double ) noexcept;
      double         Component( uint32_t ) const noexcept;

      /// returns number of array elements
      uint32_t       Size() const noexcept;
      void           Resize( uint32_t newSize, double newValue = std::numeric_limits<double>::quiet_NaN() ) noexcept;
      VARIABLE_FLAG  Flag(  ) const noexcept;
      VARIABLE_FLAG& Flag(  ) noexcept;
      void           Flag( VARIABLE_FLAG flag ) noexcept;
      bool           IsWithinRange( double min, double max ) const noexcept;
      bool           Has_NaN_Values() const noexcept;
      void           MinMax( double& min, double& max ) const noexcept;
      void           Sort() noexcept;

      /// printing array values to screen
      void           Out( long digits=3 ) const;
      bool           Out( const char* filename, long precision = 9 ) const;
      
      /// filestream I/O used for domain variables (not recommended because padding creates large storage overhead; use PropertyData instead)
      bool           Out( std::fstream& ) const;
      bool           In( std::fstream& );

      ArrayContainer::const_iterator Begin() const noexcept;
      ArrayContainer::const_iterator End()   const noexcept;

      // Standard range interface — enables range-based for and std:: algorithms.
      ArrayContainer::iterator       begin()        noexcept { return data_.begin(); }
      ArrayContainer::iterator       end()          noexcept { return data_.end();   }
      ArrayContainer::const_iterator begin()  const noexcept { return data_.begin(); }
      ArrayContainer::const_iterator end()    const noexcept { return data_.end();   }
      ArrayContainer::const_iterator cbegin() const noexcept { return data_.cbegin(); }
      ArrayContainer::const_iterator cend()   const noexcept { return data_.cend();   }

      /// To allow passing of array variables to GEMS.
      double& Front() noexcept { return data_.front(); }
      const double& Front() const noexcept { return data_.front(); }

      double NextLargestEntry(double) const noexcept;
      bool   HasLargerEntry(double) const noexcept;

    private:
      VARIABLE_FLAG  flag_;
      ArrayContainer data_;
    };

  std::ostream&  operator<<( std::ostream& stream, const ArrayVariable& o );


  // INLINE METHODS

inline ArrayVariable::ArrayVariable( ArrayVariable&& other ) noexcept
    : flag_( other.flag_ ),
      data_( std::move( other.data_ ) )
{}

inline ArrayVariable& ArrayVariable::operator=( ArrayVariable&& other ) noexcept
{
    if ( this != &other ) {
        flag_ = other.flag_;
        data_ = std::move( other.data_ );
    }
    return *this;
}

inline ArrayVariable::ArrayVariable() noexcept
    : flag_(ANY), data_()
    {
    }

  /// As PropertyDatabase ctor, but using index right away
inline ArrayVariable::ArrayVariable( const Index& arrayKey, double defaultValue, VARIABLE_FLAG flag ) noexcept
    : data_( arrayKey.dataDepth, defaultValue ), flag_(flag)
    {
    }  

inline ArrayVariable::ArrayVariable( uint32_t arraySize, double defaultValue, VARIABLE_FLAG flag ) noexcept
    : flag_(flag), data_( arraySize, defaultValue )
    {
    }

inline ArrayVariable::ArrayVariable( const ArrayVariable& av ) noexcept
    : flag_( av.flag_ ), data_( av.data_ )
    {
    }

inline double& ArrayVariable::operator()( uint32_t i ) noexcept
    {
      assert( i < Size() );
      return data_[i];
    }

inline double ArrayVariable::operator[]( uint32_t i ) const noexcept
    {
      assert( i < Size() );
      return data_[i];
    }

inline void  ArrayVariable::Component( uint32_t i, double val ) noexcept
   {
      assert( i < Size() );
      data_[i] = val;
   }

inline double  ArrayVariable::Component( uint32_t i ) const noexcept
   {
      assert( i < Size() );
      return data_[i];
   }

inline uint32_t ArrayVariable::Size() const noexcept
    {
      return static_cast<uint32_t>(data_.size());
    }

inline void ArrayVariable::Resize( uint32_t newSize, double newValue ) noexcept
    {
      data_.resize( newSize, newValue );
    }


inline VARIABLE_FLAG ArrayVariable::Flag() const noexcept
    {
      return flag_;
    }

inline VARIABLE_FLAG& ArrayVariable::Flag() noexcept
    {
      return flag_;
    }

inline void ArrayVariable::Flag( VARIABLE_FLAG flag ) noexcept
    {
      flag_ = flag;
    }


inline bool ArrayVariable::operator<( const ArrayVariable& other ) const noexcept
{
    return data_ < other.data_;
}


inline bool ArrayVariable::operator==( const ArrayVariable& other ) const noexcept
{
    return flag_ == other.flag_ && data_ == other.data_;
}


inline bool ArrayVariable::operator!=( const ArrayVariable& av ) const noexcept
  {
    return !(*this == av);
  }

inline ArrayVariable& ArrayVariable::operator=( const ArrayVariable& av ) noexcept
  {
    if( this != &av ) {
        flag_ = av.flag_;
        data_ = av.data_;
      }
    return *this;
  }


inline ArrayVariable& ArrayVariable::operator=( double val ) noexcept
  {
    for( auto& i : data_ ) i = val;
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator=( const ScalarVariable& val ) noexcept
  {
    for( auto& i : data_ ) i = val();
    return *this;
  }

/// Standart Operations with temporary object

inline ArrayVariable ArrayVariable::operator+( double val ) const noexcept
  {
    ArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) += val;
    return temp_arr;
}

inline ArrayVariable ArrayVariable::operator-( double val ) const noexcept
  {
    ArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) -= val;
    return temp_arr;
  }

inline ArrayVariable ArrayVariable::operator*( double val ) const noexcept
  {
    ArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) *= val;
    return temp_arr;
  }

inline ArrayVariable ArrayVariable::operator/( double val ) const noexcept
  {
    ArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) /= val;
    return temp_arr;
  }

inline ArrayVariable ArrayVariable::Pow( double exponent ) const noexcept
  {
    ArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) = std::pow(temp_arr(i), exponent );
    return temp_arr;
  }


inline ArrayVariable& ArrayVariable::PowInPlace( double exponent ) noexcept
{
    for ( auto& v : data_ )
        v = std::pow( v, exponent );
    return *this;
}

/// Standart Operations with current object

inline ArrayVariable& ArrayVariable::operator+=( double val ) noexcept
  {
    for( auto& i : data_ ) i += val;
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator-=( double val ) noexcept
  {
    for( auto& i : data_ ) i -= val;
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator*=( double val ) noexcept
  {
    for( auto& i : data_ ) i *= val;
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator/=( double val ) noexcept
  {
    for( auto& i : data_ ) i /= val;
    return *this;
  }

/// Standart Operations with current object

/// Operations with ScalarVariables
inline ArrayVariable& ArrayVariable::operator+=( const ScalarVariable& sc ) noexcept
  {
    for( auto& i : data_ ) i += sc();
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator-=( const ScalarVariable& sc ) noexcept
  {
    for( auto& i : data_ ) i -= sc();
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator*=( const ScalarVariable& sc ) noexcept
  {
    for( auto& i : data_ ) i *= sc();
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator/=( const ScalarVariable& sc ) noexcept
  {
    for( auto& i : data_ ) i /= sc();
    return *this;
  }

/// Operations with ArrayVariables
inline ArrayVariable& ArrayVariable::operator+=( const ArrayVariable& arr ) noexcept
  {
    for ( uint32_t i{0U}; i < Size(); ++i ) data_[i] += arr[i];
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator-=( const ArrayVariable& arr ) noexcept
  {
    for ( uint32_t i{0U}; i < Size(); ++i ) data_[i] -= arr[i];
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator*=( const ArrayVariable& arr ) noexcept
  {
    for ( uint32_t i{0U}; i < Size(); ++i ) data_[i] *= arr[i];
    return *this;
  }

inline ArrayVariable& ArrayVariable::operator/=( const ArrayVariable& arr ) noexcept
  {
    for ( uint32_t i{0U}; i < Size(); ++i ) data_[i] /= arr[i];
    return *this;
  }
  
  
inline ArrayVariable ArrayVariable::operator-( const ArrayVariable& av ) const noexcept
  {
    ArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
      returnArray(i) -=  av[i];
    return returnArray;
  }

inline ArrayVariable ArrayVariable::operator+( const ArrayVariable& av ) const noexcept
  {
    ArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
      returnArray(i) +=  av[i];
    return returnArray;
  }

inline ArrayVariable ArrayVariable::operator*( const ArrayVariable& av ) const noexcept
  {
    ArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
      returnArray(i) *=  av[i];
    return returnArray;
  }


inline ArrayVariable ArrayVariable::operator/( const ArrayVariable& av ) const noexcept
  {
    ArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
      returnArray(i) /=  av[i];
    return returnArray;
  }


inline bool ArrayVariable::IsWithinRange( double min, double max ) const noexcept
{
    if ( min >= max ) return false;  // degenerate range — nothing can be in it
    for ( const double v : data_ )
        if ( v < min || v > max ) return false;
    return true;
}


inline bool ArrayVariable::Has_NaN_Values() const noexcept
{
    return std::any_of( data_.cbegin(), data_.cend(),
                        []( double v ) { return std::isnan(v); } );
}


inline void ArrayVariable::MinMax( double& min, double& max ) const noexcept
{
    if ( data_.empty() )
    {
        min = std::numeric_limits<double>::quiet_NaN();
        max = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const auto [minIt, maxIt] =
        std::minmax_element( data_.cbegin(), data_.cend() );
    min = *minIt;
    max = *maxIt;
}


inline double ArrayVariable::NextLargestEntry( double fromValue ) const noexcept
{
    // Linear scan for the smallest value strictly greater than fromValue.
    // O(n), no allocation, no sort.
    // Returns numeric_limits<double>::max() if no larger entry exists.
    double result = std::numeric_limits<double>::max();
    for ( const double v : data_ )
        if ( v > fromValue && v < result )
            result = v;
    return result;
}

inline bool ArrayVariable::HasLargerEntry( double fromValue ) const noexcept
{
    // Linear scan — stops at first element greater than fromValue.
    // O(n) worst case, O(1) best case, no allocation, no sort.
    return std::any_of( data_.cbegin(), data_.cend(),
                        [fromValue]( double v ) { return v > fromValue; } );
}


inline void ArrayVariable::Sort() noexcept
{
  sort( data_.begin(), data_.end() );
}


inline ArrayVariable::ArrayContainer::const_iterator ArrayVariable::Begin() const noexcept
{
  return data_.begin();
}


inline ArrayVariable::ArrayContainer::const_iterator ArrayVariable::End() const noexcept
{
  return data_.end();
}


  } // csmp

#endif
