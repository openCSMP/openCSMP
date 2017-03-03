#ifndef FLAGGED_ARRAY_VARIABLE_H
#define FLAGGED_ARRAY_VARIABLE_H

#include "ArrayVariable.h"
#include "PropertyDatabase.h"


namespace csmp {

class ScalarVariable;
class ArrayVariable;

/**

 @brief Storage for an arbitrary number of variables and flags with a size 
  that is not predefined but can be determined at runtime.

  @author Roman Manasipov, Philipp Lang
  @date 2013

  In order to support arbitrary sized variables to be discretized in CSMP.
  Based on implementation of ArrayVariable by Philipp Lang but has additional
  container of flags for each data component

  @section usage Usage
  
  To instantiate an FlaggedArrayVariable of arbitrary size:

  @code
  FlaggedArrayVariable fav;                  // 0 entries, flag set to Any
  FlaggedArrayVariable fav(4, 10., DIRICH ); // 4 entries with value 10., flag set to Dirichlet
  @endcode

  When using FlaggedArrayVariable in Store/Read operations, we need to first create an instance of
  appropriate size, i.e. the same size as provided in the PropertyDatabase. A constructor which
  does exactly this:

  @code
  FlaggedArrayVariable fav( "element array", model.Database() ); // size as specified in PropertyDatabase, values = 0., flag = any
  @endcode

  Again, the alternative would be to manually initiate the correct size.
  With an FlaggedArrayVariable of appropriate size, you may use functions like

  @code
  model.InputPropertyValue( "element flagged array", fav );
  @endcode

  Or something low level in the likes of

  @code
  Index efaKey( model.Database().StorageKey("element flagged array") );
  for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
    (*it)->Read( efaKey, fav );
  @endcode

  FlaggedArrayVariable operations are possible as defined through the interface

  */
class FlaggedArrayVariable
    {
    public:
      typedef std::vector<double64> FlaggedArrayContainer;

      FlaggedArrayVariable();

      FlaggedArrayVariable( const FlaggedArrayVariable& );

      template<size_t dim>
      FlaggedArrayVariable( const char* arrayPropertyName,
                            const PropertyDatabase<dim>& pd,
                            double64 defaultValue=std::numeric_limits<double64>::quiet_NaN(),
                            VARIABLE_FLAG flag = ANY );
      
      explicit FlaggedArrayVariable( size_t arraySize,
                                     double64 defaultValue = 0.,
                                     VARIABLE_FLAG flag = ANY );
      
      explicit FlaggedArrayVariable( const Index& arrayKey,
                                     double64 defaultValue = 0.,
                                     VARIABLE_FLAG flag = ANY );
                                     
      /// assignments
      FlaggedArrayVariable& operator=( const FlaggedArrayVariable& );
      void CopyValuesOnly( ArrayVariable& av);
      FlaggedArrayVariable& operator=( double64 );
      FlaggedArrayVariable& operator=( const ScalarVariable& );

      /// standard operations
      FlaggedArrayVariable operator+( double64 ) const;
      FlaggedArrayVariable operator-( double64 ) const;
      FlaggedArrayVariable operator*( double64 ) const;
      FlaggedArrayVariable operator/( double64 ) const;
      FlaggedArrayVariable operator^( double64 ) const;

      FlaggedArrayVariable& operator+=( double64 );
      FlaggedArrayVariable& operator-=( double64 );
      FlaggedArrayVariable& operator*=( double64 );
      FlaggedArrayVariable& operator/=( double64 );

      FlaggedArrayVariable& operator+=( const ScalarVariable& );
      FlaggedArrayVariable& operator-=( const ScalarVariable& );
      FlaggedArrayVariable& operator*=( const ScalarVariable& );
      FlaggedArrayVariable& operator/=( const ScalarVariable& );

      /// value by value operations
      FlaggedArrayVariable& operator+=( const FlaggedArrayVariable& );
      FlaggedArrayVariable& operator-=( const FlaggedArrayVariable& );
      FlaggedArrayVariable& operator*=( const FlaggedArrayVariable& );
      FlaggedArrayVariable& operator/=( const FlaggedArrayVariable& );

      /// element-by-element operations
      FlaggedArrayVariable  operator-( const FlaggedArrayVariable& ) const;
      FlaggedArrayVariable  operator*( const FlaggedArrayVariable& ) const;
      FlaggedArrayVariable  operator+( const FlaggedArrayVariable& ) const;
      FlaggedArrayVariable  operator/( const FlaggedArrayVariable& ) const;

      /// arithmetical comparison of array values
      bool           operator==( const FlaggedArrayVariable& ) const;
      bool           operator!=( const FlaggedArrayVariable& ) const;
      bool           operator<( const FlaggedArrayVariable& ) const;
      bool           operator<=( const FlaggedArrayVariable& ) const;
      bool           operator>( const FlaggedArrayVariable& ) const;
      bool           operator>=( const FlaggedArrayVariable& ) const;

      /// accessors to the data
      double64       operator[]( size_t ) const;
      double64&      operator()( size_t );
      void           Component( size_t, double64 );
      double64       Component( size_t i ) const;

      size_t         Size() const;
      void           Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
      VARIABLE_FLAG  Flag( const size_t&) const;
      VARIABLE_FLAG& Flag( const size_t&);
      void           Flag( const size_t&, VARIABLE_FLAG);
      void           Fabs();
      void           Ln();
      void           Log10();
      void           Sqrt();
      bool           IsWithinRange( double64 min, double64 max ) const;
      void           MinMax( double64& min, double64& max ) const;
      void           Sort();

      /// printing array values to screen
      void           Out( long digits=3 ) const;
      bool           Out( const char* filename, size_t precision = 9 ) const;
      
      /// SKM: DO NOT USE C-STYLE I/O - for internal csmp binary IO
      bool           Out( std::FILE* fp ) const;
      bool           In( std::FILE* fp );

      FlaggedArrayContainer::const_iterator Begin() const;
      FlaggedArrayContainer::const_iterator End()   const;

      double64 NextLargestEntry(double64) const;
      bool     HasLargerEntry(double64) const;

    private:
      std::vector<VARIABLE_FLAG>  flags_;
      FlaggedArrayContainer       data_;
    };

  std::ostream&  operator<<( std::ostream& stream, const FlaggedArrayVariable& o );


  /// Automatically sets size to corresponding index, flag to ANY
  template<size_t dim>
  FlaggedArrayVariable::FlaggedArrayVariable( const char* arrayPropertyName,
                                              const PropertyDatabase<dim>& pd,
                                              double64 defaultValue,
                                              VARIABLE_FLAG flag )
    : data_       ( pd.StorageKey(arrayPropertyName).dataDepth, defaultValue ),
      flags_      ( pd.StorageKey(arrayPropertyName).dataDepth, flag )
    {
    }


  } // csmp

#endif
