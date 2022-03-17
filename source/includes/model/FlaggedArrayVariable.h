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
      static constexpr VARIABLE_TYPE VariableType = FLAGGEDARRAY;

      typedef std::vector<double> FlaggedArrayContainer;

      FlaggedArrayVariable();

      FlaggedArrayVariable( const FlaggedArrayVariable& );

      template<uint32_t dim>
      FlaggedArrayVariable( const char* arrayPropertyName,
                            const PropertyDatabase<dim>& pd,
                            double defaultValue=std::numeric_limits<double>::quiet_NaN(),
                            VARIABLE_FLAG flag = ANY );
      
      explicit FlaggedArrayVariable( size_t arraySize,
                                     double defaultValue = 0.,
                                     VARIABLE_FLAG flag = ANY );
      
      explicit FlaggedArrayVariable( const Index& arrayKey,
                                     double defaultValue = 0.,
                                     VARIABLE_FLAG flag = ANY );
                                     
      /// assignments
      FlaggedArrayVariable& operator=( const FlaggedArrayVariable& );
      void CopyValuesOnly( ArrayVariable& av);
      FlaggedArrayVariable& operator=( double );
      FlaggedArrayVariable& operator=( const ScalarVariable& );

      /// standard operations
      FlaggedArrayVariable operator+( double ) const;
      FlaggedArrayVariable operator-( double ) const;
      FlaggedArrayVariable operator*( double ) const;
      FlaggedArrayVariable operator/( double ) const;
      FlaggedArrayVariable operator^( double ) const;

      FlaggedArrayVariable& operator+=( double );
      FlaggedArrayVariable& operator-=( double );
      FlaggedArrayVariable& operator*=( double );
      FlaggedArrayVariable& operator/=( double );

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
      double         operator[]( size_t ) const;
      double&        operator()( size_t );
      void           Component( size_t, double );
      double         Component( size_t i ) const;

      size_t         Size() const;
      void           Resize( size_t newSize, double newValue = std::numeric_limits<double>::quiet_NaN() );
      VARIABLE_FLAG  Flag( const size_t&) const;
      VARIABLE_FLAG& Flag( const size_t&);
      void           Flag( const size_t&, VARIABLE_FLAG);
      void           Fabs();
      void           Ln();
      void           Log10();
      void           Sqrt();
      bool           IsWithinRange( double min, double max ) const;
      void           MinMax( double& min, double& max ) const;
      void           Sort();

      /// printing array values to screen
      void           Out( long digits=3 ) const;
      bool           Out( const char* filename, size_t precision = 9 ) const;
      
      /// SKM: DO NOT USE C-STYLE I/O - for internal csmp binary IO
      bool           Out( std::fstream& fp ) const;
      bool           In( std::fstream& fp );

      FlaggedArrayContainer::const_iterator Begin() const;
      FlaggedArrayContainer::const_iterator End()   const;

      double NextLargestEntry(double) const;
      bool     HasLargerEntry(double) const;

    private:
      std::vector<VARIABLE_FLAG>  flags_;
      FlaggedArrayContainer       data_;
    };

  std::ostream&  operator<<( std::ostream& stream, const FlaggedArrayVariable& o );


  /// Automatically sets size to corresponding index, flag to ANY
  template<uint32_t dim>
  FlaggedArrayVariable::FlaggedArrayVariable( const char* arrayPropertyName,
                                              const PropertyDatabase<dim>& pd,
                                              double defaultValue,
                                              VARIABLE_FLAG flag )
    : data_       ( pd.StorageKey(arrayPropertyName).dataDepth, defaultValue ),
      flags_      ( pd.StorageKey(arrayPropertyName).dataDepth, flag )
    {
    }


  } // csmp

#endif
