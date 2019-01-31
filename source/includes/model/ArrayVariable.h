#ifndef ARRAY_VARIABLE_H
#define ARRAY_VARIABLE_H

#include "FlaggedArrayVariable.h"
#include "PropertyDatabase.h"

namespace csmp {

class ScalarVariable;

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
  for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
    (*it)->Read( eaKey, av );
  @endcode

  ArrayVariable operations are possible as defined through the interface

  
  // TODO: why can this not be of a size that is known at compile time
  // TODO: check whether <array> may be a good basis for refactoring
  // TODO: why - as this is essentially a vector - is it not inherited from a vector?

  */
  class FlaggedArrayVariable;
  
  class ArrayVariable
    {
    public:
      static constexpr VARIABLE_TYPE VariableType = ARRAY;

      typedef std::vector<double64> ArrayContainer;

      ArrayVariable();
      template<size_t dim>
      ArrayVariable( const char* arrayPropertyName,
                     const PropertyDatabase<dim>&,
                     double64 defaultValue=std::numeric_limits<double64>::quiet_NaN(),
                     VARIABLE_FLAG flag = ANY );
      
      explicit ArrayVariable( size_t arraySize,
                              double64 defaultValue = 0.,
                              VARIABLE_FLAG flag = ANY );
      
      explicit ArrayVariable( const Index& arrayKey,
                              double64 defaultValue = 0.,
                              VARIABLE_FLAG flag = ANY );
      
      ArrayVariable( const ArrayVariable& );

      /// assigments
      ArrayVariable& operator=( const ArrayVariable& );
      void CopyValuesOnly( FlaggedArrayVariable& );
      ArrayVariable& operator=( double64 );
      ArrayVariable& operator=( const ScalarVariable& );

      /// standard operations
      ArrayVariable operator+( double64 ) const;
      ArrayVariable operator-( double64 ) const;
      ArrayVariable operator*( double64 ) const;
      ArrayVariable operator/( double64 ) const;
      ArrayVariable operator^( double64 ) const;

      ArrayVariable& operator+=( double64 );
      ArrayVariable& operator-=( double64 );
      ArrayVariable& operator*=( double64 );
      ArrayVariable& operator/=( double64 );

      ArrayVariable& operator+=( const ScalarVariable& );
      ArrayVariable& operator-=( const ScalarVariable& );
      ArrayVariable& operator*=( const ScalarVariable& );
      ArrayVariable& operator/=( const ScalarVariable& );

      /// value by value operations
      ArrayVariable& operator+=( const ArrayVariable& );
      ArrayVariable& operator-=( const ArrayVariable& );
      ArrayVariable& operator*=( const ArrayVariable& );
      ArrayVariable& operator/=( const ArrayVariable& );

      /// element-by-element operations
      ArrayVariable  operator-( const ArrayVariable& ) const;
      ArrayVariable  operator*( const ArrayVariable& ) const;
      ArrayVariable  operator+( const ArrayVariable& ) const;
      ArrayVariable  operator/( const ArrayVariable& ) const;

      /// comparison of flags and values
      bool           operator==( const ArrayVariable& ) const;
      bool           operator!=( const ArrayVariable& ) const;

      /// comparison of length (to allow ordering in containers)
      bool           operator<( const ArrayVariable& ) const;

      /// accessors to the data
      double64&      operator()( size_t );
      double64       operator[]( size_t ) const;
      void           Component( size_t, double64 );
      double64       Component( size_t ) const;

      size_t         Size() const;
      void           Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
      VARIABLE_FLAG  Flag(  ) const;
      VARIABLE_FLAG& Flag(  );
      void           Flag( VARIABLE_FLAG flag );
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
      
      /// printing values to filestream
      bool           Out( std::fstream& ) const;
      bool           In( std::fstream& );

      ArrayContainer::const_iterator Begin() const;
      ArrayContainer::const_iterator End()   const;

      /// Julian 17-07-2014 - This is to allow passing of array variables to GEMS.
      double64& Front(){ return data_.front(); }

      double64 NextLargestEntry(double64) const;
      bool     HasLargerEntry(double64) const;

    private:
      VARIABLE_FLAG  flag_;
      ArrayContainer data_;
    };

  std::ostream&  operator<<( std::ostream& stream, const ArrayVariable& o );


  } // csmp

#endif
