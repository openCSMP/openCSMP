#ifndef ARRAY_VARIABLE_H
#define ARRAY_VARIABLE_H

#include "CSMP_definitions.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

class ScalarVariable;
class FlaggedArrayVariable;
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

  
  // TODO: why can this not be of a size that is known at compile time
  // TODO: check whether <array> may be a good basis for refactoring
  // TODO: why - as this is essentially a vector - is it not inherited from a vector?

  */
  class ArrayVariable
    {
    public:
      static constexpr VARIABLE_TYPE VariableType = ARRAY;

      typedef std::vector<double> ArrayContainer;

      ArrayVariable();
      template<uint32_t dim>
      ArrayVariable( const char* arrayPropertyName,
                     const PropertyDatabase<dim>&,
                     double defaultValue=std::numeric_limits<double>::quiet_NaN(),
                     VARIABLE_FLAG flag = ANY );
      
      ArrayVariable( uint32_t arraySize,
                     double defaultValue = 0.,
                     VARIABLE_FLAG flag = ANY );
      
      ArrayVariable( const Index& arrayKey,
                     double defaultValue = 0.,
                     VARIABLE_FLAG flag = ANY );
      
      ArrayVariable( const ArrayVariable& );

      /// assigments
      ArrayVariable& operator=( const ArrayVariable& );
      void CopyValuesOnly( FlaggedArrayVariable& );
      ArrayVariable& operator=( double );
      ArrayVariable& operator=( const ScalarVariable& );

      /// standard operations
      ArrayVariable operator+( double ) const;
      ArrayVariable operator-( double ) const;
      ArrayVariable operator*( double ) const;
      ArrayVariable operator/( double ) const;
      ArrayVariable operator^( double ) const;

      ArrayVariable& operator+=( double );
      ArrayVariable& operator-=( double );
      ArrayVariable& operator*=( double );
      ArrayVariable& operator/=( double );

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
      double&        operator()( uint32_t );
      double         operator[]( uint32_t ) const;
      
      /// functions needed to make CSMP variables interoperable; do not delete
      void           Component( uint32_t, double );
      double         Component( uint32_t ) const;

      /// returns number of array elements
      uint32_t       Size() const;
      void           Resize( uint32_t newSize, double newValue = std::numeric_limits<double>::quiet_NaN() );
      VARIABLE_FLAG  Flag(  ) const;
      VARIABLE_FLAG& Flag(  );
      void           Flag( VARIABLE_FLAG flag );
      bool           IsWithinRange( double min, double max ) const;
      bool           Has_NaN_Values() const;
      void           MinMax( double& min, double& max ) const;
      void           Sort();

      /// printing array values to screen
      void           Out( long digits=3 ) const;
      bool           Out( const char* filename, size_t precision = 9 ) const;
      
      /// filestream I/O used for domain variables (not recommended because padding creates large storage overhead; use PropertyData instead)
      bool           Out( std::fstream& ) const;
      bool           In( std::fstream& );

      ArrayContainer::const_iterator Begin() const;
      ArrayContainer::const_iterator End()   const;

      /// Julian 17-07-2014 - This is to allow passing of array variables to GEMS.
      double& Front(){ return data_.front(); }

      double NextLargestEntry(double) const;
      bool     HasLargerEntry(double) const;

    private:
      VARIABLE_FLAG  flag_;
      ArrayContainer data_;
    };

  std::ostream&  operator<<( std::ostream& stream, const ArrayVariable& o );


  } // csmp

#endif
