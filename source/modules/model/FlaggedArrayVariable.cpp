#include "FlaggedArrayVariable.h"
#include "ScalarVariable.h"

using namespace std;

namespace csmp{

/// Default constructor
FlaggedArrayVariable::FlaggedArrayVariable()
  : flags_(ANY),
    data_(0)
  {
  }

/// As PropertyDatabase constructor, but using index right away
FlaggedArrayVariable::FlaggedArrayVariable( const Index& arrayKey,
                                            double defaultValue,
                                            VARIABLE_FLAG flag )
  : data_ ( arrayKey.dataDepth, defaultValue ),
    flags_( arrayKey.dataDepth, flag)
  {
  }

FlaggedArrayVariable::FlaggedArrayVariable( unsigned int arraySize,
                                            double defaultValue,
                                            VARIABLE_FLAG flag )
  : data_ ( arraySize, defaultValue ),
    flags_( arraySize, flag)
  {
  }

FlaggedArrayVariable::FlaggedArrayVariable( const FlaggedArrayVariable& av )
: data_       ( av.data_ ),
  flags_      ( av.flags_)
  {
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator=( const FlaggedArrayVariable& av )
  {
    if( this != &av )
      {
        data_         = av.data_;
        flags_        = av.flags_;
      }
    return *this;
  }

void FlaggedArrayVariable::CopyValuesOnly( ArrayVariable& av )
{
    if( this->Size() == av.Size() )
        for (auto i = 0 ; i< av.Size();i++)
            data_[i]         = av(i);
}

FlaggedArrayVariable& FlaggedArrayVariable::operator=( double val )
  {
    for ( auto& i : data_ ) i = val;
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator=( const ScalarVariable& val )
  {
    for ( auto& i : data_ ) i = val();
    return *this;
  }


uint32_t FlaggedArrayVariable::Size() const
  {
    return static_cast<uint32_t>(data_.size());
  }

void FlaggedArrayVariable::Resize( uint32_t newSize, double newValue )
  {
    data_.resize ( newSize, newValue );
    flags_.resize( newSize, ANY );
  }


VARIABLE_FLAG FlaggedArrayVariable::Flag( const size_t& i ) const
  {
    return flags_[i];
  }

VARIABLE_FLAG& FlaggedArrayVariable::Flag( const size_t& i )
  {
    return flags_[i];
  }

void FlaggedArrayVariable::Flag( const size_t& i, VARIABLE_FLAG flag )
  {
    flags_[i] = flag;
  }

double FlaggedArrayVariable::operator[]( size_t i ) const
  {
    assert( i < Size() );
    return data_[i];
  }

double& FlaggedArrayVariable::operator()( size_t i )
  {
    assert( i < Size() );
    return data_[i];
  }

void  FlaggedArrayVariable::Component( uint32_t i, double val )
 {
    assert( i < Size() );
    data_[i] = val;
 }

double  FlaggedArrayVariable::Component( uint32_t i ) const
 {
    assert( i < Size() );
    return data_[i];
 }

bool FlaggedArrayVariable::operator==( const FlaggedArrayVariable& av ) const
  {
    assert( av.Size() == Size() );
    for( uint32_t i{0U}; i < Size(); ++i )
        if( (*this)[i] != av[i] )
          return false;
      return true;
  }

bool FlaggedArrayVariable::operator!=( const FlaggedArrayVariable& av ) const
  {
    return !(*this == av);
  }

bool FlaggedArrayVariable::operator<( const FlaggedArrayVariable& av ) const
  {
    assert( av.Size() == Size() );

    for( uint32_t i{0U}; i < Size(); ++i )
      if( (*this)[i] >= av[i] )
        return false;
    return true;
  }

bool FlaggedArrayVariable::operator<=( const FlaggedArrayVariable& av ) const
  {
    assert( av.Size() == Size() );

    for( uint32_t i{0U}; i < Size(); ++i )
      if( (*this)[i] > av[i] )
        return false;
    return true;
  }

bool FlaggedArrayVariable::operator>( const FlaggedArrayVariable& av ) const
  {
    assert( av.Size() == Size() );

    for( uint32_t i{0U}; i < Size(); ++i )
      if( (*this)[i] <= av[i] )
        return false;
    return true;
  }

bool FlaggedArrayVariable::operator>=( const FlaggedArrayVariable& av ) const
  {
    assert( av.Size() == Size() );

    for( uint32_t i{0U}; i < Size(); ++i )
      if( (*this)[i] < av[i] )
        return false;
    return true;
  }

/// Standard Operations with temporary object

FlaggedArrayVariable FlaggedArrayVariable::operator+( double val ) const
  {
    FlaggedArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) += val;
    return temp_arr;
}

FlaggedArrayVariable FlaggedArrayVariable::operator-( double val ) const
  {
    FlaggedArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) -= val;

    return temp_arr;
  }

FlaggedArrayVariable FlaggedArrayVariable::operator*( double val ) const
  {
    FlaggedArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) *= val;
    return temp_arr;
  }

FlaggedArrayVariable FlaggedArrayVariable::operator/( double val ) const
  {
    FlaggedArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
        temp_arr(i) /= val;
    return temp_arr;
  }

FlaggedArrayVariable FlaggedArrayVariable::operator^( double val ) const
  {
    FlaggedArrayVariable temp_arr( *this );
    for( uint32_t i{0U}; i < Size(); ++i )
    {
        temp_arr(i) = std::pow( temp_arr[i], val);
        temp_arr.Flag(i) = flags_[i];
    }
    return temp_arr;
  }

/// Standard Operations with current object

FlaggedArrayVariable& FlaggedArrayVariable::operator+=( double val )
  {
    for ( auto& i : data_ ) i += val;
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator-=( double val )
  {
    for ( auto& i : data_ ) i -= val;
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator*=( double val )
  {
    for ( auto& i : data_ ) i *= val;
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator/=( double val )
  {
    for ( auto& i : data_ ) i /= val;
    return *this;
  }

/// Standard Operations with current object

/// Operations with ScalarVariables
FlaggedArrayVariable& FlaggedArrayVariable::operator+=( const ScalarVariable& sc )
  {
    for ( auto& i : data_ ) i += sc();
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator-=( const ScalarVariable& sc )
  {
    for ( auto& i : data_ ) i -= sc();
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator*=( const ScalarVariable& sc )
  {
    for ( auto& i : data_ ) i *= sc();
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator/=( const ScalarVariable& sc )
  {
    for ( auto& i : data_ ) i /= sc();
    return *this;
  }

/// Operations with FlaggedArrayVariables
FlaggedArrayVariable& FlaggedArrayVariable::operator+=( const FlaggedArrayVariable& av )
  {
    for( uint32_t i{0U}; i < Size(); ++i )
        data_[i] += av[i];
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator-=( const FlaggedArrayVariable& av )
  {
    for( uint32_t i{0U}; i < Size(); ++i )
        data_[i] -= av[i];
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator*=( const FlaggedArrayVariable& av )
  {
    for( uint32_t i{0U}; i < Size(); ++i )
        data_[i] *= av[i];
    return *this;
  }

FlaggedArrayVariable& FlaggedArrayVariable::operator/=( const FlaggedArrayVariable& av )
  {
    for( uint32_t i{0U}; i < Size(); ++i )
        data_[i] /= av[i];
    return *this;
  }

FlaggedArrayVariable FlaggedArrayVariable::operator-( const FlaggedArrayVariable& av ) const
  {
    FlaggedArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
    {
        returnArray(i) -= av[i];
        returnArray.Flag(i) = flags_[i];
    }
    return returnArray;
  }

FlaggedArrayVariable FlaggedArrayVariable::operator+( const FlaggedArrayVariable& av ) const
  {
    FlaggedArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
    {
        returnArray(i) += av[i];
        returnArray.Flag(i) = flags_[i];
    }
    return returnArray;
  }

FlaggedArrayVariable FlaggedArrayVariable::operator*( const FlaggedArrayVariable& av ) const
  {
    FlaggedArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
    {
      returnArray(i) *= av[i];
      returnArray.Flag(i) = flags_[i];
    }
    return returnArray;
  }


FlaggedArrayVariable FlaggedArrayVariable::operator/( const FlaggedArrayVariable& av ) const
  {
    FlaggedArrayVariable returnArray( *this );
    for ( uint32_t i{0U}; i < Size(); ++i )
    {
      returnArray(i) /= av[i];
      returnArray.Flag(i) = flags_[i];
    }
    return returnArray;
  }

bool FlaggedArrayVariable::IsWithinRange( double min, double max ) const
  {
    assert( min < max );

    for ( auto& i : data_ )
      if( i < min || i > max )
        return false;
    return true;
  }


/// returns the minimum and maximum of the values stored in the array variable 
void  FlaggedArrayVariable::MinMax( double& min, double& max ) const
  {
     min = (*min_element( data_.begin(), data_.end() ));
     max = (*max_element( data_.begin(), data_.end() ));
  }






/// outputs values to screen up to 10 values per row
void FlaggedArrayVariable::Out( long digits ) const
  {
    long   prec(cout.precision(digits));
    size_t pcols(1);

    cout <<"\n\nFlaggedArrayVariable::Out: array size: "<< Size() <<" values:\n";

    if ( digits != 0 ) cout.setf(ios::scientific);

    for ( size_t i(0); i < Size(); ++i )
      {
         cout <<"("<< i <<"):  "<< (*this)[i] <<", ";
         if ( pcols == 10 ) {
              cout << endl;
              pcols = 0;
           }
         pcols++;
      }
    cout << endl;

    if ( digits != 0 ) {
         cout.unsetf( ios::scientific );
         cout.precision(prec);
      }

 } // end Out



bool FlaggedArrayVariable::Out( const char* filename, size_t precision ) const
  {
    std::ofstream f;
    f.open(filename);
    if( !f.is_open() )
      return false;
    f.precision(precision);
    for (size_t i(0); i < Size(); ++i )
      f << data_[i] << std::endl;
    return true;
  }


/**
 @fn  bool FlaggedArrayVariable::Out( std::fstream& fp ) const

 @brief Outs the array to given file pointer (binary format)

 Order:
 1. Flag
 2. Size
 3. Data

 Roman, 2013: Changed order due to multiple flags:
 1. Size
 2. Flags
 3. Data

 @author  P. Lang
 @date  9/28/2012

 @param [in,out]  fp  If non-null, the fp.

 @return  true if it succeeds, false if it fails.
 */
bool FlaggedArrayVariable::Out( std::fstream& fp ) const
  {
    if (!fp)
      {
        std::cerr <<"\nFlaggedArrayVariable::Out(): invalid file pointer."<< std::endl;
        return false;
      }    

    // size
    const size_t depth( Size() );
    fp.write( (char*) &depth, sizeof(size_t));

    // flags
    const size_t flag_size(sizeof(int8_t));  // VARIABLE_FLAG
    vector<VARIABLE_FLAG>::const_iterator flagsEnd( flags_.end() );
    for ( vector<VARIABLE_FLAG>::const_iterator it( flags_.begin() ); it != flagsEnd; ++it )
      fp.write( (char*) &(*it), flag_size);

    // data
    const size_t bytes(sizeof(double));
    vector<double>::const_iterator dataEnd( data_.end() );
    for ( vector<double>::const_iterator it( data_.begin() ); it != dataEnd; ++it )
      fp.write( (char*) &(*it), bytes);

    return true;
  }
  
  


/**
 @fn  bool FlaggedArrayVariable::In( std::fstream& fp )

 @brief Reads the array from given file pointer (binary format)

 See bool FlaggedArrayVariable::Out( std::fstream& fp ) const

 @param [in,out]  fp  If non-null, the fp.

 @return  true if it succeeds, false if it fails.
 */
bool FlaggedArrayVariable::In( std::fstream& fp )
  {
    if (!fp)
      {
        std::cout <<"\nFlaggedArrayVariable::In(): invalid file pointer."<< std::endl;
        return false;
      }

    // depth
    size_t depth(0);
    if ( !fp.read( (char*) &depth, sizeof(size_t)) ) {
      std::cerr <<"\nFlaggedArrayVariable::In(): could not read binary record depth"<< std::endl;
      return false;
      }
    Resize( static_cast<uint32_t>(depth) );

    // flags
    const size_t  flags_size( sizeof( int8_t ) );  // VARIABLE_FLAG
    for ( size_t i(0); i < depth; ++i )
      {
        if( !fp.read( (char*) &flags_[i], flags_size) )
          {
          std::cerr <<"\nFlaggedArrayVariable::In(): could not read binary flags record"<< std::endl;
          return false;
          }
      }

    // data
    const size_t  bytes( sizeof(double) );
    for ( size_t i(0); i < depth; ++i )
      {
        if( !fp.read( (char*) &data_[i], bytes) )
          {
          std::cerr <<"\nFlaggedArrayVariable::In(): could not read binary data record"<< std::endl;
          return false;
          }
      }
    return true;
  }




namespace{
  FlaggedArrayVariable::FlaggedArrayContainer sortedArray( FlaggedArrayVariable::FlaggedArrayContainer::const_iterator plainDataBegin,
                                                           FlaggedArrayVariable::FlaggedArrayContainer::const_iterator plainDataEnd )
  {
    FlaggedArrayVariable::FlaggedArrayContainer sortedData( plainDataBegin, plainDataEnd );
    sort( sortedData.begin(), sortedData.end() );
    return sortedData;
  }
} // locally restricted


double FlaggedArrayVariable::NextLargestEntry( double fromValue ) const
{
  FlaggedArrayContainer sorted( sortedArray( data_.begin(), data_.end() ) );
  return *upper_bound( sorted.begin(), sorted.end(), fromValue );
}


bool FlaggedArrayVariable::HasLargerEntry( double fromValue ) const
{
  FlaggedArrayContainer sorted( sortedArray( data_.begin(), data_.end() ) );
  return upper_bound( sorted.begin(), sorted.end(), fromValue ) != sorted.end();
}


void FlaggedArrayVariable::Sort()
{
  sort( data_.begin(), data_.end() );
}


FlaggedArrayVariable::FlaggedArrayContainer::const_iterator FlaggedArrayVariable::Begin() const
{
  return data_.begin();
}


FlaggedArrayVariable::FlaggedArrayContainer::const_iterator FlaggedArrayVariable::End() const
{
  return data_.end();
}


ostream&  operator<<( ostream& stream, const FlaggedArrayVariable& o )
  {
    stream <<" Size: " << o.Size();
    for (size_t i(0); i < o.Size(); ++i )
      stream << endl << "Flag: " << o.Flag(i) << "Data = "<<o[i];

    return stream;
  }

  } // csmp
