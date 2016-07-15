#include "ArrayVariable.h"
#include "ScalarVariable.h"

#include <cassert>

using namespace std;

namespace csmp {

  /// Automatically sets size to corresponding index, flag to ANY
template<size_t dim>
ArrayVariable::ArrayVariable( const char* arrayPropertyName, const PropertyDatabase<dim>& pd, double64 defaultValue, VARIABLE_FLAG flag )
    : data_( pd.StorageKey(arrayPropertyName).dataDepth, defaultValue ), flag_(flag)
    {
    }

template ArrayVariable::ArrayVariable( const char*, const PropertyDatabase<1U>&, double64, VARIABLE_FLAG );
template ArrayVariable::ArrayVariable( const char*, const PropertyDatabase<2U>&, double64, VARIABLE_FLAG );
template ArrayVariable::ArrayVariable( const char*, const PropertyDatabase<3U>&, double64, VARIABLE_FLAG );


  /// As PropertyDatabase ctor, but using index right away
ArrayVariable::ArrayVariable( const Index& arrayKey, double64 defaultValue, VARIABLE_FLAG flag )
    : data_( arrayKey.dataDepth, defaultValue ), flag_(flag)
    {
    }  

ArrayVariable::ArrayVariable()
    : flag_(), data_()
    {
    }

ArrayVariable::ArrayVariable( size_t arraySize, double64 defaultValue, VARIABLE_FLAG flag )
    : flag_(flag), data_( arraySize, defaultValue )
    {

    }

ArrayVariable::ArrayVariable( const ArrayVariable& av )
    : flag_( av.flag_ ), data_( av.data_ )
    {

    }

double64& ArrayVariable::operator()( size_t i )
    {
      assert( i < Size() );
      return data_[i];
    }

double64 ArrayVariable::operator[]( size_t i ) const
    {
      assert( i < Size() );
      return data_[i];
    }

void  ArrayVariable::Component( size_t i, double64 val )
   {
      assert( i < Size() );
      data_[i] = val;
   }

double64  ArrayVariable::Component( size_t i ) const
   {
      assert( i < Size() );
      return data_[i];
   }

size_t ArrayVariable::Size() const
    {
      return data_.size();
    }

void ArrayVariable::Resize( size_t newSize, double64 newValue )
    {
      data_.resize( newSize, newValue );
    }


VARIABLE_FLAG ArrayVariable::Flag(  ) const
    {
      return flag_;
    }

VARIABLE_FLAG& ArrayVariable::Flag(  )
    {
      return flag_;
    }

void ArrayVariable::Flag( VARIABLE_FLAG flag )
    {
      flag_ = flag;
    }

bool ArrayVariable::operator==( const ArrayVariable& av ) const
  {
    if( Flag() != av.Flag() || Size() != av.Size() )
      return false;
    for( size_t i(0); i < Size(); ++i )
      if( data_[i] != av[i] )
        return false;
    return true;
  }


bool ArrayVariable::operator!=( const ArrayVariable& av ) const
  {
    return !(*this == av);
  }


bool ArrayVariable::operator<( const ArrayVariable& av ) const
  {
    if( Size() < av.Size() )
      return true;
    if( Size() > av.Size() )
      return false;
    else
      {
        double64 sumThis(0.), sumParameter(0.);
        for( size_t i(0); i < Size(); ++i )
          {
            sumThis += (*this)[i];
            sumParameter += av[i];
          }
        return sumThis < sumParameter;
      }
    return false;
  }


ArrayVariable& ArrayVariable::operator=( const ArrayVariable& av )
  {
    if( this != &av )
      {
        flag_ = av.flag_;
        data_ = av.data_;
      }
    return *this;
  }

void ArrayVariable::CopyValuesOnly( FlaggedArrayVariable& fav )
{
    if( this->Size() == fav.Size() )
        for (size_t i = 0 ; i< fav.Size();i++)
            data_[i]         = fav(i);
}

ArrayVariable& ArrayVariable::operator=( double64 val )
  {
  for( size_t i(0); i < Size(); ++i )
    data_[i] = val;
  return *this;
  }

ArrayVariable& ArrayVariable::operator=( const ScalarVariable& val )
  {
  for( size_t i(0); i < Size(); ++i )
    data_[i] = val.Value();
  return *this;
  }

/// Standart Operations with temporary object

ArrayVariable ArrayVariable::operator+( double64 val ) const
  {
    ArrayVariable temp_arr( *this );
    for( size_t i(0); i < Size(); ++i )
        temp_arr(i) += val;
    return temp_arr;
}

ArrayVariable ArrayVariable::operator-( double64 val ) const
  {
    ArrayVariable temp_arr( *this );
    for( size_t i(0); i < Size(); ++i )
        temp_arr(i) -= val;
    return temp_arr;
  }

ArrayVariable ArrayVariable::operator*( double64 val ) const
  {
    ArrayVariable temp_arr( *this );
    for( size_t i(0); i < Size(); ++i )
        temp_arr(i) *= val;
    return temp_arr;
  }

ArrayVariable ArrayVariable::operator/( double64 val ) const
  {
    ArrayVariable temp_arr( *this );
    for( size_t i(0); i < Size(); ++i )
        temp_arr(i) /= val;
    return temp_arr;
  }

ArrayVariable ArrayVariable::operator^( double64 val ) const
  {
    ArrayVariable temp_arr( *this );
    for( size_t i(0); i < Size(); ++i )
        temp_arr(i) = std::pow(temp_arr(i), val);
    return temp_arr;
  }

/// Standart Operations with current object

ArrayVariable& ArrayVariable::operator+=( double64 val )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] += val;
    return *this;
  }

ArrayVariable& ArrayVariable::operator-=( double64 val )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] -= val;
    return *this;
  }

ArrayVariable& ArrayVariable::operator*=( double64 val )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] *= val;
    return *this;
  }

ArrayVariable& ArrayVariable::operator/=( double64 val )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] /= val;
    return *this;
  }

/// Standart Operations with current object

/// Operations with ScalarVariables
ArrayVariable& ArrayVariable::operator+=( const ScalarVariable& sc )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] += sc.Value();
    return *this;
  }

ArrayVariable& ArrayVariable::operator-=( const ScalarVariable& sc )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] -= sc.Value();
    return *this;
  }

ArrayVariable& ArrayVariable::operator*=( const ScalarVariable& sc )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] *= sc.Value();
    return *this;
  }

ArrayVariable& ArrayVariable::operator/=( const ScalarVariable& sc )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] /= sc.Value();
    return *this;
  }

/// Operations with ArrayVariables
ArrayVariable& ArrayVariable::operator+=( const ArrayVariable& arr )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] += arr[i];
    return *this;
  }

ArrayVariable& ArrayVariable::operator-=( const ArrayVariable& arr )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] -= arr[i];
    return *this;
  }

ArrayVariable& ArrayVariable::operator*=( const ArrayVariable& arr )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] *= arr[i];
    return *this;
  }

ArrayVariable& ArrayVariable::operator/=( const ArrayVariable& arr )
  {
    for( size_t i(0); i < Size(); ++i )
        data_[i] /= arr[i];
    return *this;
  }
ArrayVariable ArrayVariable::operator-( const ArrayVariable& av ) const
  {
    ArrayVariable returnArray( *this );
    for (size_t i(0); i < Size(); ++i )
      returnArray(i) -=  av[i];
    return returnArray;
  }

ArrayVariable ArrayVariable::operator+( const ArrayVariable& av ) const
  {
    ArrayVariable returnArray( *this );
    for (size_t i(0); i < Size(); ++i )
      returnArray(i) +=  av[i];
    return returnArray;
  }

ArrayVariable ArrayVariable::operator*( const ArrayVariable& av ) const
  {
    ArrayVariable returnArray( *this );
    for (size_t i(0); i < Size(); ++i )
      returnArray(i) *=  av[i];
    return returnArray;
  }


ArrayVariable ArrayVariable::operator/( const ArrayVariable& av ) const
  {
    ArrayVariable returnArray( *this );
    for (size_t i(0); i < Size(); ++i )
      returnArray(i) /=  av[i];
    return returnArray;
  }

bool ArrayVariable::IsWithinRange( double64 min, double64 max ) const
  {
    assert( min < max );

    for ( size_t i(0); i < Size(); ++i )
      if( ((*this)[i] < min) || ((*this)[i] > max) )
        return false;
    return true;
  }


/// returns the minimum and maximum of the values stored in the array variable 
void  ArrayVariable::MinMax( double64& min, double64& max ) const
 {
     min = (*min_element( data_.begin(), data_.end() ));
     max = (*max_element( data_.begin(), data_.end() ));
 }



void ArrayVariable::Fabs()
  {
  for ( size_t i(0); i < Size(); ++i )
      data_[i] = fabs(data_[i]);
  }

void ArrayVariable::Sqrt()
  {
  for ( size_t i(0); i < Size(); ++i )
    data_[i] = sqrt(data_[i]);
  }

void ArrayVariable::Log10()
  {
  for ( size_t i(0); i < Size(); ++i )
    data_[i] = log10(data_[i]);
  }

void ArrayVariable::Ln()
  {
  for ( size_t i(0); i < Size(); ++i )
    data_[i] = log(data_[i]);
  }





/// outputs values to screen up to 10 values per row
void ArrayVariable::Out( long digits ) const
  {
    long   prec(cout.precision(digits));
    size_t pcols(1);

    cout <<"\n\nArrayVariable::Out: array size: "<< Size() <<" values:\n";
    
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



bool ArrayVariable::Out( const char* filename, size_t precision ) const
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
 @fn  bool ArrayVariable::Out( std::FILE* fp ) const

 @brief Outs the array to given file pointer (binary format)

 Order:
 1. Flag
 2. Size
 3. Data

 @author  P. Lang
 @date  9/28/2012

 @param [in,out]  fp  If non-null, the fp.

 @return  true if it succeeds, false if it fails.
 */
bool ArrayVariable::Out( std::FILE* fp ) const
  {
    if (!fp)
      {
        std::cout <<"\nArrayVariable::Out(): invalid file pointer."<< std::endl;
        return false;
      }    

    // flag
    const VARIABLE_FLAG flag(flag_);
    fwrite( (void*) &flag, sizeof(VARIABLE_FLAG), 1, fp );

    // size
    const size_t depth( Size() );
    fwrite( (void*) &depth, sizeof(size_t), 1, fp );

    // data
    const size_t bytes(sizeof(double64));
    vector<double>::const_iterator dataEnd( data_.end() );
    for ( vector<double>::const_iterator it( data_.begin() ); it != dataEnd; ++it )
      std::fwrite( (void*) &(*it), bytes, 1, fp );

    return true;
  }


/**
 @fn  bool ArrayVariable::In( std::FILE* fp )

 @brief Reads the array from given file pointer (binary format)

 See bool ArrayVariable::Out( std::FILE* fp ) const

 @author  P. Lang
 @date  9/28/2012

 @param [in,out]  fp  If non-null, the fp.

 @return  true if it succeeds, false if it fails.
 */
bool ArrayVariable::In( std::FILE* fp )
  {
    if (!fp)
      {
        std::cout <<"\nArrayVariable::In(): invalid file pointer."<< std::endl;
        return false;
      }

    // flag
    if ( !std::fread( (void*) &flag_, sizeof(VARIABLE_FLAG), 1, fp ) ) {
        std::cout <<"\nArrayVariable::In(): Not able to read binary record flag"<< std::endl;
        return false;
    }

    // depth
    size_t depth(0);
    if ( !std::fread( (void*) &depth, sizeof(size_t), 1, fp ) ) {
      std::cout <<"\nArrayVariable::In(): could not read bindary record depth"<< std::endl;
      return false;
      }
    Resize(depth);

    // data
    const size_t  bytes( sizeof(double64) );
    for ( size_t i(0); i < depth; ++i )
      {
        if( !std::fread( (void*) &data_[i], bytes, 1, fp ) )
          {
          std::cout <<"\nArrayVariable::In(): could not read binary data record"<< std::endl;
          return false;
          }
      }
    return true;
  }


namespace{
  ArrayVariable::ArrayContainer sortedArray( ArrayVariable::ArrayContainer::const_iterator plainDataBegin,
                                             ArrayVariable::ArrayContainer::const_iterator plainDataEnd )
  {
    ArrayVariable::ArrayContainer sortedData( plainDataBegin, plainDataEnd );
    sort( sortedData.begin(), sortedData.end() );
    return sortedData;
  }
} // locally restricted


double ArrayVariable::NextLargestEntry( double64 fromValue ) const
{
  ArrayContainer sorted( sortedArray( data_.begin(), data_.end() ) );
  return *upper_bound( sorted.begin(), sorted.end(), fromValue );
}


bool ArrayVariable::HasLargerEntry( double64 fromValue ) const
{
  ArrayContainer sorted( sortedArray( data_.begin(), data_.end() ) );
  return upper_bound( sorted.begin(), sorted.end(), fromValue ) != sorted.end();
}


void ArrayVariable::Sort()
{
  sort( data_.begin(), data_.end() );
}


ArrayVariable::ArrayContainer::const_iterator ArrayVariable::Begin() const
{
  return data_.begin();
}


ArrayVariable::ArrayContainer::const_iterator ArrayVariable::End() const
{
  return data_.end();
}


ostream&  operator<<( ostream& stream, const ArrayVariable& o )
  {
    stream << "Flag: " << parseStatus(o.Flag()) << " Size: " << o.Size();
    for (size_t i(0); i < o.Size(); ++i )
      stream << endl << o[i];

    return stream;
  }

} // end csmp
