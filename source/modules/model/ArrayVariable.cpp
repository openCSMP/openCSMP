// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ArrayVariable.h"
#include "ScalarVariable.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {

  /// Automatically sets size to corresponding index, flag to ANY
template<uint32_t dim>
ArrayVariable::ArrayVariable( const char* arrayPropertyName,
                              const PropertyDatabase<dim>& pd, double defaultValue, VARIABLE_FLAG flag )
    : data_( pd.StorageKey(arrayPropertyName).dataDepth, defaultValue ), flag_(flag)
    {
    }

template ArrayVariable::ArrayVariable( const char*, const PropertyDatabase<1U>&, double, VARIABLE_FLAG );
template ArrayVariable::ArrayVariable( const char*, const PropertyDatabase<2U>&, double, VARIABLE_FLAG );
template ArrayVariable::ArrayVariable( const char*, const PropertyDatabase<3U>&, double, VARIABLE_FLAG );




void ArrayVariable::CopyValuesOnly( FlaggedArrayVariable& fav ) noexcept
{
    assert ( this->Size() == fav.Size() );
    for ( uint32_t i{0U}; i<fav.Size(); i++ )
      data_[i] = fav(i);
}






/// outputs values to screen up to 10 values per row
void ArrayVariable::Out( long digits ) const
  {
    long   prec(cout.precision(digits));
    size_t pcols(1);

    cout <<"\n\nArrayVariable::Out: array size: "<< Size() <<" values:\n";
    
    if ( digits != 0 ) cout.setf(ios::scientific);

    for ( uint32_t i{0U}; i < Size(); ++i )
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



bool ArrayVariable::Out( const char* filename, long precision ) const
{
  std::ofstream f;
  f.open(filename);
  if( !f.is_open() )
    return false;
  f.precision(precision);
  for (size_t i{0u}; i < Size(); ++i )
      f << data_[i] << std::endl;
    return true;
  }


/**
 @fn  bool ArrayVariable::Out( std::fstream& fp ) const

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
bool ArrayVariable::Out( std::fstream& fp ) const
  {
    if (!fp.is_open())
      {
        std::cerr <<"\nArrayVariable::Out(): invalid file pointer."<< std::endl;
        return false;
      }    

    // flag
    const int8_t flag(flag_);
    fp.write( (char*) &flag, sizeof(int8_t));  // VARIABLE_FLAG

    // size
    const size_t depth( Size() );
    fp.write( (char*) &depth, sizeof(size_t));

    // data
    const size_t bytes(sizeof(double));
    vector<double>::const_iterator dataEnd( data_.end() );
    for ( vector<double>::const_iterator it( data_.begin() ); it != dataEnd; ++it )
      fp.write( (char*) &(*it), bytes);

    return true;
  }



/**
 @fn  bool ArrayVariable::In( std::fstream& fp )

 @brief Reads the array from given file pointer (binary format)

 See bool ArrayVariable::Out( std::fstream& fp ) const

 @author  P. Lang
 @date  9/28/2012

 @param [in,out]  fp  If non-null, the fp.

 @return  true if it succeeds, false if it fails.
 */
bool ArrayVariable::In( std::fstream& fp )
  {
    if (!fp.is_open())
      {
        std::cerr <<"\nArrayVariable::In(): invalid file pointer."<< std::endl;
        return false;
      }

    // flag
    if ( !fp.read( (char*) &flag_, sizeof(int8_t)) ) {  // VARIABLE_FLAG
        std::cerr <<"\nArrayVariable::In(): Not able to read binary record flag"<< std::endl;
        return false;
    }

    // depth
    size_t depth(0);
    if ( !fp.read( (char*) &depth, sizeof(size_t)) ) {
      std::cerr <<"\nArrayVariable::In(): could not read bindary record depth"<< std::endl;
      return false;
      }
    Resize( static_cast<uint32_t>(depth) );

    // data
    const size_t  bytes( sizeof(double) );
    for ( size_t i(0); i < depth; ++i )
      {
        if( !fp.read( (char*) &data_[i], bytes) )
          {
          std::cerr <<"\nArrayVariable::In(): could not read binary data record"<< std::endl;
          return false;
          }
      }
    return true;
  }





ostream&  operator<<( ostream& stream, const ArrayVariable& o )
  {
    stream << "Flag: " << parseStatus(o.Flag()) << " Size: " << o.Size();
    for ( uint32_t i{0U}; i < o.Size(); ++i )
      stream << endl << o[i];

    return stream;
  }

} // end csmp
