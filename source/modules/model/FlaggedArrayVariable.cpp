#include "FlaggedArrayVariable.h"


#include "ScalarVariable.h"
#include "ArrayVariable.h"
#include "PropertyDatabase.h"

#include <numeric>

using namespace std;

namespace csmp {

// ============================================================================
//  Template constructor — explicit instantiations
// ============================================================================

template<uint32_t dim>
FlaggedArrayVariable::FlaggedArrayVariable(
    const char*                  arrayPropertyName,
    const PropertyDatabase<dim>& pd,
    double                       defaultValue,
    VARIABLE_FLAG                flag )
    : flags_( static_cast<size_t>(
                  pd.StorageKey( arrayPropertyName ).dataDepth ), flag ),
      data_(  static_cast<size_t>(
                  pd.StorageKey( arrayPropertyName ).dataDepth ), defaultValue )
{}

template FlaggedArrayVariable::FlaggedArrayVariable(
    const char*, const PropertyDatabase<1U>&, double, VARIABLE_FLAG );
template FlaggedArrayVariable::FlaggedArrayVariable(
    const char*, const PropertyDatabase<2U>&, double, VARIABLE_FLAG );
template FlaggedArrayVariable::FlaggedArrayVariable(
    const char*, const PropertyDatabase<3U>&, double, VARIABLE_FLAG );

// ============================================================================
//  Output
// ============================================================================

void FlaggedArrayVariable::Out( long digits ) const
{
    long   prec( cout.precision( digits ) );
    size_t pcols( 1 );

    cout << "\n\nFlaggedArrayVariable::Out: array size: "
         << Size() << " values:\n";

    if ( digits != 0 ) cout.setf( ios::scientific );

    for ( size_t i = 0; i < Size(); ++i )
    {
        cout << "(" << i << "):  " << data_[i]
             << " [" << static_cast<int>( flags_[i] ) << "],  ";
        if ( pcols == 10 ) { cout << "\n"; pcols = 0; }
        ++pcols;
    }
    cout << "\n";

    if ( digits != 0 )
    {
        cout.unsetf( ios::scientific );
        cout.precision( prec );
    }
}



void FlaggedArrayVariable::CopyValuesOnly( ArrayVariable& av ) noexcept
{
    assert( Size() == av.Size() );
    for ( uint32_t i = 0; i < Size(); ++i )
        data_[i] = av(i);
}



bool FlaggedArrayVariable::Out( const char* filename,
                                 size_t      precision ) const
{
    ofstream f( filename );
    if ( !f ) return false;
    f.precision( static_cast<long>(precision) );
    for ( size_t i = 0; i < Size(); ++i )
        f << data_[i] << "\n";
    return true;
}

bool FlaggedArrayVariable::Out( fstream& fp ) const
{
    if ( !fp )
    {
        cerr << "\nFlaggedArrayVariable::Out(): invalid file pointer.\n";
        return false;
    }

    const size_t depth( Size() );
    fp.write( reinterpret_cast<const char*>( &depth ), sizeof(size_t) );

    const size_t flag_size( sizeof(int8_t) );
    for ( const auto& f : flags_ )
        fp.write( reinterpret_cast<const char*>( &f ), flag_size );

    const size_t bytes( sizeof(double) );
    for ( const auto& v : data_ )
        fp.write( reinterpret_cast<const char*>( &v ), bytes );

    return true;
}

bool FlaggedArrayVariable::In( fstream& fp )
{
    if ( !fp )
    {
        cerr << "\nFlaggedArrayVariable::In(): invalid file pointer.\n";
        return false;
    }

    size_t depth( 0 );
    if ( !fp.read( reinterpret_cast<char*>( &depth ), sizeof(size_t) ) )
    {
        cerr << "\nFlaggedArrayVariable::In(): could not read depth.\n";
        return false;
    }
    Resize( static_cast<uint32_t>( depth ) );

    const size_t flag_size( sizeof(int8_t) );
    for ( size_t i = 0; i < depth; ++i )
        if ( !fp.read( reinterpret_cast<char*>( &flags_[i] ), flag_size ) )
        {
            cerr << "\nFlaggedArrayVariable::In(): could not read flags.\n";
            return false;
        }

    const size_t bytes( sizeof(double) );
    for ( size_t i = 0; i < depth; ++i )
        if ( !fp.read( reinterpret_cast<char*>( &data_[i] ), bytes ) )
        {
            cerr << "\nFlaggedArrayVariable::In(): could not read data.\n";
            return false;
        }

    return true;
}

// ============================================================================
//  Stream output
// ============================================================================

ostream& operator<<( ostream& stream, const FlaggedArrayVariable& o )
{
    stream << " Size: " << o.Size();
    for ( size_t i = 0; i < o.Size(); ++i )
        stream << "\n  [" << i << "] Flag: "
               << static_cast<int>( o.Flag(i) )
               << "  Value: " << o[i];
    return stream;
}

} // namespace csmp

