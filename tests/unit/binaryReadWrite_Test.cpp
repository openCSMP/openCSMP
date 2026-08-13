//
//  binaryReadWrite_Test.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 7/04/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "binaryReadWrite_Test.h"

using namespace std;

namespace csmp {

// ============================================================================
// Helpers
// ============================================================================

std::fstream binaryReadWrite_Test::openWrite( const std::string& name )
{
    std::fstream fp( name, std::ios::out | std::ios::binary | std::ios::trunc );
    if ( !fp.is_open() )
        throw std::runtime_error(
            "binaryReadWrite_Test: could not open '" + name + "' for writing." );
    return fp;
}

std::fstream binaryReadWrite_Test::openRead( const std::string& name )
{
    std::fstream fp( name, std::ios::in | std::ios::binary );
    if ( !fp.is_open() )
        throw std::runtime_error(
            "binaryReadWrite_Test: could not open '" + name + "' for reading." );
    return fp;
}


// ============================================================================
// RAII helper — silences std::cerr for the duration of a scope, then
// restores it.  Used to suppress expected diagnostic output from
// deliberate error-path tests so that the test log stays clean.
// ============================================================================

struct SilenceCerr {
    SilenceCerr()
        : old_buf_( std::cerr.rdbuf( sink_.rdbuf() ) )
    {}
    ~SilenceCerr()
    {
        std::cerr.rdbuf( old_buf_ );
    }
    // Non-copyable.
    SilenceCerr( const SilenceCerr& )            = delete;
    SilenceCerr& operator=( const SilenceCerr& ) = delete;
  private:
    std::ostringstream sink_;
    std::streambuf*    old_buf_;
};


// ============================================================================
// Main test runner
// ============================================================================

void binaryReadWrite_Test::run()
{
    // =========================================================================
    // BinaryFileSectionWrite / BinaryFileSectionRead
    // =========================================================================
    _info("binaryReadWrite_Test: testing BinaryFileSectionWrite / BinaryFileSectionRead:");

    // --- Normal round-trip ---
    {
        auto fp = openWrite("brw_section.bin");
        BinaryFileSectionWrite( fp, "VSETCORD" );
        BinaryFileSectionWrite( fp, "VSETPELT" );
        fp.close();

        auto fp2 = openRead("brw_section.bin");
        bool ok = true;
        try {
            BinaryFileSectionRead( fp2, "VSETCORD" );
            BinaryFileSectionRead( fp2, "VSETPELT" );
        } catch ( const std::exception& e ) {
            ok = false;
            if ( verbose_ )
                std::cerr << "  unexpected exception: " << e.what() << "\n";
        }
        _test( ok );
    }

    // --- Mismatched tag throws ---
    // Clear errno first so that the BinaryFileSectionRead diagnostic does not
    // print a spurious "No such file or directory" from a prior OS call.
    {
        errno = 0;
        auto fp = openWrite("brw_section_bad.bin");
        BinaryFileSectionWrite( fp, "VSETCORD" );
        fp.close();

        auto fp2 = openRead("brw_section_bad.bin");
        bool threw = false;
        {
            SilenceCerr sc; // suppress the errno diagnostic printed before throw
            try {
                BinaryFileSectionRead( fp2, "VSETPELT" ); // wrong tag
            } catch ( const std::runtime_error& ) {
                threw = true;
            }
        }
        _test( threw );
    }

    // --- Header too long throws on write ---
    {
        auto fp = openWrite("brw_section_long.bin");
        bool threw = false;
        try {
            BinaryFileSectionWrite( fp, "TOOLONGHEADER" );
        } catch ( const std::invalid_argument& ) {
            threw = true;
        }
        _test( threw );
    }

    // --- Header too long throws on read ---
    {
        auto fp = openWrite("brw_section_longr.bin");
        BinaryFileSectionWrite( fp, "VSETCORD" );
        fp.close();

        auto fp2 = openRead("brw_section_longr.bin");
        bool threw = false;
        try {
            BinaryFileSectionRead( fp2, "TOOLONGHEADER" );
        } catch ( const std::invalid_argument& ) {
            threw = true;
        }
        _test( threw );
    }

    // --- Write to closed stream throws ---
    {
        std::fstream fp; // deliberately not opened
        bool threw = false;
        try {
            BinaryFileSectionWrite( fp, "VSETCORD" );
        } catch ( const std::runtime_error& ) {
            threw = true;
        }
        _test( threw );
    }


    // =========================================================================
    // readContainerSize
    // =========================================================================
    _info("binaryReadWrite_Test: testing readContainerSize:");

    // --- Normal read ---
    {
        auto fp = openWrite("brw_csize.bin");
        const size_t n = 42U;
        fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
        fp.close();

        auto fp2 = openRead("brw_csize.bin");
        _test( readContainerSize( fp2 ) == 42U );
    }

    // --- Zero is valid ---
    {
        auto fp = openWrite("brw_csize_zero.bin");
        const size_t n = 0U;
        fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
        fp.close();

        auto fp2 = openRead("brw_csize_zero.bin");
        _test( readContainerSize( fp2 ) == 0U );
    }

    // --- SIZE_MAX sentinel throws ---
    {
        auto fp = openWrite("brw_csize_max.bin");
        const size_t n = std::numeric_limits<size_t>::max();
        fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
        fp.close();

        auto fp2 = openRead("brw_csize_max.bin");
        bool threw = false;
        try {
            readContainerSize( fp2 );
        } catch ( const std::runtime_error& ) {
            threw = true;
        }
        _test( threw );
    }

    // --- Truncated file (EOF before full size_t) throws ---
    {
        auto fp = openWrite("brw_csize_trunc.bin");
        const uint16_t partial = 1U;
        fp.write( reinterpret_cast<const char*>(&partial), sizeof(uint16_t) );
        fp.close();

        auto fp2 = openRead("brw_csize_trunc.bin");
        uint16_t dummy;
        fp2.read( reinterpret_cast<char*>(&dummy), sizeof(uint16_t) );
        bool threw = false;
        try {
            readContainerSize( fp2 ); // stream is now at EOF
        } catch ( const std::runtime_error& ) {
            threw = true;
        }
        _test( threw );
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — std::string
    // =========================================================================
    _info("binaryReadWrite_Test: testing string round-trips:");

    // --- Normal string ---
    {
        const std::string original = "Hello, CSMP binary world!";
        auto fp = openWrite("brw_string.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::string result;
        auto fp2 = openRead("brw_string.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Empty string ---
    {
        const std::string original;
        auto fp = openWrite("brw_string_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::string result = "should be cleared";
        auto fp2 = openRead("brw_string_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.empty() );
    }

    // --- String with embedded spaces and special characters ---
    {
        const std::string original = "permeability_xx [m^2]\t3D tensor";
        auto fp = openWrite("brw_string_special.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::string result;
        auto fp2 = openRead("brw_string_special.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Multiple strings in sequence ---
    {
        const std::vector<std::string> originals = {
            "pressure", "temperature", "velocity", "permeability"
        };
        auto fp = openWrite("brw_string_multi.bin");
        for ( const auto& s : originals )
            _test( binaryFileWrite( fp, s ) );
        fp.close();

        auto fp2 = openRead("brw_string_multi.bin");
        for ( const auto& expected : originals ) {
            std::string result;
            _test( binaryFileRead( fp2, result ) );
            _test( result == expected );
        }
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        _test( binaryFileWrite( fp, std::string("test") ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        std::string s;
        _test( binaryFileRead( fp, s ) == false );
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — const char* / char[]
    // =========================================================================
    _info("binaryReadWrite_Test: testing const char* / char[] round-trips:");

    // --- Normal C-string ---
    {
        const char* original = "VSETCORD";
        auto fp = openWrite("brw_cstr.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        char result[INFO_STRING] = {};
        auto fp2 = openRead("brw_cstr.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( std::strcmp( result, original ) == 0 );
    }

    // --- Empty C-string ---
    {
        const char* original = "";
        auto fp = openWrite("brw_cstr_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        char result[INFO_STRING] = {};
        auto fp2 = openRead("brw_cstr_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( std::strlen( result ) == 0 );
    }

    // --- nullptr returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        auto fp = openWrite("brw_cstr_null.bin");
        _test( binaryFileWrite( fp, static_cast<const char*>(nullptr) ) == false );
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        _test( binaryFileWrite( fp, "test" ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        char buf[INFO_STRING] = {};
        _test( binaryFileRead( fp, buf ) == false );
    }

    // --- Multiple C-strings in sequence ---
    {
        const std::vector<std::string> names = {
            "pressure", "saturation", "porosity"
        };
        auto fp = openWrite("brw_cstr_multi.bin");
        for ( const auto& n : names )
            _test( binaryFileWrite( fp, n.c_str() ) );
        fp.close();

        auto fp2 = openRead("brw_cstr_multi.bin");
        for ( const auto& expected : names ) {
            char result[INFO_STRING] = {};
            _test( binaryFileRead( fp2, result ) );
            _test( std::string(result) == expected );
        }
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — std::vector<T>
    // =========================================================================
    _info("binaryReadWrite_Test: testing vector<T> round-trips:");

    // --- vector<double> including infinity and NaN ---
    {
        const std::vector<double> original = {
            1.0, 2.5, -3.14, 0.0,
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::quiet_NaN()
        };
        auto fp = openWrite("brw_vec_double.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::vector<double> result;
        auto fp2 = openRead("brw_vec_double.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.size() == original.size() );
        for ( size_t i = 0; i < original.size(); ++i ) {
            if ( std::isnan( original[i] ) )
                _test( std::isnan( result[i] ) );
            else
                _equal( result[i], original[i],
                        std::numeric_limits<double>::epsilon() );
        }
    }

    // --- vector<int> ---
    {
        const std::vector<int> original = { -1000, 0, 1, 42, INT_MAX, INT_MIN };
        auto fp = openWrite("brw_vec_int.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::vector<int> result;
        auto fp2 = openRead("brw_vec_int.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- vector<size_t> ---
    {
        std::vector<size_t> original( 100 );
        std::iota( original.begin(), original.end(), 0U );
        auto fp = openWrite("brw_vec_sizet.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::vector<size_t> result;
        auto fp2 = openRead("brw_vec_sizet.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Empty vector ---
    {
        const std::vector<double> original;
        auto fp = openWrite("brw_vec_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::vector<double> result = { 99.0 }; // must be cleared
        auto fp2 = openRead("brw_vec_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.empty() );
    }

    // --- Large vector (stress test) ---
    {
        std::vector<double> original( 100000 );
        std::iota( original.begin(), original.end(), 0.0 );
        auto fp = openWrite("brw_vec_large.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::vector<double> result;
        auto fp2 = openRead("brw_vec_large.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Read clears pre-existing contents ---
    {
        const std::vector<int> original = { 7, 8, 9 };
        auto fp = openWrite("brw_vec_clear.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::vector<int> result = { 1, 2, 3, 4, 5 }; // longer than original
        auto fp2 = openRead("brw_vec_clear.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        const std::vector<double> v = { 1.0 };
        _test( binaryFileWrite( fp, v ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        std::vector<double> v;
        _test( binaryFileRead( fp, v ) == false );
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — std::deque<T>
    // =========================================================================
    _info("binaryReadWrite_Test: testing deque<T> round-trips:");

    // --- deque<double> ---
    {
        const std::deque<double> original = { 0.1, 0.2, 0.3, -1.0, 1e10 };
        auto fp = openWrite("brw_deq_double.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::deque<double> result;
        auto fp2 = openRead("brw_deq_double.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.size() == original.size() );
        for ( size_t i = 0; i < original.size(); ++i )
            _equal( result[i], original[i],
                    std::numeric_limits<double>::epsilon() );
    }

    // --- deque<int32_t> ---
    {
        const std::deque<int32_t> original = { -5, 0, 5, 100, -100 };
        auto fp = openWrite("brw_deq_int.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::deque<int32_t> result;
        auto fp2 = openRead("brw_deq_int.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Empty deque ---
    {
        const std::deque<double> original;
        auto fp = openWrite("brw_deq_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::deque<double> result = { 1.0, 2.0 }; // must be cleared
        auto fp2 = openRead("brw_deq_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.empty() );
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        const std::deque<int> d = { 1, 2 };
        _test( binaryFileWrite( fp, d ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        std::deque<int> d;
        _test( binaryFileRead( fp, d ) == false );
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — std::deque<std::vector<T>>
    // =========================================================================
    _info("binaryReadWrite_Test: testing deque<vector<T>> round-trips:");

    // --- Normal case including inner empty vector ---
    {
        std::deque<std::vector<double>> original;
        original.push_back( { 1.0, 2.0, 3.0 } );
        original.push_back( { 4.0, 5.0 } );
        original.push_back( { 6.0 } );
        original.push_back( {} ); // inner empty vector

        auto fp = openWrite("brw_deqvec.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::deque<std::vector<double>> result;
        auto fp2 = openRead("brw_deqvec.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.size() == original.size() );
        for ( size_t i = 0; i < original.size(); ++i )
            _test( result[i] == original[i] );
    }

    // --- Empty outer deque ---
    {
        const std::deque<std::vector<int>> original;
        auto fp = openWrite("brw_deqvec_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::deque<std::vector<int>> result;
        result.push_back( { 1, 2, 3 } ); // must be cleared
        auto fp2 = openRead("brw_deqvec_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.empty() );
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        const std::deque<std::vector<int>> d;
        _test( binaryFileWrite( fp, d ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        std::deque<std::vector<int>> d;
        _test( binaryFileRead( fp, d ) == false );
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — std::map<M,T>
    // =========================================================================
    _info("binaryReadWrite_Test: testing map<M,T> round-trips:");

    // --- map<int, double> ---
    {
        const std::map<int,double> original = {
            { 0, 1.0 }, { 1, 2.5 }, { 2, -3.14 }, { 100, 0.0 }
        };
        auto fp = openWrite("brw_map.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::map<int,double> result;
        auto fp2 = openRead("brw_map.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.size() == original.size() );
        for ( const auto& [k, v] : original )
            _equal( result.at(k), v,
                    std::numeric_limits<double>::epsilon() );
    }

    // --- map<size_t, int32_t> ---
    {
        const std::map<size_t, int32_t> original = {
            { 0U, -1 }, { 10U, 0 }, { 999U, 42 }
        };
        auto fp = openWrite("brw_map_sizet.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::map<size_t, int32_t> result;
        auto fp2 = openRead("brw_map_sizet.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Empty map ---
    {
        const std::map<int,double> original;
        auto fp = openWrite("brw_map_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::map<int,double> result = { { 1, 1.0 } }; // must be cleared
        auto fp2 = openRead("brw_map_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.empty() );
    }

    // --- Read clears pre-existing contents ---
    {
        const std::map<int,int> original = { { 1, 10 }, { 2, 20 } };
        auto fp = openWrite("brw_map_clear.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::map<int,int> result = { { 99, 99 }, { 100, 100 }, { 101, 101 } };
        auto fp2 = openRead("brw_map_clear.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        const std::map<int,double> m = { { 1, 1.0 } };
        _test( binaryFileWrite( fp, m ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        std::map<int,double> m;
        _test( binaryFileRead( fp, m ) == false );
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — std::unordered_map<M,T>
    // =========================================================================
    _info("binaryReadWrite_Test: testing unordered_map<M,T> round-trips:");

    // --- unordered_map<int, double> ---
    {
        const std::unordered_map<int,double> original = {
            { 0, 1.0 }, { 1, 2.5 }, { 2, -3.14 }, { 100, 0.0 }
        };
        auto fp = openWrite("brw_umap.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::unordered_map<int,double> result;
        auto fp2 = openRead("brw_umap.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.size() == original.size() );
        for ( const auto& [k, v] : original )
            _equal( result.at(k), v,
                    std::numeric_limits<double>::epsilon() );
    }

    // --- Empty unordered_map ---
    {
        const std::unordered_map<int,int> original;
        auto fp = openWrite("brw_umap_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::unordered_map<int,int> result = { { 1, 1 } }; // must be cleared
        auto fp2 = openRead("brw_umap_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.empty() );
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        const std::unordered_map<int,double> m = { { 1, 1.0 } };
        _test( binaryFileWrite( fp, m ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        std::unordered_map<int,double> m;
        _test( binaryFileRead( fp, m ) == false );
    }


    // =========================================================================
    // binaryFileWrite / binaryFileRead — std::map<M, std::vector<T>>
    // =========================================================================
    _info("binaryReadWrite_Test: testing map<M,vector<T>> round-trips:");

    // --- map<int, vector<double>> including inner empty vector ---
    {
        const std::map<int, std::vector<double>> original = {
            { 0, { 1.0, 2.0, 3.0 } },
            { 1, { 4.0, 5.0 } },
            { 2, {} },
            { 3, { -1.0 } }
        };
        auto fp = openWrite("brw_mapvec.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::map<int, std::vector<double>> result;
        auto fp2 = openRead("brw_mapvec.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.size() == original.size() );
        for ( const auto& [k, v] : original )
            _test( result.at(k) == v );
    }

    // --- map<size_t, vector<int>> ---
    {
        std::map<size_t, std::vector<int>> original;
        for ( size_t i = 0; i < 5; ++i ) {
            std::vector<int> v( i + 1 );
            std::iota( v.begin(), v.end(), static_cast<int>(i) );
            original[i] = v;
        }
        auto fp = openWrite("brw_mapvec_sizet.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::map<size_t, std::vector<int>> result;
        auto fp2 = openRead("brw_mapvec_sizet.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result == original );
    }

    // --- Empty map ---
    {
        const std::map<int, std::vector<double>> original;
        auto fp = openWrite("brw_mapvec_empty.bin");
        _test( binaryFileWrite( fp, original ) );
        fp.close();

        std::map<int, std::vector<double>> result;
        result[99] = { 1.0, 2.0 }; // must be cleared
        auto fp2 = openRead("brw_mapvec_empty.bin");
        _test( binaryFileRead( fp2, result ) );
        _test( result.empty() );
    }

    // --- Write to closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        const std::map<int, std::vector<double>> m = { { 1, { 1.0 } } };
        _test( binaryFileWrite( fp, m ) == false );
    }

    // --- Read from closed stream returns false (expected diagnostic suppressed) ---
    {
        SilenceCerr sc;
        std::fstream fp;
        std::map<int, std::vector<double>> m;
        _test( binaryFileRead( fp, m ) == false );
    }


    // =========================================================================
    // Mixed sequential writes and reads in a single file
    // =========================================================================
    _info("binaryReadWrite_Test: testing mixed sequential writes in one file:");

    {
        const std::string           s_orig = "mixed_test";
        const std::vector<double>   v_orig = { 1.1, 2.2, 3.3 };
        const std::map<int,int>     m_orig = { { 1, 10 }, { 2, 20 } };
        const std::deque<int32_t>   d_orig = { 7, 8, 9 };

        auto fp = openWrite("brw_mixed.bin");
        BinaryFileSectionWrite( fp, "MIXED001" );
        _test( binaryFileWrite( fp, s_orig ) );
        _test( binaryFileWrite( fp, v_orig ) );
        _test( binaryFileWrite( fp, m_orig ) );
        _test( binaryFileWrite( fp, d_orig ) );
        fp.close();

        std::string         s_result;
        std::vector<double> v_result;
        std::map<int,int>   m_result;
        std::deque<int32_t> d_result;

        auto fp2 = openRead("brw_mixed.bin");
        BinaryFileSectionRead( fp2, "MIXED001" );
        _test( binaryFileRead( fp2, s_result ) );
        _test( binaryFileRead( fp2, v_result ) );
        _test( binaryFileRead( fp2, m_result ) );
        _test( binaryFileRead( fp2, d_result ) );

        _test( s_result == s_orig );
        _test( v_result == v_orig );
        _test( m_result == m_orig );
        _test( d_result == d_orig );
    }


    // =========================================================================
    // Multiple sections in one file
    // =========================================================================
    _info("binaryReadWrite_Test: testing multiple sections in one file:");

    {
        const std::vector<int>    v1 = { 1, 2, 3 };
        const std::vector<double> v2 = { 0.1, 0.2, 0.3 };

        auto fp = openWrite("brw_multisec.bin");
        BinaryFileSectionWrite( fp, "VSETHEDR" );
        _test( binaryFileWrite( fp, v1 ) );
        BinaryFileSectionWrite( fp, "VSETCORD" );
        _test( binaryFileWrite( fp, v2 ) );
        BinaryFileSectionWrite( fp, "VSETFOTR" );
        fp.close();

        std::vector<int>    r1;
        std::vector<double> r2;

        auto fp2 = openRead("brw_multisec.bin");
        BinaryFileSectionRead( fp2, "VSETHEDR" );
        _test( binaryFileRead( fp2, r1 ) );
        BinaryFileSectionRead( fp2, "VSETCORD" );
        _test( binaryFileRead( fp2, r2 ) );
        BinaryFileSectionRead( fp2, "VSETFOTR" );

        _test( r1 == v1 );
        _test( r2 == v2 );
    }

} // end run

} // namespace csmp

