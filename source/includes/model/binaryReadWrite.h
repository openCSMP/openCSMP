#ifndef CSMP_READ_WRITE_H
#define CSMP_READ_WRITE_H

#include "CSMP_definitions.h"
#include "PropertyDatabase.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

/**
  @file binaryReadWrite.h

  @brief Binary file I/O utilities for CSMP STL containers and variables.

  All write functions serialise a leading `size_t` element count followed by
  the raw data.  All read functions expect the same layout and validate the
  count before populating the output container.

  Error handling policy
  ---------------------
  - A closed or invalid file stream is always a programming error: functions
    return `false` and print a diagnostic to `std::cerr`.
  - A corrupt record (unexpected element count, truncated read) throws
    `std::runtime_error` so that the caller's stack is unwound cleanly.
  - An empty container is **not** an error: zero is written/read as a valid
    element count and the function returns `true`.
*/

// ============================================================================
// Constants
// ============================================================================

/// Size of a binary section-header tag in bytes, including null terminator.
constexpr size_t CSMP_BINARY_FILE_HDR_SIZE = 8 + 1;


// ============================================================================
// Section-header helpers
// ============================================================================

/**
  @brief Writes a fixed-width section tag to a binary file.

  The tag is exactly CSMP_BINARY_FILE_HDR_SIZE bytes (8 characters + null).
  Throws `std::runtime_error` if the write fails.
*/
class BinaryFileSectionWrite {
  public:
    BinaryFileSectionWrite( std::fstream& fp, const char* header );
  private:
    std::fstream& fp_;
};


/**
  @brief Reads and validates a fixed-width section tag from a binary file.

  Throws `std::runtime_error` if the tag does not match the expected header,
  indicating a corrupt or misaligned file.
*/
class BinaryFileSectionRead {
  public:
    BinaryFileSectionRead( std::fstream& fp, const char* header );
  private:
    std::fstream& fp_;
};


// ============================================================================
// Low-level size reader
// ============================================================================

/**
  @brief Read a leading `size_t` element count from the file stream.

  Throws `std::runtime_error` if the read fails or the value equals
  `std::numeric_limits<size_t>::max()`, which is used as a sentinel for
  an uninitialised record.
*/
size_t readContainerSize( std::fstream& fp );


// ============================================================================
// Forward declarations (template bodies follow below)
// ============================================================================

template<class T>
bool binaryFileWrite( std::fstream& fp, const std::vector<T>& c );

template<class T>
bool binaryFileRead(  std::fstream& fp, std::vector<T>& c );

template<typename T>
bool binaryFileWrite( std::fstream& fp, const std::deque<T>& c );

template<typename T>
bool binaryFileRead(  std::fstream& fp, std::deque<T>& c );

template<class T>
bool binaryFileWrite( std::fstream& fp, const std::deque<std::vector<T>>& c );

template<class T>
bool binaryFileRead(  std::fstream& fp, std::deque<std::vector<T>>& c );

template<class M, class T>
bool binaryFileWrite( std::fstream& fp, const std::map<M, T>& c );

template<class M, class T>
bool binaryFileWrite( std::fstream& fp, const std::unordered_map<M, T>& c );

template<class M, class T>
bool binaryFileRead(  std::fstream& fp, std::map<M, T>& c );

template<class M, class T>
bool binaryFileRead(  std::fstream& fp, std::unordered_map<M, T>& c );

template<class M, class T>
bool binaryFileWrite( std::fstream& fp, const std::map<M, std::vector<T>>& c );

template<class M, class T>
bool binaryFileRead(  std::fstream& fp, std::map<M, std::vector<T>>& c );

bool binaryFileWrite( std::fstream& fp, const std::string& s );
bool binaryFileRead(  std::fstream& fp, std::string& s );
bool binaryFileWrite( std::fstream& fp, const char* str );
bool binaryFileRead(  std::fstream& fp, char str[] );


// ============================================================================
// Domain variable I/O dispatch helpers
// ============================================================================

namespace femDataOutputDispatch {

template<class Var>
inline void initVariable( csmp::Index, Var& ) {}

template<>
inline void initVariable( csmp::Index key, ArrayVariable& var )
{ var.Resize( key.dataDepth ); }

template<>
inline void initVariable( csmp::Index key, FlaggedArrayVariable& var )
{ var.Resize( key.dataDepth ); }

} // namespace femDataOutputDispatch


// ============================================================================
// Domain variable binary I/O
// ============================================================================

/// Write all variables of a single type from a domain to file.
template<class V, class D, uint32_t dim> requires CsmpVariable<dim, V>
bool variablesOut( std::fstream& fp, const D& domain,
                   const PropertyDatabase<dim>& pref, VARIABLE_TYPE vtype )
{
    std::set<std::string> propList;
    pref.ListVariables( domain.Placement(), vtype, propList );
    const size_t vcount = propList.size();
    fp.write( reinterpret_cast<const char*>(&vcount), sizeof(size_t) );

    for ( const auto& name : propList ) {
        V var;
        const Index key = pref.StorageKey( name.c_str() );
        femDataOutputDispatch::initVariable( key, var );
        domain.Read( key, var );
        binaryFileWrite( fp, name.c_str() );
        if ( !var.Out( fp ) ) {
            std::cerr << "\nvariablesOut: ERROR: failed to write variable '"
                      << name << "'.\n";
            return false;
        }
    }
    return true;
}


/// Read all variables of a single type from file into a domain.
template<class V, class D, uint32_t dim> requires CsmpVariable<dim, V>
bool variablesIn( std::fstream& fp, D& domain,
                  const PropertyDatabase<dim>& pref, VARIABLE_TYPE )
{
    // Fix #10: initialise to 0, not max(); read is checked by readContainerSize.
    const size_t vcount = readContainerSize( fp );

    for ( size_t i = 0; i < vcount; ++i ) {
        V    var;
        char propName[NAME_STRING];
        binaryFileRead( fp, propName );

        if ( pref.IsDefined( propName ) ) {
            const Index key = pref.StorageKey( propName );
            femDataOutputDispatch::initVariable( key, var );
            if ( !var.In( fp ) ) {
                std::cerr << "\nvariablesIn: ERROR: failed to read variable '"
                          << propName << "' (record " << i << " of " << vcount << ").\n";
                return false;
            }
            domain.Store( key, var );
        } else {
            // Variable not in database — read and discard to keep stream aligned.
            if ( !var.In( fp ) ) {
                std::cerr << "\nvariablesIn: ERROR: failed to skip unknown variable '"
                          << propName << "' (record " << i << " of " << vcount << ").\n";
                return false;
            }
        }
    }
    return true;
}


/// Write all variable types from a domain to file.
template<class D, uint32_t dim>
bool domainVariablesOut( std::fstream& fp, const D& domain,
                         const PropertyDatabase<dim>& pref )
{
    return variablesOut<ScalarVariable>      ( fp, domain, pref, SCALAR      )
        && variablesOut<VectorVariable<dim>> ( fp, domain, pref, VECTOR      )
        && variablesOut<TensorVariable<dim>> ( fp, domain, pref, TENSOR      )
        && variablesOut<ArrayVariable>       ( fp, domain, pref, ARRAY       )
        && variablesOut<FlaggedArrayVariable>( fp, domain, pref, FLAGGEDARRAY);
}


/// Read all variable types from file into a domain.
template<class D, uint32_t dim>
bool domainVariablesIn( std::fstream& fp, D& domain,
                        const PropertyDatabase<dim>& pref )
{
    return variablesIn<ScalarVariable>      ( fp, domain, pref, SCALAR      )
        && variablesIn<VectorVariable<dim>> ( fp, domain, pref, VECTOR      )
        && variablesIn<TensorVariable<dim>> ( fp, domain, pref, TENSOR      )
        && variablesIn<ArrayVariable>       ( fp, domain, pref, ARRAY       )
        && variablesIn<FlaggedArrayVariable>( fp, domain, pref, FLAGGEDARRAY);
}


/**
  @brief Selective variable reader — reads only variables whose names appear
         in @p selection.

  Variables present in the file but absent from @p selection are read and
  discarded so that the stream remains correctly positioned.

  @author SKM
  @date   6/9/2021
*/
template<class V, class D, uint32_t dim> requires CsmpVariable<dim, V>
bool selectedVariablesIn( std::fstream& fp, D& domain,
                          const PropertyDatabase<dim>& pref,
                          VARIABLE_TYPE,
                          const std::set<std::string>& selection )
{
    // Fix #10: use readContainerSize rather than initialising to max().
    const size_t vcount = readContainerSize( fp );

    for ( size_t i = 0; i < vcount; ++i ) {
        V    var;
        char propName[NAME_STRING];
        binaryFileRead( fp, propName );

        const bool wanted = selection.contains( propName )
                         && pref.IsDefined( propName );

        if ( wanted ) {
            const Index key = pref.StorageKey( propName );
            femDataOutputDispatch::initVariable( key, var );
            if ( !var.In( fp ) ) {
                std::cerr << "\nselectedVariablesIn: ERROR: failed to read variable '"
                          << propName << "' (record " << i << " of " << vcount << ").\n";
                return false;
            }
            domain.Store( key, var );
        } else {
            // Read and discard to keep stream aligned.
            if ( !var.In( fp ) ) {
                std::cerr << "\nselectedVariablesIn: ERROR: failed to skip variable '"
                          << propName << "' (record " << i << " of " << vcount << ").\n";
                return false;
            }
        }
    }
    return true;
}


/// Selective reader for all variable types.
template<class D, uint32_t dim>
bool selectedDomainVariablesIn( std::fstream& fp, D& domain,
                                const PropertyDatabase<dim>& pref,
                                const std::set<std::string>& selection )
{
    return selectedVariablesIn<ScalarVariable>      ( fp, domain, pref, SCALAR,       selection )
        && selectedVariablesIn<VectorVariable<dim>> ( fp, domain, pref, VECTOR,       selection )
        && selectedVariablesIn<TensorVariable<dim>> ( fp, domain, pref, TENSOR,       selection )
        && selectedVariablesIn<ArrayVariable>       ( fp, domain, pref, ARRAY,        selection )
        && selectedVariablesIn<FlaggedArrayVariable>( fp, domain, pref, FLAGGEDARRAY, selection );
}


// ============================================================================
// Template definitions — vector
// ============================================================================

/**
  @brief Write a `std::vector<T>` to a binary stream.

  Writes a leading `size_t` element count followed by the raw element data
  in a single contiguous block.  An empty vector writes a zero count and
  returns `true` — zero is a valid record size.

  @return `true` on success, `false` if the stream is not open.
*/
template<class T>
bool binaryFileWrite( std::fstream& fp, const std::vector<T>& c )
{
    // Fix #1/#2: is_open checked first; empty container is not an error.
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(vector<T>): ERROR: file stream is not open.\n";
        return false;
    }
    const size_t n = c.size();
    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    if ( n > 0 )
        fp.write( reinterpret_cast<const char*>( c.data() ),
                  static_cast<std::streamsize>( sizeof(T) * n ) );
    return true;
}


/**
  @brief Read a `std::vector<T>` from a binary stream.

  Clears the output container, reads the leading element count, resizes the
  vector, then reads all elements in a single block read.

  @return `true` on success, `false` if the stream is not open.
  @throws `std::runtime_error` if the block read is truncated.
*/
template<class T>
bool binaryFileRead( std::fstream& fp, std::vector<T>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(vector<T>): ERROR: file stream is not open.\n";
        return false;
    }
    // Fix #9: always clear before populating.
    c.clear();
    const size_t n = readContainerSize( fp );
    if ( n == 0 ) return true;

    c.resize( n );
    fp.read( reinterpret_cast<char*>( c.data() ),
             static_cast<std::streamsize>( sizeof(T) * n ) );

    if ( static_cast<size_t>( fp.gcount() ) != sizeof(T) * n ) {
        throw std::runtime_error(
            "binaryFileRead(vector<T>): truncated read — file may be corrupt. "
            "Expected " + std::to_string( sizeof(T) * n ) +
            " bytes, got "  + std::to_string( fp.gcount() ) + "." );
    }
    return true;
}


// ============================================================================
// Template definitions — deque
// ============================================================================

/**
  @brief Write a `std::deque<T>` to a binary stream.

  Deque storage is not contiguous, so elements are written one at a time.
  An empty deque writes a zero count and returns `true`.

  @return `true` on success, `false` if the stream is not open.
*/
template<typename T>
bool binaryFileWrite( std::fstream& fp, const std::deque<T>& c )
{
    // Fix #1/#2: is_open checked first; empty container is not an error.
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(deque<T>): ERROR: file stream is not open.\n";
        return false;
    }
    const size_t n = c.size();
    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    for ( const auto& elem : c )
        fp.write( reinterpret_cast<const char*>(&elem), sizeof(T) );
    return true;
}


/**
  @brief Read a `std::deque<T>` from a binary stream.

  @return `true` on success, `false` if the stream is not open.
  @throws `std::runtime_error` if fewer elements than expected are read.
*/
template<typename T>
bool binaryFileRead( std::fstream& fp, std::deque<T>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(deque<T>): ERROR: file stream is not open.\n";
        return false;
    }
    c.clear();
    const size_t n = readContainerSize( fp );
    if ( n == 0 ) return true;

    size_t count = 0;
    T val;
    for ( size_t i = 0; i < n; ++i ) {
        // Fix #15: comment corrected — this is reading, not writing.
        fp.read( reinterpret_cast<char*>(&val), sizeof(T) );
        if ( fp.gcount() == static_cast<std::streamsize>( sizeof(T) ) ) {
            c.push_back( val );
            ++count;
        }
    }
    if ( count != n ) {
        throw std::runtime_error(
            "binaryFileRead(deque<T>): truncated read — file may be corrupt. "
            "Expected " + std::to_string(n) +
            " elements, successfully read " + std::to_string(count) + "." );
    }
    return true;
}


// ============================================================================
// Template definitions — deque<vector<T>>
// ============================================================================

/**
  @brief Write a `std::deque<std::vector<T>>` to a binary stream.

  Writes the outer element count, then delegates each inner vector to
  `binaryFileWrite(vector<T>)`.

  @return `true` on success, `false` if the stream is not open.
*/
template<class T>
bool binaryFileWrite( std::fstream& fp, const std::deque<std::vector<T>>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(deque<vector<T>>): ERROR: file stream is not open.\n";
        return false;
    }
    const size_t n = c.size();
    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    for ( const auto& vec : c )
        binaryFileWrite( fp, vec );
    return true;
}


/**
  @brief Read a `std::deque<std::vector<T>>` from a binary stream.

  @return `true` on success, `false` if the stream is not open.
  @throws `std::runtime_error` if fewer inner vectors than expected are read.
*/
template<class T>
bool binaryFileRead( std::fstream& fp, std::deque<std::vector<T>>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(deque<vector<T>>): ERROR: file stream is not open.\n";
        return false;
    }
    c.clear();
    const size_t n = readContainerSize( fp );
    if ( n == 0 ) return true;

    size_t count = 0;
    for ( size_t i = 0; i < n; ++i ) {
        std::vector<T> val;
        // Fix #11: binaryFileRead(vector) now returns true for empty vectors,
        // so this counter correctly reflects the number of records read.
        if ( binaryFileRead( fp, val ) ) {
            c.push_back( std::move(val) );
            ++count;
        } else {
            throw std::runtime_error(
                "binaryFileRead(deque<vector<T>>): failed to read inner vector "
                + std::to_string(i) + " of " + std::to_string(n) + "." );
        }
    }
    return true;
}


// ============================================================================
// Template definitions — map<M,T>
// ============================================================================

/**
  @brief Write a `std::map<M,T>` to a binary stream.

  Writes the element count followed by each key-value pair as raw bytes.

  @note Both `M` and `T` must be trivially copyable (no heap-allocated members).

  @return `true` on success, `false` if the stream is not open.
*/
template<class M, class T>
bool binaryFileWrite( std::fstream& fp, const std::map<M, T>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(map<M,T>): ERROR: file stream is not open.\n";
        return false;
    }
    const size_t n = c.size();
    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    for ( const auto& [key, val] : c ) {
        fp.write( reinterpret_cast<const char*>(&key), sizeof(M) );
        fp.write( reinterpret_cast<const char*>(&val), sizeof(T) );
    }
    return true;
}


/**
  @brief Read a `std::map<M,T>` from a binary stream.

  @return `true` on success, `false` if the stream is not open.
  @throws `std::runtime_error` if fewer records than expected are read.
*/
template<class M, class T>
bool binaryFileRead( std::fstream& fp, std::map<M, T>& c )
{
    // Fix #8: is_open checked before any mutation of the output container.
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(map<M,T>): ERROR: file stream is not open.\n";
        return false;
    }
    c.clear();
    const size_t n = readContainerSize( fp );
    if ( n == 0 ) return true;

    size_t count = 0;
    for ( size_t i = 0; i < n; ++i ) {
        M key;
        T val;
        fp.read( reinterpret_cast<char*>(&key), sizeof(M) );
        fp.read( reinterpret_cast<char*>(&val), sizeof(T) );
        if ( fp ) {
            c[key] = val;
            ++count;
        }
    }
    if ( count != n ) {
        throw std::runtime_error(
            "binaryFileRead(map<M,T>): truncated read — file may be corrupt. "
            "Expected " + std::to_string(n) +
            " key-value pairs, successfully read " + std::to_string(count) + "." );
    }
    return true;
}


// ============================================================================
// Template definitions — unordered_map<M,T>
// ============================================================================

/**
  @brief Write a `std::unordered_map<M,T>` to a binary stream.

  @note Iteration order of `unordered_map` is not deterministic; the file
  will be valid but the order of records may differ between runs.

  @return `true` on success, `false` if the stream is not open.
*/
template<class M, class T>
bool binaryFileWrite( std::fstream& fp, const std::unordered_map<M, T>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(unordered_map<M,T>): ERROR: file stream is not open.\n";
        return false;
    }
    const size_t n = c.size();
    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    for ( const auto& [key, val] : c ) {
        fp.write( reinterpret_cast<const char*>(&key), sizeof(M) );
        fp.write( reinterpret_cast<const char*>(&val), sizeof(T) );
    }
    return true;
}


/**
  @brief Read a `std::unordered_map<M,T>` from a binary stream.

  @return `true` on success, `false` if the stream is not open.
  @throws `std::runtime_error` if fewer records than expected are read.
*/
template<class M, class T>
bool binaryFileRead( std::fstream& fp, std::unordered_map<M, T>& c )
{
    // Fix #8: is_open checked before any mutation of the output container.
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(unordered_map<M,T>): ERROR: file stream is not open.\n";
        return false;
    }
    c.clear();
    const size_t n = readContainerSize( fp );
    if ( n == 0 ) return true;

    c.reserve( n );
    size_t count = 0;
    for ( size_t i = 0; i < n; ++i ) {
        M key;
        T val;
        fp.read( reinterpret_cast<char*>(&key), sizeof(M) );
        fp.read( reinterpret_cast<char*>(&val), sizeof(T) );
        if ( fp ) {
            c[key] = val;
            ++count;
        }
    }
    if ( count != n ) {
        throw std::runtime_error(
            "binaryFileRead(unordered_map<M,T>): truncated read — file may be corrupt. "
            "Expected " + std::to_string(n) +
            " key-value pairs, successfully read " + std::to_string(count) + "." );
    }
    return true;
}


// ============================================================================
// Template definitions — map<M, vector<T>>
// ============================================================================

/**
  @brief Write a `std::map<M, std::vector<T>>` to a binary stream.

  Writes the map size, then for each entry writes the key as raw bytes
  followed by the vector via `binaryFileWrite(vector<T>)`.

  @return `true` on success, `false` if the stream is not open.
*/
template<class M, class T>
bool binaryFileWrite( std::fstream& fp, const std::map<M, std::vector<T>>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(map<M,vector<T>>): ERROR: file stream is not open.\n";
        return false;
    }
    const size_t n = c.size();
    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    for ( const auto& [key, vec] : c ) {
        fp.write( reinterpret_cast<const char*>(&key), sizeof(M) );
        binaryFileWrite( fp, vec );
    }
    return true;
}


/**
  @brief Read a `std::map<M, std::vector<T>>` from a binary stream.

  @return `true` on success, `false` if the stream is not open.
  @throws `std::runtime_error` if fewer records than expected are read.
*/
template<class M, class T>
bool binaryFileRead( std::fstream& fp, std::map<M, std::vector<T>>& c )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(map<M,vector<T>>): ERROR: file stream is not open.\n";
        return false;
    }
    c.clear();
    const size_t n = readContainerSize( fp );
    if ( n == 0 ) return true;

    size_t count = 0;
    for ( size_t i = 0; i < n; ++i ) {
        M key;
        fp.read( reinterpret_cast<char*>(&key), sizeof(M) );
        std::vector<T> val;
        if ( binaryFileRead( fp, val ) ) {
            c[key] = std::move( val );
            ++count;
        } else {
            throw std::runtime_error(
                "binaryFileRead(map<M,vector<T>>): failed to read vector for "
                "map entry " + std::to_string(i) + " of " + std::to_string(n) + "." );
        }
    }
    if ( count != n ) {
        throw std::runtime_error(
            "binaryFileRead(map<M,vector<T>>): truncated read — file may be corrupt. "
            "Expected " + std::to_string(n) +
            " entries, successfully read " + std::to_string(count) + "." );
    }
    return true;
}

} // namespace csmp

#endif // CSMP_READ_WRITE_H

