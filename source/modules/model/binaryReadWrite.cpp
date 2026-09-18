// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "binaryReadWrite.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// ============================================================================
// BinaryFileSectionWrite
// ============================================================================

BinaryFileSectionWrite::BinaryFileSectionWrite( std::fstream& fp, const char* header )
    : fp_( fp )
{
    // Fix #3: guard against headers that would overflow the fixed buffer,
    // leaving no room for the null terminator.
    const size_t hdrlen = std::strlen( header );
    if ( hdrlen >= CSMP_BINARY_FILE_HDR_SIZE ) {
        throw std::invalid_argument(
            "BinaryFileSectionWrite: header '" + std::string(header) +
            "' is " + std::to_string(hdrlen) +
            " characters, which exceeds the maximum of " +
            std::to_string(CSMP_BINARY_FILE_HDR_SIZE - 1) + "." );
    }

    // Zero-fill so the full fixed-width tag is always written cleanly.
    char hdr[CSMP_BINARY_FILE_HDR_SIZE] = {};
    std::memcpy( hdr, header, hdrlen ); // null terminator already in place via zero-fill

    fp_.write( hdr, sizeof(hdr) );

    if ( !fp_ ) {
        throw std::runtime_error(
            "BinaryFileSectionWrite: failed to write section tag '" +
            std::string(header) + "' to file." );
    }
}


// ============================================================================
// BinaryFileSectionRead
// ============================================================================

BinaryFileSectionRead::BinaryFileSectionRead( std::fstream& fp, const char* header )
    : fp_( fp )
{
    // Fix #3/#12: hdr_ is now a local variable — it is only needed for the
    // comparison and does not need to be a class member.
    const size_t hdrlen = std::strlen( header );
    if ( hdrlen >= CSMP_BINARY_FILE_HDR_SIZE ) {
        throw std::invalid_argument(
            "BinaryFileSectionRead: expected header '" + std::string(header) +
            "' is " + std::to_string(hdrlen) +
            " characters, which exceeds the maximum of " +
            std::to_string(CSMP_BINARY_FILE_HDR_SIZE - 1) + "." );
    }

    char expected[CSMP_BINARY_FILE_HDR_SIZE] = {};
    std::memcpy( expected, header, hdrlen );

    char actual[CSMP_BINARY_FILE_HDR_SIZE] = {};
    fp_.read( actual, sizeof(actual) );

    if ( !fp_ ) {
        throw std::runtime_error(
            "BinaryFileSectionRead: failed to read section tag from file "
            "(expected '" + std::string(header) + "')." );
    }

    if ( std::memcmp( expected, actual, CSMP_BINARY_FILE_HDR_SIZE ) != 0 ) {
        // Include both tags in the message to make misalignment easy to diagnose.
        const std::string actual_str( actual,
            std::min( std::strlen(actual), CSMP_BINARY_FILE_HDR_SIZE - 1 ) );
        if ( errno != 0 )
            std::cerr << "\n\tSystem error: " << std::strerror(errno) << "\n";
        throw std::runtime_error(
            "BinaryFileSectionRead: section tag mismatch — "
            "expected '" + std::string(header) +
            "', found '"  + actual_str +
            "'. File may be corrupt or the read position is misaligned." );
    }
}


// ============================================================================
// readContainerSize
// ============================================================================

size_t readContainerSize( std::fstream& fp )
{
    assert( fp.is_open() );

    size_t n = 0;
    fp.read( reinterpret_cast<char*>(&n), sizeof(size_t) );

    // Fix #7: distinguish a failed read from a legitimate zero.
    if ( !fp ) {
        throw std::runtime_error(
            "readContainerSize: failed to read element count from file — "
            "the file may be truncated or the read position is misaligned." );
    }

    // A value of max() is used as a sentinel for an uninitialised record,
    // which indicates the writing code failed to populate the size field.
    if ( n == std::numeric_limits<size_t>::max() ) {
        throw std::runtime_error(
            "readContainerSize: element count is SIZE_MAX, which is the "
            "uninitialised sentinel value — the file record was not written "
            "correctly." );
    }

    return n;
}


// ============================================================================
// binaryFileWrite / binaryFileRead — const char*
// ============================================================================

bool binaryFileWrite( fstream& fp, const char* str )
{
    // Fix #5: nullptr check must prevent further execution.
    if ( str == nullptr ) {
        std::cerr << "\nbinaryFileWrite(const char*): ERROR: string pointer is null.\n";
        return false;
    }
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(const char*): ERROR: file stream is not open.\n";
        return false;
    }

    const size_t n = std::strlen( str );

    // Fix #6/#13: validate against streamsize max before writing anything.
    if ( n > static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) ) {
        throw std::out_of_range(
            "binaryFileWrite(const char*): string length " + std::to_string(n) +
            " exceeds the maximum writable size of " +
            std::to_string( std::numeric_limits<std::streamsize>::max() ) + "." );
    }

    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    if ( n > 0 )
        fp.write( str, static_cast<std::streamsize>(n) );

    if ( !fp ) {
        throw std::runtime_error(
            "binaryFileWrite(const char*): write failed after " +
            std::to_string(n) + " characters." );
    }
    return true;
}


bool binaryFileRead( fstream& fp, char str[] )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(char[]): ERROR: file stream is not open.\n";
        return false;
    }

    size_t n = 0;
    if ( !fp.read( reinterpret_cast<char*>(&n), sizeof(size_t) ) ) {
        std::cerr << "\nbinaryFileRead(char[]): ERROR: could not read string length.\n";
        return false;
    }

    // Fix #4: use >= to leave room for the null terminator.
    if ( n >= INFO_STRING ) {
        throw std::runtime_error(
            "binaryFileRead(char[]): string length " + std::to_string(n) +
            " would overflow the fixed buffer of size " +
            std::to_string(INFO_STRING) + " (including null terminator)." );
    }

    // Fix #13: validate against streamsize max.
    if ( n > static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) ) {
        throw std::out_of_range(
            "binaryFileRead(char[]): string length " + std::to_string(n) +
            " exceeds the maximum readable size." );
    }

    char buf[INFO_STRING] = {};
    fp.read( buf, static_cast<std::streamsize>(n) );

    if ( static_cast<size_t>( fp.gcount() ) != n ) {
        throw std::runtime_error(
            "binaryFileRead(char[]): truncated read — expected " +
            std::to_string(n) + " characters, got " +
            std::to_string( fp.gcount() ) + "." );
    }

    buf[n] = '\0';
    std::strcpy( str, buf );
    return true;
}


// ============================================================================
// binaryFileWrite / binaryFileRead — std::string
// ============================================================================

bool binaryFileWrite( fstream& fp, const std::string& str )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileWrite(string): ERROR: file stream is not open.\n";
        return false;
    }

    const size_t n = str.size();

    // Fix #6/#13: validate before writing anything.
    if ( n > static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) ) {
        throw std::out_of_range(
            "binaryFileWrite(string): string length " + std::to_string(n) +
            " exceeds the maximum writable size of " +
            std::to_string( std::numeric_limits<std::streamsize>::max() ) + "." );
    }

    fp.write( reinterpret_cast<const char*>(&n), sizeof(size_t) );
    // std::string is guaranteed contiguous since C++11.
    if ( n > 0 )
        fp.write( str.data(), static_cast<std::streamsize>(n) );

    if ( !fp ) {
        throw std::runtime_error(
            "binaryFileWrite(string): write failed after " +
            std::to_string(n) + " characters." );
    }
    return true;
}


bool binaryFileRead( fstream& fp, string& str )
{
    if ( !fp.is_open() ) {
        std::cerr << "\nbinaryFileRead(string): ERROR: file stream is not open.\n";
        return false;
    }

    size_t n = 0;
    if ( !fp.read( reinterpret_cast<char*>(&n), sizeof(size_t) ) ) {
        std::cerr << "\nbinaryFileRead(string): ERROR: could not read string length.\n";
        return false;
    }

    // Fix #13: validate against streamsize max.
    if ( n > static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) ) {
        throw std::out_of_range(
            "binaryFileRead(string): string length " + std::to_string(n) +
            " exceeds the maximum readable size of " +
            std::to_string( std::numeric_limits<std::streamsize>::max() ) + "." );
    }

    str.clear();
    if ( n == 0 ) return true;

    str.resize( n );
    // std::string is guaranteed contiguous since C++11.
    fp.read( reinterpret_cast<char*>( str.data() ),
             static_cast<std::streamsize>(n) );

    if ( static_cast<size_t>( fp.gcount() ) != n ) {
        throw std::runtime_error(
            "binaryFileRead(string): truncated read — expected " +
            std::to_string(n) + " characters, got " +
            std::to_string( fp.gcount() ) + "." );
    }

    return true;
}

} // namespace csmp

