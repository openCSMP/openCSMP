#include "binaryReadWrite.h"
#include "Exception.h"

using namespace std;

namespace csmp {

BinaryFileSectionWrite::BinaryFileSectionWrite( std::fstream& fp, const char* header )
	: fp_(fp)
{
    char hdr[CSMP_BINARY_FILE_HDR_SIZE] = {};  // Zero-initialize the header array

    size_t hdrlen = std::strlen(header);
    assert(hdrlen <= sizeof(hdr));  // Assert that the header length fits

    // Copy the header into the hdr array, respecting the size
    // Use C++17's std::min for safer size comparison
    std::memcpy(hdr, header, std::min(hdrlen, sizeof(hdr)));

    // Write the entire hdr buffer to the file
    fp_.write(hdr, sizeof(hdr));

    if (!fp) {
        throw std::runtime_error("BinaryFileSectionWrite::BinaryFileSectionWrite: Error occurred while writing data tag to file");
    }
}


BinaryFileSectionRead::BinaryFileSectionRead( std::fstream& fp, const char* header )
	: fp_(fp)
{
    char readhdr[CSMP_BINARY_FILE_HDR_SIZE] = {};  // Zero-initialize the buffer
    const size_t hdrlen = std::strlen(header);
    
    // Ensure the provided header length doesn't exceed the buffer size
    assert(hdrlen <= CSMP_BINARY_FILE_HDR_SIZE);

    // Zero-initialize and copy the header using C++17 features
    std::memset(hdr_, 0, sizeof(hdr_));
    std::memcpy(hdr_, header, std::min(hdrlen, CSMP_BINARY_FILE_HDR_SIZE ) );

    // Read from the file into the buffer
    fp_.read(readhdr, CSMP_BINARY_FILE_HDR_SIZE);

    // Compare the headers, and if they don't match, throw an exception
    if (std::memcmp(hdr_, readhdr, CSMP_BINARY_FILE_HDR_SIZE) != 0) {
        if ( errno != 0 ) std::cout << "\n\t" << std::strerror(errno) << std::endl;
        throw std::runtime_error("Binary file entry appears to be corrupt");
    }
}



/**
     Reads the size record, returning how many elements the current record contains/
     
     If the size cannot be parsed or a numbe larger than can be contained in an unsigned 32-bit integer, a range error is thrown.
*/
size_t readContainerSize( std::fstream& fp )
 {
    assert( fp.is_open() );
    size_t elements(0U);
    
    if ( !fp.read(reinterpret_cast<char*>(&elements), sizeof(size_t) ) )
      {
        std::cerr << "\nreadContainerSize: ERROR: could not read record length." << std::endl;
        return elements;
      }
    if ( elements == numeric_limits<size_t>::max() ) {
         std::cerr <<"\n\t"<<"apparent size of container that shall be read: "<< elements << std::endl;
         throw csmp::Exception( ERROR, "readContainerSize:",
                               "number not initialised or too large to fit into 'size_t aka STL container::size_type.");
      }
    return elements;
 }




bool binaryFileWrite( fstream& fp, const char* str )
{
  if ( str == nullptr )
    std::cerr <<"\nbinaryFileWrite: WARNING: string is nullptr."<< std::endl;
  
	if (!fp.is_open()) {
		cout << "\nbinaryFileWrite (const char*): ERROR: invalid file pointer." << endl;
		return false;
	}
	// writing the size of the object
	const size_t  characters = strlen(str);
	fp.write( reinterpret_cast<const char*>(&characters), sizeof(size_t) );

	// writing the character string
	fp.write( reinterpret_cast<const char*>(str), characters );

	return true;
}



bool binaryFileRead( fstream& fp, char str[] )
{
	if (!fp.is_open()) {
		cerr << "\nbinaryFileRead(char[]): ERROR: invalid file pointer." << endl;
		return false;
	}
	// read size of the record and assert this 
	size_t  characters(0);
	char    buf[INFO_STRING];

	if ( !fp.read( reinterpret_cast<char*>(&characters), sizeof(size_t) ) )
    {
      cerr << "\nbinaryFileRead(char[]): ERROR: could not read string length." << endl;
      return false;
    }

	if (characters > INFO_STRING)
		throw csmp::Exception( ERROR, "binaryFileRead", "too many characters in input string" );

	// reading the character string
	fp.read( reinterpret_cast<char*>(buf), characters );	
	if (characters != fp.gcount())
    {
      cerr << "\nbinaryFileRead(char[]) ERROR: incorrect number of characters were read: ";
      cerr << "\nIndicated number: " << characters << ", actual number read: " << strlen(buf) << endl;
      return false;
    }
	// null terminate string and copy to 'str' argument
	buf[characters] = '\0';
	strcpy(str, buf);

  if (!fp) {
      throw std::runtime_error("binaryFileRead: Error occurred while reading from file");
    }

	return true;
}



bool binaryFileWrite( fstream& fp, const std::string& str )
{
  if ( str.empty() )
    std::cerr <<"\nbinaryFileWrite (string): WARNING: string is empty."<< std::endl;
  
	if (!fp.is_open()) {
		cout << "\nbinaryFileWrite (const char*): ERROR: invalid file pointer." << endl;
		return false;
	}
	// writing the size of the object
	const size_t  characters = str.size();
	fp.write( reinterpret_cast<const char*>(&characters), sizeof(size_t) );

	// writing the character string (since C++1.7 string is guaranteed to be contiguous in memory)
	fp.write( str.data(), characters );

  if (!fp) {
      throw std::runtime_error("binaryFileWrite: Error occurred while writing string to file");
  }

	return true;
}



bool binaryFileRead( fstream& fp, string& str )
{
	if (!fp.is_open()) {
		cerr << "\nbinaryFileRead(char[]): ERROR: invalid file pointer." << endl;
		return false;
	}
	// read size of the string
	size_t  characters = 0;
	if ( !fp.read( reinterpret_cast<char*>(&characters), sizeof(size_t) ) )
    {
      cerr << "\nbinaryFileRead(string): ERROR: could not read string length." << endl;
      return false;
    }
  str.clear();
  str.resize( characters, '\0' );

	// reading the character string
	fp.read( reinterpret_cast<char*>(str.data()), characters );
	if (characters != fp.gcount())
    {
      cerr << "\nbinaryFileRead(string) ERROR: incorrect number of characters were read: ";
      cerr << "\nIndicated number: " << characters << ", actual number read: " << str.size() << endl;
      return false;
    }

	 return true;
}



} // end namespace csmp
