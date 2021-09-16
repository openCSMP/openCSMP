#include "binaryReadWrite.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// JCK 2018
BinaryFileSectionRead::BinaryFileSectionRead( std::fstream& fp, const char* header )
	: fp_(fp)
{
	char  readhdr[CSMP_BINARY_FILE_HDR_SIZE];
	const size_t hdrlen = strlen(header);
	assert(hdrlen <= CSMP_BINARY_FILE_HDR_SIZE);
	memset(hdr_, 0, sizeof(hdr_));
	memcpy(hdr_, header, std::min(hdrlen, (size_t)CSMP_BINARY_FILE_HDR_SIZE));
	fp.read(readhdr, sizeof(char) * CSMP_BINARY_FILE_HDR_SIZE);
	if ( memcmp(hdr_, readhdr, CSMP_BINARY_FILE_HDR_SIZE) ) {
		throw csmp::Exception(FATAL_ERROR, "BinaryFileSectionRead", hdr_, "Binary file entry appears to be corrupt");
	}
}



// JCK 2018
BinaryFileSectionRead::~BinaryFileSectionRead()
{
}



// JCK 2018
BinaryFileSectionWrite::BinaryFileSectionWrite( std::fstream& fp, const char* header )
	: fp_(fp)
{
	char hdr[CSMP_BINARY_FILE_HDR_SIZE];
	size_t hdrlen = strlen(header);
	assert(hdrlen <= sizeof(hdr));
	memset(hdr, 0, sizeof(hdr));
	memcpy(hdr, header, std::min(hdrlen, sizeof(hdr)));
	fp.write(hdr, sizeof(hdr) / sizeof(char));
}



// JCK 2018
BinaryFileSectionWrite::~BinaryFileSectionWrite()
{
}


/**
     Reads the size record, returning how many elements the current record contains/
     
     If the size cannot be parsed or a numbe larger than can be contained in an unsigned 32-bit integer, a range error is thrown.
*/
size_t checkContainerSize( std::fstream& fp )
 {
    assert( fp.is_open() );
    size_t elements(0U);
    
    if ( !fp.read(reinterpret_cast<char*>(&elements), sizeof(size_t) ) )
      {
        std::cerr << "\nbool binaryFileRead: ERROR: could not read record length." << std::endl;
        return elements;
      }
    if ( elements >= UINT32_MAX ) {
         std::cerr <<"\n\tapparent size of container that shall be read: "<< elements << std::endl;
         throw csmp::Exception( ERROR, "checkContainerSize:",
                               "number is too large to fit into 'size_t aka STL container::size_type.");
      }
    return elements;
 }




bool binaryFileWrite( fstream& fp, const char* str )
{
  if ( str == nullptr )
    std::cerr <<"nbinaryFileWrite: WARNING: string is nullptr."<< std::endl;
  
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

	return true;
}



bool binaryFileWrite( fstream& fp, const std::string& str )
{
  if ( str.empty() )
    std::cerr <<"nbinaryFileWrite (string): WARNING: string is empty."<< std::endl;
  
	if (!fp.is_open()) {
		cout << "\nbinaryFileWrite (const char*): ERROR: invalid file pointer." << endl;
		return false;
	}
	// writing the size of the object
	const size_t  characters = str.size();
	fp.write( reinterpret_cast<const char*>(&characters), sizeof(size_t) );

	// writing the character string (since C++1.7 string is guaranteed to be contiguous in memory)
	fp.write( &str[0], characters );

	return true;
}



bool binaryFileRead( fstream& fp, string& str )
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
      cerr << "\nbinaryFileRead(string): ERROR: could not read string length." << endl;
      return false;
    }

	if (characters > INFO_STRING)
		throw csmp::Exception( ERROR, "binaryFileRead(string)", "too many characters in input string" );

	// reading the character string
	fp.read( reinterpret_cast<char*>(buf), characters );
	if (characters != fp.gcount())
    {
      cerr << "\nbinaryFileRead(string) ERROR: incorrect number of characters were read: ";
      cerr << "\nIndicated number: " << characters << ", actual number read: " << strlen(buf) << endl;
      return false;
    }
	// null terminate string and copy to 'str' argument
	buf[characters] = '\0';
	str = buf;

	return true;
}



} // end namespace csmp
