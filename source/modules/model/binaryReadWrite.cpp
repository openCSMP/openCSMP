#include "binaryReadWrite.h"
#include "Exception.h"

using namespace std;

namespace csmp {

bool skm_C_fwrite(std::fstream& fp, const char* str)
{
	if (!fp.is_open()) {
		cout << "\nskm_C_fwrite (const char*): ERROR: invalid file pointer." << endl;
		return false;
	}
	// writing the size of the object
	size_t  characters = strlen(str);
	fp.write((char*)&characters, sizeof(size_t));

	// writing the character string
	fp.write((char*)str, characters);

	return true;
}



bool skm_C_fread(std::fstream& fp, char str[])
{
	if (!fp.is_open()) {
		cout << "\nskm_C_fread(char[]): ERROR: invalid file pointer." << endl;
		return false;
	}
	// read size of the record and assert this 
	size_t  characters(0);
	char    buf[INFO_STRING];

	if (!fp.read((char*)&characters, sizeof(size_t)))
	{
		cout << "\nskm_C_fread(char[]): ERROR: could not read string length." << endl;
		return false;
	}

	if (characters > INFO_STRING) {
		throw csmp::Exception(ERROR, "skm_C_fread", "Binary file appears to be corrupt");
	}

	// reading the character string
	fp.read((char*)buf, characters);	
	if (characters != fp.gcount())
	{
		cout << "\nskm_C_fread(char[]) ERROR: incorrect number of characters were read: ";
		cout << "\nIndicated number: " << characters << ", actual number read: " << strlen(buf) << endl;
		return false;
	}
	// null terminate string and copy to 'str' argument
	buf[characters] = '\0';
	strcpy(str, buf);

	return true;
}



BinaryFileSectionRead::BinaryFileSectionRead(std::fstream& fp, const char* header)
	: fp_(fp)
{
	char readhdr[CSMP_BINARY_FILE_HDR_SIZE];
	size_t hdrlen = strlen(header);
	assert(hdrlen <= CSMP_BINARY_FILE_HDR_SIZE);
	memset(hdr_, 0, sizeof(hdr_));
	memcpy(hdr_, header, std::min(hdrlen, (size_t)CSMP_BINARY_FILE_HDR_SIZE));
	fp.read(readhdr, sizeof(char) * CSMP_BINARY_FILE_HDR_SIZE);
	fp.read((char*)&offset_, sizeof(offset_));
	if (memcmp(hdr_, readhdr, CSMP_BINARY_FILE_HDR_SIZE)) {
		throw csmp::Exception(FATAL_ERROR, "BinaryFileSectionRead", hdr_, "Binary file appears to be corrupt");
	}
	sectoffset_ = fp.tellg();
}

BinaryFileSectionRead::~BinaryFileSectionRead()
{
	ulong64 off = fp_.tellg();
	if (off != offset_) {
		throw csmp::Exception(FATAL_ERROR, "BinaryFileSectionRead", hdr_, "Binary file appears to be corrupt");
	}
}

BinaryFileSectionWrite::BinaryFileSectionWrite(std::fstream& fp, const char* header)
	: fp_(fp)
{
	char hdr[8];
	size_t hdrlen = strlen(header);
	assert(hdrlen <= sizeof(hdr));
	memset(hdr, 0, sizeof(hdr));
	memcpy(hdr, header, std::min(hdrlen, sizeof(hdr)));
	fp.write(hdr, sizeof(hdr) / sizeof(char));
	offset_ = fp.tellg();
	fp.write((char*)&offset_, sizeof(offset_));
}

BinaryFileSectionWrite::~BinaryFileSectionWrite()
{
	ulong64 off = fp_.tellg();
	fp_.seekg(offset_, fp_.beg);
	fp_.write((char*)&off, sizeof(off));
	fp_.seekg(off, fp_.beg);
}


} // end namespace csmp
