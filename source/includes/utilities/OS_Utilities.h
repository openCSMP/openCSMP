#ifndef OPERATING_SYSTEM_UTILITIES_H
#define OPERATING_SYSTEM_UTILITIES_H

#include <string>

namespace csmp {

/// create directory
std::string currentDirectorySymbol();
std::string directorySymbol();
void createDirectoryIfDoesntExist(  const std::string& directory_name );

/// to get file modification time.
std::string getFileModificationTime(const char *filePath);

std::string CompareFileModifiedTimeStamps(std::string file0_name,std::string file1_name);

} // end namespace csmp

#endif
