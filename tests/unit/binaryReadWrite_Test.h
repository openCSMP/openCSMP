//
//  binaryReadWrite_Test.h
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 7/04/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_BINARY_READ_WRITE_TEST_H
#define CSMP_BINARY_READ_WRITE_TEST_H

#include "Test.h"
#include "CSMP_definitions.h"
#include "binaryReadWrite.h"

namespace csmp {

/**
  @brief Comprehensive test for all binaryReadWrite functions.

  Covers:
    - `BinaryFileSectionWrite` / `BinaryFileSectionRead`
    - `readContainerSize`
    - `binaryFileWrite` / `binaryFileRead` for:
        - `std::vector<T>`
        - `std::deque<T>`
        - `std::deque<std::vector<T>>`
        - `std::map<M,T>`
        - `std::unordered_map<M,T>`
        - `std::map<M, std::vector<T>>`
        - `std::string`
        - `const char*` / `char[]`
    - Empty container round-trips
    - Error handling (closed stream, corrupt tag, truncated record)
*/
class binaryReadWrite_Test : public Test {
  public:
    virtual void run();
    const static bool verbose_ = false;

  private:
    /// Open a fresh binary file for writing; throws on failure.
    static std::fstream openWrite( const std::string& name );
    /// Open an existing binary file for reading; throws on failure.
    static std::fstream openRead( const std::string& name );
};

} // namespace csmp

#endif // CSMP_BINARY_READ_WRITE_TEST_H

