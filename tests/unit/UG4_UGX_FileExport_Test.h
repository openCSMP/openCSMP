//
//  UG4_UGX_FileExport_Test.h
//  CSMP_API_library2024
//
//  Created by Stephan Matthai on 4/01/2024
//  Copyright (c) 2024 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_UG4_UGX_FILE_EXPORT_TEST_H
#define CSMP_UG4_UGX_FILE_EXPORT_TEST_H

#include "Test.h"
#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;


class  UG4_UGX_FileExport_Test : public Test {
  public:
    virtual void run();

  private:
    bool ImportModelAndRunChecks( const std::string& model_name,
                                  const std::string& variables_file ) const;
                                  
    bool RunChecks( const Model<3U>& ) const;
};

} // csmp

#endif // EXPERIMENTAL_EXAMPLE_H
