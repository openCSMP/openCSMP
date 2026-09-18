// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MODELANSYS_EXAMPLE_H
#define MODELANSYS_EXAMPLE_H

#include "Example.h"

namespace csmp {

  template<uint32_t> class Model;

class  ModelANSYS_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();

private:
  void CopyInputFiles(std::string& model_name, std::string& variable_file);
};

} // csmp

#endif // MODELANSYS_EXAMPLE_H
