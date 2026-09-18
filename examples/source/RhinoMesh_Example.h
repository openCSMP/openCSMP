// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RHINOMESH_EXAMPLE_H
#define RHINOMESH_EXAMPLE_H

#include "Example.h"
#include "CSMP_definitions.h"

namespace csmp {

  template<uint32_t> class Model;

class  RhinoMesh_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
private:
  void SideBoundaryConditions( Model<3U>& );
  void ConcentrationRectangle( Model<3U>&, double );
  void CopyInputFiles(std::string& model_name, std::string& variable_file);
};

} // csmp

#endif // RHINOMESH_EXAMPLE_H
