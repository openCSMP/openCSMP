// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SKUA_EXAMPLE_H
#define SKUA_EXAMPLE_H

#include "Example.h"

#include <string>

namespace csmp {
    template<uint32_t> class Model;
}

namespace csmp {

class  SKUA_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
private:
  static void ListModels();
  static void ImportModelAndRunChecks( const std::string& model_name );
  static void RunChecks( const Model<3U>& model );
  static void RunSKUA_model();
  static void RunSKUA_box_shaped_with_boundary();
  static void RunSKUA_split_boundary_layer();
  static void RunSKUA_cross_bedded_xsmall();
  static void RunSKUA_cross_bedded_small();
  static void RunSKUA_cross_bedded();
  static void RunSnowFlake();
};

} // csmp

#endif // SKUA_EXAMPLE_H

