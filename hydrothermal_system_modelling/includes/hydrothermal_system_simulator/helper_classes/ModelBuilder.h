// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef __ModelBuilder_h__
#define __ModelBuilder_h__

#include "ANSYS_Model2D.h"

using namespace std;
namespace csmp
{

  template<size_t dim>
  class ModelBuilder
  {

    public:

      ModelBuilder();
      void CreateVariablesFile( bool with_gold, bool with_lithium,
                                bool with_flux_visitor, bool with_tracer, bool with_degassing, bool with_split_boundary );
      void CreateLimiterLogFile();
      void DeleteVariablesFile();

      ~ModelBuilder();

    private:
      std::string var_file;

  };

}

#endif

