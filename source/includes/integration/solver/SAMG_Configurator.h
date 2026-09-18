// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SAMG_CONFIGURATOR_H
#define SAMG_CONFIGURATOR_H

#include "SAMG_Settings.h"
#include <string>

/// configuration of Fraunhofer's SAMG algebraic multigrid solver by string input.
namespace csmp {
  class SAMG_Configurator {
    public:
      virtual ~SAMG_Configurator() {};
      virtual SAMG_Settings* GetSettingsByName( const std::string& name );
    private:
  };

}// end namespace csmp

#endif
