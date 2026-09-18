// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GOCAD_OBJECT_H
#define GOCAD_OBJECT_H

#include "CSMP_definitions.h"
#include "GocadHeader.h"

namespace csmp {

class GocadObject {
  protected:
    GocadHeader  header;
  
  public:
    GocadObject();
    GocadObject( const GocadObject& go );
    ~GocadObject();
};

} // csp

#endif
