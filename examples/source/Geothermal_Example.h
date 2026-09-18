// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOTHERMAL_EXAMPLE_H
#define GEOTHERMAL_EXAMPLE_H

#include "Example.h"

#include "CSMP_definitions.h"

namespace csmp {

const size_t DIM(2U);
template<uint32_t> class Model;

class Geothermal_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
    
  private:
      void ComputeMassConductivity ( Model<DIM>& );
      void OutputToVTU( Model<DIM>&, std::string model_name, const std::list<std::string>& props, size_t timestep ) const;
      bool Compare( Model<DIM>&, const std::string& file ) const;
      void ComputeMassGravityTerm (Model<DIM>&);
};


} // csmp

#endif
