// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STOKES_DISCREPANCY_MEASURE_EXAMPLE_H
#define STOKES_DISCREPANCY_MEASURE_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  StokesDiscrepancyMeasure_Example final : public Example{
  public:
    ~StokesDiscrepancyMeasure_Example() override = default;
    
    void Run() override;
    void Specifications() override;
};

} // csmp


#endif // STOKES_DISCREPANCY_MEASURE_EXAMPLE_H
