// Copyright © 2016 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_QUADRATURE_POINT_OUTPUT_TO_VTK_H
#define CSMP_QUADRATURE_POINT_OUTPUT_TO_VTK_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;

class QuadraturePointOutputToVTK  {
public:
	QuadraturePointOutputToVTK();
	~QuadraturePointOutputToVTK();
	void OutputQuadraturePointPropertiesAsDiscontinuousNodeVariablesToVTK(const Model<3U>&, const char*, const char*);
};

} // end csmp

#endif // end defined quadraturePointOutputToVTK.h
