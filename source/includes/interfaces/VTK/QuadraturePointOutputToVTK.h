//
//  quadraturePointOutputToVTK.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 9/06/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

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
