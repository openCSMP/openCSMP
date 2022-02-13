//
//  quadraturePointOutputToVTK.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 9/06/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#ifndef QUADRATURE_POINT_OUTPUT_TO_VTK_H
#define QUADRATURE_POINT_OUTPUT_TO_VTK_H

#include "ErrorHandler.h"
#include "VTK_Interface.h"
#include "CSMP_highLevelUtilities.h"

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
