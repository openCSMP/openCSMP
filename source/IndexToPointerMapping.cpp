//
//  IndexToPointerMapping.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 20/5/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "IndexToPointerMapping.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {


template<size_t dim>
void IndexToPointerMapping<dim>::Clear()
  {
     nodeConnector_.clear();
     elementConnector_.clear();
     faceConnector_.clear();
     interFaceConnector_.clear();
  }

} // end csmp

