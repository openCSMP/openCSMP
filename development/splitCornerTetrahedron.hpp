//
//  splitCornerTetrahedron.hpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 4/11/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#ifndef splitCornerTetrahedron_hpp
#define splitCornerTetrahedron_hpp

#include <iostream>
#include <cmath>

namespace csmp {

template<size_t> class VSet;

/// converts corner tetrahedron and its neighbor into 3 tetrahedral cells, each with a face on the sides of the box
void splitCornerTetrahedron( VSet<3>&, size_t cnr, size_t only_neighbor );

}

#endif /* splitCornerTetrahedron_hpp */
