//
//  TestAids.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 11/11/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#ifndef TestAids_h
#define TestAids_h

#include "CSMP_definitions.h"

namespace csmp {

template<size_t dim, template<size_t> class CELL>
void backupNeighborConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                                 typename std::vector<CELL<dim>*>::const_iterator last,
                                std::vector<std::vector<CELL<dim>*> >& nbor_pointers );




} // end 

#endif /* TestAids_hpp */
