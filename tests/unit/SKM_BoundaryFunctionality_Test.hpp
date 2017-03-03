//
//  SKM_BoundaryFunctionality_Test.hpp
//
//  Created by Stephan Matthai on 23/03/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#ifndef SKM_BOUNDARY_FUNCTIONALITY_TEST_H
#define SKM_BOUNDARY_FUNCTIONALITY_TEST_H

#include "Test.h"
#include "Model.h"

namespace csmp {

class BoundaryFunctionality_Test : public Test {
  public:
    virtual void run();
  
    /// can Box-shaped model creation be accomplished and how fast
    void TestBoxShapedModel();
  
};


/// finding substrings in boundary name; returns '\0' if unsuccessful; use model.Boundary(findBoundary(...)) to recover boundary after search
std::string  findBoundary( const Model<3U>& model, const std::set<std::string>& intersected_regions );

/// assuming a dim-1 region, label and count material juxtaposition relationships
size_t countAndLabelRegions( Model<3U>&, std::vector<std::string>& region_names );

/// discerning patches by values for the region in terms of the diagnostic element variable
size_t labelRegionPatches( Model<3U>& model, const char* dim_1_region, const char* diagnostic_elmt_variable,
                           const std::vector<std::string>& region_names );

/// converts region into Boundarie(s) of faces, decomposed into patches; region is moved from "Model" to non-unique, connectivity is updated 
// add to BoundaryInterface
size_t createInternalBoundaryFromLowerDimensionalRegion( Model<3U>& model, const char* dim_1_region, const char* diagnostic_elmt_variable,
                                                         const std::vector<std::string>& region_names );


/// establishes neighbors of lower dimensional element
void higherDimensionalNeighbors( const Element<3U>& e, const csmp::Index& mtrl_key,
                                 std::pair<size_t,size_t>& nbors, std::pair<long,long>& materials );


} // end csmp

#endif /* SKM_BOUNDARY_FUNCTIONALITY_TEST_H */
