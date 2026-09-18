//
//  BoundaryInterface_Test.h
//
//  Created by Stephan Matthai on 23/03/2016.
//

#ifndef BOUNDARY_INTERFACE_TEST_H
#define BOUNDARY_INTERFACE_TEST_H

#include "Test.h"
#include "Model.h"

namespace csmp {

/**
    @brief tests refactored functionality related to boundaries.
    
    Input models:
      - fault_boundary_test   - normal fault offsetting central reservoir layer creating complex intersection relationships
        uses the mesh stored in FaultBoundaryTest
*/
class BoundaryInterface_Test : public Test {
  public:
    /// testing ""fault_boundary_test" creation of a segmented boundary from region 'fault'
    virtual void run();
  
    /// can Box-shaped model creation be accomplished and how fast
    void TestBoxShapedModel();
  
  private:
    /// checks whether the contact surface of 2 contacting regions is recovered correctly
    bool TestRegionContactDetection( const Model<3U>& );
    
    /// visualises topotypes by converting them to nodel scalar variables output to VTU file
    void TopoTypeToVTU( Model<3U>&, const char* region );
    
    /// VTU output if topotypes for lower-dim regions
    void TestTopoTypeIdentifiers( Model<3U>& );
    
    /// checks whether TOPOTYPE  node information is correct after the new boundaries have been inserted
    void TestBoundaryAndTopoTypeIdentifiers( Model<3U>& );

  private:
    const static bool verbose_ = true; ///< turn off to supress output of files etc.
};


/// finding substrings in boundary name; returns '\0' if unsuccessful; use model.Boundary(findBoundary(...)) to recover boundary after search
std::string  findBoundary( const Model<3U>& model, const std::set<std::string>& intersected_regions );

/// assuming a dim-1 region, label and count material juxtaposition relationships
size_t countAndLabelRegions( Model<3U>&, std::vector<std::string>& region_names );

/// discerning patches by values for the region in terms of the diagnostic element variable
size_t labelRegionPatches( Model<3U>& model, const char* dim_1_region, const char* diagnostic_elmt_variable,
                           const std::vector<std::string>& region_names );

/// establishes neighbors of lower dimensional element
void higherDimensionalNeighbors( const Element<3U>& e, const csmp::Index& mtrl_key,
                                 std::pair<size_t,size_t>& nbors, std::pair<long,long>& materials );


} // end csmp

#endif /* SKM_BOUNDARY_FUNCTIONALITY_TEST_H */
