//
//  SKUA_Interface.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/20/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_SKUA_INTERFACE_H
#define CSMP_SKUA_INTERFACE_H

#include <iostream>
#include <set>

namespace csmp {

template<size_t> class Model;

class SKUA_Interface {
  public:
    /// exports element variable from selected regions to SKUA barycenter-point cloud (ASCII table) format 
    void VariableToPointCloud( const Model<3U>&,
                               const std::set<std::string>& regions_of_interest,
                               const char* filename, const char* element_var ) const;

    /// exports multiple element variables from selected regions to SKUA barycenter-point cloud (ASCII table) format 
    void VariablesToPointCloud( const Model<3U>&,
                                const std::set<std::string>& regions_of_interest,
                                const char* filename, const std::set<std::string>& element_var ) const;

    /// creates cloud of points offset along the barycenter surface normal of the element of interest using the 'thickness' attribute
    void SurfaceArrayVariableToPointCloud( const Model<3U>&,
                                           const std::set<std::string>& regions_of_interest,
                                           const char* filename, const char* element_var ) const;
};


} // end csmp


#endif /* defined(CSMP_SKUA_INTERFACE_H) */
