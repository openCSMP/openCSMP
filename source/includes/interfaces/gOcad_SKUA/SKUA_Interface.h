//
//  SKUA_Interface.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/20/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_SKUA_INTERFACE_H
#define CSMP_SKUA_INTERFACE_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;

/**

@section Some conventions regarding neighbor connectivity

In CSMP, only the neighbors of equidimensional elements are recorded: 
the volumetric neighbors of volumetric elements, the surface neighbors of surface elements and the line element neighbors of bar elements.
Now, on the outer perimeter of the model, elements that make up the boundary surfaces would normally have the elements of adjacent 
boundary surfaces as neighbors. This connectivity is ignored. Instead, missing neighbors are indicated in the VSet::VData by the BOX_BOUNDARY
flags that correspond to the model edges.

*/
class SKUA_Interface {
  public:
  
  // IMPORT INTERFACES
  
    /// region-by-region property assignment from column-based texfile (Kuncho Kurtev)
    template<uint32_t dim>
    bool ImportElementPropertyValuesFromSKUA( Model<dim>&, const std::string& data_file=std::string() );


  // OUTPUT INTERFACES FOR POINT DATA (DOIMOI FAULT MODELLING 2011-2013)

    /// output of element barycentre points transformed into SKUA (lefthanded) coordinate system for property interpolation and transfer to CSMP
    void OutputElementNumbersAndBaryCentresRegionByRegion( const Model<3U>& );

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
  
  // MISC OPERATIONS

    /// For each rock type the elements will be converted into regions with the corresponding name
    void ConvertRockTypesIntoRegions( Model<3U>&, const std::string& rocktype_info_file );

    void Erase_NO_DATA_ElementsFromModel( Model<3U>&, const std::string& target_region, std::set<uint32_t>& no_data_elmt_numbers );

  private:
  
  // HELPER METHODS
  
    /// Enlists 'element number's of elements with NO_DATA values (-9999, -99999) in the target region
    bool Detect_NO_DATA_ElementsInDatasetFromSKUA( const std::string& input_txt_file, const std::string& target_region,
                                                   std::set<uint32_t>& no_data_elmt_numbers );

};


} // end csmp


#endif /* defined(CSMP_SKUA_INTERFACE_H) */
