//
//  smoothData.hpp
//  CSMP_FECFVM_Simulator
//
//  Created by Stephan Matthai on 4/2/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_SMOOTH_ELEMENT_DATA_H
#define CSMP_SMOOTH_ELEMENT_DATA_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;

enum SMOOTHING_TYPE { MEAN, HARMONIC_MEAN, LOG10_MEAN, HORIZONTAL_MEAN };

/**
    smoothElementData - reads CSMP binary file  and performs various kinds of smoothing on the target element variable.
    
    @parameter smoothing_passes tells functions how often the same operation is to be performed
    
    @parameter log10smoothing  takes the decadic log value of the properties before the smoothing is performed.
    
    @parameter horizontal_direction_only restricts operation to the horizontal
*/
template<uint32_t dim>
void smoothElementData( Model<dim>&, 
                        const std::string& region_to_be_smoothed, 
                        const std::string& variable_name,
                        int smoothing_passes = 1,           ///< usually sufficient
                        bool log10_smoothing = false,       ///< creates a more patchy pattern with emphasis on the high values
                        bool in_plane_smoothing = false,    ///< only in the horizontal plane in as much as algorithm can resolve 
                        bool output_model_to_binary = false ///< write the model to binary
                      ); 

/// smoothes porosity and permeability distributions in current model
template<uint32_t dim>
void smoothPorosityAndPermeabilityDistribution( Model<dim>& );


/// Edoardo Pezulli's neighbor extrapolation based method, that avoids Element and Node objects located at BOX boundary
template<uint32_t dim>
void spreadPropertiesOfInitialisedCellsAcross( typename std::vector<csmp::Element<dim>*>::iterator first,
                                               typename std::vector<csmp::Element<dim>*>::iterator last );


} // end csmp

#endif /* CSMP_SMOOTH_ELEMENT_DATA_H */
