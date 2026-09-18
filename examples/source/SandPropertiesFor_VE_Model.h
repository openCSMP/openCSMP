// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  SandPropertiesFor_VE_Model.hpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 26/8/2022.
//

#ifndef SAND_PROPERTIES_FOR_VE_MODEL_H
#define SAND_PROPERTIES_FOR_VE_MODEL_H

#include "CSMP_definitions.h"

namespace csmp {

class Index;
template<uint32_t> class Model;
template<uint32_t> class Element;

/**
      Class is used to initialise uniform properties of lower-dimensional sandbodies
      so that an application of the Brooks-Corey model will yied the correct input
      parameters for an element with a certain dip etc.
      
      While the model needs 'total velocity' as an input, all other properties are calculated
      locally.
*/
class SandPropertiesFor_VE_Model {
  public:
    /// original properties that the model returns the distributed parameters to once the specific computation has been done
    SandPropertiesFor_VE_Model( const Model<2U>& model, const std::string& dimM1_region_name,
                                double bcp, double swr, double snr );
                                
    SandPropertiesFor_VE_Model( const Model<3U>& model, const std::string& dimM1_region_name,
                                double bcp, double swr, double snr ) { /* implemented yet */ }
    
    /// specific version that only works for line elements in a  two dimensional model
    void Compute2PhaseFlowPropertiesForSandLayer( Model<2U>& );
    
    bool HasLineElementRepresentation() const { return target_region_consists_of_line_elmts_; }
    
  private:
    /// finds h/H at upstream or downstream facet of inclided FV from the updip or down-dip water saturation and the flow direction
    double CO2LevelAtFacet( const Element<2U>* const eptr ) const;
  
  private:
    const std::string target_region_; ///< the lower dimensional region that this object will manage
    // diagnostic parameters
    const csmp::Index thi_key_;       ///< "thickness"
    const csmp::Index fvV_key_;       ///< "finite volume"
    const csmp::Index vt_key_;        ///< "total velocity"
    const csmp::Index dip_key_;       ///< "dip vector"
    // multiphase flow parameters that get computed
    const csmp::Index sw_key_;        ///< "saturation aqueous phase"
    const csmp::Index bcp_key_;       ///< "brooks corey parameter"
    const csmp::Index swr_key_;       ///< "residual saturation aqueous phase"
    const csmp::Index snr_key_;       ///< "residual saturation carbonic phase"
    //const csmp::Index krw_key_;       ///< "relative permeability of water"
    //const csmp::Index krn_key_;       ///< "relative permeability of CO2"
    const csmp::Index mobw_key_;      ///< "water mobility"
    const csmp::Index mobn_key_;      ///< "CO2 mobility"
  ///
    // backup parameter values
    double bcp_;                      ///< Brooks-Corey parameter, linear model if zero
    double swr_;                      ///< irreducible water saturation
    double snr_;                      ///< residual non-wetting phase saturation
    // fluid properties
    double muw_, mun_;                ///< fluid viscosities
  ///
    bool   target_region_consists_of_line_elmts_ = true;
};

} // end csmp

#endif /* SAND_PROPERTIES_FOR_VE_MODEL_H */
