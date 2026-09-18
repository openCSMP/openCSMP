// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_REGION_MONITOR_H
#define CSMP_REGION_MONITOR_H

#include "Model.h"

namespace csmp {

template<uint32_t dim> class Model; 

/** 
    @brief Monitors value ranges and integrated property values in all regions of a model.
    
    Among the non-unique (potentially overlapping) regions those with names identical to
    BOX_BOUNDARY flags and the "All Element" region are ignored.
    
    @attention if the model contains lower-dimensional regions, a thickness attribute 
    called "thickness" must be specified for these in order to get their volume right.
    For the regions with the same dimension as the model, the thickness attribute
    should have a value of 1.
    
    @author Stephan Matthai
    @date 1999
*/
template<uint32_t dim>
class RegionMonitor {
  public:
    RegionMonitor();
  
    /// construct with 2 properties to monitor
    RegionMonitor( const Model<dim>&, 
                   const std::string& first_integral_property,
                   const std::string& first_range_property,
                   bool integrate_pore_volume_only=true );

    /// construct with lists of multiple range and/or integral properties
    RegionMonitor( const Model<dim>&, 
                   const std::list<std::string>& to_integrate_over_groups,
                   const std::list<std::string>& to_find_ranges_in_groups,
                   bool integrate_pore_volume_only=true );
                  
    ~RegionMonitor();
  
    /// tell RegionMonitor which properties it should track
    void  DefineProperties( const Model<dim>&, 
                            const std::list<std::string>& to_integrate_over_groups,
                            const std::list<std::string>& to_find_ranges_in_groups,
                            bool integrate_pore_volume_only=true );
                            
    /// switches on normalisation by volume for all property value integrations
    void  DivideIntegralPropertiesByRegionVolumes( bool divide );

    /// returns true when there are already monitoring data stored in respective containers
    bool  DataExists() const { return ( !integrals_.empty() || !ranges_.empty() ); }
  
    /// calculates integrals over current property values and stores these for the current time-step and property name
    void  ScalarPropertyIntegrals( const Model<dim>&, double current_time=0. );
  
    /// monitors current ranges of the the scalar values
    void  ScalarPropertyRanges( const Model<dim>&, double current_time=0. );

    /** adds in values in a separate container, to be output to a separate file
        with the same prefix name as the other monitoring file, except with the "_ext" added to its name.
    */
    void  InsertExternallyCalculatedProperty( std::string property, std::string regionname );
    void  InsertExternallyCalculatedProperty( std::string property_regionname );

    void  InsertPreCalculatedPropertyValue( double time,
                                            std::string property,
                                            std::string regionname,
                                            double value);
  
    void  InsertPreCalculatedPropertyValue( double time,
                                            std::string property_regionname,
                                            double value);

    void  Reset(); ///< zap all recorded values and property names
    void  EraseData(); ///< zap all recorded values
    
    void  Out( const char* text_file ) const;

  private:
    /// internal helper for calculations on models with involve surface elements
    bool  HasVolumeElements( const Model<dim>& ) const;
  
    // geometric group properties
    //         groupname,          volume,   surface area
    std::map<std::string,std::pair<double,double> >  group_specs_;
    // properties which are monitored
    std::list<std::string>  integral_properties_;
    std::list<std::string>  range_properties_;
    //    model_time,          property,             groupname, property value
    std::map<double,std::map<std::string,std::map<std::string,double> > >  integrals_;

    //    model_time,          property,              groupname,          property range
    std::map<double,std::map<std::string,std::map<std::string,std::pair<double,double> > > >  ranges_;

    // for externally calculated monitored properties
    std::set<std::string>  ext_calc_properties_column_headers_;

    //    model_time,          property_groupname,   property value
    std::map<double,std::map<std::string,double> >  ext_properties_;

    bool divide_by_volume_,
         pore_volume_integral_,
         include_thickness_attribute_;
};

} // csmp

#endif
