#include "RegionMonitor.h"
#include "Node.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
RegionMonitor<dim>::RegionMonitor()
    : divide_by_volume_(false),
      pore_volume_integral_(false),
      include_thickness_attribute_(false)
{
    cout <<"\nRegionMonitor<"<< dim;
    cout <<">: building empty monitor facility."<< endl;
}


template<uint32_t dim>
RegionMonitor<dim>::RegionMonitor( const Model<dim>& sg, 
                                   const string& first_integral_property,
                                   const string& first_range_property,
                                   bool integrate_pore_volume_only)

    : divide_by_volume_(false),
      pore_volume_integral_(integrate_pore_volume_only),
      include_thickness_attribute_(sg.Database().IsDefined("thickness"))
{
    list<string>  integral_properties; integral_properties.push_back( first_integral_property );
    list<string>  range_properties;    range_properties.push_back( first_range_property );

    DefineProperties( sg, integral_properties, range_properties, integrate_pore_volume_only );
    
} // end constructor




template<uint32_t dim>
RegionMonitor<dim>::RegionMonitor( const Model<dim>& sg, 
                                   const list<string>& to_integrate_over_groups,
                                   const list<string>& to_find_ranges_in_groups,
                                   bool integrate_pore_volume_only )
    : integral_properties_(to_integrate_over_groups),
      range_properties_(to_find_ranges_in_groups),
      divide_by_volume_(false),
      pore_volume_integral_(integrate_pore_volume_only),
      include_thickness_attribute_(sg.Database().IsDefined("thickness"))
{

    // P. Lang fix for surface calculation of 1D, 2D and 3D DFN models
    const bool hasVolumeElements( HasVolumeElements( sg ) );

    // checking whether target properties are suitable for integration
    for ( typename list<string>::const_iterator
          it=integral_properties_.begin(); it!=integral_properties_.end(); it++ ) {
        if ( !sg.Database().IsDefined( (*it).c_str() ) )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor(constructor)", (*it).c_str(),
                                   "integral property is not defined in the database (file)");

        if ( sg.Database().Type( (*it).c_str() ) != SCALAR )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor(constructor)", (*it).c_str(),
                                   "integral property must be a scalar property");

        if ( sg.Database().Placement( (*it).c_str() ) == FACE )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor(constructor)", (*it).c_str(),
                                   "Boundary properties placed on the FACE cannot be integrated over the Region volume");

        if ( sg.Database().Placement( (*it).c_str() ) == INTER_FACE )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor(constructor)", (*it).c_str(),
                                   "SplitBoundary properties placed on the INTER_FACE cannot be integrated over the Region volume");
    }

    // No need for checking the externally integrated properties as these do not
    // necesarily need to exist as a variable in the region

    // checking range properties
    for ( typename list<string>::const_iterator
          it=range_properties_.begin(); it!=range_properties_.end(); it++ )
        if ( !sg.Database().IsDefined( (*it).c_str() ) )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor(constructor)", (*it).c_str(),
                                   "range property is not defined in the database (file)");

    // getting geometric data from the existing regions
    double  volume;

    // unique regions
    for ( typename map<string,Region<dim> >::const_iterator
          git=sg.UniqueRegionsBegin(); git!=sg.UniqueRegionsEnd(); git++ ) {
        // total volume or pore volume of region
        if ( integrate_pore_volume_only ) volume = (*git).second.Volume(true);
        else                              volume = (*git).second.Volume(false);
        // surface area of region                = the length of the 1D model
        const double  surface_area = ( !hasVolumeElements ) ? (*git).second.Volume(false) : (*git).second.SurfaceArea();

        // recording the geometric properties
        string regionname=(*git).first;
        group_specs_[ regionname ] = std::make_pair(volume,surface_area);
    }

    // non-unique regions
    for ( typename map<string,Region<dim> >::const_iterator
          git=sg.RegionsBegin(); git!=sg.RegionsEnd(); git++ )
      if ( parseBoundary((*git).first) == NOT and (*git).first !="Model" ) {
          // total volume or pore volume of region
          if ( integrate_pore_volume_only ) volume = (*git).second.Volume(true);
          else                              volume = (*git).second.Volume(false);
          // surface area of region
          const double  surface_area = ( !hasVolumeElements ) ? (*git).second.Volume(false) : (*git).second.SurfaceArea();

          // recording the geometric properties
          group_specs_[ (*git).first ] = make_pair(volume,surface_area);
       }

} // end constructor




template<uint32_t dim>
void RegionMonitor<dim>::DefineProperties( const Model<dim>& sg, 
                                           const list<string>& to_integrate_over_groups,
                                           const list<string>& to_find_ranges_in_groups,
                                           bool integrate_pore_volume_only )
{
    integral_properties_  = to_integrate_over_groups;
    range_properties_     = to_find_ranges_in_groups;
    pore_volume_integral_ = integrate_pore_volume_only;

    // P. Lang fix for surface calculation of 1D, 2D and 3D DFN models
    bool hasVolumeElements( HasVolumeElements( sg ) );
    
    // perform some checks on the validity of property names
    typename list<string>::const_iterator  it;
    
    // checking whether target properties are suitable for integration
    for ( it=integral_properties_.begin(); it!=integral_properties_.end(); it++ ) {
        if ( !sg.Database().IsDefined( (*it).c_str() ) )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor::DefineProperties", (*it).c_str(),
                                   "integral property is not defined in the database (file)");

        if ( sg.Database().Type( (*it).c_str() ) != SCALAR )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor::DefineProperties)", (*it).c_str(),
                                   "integral property must be a scalar property");

        if ( sg.Database().Placement( (*it).c_str() ) == INTER_FACE )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor::DefineProperties", (*it).c_str(),
                                   "properties placed on the FACE cannot be integrated over the Region volume");
    }
    
    // checking range properties
    for ( it=range_properties_.begin(); it!=range_properties_.end(); it++ )
        if ( !sg.Database().IsDefined( (*it).c_str() ) )
            throw csmp::Exception( FATAL_ERROR, "RegionMonitor::DefineProperties", (*it).c_str(),
                                   "range property is not defined in the database (file)");

    // getting geometric data from the existing groups
    double  volume;

    for ( typename map<string,Region<dim> >::const_iterator
          git=sg.UniqueRegionsBegin(); git!=sg.UniqueRegionsEnd(); git++ ) {
        // total volume or pore volume of group
        if ( integrate_pore_volume_only ) volume = (*git).second.Volume(true);
        else                              volume = (*git).second.Volume(false);
        // surface area of group
        const double  surface_area = ( !hasVolumeElements ) ? (*git).second.Volume(false) : (*git).second.SurfaceArea();

        // recording the geometric properties
        group_specs_[ (*git).first ] = make_pair(volume,surface_area);
    }
    for ( typename map<string,Region<dim> >::const_iterator
          git=sg.RegionsBegin(); git!=sg.RegionsEnd(); git++ )
      if ( parseBoundary((*git).first) == NOT and (*git).first !="Model" ) {
          // total volume or pore volume of group
          if ( integrate_pore_volume_only ) volume = (*git).second.Volume(true);
          else                              volume = (*git).second.Volume(false);
          // surface area of group
          const double  surface_area = ( !hasVolumeElements ) ? (*git).second.Volume(false) : (*git).second.SurfaceArea();

          // recording the geometric properties
          group_specs_[ (*git).first ] = make_pair(volume,surface_area);
       }

} // end DefineProperties






template<uint32_t dim>
RegionMonitor<dim>::~RegionMonitor()
{
}





template<uint32_t dim>
void RegionMonitor<dim>::DivideIntegralPropertiesByRegionVolumes( bool doit )
{
    divide_by_volume_ = doit;
    
} // end DivideIntegralPropertiesByRegionVolumes


/** This method is meant to aid those who want to add their own integrated values (done in their own functions)
 *  while keeping the use of the region monitor. NOTE: No checks are made for the existence of the properties or
 *  region names!!
 */

template<uint32_t dim>
void RegionMonitor<dim>::InsertExternallyCalculatedProperty(string property,
                                                            string regionname)
{

    string property_column_header_entry=property+"_"+regionname;
    for(int i = 0; i < property_column_header_entry.length(); i++)
    {
           if( isspace(property_column_header_entry[i]) )
               property_column_header_entry[i] = '_';
    }

    ext_calc_properties_column_headers_.insert(property_column_header_entry);
}

template<uint32_t dim>
void RegionMonitor<dim>::InsertExternallyCalculatedProperty(string property_regionname)
{

    string property_column_header_entry=property_regionname;
    for(int i = 0; i < property_column_header_entry.length(); i++)
    {
           if( isspace(property_column_header_entry[i]) )
               property_column_header_entry[i] = '_';
    }
    ext_calc_properties_column_headers_.insert(property_column_header_entry);
}

template<uint32_t dim>
void RegionMonitor<dim>::InsertPreCalculatedPropertyValue( double time,
                                                           string property,
                                                           string regionname,
                                                           double value )
{
    string property_column_entry=property+"_"+regionname;
    for(int i = 0; i < property_column_entry.length(); i++)
    {
           if( isspace(property_column_entry[i]) )
               property_column_entry[i] = '_';
    }

    // now check that this property column entry actually exists before entering it.
    set<string>::iterator entry=ext_calc_properties_column_headers_.find(property_column_entry);
    if (entry!=ext_calc_properties_column_headers_.end())
        ext_properties_[time].insert(make_pair(property_column_entry,value));
    else
        throw csmp::Exception( ERROR, "RegionMonitor:InsertPreCalculatedPropertyValue()", property_column_entry.c_str(),
                               "Property and region name combination has not been inserted. \n Call InsertExternallyCalculatedProperty() first ");
}





template<uint32_t dim>
void RegionMonitor<dim>::InsertPreCalculatedPropertyValue( double time,
                                                           string property_regionname,
                                                           double value )
{
    string property_column_entry=property_regionname;
    for(int i = 0; i < property_column_entry.length(); i++)
    {
           if( isspace(property_column_entry[i]) )
               property_column_entry[i] = '_';
    }
    // now check that this property column entry actually exists before entering it.
    set<string>::iterator entry=ext_calc_properties_column_headers_.find(property_column_entry);
    if (entry!=ext_calc_properties_column_headers_.end())
        ext_properties_[time].insert(make_pair(property_column_entry,value));
    else
        throw csmp::Exception( ERROR, "RegionMonitor:InsertPreCalculatedPropertyValue()", property_column_entry.c_str(),
                               "Property and region name combination has not been inserted. \n Call InsertExternallyCalculatedProperty() first ");
}





/**
    calculates property integrals and stores these for the current time-step and property name
    if the variable 'thickness' is defined in the PropertyDatabase, it is used to scale
    the element by element integrals.
    
    @attention non-unique regions that bear the name of model boundaries are ignored.
*/
template<uint32_t dim>
void RegionMonitor<dim>::ScalarPropertyIntegrals( const Model<dim>& sg, double time )
{
    double  integral;
    
    // for all properties which shall be integrated over the groups
    for ( typename list<string>::const_iterator
          lit=integral_properties_.begin(); lit!=integral_properties_.end(); lit++ )
    {
        // for all unique regions in the model
        for ( typename map<string,Region<dim> >::const_iterator
              git=sg.UniqueRegionsBegin(); git!=sg.UniqueRegionsEnd(); git++ )
        {
            // integrate the property over the group
            if ( include_thickness_attribute_ )
            {
                if ( pore_volume_integral_ )
                    integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), true );
                else
                    integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), false );
            }
            else {
                if ( pore_volume_integral_ )
                    integral = (*git).second.VolumeIntegral( (*lit).c_str(), true );
                else
                    integral = (*git).second.VolumeIntegral( (*lit).c_str(), false );
            }
            // normalization
            if ( divide_by_volume_ ) {
                typename map<string,std::pair<double,double> >::const_iterator
                        gr_it=group_specs_.find( (*git).first );
                // integral value gets divided by volume/area represented by group
                integral /= (*gr_it).second.first;
            }
            // global-time  property    group-name
            integrals_[time][(*lit)][(*git).first] = integral;
        }

        // for all regions in the model (unless they are boundaries)
        for ( typename map<string,Region<dim> >::const_iterator
              git=sg.RegionsBegin(); git!=sg.RegionsEnd(); git++ )
          if ( parseBoundary((*git).first) == NOT and (*git).first !="Model" )
            {
                // integrate the property over the group
                if ( include_thickness_attribute_ )
                {
                    if ( pore_volume_integral_ )
                        integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), true );
                    else
                        integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), false );
                }
                else {
                    if ( pore_volume_integral_ )
                        integral = (*git).second.VolumeIntegral( (*lit).c_str(), true );
                    else
                        integral = (*git).second.VolumeIntegral( (*lit).c_str(), false );
                }
                // normalization
                if ( divide_by_volume_ ) {
                    typename map<string,pair<double,double> >::const_iterator
                            gr_it=group_specs_.find( (*git).first );
                    // integral value gets divided by volume/area represented by group
                    integral /= (*gr_it).second.first;
                }
                // global-time  property    group-name
                integrals_[time][(*lit)][(*git).first] = integral;
            }
    }

} // end ScalarPropertyIntegrals




/**
    Records the ranges of the target variables in each region.

    @attention non-unique regions that bear the name of model boundaries are ignored.
*/
template<uint32_t dim>
void RegionMonitor<dim>::ScalarPropertyRanges( const Model<dim>& sg, double time )
{
    double  rmin, rmax;
    
    // for all properties for which ranges shall be monitored
    for ( typename list<string>::const_iterator
          lit=range_properties_.begin(); lit!=range_properties_.end(); lit++ )
    {
        // unique regions
        for ( typename map<string,Region<dim> >::const_iterator
              git=sg.UniqueRegionsBegin(); git!=sg.UniqueRegionsEnd(); git++ )
        {
            (*git).second.MinMaxOf( (*lit).c_str(), rmin, rmax );
            ranges_[time][(*lit)][(*git).first] = make_pair(rmin,rmax);
        }
        // other regions
        for ( typename map<string,Region<dim> >::const_iterator
              git=sg.RegionsBegin(); git!=sg.RegionsEnd(); git++ )
          if ( parseBoundary((*git).first) == NOT and (*git).first !="Model" ) {
               (*git).second.MinMaxOf( (*lit).c_str(), rmin, rmax );
               ranges_[time][(*lit)][(*git).first] = make_pair(rmin,rmax);
            }
    }

} // end ScalarPropertyRanges




template<uint32_t dim>
void RegionMonitor<dim>::Reset() // zap all recorded values
{
    integral_properties_.clear();
    range_properties_.clear();
    integrals_.clear();
    ranges_.clear();
    ext_properties_.clear();

} // end Reset



template<uint32_t dim>
void RegionMonitor<dim>::EraseData() // zap all recorded values
{
    integrals_.clear();
    ranges_.clear();
    ext_properties_.clear();
    ext_calc_properties_column_headers_.clear();

} // end EraseData





template<uint32_t dim>
void RegionMonitor<dim>::Out( const char* text_file ) const
{
    string  file(text_file);
    file +=".txt";
    
    ofstream  ofs( file.c_str() );

    
    if ( !integrals_.empty() )
        ofs <<"'"<< file <<"' region data recorded up to time, t = "<< (*integrals_.rbegin()).first <<" seconds.";
    else
        ofs <<"'"<< file <<"' region data recorded up to time, t = "<< (*ranges_.rbegin()).first <<" seconds.";

    // outputting geometric properties of groups
    // -----------------------------------------
    ofs <<"\nregion_name \t volume \t surface area "<< endl;
    for ( typename map<string,pair<double,double> >::const_iterator
          gr_it=group_specs_.begin(); gr_it!=group_specs_.end(); gr_it++ ) {
        ofs << (*gr_it).first <<"\t"<< (*gr_it).second.first <<"\t"<< (*gr_it).second.second;
        ofs << endl;
    }

    // writing integrated property data blocks
    // ---------------------------------------
    // for all properties integrated over the groups
    ofs <<"\nIntegrated values of properties: "<< endl;
    for ( typename list<string>::const_iterator
          iit=integral_properties_.begin(); iit!=integral_properties_.end(); iit++ )
    {
        // name of monitored property
        ofs << endl << (*iit) << endl;

        // the column headers: time, group1, group2...
        ofs <<" time\t";
        for ( typename map<string,pair<double,double> >::const_iterator
              gr_it=group_specs_.begin(); gr_it!=group_specs_.end(); gr_it++ ){
             ofs << (*gr_it).first <<"\t";
          }
        ofs << endl;

        // for all timesteps
        //         model_time,  property, groupname, property value
        for ( typename map<double,map<string,map<string,double> > >::const_iterator
              int_it=integrals_.begin(); int_it!=integrals_.end(); int_it++ ) {
            // print the time at the beginning of the row
            ofs << (*int_it).first <<"\t";
            // find property record in the map
            typename map<string,map<string,double> >::const_iterator iprop_it=(*int_it).second.find(*iit);
            assert( iprop_it != (*int_it).second.end() );
            // for all groups print the recorded property values
            for ( typename map<string,double>::const_iterator
                  ipit=(*iprop_it).second.begin(); ipit!=(*iprop_it).second.end(); ipit++ )
                // writing values to file
                ofs << (*ipit).second <<"\t";
            // move to the next line
            ofs << endl;
        }
        ofs << endl;
    }

    // writing property-range data blocks
    // ----------------------------------
    // for all property ranges determined for the groups
    ofs <<"\nValue ranges of properties: "<< endl;
    for ( typename list<string>::const_iterator
          iit=range_properties_.begin(); iit!=range_properties_.end(); iit++ )
    {
        // name of monitored property
        ofs << endl << (*iit) << endl;

        // the column headers: time, group1, group2...
        ofs <<" time\t\t";
        for ( typename map<string,pair<double,double> >::const_iterator
              gr_it=group_specs_.begin(); gr_it!=group_specs_.end(); gr_it++ )
            ofs << (*gr_it).first <<"\t\t";
        ofs << endl;

        // for all timesteps
        //         model_time,  property, groupname, property ranges
        for ( typename map<double,map<string,map<string,pair<double,double> > > >::const_iterator
              ran_it=ranges_.begin(); ran_it!=ranges_.end(); ran_it++ ) {
            // the time
            ofs << (*ran_it).first <<"\t";
            // the property data-record searched in the map
            typename map<string,map<string,pair<double,double> > >::const_iterator rprop_it=(*ran_it).second.find(*iit);
            if ( rprop_it != (*ran_it).second.end() )
                // for all property values recorded for the groups
                for ( typename map<string,pair<double,double> >::const_iterator
                      rpit=(*rprop_it).second.begin(); rpit!=(*rprop_it).second.end(); rpit++ )
                    // writing values to file
                    ofs << (*rpit).second.first <<"\t" << (*rpit).second.second <<"\t";
            ofs << endl;
        }
        ofs << endl;
    }

    // for all properties externally calculated over the groups
    // the column headers: time, property1_group1, property2_group2, property3_group6, etc...
    if (!ext_calc_properties_column_headers_.empty()){
        ofstream  ofsext( (string(text_file)+"_ext.txt").c_str() );
        ofsext <<"(1)time\t";
        size_t ci(2);
        for ( typename set<string>::const_iterator
              cit=ext_calc_properties_column_headers_.begin(); cit!=ext_calc_properties_column_headers_.end(); cit++ )
        {
            string column=*cit;
            ofsext<<"("<<ci<<")"<<column<<"\t";
            ci++;
        }
        ofsext << endl;

        ofsext.width(10);
        ofsext.setf(ios::scientific);
        double value;
        // now the time information, in ascending order (downwards in the file).
        for ( typename std::map<double,std::map<std::string,double> >::const_iterator
              extp_it=ext_properties_.begin(); extp_it!=ext_properties_.end(); extp_it++ ) {
            ofsext<<extp_it->first;
            for ( typename set<string>::const_iterator
                  cit=ext_calc_properties_column_headers_.begin(); cit!=ext_calc_properties_column_headers_.end(); cit++ )
            {
                value=extp_it->second.at(*cit);
                if (!isnan(value))
                    ofsext<<"\t"<<value;
                else
                    ofsext<<"\t"<<0.0;
            }
            ofsext << endl;
        }
    }

} // end Out



/// P. Lang: loops over all elements of all regions of model and returns true if IsVolumeElement()
template<uint32_t dim>
bool  RegionMonitor<dim>::HasVolumeElements( const Model<dim>& mref ) const
{
    if( dim == 1U || dim == 2U ) return false;

    const Region<dim>&  rref( mref.Region("Model") );

    const auto  elementsEnd = rref.CellsEnd();
    for ( auto eit=rref.CellsBegin(); eit!=elementsEnd; ++eit )
        if ( (*eit)->IsVolumeElement() )
          return true;

    return false;
}


template class RegionMonitor<1U>;
template class RegionMonitor<2U>;
template class RegionMonitor<3U>;


} // end namespace csmp














