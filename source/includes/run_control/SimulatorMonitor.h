#ifndef SIMULATOR_MONITOR_H
#define SIMULATOR_MONITOR_H
#include "Model.h"

namespace csmp {

template<size_t dim> class Model; 

template<size_t dim>
class SimulatorMonitor {
public:

    SimulatorMonitor(std::string text_file,
                     Model<dim> &model,
                     bool restart,
                     const std::list<std::string>& variables_to_monitor_integrals,
                     const std::list<std::string>& variables_to_monitor_ranges,
                     const std::list<std::string>& variables_to_monitor_single_values,
                     const std::list<std::string>& subdomains_to_integrate,
                     const std::list<std::string>& subdomains_to_measure_ranges,
                     const std::list<std::string>& subdomains_to_calculate_dimension,
                     const std::list<std::string>& subdomains_to_calculate_perimeter,
                     const std::list<std::string>& to_divide_by_region_porevolume,
                     const std::list<std::string>& to_multiply_integrand_by_porosity,
                     bool verbose = true);

    ~SimulatorMonitor();

    /// calculates property integrals and stores these for the current time-step and property name
    void ScalarPropertyIntegrals(std::vector<double64>& rowdata);
    void ScalarPropertyRanges(std::vector<double64>& rowdata);
    void ScalarPropertyValues(std::vector<double64>& rowdata);
    void CalculateDimensionsAndPerimeters(std::vector<double64>& rowdata);
    void Monitor();

    /** adds in values in a separate container, to be output to a separate file
        with the same prefix name as the other monitoring file, except with the "_ext" added to its name.
    */
    void InsertValueHeader(std::string property_regionname , std::vector<size_t>& indexes);

    void ClearAllHeadersAndValues(); // zap all recorded values and property names
    void EraseDataValuesOnly(); // zap all recorded values
    


    bool DataExists(){return !values_.empty();}
    //bool EmptyIntegrals(){return integrals_.empty();}
    //bool EmptyRanges(){return ranges_.empty();}

    //size_t RangesSize(){return ranges_.size();}
    //size_t IntegralsSize(){return integrals_.size();}

    bool Verbose() {return this->verbose_;}
    void Verbose(bool verbose) {this->verbose_=verbose;}

    void ReadOldMonitoringData(std::string file_name );
//    void ReadDualValueOldMonitoringData(std::string file_name );
    void ReadModelTime(std::vector<double64>& rowdata);

private:
    // iternal helper for surface calculation
    bool  HasVolumeElements( const Model<dim>& ) const;
    void CheckForDuplicateHeaders();
    void  Out();
    std::string UScoreForSpace(std::string text);
    bool is_number(const std::string& s);
    SimulatorMonitor();

    Model<dim>& mref_;
    ErrorHandler& error_handler_;
    std::string file_prefix_;

    std::list<std::string>  integrated_property_names_;
    std::list<std::string>  ranged_property_names_;
    std::list<std::string>  single_value_property_names_; // for externally calculated monitored properties
    std::list<std::string>  dimensional_property_names_; // for areas and volume monitoring.

    std::list<std::string> subdomains_to_calculate_perimeter_,subdomains_to_calculate_dimension_;

    std::vector<size_t> int_header_indexes_;
    std::vector<size_t> range_header_indexes_;
    std::vector<size_t> single_value_header_indexes_;
    std::vector<size_t> dimensional_header_indexes_;
    std::vector<size_t> model_time_indexes_;

    std::vector<std::string>  values_column_headers_; // single value
    std::set<std::string> divide_by_region_volume_,multiply_integrand_with_porosity_;

    bool include_thickness_attribute_;    \
    bool verbose_,first_monitor_call_;

    std::vector<double64> rowdata_;
    double64 last_requested_monitor_time_, current_requested_monitor_time_;

    //    model_time,          property_groupname,   single property value
    //std::map<double64,std::multimap<std::string,double64 > >   single_values_;
    std::vector< std::vector<double64 > >   values_;

    //    model_time,          property_groupname,   dual value pair
    //std::map<double64,std::multimap<std::string,std::pair<double64,double64> > >   dual_values_;
    //std::map<double64,pair<std::string,std::vector<pair <double64,double64 > > >   dual_values_;

    //    model_time,          groupname,   volume,surface area pair
    //std::map<double64,std::multimap<std::string,std::pair<double64,double64> > >   volume_and_area_;

};

} // csmp

#endif // SIMULATION_MONITOR_H
