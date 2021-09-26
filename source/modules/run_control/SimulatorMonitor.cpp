#include "SimulatorMonitor.h"
#include "ErrorHandler.h"
#include "Node.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "PL_Utilities.h"
#include "CSMP_highLevelUtilities.h"

#include "Boundary.h"

using namespace std;

namespace csmp {

/** @author Julian, July 2014)
 * This monitor class is different from the RegionMonitor, mainly in the way it organizes and outputs information.
 * It will output information column wise, and separated into 3 files: ranges_monitor.txt, integral_monitor.txt, and single_value_monitor.txt
 * The names of the files indicate what goes into them.
 * The monitor class, like the RegionMonitor, only outputs its own calculations of integrals and ranges to the respective files.
 * It is up to the user to pre-calculate other values outside of this class.
 * Finally, this monitor will output data in column form, ready to be processed. It will also append existing files, without having to re-output all accumulated data
 * in the monitor (the RegionMonitor does this).
 * As a default behaviour, if the file exists, it will be overwritten.
 *
 * @todo 1) implement automatic writing of monitor data into model for re-starting.
 * 2) the monitor will search through the lines of an existing monitor file until it finds a position to begin writing (a simulated time smaller than the current one).
 * 3) it will also check that the existing monitor file has the same columns
 * 4) points 2) and 3) are a secondary feature if we can save all data to a restart file directly!!!
 */
template<size_t dim>
SimulatorMonitor<dim>::SimulatorMonitor(string file_prefix,
                                        Model<dim>& model,
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
                                        bool verbose)


    : file_prefix_(file_prefix),
      mref_(model),
      integrated_property_names_(variables_to_monitor_integrals),
      ranged_property_names_(variables_to_monitor_ranges),
      single_value_property_names_(variables_to_monitor_single_values),
      subdomains_to_calculate_dimension_(subdomains_to_calculate_dimension),
      subdomains_to_calculate_perimeter_(subdomains_to_calculate_perimeter),
      first_monitor_call_(true),
      include_thickness_attribute_(model.Database().IsDefined("thickness")),
      verbose_(verbose)
{
    // check that "model time" is defined.
    if (!mref_.Database().IsDefined("model time"))
        ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", "'model time' needs to be defined as a MODEL variable",
                               "the simulator monitor needs to know the model time!");

    this->InsertValueHeader("model time",dimensional_header_indexes_);
    if (!integrated_property_names_.empty())
        this->InsertValueHeader("model time",int_header_indexes_);
    if (!ranged_property_names_.empty())
        this->InsertValueHeader("model time",range_header_indexes_);
    this->InsertValueHeader("model time",single_value_header_indexes_);

    for (auto lit = to_divide_by_region_porevolume.begin();lit != to_divide_by_region_porevolume.end();lit++)
        divide_by_region_volume_.insert(*lit);

    for (auto lit = to_multiply_integrand_by_porosity.begin();lit != to_multiply_integrand_by_porosity.end();lit++)
        multiply_integrand_with_porosity_.insert(*lit);

    // checking whether target properties are suitable for integration
    for ( typename list<string>::const_iterator
          it=integrated_property_names_.begin(); it!=integrated_property_names_.end(); it++ ) {
        if ( !mref_.Database().IsDefined( (*it).c_str() ) )
            ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", (*it).c_str(),
                                   "integral property is not defined in the database (file)");

        if ( mref_.Database().Type( (*it).c_str() ) != SCALAR )
            ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", (*it).c_str(),
                                   "integral property must be a scalar property");

        if ( mref_.Database().Placement( (*it).c_str() ) == INTER_FACE )
            ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", (*it).c_str(),
                                   "properties placed on the FACE cannot be integrated over the Region volume");
    }

    // checking range properties
    for ( typename list<string>::const_iterator
          it=ranged_property_names_.begin(); it!=ranged_property_names_.end(); it++ )
        if ( !mref_.Database().IsDefined( (*it).c_str() ) )
            ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", (*it).c_str(),
                                   "range property is not defined in the database (file)");
    // ------------------------------------------------------------
    // FOR MONITORING VOLUMES AND AREAS OF REGIONS
    // Some might stay constant, but others might change and can be monitored.
    dimensional_property_names_.push_back("_VOID_DIM");
    dimensional_property_names_.push_back("_VOID_PERI");
    dimensional_property_names_.push_back("_BULK_DIM");
    dimensional_property_names_.push_back("_BULK_PERI");

    for ( auto sbdmit = subdomains_to_calculate_dimension_.begin();sbdmit != subdomains_to_calculate_dimension_.end();sbdmit++){
        this->InsertValueHeader(*sbdmit+"_BULK_DIM",dimensional_header_indexes_);
        this->InsertValueHeader(*sbdmit+"_VOID_DIM",dimensional_header_indexes_);
    }

    for ( auto sbdmit = subdomains_to_calculate_perimeter_.begin();sbdmit != subdomains_to_calculate_perimeter_.end();sbdmit++){
        this->InsertValueHeader(*sbdmit+"_BULK_PERI",dimensional_header_indexes_);
        this->InsertValueHeader(*sbdmit+"_VOID_PERI",dimensional_header_indexes_);
    }

    // ------------------------------------------------------------
    // FOR MONITORING INTEGRALS OF SCALAR PROPERTIES ON REGIONS
    this->integrated_property_names_.sort();
    this->integrated_property_names_.unique();
    for (auto lit = integrated_property_names_.begin(); lit != integrated_property_names_.end();lit++)
        for ( auto sbdmit = subdomains_to_integrate.begin();sbdmit != subdomains_to_integrate.end();sbdmit++)
            this->InsertValueHeader((*lit)+"_"+*sbdmit+"_INT",int_header_indexes_);

    // ------------------------------------------------------------
    // FOR MONITORING RANGES OF SCALAR PROPERTIES ON REGIONS
    this->ranged_property_names_.sort();
    this->ranged_property_names_.unique();
    for (auto lit = ranged_property_names_.begin(); lit != ranged_property_names_.end();lit++){
        for ( auto sbdmit = subdomains_to_measure_ranges.begin();sbdmit != subdomains_to_measure_ranges.end();sbdmit++){
            this->InsertValueHeader((*lit)+"_"+*sbdmit+"_RMIN",range_header_indexes_);
            this->InsertValueHeader((*lit)+"_"+*sbdmit+"_RMAX",range_header_indexes_);
        }
    }

    // ------------------------------------------------------------
    // FOR MONITORING SINGLE VALUES EXTERNALLY CALCULATED ON ANY MODEL SUBDMAIN (OR THE MODEL)
    // checking whether target properties are suitable

    // run time and timestep are standard.
    this->single_value_property_names_.push_back("run time");
    this->single_value_property_names_.push_back("timestep");
    this->single_value_property_names_.push_back("interval");
    this->single_value_property_names_.sort();
    this->single_value_property_names_.unique();

    for ( typename list<string>::const_iterator
          lit=single_value_property_names_.begin(); lit!=single_value_property_names_.end(); lit++ ) {
        if ( !mref_.Database().IsDefined( (*lit).c_str() ) )
            ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", (*lit).c_str(),
                                   "integral property is not defined in the database (file)");

        if ( mref_.Database().Type( (*lit).c_str() ) != SCALAR )
            ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", (*lit).c_str(),
                                   "integral property must be a scalar property");

        set<PLACEMENT> placements_accepted;
        placements_accepted.insert(MODEL);
        //placements_accepted.insert(REGION);
        //placements_accepted.insert(BOUNDARY);
        //placements_accepted.insert(SPLIT_BOUNDARY);

        string regionname="Model";
        if ( placements_accepted.find(mref_.Database().Placement( (*lit).c_str() )) != placements_accepted.end() )
            this->InsertValueHeader((*lit)+"_"+regionname,single_value_header_indexes_);
        else
            ErrorHandler::Instance().notice( FATAL_ERROR, "SimulatorMonitor(constructor)", (*lit).c_str(),
                                   " Externally calculated properties should be placed on MODEL (for now)");
    }

    //CheckForDuplicateHeaders();

    // if restarting, read any data from them to later append
    if (restart){
        string  file_integrals(file_prefix);
        string  file_ranges(file_prefix);
        string  file_single_value(file_prefix);
        string  file_areas_volumes(file_prefix);

        file_integrals +="_integrals_monitor.txt";
        file_ranges +="_ranges_monitor.txt";
        file_single_value +="_single_value_monitor.txt";
        file_areas_volumes+="_areas_and_volumes_monitor.txt";

        this->ReadOldMonitoringData(file_integrals);
        this->ReadOldMonitoringData(file_single_value);
        this->ReadOldMonitoringData(file_ranges);
        this->ReadOldMonitoringData(file_areas_volumes);
    }
//    cin.get();

} // end constructor

template<size_t dim>
void SimulatorMonitor<dim>::CheckForDuplicateHeaders(){
    set<string> check;
    size_t i = 0;
    for (auto it = values_column_headers_.begin(); it!= values_column_headers_.end();it++){
        cout<<i<<" value header: "<<*it<<endl;
        if (check.find(*it)==check.end()){
            check.insert(*it);

        }
        i++;
//        else
//            ErrorHandler::Instance().notice(FATAL_ERROR,"SimulatorMonitor::CheckForDuplicateHeaders()",(*it).c_str(),"Found duplicate monitor header.");
    }
    for (auto it = single_value_header_indexes_.begin();it != single_value_header_indexes_.end();it++){
        cout<<" pre-calc header "<<*it<<endl;
    }
    cin.get();
}

template<size_t dim>
void SimulatorMonitor<dim>::ReadOldMonitoringData(string file_name ){

    ifstream  ifs(file_name);
    string text_line,token;
    vector<vector<double64> > data;
    vector<string> columnheaders;

    static bool first_call(true);

    if (ifs.is_open() &&  !isInputFileEmpty(ifs)){
        std::getline( ifs, text_line );
        // read the headers. (first line of every monitoring file output by this class.)

        stringstream iss_h(text_line);

        cout<<" SimulatorMonitor<dim> Reading existing headers from file '"<<file_name<<"'"<<endl;
        cout.flush();
//        cin.get();
        while (getline(iss_h,token,'\t'))
            columnheaders.push_back(token);

        bool found=false;
        vector<string>::iterator mt;
        for (vector<string>::iterator sit = columnheaders.begin(); sit!=columnheaders.end(); sit++)
        {
            if ( sit->find("model_time")!=std::string::npos){
                found=true;
                mt=sit;
            }
        }

        if (found){ // found model time data.
            size_t mtime_pos = mt - columnheaders.begin();
            std::getline( ifs, text_line  );
            // now read the data associated to each column header
            while (!ifs.eof() && (text_line[0] != '#')){

                vector<string> listoftokens;
                vector<double64> dataline;
                stringstream iss(text_line);
                while (getline(iss,token,'\t'))
                    listoftokens.push_back(token);

                // if model time is greater than current simulation time, break out.
                if (std::stod(listoftokens[mtime_pos]) >= mref_.Read(mref_.Database().StorageKey("model time"))){
                    cout<<" WARNING: model time in monitoring file exceeds model time saved in the CSMP model"<<std::stod(listoftokens[mtime_pos])<<">="<<mref_.Read(mref_.Database().StorageKey("model time"))<<endl;
                    cout<<" Data beyond this value of simulated time will be discarded."<<endl;
                    cout.flush();
                    break;
                }

                /// if the size of the data read on one line does not match the number of headers, skip to the next data line.
                if (listoftokens.size()!= columnheaders.size()){
                    cout<<" WARNING: number of read-in column headers mismatch!"<<endl;
                    cout.flush();
                    continue;
                }
                size_t n(0);
                for (vector<string>::iterator lit = listoftokens.begin(); lit!=listoftokens.end();lit++){
                    dataline.push_back(std::stod(*lit));
                    n++;
                }
                               
                // now feed the data line to the overall container only if the time value is older than current model time.
                data.push_back(dataline);
                std::getline( ifs, text_line  );
            }

            if (first_call){
                values_.resize(data.size());
                for (size_t i = 0 ; i < values_.size();i++){
                    values_[i].resize(values_column_headers_.size());
                    fill(values_[i].begin(), values_[i].end(),0.0);
                }
                first_call=false;
            }

            // this check could be done during the reading process.  But since we would like to keep the old data throughout
            // the simulation anyway, we can do it separately.
            vector<size_t> old_index_to_new_index(columnheaders.size());
            std::fill(old_index_to_new_index.begin(),old_index_to_new_index.end(),values_column_headers_.size());

            for (size_t i = 0; i < columnheaders.size() ; i++)
                for (size_t j = 0 ; j < values_column_headers_.size(); j++)
                {
                    if ( columnheaders[i].find(string(values_column_headers_[j]+"("))!=std::string::npos){
                        old_index_to_new_index[i]=j;
                        break;
                    }
                }

            std::fill(rowdata_.begin(), rowdata_.end(),0.0);
            for (size_t i = 0 ; i < data.size() ; i++){
                for (size_t j = 0; j < columnheaders.size() ; j++)
                    if (j < values_column_headers_.size())
                        values_[i][old_index_to_new_index[j]]=data[i][j];
            }

        }
    }
}

template<size_t dim>
string SimulatorMonitor<dim>::UScoreForSpace(string text)
{
    string property_column_header_entry=text;
    for(int i = 0; i < property_column_header_entry.length(); i++)
    {
        if( isspace(property_column_header_entry[i]) )
            property_column_header_entry[i] = '_';
    }
    return property_column_header_entry;
}

template<size_t dim>
void SimulatorMonitor<dim>::InsertValueHeader(string property_regionname, vector<size_t> &indexes)
{
    string property_column_header_entry=UScoreForSpace(property_regionname);

    vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),property_column_header_entry);
    if (ith==values_column_headers_.end() || property_column_header_entry=="model_time"){
        size_t j = ith-values_column_headers_.begin();
        indexes.push_back(j);
        values_column_headers_.push_back(property_column_header_entry);
        if (this->Verbose()) cout<<" inserting value header : "<<property_column_header_entry<<endl;
        if (property_column_header_entry=="model_time")
            this->model_time_indexes_.push_back(j);
    }
    else
        ErrorHandler::Instance().notice(FATAL_ERROR,"SimulatorMonitor::InsertValueHeader()",property_column_header_entry.c_str(),"Found duplicate monitor header.");
}

/**
    calculates property integrals and stores these for the current time-step and property name
    if the variable 'thickness' is defined in the PropertyDatabase, it is used to scale
    the element by element integrals.
*/
template<size_t dim>
void SimulatorMonitor<dim>::ScalarPropertyIntegrals(vector<double64>& rowdata)
{
    double64  integral;
    
    // for all properties which shall be integrated over the groups
    for ( typename list<string>::const_iterator
          lit=integrated_property_names_.begin(); lit!=integrated_property_names_.end(); lit++ )
    {
        string read_value_key;
        // for all unique regions in the model
        for ( typename map<string,Region<dim> >::const_iterator
              git=mref_.UniqueRegionsBegin(); git!=mref_.UniqueRegionsEnd(); git++ )
        {
            string regionname=(*git).first;
            if ( include_thickness_attribute_ ) {
                if ( multiply_integrand_with_porosity_.find(*lit)!=multiply_integrand_with_porosity_.end() )
                    integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), true );
                else
                    integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), false );
            }
            else {
                if ( multiply_integrand_with_porosity_.find(*lit)!=multiply_integrand_with_porosity_.end() )
                    integral = (*git).second.VolumeIntegral( (*lit).c_str(), true );
                else
                    integral = (*git).second.VolumeIntegral( (*lit).c_str(), false );
            }

            if (divide_by_region_volume_.find(*lit)!=divide_by_region_volume_.end())
            {
                read_value_key=UScoreForSpace(regionname+"_BULK_DIM");
                vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
                if (ith!=values_column_headers_.end()){
                    size_t j = ith - values_column_headers_.begin();
                    integral /= rowdata[j];
                }
                else { // slower, but gets the job done
                    integral /= (*git).second.Volume(false);
                }
            }

            read_value_key=UScoreForSpace((*lit)+"_"+regionname+"_INT");
            vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j]=integral;
            }
        }

        // for all regions in the model
        for ( typename map<string,Region<dim> >::const_iterator
              git=mref_.RegionsBegin(); git!=mref_.RegionsEnd(); git++ )
        {
            string regionname=(*git).first;

            if ( include_thickness_attribute_ ) {
                if ( multiply_integrand_with_porosity_.find(*lit)!=multiply_integrand_with_porosity_.end() )
                    integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), true );
                else
                    integral = (*git).second.VolumeIntegral_x_Thickness( (*lit).c_str(), false );
            }
            else {
                if ( multiply_integrand_with_porosity_.find(*lit)!=multiply_integrand_with_porosity_.end() )
                    integral = (*git).second.VolumeIntegral( (*lit).c_str(), true );
                else
                    integral = (*git).second.VolumeIntegral( (*lit).c_str(), false );
            }

            if (divide_by_region_volume_.find(*lit)!=divide_by_region_volume_.end())
            {
                read_value_key=UScoreForSpace(regionname+"_BULK_DIM");
                vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
                if (ith!=values_column_headers_.end()){
                    size_t j = ith - values_column_headers_.begin();
                    integral /= rowdata[j];
                }
                else { // slower, but gets the job done
                    integral /= (*git).second.Volume(false);
                }
            }

            read_value_key=UScoreForSpace((*lit)+"_"+regionname+"_INT");
            vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j]=integral;
            }
        }
    }
} // end ScalarPropertyIntegrals

template<size_t dim>
void SimulatorMonitor<dim>::ScalarPropertyRanges( vector<double64>& rowdata)
{
    double64  rmin, rmax;
    // for all properties for which ranges shall be monitored
    for ( typename list<string>::const_iterator
          lit=ranged_property_names_.begin(); lit!=ranged_property_names_.end(); lit++ )
    {
        string read_value_key;
        // unique regions
        for ( typename map<string,Region<dim> >::const_iterator
              git=mref_.UniqueRegionsBegin(); git!=mref_.UniqueRegionsEnd(); git++ )
        {
            string regionname=(*git).first;
            (*git).second.MinMaxOf( (*lit).c_str(), rmin, rmax );
            read_value_key=UScoreForSpace((*lit)+"_"+regionname+"_RMIN");
            vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j] = rmin;
            }
            read_value_key=UScoreForSpace((*lit)+"_"+regionname+"_RMAX");
            ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j] = rmax;
            }
        }
        // other regions
        for ( typename map<string,Region<dim> >::const_iterator
              git=mref_.RegionsBegin(); git!=mref_.RegionsEnd(); git++ )
        {
            string regionname=(*git).first;
            (*git).second.MinMaxOf( (*lit).c_str(), rmin, rmax );
            read_value_key=UScoreForSpace((*lit)+"_"+regionname+"_RMIN");
            vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j] = rmin;
            }
            read_value_key=UScoreForSpace((*lit)+"_"+regionname+"_RMAX");
            ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j] = rmax;
            }
        }
    }
} // end ScalarPropertyRanges

/// this method reads values directly from the model
///
template<size_t dim>
void SimulatorMonitor<dim>::ScalarPropertyValues(vector<double64>& rowdata)
{
    string read_value_key;
    string regionname="Model";
    for (auto lit = single_value_property_names_.begin(); lit!=single_value_property_names_.end();lit++){
        read_value_key=UScoreForSpace((*lit)+"_"+regionname);

        vector<string>::iterator ith=find(values_column_headers_.begin(),values_column_headers_.end(),read_value_key);
        if (ith!=values_column_headers_.end()){
            size_t j = ith - values_column_headers_.begin();
            double64 value=mref_.Read(mref_.Database().StorageKey((*lit).c_str()));
            rowdata[j]=value;
        }
    }
}

template<size_t dim>
void SimulatorMonitor<dim>::CalculateDimensionsAndPerimeters(vector<double64>& rowdata){

    bool hasVolumeElements( HasVolumeElements( mref_ ) );

    for (auto sbdmit = subdomains_to_calculate_dimension_.begin() ; sbdmit != subdomains_to_calculate_dimension_.end();sbdmit++){

        double64 void_dim(0.0);
        double64 bulk_dim(0.0);

        if (mref_.ContainsRegion(sbdmit->c_str())){
            Region<dim> & rref = mref_.Region(sbdmit->c_str());
            void_dim = rref.Volume(true);
            bulk_dim = rref.Volume(false);
        }

        if (mref_.ContainsBoundary(sbdmit->c_str())){
            Boundary<dim> & bref = mref_.Boundary(sbdmit->c_str());
            void_dim = bref.Area(); /// @todo multiplication with porosity is not supported as with regions. To be resolved later on.
            bulk_dim = bref.Area(); /// @todo multiplication with porosity is not supported as with regions. To be resolved later on.
        }

        map<string,double64> read_value_keys;
        read_value_keys.insert(make_pair(UScoreForSpace(*sbdmit+"_VOID_DIM"),void_dim));
        read_value_keys.insert(make_pair(UScoreForSpace(*sbdmit+"_BULK_DIM"),bulk_dim));

        for (auto vit = read_value_keys.begin() ; vit != read_value_keys.end() ; vit++){
            auto ith=find(values_column_headers_.begin(),values_column_headers_.end(),vit->first);
            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j] = vit->second;
            }
        }
    }

    for (auto sbdmit = subdomains_to_calculate_perimeter_.begin() ; sbdmit != subdomains_to_calculate_perimeter_.end();sbdmit++){

        double64 void_peri(0.0);
        double64 bulk_peri(0.0);

        if (mref_.ContainsRegion(sbdmit->c_str())){
            Region<dim> & rref = mref_.Region(sbdmit->c_str());
            void_peri = ( !hasVolumeElements ) ? rref.Volume(true) : rref.SurfaceArea();
            bulk_peri = ( !hasVolumeElements ) ? rref.Volume(false) : rref.SurfaceArea();
        }

        if (mref_.ContainsBoundary(sbdmit->c_str())){
            Boundary<dim> & bref = mref_.Boundary(sbdmit->c_str());
            void_peri = ( !hasVolumeElements ) ? bref.Area() : bref.Area();///@todo multiplication with porosity is not supported as with regions. To be resolved later on.
            bulk_peri = ( !hasVolumeElements ) ? bref.Area() : bref.Area();///@todo multiplication with porosity is not supported as with regions. To be resolved later on.
        }

        map<string,double64> read_value_keys;
        read_value_keys.insert(make_pair(UScoreForSpace(*sbdmit+"_VOID_PERI"),void_peri));
        read_value_keys.insert(make_pair(UScoreForSpace(*sbdmit+"_BULK_PERI"),bulk_peri));

        for (auto vit = read_value_keys.begin() ; vit != read_value_keys.end() ; vit++){
            auto ith=find(values_column_headers_.begin(),values_column_headers_.end(),vit->first);

            if (ith!=values_column_headers_.end()){
                size_t j = ith - values_column_headers_.begin();
                rowdata[j] = vit->second;
            }
        }
    }
}

///
template<size_t dim>
void SimulatorMonitor<dim>::Monitor()
{
    if (first_monitor_call_){
        rowdata_.resize(values_column_headers_.size());
        for (size_t i = 0 ; i< rowdata_.size();i++)
            rowdata_[i]=0.0;
        this->CalculateDimensionsAndPerimeters(rowdata_); // for now, volumes and areas are calculated only once.
    }
    this->ReadModelTime(rowdata_);
    if (current_requested_monitor_time_ > last_requested_monitor_time_ || first_monitor_call_){
        this->ScalarPropertyIntegrals(rowdata_);
        this->ScalarPropertyRanges(rowdata_);
        this->ScalarPropertyValues(rowdata_);
        this->values_.push_back(rowdata_);
        this->Out();
        last_requested_monitor_time_=current_requested_monitor_time_;
    }
    first_monitor_call_=false;
}

template<size_t dim>
void SimulatorMonitor<dim>::ReadModelTime(vector<double64>& rowdata)
{
    double64 time = mref_.Read(mref_.Database().StorageKey("model time"));

    current_requested_monitor_time_=time;
    for (size_t i = 0 ; i < model_time_indexes_.size();i++) {
        rowdata[model_time_indexes_[i]]=time;
    }
    if (first_monitor_call_) last_requested_monitor_time_=time;
}

template<size_t dim>
void SimulatorMonitor<dim>::ClearAllHeadersAndValues() // zap all recorded values
{
    integrated_property_names_.clear();
    ranged_property_names_.clear();
    single_value_property_names_.clear();
    values_column_headers_.clear();
    values_.clear();

} // end ClearAllHeadersAndValues

template<size_t dim>
void SimulatorMonitor<dim>::EraseDataValuesOnly()
{
    values_.clear();
} // end EraseDataValuesOnly

template<size_t dim>
void SimulatorMonitor<dim>::Out()
{
    string  file_integrals(file_prefix_);
    string  file_ranges(file_prefix_);
    string  file_single_value(file_prefix_);
    string  file_areas_volumes(file_prefix_);

    file_integrals +="_integrals_monitor.txt";
    file_ranges +="_ranges_monitor.txt";
    file_single_value +="_single_value_monitor.txt";
    file_areas_volumes+="_areas_and_volumes_monitor.txt";

    // if on first call, output all headers on every file
    if (first_monitor_call_){
        if (!integrated_property_names_.empty()){
            ofstream  ofs_integrals( file_integrals, std::ofstream::out );
            size_t counter(1);
            for (auto tokit = int_header_indexes_.begin();tokit != int_header_indexes_.end();tokit++){
                ofs_integrals<<values_column_headers_[*tokit]<<"("<<counter<<")"<<"\t";
                counter++;
            }
            ofs_integrals<<endl;
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for  (size_t i = 0 ; i < values_.size();i++){
                for (auto hit = int_header_indexes_.begin();hit != int_header_indexes_.end();hit++)
                    ofs_integrals<<values_[i][*hit]<<"\t";
                ofs_integrals<<endl;
            }
        }

        if (!ranged_property_names_.empty()){
            ofstream  ofs_ranges( file_ranges ,std::ofstream::out );
            size_t counter(1);
            for (auto tokit = range_header_indexes_.begin();tokit != range_header_indexes_.end();tokit++){
                ofs_ranges<<values_column_headers_[*tokit]<<"("<<counter<<")"<<"\t";
                counter++;
            }
            ofs_ranges<<endl;
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for  (size_t i = 0 ; i < values_.size();i++){
                for (auto hit = range_header_indexes_.begin();hit != range_header_indexes_.end();hit++)
                    ofs_ranges<<values_[i][*hit]<<"\t";
                ofs_ranges<<endl;
            }
        }


        if (!single_value_property_names_.empty()){
            ofstream  ofs_single_value( file_single_value ,std::ofstream::out );
            size_t counter(1);
            for (auto tokit = single_value_header_indexes_.begin();tokit != single_value_header_indexes_.end();tokit++){
                ofs_single_value<<values_column_headers_[*tokit]<<"("<<counter<<")"<<"\t";
                counter++;
            }
            ofs_single_value<<endl;
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for  (size_t i = 0 ; i < values_.size();i++){
                for (auto hit = single_value_header_indexes_.begin();hit != single_value_header_indexes_.end();hit++)
                    ofs_single_value<<values_[i][*hit]<<"\t";

                ofs_single_value<<endl;

            }
        }

        if (!dimensional_property_names_.empty()){
            ofstream  ofs_areas_volumes( file_areas_volumes ,std::ofstream::out );
            size_t counter(1);
            for (auto tokit = dimensional_header_indexes_.begin();tokit != dimensional_header_indexes_.end();tokit++){
                ofs_areas_volumes<<values_column_headers_[*tokit]<<"("<<counter<<")"<<"\t";
                counter++;
            }
            ofs_areas_volumes<<endl;
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for  (size_t i = 0 ; i < values_.size();i++){
                for (auto hit = dimensional_header_indexes_.begin();hit != dimensional_header_indexes_.end();hit++)
                    ofs_areas_volumes<<values_[i][*hit]<<"\t";
                ofs_areas_volumes<<endl;
            }
        }

        first_monitor_call_=false;
    }
    else {
        if (!integrated_property_names_.empty()){
            ofstream  ofs_integrals( file_integrals, std::ofstream::out | std::ofstream::app);
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for (auto hit = int_header_indexes_.begin();hit != int_header_indexes_.end();hit++)
                ofs_integrals<<rowdata_[*hit]<<"\t";
            ofs_integrals<<endl;
        }

        if (!ranged_property_names_.empty()){
            ofstream  ofs_ranges( file_ranges ,std::ofstream::out | std::ofstream::app);
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for (auto hit = range_header_indexes_.begin();hit != range_header_indexes_.end();hit++)
                ofs_ranges<<rowdata_[*hit]<<"\t";
            ofs_ranges<<endl;
        }

        if (!single_value_property_names_.empty()){
            ofstream ofs_single_value( file_single_value ,std::ofstream::out | std::ofstream::app);
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for (auto hit = single_value_header_indexes_.begin();hit != single_value_header_indexes_.end();hit++)
                ofs_single_value<<rowdata_[*hit]<<"\t";
            ofs_single_value<<endl;
        }

        if (!dimensional_property_names_.empty()){
            ofstream  ofs_areas_volumes( file_areas_volumes ,std::ofstream::out | std::ofstream::app);
            // now write in any possible read-in data from previous runs. These are values loaded into the database during construction of this class!
            for (auto hit = dimensional_header_indexes_.begin();hit != dimensional_header_indexes_.end();hit++)
                ofs_areas_volumes<<rowdata_[*hit]<<"\t";
            ofs_areas_volumes<<endl;
        }
    }
} // end Out



/// P. Lang: loops over all elements of all regions of model and returns true if IsVolumeElement()
template<size_t dim>
bool  SimulatorMonitor<dim>::HasVolumeElements( const Model<dim>& mref ) const
{
    if( dim == 1U || dim == 2U )
        return false;

    const Region<dim>&  rref( mref.Region("Model") );

    typename vector<Element<dim>*>::const_iterator  elementsEnd = rref.ElementsEnd();
    for ( typename vector<Element<dim>*>::const_iterator eit=rref.ElementsBegin(); eit!=elementsEnd; ++eit )
        if ( (*(*eit)->FE()).IsVolumeElement() )
            return true;

    return false;
}


template<size_t dim>
SimulatorMonitor<dim>::~SimulatorMonitor(){}



template<size_t dim>
bool SimulatorMonitor<dim>::is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

template class SimulatorMonitor<1U>;
template class SimulatorMonitor<2U>;
template class SimulatorMonitor<3U>;


} // end namespace csmp














