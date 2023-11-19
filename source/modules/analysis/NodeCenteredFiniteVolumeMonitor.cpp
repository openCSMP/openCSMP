#include "NodeCenteredFiniteVolumeMonitor.h"
#include "Region.h"
#include "Model.h"
#include "ErrorHandler.h"
#include "StencilProcessor.h"
#include "ModelTime.h"

using namespace std;

namespace csmp {

/// text-file based monitoring of user-specified target variable on either, entire model or its unique regions
template<uint32_t dim>
NodeCenteredFiniteVolumeMonitor<dim>::NodeCenteredFiniteVolumeMonitor( const char* out_file, 
                                                                       const char* prop,
                                                                       bool monitor_groups )
 : output_file(out_file),
   output_variable(prop),
   group_by_group(monitor_groups)
 {
    output_file += "-";
    output_file += prop;
    output_file += ".txt";
    for ( string::iterator it=output_file.begin(); it!=output_file.end(); it++ )
      if ( *it == ' ' ) *it = '_';
    
    cout <<"\nNodeCenteredFiniteVolumeMonitor<"<< dim;
    cout <<">(constructor): Results will be written to text file '";
    cout << output_file <<"'"<< endl;
 }




/**
    Integrating current values of the target variables over the entire model domain.
    If region by region approach is chosen, only the unique regions are considered.
*/
template<uint32_t dim>
void NodeCenteredFiniteVolumeMonitor<dim>::MonitorPropertyIntegrals( const Model<dim>& sg, 
                                                                     bool consider_porosity,
                                                                     bool normalize_by_initial_integral,
                                                                     bool write_output )
 {
    double& model_time( ModelTime::Instance().modelTime );
    
    // 0. integrating the property
    const Region<dim>&  domain(sg.Region("Model"));
    double integratedPropertyValue = domain.VolumeIntegral_x_Thickness( output_variable.c_str(), consider_porosity );
    assert( !isnan( integratedPropertyValue) );
    integrals.push_back( make_pair( model_time, integratedPropertyValue ) );

    // 1. simple case: property is integrated and output for the entire model
    if ( write_output ) SaveToFile( sg, normalize_by_initial_integral );

 } // end MonitorNodePropertyIntegrals



/// updates the values one more time and then writes output to file specified in constructor
template<uint32_t dim>
void NodeCenteredFiniteVolumeMonitor<dim>::Out( const Model<dim>& model, bool normalize_values ) const
 {
    SaveToFile( model, normalize_values );
   
 } // end Out




template<uint32_t dim>
void NodeCenteredFiniteVolumeMonitor<dim>::SaveToFile( const Model<dim>& sg,
                                                       bool normalize_by_initial_integral ) const
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
     if ( integrals.empty() && group_integrals.empty() ) {
          csmp_error.Note( WARNING, "NodeCenteredFiniteVolumeMonitor::SaveToFile", "No integral data yet; no output.");
          return;
       }
 
     ofstream  ofs( output_file.c_str() ); 

     ofs <<"NodeCenteredFiniteVolumeMonitor<"<< dim;
     ofs <<">::SaveToFile: '"<< output_file <<"' variable '"<< output_variable;
     ofs <<"' recorded up to time, t = "<< (*integrals.rbegin()).first <<" seconds.";
     
     if ( !group_by_group ) {
          double  initial(1.); 
          if ( normalize_by_initial_integral ) initial = (*integrals.begin()).second;
          ofs <<"\ntime (s)\t"<< output_variable << endl;
          for ( typename list<pair<double,double> >::const_iterator 
                it=integrals.begin(); it!=integrals.end(); it++ )
            ofs << (*it).first <<"\t"<< (*it).second / initial << endl;
          cout <<"\nNodeCenteredFiniteVolumeMonitor<"<<  dim;
          cout <<">::SaveToFile: '"<< output_file <<"' written successfully."<< endl;
          return;
       }
       
     if ( group_integrals.empty() ) {
          csmp_error.Note( WARNING, "NodeCenteredFiniteVolumeMonitor::SaveToFile", "No integral data for regions yet; no output.");
          return;
       }

     ofs <<"\ntime (s)"; // and the names of the groups delimited by tabs
     for ( auto it=sg.UniqueRegionsBegin(); it!=sg.UniqueRegionsEnd(); it++ )
       ofs <<"\t"<< (*it).first;
     ofs << endl;
     
     // writing out the regional integrals for each monitored timestep
     for ( auto it=group_integrals.begin(); it!=group_integrals.end(); it++ ) {
          // writing the time
          ofs << (*it).first;
          // tracking the initial values if necessary
          typename list<double>::const_iterator  iit=(*group_integrals.begin()).second.begin();
          // writing the integrals over each group
          ofs << scientific << setprecision(numeric_limits<double>::digits10);
          for ( auto lit=(*it).second.begin(); lit!=(*it).second.end(); lit++, iit++ )
            // normalizing the integral value if so requested
            if ( normalize_by_initial_integral ) ofs <<"\t"<< (*lit) / (*iit);
            else                                 ofs <<"\t"<< (*lit);
       }
     ofs << endl;

     cout <<"\nNodeCenteredFiniteVolumeMonitor<"<<  dim;
     cout <<">::SaveToFile: '"<< output_file <<"' written successfully."<< endl;

 } // end SaveToFile

template class NodeCenteredFiniteVolumeMonitor<1U>;
template class NodeCenteredFiniteVolumeMonitor<2U>;
template class NodeCenteredFiniteVolumeMonitor<3U>;

} // end namespace csmp
