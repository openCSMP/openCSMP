// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "XGRAPH_Interface.h"
#include "Point.h"
#include "Region.h"
#include "Model.h"

using namespace std;

namespace csmp {

// use  xgraph -M   to view in X-term window
void outputToXGraph( const Model<1U>& model, const char* file_name, double model_time, bool append_to_files )
 {
    string  simulation_name(file_name);
    
    // output properties
    csmp::Index  sato_key = model.Database().StorageKey("saturation oil");
    csmp::Index  flup_key = model.Database().StorageKey("fluid pressure");
    csmp::Index  absp_key = model.Database().StorageKey("absolute fluid pressure");
    csmp::Index  flux_key = model.Database().StorageKey("volume flux");
    csmp::Index  mobt_key = model.Database().StorageKey("total mobility");

    ofstream  ofs_so, ofs_pf,  ofs_ap, ofs_lt, ofs_fx;
    
    if ( !append_to_files ) {
         // opening files in overwrite mode
         ofs_so.open((simulation_name+"_saturation-oil.txt").c_str());
         ofs_pf.open((simulation_name+"_fluid-pressure.txt").c_str());
         ofs_ap.open((simulation_name+"_absolute-fluid-pressure.txt").c_str());
         ofs_lt.open((simulation_name+"_total-mobility.txt").c_str()); 
         ofs_fx.open((simulation_name+"_volume-flux.txt").c_str());
         // writing file headers
         ofs_so <<"TitleText: oil saturation (SI units)"<< endl;
         ofs_pf <<"TitleText: fluid pressure (SI units)"<< endl;
         ofs_ap <<"TitleText: absolute fluid pressure (SI units)"<< endl;
         ofs_lt <<"TitleText: total mobility (SI units)"<< endl;
         ofs_fx <<"TitleText: volume flux (SI units)"<< endl;
      }
    else {
         // opening files in append mode
         ofs_so.open((simulation_name+"_saturation-oil.txt").c_str(), ios::out | ios::app );
         ofs_pf.open((simulation_name+"_fluid-pressure.txt").c_str(), ios::out | ios::app);
         ofs_ap.open((simulation_name+"_absolute-fluid-pressure.txt").c_str(), ios::out | ios::app );
         ofs_lt.open((simulation_name+"_total-mobility.txt").c_str(), ios::out | ios::app ); 
         ofs_fx.open((simulation_name+"_volume-flux.txt").c_str(), ios::out | ios::app);
      }

    // writing the dataset titles
    ofs_so <<"\n\"timestep = "<< model_time <<" secs\""<< endl;
    ofs_pf <<"\n\"timestep = "<< model_time <<" secs\""<< endl;
    ofs_ap <<"\n\"timestep = "<< model_time <<" secs\""<< endl;
    ofs_lt <<"\n\"timestep = "<< model_time <<" secs\""<< endl;
    ofs_fx <<"\n\"timestep = "<< model_time <<" secs\""<< endl;
    
    const Region<1>&  sgref(model.Region("Model"));
    
    // organising node output data by x-coordinate
    map<double,vector<double> >  nodeprop_map;
    vector<double>               nodeprops(3U);

    for ( auto it=sgref.NodesBegin(); it!=sgref.NodesEnd(); it++ ) {
          nodeprops[0] = (*it)->Read( sato_key );
          nodeprops[1] = (*it)->Read( flup_key );
          nodeprops[2] = (*it)->Read( absp_key );
          nodeprop_map.insert( make_pair( (*it)->x(), nodeprops ) );
      }
       
    // writing nodal variable sets
    for ( map<double,vector<double> >::const_iterator
          it=nodeprop_map.begin(); it!=nodeprop_map.end(); it++ ) 
      {
         ofs_so << (*it).first <<" "<< (*it).second[0] << endl;
         ofs_pf << (*it).first <<" "<< (*it).second[1] << endl;
         ofs_ap << (*it).first <<" "<< (*it).second[2] << endl;
      }
    nodeprop_map.clear();  
    
      
    // organising element output data by x-coordinate
    map<double,vector<double> >  elmtprop_map;
    vector<double>               elmtprops(2U);

    for ( auto it=sgref.CellsBegin(); it!=sgref.CellsEnd(); it++ ) {
          elmtprops[0] = (*it)->Read( mobt_key );
          elmtprops[1] = (*it)->Read( flux_key );
          Point<1U> p  = (*it)->BaryCenter();
          elmtprop_map.insert( make_pair( p[0U], elmtprops ) );
      }
    // writing element variable sets with barycentric property locations
    for ( auto it=elmtprop_map.begin(); it!=elmtprop_map.end(); it++ )
      {
         ofs_lt << (*it).first <<" "<< (*it).second[0] << endl;
         ofs_fx << (*it).first <<" "<< (*it).second[1] << endl;
      }
    elmtprop_map.clear();
      
   // closing the output files   
   ofs_so.close();
   ofs_pf.close();
   ofs_ap.close();
   ofs_lt.close();
   ofs_fx.close();
   
   cout <<"\noutputToXGraph: file set '"<< simulation_name <<"' written successfully."<< endl;
    
 } // end OutputToXGraph

} // end csmp
