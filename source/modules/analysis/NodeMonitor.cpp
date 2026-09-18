// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NodeMonitor.h"
#include "Node.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
NodeMonitor<dim>::NodeMonitor( const Model<dim>& sg, 
                                  const std::map<std::string,size_t>& monitored_nodes,
                                  const std::list<std::string>& properties_to_monitor )
                 
 : nodes(monitored_nodes),
   properties(properties_to_monitor)
 {
    // checking properties
    for ( typename list<string>::const_iterator
          it=properties.begin(); it!=properties.end(); it++ )
      if ( !sg.Database().IsDefined( (*it).c_str() ) or 
            sg.Database().Placement( (*it).c_str() ) != NODE ) 
        throw csmp::Exception( FATAL_ERROR, "NodeMonitor(constructor)", (*it).c_str(),
                       "range property is not a node property or is not defined in the database (file)");
                       
 } // end constructor
 
 
 
 
template<uint32_t dim>
NodeMonitor<dim>::~NodeMonitor()
 {
 }

    
    
template<uint32_t dim>
void NodeMonitor<dim>::ScalarPropertyValues( const Model<dim>& sg, double time )
 {
    list<string>::const_iterator        lit;
    map<string,size_t>::const_iterator  nit;
    double                            value;
    const Region<dim>&          sgroup(sg.Region("Model"));
    
    // for all properties for which values shall be monitored
    for ( lit=properties.begin(); lit!=properties.end(); lit++ ) {
        csmp::Index prop_key(sg.Database().StorageKey( (*lit).c_str() )); 
        for ( nit=nodes.begin(); nit!=nodes.end(); nit++ ) {
             // reading scalar property values
             value = sgroup.N( (*nit).second )->Read( prop_key ); 
             values[time][(*lit)][(*nit).first] = value;
          }
     }
      
 } // end ScalarPropertyRanges
    

 
 
 
template<uint32_t dim>
void NodeMonitor<dim>::Reset() // zap all recorded values
 {
    // properties which are monitored
    properties.erase( properties.begin(), properties.end() );
    //    model_time,       property,      groupname, property ranges
    values.erase( values.begin(), values.end() );

 } // end Reset



template<uint32_t dim>
void NodeMonitor<dim>::EraseData() // zap all recorded values
 {
    values.erase( values.begin(), values.end() );

 } // end EraseData
 


    
    
template<uint32_t dim>
void NodeMonitor<dim>::Out( const char* text_file ) const
 {
    string  file(text_file);
    file +=".txt";
    
    ofstream  ofs( file.c_str() ); 

    ofs <<"'"<< file <<"' node data recorded up to time, t = "<< (*values.rbegin()).first <<" seconds.";

    // writing property data blocks
    // ----------------------------
    // for all property ranges determined for the groups
    ofs <<"\nValue ranges of properties: "<< endl;
    for ( list<string>::const_iterator iit=properties.begin(); iit!=properties.end(); iit++ ) 
      {
         // name of monitored property
         ofs << endl << (*iit) << endl;

         // the column headers: time, node1, node2...
         ofs <<" time\t\t";
         for ( map<string,size_t>::const_iterator
               nit=nodes.begin(); nit!=nodes.end(); nit++ )
           ofs << (*nit).first <<"\t\t";
         ofs << endl;

         // for all timesteps
         // model_time       property      groupname property ranges
         for ( typename map<double,map<std::string,map<std::string,double> > >::const_iterator
               it=values.begin(); it!=values.end(); it++ ) {
              // the time
              ofs << (*it).first <<"\t";
              // the property data-record searched in the map
              typename map<string,map<std::string,double> >::const_iterator prop_it=(*it).second.find(*iit);
              if ( prop_it != (*it).second.end() )
                // for all property values recorded for the groups
                for ( typename map<std::string,double>::const_iterator
                      pit=(*prop_it).second.begin(); pit!=(*prop_it).second.end(); pit++ )
                  // writing values to file
                  ofs << (*pit).second <<"\t";
                ofs << endl;
           }
         ofs << endl;
      }

 } // end Out


template class NodeMonitor<1U>;
template class NodeMonitor<2U>;
template class NodeMonitor<3U>;


} // end namespace csmp













    
