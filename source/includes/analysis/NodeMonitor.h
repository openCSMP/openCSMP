#ifndef CSMP_NODE_MONITOR_H
#define CSMP_NODE_MONITOR_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t dim> class Model;

/**
    @brief Monitors node / finite volume properties over the course of a run
    in a cluster of nodes of interest.
    
    @author Stephan Matthai
    @date 2001
*/
template<size_t dim>
class NodeMonitor {
  public:
    /// defines nodes and scalar variables stored there that shall be monitored
    NodeMonitor( const Model<dim>&, 
                 const std::map<std::string,size_t>& monitored_nodes,
                 const std::list<std::string>& properties_to_monitor );
                  
    ~NodeMonitor();
  
    /// prompts monitoring of the target property values at the current time
    void  ScalarPropertyValues( const Model<dim>&, double current_time );
  
    /// zap all recorded values and property names
    void  Reset();
  
    /// zap all recorded values
    void  EraseData(); 
  
    /// writes results to text file
    void  Out( const char* text_file ) const;
  
  private:
    // geometric group properties
    //         nodename, idx
    std::map<std::string,size_t>  nodes;
    // properties which are monitored
    std::list<std::string>        properties;
    //     model_time,            property,            nodename, property value
    std::map<double,std::map<std::string,std::map<std::string,double> > >  values;
};

} // csmp

#endif


