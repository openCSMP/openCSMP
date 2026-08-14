#ifndef NODE_CENTERED_FINITE_VOLUME_MONITOR_H
#define NODE_CENTERED_FINITE_VOLUME_MONITOR_H

#include "Model.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class StencilProcessor;

/// monitoring of variables stored on the nodes, integrating them over the finite volumes
template<uint32_t dim>
class NodeCenteredFiniteVolumeMonitor {
  public:
    NodeCenteredFiniteVolumeMonitor( const char* output_file,
                                     const char* prop, 
                                     bool monitor_regions );
    
    /// triggering monitoring with the option to write to the file specified during class construction
    void MonitorPropertyIntegrals( const Model<dim>&,
                                   bool consider_porosity,
                                   bool normalize_by_initial_integral, 
                                   bool write_output );                                   
  
    /// file output of variable values monitored up to current point of time and without normalization
    void Out( const Model<dim>&, bool normalize_values=false ) const;
                                       
  private:
    void SaveToFile( const Model<dim>&, bool normalize_by_initial_integral ) const;
  
  private:
    std::string  output_file;
    std::string  output_variable;
    bool         group_by_group;
    
    std::list<std::pair<double,double> >              integrals;
    std::list<std::pair<double,std::list<double> > >  group_integrals;
};

} // end namespace csmp

#endif
