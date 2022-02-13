#ifndef CSMP_TECPLOT_INTERFACE_H
#define CSMP_TECPLOT_INTERFACE_H

#include "VSet.h"

namespace csmp {

template<uint32_t>  class Model;
template<uint32_t>  class Region;
    
template<uint32_t dim>
class Tecplot_Interface {
  public:
    void OutputDataToTecplotFile( const Model<dim>&,
                                  const char* file_name, 
                                  const char* var_name, 
                                  long timestep );

    void OutputDataToTecplotFile( const Model<dim>&,
                                  const char* group_name, 
                                  const char* file_name, 
                                  const char* var_name, 
                                  long timestep );
    
    void OutputGridToTecplotFile( const Model<dim>&,
                                  const char* file_name,
                                  long timestep );
    
    void OutputGridToTecplotFile( const Model<dim>&,
                                  const char* group_name,
                                  const char* file_name,
                                  long timestep );
    
  protected:
    
    void OutputDataToTecplotFile( const Model<dim>&,
                                  const Region<dim>&,
                                  const char* file_name,
                                  const char* var_name,
                                  long timestep );
    
    void OutputGridToTecplotFile( const Model<dim>&,
                                  const Region<dim>&,
                                  const char* file_name,
                                  long timestep );

  private:
    VSet<dim>  vset;
    std::map<size_t,std::vector<size_t> >  plist, transformed_plist;
    std::map<size_t,size_t>                node_mapping;
    std::map<size_t,std::vector<double> >  pxyz_data;

    void NodeBasedTopology( const Model<dim>&,
                            const std::vector<size_t>& elmt_ids,
                            std::map<size_t,std::vector<size_t> >&   plist,
                            std::map<size_t,size_t>&                 node_nums );
                            
    void PointBasedTopology( const Model<dim>&,
                            const std::vector<size_t>& elmt_ids,
                            std::map<size_t,std::vector<size_t> >&   plist,
                            std::map<size_t,size_t>&                 node_nums );

    void TransformPlist( const Model<dim>&,
                         const std::map<size_t,std::vector<size_t> >& plist,
                         std::map<size_t,std::vector<size_t> >& tplist );

    void XyzData( const Model<dim>&,
                  const std::map<size_t,size_t>& node_nums,
                  std::map<size_t,std::vector<double> >& pxyz_data );
    
    void NodeData( const Model<dim>&,
                   const Index&     prop_key,
                   const std::map<size_t,size_t>& node_nums,
                   std::map<size_t,std::vector<double> >& pxyz_data );

    void ElementPointData( const Model<dim>&,
                           const Index&     prop_key,
                           const std::map<size_t,size_t>& node_nums,
                           std::map<size_t,std::vector<double> >& pxyz_data );

 };

}

#endif
