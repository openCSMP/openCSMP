#ifndef CSMP_VSET_H
#define CSMP_VSET_H

#include "VData.h"
#include "PropertyData.h"

namespace csmp {

/**

@brief Container for polygonal data and associated properties
implemented on the basis of NCSA (Urbana, Champagne, Illinois, USA) VSet 
that forms part of their HDF data file storage.

@author Stephan K. Matthai
@author Stephen G. Roberts
@date 1996

*/
template<uint32_t dim>
class VSet : public VData {

  public:

    VSet();
    VSet( const VSet& );
    VSet( VSet&& ) = default;
    
    /// single element type constructor
    VSet( uint32_t nodes_per_element,
          uint32_t nbors_per_element,
          int8_t   csmp_etype,
          size_t nodes, size_t elmts );

    /// multiple element type constructor
    VSet( const std::vector<int8_t>&  fem_types,
          const std::deque<uint32_t>& npes,
          const std::deque<uint32_t>& epes,
          size_t nodes );

    virtual ~VSet();
    
    VSet& operator=( const VSet& );
    VSet& operator=( VSet&& ) = default;

    /// resizing meshes with only a single element type
    void Resize( uint32_t nodes_per_element,
                 uint32_t nbors_per_element,
                 int8_t csmp_etype, 
                 size_t nodes, size_t elmts );

    /// use for meshes with multiple element sizes
    void Resize( const std::deque<int8_t>& etypes,
                 const std::deque<uint32_t>& npes, 
                 const std::deque<uint32_t>& epes,
                 size_t nodes, size_t faces, size_t interfaces );

    /// node coordinates
    void AddXYZ( const std::deque<double>& x,
                 const std::deque<double>& y,
                 const std::deque<double>& z );
      
    /// the IDs of the nodes that make up each element 
    void AddPlist( typename std::map<size_t,std::vector<int64_t> >::const_iterator first,
                   typename std::map<size_t,std::vector<int64_t> >::const_iterator last );

    void AddPlist( typename std::deque<std::vector<int64_t> >::const_iterator first,
                   typename std::deque<std::vector<int64_t> >::const_iterator last );

    /// the equi-dimensional neighbors adjacent to the numbered element faces plus boundary identifiers where there is no neighbor
    void AddPfverts( typename std::map<size_t,std::vector<int64_t> >::const_iterator first,
                     typename std::map<size_t,std::vector<int64_t> >::const_iterator last );

    void AddPfverts( typename std::deque<std::vector<int64_t> >::const_iterator first,
                     typename std::deque<std::vector<int64_t> >::const_iterator last );
  
    /// adds BOX_BOUNDARY flag values to VData
    void AddBFlags( typename std::vector<std::int8_t>::const_iterator first,
                    typename std::vector<std::int8_t>::const_iterator last );

    /// rocktype identifiers for elements only
    void AddPmtrl( typename std::vector<int32_t>::const_iterator first,
                   typename std::vector<int32_t>::const_iterator last );
  
    /// checks whether the VSet contains any distributed variable values stored in PropertyData objects
    bool  DataEmpty() const;

    /// adds variable dataset to 'property_map_' data member of VSet
    bool AddData( const char* s, const PropertyData& );
  
    /// retrieves property data from VSet (if any)
    PropertyData  Data( const char* s ) const;
  
    /// if the property records contain FV data
    bool ContainsFiniteVolumeIntegrationPointData() const;
  
    /// const iterators for the material ID record (one for each element; none for face and interface objects)
    std::vector<int32_t>::const_iterator PmtrlBegin() const;
    std::vector<int32_t>::const_iterator PmtrlEnd() const;
    
    int32_t Pmtrl( size_t elmt ) const { return pmtrl_.at(elmt); }

    /// const iterators for the propery collection
    std::map<std::string,PropertyData>::const_iterator PropertyValuesBegin() const;
    std::map<std::string,PropertyData>::const_iterator PropertyValuesEnd() const;

    void RemoveData( const char* s );

    /// writes complete VSet to binary file with the given time stamp
    bool  OutputTo( const char* bin_file, double time ) const;
  
    /// reads binary files written with OutputTo() and initialises the VSet with it; the time stamp is returned in second argument
    bool  InputFrom( const char* bin_file, double& time );

    /// reads binary files written with OutputTo() and initialises the VSet with it; the time stamp is returned in second argument; it can read only a subset of variables by using third argument
    bool  InputFrom( const char* bin_file, double& time, const std::set<std::string>& subset_variables );
  
    /// permits to initialise a VSet from a text file (see detailed DOxygen documentation for format)
    bool  InputFromTextFile( const char* file_dot_txt );
  
    /// output VSet to binary file which has partitions and halos defining their overlap
    bool  ParallelOutputTo( const char* bin_file, double time, size_t first_outerhalo ) const;
    bool  ParallelInputFrom( const char* bin_file, double& time, size_t& first_outerhalo );

    /// select elements specific elements from the VSet that shall be retained while all others are deleted (including nodes)
    void  ReduceTo( const std::map<size_t,size_t>& old_and_new_consecutive_element_ids );
    
    /// updates pmtrl and property storage to size changes in VData
    void UpdatePropertyStorage();
  
    /// deletes all content of the VSet
    void  Erase();
    
    /// uses the node coordinates to infer the model dimension: if Z-range=zero, dim=2, if Y-range=2, dim=1, else dim=3
    int32_t MeshDimension( bool check_coordinates = false ) const;
  
    /// prints the VSet to the console
    void  Out( bool data_as_well=true ) const;
    
    /// writes C++17  code that reproduces a hardwired version of the current VSet
    void OutCPP17( const char* cpp_file ) const;
 

  protected:

    std::map<std::string,PropertyData>  property_map_;  ///< container for variables distributed on the mesh
    std::vector<int32_t>                pmtrl_;         ///<  rocktype identifiers; one per element; for transfer to 'material_id_'

    friend class VSet_Test;
};

/// extrapolates element property values stored in VSet to its vertices
template<uint32_t dim, class VarType>
void extrapolateElementToNodeProperty( VSet<dim>&,
                                       const std::vector<VarType>&  elmnt_values,
                                       std::vector<VarType>&        nodal_values );


} // csmp

#endif

