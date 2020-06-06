#ifndef CSMP_VSET_H
#define CSMP_VSET_H

#include "VData.h"
#include "PropertyData.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

/**

@brief Container for polygonal data and associated properties
implemented on the basis of NCSA (Urbana, Champagne, Illinois, USA) VSet 
that forms part of their HDF data file storage.

@author Stephan K. Matthai
@author Stephen G. Roberts
@date 1996

*/
template<size_t dim>
class VSet : public VData {

  public:

    VSet();
    VSet( const VSet& );
    VSet( VSet&& ) = default;
    VSet( size_t nodes_per_element, 
          size_t nbors_per_element,
          int32  csmp_etype, 
          size_t nodes, size_t elmts );

    // setting sizes without transfer of data
    VSet( const std::deque<size_t>& npes, 
          const std::deque<size_t>& epes,
          size_t nodes );

    ~VSet();
    VSet& operator=( const VSet& );
    VSet& operator=( VSet&& ) = default;

    void Resize( size_t nodes_per_element, 
                 size_t nbors_per_element,
                 int32  csmp_etype, 
                 size_t nodes, size_t elmts );

    void Resize( const std::deque<int32>& etypes,
                 const std::deque<size_t>& npes, 
                 const std::deque<size_t>& epes,
                 size_t nodes, size_t faces, size_t interfaces );

    void AddXYZ( const std::deque<double64>& x,
                 const std::deque<double64>& y,
                 const std::deque<double64>& z );
      
    void AddPlist( typename std::map<size_t,std::vector<size_t> >::const_iterator first,
                   typename std::map<size_t,std::vector<size_t> >::const_iterator last );

    void AddPlist( typename std::deque<std::vector<size_t> >::const_iterator first,
                   typename std::deque<std::vector<size_t> >::const_iterator last );

    void AddPfverts( typename std::map<size_t,std::vector<long64> >::const_iterator first,
                     typename std::map<size_t,std::vector<long64> >::const_iterator last );

    void AddPfverts( typename std::deque<std::vector<long64> >::const_iterator first,
                     typename std::deque<std::vector<long64> >::const_iterator last );
  
  
    /// checks whether the VSet contains any distributed variable values stored in PropertyData objects
    bool  DataEmpty() const;

    /// adds variable dataset to 'property_map_' data member of VSet
    bool AddData( const char* s, const PropertyData& );
  
    /// retrieves property data from VSet (if any)
    PropertyData  Data( const char* s ) const;
  
    /// if the property records contain FV data
    bool ContainsFiniteVolumeIntegrationPointData() const;
  
    /// const iterators for the propery collection
    std::map<std::string,PropertyData>::const_iterator PropertyValuesBegin() const;;
    std::map<std::string,PropertyData>::const_iterator PropertyValuesEnd() const;;

    void RemoveData( const char* s );

    void AddBFlags( typename std::unordered_map<size_t,long64>::const_iterator first,
                    typename std::unordered_map<size_t,long64>::const_iterator last );

    /// writes complete VSet to binary file with the given time stamp
    bool  OutputTo( const char* bin_file, double64 time ) const;
  
    /// reads binary files written with OutputTo() and initialises the VSet with it; the time stamp is returned in second argument; it can read only a subset of variables by using third argument
    bool  InputFrom( const char* bin_file, double64& time, const std::set<std::string>* subset_variables = nullptr );
  
    /// permits to initialise a VSet from a text file (see detailed DOxygen documentation for format)
    bool  InputFromTextFile( const char* file_dot_txt );
  
    /// output VSet to binary file which has partitions and halos defining their overlap
    bool  ParallelOutputTo( const char* bin_file, double64 time, size_t first_outerhalo ) const;
    bool  ParallelInputFrom( const char* bin_file, double64& time, size_t& first_outerhalo );

    /// select elements specific elements from the VSet that shall be retained while all others are deleted (including nodes)
    void  ReduceTo( const std::map<size_t,size_t>& old_and_new_consecutive_element_ids );
  
    /// deletes all content of the VSet
    void  Erase();
    
    /// uses the node coordinates to infer the model dimension: if Z-range=zero, dim=2, if Y-range=2, dim=1, else dim=3
    int32 MeshDimension( bool check_coordinates = false ) const;
  
    /// prints the VSet to the console
    void  Out( bool data_as_well=true ) const;

  protected:

    std::map<std::string,PropertyData>  property_map_;  ///< container for variables distributed on the mesh

    friend class VSet_Test;
};

} // csmp

#endif

