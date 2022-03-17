#ifndef CSMP_GOCAD_INTERFACE_H
#define CSMP_GOCAD_INTERFACE_H

#include "CSMP_definitions.h"
#include "GocadHeader.h"
#include "GocadPropertyClassHeader.h"
#include "Box.h"

namespace csmp {

template<uint32_t> class VSet;
template<uint32_t> class Model;

template<uint32_t dim>
class GoCadInterface {
  public:
    GoCadInterface();
    ~GoCadInterface();

    /// input of data from GoCad to CSMP
    void  ReadTetrahedralGocad3DSurface( const char* fname, GocadHeader& header, VSet<dim>& vset );

    void  ReadTetrahedralGocad3DMesh( const char* fname, GocadHeader& header, VSet<dim>& vset );
    
    /// output of data from CSMP to GoCad
    void OutputVariableToTSurface( const Model<dim>&, const char* fname, const char* var, 
                                   long timestep=0 ) const;  

    void OutputVariableToTSurface( const Model<dim>&, const char* fname, const char* group_name, 
                                   const char* var, long timestep=0 ) const;  

    void OutputVariablesToTSurface( const Model<dim>&, const char* fname, long timestep=0 ) const;  

    void OutputVariableToTSolid( const Model<dim>& sg, const char* fname, const char* var, 
                                 long timestep=0 ) const;

    void OutputVariableToTSolid( const Model<dim>&, const char* fname, const char* group_name, 
                                 const char* var, long timestep=0 ) const;  

    void   OutputRegionsToGoCadFiles( const Model<dim>&, const char* property, long timestep=0 ) const;

  private:    
    bool IsInNextLine( std::ifstream& ifs, const char* search_string ) const;
    
    BOX_BOUNDARY  IdentifyTetrahedronBoundary( const std::vector<double>& nd1,
                                               const std::vector<double>& nd2,
                                               const std::vector<double>& nd3 );
    bool ReadTSurface( std::ifstream& ifs, 
                       const std::list<GocadPropertyClassHeader>& properties,
                       VSet<dim>& vset );

    void ReadTSolid( std::ifstream& ifs, VSet<dim>& vset, bool with_node_prop=false );
    
    void FlagEdgeNodesOfBoxShapedModel( VSet<dim>& vset );
 
    bool VerifyConsecutiveNodeNumbering( std::map<size_t,std::vector<int64_t> >& plist ) const; 
};

} // csmp

#endif

