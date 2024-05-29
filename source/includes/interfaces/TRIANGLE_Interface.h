#ifndef CSMP_TRIANGLE_INTERFACE_H
#define CSMP_TRIANGLE_INTERFACE_H

#include "CSMP_definitions.h"
#include "VSet.h"

namespace csmp {

/// interface to J. Setchuk's (UC Berkeley's) Triangle mesh generator
class TRIANGLE_Interface {
  public:
    template<uint32_t dim>                        
    void ReadTriangle2DMesh( const char* fname, VSet<dim>&, bool isoparametric=false, bool perform_extra_checks=false );
    
    template<uint32_t dim>                        
    void ReadTriangle2DMeshAndCreateDiscreteFractures( const char* fname, VSet<dim>& vset, double kfrac,
                                                       bool perform_extra_checks=false );
    
    void CheckTriangleOutput( const char*  file, 
                              size_t&   nodes_per_element,
                              size_t&   nodes, 
                              size_t&   elements,
                              int&   node_attributes,
                              int&   elmt_attributes );
    
    void ReadNodeDataFile( const char* file, 
                           std::deque<double>& x, std::deque<double>& y, std::deque<double>& z,
                           std::map<size_t,int8_t>& bflags,
                           std::map<size_t,double>& bvalues );

    void ReadElementDataFile( const char* file, 
                              std::map<size_t,std::vector<size_t> >& plist,
                              std::vector<double>& evalues );

    void ReadPolyDataFile( const char* file, 
                           std::map<size_t,std::vector<size_t> >& plist );

    void ReadNeighborDataFile( const char* file, std::map<size_t,std::vector<int64_t> >& pfverts );
    
    void FindNeighborsForFractureElements( std::map<size_t,std::vector<size_t> >& plist,
                                           std::map<size_t,std::vector<int64_t> >& pfverts );

    void FlagBoundaryElements( std::map<size_t,std::vector<int64_t> >& pfverts,
                               std::map<size_t,std::vector<size_t> >& plist,
                               std::deque<double>& x, std::deque<double>& y );

    void FlagBoundaryNodes( std::map<size_t,std::vector<int64_t> >& pfverts,
                            std::map<size_t,std::vector<size_t> >& plist,
                            std::map<size_t,int8_t>& bflags );
                     
    void SplitSingleCornerElements( std::map<size_t,std::vector<size_t> >& plist,
                                    std::map<size_t,std::vector<int64_t> >& pfverts, 
                                    std::vector<double>& evalues );
  private: 
    int32_t  FindFaceBoundary( std::vector<double>& face_node1,
                               std::vector<double>& face_node2,
                               std::vector<double>& opposite_node,
                               bool verbose=false );

    bool VerifyConsecutiveNodeNumbering( std::map<size_t,std::vector<size_t> >& plist ) const; 
    
    void ListZeroPropertyValueElements( std::vector<double>& evalues,
                                        std::map<size_t,std::vector<size_t> >& plist, 
                                        std::deque<double>& x, 
                                        std::deque<double>& y ) const;

    void FlagBoundaryNodesAccordingTo( unsigned int fvert, std::int8_t bflag,
                                       const std::vector<size_t>& nds, 
                                       std::map<size_t,int8_t>& bflags );
    
    void FlagCornerNodes( std::map<size_t,int8_t>& bflags,
                          std::deque<double>& x, std::deque<double>& y, std::deque<double>& z );  
 };

 }

#endif
