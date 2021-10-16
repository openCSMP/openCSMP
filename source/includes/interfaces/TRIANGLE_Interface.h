#ifndef CSMP_TRIANGLE_INTERFACE_H
#define CSMP_TRIANGLE_INTERFACE_H

#include "CSMP_definitions.h"
#include "VSet.h"

namespace csmp {

/// interface to J. Setchuk's (UC Berkeley's) Triangle mesh generator
class TRIANGLE_Interface {
  public:
    template<size_t dim>                        
    void ReadTriangle2DMesh( const char* fname, VSet<dim>&, bool isoparametric=false, bool perform_extra_checks=false );
    
    template<size_t dim>                        
    void ReadTriangle2DMeshAndCreateDiscreteFractures( const char* fname, VSet<dim>& vset, double64 kfrac,
                                                       bool perform_extra_checks=false );
    
    void CheckTriangleOutput( const char*  file, 
                              size_t&   nodes_per_element,
                              size_t&   nodes, 
                              size_t&   elements,
                              int&   node_attributes,
                              int&   elmt_attributes );
    
    void ReadNodeDataFile( const char* file, 
                           std::deque<double64>& x, std::deque<double64>& y, std::deque<double64>& z,
                           std::vector<std::int8_t>&  bflags, std::map<size_t,double64>& bvalues );

    void ReadElementDataFile( const char* file, 
                              std::map<size_t,std::vector<long64> >& plist,
                              std::vector<double64>& evalues );

    void ReadPolyDataFile( const char* file, 
                           std::map<size_t,std::vector<long64> >& plist );

    void ReadNeighborDataFile( const char* file, std::map<size_t,std::vector<long64> >& pfverts );
    
    void FindNeighborsForFractureElements( std::map<size_t,std::vector<long64> >& plist,
                                           std::map<size_t,std::vector<long64> >& pfverts );

    void FlagBoundaryElements( std::map<size_t,std::vector<long64> >& pfverts,
                               std::map<size_t,std::vector<long64> >& plist,
                               std::deque<double64>& x, std::deque<double64>& y );

    void FlagBoundaryNodes( std::map<size_t,std::vector<long64> >& pfverts,
                            std::map<size_t,std::vector<long64> >& plist,
                            std::vector<std::int8_t>& bflags );
                     
    void SplitSingleCornerElements( std::map<size_t,std::vector<long64> >& plist,
                                    std::map<size_t,std::vector<long64> >& pfverts, 
                                    std::vector<double64>& evalues );
  private: 
    int32  FindFaceBoundary( std::vector<double64>& face_node1,
                                   std::vector<double64>& face_node2,
                                   std::vector<double64>& opposite_node,
                                   bool verbose=false ); 

    bool VerifyConsecutiveNodeNumbering( std::map<size_t,std::vector<long64> >& plist ) const; 
    
    void ListZeroPropertyValueElements( std::vector<double64>& evalues,
                                        std::map<size_t,std::vector<long64> >& plist, 
                                        std::deque<double64>& x, 
                                        std::deque<double64>& y ) const;

    void FlagBoundaryNodesAccordingTo( size_t fvert, long64 bflag,
                                       const std::vector<long64>& nds, 
                                       std::vector<std::int8_t>& bflags );
    
    void FlagCornerNodes( std::vector<std::int8_t>& bflags,
                          std::deque<double64>& x, std::deque<double64>& y, std::deque<double64>& z );  
 };

 }

#endif
