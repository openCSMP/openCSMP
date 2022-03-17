#ifndef CSMP_FACE_DATA_H
#define CSMP_FACE_DATA_H

#include <set>
#include <vector>
#include "FiniteElement.h"

namespace csmp {

class FaceData {
  public:
    FaceData() {};
    ~FaceData() {};
    FaceData( const FaceData& );
    FaceData& operator=( const FaceData& );
    bool      operator<( const FaceData& ) const;
    size_t    operator[]( size_t ) const; // gives ordered node IDs
    
    void      AssignKey( const std:: vector<uint32_t>& );
    void      AssignNodeIDs( const std::vector<uint32_t>& );
    void      AssignFirstNeighbor( int32_t );
    void      AssignSecondNeighbor( int32_t );
    void      AssignFirstNeighborFaceNumber( size_t );
    void      AssignSecondNeighborFaceNumber( size_t );
    
    int32_t     FirstNeighbor() const;
    int32_t     SecondNeighbor() const;
    size_t    FirstNeighborFaceNumber() const;
    size_t    SecondNeighborFaceNumber( ) const;

    void           SetFiniteElementType( const CSMP_FEM_TYPE& fetype );
    CSMP_FEM_TYPE  FiniteElementType() const;
    
    void      Reset();
    
  private:
    std::set<uint32_t>          key;      // ordered face-node IDs
    std::vector<uint32_t>       nodes;    // node IDs in CCW order
    std::pair<size_t,size_t>  efnumber; // face-number of connected elements
    std::pair<int32_t,int32_t>    nbors;    // ID numbers of connected elements
    CSMP_FEM_TYPE             etype;    // type of finite-element for the face
};





inline FaceData::FaceData( const FaceData& fd )
  : key(fd.key),      
    nodes(fd.nodes),    
    efnumber(fd.efnumber), 
    nbors(fd.nbors),    
    etype(fd.etype)    
  {
  }

inline bool FaceData::operator<( const FaceData& fd ) const
  {
     return key < fd.key;
  }

// gives ordered node IDs
inline size_t FaceData::operator[]( size_t i ) const
 {
    return nodes[i];
 }


inline void FaceData::AssignKey( const std::vector<uint32_t>& nids )
 {
    for ( std::vector<uint32_t>::const_iterator it=nids.begin();
          it!=nids.end(); it++ )
      key.insert( *it );  
 }
 
 
inline void FaceData::AssignNodeIDs( const std::vector<uint32_t>& nids )
 {
    nodes.assign( nids.begin(), nids.end() );  
 }

inline void FaceData::AssignFirstNeighborFaceNumber( size_t e1face_n )
 {
    efnumber.first = e1face_n;
 }
 
 
inline void FaceData::AssignSecondNeighborFaceNumber( size_t e2face_n )
 {
    efnumber.second = e2face_n;
 }



inline void FaceData::AssignFirstNeighbor( int32_t nbor1 )
 {
    nbors.first = nbor1;
 }
 
 
inline void FaceData::AssignSecondNeighbor( int32_t nbor2 )
 {
    nbors.second = nbor2;
 }


inline void FaceData::SetFiniteElementType( const CSMP_FEM_TYPE& fetype )
 {
    etype = fetype;
 }
 
 
inline CSMP_FEM_TYPE  FaceData::FiniteElementType() const
 {
    return etype;
 }


inline int32_t   FaceData::FirstNeighbor() const { return nbors.first; }
inline int32_t   FaceData::SecondNeighbor() const { return nbors.second; }
inline size_t  FaceData::FirstNeighborFaceNumber() const { return efnumber.first; }
inline size_t  FaceData::SecondNeighborFaceNumber( ) const { return efnumber.second; }



} // end namespace csp


#endif


