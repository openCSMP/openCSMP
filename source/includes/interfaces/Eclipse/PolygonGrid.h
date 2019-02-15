#ifndef POLYGON_GRID_UOM_H
#define POLYGON_GRID_UOM_H

#include "Point.h"
#include "FiniteElement.h"

namespace csmp {

namespace eclipse {

/**
    refactored - strictly 3D
    
    @todo check whether this class is needed at all, 
    if so, remove all un-used interfaces.
*/
class GridNode : public csmp::Point<3U> {
  public:
    GridNode() : idx_(UINT_MAX) {}
    GridNode( size_t index, const Point<3U>& );
    GridNode( const GridNode& );
    // add move constructor
    GridNode& operator=( const GridNode& );

    bool operator==( const GridNode& );

    size_t GetIdx() const;
    bool IsNullIdx() const;
    const csmp::Point<3U>& GetPoint() const;
    csmp::Point<3U>& GetPoint();

    void AssignNullIdx();
    void AssignIdx( size_t );
    void AssignPoint( const csmp::Point<3U>& );

private:
    size_t idx_;
};

std::ostream&  operator<<( std::ostream&, const GridNode& );




class GridFace {
  public:
    GridFace();
    GridFace( const GridFace& );
    GridFace& operator=( const GridFace& );
    ~GridFace();

    bool operator==( const GridFace& );

    // TODO: perhaps not wise to have this dependence on FiniteElement in this class
    const csmp::CSMP_FEM_TYPE& GetType() const;
    size_t GetNumNodes() const;
    GridNode* GetNode( size_t nid );
    GridNode* const GetNode( size_t nid ) const;

    void AssignIdx( size_t );
    void AssignType( const csmp::CSMP_FEM_TYPE& );
    void AddNode( GridNode* );

  private:
    /// face data
    size_t idx_;
    csmp::CSMP_FEM_TYPE     type_;
    std::vector<GridNode*>  nodes_;
};



class GridElement {
  public:
    GridElement();
    GridElement( const GridElement& );
    GridElement& operator=( const GridElement& );
    ~GridElement();

    bool operator==( const GridElement& );

    const csmp::CSMP_FEM_TYPE& GetType() const;
    size_t GetNumNodes() const;
    GridNode* GetNode( size_t nid );
    GridNode* const GetNode( size_t nid ) const;

    void AssignIdx( size_t );
    void AssignType( const csmp::CSMP_FEM_TYPE& );
    void AddNode( GridNode* );

private:
    /// element data
    size_t idx_;
    csmp::CSMP_FEM_TYPE     type_;
    std::vector<GridNode*>  nodes_;
};





/// @attention does not seem to be used anywhere! - remove?
class PolygonGrid {
  public:

    PolygonGrid();
    PolygonGrid( const PolygonGrid& );
    PolygonGrid& operator=( const PolygonGrid& );
    ~PolygonGrid();

 //   const CSMP_ElementSpecifications& GetFemSpecs() const;

    GridNode*    AddNode( GridNode& );
    GridFace*    AddFace( GridFace& );
    GridElement* AddElement( GridElement& );

  private:

    std::deque<GridNode*>     nodes_;
    std::deque<GridFace*>     faces_;
    std::deque<GridElement*>  elements_;
};




// TODO: as it stands, this is not a manager and it does not even contain a PolygonGrid ! - refactor
// TODO: if this is not used, may be it can be eliminated
/**

@class PolygonGrid  PolygonGrid "PolygonGrid.h"
@author R. Manasipov
@date 2015

The PolygonGrid contains Points and GridNodes.

*/
class PolygonGridManager {
  public:
    PolygonGridManager();
    PolygonGridManager( const PolygonGridManager& );
    PolygonGridManager& operator=( const PolygonGridManager& );
    ~PolygonGridManager();
  
    void Clear();

    /// nodes
    size_t GetNumNodes() const;
    // NOT USED: size_t GetNodeId( const csmp::Point<3U>& ) const;
    const csmp::Point<3U>& GetPoint( size_t nid ) const;
    GridNode* GetNode( size_t nid );
    GridNode* GetNode( const csmp::Point<3U>& );
    GridNode* AddNode( const csmp::Point<3U>& );
  
    /// prints object state to screen
    void Out() const;

  private:
    std::map<csmp::Point<3U>,size_t>  points_; // duplicates storage
    std::vector<GridNode*>            grid_nodes_;
};

} // eclipse

}// end namespace csmp

#endif

