#ifndef POLYGON_GRID_H
#define POLYGON_GRID_H

#include "CSMP_definitions.h"
#include "CSMP_ElementSpecifications.h"
#include "Point.h"

namespace csmp {

template<size_t dim>
class GridNode
{
public:
    GridNode(); // TODO: make private
    ~GridNode();
    GridNode( size_t index, const Point<dim>& );
    GridNode( const GridNode& );
    GridNode& operator=( const GridNode& );

    bool operator==( const GridNode& );

    size_t GetIdx() const;
    bool IsNullIdx() const;
    const csmp::Point<dim>& GetPoint() const;
    csmp::Point<dim>& GetPoint();

    void AssignNullIdx();
    void AssignIdx( size_t );
    void AssignPoint( const csmp::Point<dim>& pt );

private:
    /// node data
    size_t idx_;
    csmp::Point<dim> pt_;
};

template<size_t dim>
std::ostream&  operator<<( std::ostream&, const GridNode<dim>& );




template<size_t dim>
class GridFace
{
public:

    GridFace();
    GridFace( const GridFace& );
    GridFace& operator=( const GridFace& );
    ~GridFace();

    bool operator==( const GridFace& );

    const csmp::CSMP_FEM_TYPE& GetType() const;
    size_t GetNumNodes() const;
    GridNode<dim>* GetNode( size_t nid );
    GridNode<dim>* const GetNode( size_t nid ) const;

    void AssignIdx( size_t );
    void AssignType( const csmp::CSMP_FEM_TYPE& type );
    void AddNode( GridNode<dim>* gn );

private:
    /// face data
    size_t idx_;
    csmp::CSMP_FEM_TYPE           type_;
    std::vector<GridNode<dim>*>   nodes_;
};



template<size_t dim>
class GridElement
{
public:

    GridElement();
    GridElement( const GridElement& );
    GridElement& operator=( const GridElement& );
    ~GridElement();

    bool operator==( const GridElement& );

    const csmp::CSMP_FEM_TYPE& GetType() const;
    size_t GetNumNodes() const;
    GridNode<dim>* GetNode( size_t nid );
    GridNode<dim>* const GetNode( size_t nid ) const;

    void AssignIdx( size_t );
    void AssignType( const csmp::CSMP_FEM_TYPE& type );
    void AddNode( GridNode<dim>* pt );

private:
    /// element data
    size_t idx_;
    csmp::CSMP_FEM_TYPE         type_;
    std::vector<GridNode<dim>*> nodes_;
};


/// @attention does not seem to be used anywhere! - remove?
template<size_t dim>
class PolygonGrid
{
  public:

    PolygonGrid();
    PolygonGrid( const PolygonGrid& );
    PolygonGrid& operator=( const PolygonGrid& );
    ~PolygonGrid();

    const CSMP_ElementSpecifications& GetFemSpecs() const;

    GridNode<dim>*    AddNode( csmp::GridNode<dim>& gn );
    GridFace<dim>*    AddFace( csmp::GridFace<dim>& gf );
    GridElement<dim>* AddElement( csmp::GridElement<dim>& ge );

  private:

    std::deque<GridNode<dim>*>     nodes_;
    std::deque<GridFace<dim>*>     faces_;
    std::deque<GridElement<dim>*>  elements_;
    CSMP_ElementSpecifications     fem_specs_;
};



// TODO: as it stands, this is not a manager
template<size_t dim>
class PolygonGridManager
{
  public:

    PolygonGridManager();
    PolygonGridManager( const PolygonGridManager& );
    PolygonGridManager& operator=( const PolygonGridManager& );
    ~PolygonGridManager();
    void Clear();

    /// fem specs
    const CSMP_ElementSpecifications& GetFemSpecs() const;

    /// nodes
    size_t GetNumNodes() const;
    size_t GetNodeId( const csmp::Point<dim>& pt ) const;
    const csmp::Point<dim>& GetPoint( size_t nid ) const;
    GridNode<dim>* GetNode( size_t nid );
    GridNode<dim>* GetNode( const csmp::Point<dim>& pt );
    GridNode<dim>* AddNode( const csmp::Point<dim>& pt );
  
    /// prints object state to screen
    void Out() const;

  private:

    CSMP_ElementSpecifications          fem_specs_;
    std::map<csmp::Point<dim>,size_t>   points_;
    std::vector<csmp::GridNode<dim>*>   grid_nodes_;
};

/**

@class PolygonGrid  PolygonGrid "PolygonGrid.h"
@author R. Manasipov
@date 2015

The PolygonGrid contains Points and GridNodes.

*/

}// end namespace csmp

#endif

