#ifndef CORNER_POINT_GRID_H
#define CORNER_POINT_GRID_H

#include "CSMP_definitions.h"
#include "CSMP_highLevelUtilities.h"
#include "CSMP_ElementSpecifications.h"
#include "VSet.h"
#include "ModelTopology.h"
#include "CornerPointCell.h"

namespace csmp {

/// cell-centered grid container
template<size_t dim>
class  BlockCenteredGrid // TODO: call this cell-centered
{
public:

    BlockCenteredGrid();

    BlockCenteredGrid( const BlockCenteredGrid& );
    BlockCenteredGrid& operator=( const BlockCenteredGrid& );

    ~BlockCenteredGrid();
    void Clear();

    void AssignDimensionX( size_t NX );
    void AssignDimensionY( size_t NY );
    void AssignDimensionZ( size_t NZ );
    void AssignCellCoordinatesToPillars( std::vector<std::vector<Pillar> >& pillars );

    size_t GetNumCells() const;
    std::vector<csmp::ScalarVariable>& GetCellDepths();
    std::vector<csmp::ScalarVariable>& GetCellSizes( size_t direction );

private:

    size_t NX_; // max index of cell in x direction
    size_t NY_; // max index of cell in y direction
    size_t NZ_; // max index of cell in z direction
    std::vector<csmp::ScalarVariable> dx_;   // cell sizes  ( dx )
    std::vector<csmp::ScalarVariable> dy_;   // cell sizes  ( dy )
    std::vector<csmp::ScalarVariable> dz_;   // cell sizes  ( dz )
    std::vector<csmp::ScalarVariable> tops_; // cell depths ( ztop )
};


/// corner-point based grid container
template<size_t dim>
class  CornerPointGrid
{
  public:

    CornerPointGrid();

    CornerPointGrid( const CornerPointGrid& );
    CornerPointGrid& operator=( const CornerPointGrid& );

    ~CornerPointGrid();
    void Clear();

    /// initialise the pillar structure
    void Resize( size_t i_pillar_max, size_t j_pillar_max );
  
    /// accessing the contained pillars
    Pillar& operator()( size_t, size_t t );

    void CreateModel( const std::string&     model_name,
                      csmp::VSet<dim>&       vset,
                      csmp::ModelTopology&   model_topology,
                      std::set<std::string>& regions,
                      std::set<std::string>& faults,
                      std::set<std::string>& wells
                    );

    template<class VarType>
    void WritePropertyToVSet( csmp::VSet<dim>&            vset,
                              const std::vector<VarType>& prop_data,
                              const std::string&          prop_name,
                              const csmp::PLACEMENT&      prop_place ) const;
  
    void TetraMesh( bool );
    void ExcludeInactiveCells( bool );

    void AddWellFacePath( const std::string& well_name, const std::vector<size_t>& cell_ids );
    void AddWellFacePath( const std::string& well_name, const std::vector<size_t>& cell_ids, const std::vector<std::pair<size_t,size_t> >& face_ids );
    void AddWellEdgePath( const std::string& well_name, const std::vector<size_t>& cell_ids, const std::vector<std::pair<size_t,size_t> >& edge_ids );


    const csmp::Point<3U>& GetCellNode( size_t i, size_t j, size_t k, size_t nid ) const;
    std::vector<size_t>&   GetCellActivity();
    std::vector<std::vector<Pillar> >& GetPillars();
    std::map<std::string,std::vector<std::pair<size_t,size_t> > >& GetFaultData();
    std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& GetWellFacePath();
    std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& GetWellEdgePath();

    size_t GetNumCells() const;
    void AssignDimensionX( size_t NX );
    void AssignDimensionY( size_t NY );
    void AssignDimensionZ( size_t NZ );
  
    /// return the number of blocks the grid has in the given direction
    size_t DimensionI() const { return NX_; }
    size_t DimensionJ() const { return NY_; }
    size_t DimensionK() const { return NZ_; }

    /// find the cell number from its index
    size_t CellIndex( size_t i, size_t j, size_t k ) const;
  
    /// checking whether a cell is active
    bool IsActiveCell( size_t i, size_t j, size_t k ) const
      { return (cell_activity_[ CellIndex(i,j,k) ] == 1) ? true : false; }

protected:

    void InitializeGridSpecs();

    void ProcessPillars();

    void EstablishActiveDomain( csmp::VSet<dim>& vset,
                                csmp::ModelTopology& mesh_topology,
                                std::set<std::string>& regions,
                                const std::vector<size_t>& cell_elmts,
                                const std::set<std::string>& cell_fem_types
                              );

    /// adds surface elements to VSet and fault regions to ModelTopology
    void EstablishFaultRegions( csmp::VSet<dim>& vset,
                                csmp::ModelTopology& model_topology,
                                std::set<std::string>& faults,
                                const std::vector<CornerPointCell<3U> >& poly
                                );

    void EstablishWellRegions(  csmp::VSet<dim>& vset,
                                csmp::ModelTopology& model_topology,
                                std::set<std::string>& wells,
                                const std::vector<CornerPointCell<3U> >& poly
                                );

    void AssignCellNodeToPillar( size_t i, size_t j, size_t k, size_t nid, const csmp::Point<3U>& pt );
    void ConvertFromReservoirToCSMPcoordinateSystem( csmp::Point<3U>& pt );
    void MinMaxCoordinates( PolygonGridManager<3U>& pgm, Point<3U>& xyz_min, Point<3U>& xyz_max ) const;
    void DefineAxes( PolygonGridManager<3U>& pgm );

    /// @todo should be an operator to access the cells of the structured grid
    size_t CellIndexI( size_t cell_id ) const;
    size_t CellIndexJ( size_t cell_id ) const;
    size_t CellIndexK( size_t cell_id ) const;
    size_t NodeIndexI( size_t cell_id, size_t nid ) const;
    size_t NodeIndexJ( size_t cell_id, size_t nid ) const;
    size_t NodeIndexK( size_t cell_id, size_t nid ) const;
    size_t NodeIndexIncrementI( size_t nid ) const;
    size_t NodeIndexIncrementJ( size_t nid ) const;
    size_t NodeIndexIncrementK( size_t nid ) const;

private:

    size_t NX_; ///< max index of cell in x direction
    size_t NY_; ///< max index of cell in y direction
    size_t NZ_; ///< max index of cell in z direction
    std::vector<size_t> cell_activity_;           ///< active/inactive cells
    std::vector<std::vector<Pillar> > pillars_;   ///< (subvertical) coordinate lines defining the pillars
    std::map<size_t,std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> > >     embedded_cells_;          ///< subcells
    std::map<std::string,std::vector<std::pair<size_t,size_t> > >             faults_data_;             ///< fault representations
    std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > > well_face_path_;   ///< wells through center of cell faces
    std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > > well_edge_path_;   ///< wells through center of edges
    bool   tetra_mesh_;
    bool   exclude_inactive_cells_;

    size_t NX_x_NY_; ///< specifications of mesh derived from grid
    size_t elements_;
    std::vector<csmp::Point<3U> >  axes_; ///< reservoir coordinate system

};

void addWellPath( size_t NX, size_t NY, size_t NZ,
                  const std::string& well_name,
                  const std::vector<size_t>& cell_ids,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path );

void addWellPath( const std::string& well_name,
                  const std::vector<size_t>& cell_ids,
                  const std::vector<std::pair<size_t,size_t> >& face_ids,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path );

void addWellPath( const std::string& well_name,
                  const std::vector<size_t>& cell_ids,
                  std::pair<size_t,size_t> face_ids,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path );

void addWellPath( const std::string& well_name,
                  size_t cell_ids,
                  std::pair<size_t,size_t> face_id,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path );


/**

@class CornerPointGrid  CornerPointGrid "CornerPointGrid.h"
@author R. Manasipov
@date 2015

Treating degenerate cells ( deformed hexahedron cells or cells with less than 8 nodes ).
Each degenerate cell's subdivision is based on division of quadrilateral face.
'Tetra mesh' option will lead to tetrahedral + triangular mesh, otherwise
depending on the complexity of model resulting mixed mesh can contain
hexahedrons, pyramids, prims, tetrahedrons, quadrilaterals and triangles.

*/

}// end namespace csmp

#endif

