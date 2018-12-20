#ifndef CORNER_POINT_GRID_UOM_H
#define CORNER_POINT_GRID_UOM_H

#include "VSet.h"
#include "ModelTopology.h"
#include "CornerPointCell_UoM.h"
#include "CSMP_highLevelUtilities.h"

namespace csmp {

namespace eclipse {

class CellGenerator;

class  CellCenteredGrid
{
public:

  CellCenteredGrid();
  ~CellCenteredGrid();

  void Clear();

  void AssignDimensionX( size_t NX );
  void AssignDimensionY( size_t NY );
  void AssignDimensionZ( size_t NZ );

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


/**

@class CornerPointGrid_UoM  CornerPointGrid_UoM
@author A.J. Bromage
@date 2018

*/
class  CornerPointGrid_UoM {
public:
  CornerPointGrid_UoM();
  ~CornerPointGrid_UoM();

  void CreateModel( const std::string&     model_name,
                    csmp::VSet<3U>&       vset,
                    csmp::ModelTopology&   model_topology,
                    const std::vector<double64>& zcorn,
                    std::set<std::string>& regions,
                    std::set<std::string>& faults,
                    std::set<std::string>& wells,
                    bool tetra_mesh,
                    bool exclude_inactive_cells
                    );

  template<class VarType>
  void WritePropertyToVSet( csmp::VSet<3U>&            vset,
                            const std::vector<VarType>& prop_data,
                            const std::string&          prop_name,
                            const csmp::PLACEMENT&      prop_place ) const;

  /// return the number of blocks the grid has in the given direction
  size_t DimensionI() const { return NX_; }
  size_t DimensionJ() const { return NY_; }
  size_t DimensionK() const { return NZ_; }

  void AssignDimensions( size_t nx, size_t ny, size_t nz );

  void Resize( size_t nx, size_t ny );
  Pillar& operator()( size_t i, size_t j );

  std::vector<uint8_t>& GetCellActivity();

  size_t OrdinaryNodes() const { return ordinaryNodes_; }

  std::multimap<ijk, size_t>& IJKMap() { return elementMap; }

  void ConvertFromReservoirToCSMPcoordinateSystem( csmp::Point<3U>& pt ) const;

  bool ConstructLineElement( size_t& i, size_t& j, size_t& k, size_t& new_elemt_idx );
  bool ConstructLineElement( ColumnCell& cell, size_t& i, size_t& j, size_t& k );
  CellGenerator* GetCellGenerator() { return generator_; }

private:
  size_t NX_; ///< max index of cell in x direction
  size_t NY_; ///< max index of cell in y direction
  size_t NZ_; ///< max index of cell in z direction

  std::vector<Pillar> pillars_;
  std::map<std::pair<size_t, size_t>, Column> columns_;
  std::multimap<ijk, size_t> elementMap;
  CellGenerator* generator_;

  size_t badHexahedra_;
  size_t badPyramids_;
  size_t badTetrahedra_;
  size_t badPrisms_;
  size_t badLines_;

  // These are the model-building steps in order
  void InitializeGridSpecs();
  void ConstructPillarsAndColumns( const std::vector<double64>& zcorn );
  void ConstructFiniteElementsFromColumns( VSet<3U>& vset );
  bool ConstructEclipseCell0000( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_HEXAHEDRON
  bool ConstructEclipseCell0001( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMIDS_310_312
  bool ConstructEclipseCell0010( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMIDS_201_203
  bool ConstructEclipseCell0011( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PRISM_23
  bool ConstructEclipseCell0100( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMIDS_130_132
  bool ConstructEclipseCell0101( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_TETRAHEDRONS_130_132
  bool ConstructEclipseCell0110( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PRISM_12
  bool ConstructEclipseCell0111( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMID_0
  bool ConstructEclipseCell1000( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMIDS_021_023
  bool ConstructEclipseCell1001( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PRISM_03
  bool ConstructEclipseCell1010( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_TETRAHEDRONS_021_023
  bool ConstructEclipseCell1011( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMID_1
  bool ConstructEclipseCell1100( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PRISM_01
  bool ConstructEclipseCell1101( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMID_2
  bool ConstructEclipseCell1110( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ); //Eclipse cell type: ECLIPSE_CELL_PYRAMID_3

                                                                                                                                       // To add the degenerated elements into the map storage																															    
  void addElementToMap( size_t i, size_t j, size_t k, size_t elementID );

  size_t active_elements_;
  std::vector<uint8_t> cell_activity_;           ///< active/inactive cells

  uint8_t& CellActivity( size_t i, size_t j, size_t k );

  size_t NX_x_NY_; ///< specifications of mesh derived from grid
  size_t elements_, ordinaryNodes_;
  std::vector<csmp::Point<3U> >  axes_; ///< reservoir coordinate system
};
} // eclipse

} // end namespace csmp

#endif

