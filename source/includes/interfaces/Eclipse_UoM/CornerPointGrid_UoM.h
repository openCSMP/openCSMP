#ifndef CORNER_POINT_GRID_UOM_H
#define CORNER_POINT_GRID_UOM_H

#include "VSet.h"
#include "ModelTopology.h"
#include "CornerPointCell_UoM.h"
#include "CSMP_highLevelUtilities.h"

namespace csmp {

namespace eclipse {

class  CellCenteredGrid
{
public:

    CellCenteredGrid();
    ~CellCenteredGrid();

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
	

/**

@class CornerPointGrid_UoM  CornerPointGrid_UoM
@author A.J. Bromage
@date 2018

*/
class  CornerPointGrid_UoM
{
  public:
    CornerPointGrid_UoM();
    ~CornerPointGrid_UoM();

    void Clear();

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

    /// initialise the pillar structure
    void Initialise( size_t NX, size_t NY, size_t NZ );
  
    template<class VarType>
    void WritePropertyToVSet( csmp::VSet<3U>&            vset,
                              const std::vector<VarType>& prop_data,
                              const std::string&          prop_name,
                              const csmp::PLACEMENT&      prop_place ) const;
  
  /// return the number of blocks the grid has in the given direction
  size_t DimensionI() const { return NX_; }
  size_t DimensionJ() const { return NY_; }
  size_t DimensionK() const { return NZ_; }
  
  void AssignDimensions(size_t nx, size_t ny, size_t nz);

  void Resize( size_t nx, size_t ny );
  Pillar& operator()( size_t i, size_t j );
  
  std::vector<uint8_t>& GetCellActivity();
  
  std::multimap<ijk,size_t>& IJKMap()
  {
    return elementMap;
  }


protected:
    void ConvertFromReservoirToCSMPcoordinateSystem( csmp::Point<3U>& pt );
	
  
private:
  size_t NX_; ///< max index of cell in x direction
  size_t NY_; ///< max index of cell in y direction
  size_t NZ_; ///< max index of cell in z direction
  
  std::vector<Pillar> pillars_;
  std::map<std::pair<size_t,size_t>,Column> columns_;
  std::multimap<ijk, size_t> elementMap;

  // These are the model-building steps in order
  void InitializeGridSpecs();
  void ConstructPillarsAndColumns(const std::vector<double64>& zcorn);
  void ConstructFiniteElementsFromColumns(VSet<3U>& vset);

  // To treat the triangular element above and beneath 
  void splitElementToPyramids();		// TODO: split irregular element to pyramids or tetrahedrons
  void splitElementToTetrahedrons();	// TODO: split irregular element to pyramids or tetrahedrons
  void addElementToMap(size_t i, size_t j, size_t k, size_t elementID);
  //void addNodeListToPList(std::map<size_t, std::vector<size_t>>& plist, std::vector<std::vector<size_t>> nodeLists, size_t& elementID);
  //void addNodeListToPList(std::map<size_t, std::vector<size_t>>& plist, std::vector<size_t> nodeList, size_t& elementID);

  size_t active_elements_;
  std::vector<uint8_t> cell_activity_;           ///< active/inactive cells

  uint8_t& CellActivity(size_t i, size_t j, size_t k);

  size_t NX_x_NY_; ///< specifications of mesh derived from grid
  size_t elements_;
  std::vector<csmp::Point<3U> >  axes_; ///< reservoir coordinate system
};

} // eclipse

}// end namespace csmp

#endif

