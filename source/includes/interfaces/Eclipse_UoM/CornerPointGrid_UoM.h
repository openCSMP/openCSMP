#ifndef CORNER_POINT_GRID_UOM_H
#define CORNER_POINT_GRID_UOM_H

#include "VSet.h"
#include "ModelTopology.h"
#include "CornerPointCell_UoM.h"

namespace csmp {

namespace eclipse {

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

protected:
    void ConvertFromReservoirToCSMPcoordinateSystem( csmp::Point<3U>& pt );
  
private:
  size_t NX_; ///< max index of cell in x direction
  size_t NY_; ///< max index of cell in y direction
  size_t NZ_; ///< max index of cell in z direction
  
  std::vector<Pillar> pillars_;
  std::map<std::pair<size_t,size_t>,Column> columns_;

  // These are the model-building steps in order
  void InitializeGridSpecs();
  void ConstructPillarsAndColumns(const std::vector<double64>& zcorn);
  void ConstructFiniteElementsFromColumns();

  size_t active_elements_;
  std::vector<uint8_t> cell_activity_;           ///< active/inactive cells

  uint8_t& CellActivity(size_t i, size_t j, size_t k);

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

} // eclipse

}// end namespace csmp

#endif

