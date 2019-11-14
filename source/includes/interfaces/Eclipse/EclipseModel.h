#ifndef ECLIPSE_MODEL_UOM_H
#define ECLIPSE_MODEL_UOM_H

#include <unordered_map>

#include "CSMP_highLevelUtilities.h"
#include "EclipseInterface.h"
#include "Model.h"


namespace csmp {

/**

@brief Input interface which converts a corner-point (hexahedral cell) grid from 
Schlumberger's reservoir simulator Eclipse into a CSMP model object in which
the cells are represented by hexahedra. Degenerate cells are converted to other
element types like prisms and tetrahedra.

@author R. Manasipov
@author Stephan Matthai (refactoring in progress)
@date 2014, 2016
@revised SKM 12/04/2018

*/
class EclipseModel : public csmp::Model<3U> {
  public:
      /// CSMP model construction from Eclipse files; property database is created from "variables_file.txt"
      EclipseModel( EclipseModelSettings& settings,
                    const std::string& model_name,
                    const std::string& variables_file );

      /// CSMP model construction from mesh files; empty property database
      EclipseModel( EclipseModelSettings& settings, const std::string& model_name );

      virtual ~EclipseModel();

      /// volumetric model subdomains if any
      template<class Container>  void GetRegions( Container& data );
  
      /// lower-dimensional element regions representing geological faults
      template<class Container>  void GetFaults( Container& data );
  
      /// wells which do not have a discrete representation inside of Eclipse
      template<class Container>  void GetWells( Container& data );

      /// JC: remove it later - processing special regions
      void CreateBoundariesAroundFaults( bool keep_fault_regions = false );
  
      /// JC: remove it later - TODO: not implemented yet
      void CreateSplitBoundariesAroundFaults( bool delete_fault_regions = false );
  
	  /// JC: remove it later -  BOX flag nodes and elements of volumetric target region
      void AssignBoxBoundaryFlagsWherePossible( const char* target_region );
  
      /// Get ijk coordinates of corner-point grid cell that corresponds to csmp::Element
      ijk EclipseCoordinates( Element<3u>* e ) const { return elmt_to_ijk_.find(e)->second; }
  
	  /// inserts a well path that penetrates the centers of the faces of the supplied cells
	  void AddWell(const std::string& well_name, const Point<3U>& well_start_point, const Point<3U>& well_end_point);

  public: // SKM accessors
      /// the dimensions of the original corner-point grid
      void GridDimensions( size_t& max_I, size_t& max_J, size_t& max_K ) const
        { max_I=grid_dim_I_; max_J=grid_dim_J_; max_K=grid_dim_K_; }

  private:
      /// Master method to build the model
      void Initialize();

      // PROPS and other specs from RUNSPECS file
	  EclipseInterface			mesh_interface;
      EclipseModelSettings      eclipse_model_settings_;
      std::set<std::string>     regions_;
      std::set<std::string>     faults_;
      std::set<std::string>     wells_;
  
      // Eclipse grid dimensions
      size_t  grid_dim_I_, grid_dim_J_, grid_dim_K_;
  
      std::unordered_multimap<ijk,Element<3u>*>  ijk_to_elmt_; ///< stores mapping from i,j,k to elements
      std::unordered_map<Element<3u>*,ijk>       elmt_to_ijk_; ///< stores mapping from elements to i,j,k
};

} // csmp
#endif
