#ifndef ECLIPSE_MODEL_UOM_H
#define ECLIPSE_MODEL_UOM_H

#include <unordered_map>

#include "CSMP_highLevelUtilities.h"
#include "EclipseInterface_UoM.h"
#include "Model.h"


namespace csmp {

namespace eclipse {

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

      /// access of elements generated from corner-point cells by their i(W->E),j(S->N),k(top->bottom) grid indices
      const Element<3U>* operator()( size_t i, size_t j, size_t k ) const;
      Element<3U>*       operator()( size_t i, size_t j, size_t k );

      /// existing special regions
      template<class Container>  void GetRegions( Container& data );
  
      template<class Container>  void GetFaults( Container& data );
  
      template<class Container>  void GetWells( Container& data );

      /// processing special regions
      void CreateBoundariesAroundFaults( bool keep_fault_regions = false );
      void CreateSplitBoundariesAroundFaults( bool delete_fault_regions = false );
  
      /// BOX flag nodes and elements of volumetric target region
      void AssignBoxBoundaryFlagsWherePossible( const char* target_region );
  
  public: // SKM accessors
      /// the dimensions of the original corner-point grid
      void GridDimensions( size_t& max_I, size_t& max_J, size_t& max_K ) const
        { max_I=grid_dim_I_; max_J=grid_dim_J_; max_K=grid_dim_K_; }

  private:
      /// Master method to build the model
      void Initialize();
  
      // PROPS and other specs from RUNSPECS file
      EclipseModelSettings      eclipse_model_settings_;
      std::set<std::string>     regions_;
      std::set<std::string>     faults_;
      std::set<std::string>     wells_;
  
      // Eclipse grid dimensions
      size_t  grid_dim_I_, grid_dim_J_, grid_dim_K_;
  
      std::unordered_map<ijk,size_t>  IJK_map_; ///< stores mapping from i,j,k to elements
};

} // eclipse

} // csmp



#endif
