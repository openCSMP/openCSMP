//
//  GridIJK.h
//  Open ACGSS
//
//  Created by Stephan Matthai on 2/8/2024.
//

#ifndef CSMP_GRID_IJK_H
#define CSMP_GRID_IJK_H

#include "CornerPointGrid.h"

namespace csmp {

template<uint32_t> class Element;

// TODO: work in progress; a suitable constructor is missing
/**
      Index mapper for the retrieval of cells by the I,J,K coordinates of the original model.
      
      @attention deadcells are flagged allowing for testing.
      
      @note contained is based on hashmap and best suited for single creation followed by lookup.
*/
class GridIJK {
   public:

      /// the dimensions of the original corner-point grid
      void GridDimensions( size_t& max_I, size_t& max_J, size_t& max_K ) const
        { max_I=grid_dim_I_; max_J=grid_dim_J_; max_K=grid_dim_K_; }

      /// Get ijk coordinates of corner-point grid cell that corresponds to csmp::Element
      ijk EclipseCoordinates( const Element<3u>* e ) const { return elmt_to_ijk_.find(e)->second; }

      /// activitate or deactivate cells
      uint8_t& CellActivity( uint32_t i, uint32_t j, uint32_t k ) {
            assert( i < grid_dim_I_ && j < grid_dim_J_ && k < grid_dim_K_ );
            return cell_activity_[i + j * grid_dim_I_ + k * grid_dim_I_ * grid_dim_J_ ];
         }

   private:
   
     size_t  grid_dim_I_, grid_dim_J_, grid_dim_K_;           ///< corner-point grid dimensions
     std::unordered_multimap<ijk,Element<3u>*>  ijk_to_elmt_; ///< stores mapping from i,j,k to elements
     std::unordered_map<Element<3u>*,ijk>       elmt_to_ijk_; ///< stores mapping from elements to i,j,k
     std::vector<uint8_t> cell_activity_;                     ///< active/inactive cells
   
 };

} // end csmp

#endif /* CSMP_GRID_IJK_H */
