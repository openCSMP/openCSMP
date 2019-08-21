//
//  RectangleToCrossSectionConverter.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 21/8/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#ifndef RectangleToCrossSectionConverter_hpp
#define RectangleToCrossSectionConverter_hpp

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

/**
    Created to create a regular grid of points for the sampling of a geomodel.
    The points are numbered row by row from the top down to decreasing Y coordinates,
    Column indices increase in the X direction.
*/
class RectangularGrid {
  public:
    /// specifies the dimensions of the target grid
    RectangularGrid( size_t n_cells_x, size_t n_cells_y );

    /// transfers regular grid in local coordinates to a potentially distorted grid in physical space coordinates
    void ProjectCellCentersToCrossSection( const Point<3U>& lower_left,
                                           const Point<3U>& lower_right,
                                           const Point<3U>& upper_right,
                                           const Point<3U>& upper_left,
                                           std::vector<Point<3U> >& global_coordinates ) const;
  
    /// writes the supplied cell center coordinates to a comma-delimited ascii file
    void Out( const char* file_name, const std::vector<Point<3U> >& cell_center_coordinates ) const;
  
  private:
    /// accessor of the cell-center locations in the grid; @return r,s,0
    Point<3U> operator()( size_t i, size_t j ) const;
  
    /// as numbered starting in the upper left corner, finishing in the lower right one
    size_t CellNumber( size_t i, size_t j ) const;

  private:
    const size_t   n_cells_x_, n_cells_y_;                     ///< target dimensions of grid
    const double64 dx_, dy_;                                   ///< local grid spacing resulting from cell number for fixed local dimensions
//    std::vector<pair<double64,double64> > local_coordinates_;  ///< implicit indexing, x,y coordinates ranging from -1 to 1
};

} // end csmp

#endif /* RectangleToCrossSectionConverter_hpp */
