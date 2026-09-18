// Copyright © 2021 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_POINT_PROPERTY_TO_CELL_MAPPER_2D_H
#define CSMP_POINT_PROPERTY_TO_CELL_MAPPER_2D_H

#include "CSMP_definitions.h"
#include "DataTable.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;
template<uint32_t> class Node;
template<uint32_t> class Point;

/**
       Using the point data in the supplied spreadsheet, PointPropertyToCellMapper2D
       creates a vector of cell pointers that are non-null if they contain points are are
       part of the target region. This set is then used for the mapping of values to the mesh
       or back to the point set.
       
       @author SKM
       @date 5/1/2022
*/
class PointPropertyToCellMapper2D {
  public:
    /// reads property data from CSV text file
    explicit PointPropertyToCellMapper2D( std::string point_property_csv_file );
    
    /// initialises vector determining point containment for cells of the target region; points not found become null cells; returns whether all region cells have points
    bool MapPointsToCells( Model<2>&, std::string target_region, std::vector<Element<2>*>& );
     
    /// Assigns property data to the elements of the target region, and interpolates them across it if some points do not have corresponding cells
    void MapPointDataToElements( Model<2>&, std::string target_region, std::string target_variable );
     
    /// Interpolates the values of the target variable, to the points stored in the DataTable, storing them in there provided that the target variable is defined in the table
    // void MapNodeToPointData( Model<2>&, std::string target_region, std::string target_variable );
    void MapNodeToPointData( Model<2>&, std::string target_variable );
    
    /// writes point dataset to a CSV file
    void OutputPointDataToCSV_File( std::string file_name, std::string target_variable ) const;
    
    /// all points stored in the Mapper irrespective of their association with regions
    size_t TotalPoints() const;
    
    /// current count of cells which contain one or more of the points stored by the Mapper
    size_t CellsWithPoints() const;
    
    /// prints datatable indicating for which points mesh cells were identfied
    void Out() const;
     
  private:
    DataTable                 property_data_;     ///<  row x col table with point coordinates and property values
    std::vector<Element<2>*>  cells_with_points_; ///< pointers to the cells that are non-null if they contain points and are part of target region
};


  /// returns a vector of pointers to the cells that contain the points or NULL if there cells cannot be found
  template<uint32_t dim>
  std::vector<Element<dim>*>  findCellsEnclosingPoints( const Model<dim>&,
                                                        std::string target_region,
                                                        const std::vector<Point<dim> >& points_to_search );
   

  /// returns true if point is contained in the element or - if not - a pointer to the neighbor cell that is closest to the point of interest
  template<uint32_t dim>
  std::pair<Element<dim>*,bool>  containsPoint( Element<dim>*, const Point<dim>&, size_t region_idx );
  
  /// returns the Element's node that is closest to the supplied point
  template<uint32_t dim>
  Node<dim>*  nearestNode( Element<dim>*, const Point<dim>& );
    


} // end csmp

#endif /* CSMP_POINT_PROPERTY_TO_CELL_MAPPER_2D_H */
