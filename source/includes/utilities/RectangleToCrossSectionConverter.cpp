//
//  RectangleToCrossSectionConverter.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 21/8/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#include "RectangleToCrossSectionConverter.h"
#include "IsoParametricLinearQuadrilateral.h"
#include "Element.h"

using namespace std;

namespace csmp {

    /// specifies the dimensions of the target grid
RectangularGrid::RectangularGrid( size_t n_cells_x, size_t n_cells_y )
 : n_cells_x_(n_cells_x), n_cells_y_(n_cells_y),
   dx_(2./static_cast<double64>(n_cells_x)), dy_(2./static_cast<double64>(n_cells_y))
 {
 }



/// accessor of the cell-center locations in the grid
Point<3U> RectangularGrid::operator()( size_t i, size_t j ) const
 {
    assert( i < n_cells_y_ ); // rows
    assert( j < n_cells_x_ ); // columns

    const double64 offset_x(dx_/2.), offset_y(dy_/2.);
    const double64 r = j * dx_ + offset_x;
    const double64 s = i * dy_ + offset_y;
 
    return Point<3U>(r,s,0.);
 }




size_t RectangularGrid::CellNumber( size_t i, size_t j ) const
 {
    assert( i < n_cells_y_ ); // rows
    assert( j < n_cells_x_ ); // columns

    return i * n_cells_y_ + j;
 }
 


/**
   transfers regular grid in local coordinates to a potentially distorted grid in physical space coordinates
*/
void RectangularGrid::ProjectCellCentersToCrossSection( const Point<3U>& lower_left,
                                                        const Point<3U>& lower_right,
                                                        const Point<3U>& upper_right,
                                                        const Point<3U>& upper_left,
                                                        vector<Point<3U> >& global_coordinates ) const
 {
    // Creating an isoparametric linear quadrilateral element for the extrapolation of the point locations
    const size_t dim(3U);
    IsoparametricLinearQuadrilateral quadrilateral( dim );
    Element<dim> rectangle( &quadrilateral );
    rectangle.Idx(0);

    // Assigning the corner-point locations of the cross section to the nodes of the element
    Node<dim> n1, n2, n3, n4;
    n1.Idx(0U); n2.Idx(1); n3.Idx(2); n4.Idx(3);
    n1.Coordinate( lower_left );
    n2.Coordinate( lower_right );
    n3.Coordinate( upper_right );
    n4.Coordinate( upper_left );
    rectangle.Assign( 0, &n1 );
    rectangle.Assign( 1, &n2 );
    rectangle.Assign( 2, &n3 );
    rectangle.Assign( 3, &n4 );
 
    // computing the cell center locations for the output vector
    global_coordinates.clear();
    global_coordinates.reserve( n_cells_x_ * n_cells_y_ );
 
    for ( size_t i=0U; i<n_cells_y_; ++i )
      for ( size_t j=0U; j<n_cells_x_; ++j ) {
           // performing the coordinate transformation
           Point<3U> rs = (*this)(i,j);
           Point<3U> cell_center_xyz = rectangle.RstToXYZ( rs );
           global_coordinates.push_back( cell_center_xyz );
        }

 } // end ProjectCellCentersToCrossSection
 
 
 


/**
    Writes file with cell barycentre coordinates, preceded by element-id(number).

    The file format is

    header line
    axes titles: region-name(string), region-id(number 0..n-1), element-id(number 0..n-1), x, y, z
    comma-seperated data...
    EOF
 
    @author SKM
*/
void RectangularGrid::Out( const char* file_name, const std::vector<Point<3U> >& cell_center_coordinates ) const
 {
    ofstream ofs( string( file_name ) + "-element_barycentres.txt" );
    ofs << file_name << "-element_barycentres.txt  textfile listing the barcyentre coordinates of all elements in the rectangular grid.\n";
    ofs << "element number,x,y,z\n";
 
    size_t cell_number(0U);
    for ( auto it = cell_center_coordinates.begin(); it != cell_center_coordinates.end(); ++it, ++cell_number )
      {
         // element number
         ofs << cell_number <<",";
         // point coordinates (ordered and flipped to reflect SKUA's UTM lefthandrule coordinate system)
         //                   SKUA x             SKUA -y          SKUA -z
         //ofs << scientific << (*it)[0] << "," << (*it)[2] <<","<< (*it)[1] << endl;
         ofs << scientific << (*it)[0] << "," << (*it)[1] <<","<< (*it)[2] << endl;
      }
 
    ofs.close();
    cout << "\noRectangularGrid::Out: " << file_name << "-element_barycentres.txt written successfully.\n";
 
 } // end Out



} // end csmp
