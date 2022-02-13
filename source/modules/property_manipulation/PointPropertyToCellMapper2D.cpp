//
//  PointPropertyToCellMapper2D.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 30/12/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "PointPropertyToCellMapper2D.h"
//#include "TextFileIO.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "Point.h"
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#include "PDE_Integrator_UoM.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"

using namespace std;

namespace csmp {

// non member functions

/**
    Using the supplied cell as a starting location,
    traverses the mesh until the target point is found or a model boundary is encountered
    recording nullptr
*/
template<uint32_t dim>
std::vector<Element<dim>*>  findCellsEnclosingPoints( const Model<dim>& model,
                                                      string target_region,
                                                      const vector<Point<dim> >& points_to_search )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
          
    if ( points_to_search.empty() )
      csmp_error.notice( ERROR, "WhichCellContainsPoint<dim>::WhichCellContainsPoint", "no points were supplied, visits will have no effect");


    // marking the elements of the region with the unique key UINT_MAX so that containment is readily detected
    const Region<dim>& target_domain(model.Region(target_region.c_str()));
    
    for ( auto it=target_domain.ElementsBegin(); it!=target_domain.ElementsEnd(); ++it )
      (*it)->Idx( UINT_MAX );

   // searching for the points
   vector<Element<dim>*>  cells_containing_points;
   Element<dim>*          last_cell_found(nullptr);
   
   cells_containing_points.reserve( points_to_search.size() );
   
   for ( auto& pt : points_to_search )
     {
        // the last cell found is the starting location for the next search
        Element<dim>* cell = (last_cell_found) ? last_cell_found : target_domain.E(0);
         
        // searching
        bool target_cell_found{false};
        while (  !target_cell_found ) {
             pair<Element<dim>*,bool>  result = containsPoint( cell, pt, cell->Idx() );
             // assigning the element found
             if ( result.second ) {
                  last_cell_found   = result.first;
                  target_cell_found = true;
               }
              cell = result.first;
          }
        if ( target_cell_found )
          cells_containing_points.push_back( last_cell_found );
        else
          cells_containing_points.push_back( nullptr );
     }
     
   return cells_containing_points;

 } // end findCellsEnclosingPoints

template vector<Element<2>*>  findCellsEnclosingPoints( const Model<2>&, string, const vector<Point<2>>& );





/** returns true if point is contained in the element or - if not - a pointer to the neighbor cell that is closest to the point of interest
 
 @todo make this method general
 */
template<uint32_t dim>
pair<Element<dim>*,bool>  containsPoint( Element<dim>* elmt, const Point<dim>& point, size_t region_idx )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    // only helps if the pointer was indeed set to zero
    assert( elmt != nullptr );
    assert( isTriangularElement(elmt->FE_Type()) ||
            elmt->FE_Type() == ISOPARAMETRIC_LINEAR_TETRAHEDRON ); // equidimensional simplex elements only
    
    // getting the element interpolation function values; if all are positive, the point is contained in this element
    vector<double> IPOL;
    if constexpr( dim == 3 ) elmt->N_AtGlobalPoint( IPOL, vector<double>{point[0],point[1],point[2]} );
    if constexpr( dim == 2 ) elmt->N_AtGlobalPoint( IPOL, vector<double>{point[0],point[1]} );
    map<double,size_t>  nearest_nbor;
    bool                point_is_inside{true};

    const size_t n_nodes{ IPOL.size() };
    for ( auto i{0}; i<n_nodes; ++i ) {
         nearest_nbor.insert( make_pair( IPOL[i], i ) );
         if ( IPOL[i] < 0. ) point_is_inside = false;
      }
      
    // if the point is contained in this element, 'true' and the element pointer are returned
    if ( point_is_inside ) {
         if ( region_idx != elmt->Idx() ) {
              csmp_error.notice( ERROR, "containsPoint", "point in cell that is not part of region is ignored");
           }
         else return make_pair( elmt, true );
      }
    
    // if one of the interpolation functions is negative, the point is not contained
    // and the neighbor element opposite of the node with the most negative interpolation function value is returned if any (else nullptr)
    return make_pair( elmt->Neighbor( (*nearest_nbor.begin()).second ), false );
 
 } // end
  
template pair<Element<3>*,bool>  containsPoint( Element<3>*, const Point<3>&, size_t );
template pair<Element<2>*,bool>  containsPoint( Element<2>*, const Point<2>&, size_t );
  
  
  
  
  
/// returns the Element's node that is closest to the supplied point
template<uint32_t dim>
Node<dim>*  nearestNode( Element<dim>* e, const Point<dim>& point )
 {
    // computes distances between point and the elements node points returning the closest node
    // distance, local node number
    map<double,size_t>  node_distances;
    const size_t n_nodes{ e->Nodes() };
    for ( auto i{0}; i<n_nodes; ++i )
      node_distances.insert( make_pair( point.DistanceTo( e->N(i)->Coordinate()),i) );
    
    return e->N( (*node_distances.begin()).second );
    
 } // end

template Node<3>*  nearestNode( Element<3>*, const Point<3>& );
template Node<2>*  nearestNode( Element<2>*, const Point<2>& );






PointPropertyToCellMapper2D::PointPropertyToCellMapper2D( string csv_file )
 : property_data_( csv_file, vector<string>{"x","y"}, vector<vector<double>>{ 1, vector<double>(2) } )
 {
    vector<string>          row_labels, col_headers;
    vector<vector<double>>  rows_of_columns;
    
    size_t n_entries = read_CSV_File( csv_file, row_labels, col_headers, rows_of_columns );
    
    property_data_ = move( DataTable{ csv_file, row_labels, col_headers, rows_of_columns } );
    
    // lazy evaluation of: cells_with_points_
    cells_with_points_.resize( n_entries, nullptr );
    
    cout <<"\nPointPropertyToCellMapper2D: read "<< n_entries <<" records from file '"<< csv_file <<"'\n'";
 }



/**
   Initialises cell-pointer vector determining point containment for cells of the target region;
   points not found become null cells;
   
   @return returns true when all cells of the target region contain points.

*/
bool PointPropertyToCellMapper2D::MapPointsToCells( Model<2>& model, string target_region, vector<Element<2>*>& cells )
 {
    const Region<2>& target = model.Region(target_region);
 
    // 1. creating a search list of Point objects from the data in the table
    // ---------------------------------------------------------------------
    const size_t cx = property_data_.ColumnIndex("x");
    vector<Point<2> >  points_to_search;
    points_to_search.reserve( property_data_.Rows() );
    for ( auto i{0}; i<property_data_.Rows(); ++i )
      points_to_search.emplace_back( Point<2>( property_data_(i,cx), property_data_(i,cx+1) ) );
      
    
    // 2. searching the points in the model
    // ------------------------------------
    cells_with_points_ = findCellsEnclosingPoints( model, target_region, points_to_search );

    return ( CellsWithPoints() == target.Elements() );

 } // end MapPointsToCells
     






/// all points stored in the Mapper irrespective of their association with regions
size_t PointPropertyToCellMapper2D::TotalPoints() const
 {
    return property_data_.Rows();
 }


/// current count of cells which contain one or more of the points stored by the Mapper
size_t PointPropertyToCellMapper2D::CellsWithPoints() const
 {
    size_t non_null_cells{0};
    for ( const auto it : cells_with_points_ )
      if ( it != nullptr ) non_null_cells++;
      
    return non_null_cells;
 }
     



/**
     Attempts to assign the property data as element properties of the target region, in the process 'cells_with_points_' gets initialised
*/
void PointPropertyToCellMapper2D::MapPointDataToElements( Model<2>& model, string target_region, string target_variable )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    const csmp::Index prop_key = model.Database().StorageKey( target_variable.c_str() );
    if ( prop_key.place != ELEMENT || prop_key.type != SCALAR )
      csmp_error.notice( ERROR, "PointPropertyToCellMapper2D::MapPointDataToElements", target_variable,
                        "should be a scalar placed on the Element");
    
    // 1. searching the points in the model
    // ------------------------------------
    bool extrapolation_needed = MapPointsToCells( model, target_region, cells_with_points_ );
    
  
   // 2. assigning the property values extracted from the CSV file
   // ------------------------------------------------------------
   // (for the elements that contain the points, the properties are assigned directly, to get props on nodes,
   //  for computations, the property values are extrapolated to the nodes)

    // creating a search list of Point objects from the data in the table
    vector<Point<2> >  points_to_search;
    vector<double>     point_values;
    // finding the column index of the target variable in the table
    const size_t cx = property_data_.ColumnIndex("x");
    const size_t cy = property_data_.ColumnIndex("y");
    const size_t prop_idx = property_data_.ColumnIndex( target_variable );
    // creating the point search data
    points_to_search.reserve( property_data_.Rows() );
    point_values.reserve( property_data_.Rows() );
    for ( auto i{0}; i<property_data_.Rows(); ++i ) {
         // assuming that the x, y coordinates reside in column 0 and 1
         points_to_search.emplace_back( Point<2>( property_data_(i,cx), property_data_(i,cy) ) );
         point_values.emplace_back( property_data_(i,prop_idx) );
      }

   if ( !extrapolation_needed ) {
        for ( size_t elmt{0}; elmt < cells_with_points_.size(); ++elmt )
          if ( cells_with_points_[elmt] != nullptr )
            cells_with_points_[elmt]->Store( prop_key, makeScalar( FIELD_DATA, point_values[elmt] ) );
     }

   // 3. extrapolation to nodes, followed by computation of smooth property distribution in Region
   // --------------------------------------------------------------------------------------------
   // (if dof too small, nearest neighbor extrapolation is used)
   else {
        // for the elements that have data, these are extrapolated to the nodes
        const string target_variable_node = target_variable + " node";
        model.CreateProperty( target_variable_node.c_str(), "extrapolated", SCALAR, NODE );
        const csmp::Index prop_key_node = model.Database().StorageKey( target_variable_node.c_str() );
        // assigning Dirichlet constraints to the nearest nodes
        for ( size_t elmt{0}; elmt < cells_with_points_.size(); ++elmt )
          if ( cells_with_points_[elmt] != nullptr ) {
               Node<2>* node = nearestNode( cells_with_points_[elmt], points_to_search[elmt] );
               node->Store( prop_key_node, makeScalar( DIRICH, point_values[elmt] ) );
            }
        // extrapolating the values solving Laplace equation
        SAMG_Settings                    settings;
        SAMG_Solver                      samg_solver( &settings );
        PDE_Integrator_UoM<2,Region>     extrapolator( samg_solver );
        NumIntegral_dNT_dN_dV<2>      lhs( model.Database(), target_variable_node.c_str(), target_variable_node.c_str() );
        NumIntegral_SetRHS_to_Zero<2> rhs( model.Database(), target_variable_node.c_str() );
        extrapolator.Add( &lhs );
        extrapolator.Add( &rhs );
        // computation
        Region<2> target_domain(model.Region(target_region));
        extrapolator.IntegrateOver( target_domain );

        // interpolating nodal values back to the element barycentres in the region
        for ( auto it=target_domain.ElementsBegin(); it!=target_domain.ElementsEnd(); ++it )
          if ( (*it)->Status(prop_key) != FIELD_DATA )
            {
               double prop_val = (*it)->PropertyValueAtBaryCenter( prop_key_node );
               (*it)->Store( prop_key, makeScalar( FIELD_DATA, prop_val ) );
            }
     }
   
 } // end MapPointDataToElements

 
 
 
 
/** Interpolates the values of the target node variable, to the points specified in the CSV file, writing another CSV file with these values
*/
void PointPropertyToCellMapper2D::MapNodeToPointData( Model<2>& model, string target_region, string target_variable )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    Region<2U>    target_domain(model.Region(target_region.c_str()));

    // assuming that which cells contain the points has already been established
    if ( CellsWithPoints() < target_domain.Elements() )
      csmp_error.notice( ERROR, "PointPropertyToCellMapper2D::MapNodeToPointData",
                         "Not all cells of the target domain contain points");
        
    const csmp::Index eprop_key = model.Database().StorageKey(target_variable.c_str());
    const csmp::Index nprop_key = model.Database().StorageKey((target_variable + " node").c_str());

    if ( eprop_key.type != SCALAR || nprop_key.type != SCALAR )
      csmp_error.notice( ERROR, "PointPropertyToCellMapper2D::MapNodeToPointData",
                         "node and element variables must be SCALAR types");
       
    if ( eprop_key.place != ELEMENT || nprop_key.place != NODE )
      csmp_error.notice( ERROR, "PointPropertyToCellMapper2D::MapNodeToPointData",
                         "the target variable must be placed on the ELEMENT with a node version in addition");
       
    // interpolating the node values to the barycentres of the elements
    if ( !property_data_.ContainsColumn( target_variable) ) property_data_.AppendColumn( target_variable );
    const size_t col_idx = property_data_.ColumnIndex( target_variable );
    
    for ( size_t point{0}; point < property_data_.Rows(); ++point )
      if ( cells_with_points_[point] != nullptr ) {
           double prop_val = cells_with_points_[point]->PropertyValueAtBaryCenter( nprop_key );
           property_data_(point,col_idx) = prop_val;
        }
 
 } // end MapNodeToPointData
 
 
 
 
    
    /// writes point dataset to a CSV file
void PointPropertyToCellMapper2D::OutputPointDataToCSV_File( std::string file_name, std::string target_variable ) const
 {
    DataTable  output = property_data_.ExtractSubTable( set<string>{}, set<string>{target_variable} );
    output.Write_CSV_File( file_name );
    
 } // end OutputPointDataToCSV_File




    /// prints datatable indicating for which points mesh cells were identfied
void PointPropertyToCellMapper2D::Out() const
 {
    cout <<"\n\nPointPropertyToCellMapper2D::Out: ";
    property_data_.Out();
    
    if ( cells_with_points_.empty() ) {
         cout <<"\n\tpoint containment in cells not evaluated yet.\n";
         return;
      }
    
    cout <<"\nRows = points that are contained in cells of the model:\n";
    const size_t cx = property_data_.ColumnIndex("x");
    const size_t cy = property_data_.ColumnIndex("y");
    const size_t cz = property_data_.ColumnIndex("z");
    for ( auto i{0}; i<property_data_.Rows(); ++i ) {
         cout <<"("<< property_data_(i,cx) <<",";
         cout << property_data_(i,cy) <<",";
         cout << property_data_(i,cz) <<"): ";
         if ( cells_with_points_[i] == nullptr ) cout <<" ---\n";
         else cout << parseFiniteElementType( cells_with_points_[i]->FE_Type() ) << cells_with_points_[i]->Idx() <<"\n";
      }
 }
     



} // end csmp
