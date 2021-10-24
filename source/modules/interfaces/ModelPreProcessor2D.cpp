//
//  ModelPreProcessor2D.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 24/5/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "ModelPreProcessor2D.h"
#include "Model.h"
#include "Region.h"
#include "FiniteDifferenceGrid.h"
#include "FemToGridVisitor.h"

using namespace std;

namespace csmp {

ModelPreProcessor2D::ModelPreProcessor2D( const string& csmp_model )
 : model_ptr_(new Model<2U>( csmp_model ))
 {
     cout <<"\nModelPreProcessor2D: constructed model' "<< csmp_model <<"' with specifications:\n";
     printModelDimensions( *model_ptr_ );
     model_ptr_->RegionsOut();
     model_ptr_->BoundariesOut();
     model_ptr_->SplitBoundariesOut();

 } // end csmp model constructior
    
    
    
    
    
ModelPreProcessor2D::~ModelPreProcessor2D()
 {
    delete model_ptr_;
 }



/**
    writes log10 of permeability to regular grid; works only or quadrilateral (regular) meshes
    
        Uses file format that can be read by Quadrilaterator or Triangulator or
        the function VSet<2U> readTextPixelData();
        
        TODO: what checks can be imposed to verify that the grid is indeed regular?
*/
void ModelPreProcessor2D::ModelToMatrix( const string& file_name, size_t rows, size_t columns ) const
 {
    const size_t n_cells_y(rows), n_cells_x(columns);
    cout <<"\nModelPreProcessor2D::ModelToMatrix: creating: "<< file_name << " textfile, listing an integer for each of ";
    cout << n_cells_y <<"x"<< n_cells_x<<"="<< n_cells_y * n_cells_x <<" cells forming the regular grid.\n";

    Point<2U>  xyz_min, xyz_max;
    model_ptr_->MinMaxCoordinates( xyz_min, xyz_max );
         //                             xmin        xmax        ymin        ymax       
    FiniteDifferenceGrid  regular_grid( xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], 
                                       (xyz_max[0]-xyz_min[0]) / static_cast<double64>(n_cells_x), 
                                       (xyz_max[1]-xyz_min[1]) / static_cast<double64>(n_cells_y) );  

    // 3. setting up static visitor to write repeatedly write data to grid
    // -------------------------------------------------------------------
    FemToGridVisitor<2U> writer( model_ptr_->Database(),  regular_grid, 
                                 "permeability", model_ptr_->Region("Model").Elements() );
    writer.OverWrite( true );
    writer.OutputProperty( "permeability" );
     
    // 4. Writing data onto grid 
    // -------------------------
    model_ptr_->Accept( writer );
    
    // 5. Writing grid to JPG file
    // ---------------------------
    const bool gray_scale(false);
    const bool sqrt_of_value(true);

    regular_grid.SaveToJPG( file_name.c_str(), static_cast<int32>(0), gray_scale, sqrt_of_value );
    
    const bool with_frame(false);
    regular_grid.Out( file_name.c_str(), static_cast<int32>(0), with_frame );

/*
    ofstream ofs( file_name );
 
    for (  size_t i=0U; i<n_cells_y; ++i ) {
         for ( size_t j=0U; j<n_cells_x; ++j ) ofs <<"1\t";
         ofs <<"\n";
      }

    ofs.close();
*/    
    cout << "\nModelPreProcessor2D::ModelToMatrix '" << file_name << "' written successfully.\n";

 } // end







} // end csmp
