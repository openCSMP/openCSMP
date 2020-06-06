//
//  ModelPreProcessor2D.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 24/5/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "ModelPreProcessor2D.h"
#include "Model.h"

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
*/
void ModelPreProcessor2D::ModelToMatrix( const string& file_name ) const
 {
    size_t n_cells_y, n_cells_x;
    cout <<"\nModelPreProcessor2D::ModelToMatrix: creating: "<< file_name << " textfile, listing an integer for each of ";
    cout << n_cells_y <<"x"<< n_cells_x<<"="<< n_cells_y * n_cells_x <<" cells forming the regular grid.\n";

    ofstream ofs( file_name );
 
    for (  size_t i=0U; i<n_cells_y; ++i ) {
         for ( size_t j=0U; j<n_cells_x; ++j ) ofs <<"1\t";
         ofs <<"\n";
      }

    ofs.close();
    cout << "\nModelPreProcessor2D::ModelToMatrix '" << file_name << "' written successfully.\n";

 } // end







} // end csmp
