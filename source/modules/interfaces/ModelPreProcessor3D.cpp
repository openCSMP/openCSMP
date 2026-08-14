//
//  ModelPreProcessor3D.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 24/5/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "ModelPreProcessor3D.h"
#include "Model.h"
#include "VTU_Interface.h"

using namespace std;

namespace csmp {

ModelPreProcessor3D::ModelPreProcessor3D( const string& csmp_model )
 : model_ptr_(new Model<3U>( csmp_model ))
 {
     cout <<"\nModelPreProcessor3D: constructed model' "<< csmp_model <<"' with specifications:\n";
     printModelDimensions( *model_ptr_ );
     model_ptr_->RegionsOut();
     model_ptr_->BoundariesOut();
     model_ptr_->SplitBoundariesOut();

 } // end csmp model constructior
    

ModelPreProcessor3D::~ModelPreProcessor3D()
 {
    delete model_ptr_;
 }



/**
*/
void ModelPreProcessor3D::CheckIntegrity( const std::string& report_file /* default name: model.Name()-integrity_check.txt" */ ) const
 {
     // output to VTU for visual examination
     string property_label_file = (!report_file.empty()) ? report_file : string(model_ptr_->Name()) + "-integrity_check.txt";
     
     const bool write_VTU_output(true);
     if ( write_VTU_output ) {
          const list<string> output_vars = { "permeability", "porosity" };
          VTU_Interface<3U>  vtu_output( *model_ptr_ );
          vtu_output.OutputDataToVTU( property_label_file.c_str(), output_vars, "Model", 0 );
       }

 } // end CheckIntegrity
    

} // end csmp
