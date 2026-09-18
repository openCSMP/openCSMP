#include "SplitBoundaryTest_JK.h"

#include "ANSYS_Model2D.h"
#include "ComputationalSettings.h"
#include "InputDataManager.h"
#include "VTU_Interface.h"

using namespace std;

namespace csmp


/** 2D Test to create a split boundary using CreateSplitBoundaryFrom

   Regions File: rectangle_2SBs

**/
{
SplitBoundaryTest_JK::SplitBoundaryTest_JK()
    : geometry_name_("rectangle_2SBs"),
    model_ID_path_("./"),
    output_name_("Test_CreateSBFrom"),
    vars_name_("SplitBoundaryTest_JK-variables.txt")
{
    output_path = model_ID_path_;

    //! Create model and reference to property database
    model_ = new ANSYS_Model2D(geometry_name_, geometry_name_, vars_name_);
    model_->InstantiateFiniteVolumes();

    //! output variables
    output_props_.push_back("porosity"); // dummy0
    }
    
    

SplitBoundaryTest_JK::~SplitBoundaryTest_JK()
{
    //! clean memory
    delete model_;
}


void SplitBoundaryTest_JK::run()
{
//   Test_CreateSplitBoundaryFrom();
   Test_CreateSplitBoundaryFrom_Simplified();
}


//! Not really - I did not have to change anything except fore some cosmetics
void SplitBoundaryTest_JK::Test_CreateSplitBoundaryFrom()
{
    std::vector<std::string> fractureNames;
    fractureNames = {"FRACTURE1", "FRACTURE2"};

    bool retain_elmts_as_intervening_elements = false;

    //* Following SplitBoundaryInterface_Test workflow
    ////////////////////////////////////
    //! Splitboundary creation
    set<string> sb_names;
    // Loop through each fracture name and call CreateSplitBoundaryFrom
    for (const auto &fractureName : fractureNames) {
        auto result = model_->CreateSplitBoundaryFrom( fractureName.c_str(), retain_elmts_as_intervening_elements );
        if ( !result.second ) {
            throw std::runtime_error("CreateSplitBoundaryFrom failed for fracture: " + fractureName);
        }
        sb_names.insert( result.first.begin(), result.first.end() );  // accumulate
    }

    //cerr << endl << "created split boundaries: ";
    //for (auto &n : sb_names) cerr << endl << n;
    model_->SplitBoundariesOut();

    // NOT NEEDED because this gets done in CreateSplitBoundaryFrom:
    //model_->Mesh().UpdateConnectivity(); // this would be a costly GLOBAL mesh connectivity update

    ////////////////////////////////////////////////////////////
    //! Back inserting lower-dimensional region in each patch
    constexpr int32_t material_id = 1;
    auto new_regions = model_->InsertLowerDimensionalRegionsIntoSplitBoundaries( material_id );

    set<string> fractures_mid_regions;
    cerr << endl << "CREATED THE MID REGIONS: ";
    for (auto &item : new_regions) {
        cout << endl << item;
        cout.flush(); // if you don't trust cout
        fractures_mid_regions.insert(item.c_str());
    }

    //////////////////////////////////////////////////////
    //! Pull apart each split boundary individually
    for (const auto &sb_name : sb_names) {
        model_->SplitBoundary(sb_name).PullApartSplitBoundary(20.);
    }

    model_->MergeRegions(fractures_mid_regions, "FRACTURES MID REGIONS");
    cerr << endl << "FINISHED SPLIT BOUNDARY INITIALIZATION"<<endl<<endl;

    //! Create initial output file
    VTU_Interface<dim> vtu( *model_ );
    vtu.OutputDataToVTU( string(output_path+output_name_).c_str(), output_props_, model_->Region("Model"), 0);

} // end Test_CreateSplitBoundaryFrom



/*
    No need to delete lower-dim elements or re-create the intervening elements
*/
void SplitBoundaryTest_JK::Test_CreateSplitBoundaryFrom_Simplified()
{
    constexpr std::array fractureNames = {"FRACTURE1", "FRACTURE2"};
    const size_t         n_element_original{ model_->Mesh().Elements() };

    bool retain_elmts_as_intervening_elements = true; // KEEP THE LOWER-DIM REGIONS

    //* Following SplitBoundaryInterface_Test workflow
    ////////////////////////////////////
    set<string> sb_names;
    // Loop through each fracture name and call CreateSplitBoundaryFrom
    for (const auto &fractureName : fractureNames) {
        auto result = model_->CreateSplitBoundaryFrom( fractureName, retain_elmts_as_intervening_elements );
        if ( !result.second ) {
            throw std::runtime_error("CreateSplitBoundaryFrom failed for fracture: " + string(fractureName) );
        }
        sb_names.insert( result.first.begin(), result.first.end() );  // accumulate
    }
    // ARE THE FRACTURES STILL THERE? - yes
    model_->RegionsOut();
    _test( model_->Mesh().Elements() == n_element_original );

    //cerr << endl << "created split boundaries: ";
    //for (auto &n : sb_names) cerr << endl << n;
    // EASIER
    model_->SplitBoundariesOut();

    //model_->Mesh().UpdateConnectivity(); // this would be a costly GLOBAL mesh connectivity update
    // NOT NEEDED because this gets done in CreateSplitBoundaryFrom:
    // NOTE: mesh connectivity is a prerequisite for the correct creation of ModelSubDomain objects

    ////////////////////////////////////////////////////////////
    // THE MID REGIONS ARE ALREADY CONNECTED, NOTHING NEEDS TO BE DONE HERE
  
    // setting material ID if so desired
    constexpr int32_t material_id = 1;
    for ( auto it=model_->RegionsBegin(); it!= model_->RegionsEnd(); ++it )
      if ( (*it).second.Name().contains("FRACTURE") )
     for ( auto& mit : (*it).second.CellVector() )
       mit->Material_ID( material_id );

    //////////////////////////////////////////////////////
    // PULL-APART NOT NEEDED ANYMORE

    model_->MergeRegions( {"FRACTURE1", "FRACTURE2"}, "FRACTURES_MID_REGIONS" );
    cerr << endl << "FINISHED SPLIT BOUNDARY INITIALIZATION"<<endl<<endl;

    //! Create initial output file
    VTU_Interface<dim> vtu( *model_ );
    vtu.OutputDataToVTU( string(output_path+output_name_).c_str(), output_props_, model_->Region("Model"), 0);
    vtu.OutputDataToVTU( string(output_path+output_name_).c_str(), output_props_, model_->Region("FRACTURE1"), 0);
    vtu.OutputDataToVTU( string(output_path+output_name_).c_str(), output_props_, model_->Region("FRACTURE2"), 0);
    
    cout <<"\n"<<"That's it"<< endl;

} // end Test_CreateSplitBoundaryFrom


}//end csmp
