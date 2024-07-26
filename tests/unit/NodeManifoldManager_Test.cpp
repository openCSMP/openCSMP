#include "NodeManifoldManager_Test.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "NodeManifoldManager.h"

#include "VTU_Interface.h"


using namespace std;


namespace csmp
{

// TODO: test should not depend on ANSYS_Models
// TODO: run 2 and 3D versions
void NodeManifoldManager_Test::run()
{
  // NB: test is written in such a way that the model must contain a volumetric region "FAULT"
  //     there is currently no 3D model that matches this requirement
  Test_nodemanifolds_created_from_splitboundaries_between_regions<2U>( "box2d_fault" );
  //Test_nodemanifolds_created_from_splitboundaries_between_regions<3U>( "..." );
}




template<uint32_t dim>
void NodeManifoldManager_Test::Test_nodemanifolds_created_from_splitboundaries_between_regions( const std::string& model_name )
{
  std::ostringstream ostr;
  std::string dimension( "" );
  ostr << dim;
  dimension += ostr.str();
  dimension += "D";

  if ( verbose_ ) std::cerr << "\nStart " << dimension << " NodeManifoldManager_Test: Test_nodemanifolds_created_from_splitboundaries_between_regions \n"; 

  // load Model
  const std::string variables_file( "NodeManifoldManager_Test-variables.txt" );

  Model<dim>* modelIN(nullptr);

  const bool with_regions_file{true};
  if ( dim == 2U )
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), true, with_regions_file ));
  else if ( dim == 3U )
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), true ));

  // check the model
  if ( verbose_ ) {
    cout << "\nNodeManifoldManager_Test<" << dim << ">::Test_nodemanifolds_created_from_splitboundaries_around_regions:";
    cout <<"\n\t"<<"model '"<< modelIN->Name();
    for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ )
      cout << "\n\t\t"<<"unique region: " << (*it).first;
    for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ )
      cout << "\n\t\t"<<"non-unique region: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->BoundariesBegin(); it != modelIN->BoundariesEnd(); it++ )
      cout << "\n\t\t"<<"boundary: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->SplitBoundariesBegin(); it != modelIN->SplitBoundariesEnd(); it++ )
      cout << "\n\t\t"<<"split boundary: " << (*it).first;
    cout << endl;
  }  

  // create new element property 'region ID' for the split boundary creation
  if(!modelIN->Database().IsDefined("region ID")) modelIN->CreateProperty( "region ID", "rid", "uint", SCALAR, ELEMENT, 1, 0, 1000);
  modelIN->Region("Model").InputPropertyValue("region ID", makeScalar(PLAIN, 0));
  size_t region_id(1);
  for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      modelIN->Region(rg_name).InputPropertyValue("region ID", makeScalar(PLAIN, region_id));
      cout<<"region ID = "<<region_id<<" assign to unique region: "<<rg_name<<endl;
      region_id ++;
  }   

  VTU_Interface<dim>  vtu(*modelIN);
  vtu.OutputDataToVTU( "initial_region_ID", "region ID",  "Model", 0 );  

  cout<<"before split, region FAULT has "<<modelIN->Region("FAULT").PerimeterNodes()<<" perimeter nodes"<<endl;
  if(!modelIN->Database().IsDefined("node indicator")) modelIN->CreateProperty( "node indicator", "nind", "none", SCALAR, NODE, 1, 0, 1000);
  modelIN->Region("Model").InputPropertyValue("node indicator", makeScalar(PLAIN, 0));
  csmp::INDEX<SCALAR,NODE> key_nd = csmp::INDEX<SCALAR,NODE>( modelIN->Database().StorageKey("node indicator") );
  double min_x(1e10);
  Node<dim>* left_bt_nd(nullptr);
  double min_y(1e10);
  for(auto it = modelIN->Region("FAULT").PerimeterNodesBegin(); it!=modelIN->Region("FAULT").NodesEnd(); it++) {
      (*it)->Store(key_nd, makeScalar((*it)->Status(key_nd), 1));
      double x = (*it)->Coordinate()[0];
      double y = (*it)->Coordinate()[1];
      if(x<=min_x && y<=min_y) {min_x = x; min_y = y; left_bt_nd = (*it);};
  }
  cout<<"before split, left bottom node id = "<<left_bt_nd->Idx()<<", x = "<<min_x<<", y = "<<min_y<<", parent elements = "<<left_bt_nd->Parents()<<", node neighbors = "<<left_bt_nd->Neighbors()<<endl;
  
  if ( verbose_ ) vtu.OutputDataToVTU( "node_indicator", "node indicator",  "Model", 0 );  

  //Create splitboundaries
  Create_splitboundary_between_regions<dim>(*modelIN);

  // check the model again
  if ( verbose_ ) {
    cout << "\nNodeManifoldManager_Test<" << dim << ">::Test_nodemanifolds_created_from_splitboundaries_around_regions:";
    cout <<"\n\t"<<"model '"<< modelIN->Name() <<"' after creation of splitboundries:";
    for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ )
      cout << "\n\t\t"<<"unique region: " << (*it).first;
    for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ )
      cout << "\n\t\t"<<"non-unique region: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->BoundariesBegin(); it != modelIN->BoundariesEnd(); it++ )
      cout << "\n\t\t"<<"boundary: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->SplitBoundariesBegin(); it != modelIN->SplitBoundariesEnd(); it++ )
      cout << "\n\t\t"<<"split boundary: " << (*it).first;
    cout << endl;
  }  

  //assign different entry pressures to different regions
  if ( !modelIN->Database().IsDefined("entry pressure") )
    modelIN->CreateProperty( "entry pressure", "pd", "Pa", SCALAR, ELEMENT, 1, 0 ,50000000);
  modelIN->Region("Model").InputPropertyValue("entry pressure", makeScalar(PLAIN, 0.));
  double entry_pressure = 100.;
  for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      if ( rg_name != "Model" ) {
          modelIN->Region(rg_name).InputPropertyValue("entry pressure", makeScalar(PLAIN, entry_pressure));
          cout<<"entry pressure = "<<entry_pressure<<" assign to non-unique region: "<<rg_name<<endl;
          entry_pressure += 100.;
      }      
  }

  for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      modelIN->Region(rg_name).InputPropertyValue("entry pressure", makeScalar(PLAIN, entry_pressure));
      cout<<"entry pressure = "<<entry_pressure<<" assign to unique region: "<<rg_name<<endl;
      entry_pressure += 100.;
  }      

  if ( verbose_ ) vtu.OutputDataToVTU( "entry_pressure", "entry pressure", "Model", 0 );

  // extra variable 'element indicator'
  if ( !modelIN->Database().IsDefined("element indicator") )
    modelIN->CreateProperty( "element indicator", "eind", "none", SCALAR, ELEMENT, 1, 0 ,50000000 );
  modelIN->Region("Model").InputPropertyValue("element indicator", makeScalar(PLAIN, 0.));
  const csmp::INDEX<SCALAR,ELEMENT> key_element = csmp::INDEX<SCALAR,ELEMENT>( modelIN->Database().StorageKey("element indicator") );
  const csmp::INDEX<SCALAR,ELEMENT> key_pd = csmp::INDEX<SCALAR,ELEMENT>( modelIN->Database().StorageKey("entry pressure") );
  size_t count(0ul);
  // sorting node manifolds by 'entry pressure' value
  for( auto mit = modelIN->Mesh().NodeManifoldsBegin(); mit!=modelIN->Mesh().NodeManifoldsEnd(); mit++, count++ ) {
    (*mit).SortByVariableValue(key_pd);
    for(uint32_t n1(0);n1<(*mit).Branches();n1++){
      auto nd1 = (*mit).N(n1);
      nd1->Store(key_nd, makeScalar(nd1->Status(key_nd), 2));
      double x = nd1->Coordinate()[0];
      double y = nd1->Coordinate()[1];
      if(x==min_x && y==min_y) {
        cout<<"after split, left bottom node found, id = "<<nd1->Idx()<<", parent elements = "<<nd1->Parents()<<", node neighbors = "<<nd1->Neighbors()<<endl;
        cout<<"nd position = "<<nd1->Coordinate()<<endl;
        for(auto e(0);e<nd1->Parents();e++){
          cout<<"  parent element id = "<<nd1->Parent(e)->Idx()<<", position = "<<nd1->Parent(e)->BaryCenter()<<endl;
          if(nd1==left_bt_nd) //old node
              nd1->Parent(e)->Store(key_element, makeScalar(PLAIN, 1));
          //else //slave node
          else //new node
              nd1->Parent(e)->Store(key_element, makeScalar(PLAIN, 2));
        }
        for(auto n(0);n<nd1->Neighbors();n++){
          cout<<"  neighbor node id = "<<nd1->Neighbor(n)->Idx()<<", position = "<<nd1->Neighbor(n)->Coordinate()<<endl;
        }
        size_t parent_regions(0);
        for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ ) {
            auto rg_name = (*it).first;
            if(modelIN->Region(rg_name).Contains(nd1)) {
              cout<<"    this node is contained in region "<<rg_name<<endl;
              parent_regions++;
            }
          }
        cout<<"    this node is within "<<parent_regions<<" regions"<<endl;
      }

      if(count==5) {  
        for( uint32_t e(0);e<nd1->Parents();e++){
          if(n1==0) //master node
              nd1->Parent(e)->Store(key_element, makeScalar(PLAIN, 1));
          else //slave node
              nd1->Parent(e)->Store(key_element, makeScalar(PLAIN, 2));
        }
      }

    }
  }
  if ( verbose_ ) {
      vtu.OutputDataToVTU( "element_indicator", "element indicator",  "Model", 0 );
      vtu.OutputDataToVTU( "node_indicator", "node indicator",  "Model", 0 );
    }

  // TESTING the sorted manifolds
  Test_created_manifolds<dim>( *modelIN );
  
} // end Test_nodemanifolds_created_from_splitboundaries_between_regions



template<uint32_t dim>
void NodeManifoldManager_Test::Test_nodemanifolds_created_from_splitboundaries_around_regions( const std::string& model_name )
{
  std::ostringstream ostr;
  std::string dimension( "" );
  ostr << dim;
  dimension += ostr.str();
  dimension += "D";

  if ( verbose_ ) std::cerr << "\nStart " << dimension << " NodeManifoldManager Test: Test_nodemanifolds_created_from_splitboundaries_around_regions \n"; 

  // load Model
  const std::string variables_file( "NodeManifoldManager_Test-variables.txt" );

  Model<dim>* modelIN(nullptr);

  if ( dim == 2U )
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), true ));
  else if ( dim == 3U )
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), true ));

  // check the model
  if ( verbose_ ) {
    cout << "\nNodeManifoldManager_Test<" << dim << ">::Test_nodemanifolds_created_from_splitboundaries_around_regions:";
    for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ )
      cout << "\n\tunique region: " << (*it).first;      
    for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ )
      cout << "\n\tnon-unique region: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->BoundariesBegin(); it != modelIN->BoundariesEnd(); it++ )
      cout << "\n\tboundary: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->SplitBoundariesBegin(); it != modelIN->SplitBoundariesEnd(); it++ )
      cout << "\n\tsplit boundary: " << (*it).first;
    cout << endl;    
  }  

 
  if(!modelIN->Database().IsDefined("region ID")) modelIN->CreateProperty( "region ID", "none", SCALAR, ELEMENT, 1, 0, 1000); 
  size_t region_id(0);
  for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      modelIN->Region(rg_name).InputPropertyValue("region ID", makeScalar(PLAIN, region_id));
      cout<<"region ID = "<<region_id<<" assign to non-unique region: "<<rg_name<<endl;
      region_id ++;
  }   

  VTU_Interface<dim>  vtu(*modelIN);
  vtu.OutputDataToVTU( "initial_region_ID", "region ID",  "Model", 0 );  
  
  /*
  //assign different entry pressures to different regions
  if(!modelIN->Database().IsDefined("entry pressure")) modelIN->CreateProperty( "entry pressure", "Pa", SCALAR, ELEMENT, 1, 0 ,50000000); 
  double entry_pressure = 100.;   
  for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      if ( rg_name != "Model" ) {
          modelIN->Region(rg_name).InputPropertyValue("entry pressure", makeScalar(PLAIN, entry_pressure));
          cout<<"entry pressure = "<<entry_pressure<<" assign to non-unique region: "<<rg_name<<endl;
          entry_pressure += 100.;
      }      
  }

  for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      modelIN->Region(rg_name).InputPropertyValue("entry pressure", makeScalar(PLAIN, entry_pressure));
      cout<<"entry pressure = "<<entry_pressure<<" assign to unique region: "<<rg_name<<endl;
      entry_pressure += 100.;
  }       
  */

  //Create splitboundaries
  //node manifolds already created by function SplitBoundaryInterface::SplitNodes() 
  //during the process of creating splitboundaries
  std::vector<std::string> interfaces;
  Create_splitboundary_around_regions<dim>( *modelIN, interfaces );

  // check the model again
  if ( verbose_ ) {
    cout << "\nNodeManifoldManager_Test<" << dim << ">::Test_nodemanifolds_created_from_splitboundaries_around_regions:after creation of splitboundries:";
    for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ )
      cout << "\n\tunique region: " << (*it).first;      
    for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ )
      cout << "\n\tnon-unique region: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->BoundariesBegin(); it != modelIN->BoundariesEnd(); it++ )
      cout << "\n\tboundary: " << (*it).first;
    cout << endl;
    for ( auto it = modelIN->SplitBoundariesBegin(); it != modelIN->SplitBoundariesEnd(); it++ )
      cout << "\n\tsplit boundary: " << (*it).first;
    cout << endl;    
  }  

  //assign different entry pressures to different regions
  if(!modelIN->Database().IsDefined("entry pressure")) modelIN->CreateProperty( "entry pressure", "Pa", SCALAR, ELEMENT, 1, 0 ,50000000); 
  modelIN->Region("Model").InputPropertyValue("entry pressure", makeScalar(PLAIN, 0.));
  double entry_pressure = 100.;   
  for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      if ( rg_name != "Model" ) {
          modelIN->Region(rg_name).InputPropertyValue("entry pressure", makeScalar(PLAIN, entry_pressure));
          cout<<"entry pressure = "<<entry_pressure<<" assign to non-unique region: "<<rg_name<<endl;
          entry_pressure += 100.;
      }      
  }

  for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      modelIN->Region(rg_name).InputPropertyValue("entry pressure", makeScalar(PLAIN, entry_pressure));
      cout<<"entry pressure = "<<entry_pressure<<" assign to unique region: "<<rg_name<<endl;
      entry_pressure += 100.;
  }      
  
  region_id = 0;
  for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); it++ ) {
      auto rg_name = (*it).first;
      modelIN->Region(rg_name).InputPropertyValue("region ID", makeScalar(PLAIN, region_id));
      cout<<"region ID = "<<region_id<<" assign to non-unique region: "<<rg_name<<endl;
      region_id ++;
  }   

  vtu.OutputDataToVTU( "new_region_ID", "region ID",  "Model", 0 );  
  vtu.OutputDataToVTU( "entry_pressure", "entry pressure",  "Model", 0 );

  Test_created_manifolds<dim>(*modelIN);

  
}



template<uint32_t dim>
void NodeManifoldManager_Test::Create_splitboundary_around_regions( Model<dim>& model,
                                                       std::vector<std::string>& interfaces )
{
  /*
  // 0. Model initialization
  const std::string variables_file( "NodeManifoldManager_Test-variables.txt" );
  Model<dim>* model(nullptr);
  if ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));
  else if ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));
  */
  
  
  string  spliboundary_regions_file( model.Name() );

  // 1. InterFace sets
  std::set<std::string> interface_basic_sets;
  std::set<std::string> interface_sets;
  InputFromFile( std::string( spliboundary_regions_file + "-disconnected-interface-regions.txt" ).c_str(), interface_basic_sets );

  // 2. splitting input regions if they are discontigouos
  bool discontiguous_regions( false );
  for ( set<string>::const_iterator it = interface_basic_sets.begin(); it != interface_basic_sets.end(); ++it )
    if ( !model.IsContiguous( (*it).c_str() ) ) {
      if ( verbose_ ) cerr << "\n\tNodeManifoldManager_Test::Create_splitboundary_around_regions: discovered discontiguous region: " << (*it);
      discontiguous_regions = true;
    }
  if ( discontiguous_regions ) {
    set<string>  original_region_names;
    for ( typename RegionInterface<dim, Region>::regionIterator it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); ++it )
      if ( (*it).first != "Model" )
        original_region_names.insert( (*it).first.c_str() );

    // partitioning regions without revisiting new partitions that can inserted into region map
    for ( set<string>::const_iterator it = original_region_names.begin(); it != original_region_names.end(); ++it )
      model.PartitionRegionIntoContiguousSubRegions( (*it).c_str() );
  } else {
    if ( verbose_ ) cerr << "\n\tNodeManifoldManager_Test::Create_splitboundary_around_regions: no discontiguous region is discovered: "<<endl;
  }

 if ( discontiguous_regions ) {
  // 2. preparing lower dimensional regions for making SplitBoundaries around
  EstablishContiguousRegionsList( model, interface_basic_sets, interface_sets );
  
  /*
  model.MergeRegions( interface_sets, "interfaces" );

  cout<<"\nbefore remove regions:"<<endl;
  for ( auto it = model.RegionsBegin(); it != model.RegionsEnd(); it++ )
      cout << "\n\tregion: " << (*it).first; 
  for ( auto it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); it++ )
      cout << "\n\tunique region: " << (*it).first;        

  for ( auto name : interface_sets )
    model.RemoveRegion( name.c_str(), false );

  cout<<"\nafter remove regions:"<<endl;
  for ( auto it = model.RegionsBegin(); it != model.RegionsEnd(); it++ )
      cout << "\n\tregion: " << (*it).first; 
  for ( auto it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); it++ )
      cout << "\n\tunique region: " << (*it).first;     

  interfaces.clear();
  interfaces.push_back( "interfaces" );
  OutputToFile( std::string( spliboundary_regions_file + "-connected-interface-regions.txt" ).c_str(), interfaces );
  */

  if ( verbose_ ) cout << "\n\n\nNodeManifoldManager_Test::Create_splitboundary_around_regions: the following interface(s) / interface sets will be considered:\n\n";

  // 3. creating SplitBoundaries
  //for ( std::vector<string>::const_iterator it = interfaces.begin(); it != interfaces.end(); it++ ) {
  for ( std::set<string>::const_iterator it = interface_sets.begin(); it != interface_sets.end(); it++ ) {
       //cerr <<"\nRecode this so that it does the right thing!\n";
        /*
        model.CreateInternalBoundaryFrom( (*it).c_str() );
        //model.CreateInternalBoundaryFrom( (*it).c_str(), false);
        Boundary<dim>& bdry = model.Boundary( (*it).c_str() );
        model.CreateSplitBoundaryFrom( bdry );
        */
        model.CreateSplitBoundaryFrom( (*it).c_str());
    }

  // 4. creating lower-dimensional stand-alone meshes from SplitBoundary objects, and 
  //    insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  set<string> new_regions = model.RegionsFromSplitBoundaries();

 } else {
   // creating SplitBoundaries directly from original interfaces
   for ( std::set<string>::const_iterator it = interface_basic_sets.begin(); it != interface_basic_sets.end(); it++ )
       model.CreateSplitBoundaryFrom( (*it).c_str());
 }
// TODO: build this in
model.UpdateConnectivity();
 
  /*
  const bool remove_dim_minus1_region(false);
  pair<set<string>,bool> boundary_patches = model.CreateInternalBoundaryFrom( "FRACTURE", remove_dim_minus1_region );  
  assert( boundary_patches.second == true );
  assert( boundary_patches.first.size() == 1 );
  const string boundary_name = (*boundary_patches.first.begin());
  Boundary<dim>& fractureBoundary = model.Boundary( boundary_name.c_str() ); 
  pair<string,bool> split_boundary = model.CreateSplitBoundaryFrom( fractureBoundary );
  assert( split_boundary.second == true );
  */

  //model.CreateSplitBoundaryFrom( "FRACS" );

    // create lower-dimensional stand-alone meshes from SplitBoundary objects, and 
  // insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  set<string> new_regions = model.RegionsFromSplitBoundaries();

  // 5. Writing it 
  model.OutputToBinaryFile( model.Name() );

  //csmp::Model<dim> model_out( model.Name() );
}



/// Input name of regions to split
void NodeManifoldManager_Test::InputFromFile( const char* file_name,
                                                 std::set<string>& interface_basic_set )
{
  assert( file_name != NULL );

  ifstream ifs( file_name );
  if ( !ifs.is_open() )
    throw csmp::Exception( ERROR, "NodeManifoldManager_Test::inputFromFile",
                           file_name,
                           "file specifying regions is missing" );
  if ( !interface_basic_set.empty() )
    interface_basic_set.clear();

  // reading header line printing it to screen and swallowing empty line thereafter
  char text[256];
  ifs.getline( text, 256 );
  if ( verbose_ ) cout << "\nNodeManifoldManager_Test::inputFromFile: file header: " << text << endl;
  ifs.getline( text, 256 );

  int n_interfaces( 0 );
  ifs >> n_interfaces;
  assert( n_interfaces > 0 );
  assert( n_interfaces < 10000 );

  string  interface_name;
  for ( int n = 0; n<n_interfaces; n++ ) {
    ifs >> interface_name;
    if ( !interface_name.empty() )
      interface_basic_set.insert( interface_name );
    else
      throw csmp::Exception( ERROR, "NodeManifoldManager_Test::inputFromFile:", "encountered empty region name." );
  }

  ifs.close();

  if ( verbose_ ) cout << "\nNodeManifoldManager_Test::inputFromFile: region names stored in '" << file_name << "' read successfully." << endl;

} // end inputFromFile



/// Input name of regions to split
template<uint32_t dim>
void NodeManifoldManager_Test::EstablishContiguousRegionsList( Model<dim>& model,
                                                                  const std::set<string>& interface_basic_set,
                                                                  std::set<string>& interface_sets )
{
  if ( verbose_ ) cout << "\nNodeManifoldManager_Test::EstablishContiguosRegionsList:" << endl;

  if ( !interface_sets.empty() )
    interface_sets.clear();
  bool first_call( true );
  for ( typename map<string, csmp::Region<dim> >::const_iterator
        git = model.UniqueRegionsBegin(); git != model.UniqueRegionsEnd(); git++ )
  {
    for ( set<string>::const_iterator it = interface_basic_set.begin(); it != interface_basic_set.end(); ++it )
      if ( (*git).first.find( *it ) != string::npos )
      {
        interface_sets.insert( (*git).first );
        if ( first_call ) {
          cout << "\n\tregion(s) incorporated into the input list: ";
          first_call = false;
        }
        cout << (*git).first << " ";
        break;
      }
  }
  if ( verbose_ ) cout << "\n\n";

}



/// Write contiguous regions
void NodeManifoldManager_Test::OutputToFile( const char* file_name,
                                                const std::vector<string>& interfaces )
{
  assert( file_name != NULL );

  ofstream ofs( file_name );
  if ( !ofs.is_open() )
    throw csmp::Exception( ERROR, "NodeManifoldManager_Test::outputToFile",
                           file_name, "file specifying permeability-model input variables is missing" );

  ofs << "'" << file_name << "' interface regions to be included.\n\n";

  ofs << interfaces.size() << " ";
  for ( std::vector<string>::const_iterator it = interfaces.begin(); it != interfaces.end(); ++it )
    ofs << (*it) << " ";

  ofs << "\n";
  ofs.close();

  if ( verbose_ ) cout << "\nNodeManifoldManager_Test::outputToFile: '" << file_name << "' written successfully." << endl;

} // end outputToFile



template<uint32_t dim>
void NodeManifoldManager_Test::Test_created_manifolds( Model<dim>& modelIN )
{
  // 1. test created node manifolds
  size_t n_manifolds = modelIN.Mesh().Nodes();
  cout<<"total node manifolds = "<<n_manifolds<<endl;
  _test(n_manifolds>0);
  
  // 2. test whether each node in a manifold is unique
  size_t count(0);
  for(auto mit = modelIN.Mesh().NodeManifoldsBegin();mit!=modelIN.Mesh().NodeManifoldsEnd();mit++,count++) {
      auto md = (*mit);
      if(md.Branches()>2) cout<<"  manifold "<<count<<":"<<endl;
      if((*mit).Branches()>2) cout<<"    number of collocated nodes = "<<(*mit).Branches()<<endl;
      _test((*mit).Branches()>=2);

      for(size_t n1(0);n1<(*mit).Branches();n1++){
        auto nd1 = (*mit).N(n1);
        if((*mit).Branches()>2) cout<<"    node "<<n1<<": index = "<<nd1->Idx()<<endl;
        for(size_t n2(n1+1);n2<(*mit).Branches();n2++){
          auto nd2 = (*mit).N(n2);
          _test(nd1->Idx()!=nd2->Idx());
        }
      }
    }

  // 3. test whether each node is contained in one manifold only
  const Region<dim>& model_domain(modelIN.Region( "Model" ));
  for( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ ) {
        if ( (*nit)->IsManifold() ) {
           for ( auto mit = modelIN.Mesh().NodeManifoldsBegin(); mit!=modelIN.Mesh().NodeManifoldsEnd(); mit++ )
             // if this is not the same manifold (if both have the same address)
             if ( (*nit)->Manifold() != &(*mit) ) {
                 for ( uint32_t n{0u}; n<(*mit).Branches(); n++ )
                   // we test that the Node has not the same address as nodes contained in other manifolds
                   _test( (*nit) != (*mit).N(n) );
               }
         }
    }

  // 4. testing the sorting function
  //    (nodes were sorted by entry pressure)
  const csmp::INDEX<SCALAR,ELEMENT> key_pd = csmp::INDEX<SCALAR,ELEMENT>( modelIN.Database().StorageKey("entry pressure") );
  
  size_t id{0ul};
  for(auto mit = modelIN.Mesh().NodeManifoldsBegin();mit!=modelIN.Mesh().NodeManifoldsEnd();mit++, id++) {
      (*mit).SortByVariableValue(key_pd);
      if ( verbose_ ) cout<<"\nManifold id = "<<id<<":"<<endl;
      for(size_t n1(0);n1<(*mit).Branches();n1++){
        //double pd1 = (*mit).N(n1)->Parent(0)->Read(key_pd);
        double pd1(0.);
        size_t pd1_count(0);
        for( uint32_t elmt = 0u; elmt < (*mit).N(n1)->Parents(); elmt++) {
            if ( verbose_ ) {
                 if ( (dim==2 && (*mit).N(n1)->Parent(elmt)->IsSurface()) ||
                      (dim==3 && (*mit).N(n1)->Parent(elmt)->IsVolume()) ) {
                    bool surface_element(false), volume_element(false);
                    if ( dim==2 && (*mit).N(n1)->Parent(elmt)->IsSurface() ) {
                        surface_element = true;
                        cout<<"  n1 = "<<n1<<" e = "<< elmt <<" pd = "<<(*mit).N(n1)->Parent(elmt)->Read(key_pd)<<" IsSurfaceElement = "<<surface_element <<endl;
                      }
                    if ( dim==3 && (*mit).N(n1)->Parent(elmt)->IsVolume() ) {
                        volume_element = true;
                        cout<<"  n1 = "<<n1<<"  e = "<< elmt <<" pd = "<<(*mit).N(n1)->Parent(elmt)->Read(key_pd)<<" IsVolumeElement = "<<volume_element <<endl;
                      }
                  }
              }
            if(!isnan((*mit).N(n1)->Parent(elmt)->Read(key_pd))) {
                pd1 += (*mit).N(n1)->Parent(elmt)->Read(key_pd);
                pd1_count ++;
            }
         }
        pd1 /= pd1_count;
        
        if ( verbose_ ) {
             cout<<"  Node id = "<<(*mit).N(n1)->Idx()<<":"<<endl;
             for(uint32_t e1(0);e1<(*mit).N(n1)->Parents();e1++)
               cout<<"    Parent element id = "<<(*mit).N(n1)->Parent(e1)->Idx()<<" pd = "<<(*mit).N(n1)->Parent(e1)->Read(key_pd)<<endl;
             cout<<"pd = "<<pd1;
          }
        for(size_t n2(n1+1);n2<(*mit).Branches();n2++){
            double pd2(0.);
            size_t pd2_count(0);
            for( uint32_t e3 = 0; e3 < (*mit).N(n2)->Parents(); e3++) {
                if(!isnan((*mit).N(n2)->Parent(e3)->Read(key_pd))) {
                    pd2 += (*mit).N(n2)->Parent(e3)->Read(key_pd);
                    pd2_count ++;
                }
            }
          pd2 /= pd2_count;

          if(pd1>pd2) cout<<"pd1 = "<<pd1<<" pd2 = "<<pd2<<" Manifold id = "<<id<<":"<<endl; 

          _test(pd1<=pd2);  
          if(pd1==pd2) cout<<"pd1 = pd2 = "<<pd1<<" Manifold id = "<<id<<":"<<endl;     
          //cout<<", "<<pd2;
        }
        //cout<<endl;
      }
  }
} // end Test_created_manifolds



template<uint32_t dim>
void NodeManifoldManager_Test::Create_splitboundary_between_regions( Model<dim>& modelIN )
{
  std::ostringstream ostr;
  std::string dimension( "" );
  ostr << dim;
  dimension += ostr.str();
  dimension += "D";

  if ( verbose_ ) std::cerr << "\nStart " << dimension << " NodeManifoldManager Test: Create_splitboundary_between_regions\n";

  // create SplitBoundaries
  std::vector<std::string> regions;
  regions.reserve( modelIN.UniqueRegions() );

  const auto  model_dim = modelIN.Region( "Model" ).SpatialDimensions();
  for ( auto it = modelIN.UniqueRegionsBegin(); it != modelIN.UniqueRegionsEnd(); ++it ) {
    const auto  sub_dim = (*it).second.SpatialDimensions();
    if ( sub_dim.second == model_dim.second ) // check whether the highest dimension of the region is equal to the highest dimension of the model
      regions.push_back( (*it).second.Name() );
  }

  // search the existing unique sub-regions
  set<pair<string, string>>	discovered;
  deque<string>	current_regions;
  vector<pair<string, string>>  region_final_pairs;
  for ( size_t i = 0U; i < regions.size(); i++ ) {
    string root = regions[i];
    // starting at the first region
    current_regions.push_back( root );
    while ( !current_regions.empty() ) {
      std::string current_region( *current_regions.begin() );
      set<string> neighbors;
      // for all neighbor sub-regions of the current region
      const csmp::Region<dim>&  gref1( modelIN.Region( current_region ) );
      for ( size_t j = 0; j < regions.size(); j++ ) {
        if ( current_region.compare( regions[j] ) == 0 ) continue;
        const csmp::Region<dim>&  gref2( modelIN.Region( regions[j] ) );
        const size_t  shared_nodes( sharedNodes( gref1, gref2 ) );
        if ( shared_nodes > 1 )
          neighbors.insert( regions[j] );
      }

      for ( auto neighbour_region : neighbors ) {
        // if this neighbor is new one        
        auto new_region = discovered.insert( make_pair( current_region, neighbour_region ) );
        if ( new_region.second ) {
          discovered.insert( make_pair( neighbour_region, current_region ) );
          current_regions.push_back( neighbour_region );
          region_final_pairs.push_back( make_pair( current_region, neighbour_region ) );
        }
      }
      // removing the sub-region from the discovered (but not yet explored) deque
      current_regions.pop_front();
    }
  }

  /*
  // assigin different entry pressures to different regions
  std::set<std::string> regions_to_assign;
  for ( auto it : region_final_pairs ) {
    regions_to_assign.insert(it.first);
    regions_to_assign.insert(it.second);
  }
  if(!modelIN.Database().IsDefined("entry pressure")) modelIN.CreateProperty( "entry pressure", "Pa", SCALAR, ELEMENT, 1, 0 ,50000000);  
  double entry_pressure = 100.;
  for (auto rg : regions_to_assign) {
    modelIN.Region(rg).InputPropertyValue("entry pressure", makeScalar(PLAIN, entry_pressure));
    cout<<"entry pressure = "<<entry_pressure<<" assign to region: "<<rg<<endl;
    entry_pressure += 100.;
  }  
  */

  // CREATION OF SPLITBOUNDARY BETWEEN 2 REGIONS
  // -------------------------------------------
  // do this sequentially according to neighbours, otherwise boundaries are not assgiend properly
  for ( auto it : region_final_pairs ) {// for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
     modelIN.CreateSplitBoundaryBetween( it.first.c_str(), it.second.c_str() );
  }
// TODO: build this in
modelIN.Mesh().UpdateConnectivity();
  
  // create lower-dimensional stand-alone meshes from SplitBoundary objects, and 
  // insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
// TODO: method does not exist:  set<string> new_regions = modelIN.RegionsFromSplitBoundaries();

  // output Model to Binary
  modelIN.OutputToBinaryFile(modelIN.Name());

  // screen output
  modelIN.RegionsOut();
  modelIN.BoundariesOut();
  modelIN.SplitBoundariesOut();
  /*
  csmp::Model<dim> model_out( modelIN.Name() );
  cout << "\nNodes: " << model_out.Mesh().Nodes() << "\n";
  cout << "\nNode Groups: " << model_out.Mesh().NodeGroups() << "\n";
  cout << "\nElements: " << model_out.Mesh().Elements() << "\n";
  cout << "\nElement Groups: " << model_out.Mesh().ElementGroups() << "\n";
  cout << "\nFaces: " << model_out.Mesh().Faces() << "\n";
  cout << "\nFace Groups: " << model_out.Mesh().FaceGroups() << "\n";
  cout << "\nInterfaces: " << model_out.Mesh().InterFaces() << "\n";
  cout << "\nInterface Groups: " << model_out.Mesh().InterFaceGroups() << "\n";  
  */

  return;
}


} // csmp
