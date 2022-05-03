#include "ANSYS_SplitBoundaryMatch_Test.h"
#include "ErrorHandler.h"

#include "ANSYS_Model3D.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "PL_Utilities.h"
#include "VTK_Interface.h"
#include "variableOperations.h"
#include "CSMP_mathUtilities.h"
#include "MeshManagementUtilities.h"

#include "Model.h"
#include "VTU_Interface.h"

using namespace std;

namespace csmp {

/**
    tests findSplitInterfaceElements() and the creation of a SplitBoundary via DetectAndCreateSplitBoundaries()
    using a node-matched model with internal split boundaries built in ANSYS.
 
    SKM 15/01/18
*/
void ANSYS_SplitBoundaryMatch_Test::run()
  {
      _test( TestForContiguousModel() );
      _test( TestForDiscontiguousModel() );
   
} // end run




/**
      Model 'Dyke_Split' contains several volumetric domains disconnected along node-matched split boundaries.
      These domains are called:
      
      GEOL1_VOL
      GEOL2_VOL
      GEOL3_VOL
      GEOL4_VOL
      DYKE_VOL
*/
bool ANSYS_SplitBoundaryMatch_Test::TestForDiscontiguousModel()
 {
      ErrorHandler& csmp_error{ ErrorHandler::Instance() };
      
      // BINARY IO
      // ---------
      const string model_name("Dyke_Split");
      if ( verbose_ ) {
           cout <<"\nNSYS_SplitBoundaryMatch_Test::run: "<<this->getName()<<endl<<endl;
           cout <<"Building Model: '"<< model_name <<"''"<< endl;
        }
      const string variablesFile("ANSYS_SplitBoundaryMatch_Test-variables.txt");
      ANSYS_Model3D model( model_name.c_str(), model_name.c_str(), variablesFile.c_str(), true, true );
      Region<3U>    model_domain(model.Region("Model"));
      model_domain.UpdateMemberIndexes();
      VTK_Interface<3U>  vtk_output;
      VTU_Interface<3U>  vtu_output( model );
      
      return true;
      
 } // end TestForDiscontiguousModel


    // TESTING WHETHER PRECISION HAS IMPACT ON NODE MATCHING
    /*
        const double n_dec_places(1.0e-3);
        // regularising node positions
        for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
             (*nit)->x( quantiseToScale( (*nit)->x(), n_dec_places ) );
             (*nit)->y( quantiseToScale( (*nit)->y(), n_dec_places ) );
             (*nit)->z( quantiseToScale( (*nit)->z(), n_dec_places ) );
          }
     */
     



/**
      Model Split_Edges contains a half-open box in the middle some of the walls of which are
      ANSYS internal boundaries which have been split while node-matching is retained (NOTE: make sure to enforce this when
      trying to use this functionality).
      
      To avoid self-intersection (for which the   findSplitInterFace() mesh utility function does not work,
      a model subregion other than "Model" has to be specified as input.
*/
bool ANSYS_SplitBoundaryMatch_Test::TestForContiguousModel()
 {
      ErrorHandler& csmp_error{ ErrorHandler::Instance() };
      bool test_was_successful{true};
      
      // BINARY IO
      // ---------
      const string model_name("Split_Edges");
      if ( verbose_ ) {
           cout <<"\nNSYS_SplitBoundaryMatch_Test::run: "<<this->getName()<<endl<<endl;
           cout <<"Building Model: '"<< model_name <<"'"<< endl;
        }
      const string variablesFile("ANSYS_SplitBoundaryMatch_Test-variables.txt");
      //                   fileset             regions-file
      ANSYS_Model3D model( model_name.c_str(), model_name.c_str(), variablesFile.c_str(), true, true );
      Region<3U>    model_domain(model.Region("Model"));
      model_domain.UpdateMemberIndexes();
      VTK_Interface<3U>  vtk_output;
      VTU_Interface<3U>  vtu_output( model );
      
    // 0. assigning unique values to visualise the model regions
    // ---------------------------------------------------------
     double evar(1.);
     cerr <<"\nrun: assigning values to unique regions:\n";
     for ( auto it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); ++it ) {
          cerr <<"\t'"<< (*it).first <<"'";
          (*it).second.InputPropertyValue( "element variable", makeScalar(PLAIN,evar) );
          evar += 1.;
       }
     cerr << endl;
     // painting internal boundary elements (they are flagged neither INTERNAL nor IRREGULAR
     const double color(0.);
     const csmp::Index evar_key(model.Database().StorageKey("element variable"));
     for ( auto eit=model_domain.PerimeterCellsBegin(); eit!=model_domain.CellsEnd(); ++eit )
       (*eit)->Store( evar_key, makeScalar(PLAIN,color) );

     // "nodal variable" is assigned 0 in interior and 2 on perimeter; if perimeter is not
     // flagged correctly a value of 5 is assigned to the perimeter node
     const csmp::Index    nvar_key(model.Database().StorageKey("nodal variable"));
     multiset<Point<3U> > split_nodes;
     model_domain.InputPropertyValue( "nodal variable", makeScalar(PLAIN,0.) );
     for ( auto nit=model_domain.PerimeterNodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
          // external model boundary
          if ( (*nit)->AtBoundary() != NOT && (*nit)->AtBoundary() != INTERNAL )
            (*nit)->Store( nvar_key, makeScalar(PLAIN,5.) );
          // internal (split) boundary
          else {
               (*nit)->Store( nvar_key, makeScalar(PLAIN,2.) );
               split_nodes.insert( (*nit)->Coordinate() );
            }
       }
  
     if ( verbose_ ) {
          vtk_output.OutputDataToVTK( model, "test-evar", "element variable", 0U );
          vtk_output.OutputDataToVTK( model, "test-nvar", "nodal variable", 0U );
          // output of perimeter elements to VTK (whole perimenter gets identified)
          outputRegionBoundaryToVTK( model, "Model", "model-boundary" );
       }
  

    // 1. create element correspondance set
    // ---------------------------------------------------------------------------------------
    // creating region labels and tagging the regions with unique integer indentifiers
    const string region_tag("region identifier");
    if ( !model.Database().IsDefined(region_tag.c_str()) )
      model.CreateProperty( region_tag.c_str(), "X", SCALAR, ELEMENT );
    const csmp::Index mtrl_key = model.Database().StorageKey(region_tag.c_str());
    vector<string>  region_names;
    const size_t model_regions = model.CountAndLabelRegions( region_tag.c_str(), region_names );
    assert( model_regions > 1U );

      // reports local ids of shared face on either side of the split boundary
      // std::pair<std::pair<Element<3U>*,size_t>,std::pair<Element<3U>*,size_t> >
      vector<OppositeElements>  interface_elmt_pairs;
      _test( findSplitInterfaceElements( model, region_tag, interface_elmt_pairs ) );
  
      // changing element property values in the discovered elements
      for ( auto eit=interface_elmt_pairs.begin(); eit!=interface_elmt_pairs.end(); ++eit ) {
           (*eit).first.first->Store( evar_key, makeScalar(PLAIN,6.) );
           (*eit).second.first->Store( evar_key, makeScalar(PLAIN,7.) );
        }
// ERROR: some INSIDE elements appear on the outside
      vtk_output.OutputDataToVTK(model, "high-vals-next-to-interface", "element variable", 1U );


   // 2. build CSMP SplitBoundary
   // ---------------------------------------------------------------------------------------
   model.DetectAndCreateSplitBoundaries();
   // checking which boundaries were created
   model.SplitBoundariesOut();

   // 3. change some property along split boundary to verify that assignments are made correctly
   // ------------------------------------------------------------------------------------------
   if ( model.SplitBoundariesBegin() != model.SplitBoundariesEnd() )
     {
       SplitBoundary<3U>& splitdomain = (*model.SplitBoundariesBegin()).second;
       // node variable "nodal variable"
       if ( verbose_ ) cout <<"\nrun: parameterising region: "<< splitdomain.Name() <<"\n";
       model.InputPropertyValue( "nodal variable", makeScalar(ANY,0.) );
       // inside of boundary: OK
       splitdomain.InputPropertyValue( "nodal variable", makeScalar(ANY,10.), INTERIOR, INSIDE );
       splitdomain.InputPropertyValue( "nodal variable", makeScalar(ANY,11.), PERIMETER, INSIDE );
       // outside of boundary: OK
       splitdomain.InputPropertyValue( "nodal variable", makeScalar(ANY,-10.), INTERIOR, OUTSIDE );
       splitdomain.InputPropertyValue( "nodal variable", makeScalar(ANY,-11.), PERIMETER, OUTSIDE );

       // split boundary variable "interface flux" OK
       const csmp::Index key = model.Database().StorageKey("split boundary flux");
       splitdomain.Store( key, makeScalar(ANY,1.0e-5) );
       _equal( splitdomain.Read(key), 1.0e-5, numeric_limits<double>::epsilon() );

       if ( verbose_ ) {
            // output the split boundary here
            SplitBoundary<3U>& splitBoundary1 = (*model.SplitBoundariesBegin()).second;
            splitBoundary1.Out(); // too many nodes in manifolds!
            TestThatManifoldNodesAreColocated( splitBoundary1 );
            const csmp::Index var_key = model.Database().StorageKey("nodal variable");
            PrintSplitBoundaryNodeVariableVector( var_key, splitBoundary1.NodesBegin(), splitBoundary1.NodesEnd(), INSIDE );
            set<string> var_names{"nodal variable"};
            // TODO: this method is broken
            // vtu_output.OutputDataToVTU( splitBoundary1.Name(), var_names, splitBoundary1, 0L );
         }
     }
   else {
        csmp_error.Note( ERROR, "ANSYS_SplitBoundaryMatch_Test::run", "split boundary could not be formed" );
        return false;
     }
   

   // 4. write altered element properties on either side to VTK
   // ---------------------------------------------------------
   string filename = model_name + "-nodal-variable";
   // writing all volumetric regions to files
   if ( verbose_ ) {
        for ( auto rit=model.UniqueRegionsBegin(); rit!= model.UniqueRegionsEnd(); ++rit ) {
          cout <<"\nrun: saving region: "<< (*rit).first <<"\n";
          pair<int32_t,int32_t> region_shape = (*rit).second.SpatialDimensions();
          if ( region_shape.first == 1 and region_shape.second == 3 )
            vtk_output.OutputDataToVTK( model, (*rit).first, filename, string("nodal variable"), 2U );
          }
     }
  
   // 5. see whether the split boundary survives being writting to and recovered from file
   // ------------------------------------------------------------------------------------
   model.OutputToBinaryFile( model_name.c_str() );
   // reconstructing model from file
   set<string> test_variables{ "nodal variable", "split boundary flux" };
   Model<3U>   split_model( model_name, test_variables );
   // writing all volumetric regions to files
   if ( verbose_ ) {
        filename += "_restored";
        for ( auto rit=split_model.UniqueRegionsBegin(); rit!= split_model.UniqueRegionsEnd(); ++rit ) {
             cout <<"\nrun: saving region: "<< (*rit).first <<"\n";
             pair<int32_t,int32_t> region_shape = (*rit).second.SpatialDimensions();
             if ( region_shape.first == 1 and region_shape.second == 3 )
               vtk_output.OutputDataToVTK( model, (*rit).first, filename, string("nodal variable"), 2U );
          }
     }
     
   return test_was_successful;

 } // end TestForContiguousModel
 
 
 
     /// for testing how the nodes are organised in the node vector
void ANSYS_SplitBoundaryMatch_Test::PrintSplitBoundaryNodeVariableVector( const csmp::Index& var_key,
                                                                          vector<Node<3U>*>::const_iterator begin,
                                                                          vector<Node<3U>*>::const_iterator end,
                                                                          INTERFACE_SIDE iside )
 {
    assert( distance(begin,end) > 0U );
    assert( var_key.place == NODE );
    assert( var_key.type == SCALAR );
    assert( iside != MIDDLE ); // there should be no access to the nodes of another region, other than through the node manifold
    
    if ( iside == INSIDE ) {
         cout <<"\n\nPrintSplitBoundaryNodeVariableVector: "<< distance(begin,end);
         cout <<" INSIDE values of scalar node variable in the order stored in 'node_vec_':\n";
         while( begin != end ) {
              cout << (*begin)->Idx() <<": "<< (*begin)->Read( var_key ) <<", ";
              begin++;
           }
         cout << endl << endl;
         return;
      }
      
    // OUTSIDE
    cout <<"\n\nPrintSplitBoundaryNodeVariableVector: use ManifoldNode(size_t) method to get access to outside.";
    cout << endl << endl;
     
 } // end PrintSplitBoundaryNodeVariableVector




void ANSYS_SplitBoundaryMatch_Test::TestThatManifoldNodesAreColocated( const SplitBoundary<3U>& split_boundary )
 {
    // getting the coordinate range for scaling
    Point<3U> minCoord{1e9,1e9,1e9}, maxCoord{-1e9,-1e9,-1e9};
    for ( size_t i{0U}; i<split_boundary.Nodes(); i++ )
      {
         Point<3U> ref_point = split_boundary.N(i)->Coordinate();
         minCoord = min( minCoord, ref_point );
         maxCoord = max( maxCoord, ref_point );
      }
    const double scale_factor = minCoord.DistanceTo( maxCoord );

    for ( size_t i{0U}; i<split_boundary.Nodes(); i++ )
      {
         // adding all manifold points and then dividing by their number before comparing to inside node location
         Point<3U> ref_point = split_boundary.N(i)->Coordinate();
         minCoord = min( minCoord, ref_point );
         Point<3U> avg_point{0.,0.,0.};
         for ( auto j{0U}; j<split_boundary.ManifoldNode(i)->Branches(); j++ )
           avg_point += split_boundary.ManifoldNode(i)->N(j)->Coordinate();
         // averaging
         avg_point /= static_cast<double>(split_boundary.ManifoldNode(i)->Branches());
         // testing
         if ( ref_point.DistanceTo(avg_point) > numeric_limits<double>::epsilon()*scale_factor )
           cout <<"\nnode "<< split_boundary.N(i)->Idx() <<": offset by: "<< ref_point.DistanceTo(avg_point) << endl;
         _test( ref_point.DistanceTo(avg_point) <= numeric_limits<double>::epsilon()*scale_factor );
      }
      
 } // end TestThatManifoldNodesAreColocated



} // csmp
