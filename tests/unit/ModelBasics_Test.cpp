/*
 *  ModelBasics_Test.h
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *
 */

#include "ModelBasics_Test.h"

#include "vsetMakers.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyConstraints.h"
#include "VSet.h"
#include "VTK_Interface.h"
#include "ANSYS_Model3D.h"
#include "meshManagementUtilities.h"
#include "compareFloats.h"
#include "ModelComparator.h"

using namespace std;

namespace csmp {

/**
    Reading of all basic variable types and placements.
*/
void ModelBasics_Test::run()
  {
     TestModelConstructionFromVSet();
     
     // builds model from regions, converting lower-dimensional elements into faces
     _test( TestWriteModelToDiskAndReadBack() );
     
     // 2D model with split boundaries
     _test( TestWriteModelToDiskAndReadBackWithInterfaces() );
 
   } // end run



bool ModelBasics_Test::TestModelConstructionFromVSet()
 {
    const  bool bSkewed{false};
    VSet<3U>    vset, vset1;
    create_Pyramid_Hexa_VSet( vset, bSkewed ); // no boundaries!
    
    // does the VSet write/reads correctly?
    double time0{3600.123}, time1;
    vset.OutputTo( "ModelBasics_Test", time0 );
    vset1.InputFrom( "ModelBasics_Test", time1 );
    // testing
    _test( vset1 == vset );
    _test( approximatelyEqual(time0,time1) );
    
    // test: basic constructor
    Model<3U>   model( vset, "ModelBasics_Test-variables.txt" );
    // save to native binary
    model.OutputToBinaryFile( "ModelBasics_Test" );
    // bring back from binary
    set<string>  subset_variables; // all variables
    Model<3U>    restored_model( string{"ModelBasics_Test"}, subset_variables );
     
    // compare = test
    return true;
    
 } // end TestModelConstructionFromVSet




bool ModelBasics_Test::TestWriteModelToDiskAndReadBack()
 {
     VSet<3U>      vset;
     ModelTopology topology = create_FracBox( vset );     

     // creating model
     const bool do_not_use_regions_file{true};
     Model<3U>  model( topology, vset, "ModelBasics_Test-variables.txt", do_not_use_regions_file );
     model.Name("FracBox");
     
     // assigning some dummy values to verify functionality
     Region<3U>& model_domain = model.Region("Model");
     Region<3U>& matrix_domain = model.Region("MATRIX");
     model_domain.InputPropertyValue( "fluid pressure", makeScalar(ANY,1e5) );
     model_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-15) );
     matrix_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );
     
     // a random fluid pressure distribution to verify that nodal values are read back correctly
     const size_t n = model_domain.Nodes();
     constexpr double p_min = 100325.0;
     constexpr double p_max = 2.7e7;
     vector<double>   random_pressure(n);

     // Random engine (seeded properly)
     std::random_device rd;
     std::mt19937_64 gen(rd());  // 64-bit Mersenne Twister
     std::uniform_real_distribution<double> dist(p_min, p_max);
     // Fill using C++20 ranges and lambda
     std::ranges::generate( random_pressure, [&] { return dist(gen); } );
     
     // assigning the values to the model, and mapping them to node-coordinates
     map<Point<3>,double>  pressure_map; // (needed because the nodes will not necessarily be in the same order in the new mode)
     csmp::Index pkey = model.Database().StorageKey("fluid pressure");
     size_t count{0};
     for ( auto& nit : model_domain.NodeVector() ) {
           nit->Store( pkey, makeScalar(ANY,random_pressure[count++]) );
           pressure_map.insert( make_pair( nit->Coordinate(), nit->Read(pkey) ) );
       }

     if ( verbose_ ) model.Out();
       
     // saving model to disk
     model.OutputToBinaryFile( "ModelBasics_Test" );
     
     // bringing model back
     set<string>  subset_variables; // all variables
     Model<3U>    restored_model( string{"ModelBasics_Test"}, subset_variables );
     double       pmin, pmax;
     restored_model.MinMaxOf( "permeability", pmin, pmax );
     
     _test( approximatelyEqual( pmin, 1.0e-15 ) );
     _test( approximatelyEqual( pmax, 1.0e-12 ) );
     
     // testing fluid pressure values
     const Region<3>& restored_model_domain = restored_model.Region("Model");
     map<Point<3>,double>  restored_pressure_map;
     bool first_call{true};
     count = 0;
     for ( auto& nit : restored_model_domain.NodeVector() )
       restored_pressure_map.insert( make_pair( nit->Coordinate(), nit->Read(pkey) ) );

     // lexicographical comparison
     _test( pressure_map.size() == restored_pressure_map.size() );
     _test( pressure_map == restored_pressure_map );

     // tolerance based comparison
     // (co-iterating the maps and comparing the pressures)
     {
        auto it1 = pressure_map.begin();
        auto it2 = restored_pressure_map.begin();
        for (; it1 != pressure_map.end() && it2 != restored_pressure_map.end(); ++it1, ++it2) {
            // keys must match
            const auto& [k1, v1] = *it1;
            const auto& [k2, v2] = *it2;
            _test(k1 == k2);
            // comparing values v1 and v2
            _test( approximatelyEqual( v1, v2 ) );
            if ( verbose_ ) {
                 if ( !approximatelyEqual( v1, v2 ) && first_call ) {
                      cout <<"\n\n"<<"pressure mismatch between original and restored model:"<< endl;
                      cout <<" "<< fabs(v1 - v2) <<",";
                      first_call = false;
                   }
                 if ( !approximatelyEqual( v1, v2 ) && !first_call ) cout <<" "<< fabs(v1 - v2) <<",";
                 cout << std::setprecision(6);
             }
         }
     }
      
     // testing fundamental assumption made working with default initialisations of 'size_t'
     uint32_t default_uint = std::numeric_limits<uint32_t>::max();
     _test( std::numeric_limits<size_t>::max() != UINT_MAX ); // false because UINT_MAX is not for size_t
     _test( isUninitialisedInteger( default_uint ) );
     _test( default_uint == UINT_MAX );
     
     // some visual QC, using VTK
     if ( verbose_ ) restored_model.Out();
     
     return true;
 }




// 2D model SPLIT_22 with crossing split boundaries
bool ModelBasics_Test::TestWriteModelToDiskAndReadBackWithInterfaces()
 {
     VSet<2U> vset;
     ModelTopology topology = create_BoundarySplitBoundaryPatch( vset );

     // creating model
     Model<2U>  model( topology, vset, "CSMP-1phase-variables.txt", false );
     model.Name("FracBox");
     
     // assigning some dummy values to verify functionality
     Region<2U>& model_domain = model.Region("Model");
     Region<2U>& matrix_domain = model.Region("lower");
     model_domain.InputPropertyValue( "fluid pressure", makeScalar(ANY,1e5) );
     model_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-15) );
     matrix_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );

     if ( verbose_ ) model.Out(); // crashes when trying to print normals to FV Stencil for ISO_LIN_HEX
       
     // saving model to disk
     model.OutputToBinaryFile( "ModelBasics_Test" );
     
     // bringing model back
     set<string>  subset_variables; // all variables
     Model<2U>    restored_model( string{"ModelBasics_Test"}, subset_variables );
     double       pmin, pmax;
     restored_model.MinMaxOf( "permeability", pmin, pmax );
     
     _test( approximatelyEqual( pmin, 1.0e-15 ) );
     _test( approximatelyEqual( pmax, 1.0e-12 ) );
     
     // testing fundamental assumption made working with default initialisations of 'size_t'
     uint32_t default_uint = std::numeric_limits<uint32_t>::max();
     _test( std::numeric_limits<size_t>::max() != UINT_MAX ); // false because UINT_MAX is not for size_t
     _test( isUninitialisedInteger( default_uint ) );
     _test( default_uint == UINT_MAX );
     
     // some visual QC, using VTK
     if ( verbose_ ) restored_model.Out();
     
     return true;
     
 } // TestWriteModelToDiskAndReadBackWithInterfaces


  
} // end csmp



