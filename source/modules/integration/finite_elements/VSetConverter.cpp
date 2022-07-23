#include "VSetConverter.h"
#include "CSMP_ElementSpecifications.h"
#include "VSet.h"
#include "FiniteElementManager.h"
#include "MJL_Point.h"
#include "MJL_Triangle3D.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "compareFloats.h"
#include "Box.h"

using namespace std; 

namespace csmp {


    /// replaces straight-sided global element types with isoparametric ones
template<uint32_t dim>
void VSetConverter<dim>::ConvertElementTypesToOnesUsingLocalCoordinateSystem( VSet<dim>& vset )
 {
    // if the VSet already consists of an isoparametric element family, no conversion is needed
    if ( vset.IsoparametricElementMesh() ) return;
    
    // looping over the element types converting them
    for ( auto n{0}; n<vset.Vertices(); ++n )
      vset.BFlag( n, CSMP_ElementSpecifications::CSMP_TypeUsingLocalCoordinates( vset.BFlag(n)) );
    
 } // end ConvertElementTypesToOnesUsingLocalCoordinateSystem



// inline function definitions

/**
    Interpolation of variable values on the boundary, assuming that it lies in one
    of the coordinate planes.
*/
template<uint32_t dim>
double VSetConverter<dim>::BoundaryValue( const std::map<size_t,double>& bvals,
                                          size_t nID1, size_t nID2 )
 const
  {
      typename std::map<size_t,double>::const_iterator  bvit1(bvals.find(nID1)),
                                                          bvit2(bvals.find(nID2));
      assert( bvit1 != bvals.end() );
      assert( bvit2 != bvals.end() );

      return ((*bvit1).second + (*bvit2).second) / 2.;

  } // end BoundaryValue
                                         



/**

Takes the plist, px, py, pz etc. data from the VSet<dim> and changes the
order of the interpolation functions of linear triangular elements to
quadratic.

The Vset is modified according to the new specifications and the 
property values are interpolated onto the newly generated points just as the 
node coordinates are. 

In the process the 'plist', 'px', 'py', 'pz', 'bflags' arrays
in the Vdata are changed to adapt to the quadratic element. 

@param vset the initialized VSet with the data described above as needed to carry out
the conversion.  

@section implementation  Implementation 

The numbering of the new nodes starts from previous-nodes+1 and goes to
new-total nodes of the mesh. 

@attention This approach
implies that, if node numbers are used to determine the placement in
global solution matrices, nodes from a single triangle will be separated
by many rows/cols in the global solution matrix. If a solver is used,
which applies no pre-conditioning of the solution matrix, before inverting
the solution matrix, this numbering leads to high storage requirements and 
long convergence times. In this case, an algorithm like the Cuthill-McKhee 
node-numbering graph-tree traversal should be used to improve the matrix 
occupancy and to reduce the number of nodes which are far off the 
diagonal.  
*/
template<uint32_t dim>
void VSetConverter<dim>::ConvertLinearToQuadraticTriangles( VSet<dim>& vset ) 
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToQuadraticTriangles", 
                                     "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != LINEAR_TRIANGLE and vset.ElementType(0U) != ISOPARAMETRIC_LINEAR_TRIANGLE) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToQuadraticTriangles", 
                                     "This method only works for linear triangle elements, vs.", 
                          parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }
     const bool  debug(false);
     map<pair<double,double>,int64_t>  nodeIDs;
 //    typename map<pair<double,double>,int64_t>::iterator  ndit;
     double x, y;

     xmin = xmax = vset.Px(0);
     ymin = ymax = vset.Py(0);

     // 1. mapping already existing node points O.K.
     // ---------------------------------------
     const size_t n_nodes{vset.Vertices()};
     for ( size_t i{0U}; i<n_nodes; i++ )
       {
          x = vset.Px(i);
          y = vset.Py(i);
          if ( x < xmin ) xmin = x;
          if ( x > xmax ) xmax = x;
          if ( y < ymin ) ymin = y;
          if ( y > ymax ) ymax = y;
          nodeIDs[ pair<double,double>(x,y) ] = static_cast<int64_t>(i);
       }
         
     // 2. copying already existing boundary flags and values
     //    after flagging the corner nodes
     // -----------------------------------------------------
     FlagCornerNodes( vset );
     vector<int8_t>  bflags( vset.BFlagsBegin(), vset.BFlagsEnd() );
     
     // 3. looping through plist:
     // -------------------------
     //  - expanding node-ID vectors for each element
     //  - adding new node points to px, py, pz vectors
     //  - assigning boundary flags to new node points
     //  - interpolating properties to new node points
     //    (faceverts stay exactly as they were before)
     // to test whether new node point is already part of the mesh or whether
     // it must be created     
     size_t              nID(nodeIDs.size()); // new node ID tracker
     int8_t              bflag;
     map<size_t,int8_t>  new_bflags;

     for ( auto pit=vset.PlistBegin(); pit!=vset.PlistEnd(); pit++ )
      {
         // The Plist node ID vector is resized and the new node coordinates are entered
         (*pit).reserve(6);

         // node 4 is initialized as middle node of triangle face 1
         x = (vset.Px( (*pit)[0] ) + vset.Px( (*pit)[1] )) / 2.;
         y = (vset.Py( (*pit)[0] ) + vset.Py( (*pit)[1] )) / 2.;
 
         if ( debug ) cout <<"\nmidpoint face 1 (nodes:"<< (*pit)[0] <<","<< (*pit)[1] <<"): "<< x <<", "<< y;
 
         auto test_it = nodeIDs.insert( make_pair(make_pair(x,y),nID) );
         // if the middle-node point already exists in the list
         // the node ID which was found for it inside the map is used
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 4("<< nID <<"): "<< x <<", "<< y << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              new_bflags.insert( make_pair( nID, TestForBoundaryFlags( bflags, (*pit)[0], (*pit)[1] ) ) );
             (*pit).push_back( nID++ );
          }
         // node 5 is initialized as middle node of triangle face 2
         x = (vset.Px( (*pit)[1] ) + vset.Px( (*pit)[2] )) / 2.;
         y = (vset.Py( (*pit)[1] ) + vset.Py( (*pit)[2] )) / 2.;
 
         if ( debug ) cout <<"\nmidpoint face 2 (nodes:"<< (*pit)[1] <<","<< (*pit)[2] <<"): "<< x <<", "<< y;
 
         test_it = nodeIDs.insert( make_pair(make_pair(x,y),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 5("<< nID <<"): "<< x <<", "<< y << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              new_bflags.insert( make_pair( nID, TestForBoundaryFlags( bflags, (*pit)[1], (*pit)[2] ) ) );
             (*pit).push_back( nID++ );
          }
         // node 6 is initialized as middle node of triangle face 3
         x = (vset.Px( (*pit)[2] ) + vset.Px( (*pit)[0] )) / 2.;
         y = (vset.Py( (*pit)[2] ) + vset.Py( (*pit)[0] )) / 2.;
 
         if ( debug ) cout <<"\nmidpoint face 3 (nodes:"<< (*pit)[2] <<","<< (*pit)[0] <<"): "<< x <<", "<< y;
 
         test_it = nodeIDs.insert( make_pair(make_pair(x,y),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 6("<< nID <<"): "<< x <<", "<< y << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              new_bflags.insert( make_pair( nID, TestForBoundaryFlags( bflags, (*pit)[2], (*pit)[0] )) );
             (*pit).push_back( nID++ );
          }
      }
 
    // testing whether the new numbers are O.K.
    if ( debug ) {
         cout <<"\nListing old and new nodes and their coordinates:";
         for ( auto ndit=nodeIDs.begin(); ndit!=nodeIDs.end(); ndit++ )   
           cout <<"\nx,y,id: "<< (*ndit).first.first <<", "<< (*ndit).first.second <<": "<< (*ndit).second;
         cout << endl << endl;
      }
      
    // adding the new unique and ordered bflags to the vector and re-assigning it to VSet
    assert( (*new_bflags.begin()).first == bflags.size() );
    bflags.resize( nID, NOT );
    for ( auto f : new_bflags )
      bflags[ f.first ] = f.second;
       
       
    // 4. Creating new 'px', 'py', and 'bflag' arrays and assigning them to VSet<dim>
    // ------------------------------------------------------------------------------
    deque<double>  px( nodeIDs.size() ),
                   py( nodeIDs.size() ),  // new node-point coordinates
                   pz( nodeIDs.size() );
    
    for ( auto nit=nodeIDs.begin(); nit!=nodeIDs.end(); ++nit )
      {
         px[ (*nit).second ] = (*nit).first.first;
         py[ (*nit).second ] = (*nit).first.second;
         pz[ (*nit).second ] = 0.;
      }  
    vset.AddXYZ( px, py, pz );
    vset.AddBFlags( bflags.begin(), bflags.end() );
    
    // 5. Converting CSMP finite element types 
    // ------------------------------------------------------------------------
    vset.ElementType( 0U, ISOPARAMETRIC_QUADRATIC_TRIANGLE ); 

    // 6. Interpolating Node properties to the new node points if there are any
    // ------------------------------------------------------------------------
    // TODO: check whether this is ever used
    //InterpolateNodeProperties( px.size(), 6, vset );
    
    // 7. Converting triangle node order such that first node is in lower-left corner of bounding box
    // ----------------------------------------------------------------------------------------------
    OrderQuadraticTriangleCoordinateOrigins( vset );

 } // end ConvertLinearToQuadraticTriangles









/// as ConvertLinearToQuadraticTriangles() but for surface elements in 3D space
template<uint32_t dim>
void VSetConverter<dim>::ConvertLinearToQuadraticTriangles3D( VSet<dim>& vset ) 
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( dim != 3U ) {
          cout <<"\nVSetConverter<dim>::ConvertLinearToQuadraticTriangles3D: Use only for 3D models."<< endl;
          throw domain_error("VSetConverter<dim>::ConvertLinearToQuadraticTriangles3D");
       }
     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToQuadraticTriangles3D", 
                                     "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != LINEAR_TRIANGLE3D and vset.ElementType(0U) != ISOPARAMETRIC_LINEAR_TRIANGLE) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToQuadraticTriangles3D", 
                                     "This method only works for linear triangle elements, vs.", 
                          parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }
      
     bool  debug(false);
     map<mjl::Point3D,int64_t>                     nodeIDs;
     typename map<mjl::Point3D,size_t>::iterator  ndit;
     double x, y, z;

     xmin = xmax = vset.Px(0);
     ymin = ymax = vset.Py(0);
     zmin = zmax = vset.Pz(0);

     // 1. mapping already existing node points O.K.
     // ---------------------------------------
     for ( size_t i{0U}; i<vset.Vertices(); i++ )
       {
          x = vset.Px(i);
          y = vset.Py(i);
          z = vset.Pz(i);
          if ( x < xmin ) xmin = x;
          if ( x > xmax ) xmax = x;
          if ( y < ymin ) ymin = y;
          if ( y > ymax ) ymax = y;
          if ( z < zmin ) zmin = z;
          if ( z > zmax ) zmax = z;
          // nodes are numbered from 1...n-nodes
          nodeIDs[ mjl::Point3D(x,y,z) ] = i + 1;
       }
         
     // 2. copying already existing boundary flags and values
     // -----------------------------------------------------
     vector<std::int8_t>  bflags( vset.BFlagsBegin(), vset.BFlagsEnd() );
    
     // 3. looping through plist:
     // -------------------------
     //  - expanding node-ID vectors for each element
     //  - adding new node points to px, py, pz vectors
     //  - assigning boundary flags to new node points
     //  - interpolating properties to new node points
     //    (faceverts stay exactly as they were before)
     // to test whether new node point is already part of the mesh or whether
     // it must be created     
     int64_t   nID(nodeIDs.size()); // new node ID tracker
     int8_t  bflag;

     for ( typename deque<vector<int64_t> >::iterator
           pit=vset.PlistBegin(); pit!=vset.PlistEnd(); pit++ )
      {
         // The Plist node ID vector is resized and the new node coordinates are entered
         (*pit).reserve(6);

         // node 4 is initialized as middle node of triangle face 1
         // -------------------------------------------------------
         x = (vset.Px( (*pit)[0] ) + vset.Px( (*pit)[1] )) / 2.;
         y = (vset.Py( (*pit)[0] ) + vset.Py( (*pit)[1] )) / 2.;
         z = (vset.Pz( (*pit)[0] ) + vset.Pz( (*pit)[1] )) / 2.;
 
         if ( debug ) 
           cout <<"\nmidpoint face 1 (nodes:"<< (*pit)[0] <<","<< (*pit)[1] <<"): "<< x <<", "<< y <<", "<< z;
 
         pair<map<mjl::Point3D,int64_t>::iterator, bool>
           test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         // if the middle-node point already exists in the list
         // the node ID which was found for it inside the map is used
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 4("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              if ( (bflag=TestForBoundaryFlags( bflags, (*pit)[0], (*pit)[1] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
         // node 5 is initialized as middle node of triangle face 2
         // -------------------------------------------------------
         x = (vset.Px( (*pit)[1] ) + vset.Px( (*pit)[2] )) / 2.0;
         y = (vset.Py( (*pit)[1] ) + vset.Py( (*pit)[2] )) / 2.0;
         z = (vset.Pz( (*pit)[1] ) + vset.Pz( (*pit)[2] )) / 2.0;
 
         if ( debug ) 
           cout <<"\nmidpoint face 2 (nodes:"<< (*pit)[1] <<","<< (*pit)[2] <<"): "<< x <<", "<< y <<", "<< z;
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 5("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              if ( (bflag=TestForBoundaryFlags( bflags, (*pit)[1], (*pit)[2] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
         // node 6 is initialized as middle node of triangle face 3
         // -------------------------------------------------------
         x = (vset.Px( (*pit)[2] ) + vset.Px( (*pit)[0] )) / 2.0;
         y = (vset.Py( (*pit)[2] ) + vset.Py( (*pit)[0] )) / 2.0;
         z = (vset.Pz( (*pit)[2] ) + vset.Pz( (*pit)[0] )) / 2.0;
 
         if ( debug ) 
           cout <<"\nmidpoint face 3 (nodes:"<< (*pit)[2] <<","<< (*pit)[0] <<"): "<< x <<", "<< y <<", "<< z;
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 6("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              if ( (bflag=TestForBoundaryFlags( bflags, (*pit)[2], (*pit)[0] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
      }
 
    // testing whether the new numbers are O.K.
    if ( debug ) 
      {
         cout <<"\nListing old and new nodes and their coordinates:";
         for ( auto ndit2=nodeIDs.begin(); ndit2!=nodeIDs.end(); ndit2++ ) {
              cout <<"\nx,y,z,id: "<< (*ndit2).first.X() <<", "<< (*ndit2).first.Y() <<", "<< (*ndit2).first.Z();
              cout  <<": "<< (*ndit2).second;
           }
         cout << endl << endl;
      }
 
       
    // 4. Creating new 'px' and 'py' arrays and assigning them to VSet<dim>
    // -------------------------------------------------------------------
    deque<double>  px( nodeIDs.size() ),
                      py( nodeIDs.size() ),  // new node-point coordinates
                      pz( nodeIDs.size() );
                    
    for ( typename map<mjl::Point3D,int64_t>::const_iterator
          nit=nodeIDs.begin(); nit!=nodeIDs.end(); nit++ )
      {
         px[ (*nit).second ] = (*nit).first.X();
         py[ (*nit).second ] = (*nit).first.Y();
         pz[ (*nit).second ] = (*nit).first.Z();
      }  
    vset.AddXYZ( px, py, pz );
    
    // 5. Converting CSMP finite element types 
    // ------------------------------------------------------------------------
    vset.ElementType( 0U, ISOPARAMETRIC_QUADRATIC_TRIANGLE ); 

    // 6. Interpolating Node properties to the new node points if there are any
    // ------------------------------------------------------------------------
    //InterpolateNodeProperties( px.size(), 6, vset );

    // 7. Converting triangle node order such that first node is in lower-left corner of bounding box
    // ----------------------------------------------------------------------------------------------
    OrderQuadraticTriangleCoordinateOrigins( vset );

 } // end ConvertLinearToQuadraticTriangles3D










/**
 
Takes the plist, px, py, pz etc. data from the VSet<dim> and changes the
the linear triagular element to a 7-noded quadratic barycentric triangle
according to the element specifications 
supplied by the barycentric 7-noded triangular finite element.  

The Vset is modified according to the new specifications and the 
property values are interpolated onto the newly generated points just as the 
node coordinates are. 

In the process the 'plist', 'px', 'py', 'pz', 'bflags' arrays
in the Vdata are changed to adapt to the higher-order element. 

@section implementation  Implementation 

The numbering of the new nodes starts from previous-nodes+1 and goes to
new-total nodes of the mesh. This approach
implies that, if node numbers are used to determine the placement in
global solution matrices, nodes from a single triangle will be separated
by many rows/cols in the global solution matrix. If a solver is used,
which applies no pre-conditioning of the solution matrix, before inverting
the solution matrix, this numbering leads to high storage requirements and 
long convergence times. In this case, an algorithm like the Cuthill-McKhee 
node-numbering graph-tree traversal should be used to improve the matrix 
occupancy and to reduce the number of nodes which are far off the 
diagonal.   */
template<uint32_t dim>
void VSetConverter<dim>::ConvertLinearToBarycentricTriangles( VSet<dim>& vset ) 
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToBarycentricTriangles", 
                                     "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != LINEAR_TRIANGLE and vset.ElementType(0U) != ISOPARAMETRIC_LINEAR_TRIANGLE ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToBarycentricTriangles", 
                                     "This method only works for linear triangle elements, vs.", 
                          parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }

    assert( dim == 2U );

     bool  debug(true);
     map<pair<double,double>,size_t>        nodeIDs;
     typename map<pair<double,double>,size_t>::iterator  ndit;
     double  x, y;
      
     xmin = xmax = vset.Px(0);
     ymin = ymax = vset.Py(0);

     // 1. mapping already existing node points
     // ---------------------------------------
     for ( typename deque<vector<int64_t> >::iterator
           pit=vset.PlistBegin(); pit!=vset.PlistEnd(); pit++ )
       for ( typename vector<int64_t>::iterator 
             it=(*pit).begin(); it!=(*pit).end(); it++ )
         {
            x = vset.Px( (*it) );
            y = vset.Py( (*it) );
            if ( x < xmin ) xmin = x;
            if ( x > xmax ) xmax = x;
            if ( y < ymin ) ymin = y;
            if ( y > ymax ) ymax = y;
            nodeIDs[ make_pair(x,y) ] = *it;
         }
         
     // 2. copying already existing boundary flags and values
     //    after flagging the corner nodes
     // -----------------------------------------------------
     FlagCornerNodes( vset );
     
     vector<std::int8_t>  bflags( vset.BFlagsBegin(), vset.BFlagsEnd() );

     // 3. looping through plist:
     // -------------------------
     //  - expanding node-ID vectors for each element
     //  - adding new node points to px, py, pz vectors
     //  - assigning boundary flags to new node points
     //  - interpolating properties to new node points
     //    (faceverts stay exactly as they were before)
     // to test whether new node point is already part of the mesh or whether
     // it must be created     
     double  cx, cy;
     size_t  nID(nodeIDs.size()); // new node ID tracker
     int8_t  bflag;

     for ( typename deque<vector<int64_t> >::iterator
           pit=vset.PlistBegin(); pit!=vset.PlistEnd(); pit++ )
      {
         // The Plist node ID vector is resized and the new node coordinates are entered
         (*pit).reserve(7);

         // node 4 is initialized as middle node of triangle face 1
         // -------------------------------------------------------
         x = (vset.Px( (*pit)[0] ) + vset.Px( (*pit)[1] )) / 2.;
         y = (vset.Py( (*pit)[0] ) + vset.Py( (*pit)[1] )) / 2.;
 
         if ( debug ) cout <<"\nmidpoint face 1 (nodes:"<< (*pit)[0] <<","<< (*pit)[1] <<"): "<< x <<", "<< y;
 
         pair<map<pair<double,double>,size_t>::iterator, bool>
           test_it = nodeIDs.insert( make_pair(make_pair(x,y),nID) );
         // if the middle-node point already exists in the list
         // the node ID which was found for it inside the map is used
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 4("<< nID <<"): "<< x <<", "<< y << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              if ( (bflag=TestForBoundaryFlags( bflags, (*pit)[0], (*pit)[1] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
         // node 5 is initialized as middle node of triangle face 2
         // -------------------------------------------------------
         x = (vset.Px( (*pit)[1] ) + vset.Px( (*pit)[2] )) / 2.;
         y = (vset.Py( (*pit)[1] ) + vset.Py( (*pit)[2] )) / 2.;
 
         if ( debug ) cout <<"\nmidpoint face 2 (nodes:"<< (*pit)[1] <<","<< (*pit)[2] <<"): "<< x <<", "<< y;
 
         test_it = nodeIDs.insert( make_pair(make_pair(x,y),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 5("<< nID <<"): "<< x <<", "<< y << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              if ( (bflag=TestForBoundaryFlags( bflags, (*pit)[1], (*pit)[2] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
         // node 6 is initialized as middle node of triangle face 3
         // -------------------------------------------------------
         x = (vset.Px( (*pit)[2] ) + vset.Px( (*pit)[0] )) / 2.;
         y = (vset.Py( (*pit)[2] ) + vset.Py( (*pit)[0] )) / 2.;
 
         if ( debug ) cout <<"\nmidpoint face 3 (nodes:"<< (*pit)[2] <<","<< (*pit)[0] <<"): "<< x <<", "<< y;
 
         test_it = nodeIDs.insert( make_pair(make_pair(x,y),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) cout <<"\nNew coordinates, Node 6("<< nID <<"): "<< x <<", "<< y << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              if ( (bflag=TestForBoundaryFlags( bflags, (*pit)[2], (*pit)[0] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }

         // node 7 is initialized as node in the element center
         // ---------------------------------------------------
         cx  = vset.Px( (*pit)[0] );
         cy  = vset.Py( (*pit)[0] );
         cx += vset.Px( (*pit)[1] );
         cy += vset.Py( (*pit)[1] );
         cx += vset.Px( (*pit)[2] );
         cy += vset.Py( (*pit)[2] );
         cx /= 3.;
         cy /= 3.;

         if ( debug ) cout <<"\ncoordinates of central node: "<< cx <<","<< cy << endl;
 
         test_it = nodeIDs.insert( make_pair(make_pair(cx,cy),nID) );
         if ( !test_it.second ) 
           {
              csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToBarycentricTriangles", 
                                   "Bary-centre node could not be added in conversion process");
              cout <<"\nCoordinates of central node which could not be added: "<< cx <<","<< cy << endl;
           }
         else (*pit).push_back( nID++ );
      }
 
    // testing whether the new numbers are O.K.
    if ( debug ) 
      {
         cout <<"\nListing old and new nodes and their coordinates:";
         for ( ndit=nodeIDs.begin(); ndit!=nodeIDs.end(); ndit++ )   
           cout <<"\nx,y,id: "<< (*ndit).first.first <<", "<< (*ndit).first.second <<": "<< (*ndit).second;
         cout << endl << endl;
      }
 
       
    // 4. Creating new 'px' and 'py' arrays and assigning them to VSet<dim>
    // -------------------------------------------------------------------
    deque<double>  px( nodeIDs.size() ),
                     py( nodeIDs.size() ),  // new node-point coordinates
                     pz( nodeIDs.size() );

    for ( typename map<pair<double,double>,size_t>::const_iterator
          nit=nodeIDs.begin(); nit!=nodeIDs.end(); nit++ )
      {
         px[ (*nit).second ] = (*nit).first.first;
         py[ (*nit).second ] = (*nit).first.second;
         pz[ (*nit).second ] = 0.;
      }  
    vset.AddXYZ( px, py, pz );
    
    // 5. Interpolating Node properties to the new node points if there are any
    // ------------------------------------------------------------------------
    //InterpolateNodeProperties( px.size(), 7, vset );
    
    // 6. Converting triangle node order such that first node is in lower-left corner of bounding box
    // ----------------------------------------------------------------------------------------------
    OrderBarycentricQuadraticTriangleCoordinateOrigins( vset );
  
    // 7. Converting CSMP finite element types 
    // ------------------------------------------------------------------------
    vset.ElementType( 0U, ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ); 

 } // end ConvertLinearToBarycentricTriangles












/**

For input data from linear triangular elements the function get the 
boundary flags of the face nodes and determines from it the face on which 
the midpoint lies. An error is reported if one of the nodes does not lie 
on the model boundary. 

At this stage the corner nodes have not been identified yet !.
 
tested: O.K. */
template<uint32_t dim>
int8_t  VSetConverter<dim>::TestForBoundaryFlags( const vector<std::int8_t>& bflags,
                                                  size_t nID1, size_t nID2 ) const
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

      assert( !bflags.empty() );
      assert( nID1 < bflags.size() );
      assert( nID2 < bflags.size() );
      
      // since this is about the insertion of midside nodes
      std::int8_t  flag1, flag2;
      if ( bflags[nID1] < 0 ) flag1 = bflags[nID1];
      else 
      return NOT;

      if ( bflags[nID2] < 0 ) flag2 = bflags[nID2];
      else 
      return NOT;
      
      // TWO DIMENSIONAL CASE
      // side boundaries
      if ( flag1 == TOP_OUTSIDE    && flag2 == TOP_OUTSIDE )    return TOP;
      if ( flag1 == BOTTOM_OUTSIDE && flag2 == BOTTOM_OUTSIDE ) return BOTTOM;
      if ( flag1 == LEFT_OUTSIDE   && flag2 == LEFT_OUTSIDE )   return LEFT;
      if ( flag1 == RIGHT_OUTSIDE  && flag2 == RIGHT_OUTSIDE )  return RIGHT;
      // corner boundaries
      if ( flag1 == TOP_OUTSIDE    && flag2 == CNR_MAX_MAXX )   return TOP;
      if ( flag1 == TOP_OUTSIDE    && flag2 == CNR_MAX_MINXZ )  return TOP;
      if ( flag1 == BOTTOM_OUTSIDE && flag2 == CNR_MIN )        return BOTTOM;
      if ( flag1 == BOTTOM_OUTSIDE && flag2 == CNR_MIN_MAXX )   return BOTTOM;
      if ( flag1 == LEFT_OUTSIDE   && flag2 == CNR_MIN )        return LEFT;
      if ( flag1 == LEFT_OUTSIDE   && flag2 == CNR_MAX_MINXZ )  return LEFT;
      if ( flag1 == RIGHT_OUTSIDE  && flag2 == CNR_MIN_MAXX )   return RIGHT;
      if ( flag1 == RIGHT_OUTSIDE  && flag2 == CNR_MAX_MAXX )   return RIGHT;
      // permutations
      if ( flag2 == TOP_OUTSIDE    && flag1 == CNR_MAX_MAXX )   return TOP;
      if ( flag2 == TOP_OUTSIDE    && flag1 == CNR_MAX_MINXZ )  return TOP;
      if ( flag2 == BOTTOM_OUTSIDE && flag1 == CNR_MIN )        return BOTTOM;
      if ( flag2 == BOTTOM_OUTSIDE && flag1 == CNR_MIN_MAXX )   return BOTTOM;
      if ( flag2 == LEFT_OUTSIDE   && flag1 == CNR_MIN )        return LEFT;
      if ( flag2 == LEFT_OUTSIDE   && flag1 == CNR_MAX_MINXZ )  return LEFT;
      if ( flag2 == RIGHT_OUTSIDE  && flag1 == CNR_MIN_MAXX )   return RIGHT;
      if ( flag2 == RIGHT_OUTSIDE  && flag1 == CNR_MAX_MAXX )   return RIGHT;
      // special cases for boundaries that are only one element long 
      if ( flag1 == CNR_MAX_MAXX && flag2 == CNR_MAX_MINXZ )    return TOP;
      if ( flag1 == CNR_MIN      && flag2 == CNR_MIN_MAXX )     return BOTTOM;
      if ( flag1 == CNR_MIN      && flag2 == CNR_MAX_MINXZ )    return LEFT;
      if ( flag1 == CNR_MIN_MAXX && flag2 == CNR_MAX_MAXX )     return RIGHT;
      if ( flag2 == CNR_MAX_MAXX && flag1 == CNR_MAX_MINXZ )    return TOP;
      if ( flag2 == CNR_MIN      && flag1 == CNR_MIN_MAXX )     return BOTTOM;
      if ( flag2 == CNR_MIN      && flag1 == CNR_MAX_MINXZ )    return LEFT;
      if ( flag2 == CNR_MIN_MAXX && flag1 == CNR_MAX_MAXX )     return RIGHT;
      // irregular boundaries
      if ( flag1 == IRREGULAR_OUTSIDE  && flag2 == IRREGULAR_OUTSIDE ) return IRREGULAR;

if constexpr ( dim == 3 )
  throw csmp::Exception( ERROR, "VSetConverter<dim>::TestForBoundaryFlags", "3D case not implemented yet ");
  
      // degenerate cases that should have been picked up by the boundary flagger
      // beforehand
      csmp_error.Note( WARNING, "VSetConverter<dim>::TestForBoundaryFlags:",
                                  "Unable to parse boundary flags.");
      
      cout <<"\nflags: flag1="<< parseBoundary( static_cast<BOX_BOUNDARY>(flag1) ) <<", flag2="<< parseBoundary( static_cast<BOX_BOUNDARY>(flag2) ) << endl;

      return NOT;
      
  } // end TestForBoundaryFlags
     


// tested: O.K. SKM 4/3/02
template<uint32_t dim>
int8_t  VSetConverter<dim>::BoundaryFlags3D( const vector<std::int8_t>& bflags,
                                             size_t nID1, size_t nID2 ) const
  {
      assert( dim == 3U );
      std::int8_t flag1, flag2;

      assert( !bflags.empty() );
      assert( nID1 < bflags.size() );
      assert( nID2 < bflags.size() );
      
      if ( bflags[nID1] < 0 ) flag1 = bflags[nID1];
      else
      return 0;

      if ( bflags[nID2] < 0 ) flag2 = bflags[nID2];
      else
      return 0;
     
      // if both boundary flags are the same
      if ( flag1 == flag2 )
        {
	       // side boundaries
	       if ( flag1 == TOP_OUTSIDE )     return TOP_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE )  return BOTTOM_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE )    return LEFT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE )   return RIGHT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE )   return FRONT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE )    return RIGHT_OUTSIDE;
	       if ( flag1 == REGION_BOUNDARY ) return REGION_BOUNDARY;
	       // edges
	       if ( flag1 == BACK_BOTTOM )  return BACK_BOTTOM;
	       if ( flag1 == BACK_RIGHT )   return BACK_RIGHT;
	       if ( flag1 == BACK_TOP )     return BACK_TOP;
	       if ( flag1 == BACK_LEFT )    return BACK_LEFT;
	       if ( flag1 == BOTTOM_LEFT )  return BOTTOM_LEFT;
	       if ( flag1 == BOTTOM_RIGHT ) return BOTTOM_RIGHT;
	       if ( flag1 == TOP_RIGHT )    return TOP_RIGHT;
	       if ( flag1 == TOP_LEFT )     return TOP_LEFT;
	       if ( flag1 == FRONT_BOTTOM ) return FRONT_BOTTOM;
	       if ( flag1 == FRONT_RIGHT )  return FRONT_RIGHT;
	       if ( flag1 == FRONT_TOP )    return FRONT_TOP;
	       if ( flag1 == FRONT_LEFT )   return FRONT_LEFT;
	       // corners
	       if ( flag1 == CNR_MIN  )      return CNR_MIN;
	       if ( flag1 == CNR_MAX )       return CNR_MAX;
	       if ( flag1 == CNR_MIN_MAXX )  return CNR_MIN_MAXX;
	       if ( flag1 == CNR_MIN_MAXXZ ) return CNR_MIN_MAXXZ;
	       if ( flag1 == CNR_MIN_MAXZ )  return CNR_MIN_MAXZ;
	       if ( flag1 == CNR_MAX_MINXZ ) return CNR_MAX_MINXZ;
	       if ( flag1 == CNR_MAX_MAXX )  return CNR_MAX_MAXX;
	       if ( flag1 == CNR_MAX_MAXZ )  return CNR_MAX_MAXZ;
        }
      else
        {
	       // TOP
	       if ( flag1 == TOP_OUTSIDE && flag2 == BACK_TOP )      return TOP_OUTSIDE;
	       if ( flag1 == TOP_OUTSIDE && flag2 == TOP_LEFT )      return TOP_OUTSIDE;
	       if ( flag1 == TOP_OUTSIDE && flag2 == TOP_RIGHT )     return TOP_OUTSIDE;
	       if ( flag1 == TOP_OUTSIDE && flag2 == FRONT_TOP )     return TOP_OUTSIDE;
	       if ( flag1 == TOP_OUTSIDE && flag2 == CNR_MAX_MAXX )  return TOP_OUTSIDE;
	       if ( flag1 == TOP_OUTSIDE && flag2 == CNR_MAX_MINXZ ) return TOP_OUTSIDE;
	       if ( flag1 == TOP_OUTSIDE && flag2 == CNR_MAX )       return TOP_OUTSIDE;
	       if ( flag1 == TOP_OUTSIDE && flag2 == CNR_MAX_MAXZ )  return TOP_OUTSIDE;
	       
	       // BOTTOM
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == BACK_BOTTOM )   return BOTTOM_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == BOTTOM_LEFT )   return BOTTOM_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == BOTTOM_RIGHT )  return BOTTOM_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == FRONT_BOTTOM )  return BOTTOM_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == CNR_MIN )       return BOTTOM_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == CNR_MIN_MAXX )  return BOTTOM_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == CNR_MIN_MAXZ )  return BOTTOM_OUTSIDE;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == CNR_MIN_MAXXZ ) return BOTTOM_OUTSIDE;
           
           // LEFT
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == BOTTOM_LEFT )   return LEFT_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == BACK_LEFT )     return LEFT_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == TOP_LEFT )      return LEFT_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == FRONT_LEFT )    return LEFT_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == CNR_MIN )       return LEFT_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == CNR_MAX_MINXZ ) return LEFT_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == CNR_MIN_MAXZ )  return LEFT_OUTSIDE;
	       if ( flag1 == LEFT_OUTSIDE   && flag2 == CNR_MAX_MAXZ )  return LEFT_OUTSIDE;
	       
	       // RIGHT
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == BACK_RIGHT )    return RIGHT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == BOTTOM_RIGHT )  return RIGHT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == TOP_RIGHT )     return RIGHT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == FRONT_RIGHT )   return RIGHT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == CNR_MIN_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == CNR_MAX_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == CNR_MIN_MAXXZ ) return RIGHT_OUTSIDE;
	       if ( flag1 == RIGHT_OUTSIDE  && flag2 == CNR_MAX )       return RIGHT_OUTSIDE;

           // FRONT
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == FRONT_BOTTOM )  return FRONT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == FRONT_RIGHT )   return FRONT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == FRONT_TOP )     return FRONT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == FRONT_LEFT )    return FRONT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == CNR_MIN_MAXZ )  return FRONT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == CNR_MIN_MAXXZ ) return FRONT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == CNR_MAX )       return FRONT_OUTSIDE;
	       if ( flag1 == FRONT_OUTSIDE  && flag2 == CNR_MAX_MAXZ )  return FRONT_OUTSIDE;
	       
	       // BACK
	       if ( flag1 == BACK_OUTSIDE   && flag2 == BACK_BOTTOM )   return RIGHT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE   && flag2 == BACK_RIGHT )    return RIGHT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE   && flag2 == BACK_TOP )      return RIGHT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE   && flag2 == BACK_LEFT )     return RIGHT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE   && flag2 == CNR_MIN )       return RIGHT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE   && flag2 == CNR_MIN_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE   && flag2 == CNR_MAX_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag1 == BACK_OUTSIDE   && flag2 == CNR_MAX_MINXZ ) return RIGHT_OUTSIDE;

           // flag1 flipped with flag2

	       // TOP
	       if ( flag2 == TOP_OUTSIDE && flag1 == BACK_TOP )      return TOP_OUTSIDE;
	       if ( flag2 == TOP_OUTSIDE && flag1 == TOP_LEFT )      return TOP_OUTSIDE;
	       if ( flag2 == TOP_OUTSIDE && flag1 == TOP_RIGHT )     return TOP_OUTSIDE;
	       if ( flag2 == TOP_OUTSIDE && flag1 == FRONT_TOP )     return TOP_OUTSIDE;
	       if ( flag2 == TOP_OUTSIDE && flag1 == CNR_MAX_MAXX )  return TOP_OUTSIDE;
	       if ( flag2 == TOP_OUTSIDE && flag1 == CNR_MAX_MINXZ ) return TOP_OUTSIDE;
	       if ( flag2 == TOP_OUTSIDE && flag1 == CNR_MAX )       return TOP_OUTSIDE;
	       if ( flag2 == TOP_OUTSIDE && flag1 == CNR_MAX_MAXZ )  return TOP_OUTSIDE;
	       
	       // BOTTOM
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == BACK_BOTTOM )   return BOTTOM_OUTSIDE;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == BOTTOM_LEFT )   return BOTTOM_OUTSIDE;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == BOTTOM_RIGHT )  return BOTTOM_OUTSIDE;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == FRONT_BOTTOM )  return BOTTOM_OUTSIDE;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == CNR_MIN )       return BOTTOM_OUTSIDE;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == CNR_MIN_MAXX )  return BOTTOM_OUTSIDE;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == CNR_MIN_MAXZ )  return BOTTOM_OUTSIDE;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == CNR_MIN_MAXXZ ) return BOTTOM_OUTSIDE;
           
           // LEFT
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == BOTTOM_LEFT )   return LEFT_OUTSIDE;
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == BACK_LEFT )     return LEFT_OUTSIDE;
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == TOP_LEFT )      return LEFT_OUTSIDE;
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == FRONT_LEFT )    return LEFT_OUTSIDE;
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == CNR_MIN )       return LEFT_OUTSIDE;
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == CNR_MAX_MINXZ ) return LEFT_OUTSIDE;
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == CNR_MIN_MAXZ )  return LEFT_OUTSIDE;
	       if ( flag2 == LEFT_OUTSIDE   && flag1 == CNR_MAX_MAXZ )  return LEFT_OUTSIDE;
	       
	       // RIGHT
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == BACK_RIGHT )    return RIGHT_OUTSIDE;
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == BOTTOM_RIGHT )  return RIGHT_OUTSIDE;
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == TOP_RIGHT )     return RIGHT_OUTSIDE;
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == FRONT_RIGHT )   return RIGHT_OUTSIDE;
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == CNR_MIN_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == CNR_MAX_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == CNR_MIN_MAXXZ ) return RIGHT_OUTSIDE;
	       if ( flag2 == RIGHT_OUTSIDE  && flag1 == CNR_MAX )       return RIGHT_OUTSIDE;

           // FRONT
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == FRONT_BOTTOM )  return FRONT_OUTSIDE;
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == FRONT_RIGHT )   return FRONT_OUTSIDE;
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == FRONT_TOP )     return FRONT_OUTSIDE;
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == FRONT_LEFT )    return FRONT_OUTSIDE;
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == CNR_MIN_MAXZ )  return FRONT_OUTSIDE;
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == CNR_MIN_MAXXZ ) return FRONT_OUTSIDE;
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == CNR_MAX )       return FRONT_OUTSIDE;
	       if ( flag2 == FRONT_OUTSIDE  && flag1 == CNR_MAX_MAXZ )  return FRONT_OUTSIDE;
	       
	       // BACK
	       if ( flag2 == BACK_OUTSIDE   && flag1 == BACK_BOTTOM )   return RIGHT_OUTSIDE;
	       if ( flag2 == BACK_OUTSIDE   && flag1 == BACK_RIGHT )    return RIGHT_OUTSIDE;
	       if ( flag2 == BACK_OUTSIDE   && flag1 == BACK_TOP )      return RIGHT_OUTSIDE;
	       if ( flag2 == BACK_OUTSIDE   && flag1 == BACK_LEFT )     return RIGHT_OUTSIDE;
	       if ( flag2 == BACK_OUTSIDE   && flag1 == CNR_MIN )       return RIGHT_OUTSIDE;
	       if ( flag2 == BACK_OUTSIDE   && flag1 == CNR_MIN_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag2 == BACK_OUTSIDE   && flag1 == CNR_MAX_MAXX )  return RIGHT_OUTSIDE;
	       if ( flag2 == BACK_OUTSIDE   && flag1 == CNR_MAX_MINXZ ) return RIGHT_OUTSIDE;
	       
	       // If unfortunate elements across 2 boundaries, the midside nodes on
	       // such edges must not lie on the boundary !
	       if ( flag1 == BACK_OUTSIDE && flag2 == BOTTOM_OUTSIDE )  return NOT;
	       if ( flag1 == BACK_OUTSIDE && flag2 == RIGHT_OUTSIDE)    return NOT;
	       if ( flag1 == BACK_OUTSIDE && flag2 == TOP_OUTSIDE )     return NOT;
	       if ( flag1 == BACK_OUTSIDE && flag2 == LEFT_OUTSIDE )    return NOT;

	       if ( flag1 == FRONT_OUTSIDE && flag2 == BOTTOM_OUTSIDE ) return NOT;
	       if ( flag1 == FRONT_OUTSIDE && flag2 == RIGHT_OUTSIDE )  return NOT;
	       if ( flag1 == FRONT_OUTSIDE && flag2 == TOP_OUTSIDE )    return NOT;
	       if ( flag1 == FRONT_OUTSIDE && flag2 == LEFT_OUTSIDE )   return NOT;

	       if ( flag1 == TOP_OUTSIDE && flag2 == BACK_OUTSIDE )     return NOT;
	       if ( flag1 == TOP_OUTSIDE && flag2 == LEFT_OUTSIDE )     return NOT;
	       if ( flag1 == TOP_OUTSIDE && flag2 == RIGHT_OUTSIDE)     return NOT;
	       if ( flag1 == TOP_OUTSIDE && flag2 == FRONT_OUTSIDE )    return NOT;

	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == BACK_OUTSIDE )  return NOT;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == LEFT_OUTSIDE )  return NOT;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == RIGHT_OUTSIDE)  return NOT;
	       if ( flag1 == BOTTOM_OUTSIDE && flag2 == FRONT_OUTSIDE ) return NOT;
	       
	       if ( flag1 == LEFT_OUTSIDE && flag2 == BOTTOM_OUTSIDE )  return NOT;
	       if ( flag1 == LEFT_OUTSIDE && flag2 == BACK_OUTSIDE )    return NOT;
	       if ( flag1 == LEFT_OUTSIDE && flag2 == TOP_OUTSIDE )     return NOT;
	       if ( flag1 == LEFT_OUTSIDE && flag2 == FRONT_OUTSIDE )   return NOT;
	       // flag2 = flag1
	       if ( flag2 == BACK_OUTSIDE && flag1 == BOTTOM_OUTSIDE )  return NOT;
	       if ( flag2 == BACK_OUTSIDE && flag1 == RIGHT_OUTSIDE)    return NOT;
	       if ( flag2 == BACK_OUTSIDE && flag1 == TOP_OUTSIDE )     return NOT;
	       if ( flag2 == BACK_OUTSIDE && flag1 == LEFT_OUTSIDE )    return NOT;

	       if ( flag2 == FRONT_OUTSIDE && flag1 == BOTTOM_OUTSIDE ) return NOT;
	       if ( flag2 == FRONT_OUTSIDE && flag1 == RIGHT_OUTSIDE )  return NOT;
	       if ( flag2 == FRONT_OUTSIDE && flag1 == TOP_OUTSIDE )    return NOT;
	       if ( flag2 == FRONT_OUTSIDE && flag1 == LEFT_OUTSIDE )   return NOT;

	       if ( flag2 == TOP_OUTSIDE && flag1 == BACK_OUTSIDE )     return NOT;
	       if ( flag2 == TOP_OUTSIDE && flag1 == LEFT_OUTSIDE )     return NOT;
	       if ( flag2 == TOP_OUTSIDE && flag1 == RIGHT_OUTSIDE)     return NOT;
	       if ( flag2 == TOP_OUTSIDE && flag1 == FRONT_OUTSIDE )    return NOT;

	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == BACK_OUTSIDE )  return NOT;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == LEFT_OUTSIDE )  return NOT;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == RIGHT_OUTSIDE)  return NOT;
	       if ( flag2 == BOTTOM_OUTSIDE && flag1 == FRONT_OUTSIDE ) return NOT;
	       
	       if ( flag2 == LEFT_OUTSIDE && flag1 == BOTTOM_OUTSIDE )  return NOT;
	       if ( flag2 == LEFT_OUTSIDE && flag1 == BACK_OUTSIDE )    return NOT;
	       if ( flag2 == LEFT_OUTSIDE && flag1 == TOP_OUTSIDE )     return NOT;
	       if ( flag2 == LEFT_OUTSIDE && flag1 == FRONT_OUTSIDE )   return NOT;

           // edge and corner cases
	       if ( flag1 == BACK_BOTTOM && flag2 == CNR_MIN )        return BACK_BOTTOM;
	       if ( flag1 == BACK_BOTTOM && flag2 == CNR_MIN_MAXX )   return BACK_BOTTOM;
	       if ( flag1 == BACK_RIGHT && flag2 == CNR_MIN_MAXX )    return BACK_RIGHT;
	       if ( flag1 == BACK_RIGHT && flag2 == CNR_MAX_MAXX )    return BACK_RIGHT;
	       if ( flag1 == BACK_TOP && flag2 == CNR_MAX_MAXX )      return BACK_TOP;
	       if ( flag1 == BACK_TOP && flag2 == CNR_MAX_MINXZ )     return BACK_TOP;
	       if ( flag1 == BACK_LEFT && flag2 == CNR_MAX_MINXZ )    return BACK_LEFT;
	       if ( flag1 == BACK_LEFT && flag2 == CNR_MIN )          return BACK_LEFT;
	       if ( flag1 == BOTTOM_LEFT && flag2 == CNR_MIN )        return BOTTOM_LEFT;
	       if ( flag1 == BOTTOM_LEFT && flag2 == CNR_MIN_MAXZ )   return BOTTOM_LEFT;
	       if ( flag1 == BOTTOM_RIGHT && flag2 == CNR_MIN_MAXX )  return BOTTOM_RIGHT;
	       if ( flag1 == BOTTOM_RIGHT && flag2 == CNR_MIN_MAXXZ ) return BOTTOM_RIGHT;
	       if ( flag1 == TOP_RIGHT && flag2 == CNR_MAX_MAXX )     return TOP_RIGHT;
	       if ( flag1 == TOP_RIGHT && flag2 == CNR_MAX )          return TOP_RIGHT;
	       if ( flag1 == TOP_LEFT && flag2 == CNR_MAX_MINXZ )     return TOP_LEFT;
	       if ( flag1 == TOP_LEFT && flag2 == CNR_MAX_MAXZ )      return TOP_LEFT;
	       if ( flag1 == FRONT_BOTTOM && flag2 == CNR_MIN_MAXZ )  return FRONT_BOTTOM;
	       if ( flag1 == FRONT_BOTTOM && flag2 == CNR_MIN_MAXXZ ) return FRONT_BOTTOM;
	       if ( flag1 == FRONT_RIGHT && flag2 == CNR_MIN_MAXXZ )  return FRONT_RIGHT;
	       if ( flag1 == FRONT_RIGHT && flag2 == CNR_MAX )        return FRONT_RIGHT;
	       if ( flag1 == FRONT_TOP && flag2 == CNR_MAX )          return FRONT_TOP;
	       if ( flag1 == FRONT_TOP && flag2 == CNR_MAX_MAXZ )     return FRONT_TOP;
	       if ( flag1 == FRONT_LEFT && flag2 == CNR_MAX_MAXZ )    return FRONT_LEFT;
	       if ( flag1 == FRONT_LEFT && flag2 == CNR_MIN_MAXZ )    return FRONT_LEFT;
	       
	       if ( flag2 == BACK_BOTTOM && flag1 == CNR_MIN )        return BACK_BOTTOM;
	       if ( flag2 == BACK_BOTTOM && flag1 == CNR_MIN_MAXX )   return BACK_BOTTOM;
	       if ( flag2 == BACK_RIGHT && flag1 == CNR_MIN_MAXX )    return BACK_RIGHT;
	       if ( flag2 == BACK_RIGHT && flag1 == CNR_MAX_MAXX )    return BACK_RIGHT;
	       if ( flag2 == BACK_TOP && flag1 == CNR_MAX_MAXX )      return BACK_TOP;
	       if ( flag2 == BACK_TOP && flag1 == CNR_MAX_MINXZ )     return BACK_TOP;
	       if ( flag2 == BACK_LEFT && flag1 == CNR_MAX_MINXZ )    return BACK_LEFT;
	       if ( flag2 == BACK_LEFT && flag1 == CNR_MIN )          return BACK_LEFT;
	       if ( flag2 == BOTTOM_LEFT && flag1 == CNR_MIN )        return BOTTOM_LEFT;
	       if ( flag2 == BOTTOM_LEFT && flag1 == CNR_MIN_MAXZ )   return BOTTOM_LEFT;
	       if ( flag2 == BOTTOM_RIGHT && flag1 == CNR_MIN_MAXX )  return BOTTOM_RIGHT;
	       if ( flag2 == BOTTOM_RIGHT && flag1 == CNR_MIN_MAXXZ ) return BOTTOM_RIGHT;
	       if ( flag2 == TOP_RIGHT && flag1 == CNR_MAX_MAXX )     return TOP_RIGHT;
	       if ( flag2 == TOP_RIGHT && flag1 == CNR_MAX )          return TOP_RIGHT;
	       if ( flag2 == TOP_LEFT && flag1 == CNR_MAX_MINXZ )     return TOP_LEFT;
	       if ( flag2 == TOP_LEFT && flag1 == CNR_MAX_MAXZ )      return TOP_LEFT;
	       if ( flag2 == FRONT_BOTTOM && flag1 == CNR_MIN_MAXZ )  return FRONT_BOTTOM;
	       if ( flag2 == FRONT_BOTTOM && flag1 == CNR_MIN_MAXXZ ) return FRONT_BOTTOM;
	       if ( flag2 == FRONT_RIGHT && flag1 == CNR_MIN_MAXXZ )  return FRONT_RIGHT;
	       if ( flag2 == FRONT_RIGHT && flag1 == CNR_MAX )        return FRONT_RIGHT;
	       if ( flag2 == FRONT_TOP && flag1 == CNR_MAX )          return FRONT_TOP;
	       if ( flag2 == FRONT_TOP && flag1 == CNR_MAX_MAXZ )     return FRONT_TOP;
	       if ( flag2 == FRONT_LEFT && flag1 == CNR_MAX_MAXZ )    return FRONT_LEFT;
	       if ( flag2 == FRONT_LEFT && flag1 == CNR_MIN_MAXZ )    return FRONT_LEFT;
        }
      // special cases for boundaries that are only one element long 
      if ( flag1 == IRREGULAR_OUTSIDE  && flag2 == IRREGULAR_OUTSIDE ) return IRREGULAR_OUTSIDE;

      // degenerate cases that should have been picked up by the boundary flagger
      // beforehand
      cout <<"\nINFO, VSetConverter<dim>::BoundaryFlags3D: Unable to parse boundary flags"<< endl;
      cout <<"\nflags: flag1="<< flag1 <<", flag2="<< flag2 << endl;

      return NOT;
      
  } // end BoundaryFlags3D




template<uint32_t dim>
int8_t  VSetConverter<dim>::TestForBoundaryFlags3D( double x, double y, double z ) const
  {
     assert( dim == 3U );
  
     // Sides and Edges spanned by Corners 1, 2, 3, 4
     if ( z == zmin )
       {
          // BACK  
          if ( y > ymin && y < ymax && x > xmin && x < xmax ) return BACK_OUTSIDE;
          // EDGE1
          if ( y == ymin && x > xmin && x < xmax ) return BACK_BOTTOM;
          // EDGE2
          if ( x == xmax && y > ymin && y < ymax ) return BACK_RIGHT;
          // EDGE3
          if ( y == ymax && x > xmin && x < xmax ) return BACK_TOP;
          // EDGE4
          if ( x == xmin && y > ymin && y < ymax ) return BACK_LEFT;
       }
       
     // Sides and Edges spanned by Corners 5, 6, 7, 8
     if ( z == zmax )
       {
          // FRONT 
          if ( y > ymin && y < ymax && x > xmin && x < xmax ) return FRONT_OUTSIDE;
          // EDGE9
          if ( y == ymin && x > xmin && x < xmax ) return FRONT_BOTTOM;
          // EDGE10
          if ( x == xmax && y > ymin && y < ymax ) return FRONT_RIGHT;
          // EDGE11
          if ( y == ymax && x > xmin && x < xmax ) return FRONT_TOP;
          // EDGE12
          if ( x == xmin && y > ymin && y < ymax ) return FRONT_LEFT;
       }
     
     // Sides and Edges spanned by Corners 1, 4, 5, 8
     if ( x == xmin )
       {
          // LEFT 
          if ( y > ymin && y < ymax && z > zmin && z < zmax ) return LEFT_OUTSIDE;
          // EDGE5
          if ( y == ymin && z > zmin && z < zmax ) return BOTTOM_LEFT;
          // EDGE8
          if ( y == ymax && z > zmin && z < zmax ) return TOP_LEFT;
       } 
        
     // Sides and Edges spanned by Corners 2, 3, 6, 7
     if ( x == xmax )
       {
          // RIGHT 
          if ( y > ymin && y < ymax && z > zmin && z < zmax ) return RIGHT_OUTSIDE;
          // EDGE6
          if ( y == ymin && z > zmin && z < zmax ) return BOTTOM_RIGHT;
          // EDGE7
          if ( y == ymax && z > zmin && z < zmax ) return TOP_RIGHT;
       } 
     
     // The remaining two Sides (Bottom and Top)
     
     // BOTTOM
     if ( y == ymin && x > xmin && x < xmax && z > zmin && z < zmax ) return BOTTOM_OUTSIDE;
     // TOP
     if ( y == ymax && x > xmin && x < xmax && z > zmin && z < zmax ) return TOP_OUTSIDE;
        
     // nodes to which no boundaries could be assigned
     return 0;
  
  } // end TestForBoundaryFlags3D						 













// tested: O.K.
template<uint32_t dim>
void VSetConverter<dim>::FlagCornerNodes( VSet<dim>& vset, bool three_dimensional ) const
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     // finding min-max of x, and y coordinates of boundary nodes
     map<size_t,int8_t>  bflags;
     size_t              flagging_count(0U);
     
     // 2D CASE finding corners nodes on the basis of their coordinates
     // and flagging them accordingly
     if ( !three_dimensional ) 
       {
         size_t n_node(0U);
         for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
           {
              if ( (*bit) < INTERNAL ) {
                   csmp_error.Note( WARNING, "VSetConverter<dim>::FlagCornerNodes (2D case)",
                                     "Skipping node, the BOX_boundary flag of which could not be identified" );
                   cerr <<"\nNode "<< n_node <<", flagged: "<< parseBoundary( static_cast<BOX_BOUNDARY>(*bit) ) << endl;
                }
              else if ( (*bit) != NOT )
                {
                   // CNR1
                   if ( approximatelyEqual( vset.Px( n_node ), xmin ) && approximatelyEqual( vset.Py( n_node ), ymin ) )
                     { (*bit) = CNR_MIN; flagging_count++; }
                   // CNR2
                   else if ( approximatelyEqual( vset.Px( n_node ), xmax ) && approximatelyEqual( vset.Py( n_node ), ymin ) )
                     { (*bit) = CNR_MIN_MAXX; flagging_count++; }
                   // CNR3
                   else if ( approximatelyEqual( vset.Px( n_node ), xmax ) && approximatelyEqual( vset.Py( n_node ), ymax ) )
                     { (*bit) = CNR_MAX_MAXX; flagging_count++; }
                   // CNR4
                   else if ( approximatelyEqual( vset.Px( n_node ), xmin ) && approximatelyEqual( vset.Py( n_node ), ymax ) )
                     { (*bit) = CNR_MAX_MINXZ; flagging_count++; }
                }
           }
         if ( flagging_count < 4U ) {
              csmp_error.Note( WARNING, "VSetConverter<dim>::FlagCornerNodes (2D case)",
                                          "Less than 4 corner nodes could be flagged" );
              cout <<"\nNumber of flagged nodes: "<< flagging_count << endl;
           }
         return;
       } // end 2D case

     // 3D CASE finding corners nodes on the basis of their coordinates
     // and flagging them accordingly
     const double tol(5.0e-3); // 5 mm to deal with potential imprecision of ANSYS Tetra

     size_t n_node(0U);
     for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
       {
          if ( (*bit) < INTERNAL ) {
               csmp_error.Note( WARNING, "VSetConverter<dim>::FlagCornerNodes (3D case)",
                                         "Skipping node, the boundary flag of which could not be identified" );
               cout <<"\nNode "<< n_node <<", flagged: "<< (*bit) << endl;
            }
          else if ( (*bit) != NOT )
            {
               if ( approximatelyEqual( vset.Pz( n_node ), zmin, tol ) )
                 { 
                    // CNR1
                    if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                      { (*bit) = CNR_MIN; flagging_count++; }
                    // CNR2
                    else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                      { (*bit) = CNR_MIN_MAXX; flagging_count++; }
                    // CNR3
                    else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                      { (*bit) = CNR_MAX_MAXX; flagging_count++; }
                    // CNR4
                    else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                      { (*bit) = CNR_MAX_MINXZ; flagging_count++; }
                 }
               else if ( approximatelyEqual( vset.Pz( n_node ), zmax, tol ) )
                 {
                    // CNR5
                    if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                      { (*bit) = CNR_MIN_MAXZ; flagging_count++; }
                    // CNR6
                    else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                      { (*bit) = CNR_MIN_MAXXZ; flagging_count++; }
                    // CNR7
                    else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                      { (*bit) = CNR_MAX; flagging_count++; }
                    // CNR8
                    else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                      { (*bit) = CNR_MAX_MAXZ; flagging_count++; }
                 }
            }
       }     
     if ( flagging_count < 8U ) {
          csmp_error.Note( WARNING, "VSetConverter<dim>::FlagCornerNodes (3D case)",
                                     "Less than 8 nodes were identified as model corners" );
          cout <<"\nNumber of flagged nodes: "<< flagging_count << endl;
       }
 
 } // end FlagCornerNodes





template<uint32_t dim>
void VSetConverter<dim>::OrderQuadraticTriangleCoordinateOrigins( VSet<dim>& vset ) const
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::OrderQuadraticTriangleCoordinateOrigins", 
                                         "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != ISOPARAMETRIC_QUADRATIC_TRIANGLE ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::OrderQuadraticTriangleCoordinateOrigins", 
                            "This method only works for quadratic triangle elements, vs.", 
                             parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }

    // 1. while looping over all triangles build ordered triangles using the first three
    //    nodes
    typename deque<vector<int64_t> >::iterator  pit(vset.PlistBegin());
    vector<size_t>    pdata(6U);
    typename deque<vector<int64_t> >::iterator  fit(vset.PfvertsBegin());
    vector<int64_t>   pfvert(3U);
    map<mjl::Point,size_t>  ordered_nodes;

    while( pit!=vset.PlistEnd()  and  fit!=vset.PfvertsEnd() )
      {
         // finding node of triangle origin
         ordered_nodes.clear();
         for ( auto i{0U}; i<3U; i++ )
           ordered_nodes[ mjl::Point(vset.Px((*pit)[i]), vset.Py((*pit)[i])) ] = i;
         size_t  offset = (*ordered_nodes.begin()).second;
                                                   
         // resetting the plist entry, remembering that the nodes were ordered
         // counterclockwise initially (only needed if the first node in the 
         // triangle is not the first node in the plist entry)
         if ( offset != 0U )
           { 
              if ( offset == 1U )
                {
                   pdata[0] = (*pit)[1];
                   pdata[1] = (*pit)[2];
                   pdata[2] = (*pit)[0];
                   pdata[3] = (*pit)[4];
                   pdata[4] = (*pit)[5];
                   pdata[5] = (*pit)[3];
                   pfvert[0] = (*fit)[1];
                   pfvert[1] = (*fit)[2];
                   pfvert[2] = (*fit)[0];
                } 
              else if ( offset == 2U )
                {
                   pdata[0] = (*pit)[2];
                   pdata[1] = (*pit)[0];
                   pdata[2] = (*pit)[1];
                   pdata[3] = (*pit)[5];
                   pdata[4] = (*pit)[3];
                   pdata[5] = (*pit)[4];
                   pfvert[0] = (*fit)[2];
                   pfvert[1] = (*fit)[0];
                   pfvert[2] = (*fit)[1];
                } 
         
              // reassigning the new node list to the plist
              for ( auto i{0U}; i<(*pit).size(); i++ ) (*pit)[i] = pdata[i];
         
              // reorganizing the neighbor element list 'pfverts' as well
              for ( auto i{0U}; i<(*fit).size(); i++ ) (*fit)[i] = pfvert[i];
           }
         pit++;
         fit++;
      }
      
 } // end OrderQuadraticTriangleCoordinateOrigins





/**
     Finds edges in boxed shaped model and gibes them a BOX_BOUDARY_FLAG
     dependent on their location.
*/
template<uint32_t dim>
void VSetConverter<dim>::FlagEdges( VSet<dim>& vset, double tol ) const
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::FlagEgdes", 
                                     "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != ISOPARAMETRIC_QUADRATIC_TRIANGLE) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::FlagEgdes", 
                                        "This method only works for linear tetrahedral elements, vs.", 
                          parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }
     // finding min-max of x, and y coordinates of boundary nodes
     map<size_t,double>  bflags;

     // 3D CASE finding the edges on the basis of their coordinates
     // and flagging them accordingly. note there is no 2D case since in 2D we
     // are splitting up triangles that belong to two boundaries so that
     // they belong to only one side
     size_t n_node(0U);
     for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
       {
          // no flagging at all
          if ( (*bit) == 0 || (*bit) == IRREGULAR )
            {
               csmp_error.Note( ERROR, "VSetConverter<dim>::FlagEdges",
                               "Skipping node, the boundary flag of which could not be identified" );
               cout <<"\nNode "<< n_node <<", flagged: "<< (*bit) << endl;
            }
          
          // Edges along the back side (Edges 1 to 4)
          if ( approximatelyEqual(vset.Pz( n_node ),zmin,tol) )
             {
                // EDGE1
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = BACK_BOTTOM;
                // EDGE2
                if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = BACK_RIGHT;
                // EDGE3
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = BACK_TOP;
                // EDGE4
                if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = BACK_LEFT;
             }
             
          // Edges along the front side (Edges 9 to 12)
          if ( approximatelyEqual(vset.Pz( n_node ),zmax) )
             {
                // EDGE9
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = FRONT_BOTTOM;
                // EDGE10
                if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = FRONT_RIGHT;
                // EDGE11
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = FRONT_TOP;
                // EDGE12
                if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = FRONT_LEFT;
             }
          
          // Edges along the right side (Edges 6 and 7)
          if ( approximatelyEqual(vset.Px( n_node ),xmax)  )
             {
                // EDGE6
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = BOTTOM_RIGHT;
                // EDGE7
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = TOP_RIGHT;
             }
          
          // Edges along the left side (Edges 5 and 8)
          if ( approximatelyEqual(vset.Px( n_node ),xmin)  )
             {
                // EDGE5
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = BOTTOM_LEFT;
                // EDGE8
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = TOP_LEFT;
             }
                                        
       }
       
 } // end FlagEdges



template<uint32_t dim>
void VSetConverter<dim>::OrderBarycentricQuadraticTriangleCoordinateOrigins( VSet<dim>& vset ) const
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::OrderBarycentricQuadraticTriangleCoordinateOrigins", 
                                     "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::OrderBarycentricQuadraticTriangleCoordinateOrigins", 
                                     "This method only works for barycentric quadratic triangle elements, vs.", 
                          parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }
       
    assert( dim == 2U );
    // 1. while looping over all triangles build ordered triangles using the first three
    //    nodes
    typename deque<vector<int64_t> >::iterator  pit(vset.PlistBegin());
    vector<size_t>                              pdata(7);
    typename deque<vector<int64_t> >::iterator  fit(vset.PfvertsBegin());
    vector<int64_t>                             pfvert(3);
    map<mjl::Point,size_t>                     ordered_nodes;
    typename map<mjl::Point,size_t>::const_iterator  it;

    while ( pit!=vset.PlistEnd() and fit!=vset.PfvertsEnd() )
      {
         // finding node of triangle origin
         ordered_nodes.erase( ordered_nodes.begin(), ordered_nodes.end() );
         for ( auto i{0U}; i<3U; i++ )
           ordered_nodes[ mjl::Point(vset.Px((*pit)[i]), vset.Py((*pit)[i])) ] = i;
         size_t offset = (*ordered_nodes.begin()).second;
                                                   
         // resetting the plist entry, rememebering that the nodes were ordered
         // counterclockwise initially (only needed if the first node in the 
         // triangle is not the first node in the plist entry)
         if ( offset != 0 )
           { 
              if ( offset == 1 )
                {
                   pdata[0] = (*pit)[1];
                   pdata[1] = (*pit)[2];
                   pdata[2] = (*pit)[0];
                   pdata[3] = (*pit)[4];
                   pdata[4] = (*pit)[5];
                   pdata[5] = (*pit)[3];
                   pdata[6] = (*pit)[6];
                   pfvert[0] = (*fit)[1];
                   pfvert[1] = (*fit)[2];
                   pfvert[2] = (*fit)[0];
                } 
              else if ( offset == 2 )
                {
                   pdata[0] = (*pit)[2];
                   pdata[1] = (*pit)[0];
                   pdata[2] = (*pit)[1];
                   pdata[3] = (*pit)[5];
                   pdata[4] = (*pit)[3];
                   pdata[5] = (*pit)[4];
                   pdata[6] = (*pit)[6];
                   pfvert[0] = (*fit)[2];
                   pfvert[1] = (*fit)[0];
                   pfvert[2] = (*fit)[1];
                } 
         
              // reassigning the new node list to the plist
              for ( size_t i{0U}; i<(*pit).size(); i++ ) (*pit)[i] = pdata[i];
         
              // reorganizing the neighbor element list 'pfverts' as well
              for ( size_t i{0U}; i<(*fit).size(); i++ ) (*fit)[i] = pfvert[i];
           }
         pit++; 
         fit++;
      }
      
 } // end 



// *******************************************************************************
//
//
// THREE DIMENSIONAL METHODS
//
//
// *******************************************************************************


/**
 
 void VSetConverter<dim>::ConvertLinearToBarycentricTetrahedra( VSet<dim>& vset ) 
 
 

\n \b Description 

ConvertLinearToBarycentricTetrahedra() transforms a 'plist' with linear
tetrahedral elements into a corresponding Vdata for 11-noded barycentric
tetrahedral elements numbered as follows. 

  
        
   ^ t 
   |             / s
   |    9      /
 3 o -- o --- o 2
   |  \     / |
   |    \ /   |
 7 o  6 o \   o 5
   |   /  8 o |
   | /       \|
 0 o -- o --- o 1  ----> r
        4       

 
The center node of the barycentric tetrahedron has the number 10 
(numbering 0...n). The missing node coordinates are calculated from the 
node coordinates of the linear tetrahedra and the lists px, py, pz'
are extended correspondingly.  

The property values are interpolated onto the newly generated points 
just as the node coordinates are.  

\n \b Input \b Arguments 

The supplied VSet must hold the "P-data" 'plist', 'pfverts' 'px', 'py', 'pz',
in order to be processable.  

\n \b Implementation 

The numbering of the new nodes starts from previous-nodes+1 and goes to
new-total nodes of the mesh. This approach
implies that, if node numbers are used to determine the placement in
global solution matrices, nodes from a single triangle will be separated
by many rows/cols in the global solution matrix. If a solver is used,
which applies no pre-conditioning of the solution matrix, before inverting
the solution matrix, this numbering leads to high storage requirements and 
long convergence times. In this case, an algorithm like the Cuthill-McKhee 
node-numbering graph-tree traversal should be used to improve the matrix 
occupancy and to reduce the number of nodes which are far off the 
diagonal.  

\n \b Application 


\n \b Messages 

 
tested: */
template<uint32_t dim>
void VSetConverter<dim>::ConvertLinearToBarycentricTetrahedra( VSet<dim>& vset ) 
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToBarycentricTetrahedra", 
                                     "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != LINEAR_TETRAHEDRON and vset.ElementType(0U) != ISOPARAMETRIC_LINEAR_TETRAHEDRON ) {
          csmp_error.Note( ERROR, "VSetConverter<dim>::ConvertLinearToBarycentricTetrahedra", 
                                     "This method only works for linear tetrahedral elements, vs.", 
                          parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }
     
     bool                      debug(false);
     map<mjl::Point3D,int64_t>  nodeIDs;
     double                  x, y, z;

     xmin = xmax = vset.Px(0);
     ymin = ymax = vset.Py(0);
     zmin = zmax = vset.Pz(0);

     // 1. mapping already existing node points O.K.
     // ---------------------------------------
     for ( size_t i{0U}; i<vset.Vertices(); i++ )
       {
          x = vset.Px(i);
          y = vset.Py(i);
          z = vset.Pz(i);
          if ( x < xmin ) xmin = x;
          if ( x > xmax ) xmax = x;
          if ( y < ymin ) ymin = y;
          if ( y > ymax ) ymax = y;
          if ( z < zmin ) zmin = z;
          if ( z > zmax ) zmax = z;
          // nodes are numbered from 1...n-nodes
          nodeIDs[ mjl::Point3D(x,y,z) ] = static_cast<int64_t>(i + 1);
       }
               
     // 2. copying already existing boundary flags and values
     //    after flagging the corner nodes and the edges
     // -----------------------------------------------------
     FlagCornerNodes( vset, true );
     FlagEdges( vset );
     
     vector<std::int8_t>  bflags( vset.BFlagsBegin(), vset.BFlagsEnd() );

     // 3. looping through plist:
     // -------------------------
     //  - expanding node-ID vectors for each element
     //  - adding new node points to px, py, pz vectors
     //  - assigning boundary flags to new node points
     //  - interpolating properties to new node points
     //    (faceverts stay exactly as they were before)
     
     // iterator to test whether new node point is already part of the mesh or whether
     // it must be created     
     int64_t    nID(nodeIDs.size()+1); // new node ID tracker

     for ( auto pit=vset.PlistBegin(); pit!=vset.PlistEnd(); pit++ )
      {
         // The Plist node ID vector is resized and the new node coordinates are entered
         (*pit).reserve(11);

         // node 4 is initialized as midside node of segment 1 (nodes 0 & 1)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[0] ) + vset.Px( (*pit)[1] )) / 2.0;
         y = (vset.Py( (*pit)[0] ) + vset.Py( (*pit)[1] )) / 2.0;
         z = (vset.Pz( (*pit)[0] ) + vset.Pz( (*pit)[1] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 1 (nodes:"<< (*pit)[0] <<","<< (*pit)[1] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         pair<typename map<mjl::Point3D,int64_t>::iterator,bool>
           test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         // if the middle-node point already exists in the list
         // the node ID which was found for it inside the map is used
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 4("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              // (should also work in 3D)
              const int8_t bflag{BoundaryFlags3D( bflags, (*pit)[0], (*pit)[1] )};
              if ( bflag != 0 )
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 5 is initialized as midside node of segment 2 (nodes 1 & 2)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[1] ) + vset.Px( (*pit)[2] )) / 2.0;
         y = (vset.Py( (*pit)[1] ) + vset.Py( (*pit)[2] )) / 2.0;
         z = (vset.Pz( (*pit)[1] ) + vset.Pz( (*pit)[2] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 2 (nodes:"<< (*pit)[1] <<","<< (*pit)[2] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 5("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              const int8_t bflag{BoundaryFlags3D( bflags, (*pit)[1], (*pit)[2] )};
              if ( bflag != 0 )
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 6 is initialized as midside node of segment 3 (nodes 2 & 0)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[2] ) + vset.Px( (*pit)[0] )) / 2.0;
         y = (vset.Py( (*pit)[2] ) + vset.Py( (*pit)[0] )) / 2.0;
         z = (vset.Pz( (*pit)[2] ) + vset.Pz( (*pit)[0] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 3 (nodes:"<< (*pit)[2] <<","<< (*pit)[0] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 6("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              const int8_t bflag{BoundaryFlags3D( bflags, (*pit)[2], (*pit)[0] )};
              if ( bflag != 0 )
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 7 is initialized as midside node of segment 4 (nodes 0 & 3)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[0] ) + vset.Px( (*pit)[3] )) / 2.0;
         y = (vset.Py( (*pit)[0] ) + vset.Py( (*pit)[3] )) / 2.0;
         z = (vset.Pz( (*pit)[0] ) + vset.Pz( (*pit)[3] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 4 (nodes:"<< (*pit)[0] <<","<< (*pit)[3] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 7("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              const int8_t bflag{BoundaryFlags3D( bflags, (*pit)[0], (*pit)[3] )};
              if ( bflag != 0 )
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 8 is initialized as midside node of segment 5 (nodes 1 & 3)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[1] ) + vset.Px( (*pit)[3] )) / 2.0;
         y = (vset.Py( (*pit)[1] ) + vset.Py( (*pit)[3] )) / 2.0;
         z = (vset.Pz( (*pit)[1] ) + vset.Pz( (*pit)[3] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 5 (nodes:"<< (*pit)[1] <<","<< (*pit)[3] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 8("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              const int8_t bflag{BoundaryFlags3D( bflags, (*pit)[1], (*pit)[3] )};
              if ( bflag != 0 )
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 9 is initialized as midside node of segment 6 (nodes 0 & 1)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[2] ) + vset.Px( (*pit)[3] )) / 2.0;
         y = (vset.Py( (*pit)[2] ) + vset.Py( (*pit)[3] )) / 2.0;
         z = (vset.Pz( (*pit)[2] ) + vset.Pz( (*pit)[3] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 6 (nodes:"<< (*pit)[2] <<","<< (*pit)[3] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 9("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              const int8_t bflag{BoundaryFlags3D( bflags, (*pit)[2], (*pit)[3] )};
              if ( bflag != 0 )
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 10 is positioned at the center of the tetrahedron
         // ------------------------------------------------------
         x = (vset.Px( (*pit)[0] )+vset.Px( (*pit)[1] )+vset.Px( (*pit)[2] )+vset.Px( (*pit)[3] )) / 4.0;
         y = (vset.Py( (*pit)[0] )+vset.Py( (*pit)[1] )+vset.Py( (*pit)[2] )+vset.Py( (*pit)[3] )) / 4.0;
         z = (vset.Pz( (*pit)[0] )+vset.Pz( (*pit)[1] )+vset.Pz( (*pit)[2] )+vset.Pz( (*pit)[3] )) / 4.0;
 
         if ( debug )
           { 
              cout <<"\nCoordinates of center node of the tetrahedron: ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         // if the middle-node point already exists in the list
         // the node ID which was found for it inside the map is used
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 10("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              (*pit).push_back( nID++ );
          }
          
      } // end plist loop (pit)
 
    // testing whether the new numbers are O.K.
    if ( debug ) 
      {
         cout <<"\nListing old and new nodes and their coordinates:";
         for ( auto ndit=nodeIDs.begin(); ndit!=nodeIDs.end(); ndit++ )
           {   
              cout <<"\nx,y,id: "<< endl;
              (*ndit).first.Out();
           }
         cout << endl << endl;
      }
 
       
    // 4. Creating new 'px' and 'py' arrays and assigning them to VSet<dim>
    // -------------------------------------------------------------------
    deque<double>  px( nodeIDs.size() ),
                   py( nodeIDs.size() ),  // new node-point coordinates
                   pz( nodeIDs.size() );

    for ( auto n : nodeIDs )
      {
         px[ n.second ] = n.first.x_;
         py[ n.second ] = n.first.y_;
         pz[ n.second ] = n.first.z_;
      }  
    vset.AddXYZ( px, py, pz );
    
    // 5. Converting CSMP finite element types 
    // ------------------------------------------------------------------------
    vset.ElementType( 0U, ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ); 
    
    // 6. Converting triangle node order such that first node is in lower-left corner of bounding box
    // ----------------------------------------------------------------------------------------------
  
 } // end ConvertLinearToBarycentricTetrahedra








/**
 
  void VSetConverter<dim>::ConvertLinearToQuadraticTetrahedra( VSet<dim>& vset ) 
 
 

\n \b Description 

Interpolates midside node locations onto the sides of the quadratic
tetrahedron and adds these points as nodes to the VSet<dim> which is 
returned.  

The node numbering of the tetrahedron is:
          
          ^ t
          |
          |
          
        4 o
          | \
          |   \
       10 o     o 9
          |      \
   8 o    |        \
        1 o -- o --- o 3  --> s
         /     7
       /         
   5 o     o 6
   /
  /
o 2   
r   

\n \b Input \b Arguments 

A reference to the VSet<dim> which shall be modified by the function.  

\n \b Output \b Arguments &amp; Return Value 

The modified VSet<dim> in returned into the first method argument.  

\n \b Implementation 


\n \b Application 

To facilitate computations with quadratic tetrahedral elements on the 
basis of simple tetrahedral element meshes.  

\n \b Messages 

 
// tested: O.K. SKM 4/3/02 */
template<uint32_t dim>
void VSetConverter<dim>::ConvertLinearToQuadraticTetrahedra( VSet<dim>& vset ) 
  {
     cout <<"\nVSetConverter<dim>::ConvertLinearToQuadraticTetrahedra: Warning: Output node numbering ";
     cout <<"is incompatible with IsoparametricQuadraticTetrahedron !"<< endl;
  
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
     if ( vset.HybridElementTypeMesh() ) {
          csmp_error.Note( FATAL_ERROR, "VSetConverter<dim>::ConvertLinearToQuadraticTetrahedra", 
                                     "This method does not work for mixed element meshes" );
          return;
       }
     if ( vset.ElementType(0U) != LINEAR_TETRAHEDRON and vset.ElementType(0U) != ISOPARAMETRIC_LINEAR_TETRAHEDRON ) {
          csmp_error.Note( FATAL_ERROR, "VSetConverter<dim>::ConvertLinearToQuadraticTetrahedra", 
                         "This method only works for linear tetrahedral elements, vs.", 
                          parseFiniteElementType(parseFiniteElementTypeEnum(vset.ElementType(0U))) );
          return;
       }
     
     map<mjl::Point3D,int64_t>  nodeIDs;
     double                  x, y, z;
     bool                      debug(false);

     xmin = xmax = vset.Px(0);
     ymin = ymax = vset.Py(0);
     zmin = zmax = vset.Pz(0);

     // 1. mapping already existing node points O.K.
     // ---------------------------------------
     for ( size_t i{0U}; i<vset.Vertices(); i++ )
       {
          x = vset.Px(i);
          y = vset.Py(i);
          z = vset.Pz(i);
          if ( x < xmin ) xmin = x;
          if ( x > xmax ) xmax = x;
          if ( y < ymin ) ymin = y;
          if ( y > ymax ) ymax = y;
          if ( z < zmin ) zmin = z;
          if ( z > zmax ) zmax = z;
          // nodes are numbered from 1...n-nodes
          nodeIDs[ mjl::Point3D(x,y,z) ] = i + 1;
       }
               
     // 2. copying already existing boundary flags and values
     //    after flagging the corner nodes and the edges
     // -----------------------------------------------------
     vector<std::int8_t>  bflags( vset.BFlagsBegin(), vset.BFlagsEnd() );

     // 3. looping through plist:
     // -------------------------
     //  - expanding node-ID vectors for each element
     //  - adding new node points to px, py, pz vectors
     //  - assigning boundary flags to new node points
     //  - interpolating properties to new node points
     //    (faceverts stay exactly as they were before)
     
     // iterator to test whether new node point is already part of the mesh or whether
     // it must be created
     int64_t    nID(nodeIDs.size()+1); // new node ID tracker
     int8_t   bflag;

     for ( auto pit=vset.PlistBegin(); pit!=vset.PlistEnd(); pit++ )
      {
         // The Plist node ID vector is resized and the new node coordinates are entered
         (*pit).reserve(10);

         // node 4 is initialized as midside node of segment 1 (nodes 0 & 1)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[0] ) + vset.Px( (*pit)[1] )) / 2.0;
         y = (vset.Py( (*pit)[0] ) + vset.Py( (*pit)[1] )) / 2.0;
         z = (vset.Pz( (*pit)[0] ) + vset.Pz( (*pit)[1] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 1 (nodes:"<< (*pit)[0] <<","<< (*pit)[1] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         pair<typename map<mjl::Point3D,int64_t>::iterator,bool>
           test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         // if the middle-node point already exists in the list
         // the node ID which was found for it inside the map is used
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 4("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              // If new point lies at the model boundary a boundary flag is assigned to the new point
              // (should also work in 3D)
              if ( (bflag=BoundaryFlags3D( bflags, (*pit)[0], (*pit)[1] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 5 is initialized as midside node of segment 2 (nodes 1 & 2)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[1] ) + vset.Px( (*pit)[2] )) / 2.0;
         y = (vset.Py( (*pit)[1] ) + vset.Py( (*pit)[2] )) / 2.0;
         z = (vset.Pz( (*pit)[1] ) + vset.Pz( (*pit)[2] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 2 (nodes:"<< (*pit)[1] <<","<< (*pit)[2] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 5("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              if ( (bflag=BoundaryFlags3D( bflags, (*pit)[1], (*pit)[2] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 6 is initialized as midside node of segment 3 (nodes 2 & 0)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[2] ) + vset.Px( (*pit)[0] )) / 2.0;
         y = (vset.Py( (*pit)[2] ) + vset.Py( (*pit)[0] )) / 2.0;
         z = (vset.Pz( (*pit)[2] ) + vset.Pz( (*pit)[0] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 3 (nodes:"<< (*pit)[2] <<","<< (*pit)[0] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 6("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              if ( (bflag=BoundaryFlags3D( bflags, (*pit)[2], (*pit)[0] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 7 is initialized as midside node of segment 4 (nodes 1 & 3)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[1] ) + vset.Px( (*pit)[3] )) / 2.0;
         y = (vset.Py( (*pit)[1] ) + vset.Py( (*pit)[3] )) / 2.0;
         z = (vset.Pz( (*pit)[1] ) + vset.Pz( (*pit)[3] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 4 (nodes:"<< (*pit)[1] <<","<< (*pit)[3] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 7("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              if ( (bflag=BoundaryFlags3D( bflags, (*pit)[1], (*pit)[3] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 8 is initialized as midside node of segment 5 (nodes 2 & 3)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[2] ) + vset.Px( (*pit)[3] )) / 2.0;
         y = (vset.Py( (*pit)[2] ) + vset.Py( (*pit)[3] )) / 2.0;
         z = (vset.Pz( (*pit)[2] ) + vset.Pz( (*pit)[3] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 5 (nodes:"<< (*pit)[2] <<","<< (*pit)[3] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 8("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              if ( (bflag=BoundaryFlags3D( bflags, (*pit)[2], (*pit)[3] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
         // node 9 is initialized as midside node of segment 6 (nodes 0 & 3)
         // ----------------------------------------------------------------
         x = (vset.Px( (*pit)[0] ) + vset.Px( (*pit)[3] )) / 2.0;
         y = (vset.Py( (*pit)[0] ) + vset.Py( (*pit)[3] )) / 2.0;
         z = (vset.Pz( (*pit)[0] ) + vset.Pz( (*pit)[3] )) / 2.0;
 
         if ( debug )
           { 
              cout <<"\nMidside node of segment 6 (nodes:"<< (*pit)[0] <<","<< (*pit)[3] <<"): ";
              cout << x <<", "<< y <<", "<< z << endl;
           }
 
         test_it = nodeIDs.insert( make_pair(mjl::Point3D(x,y,z),nID) );
         if ( !test_it.second ) (*pit).push_back( (*test_it.first).second );
         else                   
          {
              if ( debug ) 
                cout <<"\nNew coordinates, node 9("<< nID <<"): "<< x <<", "<< y <<", "<< z << endl;
              if ( (bflag=BoundaryFlags3D( bflags, (*pit)[0], (*pit)[3] )) != 0 ) 
                {
                   vset.BFlag( nID, bflag );
                }
             (*pit).push_back( nID++ );
          }
          
      } // end plist loop (pit)
 
    // testing whether the new numbers are O.K.
    if ( debug ) 
      {
         cout <<"\nListing old and new nodes and their coordinates:";
         for ( auto ndit=nodeIDs.begin(); ndit!=nodeIDs.end(); ndit++ )
           {   
              cout <<"\nx,y,id: "<< endl;
              (*ndit).first.Out();
           }
         cout << endl << endl;
      }
 
       
    // 4. Creating new 'px', 'py', and 'pz' arrays and assigning them to VSet<dim>
    // --------------------------------------------------------------------------
    deque<double>  px( nodeIDs.size() ),
                      py( nodeIDs.size() ),  // new node-point coordinates
                      pz( nodeIDs.size() );

    for ( auto nit=nodeIDs.begin(); nit!=nodeIDs.end(); nit++ )
      {
         px[ (*nit).second ] = (*nit).first.x_;
         py[ (*nit).second ] = (*nit).first.y_;
         pz[ (*nit).second ] = (*nit).first.z_;
      }  
    vset.AddXYZ( px, py, pz );
    

    // 5. Converting CSMP finite element types 
    // ------------------------------------------------------------------------
    vset.ElementType( 0U, QUADRATIC_TETRAHEDRON ); 

    // 5. Interpolating Node properties to the new node points if there are any
    // ------------------------------------------------------------------------
    
    
    // 6. Converting triangle node order such that first node is in lower-left corner of bounding box
    // ----------------------------------------------------------------------------------------------
  
 } // end ConvertLinearToQuadraticTetrahedra









// brute force approach: based on their position the nodes are flagged as boundary nodes
// tested: O.K.
template<uint32_t dim>
void VSetConverter<dim>::EstablishBoundaryFlagsForBoxModel( VSet<dim>& vset, double tol )
 {
     // 0. getting rid of the existing boundary conditions
     // --------------------------------------------------
     vset.RemoveBflags();

     // 1. finding the dimensions of the box-shaped model
     //    assuming that the boundary have been identified
     //    but the flags are incorrect thus far
     // -------------------------------------------------
     xmin = xmax = vset.Px(0);
     ymin = ymax = vset.Py(0);
     zmin = zmax = vset.Pz(0);

     double  x, y, z;

     for ( size_t i{0U}; i<vset.Vertices(); i++ )
       {
          x = vset.Px(i);
          y = vset.Py(i);
          z = vset.Pz(i);
          if ( x < xmin ) xmin = x;
          if ( x > xmax ) xmax = x;
          if ( y < ymin ) ymin = y;
          if ( y > ymax ) ymax = y;
          if ( z < zmin ) zmin = z;
          if ( z > zmax ) zmax = z;
       }

     // making new boundary conditions
     vector<std::int8_t>  new_bflags( vset.Vertices(), 0 );

     for ( size_t i{0U}; i<vset.Vertices(); i++ )  
       if ( approximatelyEqual(vset.Px(i),xmin,tol) ||
            approximatelyEqual(vset.Px(i),xmax,tol) ||
            approximatelyEqual(vset.Py(i),ymin,tol) ||
            approximatelyEqual(vset.Py(i),ymax,tol) ||
            approximatelyEqual(vset.Pz(i),zmin,tol) ||
            approximatelyEqual(vset.Pz(i),zmax,tol) )
         { 
            new_bflags[ i ] = UNSPECIFIED;
         }

     vset.AddBFlags( new_bflags.begin(), new_bflags.end() );

               
     // 2. flagging the sides of the model using the user-defined
     //    tolerances (the edges are dealt with later)
     // ---------------------------------------------------------
     size_t n_node(0U);
     for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
       {
          // bottom (y=ymin)
          if      ( approximatelyEqual(vset.Py( n_node ),ymin,tol) ) (*bit) = BOTTOM_OUTSIDE;
          // top    (y=ymax)
          else if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) ) (*bit) = TOP_OUTSIDE;
          // left   (x=xmin)
          else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) ) (*bit) = LEFT_OUTSIDE;
          // right  (x=xmax)
          else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) ) (*bit) = RIGHT_OUTSIDE;
          // back   (z=zmin)
          else if ( approximatelyEqual(vset.Pz( n_node ),zmin,tol) ) (*bit) = BACK_OUTSIDE;
          // front  (z=zmax)
          else if ( approximatelyEqual(vset.Pz( n_node ),zmax,tol) ) (*bit) = FRONT_OUTSIDE;
       }
       
     // 3. now the corners & edges are flagged
     // --------------------------------------
     FlagEdges( vset );  
     FlagCornerNodes( vset, true );
       
 } // end
    

template class VSetConverter<1U>;
template class VSetConverter<2U>;
template class VSetConverter<3U>;



} // end namespace csmp




