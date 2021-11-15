#include "TRIANGLE_Interface.h"
#include "Box.h"
#include "LinearTriangle.h"
#include "IsoparametricLinearTriangle.h"
#include "Exception.h"
#include "ErrorHandler.h"
//#include "FEM_Data.h"
#include "PropertyData.h"

using namespace std;

namespace csmp {

/**

Reads the series of mesh output files from Setchuk's triangle mesher. The
data from the text files are used to build a VSet which holds all the
necessary information to initialize a CSMP mesh, including property data
and boundary flags. 

However, Setchuk's mesher numbers nodes and elements from 1..n and 
CSMP wants them numbered 0..n-1. Therefore a conversion is performed
to adjust this numbering before the VSet is returned. 

The methods expects x files with the extensions '.edge', 'ele', '.neigh',
'.node' and '.poly'. The neighbor element file ('.neigh') is not generated
automatically by triangle, but must be prompted for by using the '-n'
command line option. It contains a list of neighboring elements for each
triangle. 

@section arguments Input Arguments 

The name of the file series without extension, and a VSet into which
the data shall be stored. 

@param vset The file data are returned into the second method argument. 

@section application Application

To read set of input files (.poly, .node, .ele, .neigh) from J. Setchuk's 2D
meshing programm 'Triangle' and puts the results into the supplied
VSet. 

@section messages Messages 

If one of the files cannot be found the method will report this as error
and terminate. */
template<size_t dim>                        
void TRIANGLE_Interface::ReadTriangle2DMesh( const char* fname, VSet<dim>& vset, 
                                             bool isoparametric, bool perform_extra_checks )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    char  file[31];
    strcpy( file, fname );
 	
	// mesh characteristics
    size_t  nodes_per_element;
    size_t  nodes;
    size_t  elements;
    int32_t   node_attributes;
    int32_t   elmt_attributes;
    bool    verbose(false);
    
    // temporary mesh storage
    deque<double>              x, y, z; 
    vector<double>             evalues;
    map<size_t,vector<int64_t> >  plist;
    map<size_t,vector<int64_t> >  pfverts;
    vector<std::int8_t>          bflags;
    map<size_t,double>         bvalues;

    // 1. Read input files 
    // ---------------------------------------------------------------
    CheckTriangleOutput( file, nodes_per_element, nodes, elements, 
                         node_attributes, elmt_attributes );
                              
    ReadNodeDataFile(     file, x, y, z, bflags, bvalues );
    ReadElementDataFile(  file, plist, evalues );
    ReadNeighborDataFile( file, pfverts ); // missing neighbors are flagged -1
    
    // 2. Checking for consecutive node numbering
    // ------------------------------------------
    if ( perform_extra_checks && !VerifyConsecutiveNodeNumbering( plist ) )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadTriangle2DMesh", 
            "Nodes in input mesh are not numbered consecutivly");
    
    // 3. If triangles span the element corners
    //    the two adjacent triangles are flipped such that
    //    the new segment splits the triangle into 2 corner elements
    // -------------------------------------------------------------
    SplitSingleCornerElements( plist, pfverts, evalues );

    // 4. BOX_BOUNDARY flags are identified 
    // -----------------------------------
    FlagBoundaryElements( pfverts, plist, x, y );
    FlagBoundaryNodes( pfverts, plist, bflags );
    FlagCornerNodes( bflags, x, y, z );

    // 5. Assigning and testing the permeability data
    // ----------------------------------------------
    if (  evalues.size() < plist.size() )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadTriangle2DMesh", 
                            "There were no element property values read");

    if ( evalues.size() != plist.size() )
      csmp_error.notice( ERROR, "TRIANGLE_Interface::ReadTriangle2DMesh", 
                            "This method assumes that there is only one property value per element." );
    
    ListZeroPropertyValueElements( evalues, plist, x, y );

    // creating the storage for the scalar permeability values
    PropertyData mesh_regions( ELEMENT, SCALAR, 2U );
    mesh_regions.Reserve( evalues.size(), evalues.size() );
    for ( size_t i=0U; i<evalues.size(); i++ )
      pushBack( mesh_regions, makeScalar(PLAIN,evalues[i]) );

    evalues.erase( evalues.begin(), evalues.end() );
  
    // 6. Putting mesh into the VSet
    // -----------------------------
    if ( !isoparametric ) vset.Resize( LinearTriangle().Nodes(),
                                       LinearTriangle().Neighbors(),
                                       LinearTriangle().ElementType(), x.size(), plist.size() );
    else                  vset.Resize( IsoparametricLinearTriangle().Nodes(),
                                       IsoparametricLinearTriangle().Neighbors(),
                                       IsoparametricLinearTriangle().ElementType(), x.size(), plist.size() );
    vset.AddXYZ( x, y, z );
    vset.AddPlist( plist.begin(), plist.end() );
    vset.AddPfverts( pfverts.begin(), pfverts.end() );
    vset.AddBFlags( bflags.begin(), bflags.end() );
    vset.AddData( "permeability", mesh_regions );
    
    // 7. Checking the VSet
    // --------------------
    vset.EstablishZeroBasedNumbering();
    if ( verbose ) vset.Out(); 
   
} // end ReadTriangle2DMesh











/// not ported to new numbering of elements and nodes yet
template<size_t dim>                        
void TRIANGLE_Interface::ReadTriangle2DMeshAndCreateDiscreteFractures( const char* fname, 
                                                                       VSet<dim>& vset, 
                                                                       double kfrac, 
                                                                       bool perform_extra_checks )
{
    char  file[31];
    strcpy( file, fname );
	
	// mesh characteristics
    size_t          nodes_per_element;
    size_t          nodes;
    size_t          elements;
    int32_t              node_attributes;
    int32_t              elmt_attributes;
    bool               verbose(false);
    
    // temporary mesh storage
    deque<double>              x, y, z; 
    vector<double>             evalues;
    vector<ScalarVariable >      sc_values;
    map<size_t,vector<int64_t> >  plist_tria;
    map<size_t,vector<int64_t> >  plist_bar;
    map<size_t,vector<int64_t> >  pfverts_tria;
    map<size_t,vector<int64_t> >  pfverts_bar;
    std::vector<std::int8_t>     bflags;
    map<size_t,double>         bvalues;

    // 1. Read input files 
    // ---------------------------------------------------------------
    CheckTriangleOutput( file, nodes_per_element, nodes, elements, 
                         node_attributes, elmt_attributes );
                              
    ReadNodeDataFile(     file, x, y, z, bflags, bvalues );
    ReadElementDataFile(  file, plist_tria, evalues );
    ReadPolyDataFile( file, plist_bar );
    FindNeighborsForFractureElements( plist_bar, pfverts_bar );
    ReadNeighborDataFile( file, pfverts_tria ); // missing neighbors are flagged -1
    
    // 2. Checking for consecutive node numbering
    // ------------------------------------------
    if ( perform_extra_checks && !VerifyConsecutiveNodeNumbering( plist_tria ) )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadTriangle2DMeshAndCreateDiscreteFractures", 
            "Nodes in input mesh are not numbered consecutivly");
    
    // 3. If triangles span the element corners
    //    the two adjacent triangles are flipped such that
    //    the new segment splits the triangle into 2 corner elements
    // -------------------------------------------------------------
    SplitSingleCornerElements( plist_tria, pfverts_tria, evalues );

    // 4. BOX_BOUNDARY flags are identified 
    // -----------------------------------
    FlagBoundaryElements( pfverts_tria, plist_tria, x, y );
    FlagBoundaryNodes( pfverts_tria, plist_tria, bflags );
    FlagCornerNodes( bflags, x, y, z );

    // 5. Assigning and testing the permeability data
    // ----------------------------------------------
    if ( evalues.size() != plist_tria.size() )
      throw csmp::Exception( ERROR, "TRIANGLE_Interface::ReadTriangle2DMeshAndCreateDiscreteFractures", "Element property array has wrong size." );
    
    if ( kfrac == 0. )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadTriangle2DMeshAndCreateDiscreteFractures", "Fracture permeability is set to zero." );

    // additional bar elements
    size_t total_elements(plist_tria.size()+plist_bar.size()), tria_elements(plist_tria.size());
    for ( size_t i=plist_tria.size(); i<total_elements; i++ ) evalues.push_back( kfrac );

    PropertyData mesh_regions( ELEMENT, SCALAR, 2U );
    mesh_regions.Reserve( evalues.size(), evalues.size() );
    for ( size_t i=0U; i<evalues.size(); i++ )
      pushBack( mesh_regions, makeScalar(PLAIN,evalues[i]) );

    evalues.erase( evalues.begin(), evalues.end() );
  
    // 6. Merge containers for TRIA and BAR elements
    // ---------------------------------------------
    map<size_t,vector<int64_t> > plist;
    map<size_t,vector<int64_t> > pfverts;
    plist = plist_tria;
    for ( size_t i=1; i<=plist_bar.size(); i++ )   plist[tria_elements+i]   = plist_bar[i];
    pfverts = pfverts_tria;
    for ( size_t i=1; i<=pfverts_bar.size(); i++ ) pfverts[tria_elements+i] = pfverts_bar[i];
    
    deque<size_t>  ndele(total_elements); // nodes per element and neighbours per element are identical for BAR and TRIA
                       
    for ( size_t i=0; i<plist_tria.size(); i++ )                  
      ndele[i] = 3U; 
    for ( size_t i=plist_tria.size(); i<total_elements; i++ )                  
      ndele[i] = 2U; 
    vset.ResizePlist( ndele );
    vset.ResizePfverts( ndele );
     
    // 7. Putting results into the VSet
    // ------------------------------------
    vset.AddXYZ( x, y, z );
    vset.AddPlist( plist.begin(), plist.end() );
    vset.AddPfverts( pfverts.begin(), pfverts.end() );
    vset.AddBFlags( bflags.begin(), bflags.end() );
    vset.AddData( "permeability", mesh_regions );
    
    vector<int8_t> elmt_types(total_elements);
    for ( size_t i=0U; i<tria_elements; i++ )
      //vset.ElementType( i, ISOPARAMETRIC_LINEAR_TRIANGLE );
      elmt_types[i] = ISOPARAMETRIC_LINEAR_TRIANGLE;
    for ( size_t i=tria_elements; i<total_elements; i++ )
      //vset.ElementType( i, ISOPARAMETRIC_LINEAR_BAR );
      elmt_types[i] = ISOPARAMETRIC_LINEAR_BAR;      
    vset.AddElementTypes( elmt_types.begin(), elmt_types.end() );

    
    // 7. Checking the VSet
    // --------------------
    vset.EstablishZeroBasedNumbering();
    if ( verbose ) vset.Out(); 
   
} // end ReadTriangle2DMesh




void TRIANGLE_Interface::CheckTriangleOutput( const char*  file, 
                                         size_t&   nodes_per_element,
                                         size_t&  nodes, 
                                         size_t&  elements,
                                         int&   node_attributes,
                                         int&   elmt_attributes )
 {
    char fnode[200], fele[200], fedge[200];
    strcpy( fnode, file );
    strcat( fnode, ".node");
    strcpy( fele, file );
    strcat( fele, ".ele");
    strcpy( fedge, file );
    strcat( fedge, ".edge"); // invoke with -e flag of Triangle
    
    // getting the size information out of the input files
    // ---------------------------------------------------
    size_t segms, bms, dimension;
    
    // reading the first line of '.node' file
    ifstream ifn( fnode );
    
    if ( !ifn.is_open() )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::CheckTriangleOutput", 
                     "'*.node' input file could not be opened. Does it exist (in the same directory as program) ?" );
    
    ifn >> nodes >> dimension >> node_attributes >> bms;
    ifn.close();

    // reading the first line of '.ele' file
    ifstream ife( fele );
    
    if ( !ife.is_open() )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::CheckTriangleOutput", 
                     "'*.ele' input file could not be opened. Does it exist (in the same directory as program) ?" );
   
    ife >> elements >> nodes_per_element >> elmt_attributes >> bms;
    ife.close();
    
    if ( nodes_per_element != 3 )
      cout <<"\nMeshInterface::CheckTriangleOutput: '.ele' file contains "<< nodes_per_element <<"-noded triangles..."<< endl;

    // reading the first line of '.edge' file
    ifstream ifd( fedge );
    
    if ( !ifd.is_open() )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::CheckTriangleOutput", 
                     "'*.edge' input file could not be opened. Does it exist (in the same directory as program) ?" );
    
    ifd >> segms >> bms;
    ifd.close();
    
    if ( elmt_attributes == 0 ) 
      cout <<"\nMeshInterface::CheckTriangleOutput: '.ele' file contains no element attributes..."<< endl;

    cout <<"\nMeshInterface::CheckTriangleOutput '"<< file;
    cout <<"' Input fileset apparently contains: " << endl;
    cout << nodes <<" nodes, "<< segms <<" segments=faces, "<< elements <<" triangles."<< endl; 
    
 } // end CheckTriangleOutput




void TRIANGLE_Interface::ReadNodeDataFile( const char* file, 
                                           deque<double>& x, deque<double>& y, 
                                           deque<double>& z,
                                           vector<std::int8_t>&  bflags,
                                           map<size_t,double>& bvalues )
 {
    char    fname[200];
    strcpy( fname, file );
    strcat( fname, ".node");
    ifstream ifs( fname );
    
    assert( ifs.is_open() );

    size_t    id, nodes;
    int32_t     pbflag(0), dim, node_attributes, boundary_markers;
    double  pbval, xval, yval;

    // erasing vectors and maps
    if ( !x.empty() ) x.erase( x.begin(), x.end() );
    if ( !y.empty() ) y.erase( y.begin(), y.end() );
    if ( !z.empty() ) z.erase( z.begin(), z.end() );
    if ( !bflags.empty() ) bflags.erase( bflags.begin(), bflags.end() );
    if ( !bvalues.empty() ) bvalues.erase( bvalues.begin(), bvalues.end() );

    // reading the header line and discarding it because the info is already there
    ifs >> nodes >> dim >> node_attributes >> boundary_markers;
    
    // reading the body of data
    for( size_t i=0U; i<nodes; i++ )
      {
         //     id    x-coordinate   y-coordinate   
         ifs >> id >> xval >> yval; // node attributes
         for ( int32_t j=0; j<node_attributes; j++ ) ifs >> pbval;
         if ( boundary_markers > 0 ) ifs >> pbflag; // boundary flags
         x.push_back( xval );
         y.push_back( yval );
         z.push_back( 0. );
         if ( id == 0 || id > nodes )
           throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadNodeDataFile", 
                     "encountered erratic node ID's (nodes should be numbered 1...n)" );

         // boundary values and flags are stored only if a node is at the model
         // boundary
         if ( pbflag == 1 )
           {
              bvalues[ id ] = pbval;
              bflags[ id ]  = pbflag;
           }
      }
    ifs.close();

    cout <<"\nMeshInterface::ReadNodeDataFile: file '"<< fname;
    cout <<"' read successfully..." << endl;
    
 } // end ReadNodeDataFile





void TRIANGLE_Interface::ReadElementDataFile( const char* file, 
                                              map<size_t,vector<int64_t> >& plist,
                                              vector<double>& evalues )
 {
    char fname[200], text_line[256], *token;
    const char* const  delims = " ,=#%<>,\n,\r,\t";
    strcpy( fname, file );
    strcat( fname, ".ele");
    ifstream ifs( fname );
    
    if ( !ifs.is_open() )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadElementDataFile", 
                     "'*.ele' input file could not be opened. Does it exist (in the same directory as program) ?" );
    
    // erasing vectors and maps
    if ( !plist.empty() ) plist.erase( plist.begin(), plist.end() );
    if ( !evalues.empty() ) evalues.erase( evalues.begin(), evalues.end() );
    
    // file information
    long      elements;
    size_t    id, pval, points_per_element, n_attributes;
    int32_t     n_boundary_markers, bmark;
    double  eval;

    // reading the header line and discarding it because the info is already there
    ifs.getline( text_line, 256 );
    token              = strtok( text_line, delims );
    elements           = atol( token );   
    token              = strtok( NULL, delims );
    points_per_element = static_cast<size_t>(atol( token ));  
    token              = strtok( NULL, delims );
    n_attributes       = static_cast<size_t>(atol( token ));
    token              = strtok( NULL, delims );
    if ( token != NULL ) n_boundary_markers = atoi( token );
    else n_boundary_markers = 0;

    // evalues will have a size that corresponds to n-elements * n-attributes  
    evalues.reserve( static_cast<size_t>(elements) * n_attributes );
    vector<int64_t>  pdata(points_per_element);
    
    // reading the pdata
    for ( size_t i=0U; i<static_cast<size_t>(elements); i++ )
       {
          ifs >> id;
          for ( size_t j=0U; j<points_per_element; j++ )
            {
               ifs     >> pval;
               pdata[j] = pval;
            }
          for ( size_t l=0U; l<n_attributes; l++ ) ifs >> eval;
          evalues.push_back( eval );
          if  ( n_boundary_markers > 0 ) ifs >> bmark;
          plist[ id ] = pdata;
       }
    ifs.close();

    cout <<"\nMeshInterface::ReadElementDataFile: file '"<< fname;
    cout <<"' read successfully..." << endl;
    
 } // end ReadElementDataFile








void TRIANGLE_Interface::ReadPolyDataFile( const char* file, 
                                           map<size_t,vector<int64_t> >& plist )
 {
    char fname[200];
    strcpy( fname, file );
    strcat( fname, ".poly");
    ifstream ifs( fname );
    
    if ( !ifs.is_open() )
      throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadPolyDataFile", 
                     "'*.poly' input file could not be opened. Does it exist (in the same directory as program) ?" );
    
    // erasing vectors and maps
    if ( !plist.empty() ) plist.erase( plist.begin(), plist.end() );
    
    // file information
    int32_t     nodes, dimensions, attributes, markers, segments, bmark;
    size_t id, pval, fid(1);

    // reading the header line and discarding it because the info is already there
    ifs >> nodes >> dimensions >> attributes >> markers;
    ifs >> segments >> markers;
    
    if ( nodes != 0 ) 
    throw csmp::Exception( FATAL_ERROR, "TRIANGLE_Interface::ReadPolyDataFile", 
                     "Aparently there are nodes listed in '*.poly' file ?" );
                     
    vector<int64_t>  pdata(2); // bar elements have 2 nodes per element.
    
    // reading the pdata
    for ( size_t i=0U; i<static_cast<size_t>(segments); i++ )
       {
          ifs >> id;
          for ( size_t j=0; j<2; j++ )
            {
               ifs     >> pval;
               pdata[j] = pval;
            }
          ifs >> bmark;
          if ( bmark != 1 ) {
              plist[fid] = pdata; // only segments inside the model are considered
              fid++;
            }
       }
    ifs.close();
    
    cout <<"\nMeshInterface::ReadElementDataFile: file '"<< fname;
    cout <<"' read successfully..." << endl;
    
 } // end ReadElementDataFile





void TRIANGLE_Interface::FindNeighborsForFractureElements( map<size_t,vector<int64_t> >& plist,
                                                           map<size_t,vector<int64_t> >& pfverts )
{
  vector<int64_t> nbors(2);
  
  size_t i(1);
  
  for ( map<size_t,vector<int64_t> >::iterator plit1 = plist.begin(); plit1 != plist.end(); plit1++ ) {
      nbors[0] = nbors[1] = -1; // flag as -1 in case BAR element has no neighbor
      for ( map<size_t,vector<int64_t> >::iterator plit2 = plist.begin(); plit2 != plist.end(); plit2++) {
          if ( ( plit1->second[0] == plit2->second[0] || plit1->second[0] == plit2->second[1] ) && plit1 != plit2 )  
              nbors[0] = (plit2->first);
          if ( ( plit1->second[1] == plit2->second[0] || plit1->second[1] == plit2->second[1] ) && plit1 != plit2 )  
              nbors[1] = (plit2->first);
        }
      pfverts[i] = nbors;
      i++;
    }
}                                                                   





void TRIANGLE_Interface::ReadNeighborDataFile( const char* file, 
                                               map<size_t,vector<int64_t> >& pfverts )
{
    char  fname[200];
    strcpy( fname, file );
    strcat( fname, ".neigh");
    ifstream ifs( fname );
    
    assert( ifs.is_open() );
    
    // erasing map
    if ( !pfverts.empty() ) pfverts.erase( pfverts.begin(), pfverts.end() );
    
    // file information
    size_t  elements, neighbors_per_element, id;
    int64_t   nid;

    // reading the header line and discarding it because the info is already there
    ifs >> elements >> neighbors_per_element;
    vector<int64_t> nbors( neighbors_per_element );
    
    // reading the faceverts as being opposite to the nodes in the triangles
    for ( size_t i=0; i<elements; i++ )
      {
         ifs >> id;
         for ( size_t j=0; j<neighbors_per_element; j++ )
           {
              ifs     >> nid;
              nbors[j] = nid; 
           }
         pfverts[ id ] = nbors;   
      }     
    ifs.close();

    cout <<"\nMeshInterface::ReadNeighborDataFile: file '"<< fname;
    cout <<"' read successfully..." << endl;
    
 } // end ReadNeighborDataFile
  




/**

Private method used by ReadTriangle2DMesh(). 
This function flags the boundaries of the Mesh using the enum BOX_BOUNDARY.
For this purpose, missing neighbor elements which are flagged -1 in 
the '.neigh' file are used as indicator for elements which lie at the model
boundary.  

Is used when a 2D mesh created with the Triangle mesher is read from file. 
*/
void TRIANGLE_Interface::FlagBoundaryElements( map<size_t,vector<int64_t> >& pfverts,
                                               map<size_t,vector<int64_t> >& plist,
                                               deque<double>& bx, deque<double>& by )
 {
    const int32_t MINUS1 = -1;  // flag of program 'triangle' for no neighbor element 
    
    // 0. Testing whether the coordination can be made
    // -----------------------------------------------
    if ( pfverts.size() != plist.size() )
      {
         cout <<"\nMeshInterface::FlagBoundaryElements: size of 'pfverts' not equal to size of ";
         cout <<"'plist'. No flagging can be done."<< endl;
         return;
      }
    
    // 4. BOX_BOUNDARY flags are identified 
    // -----------------------------------
    vector<double>  n1(3), n2(3), n3(3);
    
    map<size_t,vector<int64_t> >::iterator  pit(plist.begin());
    map<size_t,vector<int64_t> >::iterator  it(pfverts.begin());
    
    for ( ; it!=pfverts.end(); it++, pit++ )
      {
         for ( size_t n=0; n<3; n++ )
            {
               // if a boundary face has been found
               if ( (*it).second[n] == MINUS1 )
                 { 
                    // in 'Triangle' 
                    // elements are numbered counter-clockwise 
                    // the first 3 nodes give the corner points irrespective of
                    // whether the element is linear or quadratic
                    n1[0] = bx[ (*pit).second[0]-1 ];
                    n1[1] = by[ (*pit).second[0]-1 ];
                    n1[2] = 0.0;
                    n2[0] = bx[ (*pit).second[1]-1 ];
                    n2[1] = by[ (*pit).second[1]-1 ];
                    n2[2] = 0.0;
                    n3[0] = bx[ (*pit).second[2]-1 ];
                    n3[1] = by[ (*pit).second[2]-1 ];
                    n3[2] = 0.0;

                    // finding the boundary face with 2 nodes flagged as ONE
                    // (face1 contains node 1&2, face2 2&3, face3 3&1 as in TRIANGLE_Interface::ReadNastran() )
                    // facevert1 is opposite node1, facevert2 is opposite node2...
                    // -----------------------------------------------------------
                    // facevert 1 (node2 and node3)
                    if      ( n == 0 )
                      // finding on which side of the model the face lies
                       (*it).second[0] = FindFaceBoundary( n2, n3, n1 );

                    // facevert 2 (node3 and node1)
                    else if ( n == 1 )
                       (*it).second[1] = FindFaceBoundary( n3, n1, n2 );
              
                    // facevert 3 (node1 and node2)
                    else if ( n == 2 )
                       (*it).second[2] = FindFaceBoundary( n1, n2, n3 );
                 } 
            }   
      }    
 } // end   






/**

Private method used by ReadTriangle2DMesh(). 

Nodes are flagged according to their parent faces across which the face-
vert vectors point. Thus, nodes 2,3 correspond to facevert 1, nodes 3,1
to facevert 2, and nodes 1,2 to facevert 3. The method does not identify
corner points of the model. This task is left to CSMP.  
*/
void TRIANGLE_Interface::FlagBoundaryNodes( map<size_t,vector<int64_t> >& pfverts,
                                            map<size_t,vector<int64_t> >& plist,
                                            std::vector<std::int8_t>& bflags )
 {
    // flag of program 'triangle' for boundary node
    map<size_t,vector<int64_t> >::iterator  pit;
    map<size_t,vector<int64_t> >::iterator  it;
    map<size_t,int64_t>::iterator           bit;
    
    assert( !bflags.empty() );

    for ( it=pfverts.begin(), pit=plist.begin(); it!=pfverts.end(); it++, pit++ ) 
      for ( size_t n=0; n<(*it).second.size(); n++ ) 
        {
           // if an element face is at the model boundary
           // -------------------------------------------
           if ( (*it).second[n] < 0 )
             { 
                // check that elements have been flagged properly before
                // -----------------------------------------------------
                assert( (*it).second[n] != -1 ); 

                // find the boundary flags of the elements nodes;
                // if these are ONE, i.e. if the node is at the model
                // boundary, then the node is flagged equivalent to the element
                // ------------------------------------------------------------
                //                          fvert  bflag          node-ID vector<double> bflag map
                FlagBoundaryNodesAccordingTo( n, (*it).second[n], (*pit).second, bflags );
             }
         } 

 } // end FlagBoundaryNodes






/**

Given the flagging and number (0...2) of boundary face, flag boundary
nodes accordingly using the flagging enumeration BOX_BOUNDARY.
The faceverts are required to lie opposite to the nodes with the same 
numbers. 
*/
void TRIANGLE_Interface::FlagBoundaryNodesAccordingTo( size_t fvert, int64_t  bflag_int,
                                                       const vector<int64_t>& nds,
                                                       std::vector<std::int8_t>&  bflags )
 {
    assert( fvert < 3 );
    // 1. parsing the boundary identifying integer to BOX_BOUNDARY enum
    BOX_BOUNDARY  bflag(NOT);
    if      ( bflag_int == TOP_OUTSIDE )    bflag = TOP;
    else if ( bflag_int == BOTTOM_OUTSIDE ) bflag = BOTTOM;
    else if ( bflag_int == LEFT_OUTSIDE )   bflag = LEFT;
    else if ( bflag_int == RIGHT_OUTSIDE )  bflag = RIGHT;
    else if ( bflag_int == FRONT_OUTSIDE )  bflag = FRONT;
    else if ( bflag_int == BACK_OUTSIDE )   bflag = BACK;
 
    // 2. Flagging the nodes
    if ( fvert == 0U ) {
         bflags[nds[1]] = bflag;
         bflags[nds[2]] = bflag;
      }
    else if ( fvert == 1U ) {
         bflags[nds[2]] = bflag;
         bflags[nds[0]] = bflag;
      }
    else if ( fvert == 2U ) {
         bflags[nds[0]] = bflag;
         bflags[nds[1]] = bflag;
      }
    
 } // end FlagBoundaryNodesAccordingTo                                       







/**

Loops through plist making a map of node numbers. If this map has jumps 
in the numbering, these are detected when looping through it again. 
*/
bool TRIANGLE_Interface::VerifyConsecutiveNodeNumbering( map<size_t,vector<int64_t> >& plist )
 const 
 {
    set<size_t>  node_numbers;
    size_t       counter(1);
    
    for ( map<size_t,vector<int64_t> >::const_iterator
          it=plist.begin(); it!=plist.end(); it++ )
      for ( vector<int64_t>::const_iterator
            nit=(*it).second.begin(); nit!=(*it).second.end(); nit++ ) 
        node_numbers.insert( *nit );
          
    for ( set<size_t>::const_iterator 
          sit=node_numbers.begin(); sit!=node_numbers.end(); sit++ )
      // if a gap in the numbering leads to a difference in incremented counter
      // and the ordered node numbers in the list, false is returned.
      if ( *sit != counter++ ) 
        {
           cout <<"\nMeshInterface::VerifyConsecutiveNodeNumbering:"<< endl;
           cout <<"\n     Node numbering gap in 'plist': Current node: " << endl << endl;
           cout <<"     "<< *sit     <<" versus "<< counter-1 <<", last correct node: ";
           cout << *(--sit) <<" versus "<< counter-2 << endl << endl;
           // to continue tracking
           counter = *(++sit);
           return false;
        }      
    
    cout <<"\nMeshInterface::VerifyConsecutiveNodeNumbering: Yes, nodes are numbered consecutively."<< endl;
    return true;
  
 } // end VerifyConsecutiveNodeNumbering                                         






/**

If a 2D triangular finite element has two of its sides on the boundary
of the model, this poses problems when boundary conditions are 
assigned because one always first has to identify the correct model side
and secondly because elements rarely are compliant enough to satisfactorily
handle the different conditions.  

The remedy used here, is to find the elements element neighbor and to 
create a segment from the corner node and the node diametrically opposite.
This segment then separates two new elements which are created by 
flipping the two adjacent elements, example lower left corner:  

 

 o-----o        o-----o   element 2 also becomes a corner element 
 | \ 2 |        | 1 / |
 |  \  |   -->  |  /  |
 | 1 \ |        | / 2 |
 o-----o        o-----o

 

If element properties were assigned before, these are averaged between the
two flipped elements and then re-assigned.  

@section arguments Input Arguments 

The method uses the 'plist', 'pfverts', and potential element property
values in the process of flipping the nodes.  

@section implementation Implementation

1. In the pfverts array those elements are singled out which have only one 
   neighbor. These elements must lie at the mesh corners.
   
2. If such an element is found, its neighbor is found as well in the plist.
   The element ID's are remembered.

3. A segment (MJL_Edge) is constructed from the corner node to the node
   diametrically opposite in the triangle.

4. Depending which mesh corner we are in, the two remaining nodes left and 
   right of the new segment form corner nodes of the new triangles
   (MJL_Triangle).
   
5. The new triangles are reassigned to the plist and the properties, if
   present are averaged and mapped back onto the elements. If both property
   values were the same, nothing is done.
   
6. The elements pfverts are re-assigned according to the new connectivity.

7. The elements adjacent to the flipped elements must be updated with regard
   to their neighbor lists.


@section application Application

The method is applied before the VSet is output from the 2D mesh reading
method above.  

@section messages Messages 

The method reports how many edge flips were performed and whether and
which property values were averaged at which mesh corner.  
 
        @todo (3) Messes up neighbor connectivity SKM12/5/02
        O.K. for trianglePatch8 retested by SKM 27/11/02 (previously reported error must be in Triangle mesher)
*/
void  TRIANGLE_Interface::SplitSingleCornerElements( map<size_t,vector<int64_t> >& plist,
                                                     map<size_t,vector<int64_t> >& pfverts,
                                                     vector<double>& evalues )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    map<size_t,vector<int64_t> >::iterator  pfit, neigh_it, findit;
    map<size_t,vector<int64_t> >::iterator  eit1, eit2;
    size_t                                 i, neighbors, cnr_node, opp_node, corner_elements(0);
    int64_t                                  cnr_id, neighbor_id;
    vector<size_t>                         new_pl1(3), new_pl2(3);
    vector<int64_t>                         new_pf1(3), new_pf2(3);
    double                               new_prop, eprop1, eprop2;
    
    
    for ( pfit=pfverts.begin(); pfit!=pfverts.end(); pfit++ )
      {
         // 1. In the pfverts array those elements are found which only have one 
         //    neighbor. These elements must lie at the mesh corners.
         for ( neighbors=i=0; i<(*pfit).second.size(); i++ ) if ( (*pfit).second[i] > 0 ) neighbors++;

         // 2. If such an element is found, its neighbor is found as well in the plist.
         //    The element ID's are remembered.
         if ( neighbors == 1 )
           {  
              for ( cnr_node=0; cnr_node<(*pfit).second.size(); cnr_node++ ) 
                if ( (*pfit).second[cnr_node] > 0 ) break;
              cnr_id      = ((*pfit).first);
              neighbor_id = (*pfit).second[cnr_node];

              cout <<"\nTRIANGLE_Interface::SplitSingleCornerElements: "<< endl;
              cout <<"\n   Corner element found:    "<< cnr_id;
              cout <<"\n   Its neighbor element is: "<< neighbor_id << endl;
              // 3. A segment (MJL_Edge) is constructed from the corner node to the node
              //    diametrically opposite in the triangle.
              //    (the nodes opposite the corner nodes are found because they correspond
              //    to the positions where the neighbor elements have the identified ID(s).
              assert( neighbor_id > 0 );
              neigh_it=pfverts.find( static_cast<size_t>(neighbor_id) );
              assert( neigh_it != pfverts.end() );
              for ( opp_node=0; opp_node<(*neigh_it).second.size(); opp_node++ ) 
                if ( (*neigh_it).second[opp_node] == cnr_id ) break;
                
              // the nodes of the new segment have now been identified as 'cnr_id' and 'neighbor_id'

              // 4. Depending which mesh corner we are in, the two remaining nodes left and 
              //    right of the new segment form corner nodes of the new triangles
              assert( cnr_id > 0 );
              assert( neighbor_id > 0 );
              eit1=plist.find(static_cast<int64_t>(cnr_id));
              assert( eit1 != plist.end() );
              eit2=plist.find(static_cast<int64_t>(neighbor_id));
              assert( eit2 != plist.end() );

              // ----------------------------------
              // 4.1 cases for first new triangle 1
              // ----------------------------------
              if ( cnr_node == 0 ) // case A1 doc
                {
                   // since nodes sit opposite neighbors
                   new_pl1[0] = (*eit1).second[0]; 
                   new_pl1[1] = (*eit1).second[1]; 
                   new_pl1[2] = (*eit2).second[opp_node]; 
                   // new pfverts dependent on numbering of neighbor triangle
                   if      ( opp_node == 0 ) new_pf1[0] = (*neigh_it).second[1];
                   else if ( opp_node == 1 ) new_pf1[0] = (*neigh_it).second[2];
                   else if ( opp_node == 2 ) new_pf1[0] = (*neigh_it).second[0];
                   new_pf1[1] = neighbor_id; 
                   new_pf1[2] = (*pfit).second[2];
                }
              else if ( cnr_node == 2 ) // case B1 doc
                {
                   // since nodes sit opposite neighbors
                   new_pl1[0] = (*eit1).second[0]; 
                   new_pl1[1] = (*eit2).second[opp_node]; 
                   new_pl1[2] = (*eit1).second[2];
                   // new pfverts dependent on numbering of neighbor triangle
                   new_pf1[0] = neighbor_id; 
                   new_pf1[1] = (*pfit).second[1];
                   if      ( opp_node == 0 ) new_pf1[2] = (*neigh_it).second[1];
                   else if ( opp_node == 1 ) new_pf1[2] = (*neigh_it).second[2];
                   else if ( opp_node == 2 ) new_pf1[2] = (*neigh_it).second[0];
                }
              else if ( cnr_node == 1 ) // case C1 doc
                {
                   // since nodes sit opposite neighbors
                   new_pl1[0] = (*eit2).second[opp_node]; 
                   new_pl1[1] = (*eit1).second[1]; 
                   new_pl1[2] = (*eit1).second[2];
                    // new pfverts dependent on numbering of neighbor triangle
                   new_pf1[0] = (*pfit).second[0]; 
                   if      ( opp_node == 0 ) new_pf1[1] = (*neigh_it).second[1];
                   else if ( opp_node == 1 ) new_pf1[1] = (*neigh_it).second[2];
                   else if ( opp_node == 2 ) new_pf1[1] = (*neigh_it).second[0];
                   new_pf1[2] = neighbor_id;
               }
              // -----------------------------------
              // 4.2 cases for second new triangle 2
              // -----------------------------------
              if ( opp_node == 0 ) // case A1 doc
                {
                   // since nodes sit opposite neighbors
                   new_pl2[0] = (*eit2).second[0]; 
                   new_pl2[1] = (*eit2).second[1]; 
                   new_pl2[2] = (*eit1).second[cnr_node]; 
                   // new pfverts dependent on numbering of neighbor triangle
                   if      ( cnr_node == 0 ) new_pf2[0] = (*pfit).second[1];
                   else if ( cnr_node == 1 ) new_pf2[0] = (*pfit).second[2];
                   else if ( cnr_node == 2 ) new_pf2[0] = (*pfit).second[0];
                   new_pf1[1] = cnr_id; 
                   new_pf1[2] = (*neigh_it).second[2];
                }
              else if ( opp_node == 2 ) // case B1 doc
                {
                   // since nodes sit opposite neighbors
                   new_pl2[0] = (*eit2).second[0]; 
                   new_pl2[1] = (*eit1).second[cnr_node]; 
                   new_pl2[2] = (*eit2).second[2];
                   // new pfverts dependent on numbering of neighbor triangle
                   new_pf2[0] = cnr_id; 
                   new_pf2[1] = (*neigh_it).second[1];
                   if      ( cnr_node == 0 ) new_pf2[2] = (*pfit).second[1];
                   else if ( cnr_node == 1 ) new_pf2[2] = (*pfit).second[2];
                   else if ( cnr_node == 2 ) new_pf2[2] = (*pfit).second[0];
                }
              else if ( opp_node == 1 ) // case C1 doc
                {
                   // since nodes sit opposite neighbors
                   new_pl2[0] = (*eit1).second[cnr_node]; 
                   new_pl2[1] = (*eit2).second[1]; 
                   new_pl2[2] = (*eit2).second[2];
                    // new pfverts dependent on numbering of neighbor triangle
                   new_pf2[0] = (*neigh_it).second[0]; 
                   if      ( cnr_node == 0 ) new_pf2[1] = (*pfit).second[1];
                   else if ( cnr_node == 1 ) new_pf2[1] = (*pfit).second[2];
                   else if ( cnr_node == 2 ) new_pf2[1] = (*pfit).second[0];
                   new_pf2[2] = cnr_id;
                }
                
              // 5. The new triangles are reassigned to the plist.
              //    The elements pfverts are re-assigned according to the new connectivity.
              for ( i=0; i<(*eit1).second.size(); i++ )
                {
                   (*eit1).second[i]     = new_pl1[i];
                   (*pfit).second[i]     = new_pf1[i];
                   (*eit2).second[i]     = new_pl2[i];
                   (*neigh_it).second[i] = new_pf2[i];
                }
            
              // 6. If element property values are present, the properties are averaged and mapped 
              //    back onto the elements. If both property values were the same, nothing is done.
              if ( evalues.size() == plist.size() )
                {
                   // reading property values 
                   eprop1 = evalues[ (*eit1).first-1 ];
                   eprop2 = evalues[ (*eit2).first-1 ];
                   // only if these values differ, something is done
                   if ( eprop1 != eprop2 )
                     {
                        new_prop = (eprop1 + eprop2) / 2.0;
                        evalues[ (*eit1).first-1 ] = evalues[ (*eit2).first-1 ] = new_prop;
                        csmp_error.notice( INFO, "TRIANGLE_Interface::SplitSingleCornerElements", 
                                       "Element property values were averaged when corner element was flipped");
                        cout <<"\nElement: "<< (*eit1).first <<", original property value: "<< eprop1;
                        cout <<"\nElement: "<< (*eit2).first <<", original property value: "<< eprop2;
                        cout <<"\nNew property value for both new elements: "<< new_prop << endl;
                     }
                }
              
              // 7. Updating the neighbors of the flipped elements
              // -------------------------------------------------
              // only one element must be updated with regard to its neighbors
              if ( cnr_node == 0 ) {
                  assert( new_pf1[0] > 0 );
                  findit = pfverts.find( static_cast<size_t>(new_pf1[0]) );
                  assert( findit != pfverts.end() );
                }
              else if ( cnr_node == 1 ) {
                  assert( new_pf1[1] > 0 );
                  findit = pfverts.find( static_cast<size_t>(new_pf1[1]) );
                  assert( findit != pfverts.end() );
                }
              else if ( cnr_node == 2 ) {
                  assert( new_pf1[2] > 0 );
                  findit = pfverts.find( static_cast<size_t>(new_pf1[2]) );
                  assert( findit != pfverts.end() );
                }

              for ( size_t n = 0; n < (*findit).second.size(); n++ )
                if ( (*findit).second[n] == neighbor_id ) {
                     (*findit).second[n] = cnr_id;
                     break;
                  }
              
              corner_elements++;
              
           } // end if corner element
        
      } // end looping over the elements 
      
    if ( corner_elements > 0 ) {
         cout <<"\nTRIANGLE_Interface::SplitSingleCornerElements: Split ";
         cout << corner_elements <<" corner elements."<< endl;
      }

 } // end SplitSingleCornerElements





/**

 function finds out whether the face of a triangle, the coordinates of
 which are provided as arguments, faces to the TOP, BOTTOM, LEFT or
 or right in an x(horizontal increasing to the right) y(vertical, increasing
 to the left) coordinate system.
 if the boundary is 45o the function will fail, return NOT and a cout
 message.

*/
int32_t  TRIANGLE_Interface::FindFaceBoundary( vector<double>& face_node1,
                                             vector<double>& face_node2,
                                             vector<double>& opposite_node,
                                             bool verbose )
 {
    int32_t  result(0); // = NOT
    bool   left, right, top, bottom;
    
    // getting coordinate average of boundary segment
    const double avg_x = (face_node1[0]+face_node2[0])/2.0;
    const double avg_y = (face_node1[1]+face_node2[1])/2.0;

    // finding the x and y distance between segment=face center and the 
    // opposite node
    const double dx = fabs( face_node1[0] - face_node2[0] );
    const double dy = fabs( face_node1[1] - face_node2[1] );

    // evaluating face position in triangle
    left=right=top=bottom = false;
    if ( avg_y < opposite_node[1] ) bottom = true;
    if ( avg_y > opposite_node[1] ) top    = true; // else false already set
    if ( avg_x < opposite_node[0] ) left   = true;
    if ( avg_x > opposite_node[0] ) right  = true; // else false already set

    // discriminating
    if      ( top    && dy < dx ) result = TOP_OUTSIDE;
    else if ( bottom && dy < dx ) result = BOTTOM_OUTSIDE;
    else if ( left   && dy > dx ) result = LEFT_OUTSIDE;
    else if ( right  && dy > dx ) result = RIGHT_OUTSIDE;    
    else {
       cout <<"\nMeshInterface::FindFaceBoundary: Could not identify side to which triangle boundary faces"<< endl;
       cout <<"\tcoordinates of boundary nodes (x,y): "<< endl;
       cout <<"\tnode1: "<< face_node1[0] <<", "<< face_node1[1] <<", node2: "<< face_node2[0] <<", "<< face_node2[1];
       cout <<"\n\tinside (opposite) node: ";
       cout << opposite_node[0] <<", "<< opposite_node[1] << endl;
     }
    // writing out results for testing
    if ( verbose )
      switch( result )
       {
          case TOP_OUTSIDE:    cout <<"\nresult: TOP_OUTSIDE"<< endl;
            break;
          case BOTTOM_OUTSIDE: cout <<"\nresult: BOTTOM_OUTSIDE"<< endl;
            break;
          case LEFT_OUTSIDE:   cout <<"\nresult: LEFT_OUTSIDE"<< endl;
            break;
          case RIGHT_OUTSIDE:  cout <<"\nresult: RIGHT_OUTSIDE"<< endl;
       }
       
   return result;
   
 } // end





/**

Using the bflags map, the min-max x, y, z values are found in the Vdata
and, assuming, that the model is brick-shaped, the corners are flagged
such that the 2D model has, counting from the bottom up in counter-
clockwise fashion, the corners 1, 2, 3, 4. In a 3D model another rectangle
further front (z-direction) has the corresponding corners 4, 5, 6, 7, 8.  
 */
void TRIANGLE_Interface::FlagCornerNodes( vector<std::int8_t>& bflags,
                                          deque<double>& x,
                                          deque<double>& y,
                                          deque<double>& z )
 {
     deque<double>::const_iterator xmin = min_element( x.begin(), x.end() ), 
                                     xmax = max_element( x.begin(), x.end() ),
                                     ymin = min_element( y.begin(), y.end() ),
                                     ymax = max_element( y.begin(), y.end() ),
                                     zmin = min_element( z.begin(), z.end() ),
                                     zmax = max_element( z.begin(), z.end() );
      
     size_t n_node(0U);
     for ( auto it=bflags.begin(); it!=bflags.end(); it++, n_node++ )
       {
          if ( x[ n_node ] == *xmin && y[ n_node ] == *ymin && z[ n_node ] == *zmin )
            (*it) = CNR_MIN;
          if ( x[ n_node ] == *xmax && y[ n_node ] == *ymax && z[ n_node ] == *zmax )
            (*it) = CNR_MAX;
          if ( x[ n_node ] == *xmax && y[ n_node ] == *ymin && z[ n_node ] == *zmin )
            (*it) = CNR_MIN_MAXX;
          if ( x[ n_node ] == *xmax && y[ n_node ] == *ymin && z[ n_node ] == *zmax )
            (*it) = CNR_MIN_MAXXZ;
          if ( x[ n_node ] == *xmin && y[ n_node ] == *ymin && z[ n_node ] == *zmax )
            (*it) = CNR_MIN_MAXZ;
          if ( x[ n_node ] == *xmin && y[ n_node ] == *ymax && z[ n_node ] == *zmin )
            (*it) = CNR_MAX_MINXZ;
          if ( x[ n_node ] == *xmax && y[ n_node ] == *ymax && z[ n_node ] == *zmin )
            (*it) = CNR_MAX_MAXX;
          if ( x[ n_node ] == *xmin && y[ n_node ] == *ymax && z[ n_node ] == *zmax )
            (*it) = CNR_MAX_MAXZ;
       }    
 
 } // end FlagCornerNodes  





void TRIANGLE_Interface::ListZeroPropertyValueElements( vector<double>& evalues,
                                                        map<size_t,vector<int64_t> >& plist,
                                                        deque<double>& x, deque<double>& y ) const
 {
    double avgx, avgy;
 
    cout <<"\nTRIANGLE_Interface::ListZeroPropertyValueElements: Searching for elements with zero ";
    cout <<"'permeability' values..."<< endl;
    for ( size_t i=0U; i<evalues.size(); i++ )
      if ( evalues[i] == 0. )
        {
           map<size_t,vector<int64_t> >::const_iterator it=plist.find(i+1U);
           assert( it != plist.end() );
           avgx = avgy = 0.0;
           for ( size_t j=0; j<(*it).second.size(); j++ ) {
                avgx += x[ (*it).second[j]-1 ];
                avgy += y[ (*it).second[j]-1 ];
             }
           avgx /= (*it).second.size();
           avgy /= (*it).second.size();
           cout <<"\n\t\tfound Element "<< i+1 <<" with coordinates: ";
           cout << avgx <<", "<< avgy <<", " << endl << endl;
        }
        
 } // end ListZeroPropertyValueElements 


template void TRIANGLE_Interface::ReadTriangle2DMesh<2U>( const char*, 
                                                          VSet<2U>&,
                                                          bool, bool );

template void TRIANGLE_Interface::ReadTriangle2DMeshAndCreateDiscreteFractures<2U>( const char*, 
                                                                                    VSet<2U>&,
                                                                                    double, bool );


} // end namespace csmp
