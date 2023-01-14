#include "GoCadInterface.h"
#include "InterFace.h"
#include "Region.h"
#include "Model.h"
#include "HashKey.h"
#include "LinearTriangle3D.h"
#include "LinearTetrahedron.h"
#include "MJL_Triangle3D.h"
#include "Exception.h"
#include "PropertyData.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
GoCadInterface<dim>::GoCadInterface()
 {
 }


template<uint32_t dim>
GoCadInterface<dim>::~GoCadInterface()
 {
 }


/** Peaks into next file of file to see whether it contains the search string
of interest. 
*/
template<uint32_t dim>
bool  GoCadInterface<dim>::IsInNextLine( ifstream& ifs, const char* search_string ) const
 {
    char text_line[256];
   
    // 1. remembering stream positions
    // -------------------------------
    streampos pos = ifs.tellg();
    ifs.get( text_line, 256 );
    ifs.seekg( pos );

    // seeking for word in text string
    if ( strstr( text_line, search_string ) == NULL ) return false;
    return true;
}


template<uint32_t dim>
BOX_BOUNDARY  GoCadInterface<dim>::IdentifyTetrahedronBoundary( const vector<double>& nd1,
                                                                      const vector<double>& nd2,
                                                                      const vector<double>& nd3 )
 {
    assert( nd1.size() == 3U );
    assert( nd2.size() == 3U );
    assert( nd3.size() == 3U );

    mjl::Triangle3D  face( (mjl::Point3D( nd1[0],nd1[1],nd1[2] )),
                           (mjl::Point3D( nd2[0],nd2[1],nd2[2] )),
                           (mjl::Point3D( nd3[0],nd3[1],nd3[2] )), 1 );

    int32_t  facing_direction = face.FacingDirection();

    if ( facing_direction == MJL3D_BOTTOM ) return BOTTOM;
    if ( facing_direction == MJL3D_FRONT )  return FRONT;
    if ( facing_direction == MJL3D_RIGHT )  return RIGHT;
    if ( facing_direction == MJL3D_BACK )   return BACK;
    if ( facing_direction == MJL3D_LEFT )   return LEFT;
    if ( facing_direction == MJL3D_TOP )    return TOP;

    cout <<"\nGoCadInterface::IdentifyTetrahedronBoundary: facing direction was not found."<< endl;
    cout <<"Coordinates of face (triangle) nodes: "<< endl;
    for ( size_t i{0U}; i<nd1.size(); i++ ) cout << nd1[i] <<" ";
    cout << endl;
    for ( size_t q=0; q<nd2.size(); q++ ) cout << nd2[q] <<" ";
    cout << endl;
    for ( size_t r=0; r<nd3.size(); r++ ) cout << nd3[r] <<" ";
    cout << endl;

    if ( facing_direction == MJL3D_TILTED ) {
         cout <<"\nGoCadInterface::IdentifyTetrahedronBoundary: Assigned IRREGULAR flag."<< endl;
         return IRREGULAR;
      }

    return NOT; // NOT = 0

 } // end IdentifyTetrahedronBoundary




/**
 
Reads the GoCad header and property class headers from the input text
file created by GoCad. The method subsequently calls ReadTSurface() 
which reads the actual data from the input file and puts them into
the supplied VSet, including property data if present.  

@section arguments Input Arguments 

The first argument specifies the input file name without (.ts) extension
which is appended automatically.  

@param header The second method argument is used to store the header information
which is retrieved from the input file, 

@param vset the third argument will hold
the mesh which was read, including the property data if such were 
present.  

@section application Application

To read GoCad 'TSurf' surfaces with assigned properties from textfiles
output by GoCad.  

@section messages Messages 

The method reports the information given in the header section of the 
GoCad file. This includes the header with the object type and name, as
well as the property class headers for the properties associated with
the surface.  
 */
template<uint32_t dim>                        
void GoCadInterface<dim>::ReadTetrahedralGocad3DSurface( const char* fname, 
                                                            GocadHeader& header, 
                                                            VSet<dim>& vset )
 {
    bool  read_data, debug(true);
    char  file_name[200];
    
    strcpy( file_name, fname );
    strcat( file_name, ".ts");

    ifstream ifs( file_name );
    if ( !ifs ) {
         cout <<"\nGoCadInterface::ReadTetrahedralGocad3DSurface: File not found...";
         return;
      }

    // 0. variables
    // ------------
    char  text_line[256];
    char* token;
    const char* const delims = " ,=,\n,\r,\t"; 
    GocadPropertyClassHeader                  dummy;
    list<GocadPropertyClassHeader>            prop_headers;
    typename list<GocadPropertyClassHeader>::iterator  it;
   
    // 1. reading first line and checking whether a gocad file is read
    // ---------------------------------------------------------------
    ifs.getline( text_line, 256 );

    // seeking for word in text string
    if ( strstr( text_line, "GOCAD" ) == NULL ) {
         cout <<"\nGoCadInterface::ReadTetrahedralGocad3DSurface: not a GoCad input file. Now exciting..."<< endl;
         return;
      }
   
    // 2. Identifying type of GoCad object
    // -----------------------------------
    if ( strstr( text_line, "TSurf" ) != NULL ) {
         cout <<"\nGoCadInterface::ReadTetrahedralGocad3DSurface: reading 'TSurf' object..."<< endl;
         // 2.1 reading header into GocadHeader
         header.InitializeFrom( ifs );
         header.Out();
      }
    streampos pos = ifs.tellg(); // remembering position in file stream

    // 3. If necessary, search for PROPERTIES information block
    // --------------------------------------------------------
    size_t tries(0), max_tries(1000);
    if ( strstr( text_line, "PROPERTIES" ) == NULL )
      do {
            ifs.getline( text_line, 256 );
            tries++;
         }
    while ( strstr( text_line, "PROPERTIES" ) == NULL && tries <= max_tries );

    if ( tries >= max_tries )     
      throw csmp::Exception( ERROR, "GoCadInterface::ReadTetrahedralGocad3DSurface", 
                             "PROPERTIES specifier not found" );

    // 3.1 reading property class header
    // ---------------------------------
    if ( strstr( text_line, "PROPERTIES" ) != NULL )
      {
         // this is expected
         // ---------------------------------------------------------------------------------
         // PROPERTIES prop1 prop2 ...                   property name list
         // NO_DATA_VALUES 1e-16 1e-14 ...               values used as init.vals. for perm
         // PROPERTY_CLASSES permeability  thickness ... 
         // ESIZES 1 1 ...                               type of prop (scalar=1, vector=2...)
         // ---------------------------------------------------------------------------------
         // parsing the PROPERTIES line, reading the properties
         token = strtok( text_line, delims ); // PROPERTIES string
         token = strtok( NULL, delims );      // first property name string
             
         // reading all property names, creating a class header for each
         do prop_headers.push_back( GocadPropertyClassHeader(token,0,1,0) );
         while ( (token=strtok( NULL, delims )) != NULL );

         // reading NO_DATA_VALUES for each property
         ifs.getline( text_line, 256 );
         assert( strstr( text_line, "NO_DATA_VALUES" ) != NULL );
         token = strtok( text_line, delims ); 
         for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ ) {
              assert( (token=strtok( NULL, delims )) != NULL );
              (*it).NoDataValue( atof(token) );
           } 
             
         // reading PROPERTY_CLASSES for each property 
         ifs.getline( text_line, 256 );
         assert( strstr( text_line, "PROPERTY_CLASSES" ) != NULL );
         token = strtok( text_line, delims ); 
         for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ ) {
              assert( (token=strtok( NULL, delims )) != NULL );
              (*it).PropertyClass( token );
           } 
             
         // reading ESIZES information for each property
         ifs.getline( text_line, 256 );  
         assert( strstr( text_line, "ESIZES" ) != NULL );
         token = strtok( text_line, delims ); 
         for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ ) {
              assert( (token=strtok( NULL, delims )) != NULL );
              (*it).ESize( static_cast<uint32_t>(atoi(token)) );
           } 
             
         // 3.2 Reading the property class headers for each property
         // --------------------------------------------------------
         if ( !IsInNextLine( ifs, "PROPERTY_CLASS_HEADER" ) )
           throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTetrahedralGocad3DSurface", 
                                        "First PROPERTY_CLASS_HEADER not found" );
         else
           {
              it = prop_headers.begin();
              
              while ( it != prop_headers.end() ) {
                   // Since GoCad tends to list other properties which have no nodal values
                   // associated with the PVTRX lines, only the ones listed above in the 
                   // property class header above must be registered.
                   if ( !IsInNextLine( ifs, (*it).Property().c_str() ) ) 
                     dummy.InitializeFrom( ifs );
                   else {
                           (*it).InitializeFrom( ifs );
                           it++;
                        }
                }
           }
           
         if ( debug )
           for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ )
             (*it).Out();

         // 3.3 reading solid including properties
         // --------------------------------------
         read_data = ReadTSurface( ifs, prop_headers, vset );           
            
      } // end TSurf with property loop
         
     
    // 4. reading solid description only
    // ---------------------------------
    else 
      {
         cout <<"\nGoCadInterface::ReadTetrahedralGocad3DSurface: No property information found.";
         cout <<" Reading geometry only..."<< endl;
         
         // return to previous position in the file
         ifs.seekg( pos );
         read_data = ReadTSurface( ifs, prop_headers, vset );
      }
      
    // 5. Closing the output file
    // --------------------------
    ifs.close();
     
    if ( read_data ) {  
         cout <<"\nGoCadInterface::ReadTetrahedralGocad3DSurface: file '"<< file_name;
         cout <<"' read successfully..." << endl;
      }
   
} // end ReadTetrahedralGocad3DSurface








/**
 
The method reads a single TSurf from the input file stream. Associated
properties in as much as specified in the supplied property list are
red and assigned as node data to the VSet.  

@section arguments Input Arguments 

The input file stream in reading mode, the list of property class 
headers initialized with the header of the input file and the VSet
in which the read data are stored.  

@return The surface geometry and the associated properties are returned into
third method argument. If the reading of the TSurf information is 
performed successfully the boolean variable true is returned.  

@section messages Messages 

The method reports potential file reading errors.  
*/
template<uint32_t dim>                        
bool GoCadInterface<dim>::ReadTSurface( ifstream& ifs, 
                                           const list<GocadPropertyClassHeader>& properties,
                                           VSet<dim>& vset )
 {
    char               text_line[256];
    char              *token;
    const char* const  delims = " ,=,\n,\r,\t"; 
    bool               debug(false);
    
    if ( !ifs )
      throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTSurface:", 
                                   "Invalid file input stream",
                                   "Reading could not be performed.");
   
    // 1. reading first line and checking whether the right data will follow
    // ---------------------------------------------------------------------
    streampos pos = ifs.tellg();
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims ); 

    if ( strstr( text_line, "TFACE" ) == NULL ) {
         ifs.seekg( pos );
         throw csmp::Exception( ERROR, "GoCadInterface::ReadTSurface","TFACE descriptor missing");
         return false;
      }
      
    // 2. Reading node data VRTX or PVRTX
    // ---------------------------------
    typename list<GocadPropertyClassHeader>::const_iterator  propit;
    size_t                                       data_entries(0);
    
    for ( propit=properties.begin(); propit!=properties.end(); propit++ )
      data_entries += (*propit).ESize();
      
    //   id      attr   node xyz
    map<size_t,pair<int32_t,vector<double> > >                  nodes;
    typename map<size_t,pair<int32_t,vector<double> > >::const_iterator  nit;
    map<size_t,vector<double> >                             node_properties;
    typename map<size_t,vector<double> >::const_iterator             nprop_it;
    vector<double>                                                   props(data_entries);
    vector<double>          xyz(3);
    size_t                  nodeID(UNSPECIFIED);
    int32_t                 attr;
    //   attr       value
    pair<int,vector<double> >  node;

    while ( ifs.getline( text_line, 256 ) )
      {
         // tokenizing text line
         token = strtok( text_line, delims ); 
         if ( token != NULL && (!strcmp( token, "TRGL" ) ||
                                !strcmp( token, "PROPERTY_CN" )) ) break;
         // getting line
         if ( token != NULL && (!strcmp( token, "VRTX" ) || !strcmp( token, "PVRTX" )) ) 
           {
             // getting point coordinates
             token = strtok( NULL, delims );
             nodeID = static_cast<uint32_t>(atol( token ));   // node ID
             token = strtok( NULL, delims );
             xyz[0] = atof( token );   // x
             token = strtok( NULL, delims ); 
             xyz[1] = atof( token );   // y         
             token = strtok( NULL, delims );
             xyz[2] = atof( token );   // z
             // node properties
             if ( !properties.empty() ) {
                  for ( auto n=0; n<data_entries; n++ ) {
                       token    = strtok( NULL, delims );
                       props[n] = atof( token ); // property value
                    }
                  node_properties[ nodeID ] = props;
               }
             // constraint information
             token = strtok( NULL, delims ); 
             if ( token != NULL && !strcmp( token, "CNXYZ" ) )  attr = 1;  // boundary node
             else                                               attr = 0;
             node.first  = attr;
             node.second = xyz;
             
             // assigning node to node list (1...n)
             nodes[ nodeID ] = node;
           } 
      }
    assert( nodeID > 0 && nodeID <= nodes.size() );


   // 3. Reading but ignoring constraint point flags
   // ----------------------------------------------
   if ( token != NULL && !strcmp( token, "PROPERTY_CN" ) ) 
     while ( ifs.getline( text_line, 256 ) )
       {
          // tokenizing first word of text line
          token = strtok( text_line, delims ); 
          if ( token != NULL && ( !strcmp( token,"END") || !strcmp( token,"END_PROPERTY_CN") ) )
            {
                // next line is needed for function below
                ifs.getline( text_line, 256 );
                token = strtok( text_line, delims ); 
                break;
            }
          // abnormal termination
          if ( token != NULL && strcmp( token, "CNP" ) )
            throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTSurface:", 
                                         "Abnormal string read", token );
       }

    // 4. Reading triangle elements TRGL
    // ---------------------------------
    map<size_t,vector<int64_t> >  triangles;
    vector<int64_t>               node_ids(3);
    size_t                       triangleID(1);

    do
      {
         if ( triangleID != 1 )
           {
              token = strtok( text_line, delims ); 
              if ( token != NULL && !strcmp( token, "BSTONE" ) ) break;
              if ( token != NULL && !strcmp( token, "BORDER" ) ) break;
              if ( token != NULL && !strcmp( token, "END" ) )    break;
              if ( token != NULL && !strcmp( token, "TFACE" ) )  break;
           }
         if ( token != NULL && !strcmp( token, "TRGL" ) ) 
           {
             // getting node ID numbers for each triangle
             token = strtok( NULL, delims );
             node_ids[0] = static_cast<uint32_t>(atol( token ));   // node 1
             token = strtok( NULL, delims ); 
             node_ids[1] = static_cast<uint32_t>(atol( token ));   // node 2        
             token = strtok( NULL, delims );
             node_ids[2] = static_cast<uint32_t>(atol( token ));   // node 3
             
             // assigning nodes to triangles and triangles to tface
             triangles[ triangleID++ ] = node_ids;
           } 
      }
    while ( ifs.getline( text_line, 256 ) );

    cout <<"\nGoCadInterface::ReadTSurface: ";
    cout <<"completed reading first TFACE; Reading is stopped."<< endl;

    throw csmp::Exception( WARNING, "GoCadInterface::ReadTSurface:", 
                   "If there is more than 1 surface in this file, the latter are ignored" );


    // 4. Testing whether input was read correctly
    // -------------------------------------------
    typename map<size_t,vector<int64_t> >::iterator  plit;
    size_t n;
    
    if ( !properties.empty() ) assert( node_properties.size() == nodes.size() );
    
    if ( debug )
      {
         cout <<"\nNodes ("<< nodes.size() <<") and their coordinates: "<< endl;
         cout <<"Printing first 100..."<< endl;
         for ( n=0, nit=nodes.begin(); nit!=nodes.end(); nit++ )
           {
              if ( n++ >= 100 ) break;
              cout <<"Node: "<< (*nit).first <<":  ";
              for ( auto i{0U}; i<(*nit).second.second.size(); i++ )
              cout << (*nit).second.second[i] <<",\t";
              cout <<"flag: "<< (*nit).second.first << endl;
           }
    
         cout <<"\nTriangles ("<< triangles.size() <<") and their member Nodes: "<< endl;
         cout <<"Printing first 100..."<< endl;
         for ( n=0, plit=triangles.begin(); plit!=triangles.end(); plit++ )
           {
              if ( n++ >= 100 ) break;
              cout <<"Triangle: "<< (*plit).first <<":  ";
              for ( auto i{0U}; i<(*plit).second.size(); i++ )
                cout << (*plit).second[i] <<", ";
              cout << endl;
           }
      }

    // 4b. Checking for consecutive node numbering
    // -------------------------------------------
    if ( !VerifyConsecutiveNodeNumbering( triangles ) )
      throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTSurface", 
            "Nodes in input mesh are not numbered consecutivly");


    // ----------------------------------------------------------------------------------
    // 5. setup VSet
    // ---------------------------------------------------------------------------------- 
    LinearTriangle3D  triangle;
    vset.Resize( triangle.Nodes(),
                 triangle.Neighbors(),
                 triangle.ElementType(), 
                 nodes.size(), triangles.size() );

    // node coordinates
    // NOTE: Y and Z are flipped, since Z points upward in GoCad as 
    // opposed to CSMP, where Z points to the front.
    // ------------------------------------------------------------
    for ( nit=nodes.begin(); nit!=nodes.end(); nit++ )
      {
         vset.Px( (*nit).first, (*nit).second.second[1] ); // Y becomes X
         vset.Py( (*nit).first, (*nit).second.second[2] ); // Z becomes Y
         vset.Pz( (*nit).first, (*nit).second.second[0] ); // X becomes Z
       }
    nodes.erase( nodes.begin(), nodes.end() );

    // plist
    // -----
    vset.AddPlist( triangles.begin(), triangles.end() );  


   // -----------------------------------------------------
   // 6. FINDING THE NEIGHBOR ELEMENTS (opposite each node)
   // -----------------------------------------------------
   cout <<"\nGoCadInterface::ReadTSurface: building 'pfvert' (neighbor element) dataset..."<< endl;
   // 3.1 Making a list of all the faces of each triangular element:
   //     - the faces are located opposite each node
   //     - the nodes of each face are listed in counterclockwise fashion
   //    face_key  e-id
   pair<string,size_t>                              face;
   multimap<string,size_t>            face_tree;
   typename multimap<string,size_t>::iterator  face_it; 
   HashKey                                             hasher; 
   
   for( plit=triangles.begin(); plit!=triangles.end(); plit++ )
    {
       // making unique search keys from the node ids of the face nodes
       // (convention: face 1 lies opposite of node 1, face 2 opposite of node 2 etc.)
       // ------
       // face 1
       // ------
       hasher.Key( (*plit).second[1], (*plit).second[2], face.first );
       face.second = (*plit).first;
       face_tree.insert( make_pair(face.first,face.second) );
       // ------
       // face 2
       // ------
       hasher.Key( (*plit).second[2], (*plit).second[0], face.first );
       face.second = (*plit).first;
       face_tree.insert( make_pair(face.first,face.second) );
      // ------
       // face 3
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[1], face.first );
       face.second = (*plit).first;
       face_tree.insert( make_pair(face.first,face.second) );
    }

   // --------------------------------------------------------------------------------------
   // 7. BUILDING Pfverts() by searching for neighbors in hash key map.
   //    The boundary nodes at corresponding faces are flagged exactly as the 
   //    as the faceverts.
   // --------------------------------------------------------------------------------------
   deque<vector<int64_t> > pfverts( triangles.size() );
   vector<int64_t>         pfvert(  triangle.Neighbors() );
   vector<int8_t>         pbflags( vset.Vertices(), 0 );
   vector<double>       nd1(3), nd2(3), nd3(3);
   char                   face_key[30];

   cout <<"\nGoCadInterface::ReadTSurface: Identifying the neighbors of each element..." << endl;

   for( plit=triangles.begin(); plit!=triangles.end(); plit++ )
    {
       // 1. identify the coordinates of the four nodes of the element
       // ------------------------------------------------------------
       // node 1 (put into else statement as only needed there)
       nd1[0] = vset.Px( (*plit).second[0] ); 
       nd1[1] = vset.Py( (*plit).second[0] );
       nd1[2] = vset.Pz( (*plit).second[0] );
       // node 2
       nd2[0] = vset.Px( (*plit).second[1] ); 
       nd2[1] = vset.Py( (*plit).second[1] );
       nd2[2] = vset.Pz( (*plit).second[1] );
       // node 3
       nd3[0] = vset.Px( (*plit).second[2] ); 
       nd3[1] = vset.Py( (*plit).second[2] );
       nd3[2] = vset.Pz( (*plit).second[2] );
  
       if ( debug )
         {
             cout <<"\nTriangular element "<< (*plit).first <<" node coordinates:"<< endl;
             cout <<"Node 1 (x,y,z): "<< nd1[0] <<" "<< nd1[1] <<" "<< nd1[2] << endl;
             cout <<"Node 2 (x,y,z): "<< nd2[0] <<" "<< nd2[1] <<" "<< nd2[2] << endl;
             cout <<"Node 3 (x,y,z): "<< nd3[0] <<" "<< nd3[1] <<" "<< nd3[2] << endl;
         }
       // ------
       // face 1
       // ------
       hasher.Key( (*plit).second[1], (*plit).second[2], face_key );
       // finding face with the same key and checking to which element
       // it belongs
       if ( (face_it=face_tree.find( face_key )) == face_tree.end() )
         throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTSurface:", 
                                      "Failed to retrieve face hash-key from multimap for element", face_key );

       // pre-empt effects of probable bug in CW Pro 5 STL: multimap should always return iterator
       // to first element with key in the map, here it doesn't
       // -----------------------------------------------------
       // wind back one, but only if previous element has same key
       if ( face_it!=face_tree.begin() ) if ( (*(--face_it)).first != face_key ) face_it++;
         
       // if the element ID associated with the key is different from the
       // current element, the neighbor has already been found
       // ----------------------------------------------------
       if ( debug ) cout <<"\nTest1 face1: Element ID recovered for key: "<< face_key;
       if ( debug ) cout <<" from map vs. plist iterator elmt ID: "<< (*face_it).second <<" "<< (*plit).first << endl;
       if ( (*face_it).second != (*plit).first ) pfvert[0] = ((*face_it).second);
       // if the face belongs to the same element one searches for
       // the next occurrence of the key in the multimap
       // ----------------------------------------------
       else
       // advance in the multimap
       // -----------------------
         {
            face_it++;
            if ( (*face_it).first == face_key ) pfvert[0] = ((*face_it).second);
            else {
                if ( debug ) cout <<"\nElement: "<< (*plit).first <<" Test face 1: failed comparison: list-face key: ";
                if ( debug ) cout <<(*face_it).first <<" vs. hash key: "<< face_key << endl;

                // If there is no shared occurrence of the face, the face lies at the model
                // boundary. If so, the type of this boundary is identified.
                // ---------------------------------------------------------
                pfvert[0] = IRREGULAR;
                  
                // flagging the boudary nodes as such
                // ----------------------------------
                assert( pfvert[0] < 0 );
                pbflags[ (*plit).second[1] ] = static_cast<int8_t>(pfvert[0]);
                pbflags[ (*plit).second[2] ] = static_cast<int8_t>(pfvert[0]);
             }
         }
       // ------
       // face 2
       // ------
       hasher.Key( (*plit).second[2], (*plit).second[0], face_key );
       if ( (face_it=face_tree.find( face_key )) == face_tree.end() )
         throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTSurface:", 
                                      "Failed to retrieve face hash-key from multimap for element", face_key );

       if ( face_it!=face_tree.begin() ) if ( (*(--face_it)).first != face_key ) face_it++;
       if ( debug ) cout <<"\nTest1 face2: Element ID recovered for key: "<< face_key;
       if ( debug ) cout <<" from map vs. plist iterator elmt ID: "<< (*face_it).second <<" "<< (*plit).first << endl;
       if ( (*face_it).second != (*plit).first ) pfvert[1] = ((*face_it).second);
       else
         {
            face_it++;
            if ( (*face_it).first == face_key ) pfvert[1] = (*face_it).second;
            else {
                if ( debug ) cout <<"\nElement: "<< (*plit).first <<" Test face 2: failed comparison: list-face key: ";
                if ( debug ) cout <<(*face_it).first <<" vs. hash key: "<< face_key << endl;
                pfvert[1] = IRREGULAR;
                pbflags[ (*plit).second[0] ] = static_cast<int8_t>(pfvert[1]);
                pbflags[ (*plit).second[2] ] = static_cast<int8_t>(pfvert[1]);
             }
         }
       // ------
       // face 3
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[1], face_key );
       if ( (face_it=face_tree.find( face_key )) == face_tree.end() )
         throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTSurface:", 
                                      "Failed to retrieve face hash-key from multimap for element", face_key );

       if ( face_it!=face_tree.begin() ) if ( (*(--face_it)).first != face_key ) face_it++;
       if ( debug ) cout <<"\nTest1 face3: Element ID recovered for key: "<< face_key;
       if ( debug ) cout <<" from map vs. plist iterator elmt ID: "<< (*face_it).second <<" "<< (*plit).first << endl;
       if ( (*face_it).second != (*plit).first ) pfvert[2] = ((*face_it).second);
       else
         {
            face_it++;
            if ( (*face_it).first == face_key ) pfvert[2] = ((*face_it).second);
            else {
                if ( debug ) cout <<"\nElement: "<< (*plit).first <<" Test face 3: failed comparison: list-face key: ";
                if ( debug ) cout <<(*face_it).first <<" vs. hash key: "<< face_key << endl;
                pfvert[2] = IRREGULAR;
                pbflags[ (*plit).second[0] ] = static_cast<int8_t>(pfvert[2]);
                pbflags[ (*plit).second[1] ] = static_cast<int8_t>(pfvert[2]);
             }
         }
       pfverts[ (*plit).first ] = pfvert;
    }

  if ( debug )
    {
       cout <<"\nDetermined boundary flags:"<< endl;
       size_t counter(0U);
       for ( auto bflit=pbflags.begin(); bflit!=pbflags.end(); bflit++, counter++ )
         if ( (*bflit) < 0 )
           cout <<"\nNode: "<< counter <<" flagged: "<< (*bflit);
           
       cout << endl; 
    }
    
  // -------------------------------------------
  // 8. Writing remaining results to the VSet 
  // -------------------------------------------
  cout <<"\nGoCadInterface::ReadTSurface: Writing remaining finite-element data to VSet..." << endl;
  vset.AddPfverts( pfverts.begin(), pfverts.end() );
  vset.AddBFlags( pbflags.begin(), pbflags.end() );

  // -------------------------------------------
  // 9. Adding property information to Vset
  // -------------------------------------------
  // only nodal data can be transferred by GoCad
  // looping over the property storage and assigning the data to CSMP
  size_t data_entry(0);
  
  cout <<"\nGoCadInterface::ReadTSurface: Writing property data to VSet..." << endl;
  for ( propit=properties.begin(); propit!=properties.end(); propit++ )
    {
       cout <<"\n\tWriting property: '"<< (*propit).Property() <<"'"<< endl;
       // scalar node properties
       if ( (*propit).ESize() == 1 ) {
            csmp::Index  setting(SCALAR,NODE,0);
            PropertyData  data( setting.place, setting.type, dim );
            data.Reserve( vset.Vertices() );
            // map<size_t,vector<double> >::const_iterator
            for ( auto nprop_it2=node_properties.begin(); nprop_it2!=node_properties.end(); nprop_it2++ )
              pushBack( data, makeScalar( PLAIN, (*nprop_it2).second[ data_entry ] ) );
            // adding scalars to the VSet
            vset.AddData( (*propit).Property().c_str(), data );
            data_entry++;
         }
       // vector<double> node properties (3 components in 3D)
       if ( (*propit).ESize() == 3 ) {
            csmp::Index  setting(VECTOR,NODE,0);
            PropertyData  data( setting.place, setting.type, dim );
            data.Reserve( vset.Vertices() );
            VectorVariable<dim> vc;
            for ( nprop_it=node_properties.begin(); nprop_it!=node_properties.end(); nprop_it++ ) {
                 for ( size_t i{0U}; i<3; i++ ) vc(i) = (*nprop_it).second[ data_entry+i ];
                 pushBack( data, vc );
              }
            // adding scalars to the VSet
            vset.AddData( (*propit).Property().c_str(), data );
            data_entry += 3;
         }
       // tensor node properties (9 components in 3D)
       if ( (*propit).ESize() == 9 ) {
            csmp::Index  setting(TENSOR,NODE,0);
            PropertyData  data( setting.place, setting.type, dim );
            data.Reserve( vset.Vertices() );
            TensorVariable<dim> ts;
            for ( nprop_it=node_properties.begin(); nprop_it!=node_properties.end(); nprop_it++ ) {
                 for ( size_t i{0U}; i<3; i++ )
                   for ( size_t j{0U}; j<3; j++ ) ts(i,j) = (*nprop_it).second[ data_entry+i*3+j ];
                 pushBack( data, ts );
              }
            // adding scalars to the VSet
            vset.AddData( (*propit).Property().c_str(), data );
            data_entry += 9;
         }
    }

  if ( debug ) vset.Out();
   
  return true;
    
} // end GoCadInterface::ReadTSurface











/**
 
Top level method to read ASCII text files from the geometric modeling 
program GoCad which were output in the TSolid format. The method just 
reads the header of the GoCad file diagnosing its format and interpreting
which datasets are contained. If the file header can be interpreted
correctly, control is transferred to ReadTSolid() which will read the data
from the TSolid file. 

ReadGocadText() expects either no variable or the variable 'permeability'
to be specified in the GoCad output file. At the moment, the method
only supports the reading of a single property. Thus, if the file contains
multiple properties, the reading process will fail. 

@section arguments Input Arguments 

The name of the GoCaD text, ascii file with the TSolid description, 
a GoCadHeader object which will store the information about the GoCad
dataset, and the VSet into which the information from the GoCad file
will be returned. 

@param vset The mesh and associated data from the TSolid dataset will be returned 
into the VSet supplied as third method argument. 

@section implementation Implementation

The method parses the values of the GoCad variables PROPERTIES,
PROPERTY_CLASSES, and ESIZES (dimensions of each variable) and starts 
the method ReadTSolid() with the according input arguments. 

@section application Application

The methods provides an interface to the GoCaD geometric modeling tool,
facilitating that tetrahedral meshes created by this program can be used 
in CSMP computations. 

@section messages Messages 

If the GoCad input file does not exist, if more than a single property
was specified therein, or if the variable 'permeability' is missing
the method will report an error and exit.  
 */
template<uint32_t dim>                        
void GoCadInterface<dim>::ReadTetrahedralGocad3DMesh( const char* file, 
                                                         GocadHeader& header, 
                                                         VSet<dim>& vset )
 {
    char                    file_name[200];
    bool                    debug(false);
    
    strcpy( file_name, file );
    strcat( file_name, ".so");

    ifstream ifs( file_name );
    
    if ( !ifs.is_open() )
      throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTetrahedralGocad3DMesh:",
                                   "File not found.",
                                   "Remember that extension '.so' is automatically appended" );

    // 0. variables
    // ------------
    char                                      text_line[256];
    char                                     *token;
    const char* const  delims = " ,=,\n,\r,\t"; 
    GocadPropertyClassHeader                  dummy;
    list<GocadPropertyClassHeader>            prop_headers;
    typename list<GocadPropertyClassHeader>::iterator  it;
   
    // 1. reading first line and checking whether a gocad file is read
    // ---------------------------------------------------------------
    ifs.getline( text_line, 256 );

    // seeking for word in text string
    if ( strstr( text_line, "GOCAD" ) == NULL )
      throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTetrahedralGocad3DMesh:",
                                   "GOCAD header string could not be found in input file");
   
    // 2. Identifying type of GoCad object
    // -----------------------------------
    if ( strstr( text_line, "TSolid" ) != NULL )
      {
         cout <<"\nGoCadInterface::ReadTetrahedralGocad3DMesh: reading 'TSolid' object..."<< endl;
         
         // 2.1 reading header
         header.InitializeFrom( ifs );
         header.Out();

         // 2.2 reading property class header
         //     If there is property information
         streampos pos = ifs.tellg(); // remembering position in file stream

    // 3. If necessary, search for PROPERTIES information block
    // --------------------------------------------------------
    size_t tries(0), max_tries(1000);
    if ( strstr( text_line, "PROPERTIES" ) == NULL )
      do {
            ifs.getline( text_line, 256 );
            tries++;
         }
    while ( strstr( text_line, "PROPERTIES" ) == NULL && tries <= max_tries );

    if ( tries >= max_tries )     
      throw csmp::Exception( ERROR, "GoCadInterface::ReadTetrahedralGocad3DMesh", 
                             "PROPERTIES specifier not found" );

    // 3.1 reading property class header
    // ---------------------------------
    if ( strstr( text_line, "PROPERTIES" ) != NULL )
      {
         // this is expected
         // ---------------------------------------------------------------------------------
         // PROPERTIES prop1 prop2 ...                   property name list
         // NO_DATA_VALUES 1e-16 1e-14 ...               values used as init.vals. for perm
         // PROPERTY_CLASSES permeability  thickness ... 
         // ESIZES 1 1 ...                               type of prop (scalar=1, vector=2...)
         // ---------------------------------------------------------------------------------
         // parsing the PROPERTIES line, reading the properties
         token = strtok( text_line, delims ); // PROPERTIES string
         token = strtok( NULL, delims );      // first property name string
             
         // reading all property names, creating a class header for each
         do prop_headers.push_back( GocadPropertyClassHeader(token,0,1,0) );
         while ( (token=strtok( NULL, delims )) != NULL );

         // reading NO_DATA_VALUES for each property
         ifs.getline( text_line, 256 );
         assert( strstr( text_line, "NO_DATA_VALUES" ) != NULL );
         token = strtok( text_line, delims ); 
         for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ ) {
              assert( (token=strtok( NULL, delims )) != NULL );
              (*it).NoDataValue( atof(token) );
           } 
             
         // reading PROPERTY_CLASSES for each property 
         ifs.getline( text_line, 256 );
         assert( strstr( text_line, "PROPERTY_CLASSES" ) != NULL );
         token = strtok( text_line, delims ); 
         for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ ) {
              assert( (token=strtok( NULL, delims )) != NULL );
              (*it).PropertyClass( token );
           } 
             
         // reading ESIZES information for each property
         ifs.getline( text_line, 256 );  
         assert( strstr( text_line, "ESIZES" ) != NULL );
         token = strtok( text_line, delims ); 
         for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ ) {
              token=strtok( NULL, delims );
              assert( token != NULL );
              (*it).ESize( static_cast<uint32_t>(atoi(token)) );
           } 
             
         // 3.2 Reading the property class headers for each property
         // --------------------------------------------------------
         if ( !IsInNextLine( ifs, "PROPERTY_CLASS_HEADER" ) )
           throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTetrahedralGocad3DMesh", 
                                        "First PROPERTY_CLASS_HEADER not found" );
         else
           {
              it = prop_headers.begin();
              
              while ( it != prop_headers.end() ) {
                   // Since GoCad tends to list other properties which have no nodal values
                   // associated with the PVTRX lines, only the ones listed above in the 
                   // property class header above must be registered.
                   if ( !IsInNextLine( ifs, (*it).Property().c_str() ) ) 
                     dummy.InitializeFrom( ifs );
                   else {
                           (*it).InitializeFrom( ifs );
                           it++;
                        }
                }
           }
           
         if ( debug )
           for ( it=prop_headers.begin(); it!=prop_headers.end(); it++ )
             (*it).Out();

         // reading solid including properties
         // ----------------------------------
         ReadTSolid( ifs, vset, true );
      }

     // 2.3 reading solid description only
     else 
       {
          cout <<"\nGoCadInterface::ReadTetrahedralGocad3DMesh: No property information found. ";
          cout <<" Reading geometry only..."<< endl;
          // return to previous position in the file
          ifs.seekg( pos );
          ReadTSolid( ifs, vset );
        }
      }
//    ifs.close();
    
    cout <<"\nGoCadInterface::ReadTetrahedralGocad3DMesh: file '"<< file_name;
    cout <<"' read successfully..." << endl;

} // end ReadTetrahedralGocad3DMesh





/**
 
Reads a Gocad TSolid data input file and inputs it into a VSet. In this
process the neighbors of each finite element are identified and the
boundaries of the model are flagged, not discriminating yet the corners and 
edges. This is done later when the CSP Model ist constructured in the 
method MeshManager::InitializeMeshFrom(VSet).  

Material properties assigned in GoCad to the TVOLUMES are not read by 
this method since these are nodal properties in GoCad. This leads to an
ambiguity in the TSolid file where nodes at the boundary of TVOLUMES are
shared between different volumes.  

As a remedy, ReadTSolid() assigns the TVOLUME numbers to the variable
'permeability'. Using these values, groups can be build in CSP whose
properties can then be changed to the desired value.  

Currently, the method reads the nodal permeability values but they
are not used.  

Importantly, there is a difference between the Gocad and the CSP3D 
coordinate systems. In Gocad Z points upward while Z points to the front
in CSP such that 3D is merely an extension of 2D by 1 dimensions. Adhering
to this convention, makes CSP3D compatible with OpenGL and with the VTK 
toolkit. This is also shown in the sketch below: 
 

            CSMP (2D)           CSMP, OpenGl, VTK         GoCaD
 
           sc2 ^                    y |                     z |
   2D:   y     |                3D:   |                 3D:   | 
               |                      |                       |
           sc1 +------->              o------->               o------->
               sc1     sc2           /       x               /       y
                    x            z  /                    x  / 
 
Multiple TVOLUME entitities and associated data can be read from the GoCad 
input file. 

@section arguments Input Arguments 

The ASCII text-file stream from which the gocad model shall be read, a
VSet instance into which the data shall be stored, and a flag which
indicates whether node properties shall be read as well, are the input
arguments of ReadTSolid(). 

A fourth argument with the default value 'false' gives the user the 
possibility to closely examine the file reading process. This may be
useful since not all GoCad output files are flawless. 

@param vset the information contained in the GoCaD TSolid dataset extracted from the 
file is stored into the VSet instance supplied as second
argument. 

@section implementation Implementation

The method reads TVOLUME by TVOLUME from the input file and stores these
in STL maps. The neighbors of each Element are identified via string
keys which are constructed for each tetrahedron face.   

@section application Application

The methods provides an interface to the GoCaD geometric modeling tool
such that tetrahedral meshes created by this program can be used in CSMP
computations. 

@section messages Messages 

The method reports the file-reading progress and will print error messages
if elements on the sides of the model bounding box cannot be assigned 
boundary flags by "IndentifyTetrahedronBoundary()". This may be the case,
if the element boundary faces are appreciably inclined relative to the 
normals of the coordinate system (>45o). 
 */
template<uint32_t dim>                        
void GoCadInterface<dim>::ReadTSolid( ifstream& ifs, VSet<dim>& vset, 
                                             bool with_node_prop )
 {
    bool verbose(false);

    if ( !ifs )
      {
         cout <<"\nGoCadInterface::ReadTSolid: Invalid file input stream 'ifs'.";
         cout <<" Reading could not be performed."<< endl;
         return;
      }

    char text_line[256];
    char *token;
    const char* const delims = " ,=,\n,\r,\t"; 
   
    // 1. reading first line and checking whether the right data will follow
    // ---------------------------------------------------------------------
    streampos pos = ifs.tellg();
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims ); 


// -------------------------------------------------------------------------
//
//    FIRST TVOLUME OBJECT
//
//  The tetrahedra for each TVolume object are stored as groups with unique
//  id numbers. Statistics of these are printed after the reading process
//  is finished.
    map<size_t,size_t>  tvolume_map;
    size_t                  tvolume_key = 1;
//
// ------------------------------------------------------------------------- 

   // seeking for word in text string
    if ( strstr( text_line, "TVOLUME" ) == NULL )
      {
         ifs.seekg( pos );
         throw csmp::Exception( ERROR, "GoCadInterface::ReadTSolid", 
                                "TVOLUME descriptor missing.");
         return;
      }
      
   // 2. Reading node data VRTX (here nodes are still consecutively numbered)
    // -------------------------
    //   id      attr   node xyz
    map<size_t,pair<int32_t,vector<double> > >                  nodes;
    typename map<size_t,pair<int32_t,vector<double> > >::const_iterator  nit;
    map<size_t,double>                                     node_property;
    typename map<size_t,double>::const_iterator                     prit;
    vector<double>                 xyz(3);
    size_t                         nodeID(UNSPECIFIED);
    int32_t                        attr;
    pair<int32_t,vector<double> >  node;

    // reading line by line 
    while ( ifs.getline( text_line, 256 ) )
      {
         // tokenizing text line
         token = strtok( text_line, delims ); 
         if ( token != NULL && (!strcmp( token, "TETRA" ) || 
                                !strcmp( token, "PROPERTY_CN" )) ) break;

         // getting line
         if ( token != NULL && (!strcmp( token, "VRTX" ) || !strcmp( token, "PVRTX" )) ) 
           {
             // getting point coordinates
             token = strtok( NULL, delims );
             nodeID = static_cast<uint32_t>(atol( token ));   // node ID
             token = strtok( NULL, delims );
             xyz[0] = atof( token );   // x
             token = strtok( NULL, delims ); 
             xyz[1] = atof( token );   // y         
             token = strtok( NULL, delims );
             xyz[2] = atof( token );   // z
             // node properties
             if ( with_node_prop )
               {
                  token = strtok( NULL, delims );
                  node_property[ nodeID ] = atof( token );   // property value
              }
             // constraint information
             token = strtok( NULL, delims ); 
             if ( token != NULL && !strcmp( token, "CNXYZ" ) )  attr = 1;  // boundary node
             else                                               attr = 0;
             node.first  = attr;
             node.second = xyz;
             
             // assigning node to node list (1...n)
             nodes[ nodeID ] = node;
           } 
      }
    assert( nodeID > 0 && nodeID <= nodes.size() );


   // 4. Ignoring constraint point information
   // ----------------------------------------
   if ( token != NULL && !strcmp( token, "PROPERTY_CN" ) ) 
     while ( ifs.getline( text_line, 256 ) )
       {
          // tokenizing first word of text line
          token = strtok( text_line, delims ); 
          if ( token != NULL && ( !strcmp( token,"END") || !strcmp( token,"END_PROPERTY_CN") ) )
            {
                // next line is needed for function below
                ifs.getline( text_line, 256 );
                token = strtok( text_line, delims ); 
                break;
            }
          // abnormal termination
          if ( token != NULL && strcmp( token, "CNP" ) )
            {
               cout <<"\nGoCadInterface::ReadTSolid: abnormal string read: ";
               cout << token <<"  Now exciting..."<< endl;
               return;
            }
       }

    // 3. Reading tetrahedral elements TETRA
    // -------------------------------------
    //   id      attr   node xyz
    map<size_t,vector<int64_t> >  tetrahedra;
    vector<int64_t>  node_ids(4);
    size_t          tetraID = 1;

    // reading line by line 
    do
      {
         // tokenizing text line
         if ( tetraID != 1 )
           {
              token = strtok( text_line, delims ); 
              if ( token != NULL && !strcmp( token, "END" ) ) break;
              if ( token != NULL && !strcmp( token, "TVOLUME" ) ) break;
           }
         if ( token != NULL && !strcmp( token, "TETRA" ) ) 
           {
             // getting node ID numbers for each tetrahedron
             token = strtok( NULL, delims );
             node_ids[0] = static_cast<uint32_t>(atol( token ));   // node 1
             token = strtok( NULL, delims ); 
             node_ids[1] = static_cast<uint32_t>(atol( token ));   // node 2        
             token = strtok( NULL, delims );
             node_ids[2] = static_cast<uint32_t>(atol( token ));   // node 3
             token = strtok( NULL, delims );
             node_ids[3] = static_cast<uint32_t>(atol( token ));   // node 4
             
             // assigning nodes to tetrahedra and tetrahedra to tvolume
             tetrahedra[ tetraID ]    = node_ids;
             tvolume_map[ tetraID++ ] = tvolume_key;
           } 
      }
    while ( ifs.getline( text_line, 256 ) );

    cout <<"\nGoCadInterface::ReadTSolid: completed reading TVOLUME: 1"<< endl;


// -------------------------------------------------------------------------
//
//   SUBSEQUENT TVOLUME OBJECTS
//
//   - A TSolid file created from a Gocad model will contain several TVOLUMES.
//   - Each TVOLUME description consists of vertices, constraint-points, and
//     tetrahedral-element groupings.
//   - vertices can either be specified as PVTRX (vertex with properties) or 
//     as PATOM id2 id1. This gives a new id2 to an existing vertex with id1.
//   - PATOMS and PVTRX's are numbered sequentially such that new vertices are
//     not in continuation of those of the first TVOLUME.
//
//   NB: here nodes are no longer consecutively numbered ! - but are 
//     numbered interchangeably with PATOM's.
//
// ------------------------------------------------------------------------- 
   // newid oldid
   map<int64_t ,int64_t>                  shared_vtrx;
   typename map<int64_t ,int64_t>::const_iterator  shit;
   int64_t                                  newID;
   bool                                      cnp_data_read;

   // filling node ID's from first TVOLUME into 'shared_vtrx' map
   for ( nit=nodes.begin(); nit!=nodes.end(); nit++ )
     shared_vtrx[ (*nit).first ] = (*nit).first;

   while( !ifs.eof() && token != NULL && !strcmp( token, "TVOLUME" ) )
     {   
        cnp_data_read = false;

        // 1. incrementing the TVOLUME ID
        // ------------------------------ 
        tvolume_key++;
  
        // 2. Reading PATOM or PVTRX lines & skipping CNP data
        // ---------------------------------------------------
        while ( !cnp_data_read && ifs.getline( text_line, 256 )  )
          {
             // if the line contains a reference to an existing node
             // ----------------------------------------------------
             token = strtok( text_line, delims ); 
             if ( token != NULL && !strcmp( token, "PATOM" ) )
               {
                  // new vertex id
                  token  = strtok( NULL, delims );
                  newID  = static_cast<int64_t>(atol( token ));
                  // collocated old vertex id
                  token  = strtok( NULL, delims );
                  nodeID = static_cast<int64_t>(atol( token ));
                  // if there is property information, it is ignored since it
                  // duplicates already existing node data.
                  
                  // storing ID number checking whether it may reference already existing
                  // node ID numbers in the shared node list
                  if ( (shit=shared_vtrx.find(nodeID)) != shared_vtrx.end() )
                    shared_vtrx[ newID ] = (*shit).second;
                  else 
                    shared_vtrx[ newID ] = nodeID;
               }
             // or if a new node is specified
             // -----------------------------
             else if ( token != NULL && (!strcmp( token, "VRTX" ) || !strcmp( token, "PVRTX" )) ) 
               {
                  // node ID is made consecutive using existing list
                  nodeID = nodes.size() + 1;
               
                  // parsed PVTRX ID is not numbered consecutively !
                  token = strtok( NULL, delims );
                  newID = static_cast<uint32_t>(atol( token ));   
                  // getting point coordinates
                  token = strtok( NULL, delims );
                  xyz[0] = atof( token );   // x
                  token = strtok( NULL, delims ); 
                  xyz[1] = atof( token );   // y         
                  token = strtok( NULL, delims );
                  xyz[2] = atof( token );   // z
                  // node properties
                  if ( with_node_prop )
                    {
                       token = strtok( NULL, delims );
                       node_property[ nodeID ] = atof( token );   // property value
                    }
                  // constraint information
                  token = strtok( NULL, delims ); 
                  if ( token != NULL && !strcmp( token, "CNXYZ" ) )  attr = 1;  // boundary node
                  else                                               attr = 0;
                  node.first  = attr;
                  node.second = xyz;
             
                  // assigning node to node list (1...n)
                  nodes[ nodeID ] = node;
                
                  // storing PVRTX id number, since it is used later in TETRA list
                  if ( (shit=shared_vtrx.find(newID)) != shared_vtrx.end() )
                    shared_vtrx[ newID ] = (*shit).second;
                  else 
                    shared_vtrx[ newID ] = nodeID;
               }
             // any constraint point information is ignored
             // -------------------------------------------
             else if ( token != NULL && !strcmp( token, "PROPERTY_CN" ) )
               { 
                  while ( ifs.getline( text_line, 256 ) )
                    {
                       token = strtok( text_line, delims ); 
                 
                       // if end of CNP list is reached, END and tetra section must follow
                       if ( token != NULL && ( !strcmp( token,"END") || !strcmp( token,"END_PROPERTY_CN") ) )
                         {
                            cnp_data_read = true;
                            break;
                         }
                       // if line does not start with CNP or END an error occurred
                       if ( token != NULL && strcmp( token, "CNP" ) )
                         if ( token != NULL && strcmp( token, "END" ) )
                           {
                              // abnormal termination if END is not encountered after CNP list
                              cout <<"\nGoCadInterface::ReadTSolid: abnormal string read: ";
                              cout << token <<"  Now exciting..."<< endl;
                              return;
                           }
                    }
               }
          } // end Reading PATOM or PVTRX lines & skipping CNP data

        // 3. Reading tetrahedral elements TETRA translating their node id's to
        //    CSMP consistent numbering scheme
        // -------------------------------------
        tetraID = tetrahedra.size() + 1;

        while ( ifs.getline( text_line, 256 ) )
          {
             token = strtok( text_line, delims ); 
             if ( token != NULL && !strcmp( token, "TETRA" ) ) 
               {
                  // getting node ID numbers for each tetrahedron
                  token = strtok( NULL, delims );
                  node_ids[0] = static_cast<uint32_t>(atol( token ));   // node 1
                  token = strtok( NULL, delims ); 
                  node_ids[1] = static_cast<uint32_t>(atol( token ));   // node 2        
                  token = strtok( NULL, delims );
                  node_ids[2] = static_cast<uint32_t>(atol( token ));   // node 3
                  token = strtok( NULL, delims );
                  node_ids[3] = static_cast<uint32_t>(atol( token ));   // node 4
            
                  // translating node IDs to CSP ids
                  // -------------------------------
                  for ( auto n=0; n<4; n++ )
                    {
                       if ( (shit=shared_vtrx.find( node_ids[n] )) != shared_vtrx.end() )
                         node_ids[n] = (*shit).second;
                       else
                         {
                            cout <<"\nGoCadInterface::ReadTSolid: ";
                            cout <<"Error: Orphan PATOM vertex encountered: ";
                            cout << node_ids[n] << endl;
                            return;
                         }
                    }
                  // assigning nodes to tetraedra and tetrahedra to tvolumes
                  tetrahedra[ tetraID ]    = node_ids;
                  tvolume_map[ tetraID++ ] = tvolume_key;
               }
             else if ( token != NULL && !strcmp( token, "END" ) ) break;
             else if ( token != NULL && !strcmp( token, "TVOLUME" ) ) break;
          }
        cout <<"\nGoCadInterface::ReadTSolid: completed reading TVOLUME: "<< tvolume_key << endl;

   } // end of SUBSEQUENT TVOLUMES



    // 4. Testing whether input was read correctly
    // -------------------------------------------
    typename map<size_t,vector<int64_t> >::iterator  plit;
    size_t n;
    
    if ( verbose )
      {
         cout <<"\nNodes ("<< nodes.size() <<") and their coordinates: "<< endl;
         cout <<"Printing first 100..."<< endl;
         for ( n=0, nit=nodes.begin(); nit!=nodes.end(); nit++ )
           {
              if ( n++ >= 100 ) break;
              cout <<"Node: "<< (*nit).first <<":  ";
              for ( size_t i{0U}; i<(*nit).second.second.size(); i++ )
              cout << (*nit).second.second[i] <<",\t";
              cout <<"flag: "<< (*nit).second.first << endl;
           }
    
         cout <<"\nTetrahedra ("<< tetrahedra.size() <<") and their member Nodes: "<< endl;
         cout <<"Printing first 100..."<< endl;
         for ( n=0, plit=tetrahedra.begin(); plit!=tetrahedra.end(); plit++ )
           {
              if ( n++ >= 100 ) break;
              cout <<"Tetrahedron: "<< (*plit).first <<":  ";
              for ( size_t i{0U}; i<(*plit).second.size(); i++ )
                cout << (*plit).second[i] <<", ";
              cout << endl;
           }
      }

    // 4b. Checking for consecutive node numbering
    // ------------------------------------------
    if ( !VerifyConsecutiveNodeNumbering( tetrahedra ) )
      throw csmp::Exception( FATAL_ERROR, "GoCadInterface::ReadTSolid", 
            "Nodes in input mesh are not numbered consecutivly");


    // ----------------------------------------------------------------------------------
    // 5. setup VSet
    // ---------------------------------------------------------------------------------- 
    LinearTetrahedron  tetrahedron;
    vset.Resize( tetrahedron.Nodes(),
                 tetrahedron.Neighbors(),
                 tetrahedron.ElementType(), 
                 nodes.size(), tetrahedra.size() );

    // node coordinates
    // NOTE: Y and Z are flipped, since Z points upward in GoCad as 
    // opposed to CSMP, where Z points to the front.
    // ------------------------------------------------------------
    for ( nit=nodes.begin(); nit!=nodes.end(); nit++ )
      {
         vset.Px( (*nit).first, (*nit).second.second[1] ); // Y becomes X
         vset.Py( (*nit).first, (*nit).second.second[2] ); // Z becomes Y
         vset.Pz( (*nit).first, (*nit).second.second[0] ); // X becomes Z
       }
    nodes.erase( nodes.begin(), nodes.end() );

    // plist
    // -----
    vset.AddPlist( tetrahedra.begin(), tetrahedra.end() );  


   // -----------------------------------------------------
   // 6. FINDING THE NEIGHBOR ELEMENTS (opposite each node)
   // -----------------------------------------------------
   cout <<"\nGoCadInterface::ReadTSolid: building 'pfvert' (neighbor element) dataset..."<< endl;
   // 3.1 Making a list of all the faces of each element:
   //     - the faces are located opposite each node
   //     - the nodes of each face are listed in counterclockwise fashion
   //     - looking from the outside in (thus in the basal xy plane the
   //       numbers 1,2,3 are located and the peak node is node 4)
   //    face_key  e-id
   pair<string,size_t>                face;
   multimap<string,size_t>            face_tree;
   typename multimap<string,size_t>::iterator  face_it; 
   HashKey                               hasher; 
   
   // key generation:
   for( plit=tetrahedra.begin(); plit!=tetrahedra.end(); plit++ )
    {
       // making unique search keys from the node ids of the face nodes
       // (convention: face 1 lies opposite of node 1, face 2 opposite of node 2 etc.)
       // ------
       // face 1
       // ------
       hasher.Key( (*plit).second[1], (*plit).second[2], (*plit).second[3], face.first );
       face.second = (*plit).first;
       face_tree.insert( make_pair(face.first,face.second) );
       // ------
       // face 2
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[3], (*plit).second[2], face.first );
       face.second = (*plit).first;
       face_tree.insert( make_pair(face.first,face.second) );
      // ------
       // face 3
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[1], (*plit).second[3], face.first );
       face.second = (*plit).first;
       face_tree.insert( make_pair(face.first,face.second) );
       // ------
       // face 4
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[2], (*plit).second[1], face.first );
       face.second = (*plit).first;
       face_tree.insert( make_pair(face.first,face.second) );
    }
   if ( verbose )
     {
        cout <<"\nGoCadInterface::ReadTSolid: Printing multimap with face keys and corresponding Element ID's:"<< endl;
        for ( face_it=face_tree.begin(); face_it!=face_tree.end(); face_it++ )
          cout <<"\nKey: "<< (*face_it).first <<", parent element ID: "<< (*face_it).second << endl;
        cout << endl;
     }

   // --------------------------------------------------------------------------------------
   // 7. BUILDING Pfverts() by searching for neighbors in hash key map.
   //    The boundary nodes at corresponding faces are flagged exactly as the 
   //    as the faceverts.
   // --------------------------------------------------------------------------------------
   deque<vector<int64_t> > pfverts( tetrahedra.size() );
   vector<int64_t>         pfvert(  tetrahedron.Neighbors() );
   vector<int8_t>         pbflags( vset.Vertices() );
   vector<double>       nd1(3), nd2(3), nd3(3), nd4(3);
   char                   face_key[30];

   cout <<"\nGoCadInterface::ReadTSolid:Identifying the neighbors of each element..." << endl;

   for( plit=tetrahedra.begin(); plit!=tetrahedra.end(); plit++ )
    {
       // 1. identify the coordinates of the four nodes of the element
       // ------------------------------------------------------------
       // node 1 (put into else statement as only needed there)
       nd1[0] = vset.Px( (*plit).second[0] ); 
       nd1[1] = vset.Py( (*plit).second[0] );
       nd1[2] = vset.Pz( (*plit).second[0] );
       // node 2
       nd2[0] = vset.Px( (*plit).second[1] ); 
       nd2[1] = vset.Py( (*plit).second[1] );
       nd2[2] = vset.Pz( (*plit).second[1] );
       // node 3
       nd3[0] = vset.Px( (*plit).second[2] ); 
       nd3[1] = vset.Py( (*plit).second[2] );
       nd3[2] = vset.Pz( (*plit).second[2] );
       // node 4
       nd4[0] = vset.Px( (*plit).second[3] ); 
       nd4[1] = vset.Py( (*plit).second[3] );
       nd4[2] = vset.Pz( (*plit).second[3] );
  
       if ( verbose )
         {
             cout <<"\nElement "<< (*plit).first <<" node coordinates:"<< endl;
             cout <<"Node 1 (x,y,z): "<< nd1[0] <<" "<< nd1[1] <<" "<< nd1[2] << endl;
             cout <<"Node 2 (x,y,z): "<< nd2[0] <<" "<< nd2[1] <<" "<< nd2[2] << endl;
             cout <<"Node 3 (x,y,z): "<< nd3[0] <<" "<< nd3[1] <<" "<< nd3[2] << endl;
             cout <<"Node 4 (x,y,z): "<< nd4[0] <<" "<< nd4[1] <<" "<< nd4[2] << endl;
         }
       // ------
       // face 1
       // ------
       hasher.Key( (*plit).second[1], (*plit).second[2], (*plit).second[3], face_key );
       // finding face with the same key and checking to which element
       // it belongs
       if ( (face_it=face_tree.find( face_key )) == face_tree.end() )
         {
            cout<<"\nGoCadInterface::ReadTSolid:: Error retrieving face hash-key from multimap for element ";
            cout << (*plit).first <<" Erratic key: "<< face_key <<" Reading process was terminated."<< endl;
            return;
         }
       // pre-empt effects of probable bug in CW Pro 5 STL: multimap should always return iterator
       // to first element with key in the map, here it doesn't
       // -----------------------------------------------------
       // wind back one, but only if previous element has same key
       if ( face_it!=face_tree.begin() ) if ( (*(--face_it)).first != face_key ) face_it++;
         
       // if the element ID associated with the key is different from the
       // current element, the neighbor has already been found
       // ----------------------------------------------------
       if ( verbose ) cout <<"\nTest1 face1: Element ID recovered for key: "<< face_key;
       if ( verbose ) cout <<" from map vs. plist iterator elmt ID: "<< (*face_it).second <<" "<< (*plit).first << endl;
       if ( (*face_it).second != (*plit).first ) pfvert[0] = ((*face_it).second);
       // if the face belongs to the same element one searches for
       // the next occurrence of the key in the multimap
       // ----------------------------------------------
       else
       // advance in the multimap
       // -----------------------
         {
            face_it++;
            if ( (*face_it).first == face_key ) pfvert[0] = ((*face_it).second);
            else {
                if ( verbose ) cout <<"\nElement: "<< (*plit).first <<" Test face 1: failed comparison: list-face key: ";
                if ( verbose ) cout <<(*face_it).first <<" vs. hash key: "<< face_key << endl;

                // If there is no shared occurrence of the face, the face lies at the model
                // boundary. If so, the type of this boundary is identified.
                // ---------------------------------------------------------
                pfvert[0] = IdentifyTetrahedronBoundary( nd2, nd3, nd4 );
                  
                // flagging the boudary nodes as such
                // ----------------------------------
                pbflags[ (*plit).second[1] ] = static_cast<int8_t>(pfvert[0]);
                pbflags[ (*plit).second[2] ] = static_cast<int8_t>(pfvert[0]);
                pbflags[ (*plit).second[3] ] = static_cast<int8_t>(pfvert[0]);
             }
         }
       // ------
       // face 2
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[3], (*plit).second[2], face_key );
       if ( (face_it=face_tree.find( face_key )) == face_tree.end() )
         {
            cout<<"\nGoCadInterface::ReadTSolid: Error retrieving face hash-key from multimap for element ";
            cout << (*plit).first <<" Erratic key: "<< face_key <<" Reading process was terminated."<< endl;
            return;
         }
       if ( face_it!=face_tree.begin() ) if ( (*(--face_it)).first != face_key ) face_it++;
       if ( verbose ) cout <<"\nTest1 face2: Element ID recovered for key: "<< face_key;
       if ( verbose ) cout <<" from map vs. plist iterator elmt ID: "<< (*face_it).second <<" "<< (*plit).first << endl;
       if ( (*face_it).second != (*plit).first ) pfvert[1] = ((*face_it).second);
       else
         {
            face_it++;
            if ( (*face_it).first == face_key ) pfvert[1] = ((*face_it).second);
            else {
                if ( verbose ) cout <<"\nElement: "<< (*plit).first <<" Test face 2: failed comparison: list-face key: ";
                if ( verbose ) cout <<(*face_it).first <<" vs. hash key: "<< face_key << endl;
                pfvert[1] = IdentifyTetrahedronBoundary( nd1, nd4, nd3 );
                pbflags[ (*plit).second[0] ] = static_cast<int8_t>(pfvert[1]);
                pbflags[ (*plit).second[3] ] = static_cast<int8_t>(pfvert[1]);
                pbflags[ (*plit).second[2] ] = static_cast<int8_t>(pfvert[1]);
             }
         }
       // ------
       // face 3
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[1], (*plit).second[3], face_key );
       if ( (face_it=face_tree.find( face_key )) == face_tree.end() )
         {
            cout<<"\nGoCadInterface::ReadTSolid: Error retrieving face hash-key from multimap for element ";
            cout << (*plit).first <<" Erratic key: "<< face_key <<" Reading process was terminated."<< endl;
            return;
         }
       if ( face_it!=face_tree.begin() ) if ( (*(--face_it)).first != face_key ) face_it++;
       if ( verbose ) cout <<"\nTest1 face3: Element ID recovered for key: "<< face_key;
       if ( verbose ) cout <<" from map vs. plist iterator elmt ID: "<< (*face_it).second <<" "<< (*plit).first << endl;
       if ( (*face_it).second != (*plit).first ) pfvert[2] = ((*face_it).second);
       else
         {
            face_it++;
            if ( (*face_it).first == face_key ) pfvert[2] = ((*face_it).second);
            else {
                if ( verbose ) cout <<"\nElement: "<< (*plit).first <<" Test face 3: failed comparison: list-face key: ";
                if ( verbose ) cout <<(*face_it).first <<" vs. hash key: "<< face_key << endl;
                pfvert[2] = IdentifyTetrahedronBoundary( nd1, nd2, nd4 );
                pbflags[ (*plit).second[0] ] = static_cast<int8_t>(pfvert[2]);
                pbflags[ (*plit).second[1] ] = static_cast<int8_t>(pfvert[2]);
                pbflags[ (*plit).second[3] ] = static_cast<int8_t>(pfvert[2]);
             }
         }
       // ------
       // face 4
       // ------
       hasher.Key( (*plit).second[0], (*plit).second[2], (*plit).second[1], face_key );
       if ( (face_it=face_tree.find( face_key )) == face_tree.end() )
         {
            cout<<"\nGoCadInterface::ReadTSolid: Error retrieving face hash-key from multimap for element ";
            cout << (*plit).first <<" Erratic key: "<< face_key <<" Reading process was terminated."<< endl;
            return;
         }
       if ( face_it!=face_tree.begin() ) if ( (*(--face_it)).first != face_key ) face_it++;
       if ( verbose ) cout <<"\nTest1 face4: Element ID recovered for key: "<< face_key;
       if ( verbose ) cout <<" from map vs. plist iterator elmt ID: "<< (*face_it).second <<" "<< (*plit).first << endl;
       if ( (*face_it).second != (*plit).first ) pfvert[3] = ((*face_it).second);
       else
         {
            face_it++;
            if ( (*face_it).first == face_key ) pfvert[3] = ((*face_it).second);
            else {
                if ( verbose ) cout <<"\nElement: "<< (*plit).first <<" Test face 4: failed comparison: list-face key: ";
                if ( verbose ) cout <<(*face_it).first <<" vs. hash key: "<< face_key << endl;
                pfvert[3] = IdentifyTetrahedronBoundary( nd1, nd3, nd2 );
                pbflags[ (*plit).second[0] ] = static_cast<int8_t>(pfvert[3]);
                pbflags[ (*plit).second[2] ] = static_cast<int8_t>(pfvert[3]);
                pbflags[ (*plit).second[1] ] = static_cast<int8_t>(pfvert[3]);
             }
         }
       pfverts[ (*plit).first ] = pfvert;
    }

  if ( verbose )
    {
       cout <<"\nDetermined boundary flags:"<< endl;
       size_t counter(0U);
       for ( auto bflit=pbflags.begin(); bflit!=pbflags.end(); bflit++, counter++ )
         if ( (*bflit) < 0 )
           cout <<"\nNode: "<< counter <<" flagged: "<< (*bflit);
       cout << endl; 
    }
    
  // -------------------------------------------
  // 8. Writing remaining results to the VSet 
  // -------------------------------------------
  cout <<"\nGoCadInterface::ReadTSolid: Writing remaining finite-element data to VSet..." << endl;
  vset.AddPfverts( pfverts.begin(), pfverts.end() );
  vset.AddBFlags( pbflags.begin(), pbflags.end() );

  FlagEdgeNodesOfBoxShapedModel( vset );      

  // -------------------------------------------
  // 9. Adding property information to Vset
  // -------------------------------------------
  typename map<size_t,size_t>::const_iterator  tvit;

  // elemental permeability is used to convey tvolume IDs to CSMP
  assert( tvolume_map.size() == tetrahedra.size() );
  //FEM_Data<ScalarVariable >  elmt_data( ELEMENT, tetrahedra.size() );
  PropertyData  data( ELEMENT, SCALAR, dim );
  data.Reserve( tetrahedra.size() );
 
  cout <<"\nGoCadInterface::ReadTSolid: Assigning TVOLUME numbers to element variable 'permeability'" << endl;
  for ( tvit=tvolume_map.begin(); tvit!=tvolume_map.end(); tvit++ )
    pushBack( data, makeScalar( PLAIN, (*tvit).second ) );

  vset.AddData( "permeability", data );

  if ( verbose ) vset.Out();

 } // end ReadTSolid





/**
 
On the basis of the coordinate extrema, method applies edge flags to the
boundary nodes of orthogonally box-shaped models according to the edge
flags defined in 'CSMP_definitions.h'.

@section messages Messages 

Method will detect when the model is not orthogonal and issue warnings for
the nodes which could not be assigned the desired edge flags.  
*/
template<uint32_t dim>
void GoCadInterface<dim>::FlagEdgeNodesOfBoxShapedModel( VSet<dim>& vset )
 {
    // 1. finding coordinate extrema for the boundary nodes
    double xmin, xmax, ymin, ymax, zmin, zmax, xco, yco, zco;
    vset.CoordinateRange( 0, xmin, xmax );
    vset.CoordinateRange( 1, ymin, ymax );
    vset.CoordinateRange( 2, zmin, zmax );
    
    // 2. checking whether the model is actually 3-dimensional
    if ( xmin == xmax || ymin == ymax || zmin == zmax ) {
         throw csmp::Exception( WARNING, "GoCadInterface::FlagEdgeNodesOfBoxShapedModel",
                         "Model appears to be 2-dimensional. No edge-node flagging needed");
         return;
      }
      
    // 3. Flagging the edge nodes according to their coordinates
    size_t counter(0U);
    for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, counter++ )
      // ignoring nodes which were already identified as model corners 
      if ( (*bit) != CNR_MIN      && (*bit) != CNR_MAX &&
           (*bit) != CNR_MIN_MAXX && (*bit) != CNR_MIN_MAXXZ &&
           (*bit) != CNR_MIN_MAXZ && (*bit) != CNR_MAX_MINXZ &&
           (*bit) != CNR_MAX_MAXX && (*bit) != CNR_MAX_MAXZ )
        {
           xco = vset.Px( counter );
           yco = vset.Py( counter );
           zco = vset.Pz( counter );
        
           // 3.1 back wall (z=zmin)
           if ( zco == zmin )
             {
                if ( yco == ymin && (xmin < xco && xco < xmax) ) (*bit) = BACK_BOTTOM;
                if ( yco == ymax && (xmin < xco && xco < xmax) ) (*bit) = BACK_TOP;
                if ( xco == xmin && (ymin < yco && yco < ymax) ) (*bit) = BACK_LEFT;
                if ( xco == xmax && (ymin < yco && yco < ymax) ) (*bit) = BACK_RIGHT;
             }
           // 3.2 front wall (z=zmax)
           else if ( zco == zmax )
             {
                if ( yco == ymin && (xmin < xco && xco < xmax) ) (*bit) = FRONT_BOTTOM;
                if ( yco == ymax && (xmin < xco && xco < xmax) ) (*bit) = FRONT_TOP;
                if ( xco == xmin && (ymin < yco && yco < ymax) ) (*bit) = FRONT_LEFT;
                if ( xco == xmax && (ymin < yco && yco < ymax) ) (*bit) = FRONT_RIGHT;
             }
           // 3.3 side edges ( zmin < z < zmax )
           else if ( zco > zmin && zco < zmax )
             {
                if ( xco == xmin && yco == ymin ) (*bit) = BOTTOM_LEFT;
                if ( xco == xmax && yco == ymin ) (*bit) = BOTTOM_RIGHT;
                if ( xco == xmax && yco == ymax ) (*bit) = TOP_RIGHT;
                if ( xco == xmin && yco == ymax ) (*bit) = TOP_LEFT;
             }
      }
 
 } // end FlagEdgeNodesOfBoxShapedModel




/**
 
Outputs a physical variable from the Model to a GoCad TSurf object
which is stored into a textfile with the recommended extension '.ts' for 
direct import into GoCad as GoCad object. 

The GoCad TSurf format is defined in Appendix B of the GoCad manual
(p. B-10). In brief, it is a list of the following items: 
 
a GoCad header,
PVTRX information: vertex coordinates and associated property values, 
TRGL triangle connectivity information: the numbers of nodes in each triangle,
END the end of file marker.
 

@section arguments Input Arguments 

OutputVariableToTSurface() takes three arguments, the name of the output 
file to which the third argument which is the timestep and the extension
'.ts' are appended. The second argument specifies the variable which is 
output to file. Thus far, if variables  
are element variables they will still be stored on the nodes in the 
output file. This may cause boundary nodes of regions appear oddly valued,
since the element node values will overlap with eachother depending on the
node numbering in the model. 

@section implementation Implementation

OutputVariableToTSurface() uses several atomistic class objects to 
represent gocad formats, these are the GocadHeader and the 
GocadPropertyClassHeader objects. 

The method only works for two-dimensional models. 3D models must
be output as GoCad TSolid objects (see corresponding documentation). 

IntegrationPoint variables cannot be output to TSurf files. 

@section application Application

OutputVariableToTSurface() provides a direct interface to the GoCad 
visualization and modeling program. Use it to visualize scalar or vector
variables in GoCad. 

@section messages Messages 

If the output file cannot be opened, if the model is not two-dimensional,
or if an attempt is made to output a IntegrationPoint variable, the method
will report an error and exit without completing its task. 
*/
template<uint32_t dim>
void  GoCadInterface<dim>::OutputVariableToTSurface( const Model<dim>& sgroup, 
                                                        const char* fname, 
                                                        const char* var, 
                                                        long        timestep )
const
 {
   const Region<dim>&  sg(sgroup.Region("Model"));
   sg.UpdateMemberIndexes();
 
   csmp::Index prop_key = sgroup.Database().StorageKey(var);
   ofstream  ofs;
   char      file[200], step[100];
   strcpy( file, fname ); 
   if ( timestep > 0 )
     {
        snprintf( step, sizeof(step), "%lu", timestep );
        strcat( file, step );
     }
   strcat( file, ".ts" );
   
   ofs.open ( file, ios::out|ios::trunc );
   if ( !ofs )
     {
         throw csmp::Exception( FATAL_ERROR, "Model::OutputVariableToTSurface", 
                                      "Output file could not be opened");
         return;
      }
   if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) 
     {
         throw csmp::Exception( ERROR, "Model::OutputVariableToTSurface", 
                                "IntegrationPoint variables must be output to Point format");
         return;
     } 

   // 4. writing the file header
   // --------------------------
   double  pmin, pmax;
   char object_name[200];
   strcpy( object_name, "CSP_object" );
   if ( timestep > 0 ) strcat( object_name, step );
   GocadHeader  header( "TSurf", object_name, fname );
   sg.MinMaxOf( var, pmin, pmax ); 
   if ( prop_key.type == SCALAR )      header.AddProperty( var, pmin, pmax );
   else if ( prop_key.type == VECTOR ) header.AddProperty( var, pmin, pmax, 99, 3 );
   else if ( prop_key.type == TENSOR ) header.AddProperty( var, pmin, pmax, 99, 9 );
   // so far this only handles scalars in the e-size section
   header.WriteToText( ofs );

   // 4.2 Writing the vertex=Node data to file
   // ----------------------------------------
   ScalarVariable                       sc;
   VectorVariable<dim>                   vc;
   TensorVariable<dim>                   ts;

   // 5. All properties are output as NODE Properties
   // ------------------------------------------------
   // Thus, if element properties are to be output one
   // must first create an intermediate array to store these 
   // and then output from this array.
   vector<ScalarVariable >      sc_elmt_data;
   vector<VectorVariable<dim> > vc_elmt_data;
   vector<TensorVariable<dim> > ts_elmt_data;
   vector<size_t>               node_ids;
   
   if ( prop_key.place == ELEMENT ) 
     {
        if       ( prop_key.type == SCALAR ) { 
             sc_elmt_data.reserve( sg.Nodes() );
             for ( size_t i{0U}; i<sg.Nodes(); i++ ) sc_elmt_data.push_back( ScalarVariable() );
          }
        else if  ( prop_key.type == VECTOR ) {
             vc_elmt_data.reserve( sg.Nodes() );
             for ( size_t i{0U}; i<sg.Nodes(); i++ ) vc_elmt_data.push_back( VectorVariable<dim>() );
          }
        else if  ( prop_key.type == TENSOR ) {
             ts_elmt_data.reserve( sg.Nodes() );
             for ( size_t i{0U}; i<sg.Nodes(); i++ ) ts_elmt_data.push_back( TensorVariable<dim>() );
          }
        for ( auto eit=sg.CellsBegin(); eit!=sg.CellsEnd(); eit++ )
          {
              switch ( prop_key.type )
                {
                    case SCALAR:
                         (*eit)->Read(prop_key, sc );
                         // triangles only
                         for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) sc_elmt_data[ (*eit)->N(i)->Idx() ] = sc;
                      break;
                    case VECTOR:
                         (*eit)->Read(prop_key, vc );
                         for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) vc_elmt_data[ (*eit)->N(i)->Idx() ] = vc;
                      break;
                    case TENSOR:
                         (*eit)->Read(prop_key, ts );
                         for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) ts_elmt_data[ (*eit)->N(i)->Idx() ] = ts;
                    default:
                      throw out_of_range("GoCadInterface<dim>: variable type not handled");
                }
          }
     }
   
   // outputting
   // ----------
   for ( auto nit=sg.NodesBegin(); nit != sg.NodesEnd(); ++nit )
     {
        ofs <<"PVRTX "<< (*nit)->Idx() <<" ";
        ofs << (*nit)->x() <<" "<< (*nit)->y() <<" "<< 0.0 <<" ";

        // 5.1 if Node properties are output
        // ---------------------------------
        if ( prop_key.place == NODE ) switch( prop_key.type )
           {  
              case SCALAR:
                   (*nit)->Read(prop_key, sc ); 
                   ofs << sc() <<" ";
                break;
              case VECTOR:
                   (*nit)->Read(prop_key, vc );
                   // z-component was assigned 0.0 in 2D case 
                   for( auto i{0U}; i<3; i++ ) ofs << vc[i] <<" ";
                break;
              case TENSOR:
                   (*nit)->Read(prop_key, ts );
                   for( auto i{0U}; i<3; i++ )
                      for( auto j{0U}; j<3; j++ ) ofs << ts(i,j) <<" ";
                    default:
                      throw out_of_range("GoCadInterface<dim>: variable type not handled");
           }
        // 5.2 Element properties
        // ----------------------
        else if ( prop_key.place == ELEMENT ) switch( prop_key.type )
           {  
              case SCALAR:
                   ofs << sc_elmt_data[ (*nit)->Idx() ]() <<" ";
                break;
              case VECTOR:
                   // z-component was assigned 0.0 in 2D case 
                   for ( auto i{0U}; i<3; i++ )
                      ofs << vc_elmt_data[ (*nit)->Idx() ][i] <<" ";
                break;
              case TENSOR:
                   for ( auto i{0U}; i<3; i++ )
                      for ( auto j{0U}; j<3; j++ )
                          ofs << ts_elmt_data[ (*nit)->Idx() ](i,j) <<" ";
                    default:
                      throw out_of_range("GoCadInterface<dim>: variable type not handled");
           }
        if ( (*nit)->AtBoundary() != NOT ) ofs <<" CNXYZ";
        ofs << endl;
     }

    // 4.3 Writing plist (nodes that make up the tetrahedra)
    // -----------------------------------------------------
    for ( auto eit = sg.CellsBegin(); eit != sg.CellsEnd(); eit++ )
      {
         ofs <<"TRGL ";
         for ( auto j{0U}; j<(*eit)->Nodes(); j++ ) ofs << (*eit)->N(j)->Idx() <<" ";
         ofs << endl;
      }
    ofs << "END"<< endl;    

    ofs.close();
    cout <<"\nModel::OutputVariableToTSurface: file '"<< file <<"' written successfully."<< endl;

   } // end OutputVariableToTSurface

  




/**
 
Outputs a physical variable in Region subregions of the Model to a 
GoCad TSurf object which is stored into a textfile with the recommended 
extension '.ts' for direct import into GoCad as GoCad object. 

The GoCad TSurf format is defined in Appendix B of the GoCad manual
(p. B-10). In brief, it is a list of the following items: 

 
a GoCad header,
PVTRX information: vertex coordinates and associated property values, 
TRGL triangle connectivity information: the numbers of nodes in each triangle,
END the end of file marker.
 

@section arguments Input Arguments 

OutputVariableToTSurface() takes four arguments, the name of the output 
file to which the fourth argument which is the timestep and the extension
'.ts' are appended. The second argument specifies the name of the Region subregion
of the model and the third argument indicates the variable which shall be
output to file. Thus far, if variables  
are element variables they will still be stored on the nodes in the 
output file. This may cause boundary nodes of regions appear oddly valued,
since the element node values will overlap with eachother depending on the
node numbering in the model. 

@section implementation Implementation

OutputVariableToTSurface() uses several atomistic class objects to 
represent gocad formats, these are the GocadHeader and the 
GocadPropertyClassHeader objects. 

The method only works for two-dimensional models. 3D models must
be output as GoCad TSolid objects (see corresponding documentation). 

IntegrationPoint variables cannot be output to TSurf files. 

@section application Application

OutputVariableToTSurface() provides a direct interface to the GoCad 
visualization and modeling program. Use it to visualize scalar or vector
variables in GoCad. 

@section messages Messages 

If the output file cannot be opened, if the model is not two-dimensional,
if the target Region is not defined,
or if an attempt is made to output a IntegrationPoint variable, the method
will report an error and exit without completing its task. 
*/
template<uint32_t dim>
void  GoCadInterface<dim>::OutputVariableToTSurface( const Model<dim>& sgroup, 
                                                    const char* fname,
                                                    const char* group_name, 
                                                    const char* var, 
                                                    long   timestep ) 
const
 {
   const Region<dim>&  sg(sgroup.Region("Model"));
   sg.UpdateMemberIndexes();

   // 1. Property Access in the Model
   // ------------------------------------
   csmp::Index           prop_key = sgroup.Database().StorageKey(var);
   const Region<dim>&  gref     = sgroup.Region( group_name );

   if ( prop_key.type != SCALAR ) {
        throw csmp::Exception( ERROR, "Model<dim>::OutputVariableToTSurface", 
                       "With this method only scalars can be handled. Nothing is done..." ); 		       
        return;
     }

   // 3. Getting output file ready
   // ----------------------------
   ofstream  ofs;
   char      file[200], step[100];
   strcpy( file, fname ); 
   if ( timestep > 0 )
     {
        snprintf( step, sizeof(step), "%ld", timestep );
        strcat( file, step );
     }
   strcat( file, ".ts" );
   
   ofs.open ( file, ios::out|ios::trunc );
   if ( !ofs )
     {
         throw csmp::Exception( FATAL_ERROR, "Model<dim>::OutputVariableToTSurface", 
                                      "Output file could not be opened");
         return;
      }

   // 4. Writing data to file
   // -----------------------
   // 4.1 writing the file header
   // ---------------------------
   double  pmin, pmax;
   gref.MinMaxOf( var, pmin, pmax );
   char object_name[200];
   strcpy( object_name, group_name );
   if ( timestep > 0 ) strcat( object_name, step );
   GocadHeader  header( "TSurf", object_name, "rock_surface" );
   header.AddProperty( var, pmin, pmax );
   header.WriteToText( ofs );

   // 4.2 Writing the vertex=Node data to file
   // ----------------------------------------
   ScalarVariable  val;

   // 4.3 making a list of nodes that belong to the group
   // ---------------------------------------------------
   vector<size_t>  nodes(   gref.Nodes() );
   vector<size_t>  cpoints( gref.IntegrationPoints() );
   
   size_t n{0};
   for ( auto it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ ) nodes[n++] = (*it)->Idx();

   // 4.4 making a list of element properties if these are required
   // -------------------------------------------------------------
   //  All properties are output as NODE Properties.
   // Thus, if element properties are to be output one
   // must first create an intermediate array to store these 
   // and then output from this array.
   ScalarVariable              sc;
   VectorVariable<dim>          vc;
   TensorVariable<dim>          ts;
   Point<dim>                   xyz;
   vector<ScalarVariable >     sc_elmt_data;
   vector<VectorVariable<dim> > vc_elmt_data;
   vector<TensorVariable<dim> > ts_elmt_data;
   vector<size_t>               node_ids;
   // node-id's of elements
   set<size_t>                          elmts;
   typename set<size_t>::const_iterator en_it;
   
   if ( prop_key.place == ELEMENT ) 
     {
        if       ( prop_key.type == SCALAR ) { 
             sc_elmt_data.reserve( gref.Nodes() );
             for ( size_t i{0U}; i<gref.Nodes(); i++ ) sc_elmt_data.push_back( ScalarVariable() );
          }
        else if  ( prop_key.type == VECTOR ) {
             vc_elmt_data.reserve( gref.Nodes() );
             for ( size_t i{0U}; i<gref.Nodes(); i++ ) vc_elmt_data.push_back( VectorVariable<dim>() );
          }
        else if  ( prop_key.type == TENSOR ) {
             ts_elmt_data.reserve( gref.Nodes() );
             for ( size_t i{0U}; i<gref.Nodes(); i++ ) ts_elmt_data.push_back( TensorVariable<dim>() );
          }
        for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); it++ )
          {
              for ( auto j{0U}; j<(*it)->Nodes(); j++ ) elmts.insert( (*it)->N(j)->Idx() );
              switch ( prop_key.type )
                {
                    case SCALAR:
                         (*it)->Read(prop_key, sc );
                         // triangles only
                         for ( auto i{0U}; i<(*it)->Nodes(); i++ ) sc_elmt_data[ (*it)->N(i)->Idx() ] = sc;
                      break;
                    case VECTOR:
                         (*it)->Read(prop_key, vc );
                         for ( auto i{0U}; i<(*it)->Nodes(); i++ ) vc_elmt_data[ (*it)->N(i)->Idx() ] = vc;
                      break;
                    case TENSOR:
                         (*it)->Read(prop_key, ts );
                         for ( auto i{0U}; i<(*it)->Nodes(); i++ ) ts_elmt_data[ (*it)->N(i)->Idx() ] = ts;
                    default:
                      throw Exception( ERROR, "GoCadInterface<dim>::OutputVariableToTSurface",
                                      "Array or FlaggedArray variables are not handled yet.");
                }
          }
     }

   // 5. Outputting data to GoCad text file
   // -------------------------------------
   switch( prop_key.place )
     {
        case NODE:
         for ( size_t i{0U}; i<nodes.size(); i++ )
               {
                  gref.N( nodes[i] )->Read(prop_key, val );
                  ofs <<"PVRTX "<< i+1 <<" ";
                  ofs << gref.N( nodes[i] )->x() <<" "<< gref.N( nodes[i] )->y() <<" "<< 0.0 <<" ";
                  ofs << val();
                  if ( gref.N( nodes[i] )->AtBoundary() != NOT ) ofs <<" CNXYZ";
                  ofs << endl;
               } 
          break;
        case ELEMENT_INTEGRATION_POINT:
             for ( auto it= gref.CellsBegin(); it!=gref.CellsEnd(); it++ )
               for ( auto i{0}; i<(*it)->IntegrationPoints(); i++ )
               {
                  // getting constraint point coordinates
                  Point<dim>  pt((*it)->IntegrationPoint(i));
                                                                    
                  (*it)->Read( i, prop_key, val );
                  ofs <<"PVRTX "<< i+1 <<" ";
                  if ( dim == 3U ) ofs << pt[0] <<" "<< pt[1] <<" "<< pt[2] <<" ";
                  else ofs << pt[0] <<" "<< pt[1] <<" "<< 0. <<" ";
                  ofs << val();
                  if ( atBoundary(*it) != NOT ) ofs <<" CNXYZ";
                  ofs << endl;
               } 
           break;
       case ELEMENT: {
             size_t counter{0};
             for ( en_it=elmts.begin(); en_it!=elmts.end(); en_it++ )
               {
                  ofs <<"PVRTX "<< counter++ <<" ";
                  ofs << gref.N( (*en_it) )->x() <<" "<< gref.N( (*en_it) )->y() <<" "<< 0.0 <<" ";
                  switch( prop_key.type )
                    {  
                       case SCALAR:
                            ofs << sc_elmt_data[ (*en_it) ]() <<" ";
                         break;
                       case VECTOR:
                            for ( auto j{0U}; j<3; j++ )
                               ofs << vc_elmt_data[ (*en_it) ][j] <<" ";
                         break;
                       case TENSOR:
                            for ( auto k=0; k<3; k++ )
                              for ( auto j{0U}; j<3; j++ )
                                 ofs << ts_elmt_data[ (*en_it) ](k,j) <<" ";
                      default:
                         throw out_of_range("GoCadInterface<dim>: variable type not handled");
                     }
                  if ( gref.N( (*en_it) )->AtBoundary() != NOT ) ofs <<" CNXYZ";
                  ofs << endl;
               }
            }
          break;
        default:
          cout <<"\nModel<dim>::OutputVariableToTSurface: ";
          cout <<"Property placement could not be parsed. Nothing was done."<< endl;
          return;
     }

    // 4.3 Writing plist (nodes that make up the tetrahedra)
    // -----------------------------------------------------
    map<size_t,size_t> new_node_ids;
    for ( size_t i{0U}; i<nodes.size(); i++ ) new_node_ids[ nodes[i] ] = i+1;

    for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); it++ )
      {
         ofs <<"TRGL ";
         for ( auto j{0U}; j<(*it)->Nodes(); j++ )
           ofs << (*new_node_ids.find( (*it)->N(j)->Idx() )).second <<" ";
         ofs << endl;
      }
    ofs << "END"<< endl;    

    ofs.close();
    cout <<"\nModel<dim>::OutputVariableToTSurface: file '"<< file <<"' written successfully."<< endl;

 } // end OutputVariableToTSurface









/**
 
Uses the Model method OutputVariableToTSurface() to output all the 
variables of a current model to a GoCad TSurf file with the extension
'.ts'.

@section arguments Input Arguments 

The name of the output file and the model timestep for which the variables
shall be reported. The timestep will be appended to the filename before 
the extension '.ts'. 

@section implementation Implementation

Uses the Model method OutputVariableToTSurface() to output individual
variables. Please to refer to its documentation for further implementation
details. 

@section application Application

Externalize the property values of a run to visualize them in GoCad. 

@section messages Messages 

The method will report errors if the output file cannot be opened, if no 
physical variables are defined or if the model is not two-dimensional. 
*/
template<uint32_t dim>
void  GoCadInterface<dim>::OutputVariablesToTSurface( const Model<dim>& sgroup, 
                                                      const char* fname, long timestep ) const
 {
   const Region<dim>&  sg(sgroup.Region("Model"));
   sg.UpdateMemberIndexes();

   // 1. Property Access in the Model
   // ------------------------------------
   map<string,Index>  prop_list;
   sgroup.Database().ListVariables( prop_list );

   if ( prop_list.empty() )
     {
        throw csmp::Exception( WARNING, "Model<dim>::OutputVariablesToTSurface", 
                                 "No properties were found. Nothing is done..." ); 		       
        return;
     }
   if ( dim != 2 )
     {
        throw csmp::Exception( WARNING, "Model<dim>::OutputVariablesToTSurface", 
                                 "This method is only for surface (2D) data. Nothing is done..." ); 		       
        return;
     }

   // 3. Getting output file ready
   // ----------------------------
   ofstream  ofs;
   char      file[200], step[100];
   strcpy( file, fname ); 
   if ( timestep > 0 )
     {
        snprintf( step, sizeof(step), "%ld", timestep );
        strcat( file, step );
     }
   strcat( file, ".ts" );
   
   ofs.open ( file, ios::out|ios::trunc );
   if ( !ofs )
     {
         throw csmp::Exception( FATAL_ERROR, "Model<dim>::OutputVariablesToTSurface", 
                                      "Output file could not be opened");
         return;
      }

   // 4. Writing data to file
   // -----------------------
   // 4.1 writing the file header
   // ---------------------------
   double  pmin, pmax;
   char object_name[200];
   strcpy( object_name, "CSP_object" );
   strcat( object_name, step );
   GocadHeader  header( "TSurf", object_name, fname );

   for ( auto pit=prop_list.begin(); pit!=prop_list.end(); pit++ )
     {
        if ( (*pit).second.place == NODE )
          {
             sg.MinMaxOf( (*pit).first.c_str(), pmin, pmax ); 
             if      ( (*pit).second.type == SCALAR ) 
                header.AddProperty( (*pit).first.c_str(), pmin, pmax );
             else if ( (*pit).second.type == VECTOR ) 
                header.AddProperty( (*pit).first.c_str(), pmin, pmax, 99, 3 );
             else if ( (*pit).second.type == TENSOR ) 
                header.AddProperty( (*pit).first.c_str(), pmin, pmax, 99, 9 );
          }
     }
   // so far this only handles scalars in the e-size section
   header.WriteToText( ofs );

   // 4.2 Writing the vertex=Node data to file
   // ----------------------------------------
   ScalarVariable                                val;
   VectorVariable<dim>                           vc;
   TensorVariable<dim>                           ts;
   auto nit = sg.NodesBegin();

   // All properties are output as NODE Properties
   // --------------------------------------------
   while ( nit != sg.NodesEnd() )
     {
        ofs <<"PVRTX "<< (*nit)->Idx() <<" ";
        ofs << (*nit)->x() <<" "<< (*nit)->y() <<" "<< 0.0 <<" ";
                
        for ( auto pit=prop_list.begin(); pit!=prop_list.end(); pit++ )
          {
             csmp::Index  prop_key((*pit).second);
             if ( prop_key.place == NODE )
               {
                  switch( prop_key.type )
                    {  
                       case SCALAR:
                            (*nit)->Read(prop_key, val ); 
                            ofs << val() <<" ";
                         break;
                       case VECTOR:
                            (*nit)->Read(prop_key, vc );
                            // z-component was assigned 0.0 in 2D case 
                            for( auto i{0U}; i<3; i++ ) ofs << vc[i] <<" ";
                         break;
                       case TENSOR:
                           for( auto i{0U}; i<3; i++ )
                             for( auto j{0U}; j<3; j++ ) ofs << ts(i,j) <<" ";
                       default:
                         throw out_of_range("GoCadInterface<dim>::OutputVariablesToTSurface: variable type not handled");
                     }
               }
          }
        if ( (*nit)->AtBoundary() != NOT ) ofs <<" CNXYZ";
        ofs << endl;
        nit++;
     } 

    // 4.3 Writing plist (nodes that make up the tetrahedra)
    // -----------------------------------------------------
    for ( auto eit=sg.CellsBegin(); eit != sg.CellsEnd(); ++eit )
      {
         ofs <<"TRGL ";
         for ( auto j{0U}; j<(*eit)->Nodes(); j++ ) ofs << (*eit)->N(j)->Idx() <<" ";
         ofs << endl;
      }
    ofs << "END"<< endl;    

    ofs.close();
    cout <<"\nGoCadInterface<dim>::OutputVariablesToTSurface: file '"<< file <<"' written successfully."<< endl;

 } // end OutputVariablesToTSurface








/**
 
Given a three-dimensional CSMP model, OutputVariableToTSolid() creates a
GoCad TSolid dataset for the indicated variable and output it to a textfile
with the GoCad specific extension '.so'. The TSolid format is described in
the GoCad manual, Appendix B (p. B-14) and consists of: 

 
GoCad header,
PVTRX list: node coordinate and property data listed by node number,
TETRA list: nodes per tetrahedral finite element,
END end of file specifier
  

Importantly, a TSolid file written by GoCad will contain constraint point
lists preceding each TETRA block. In this fashion, the model volume is 
described object by object.  

@section arguments Input Arguments 

OutputVariableToTSolid() expects three arguments, the third of which is
the timestep in a simulation which is represented by the dataset. This
argument is optional. The first argument is the name of the file to which
the TSolid dataset is written to and to whose name the timestep is 
appended to. If the timestep is not specified, a zero will be appended to
the filename before the extension '.so'. The second argument defines the 
variable which shall be output to file.  

@section implementation Implementation

The method uses the GocadHeader object to write the file and the 
property class header for the TSolid file. 

If the method is used to output element properties, these will be mapped
onto the nodes, since GoCad only permits nodal variables. To achieve 
the mapping, intermediate variable vectors are created and all nodes
of each element are assigned the corresponding element property value. 
This strategy has the disadvantage that nodes with a higher ID number 
will overwrite the property values on nodes with lower numbers. The effect
may be a rugged appearance of element property boundaries. User
suggestions are welcome to overcome this deficiency. 

@section application Application

The method is used to output a property to a TSolid for an entire- as 
opposed to a region of a CSMP model. 

@section messages Messages 

The method will report an error and fail to write an output file, if it
cannot open the specified file with extension '.so' in text mode. 
*/
template<uint32_t dim>
void GoCadInterface<dim>::OutputVariableToTSolid( const Model<dim>& sgroup, 
                                                     const char* fname, 
                                                     const char* var, 
                                                     long timestep ) const
 {
   const Region<dim>&  sg(sgroup.Region("Model"));
   sg.UpdateMemberIndexes();

   // 1. Property Access in the Model
   // ------------------------------------
   csmp::Index prop_key = sgroup.Database().StorageKey(var);

   // 2. Getting output file ready
   // ----------------------------
   ofstream  ofs;
   char      file[200], step[100];
   strcpy( file, fname ); 
   if ( timestep > 0 )
     {
        snprintf( step, sizeof(step), "%ld", timestep );
        strcat( file, step );
     }
   strcat( file, ".so" );
   
   ofs.open ( file, ios::out|ios::trunc );
   if ( !ofs )
     {
         throw csmp::Exception( ERROR, "Model<dim>::OutputVariableToTSolid", 
                                "Output file could not be opened");
         return;
      }

   // 3. Writing data to file
   // -----------------------
   // 3.1 writing the file header
   // ---------------------------
   double  pmin, pmax;
   sg.MinMaxOf( var, pmin, pmax );
   char object_name[200];
   strcpy( object_name, "CSMP_object" );
   if ( timestep > 0 ) strcat( object_name, step );
   GocadHeader  header( "TSolid", object_name, "rock_volume" );
   if      ( prop_key.type == SCALAR ) header.AddProperty( var, pmin, pmax );
   else if ( prop_key.type == VECTOR ) header.AddProperty( var, pmin, pmax, 99, 3 );
   else if ( prop_key.type == TENSOR ) header.AddProperty( var, pmin, pmax, 99, 9 );
   header.WriteToText( ofs );

   // 3.2 Writing the vertex=Node data to file
   // ----------------------------------------
   ScalarVariable                               sc;
   VectorVariable<dim>                           vc;
   TensorVariable<dim>                           ts;
   vector<double>                                xyz;

   // 3.3 All properties are output as NODE Properties
   // ------------------------------------------------
   // Thus, if element properties are to be output one
   // must first create an intermediate array to store these 
   // and then output from this array.
   vector<ScalarVariable > sc_elmt_data;
   vector<VectorVariable<dim> > vc_elmt_data;
   vector<TensorVariable<dim> > ts_elmt_data;
   
   if ( prop_key.place == ELEMENT ) 
     {
        if ( prop_key.type == SCALAR ) { 
             sc_elmt_data.reserve( sg.Nodes() );
             for ( size_t i{0U}; i<sg.Nodes(); i++ ) sc_elmt_data.push_back( ScalarVariable() );
          }
        else if  ( prop_key.type == VECTOR ) {
             vc_elmt_data.reserve( sg.Nodes() );
             for ( size_t i{0U}; i<sg.Nodes(); i++ ) vc_elmt_data.push_back( VectorVariable<dim>() );
          }
        else if  ( prop_key.type == TENSOR ) {
             ts_elmt_data.reserve( sg.Nodes() );
             for ( size_t i{0U}; i<sg.Nodes(); i++ ) ts_elmt_data.push_back( TensorVariable<dim>() );
          }
        for ( auto eit=sg.CellsBegin(); eit!=sg.CellsEnd(); eit++ )
          {
              switch ( prop_key.type )
                {
                    case SCALAR:
                         (*eit)->Read(prop_key, sc );
                         // triangles only
                         for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) sc_elmt_data[ (*eit)->N(i)->Idx() ] = sc;
                      break;
                    case VECTOR:
                         (*eit)->Read(prop_key, vc );
                         for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) vc_elmt_data[ (*eit)->N(i)->Idx() ] = vc;
                      break;
                    case TENSOR:
                         (*eit)->Read(prop_key, ts );
                         for ( auto i{0}; i<(*eit)->Nodes(); i++ ) ts_elmt_data[ (*eit)->N(i)->Idx() ] = ts;
                    default:
                      throw out_of_range("GoCadInterface<dim>: variable type not handled");
                }
          }
     }

   // 4. Writing properties to file
   // -----------------------------
   switch( prop_key.place )
     {
        case NODE:
             for ( auto nit=sg.NodesBegin(); nit != sg.NodesEnd(); ++nit )
               {
                  (*nit)->Read(prop_key, sc );
                  ofs <<"PVRTX "<< (*nit)->Idx() <<" ";
                  ofs << (*nit)->x() <<" "<< (*nit)->y() <<" "<< (*nit)->z() <<" ";
                  ofs << sc();
                  if ( (*nit)->AtBoundary() != NOT ) ofs <<" CNXYZ";
                  ofs << endl;
               }
          break;
        case ELEMENT_INTEGRATION_POINT: {
               size_t  counter(0U);
               for ( auto eit=sg.CellsBegin(); eit != sg.CellsEnd(); ++eit )
                 {
                    for ( auto i{0}; i<(*eit)->IntegrationPoints(); i++ ) {
                        // getting constraint point coordinates
                        Point<dim>  pt((*eit)->IntegrationPoint(i));
                        (*eit)->Read( i, prop_key, sc );
                        ofs <<"PVRTX "<< counter++ <<" ";
                        ofs << pt[0] <<" "<< pt[1] <<" "<< pt[2] <<" ";
                        ofs << sc();
                        if ( atBoundary(*eit) != NOT ) ofs <<" CNXYZ";
                        ofs << endl;
                      }
                 }
             }
           break;
        case ELEMENT:
             for ( auto nit=sg.NodesBegin(); nit != sg.NodesEnd(); ++nit )
               {
                   ofs <<"PVRTX "<< (*nit)->Idx() <<" ";
                   ofs << (*nit)->x() <<" "<< (*nit)->y() <<" "<< (*nit)->z() <<" ";
                   switch( prop_key.type )
                     {  
                        case SCALAR:
                             ofs << sc_elmt_data[ (*nit)->Idx() ]() <<" ";
                          break;
                        case VECTOR:
                             for ( auto i{0U}; i<dim; i++ )
                               ofs << vc_elmt_data[ (*nit)->Idx() ][i] <<" ";
                          break;
                        case TENSOR:
                             for ( auto i{0U}; i<dim; i++ )
                               for ( auto j{0U}; j<dim; j++ )
                                 ofs << ts_elmt_data[ (*nit)->Idx() ](i,j) <<" ";
                        default:
                          throw out_of_range("GoCadInterface<dim>: variable type not handled");
                    }
                  if ( (*nit)->AtBoundary() != NOT ) ofs <<" CNXYZ";
                  ofs << endl;
              }
          break;
        default:
          cout <<"\nModel<dim>::OutputVariableToTSolid: ";
          cout <<"Property placement could not be parsed. Nothing was done."<< endl;
          return;
     }

    // 5.3 Writing plist (nodes that make up the tetrahedra)
    // -----------------------------------------------------
    for ( auto eit=sg.CellsBegin(); eit!= sg.CellsEnd(); eit++ )
      {
         ofs <<"TETRA ";
         for ( auto j{0U}; j<(*eit)->Nodes(); j++ ) ofs << (*eit)->N(j)->Idx() <<" ";
         ofs << endl;
      }
    ofs << "END"<< endl;    

    ofs.close();
    cout <<"\nModel<dim>::OutputVariableToTSolid: file '"<< file <<"' written successfully."<< endl;

 } // end OutputVariableToTSolid







/**
 
Given a three-dimensional CSMP model, OutputVariableToTSolid() creates a
GoCad TSolid dataset for the indicated variable in a subregion of the
model which was previously identified as a Region object. The variable
is to a textfile with the GoCad specific extension '.so'. The TSolid 
format is described in the GoCad manual, Appendix B (p. B-14) and consists 
of: 

 
GoCad header,
PVTRX list: node coordinate and property data listed by node number,
TETRA list: nodes per tetrahedral finite element,
END end of file specifier
  

Importantly, a TSolid file written by GoCad will contain constraint point
lists preceding each TETRA block. In constrast to CSMP jargon, a 
constraint point in GoCad is a node point which lies on an object 
boundary. In this fashion, the model volume is 
described object by object.  

@section arguments Input Arguments 

OutputVariableToTSolid() expects four arguments, the last is
the timestep in a simulation which is represented by the dataset. This
argument is optional. The first argument is the name of the file to which
the TSolid dataset is written to and to whose name the timestep is 
appended. If the timestep is not specified, a zero will be appended to
the filename before the extension '.so'. The second argument is the name 
of the group object defining the model subregion over which the variable 
(third argument) shall be output.   

@section implementation Implementation

The method uses a GocadHeader class object to write the file and the 
property class header in the TSolid file. 

If the method is used to output element properties, these will be mapped
onto the nodes, since GoCad only permits nodal variables. To achieve 
the mapping, intermediate variable vectors are created and all nodes
of each element are assigned the corresponding element property value. 
This strategy has the disadvantage that nodes with a higher ID number 
will overwrite the property values on nodes with lower numbers. The effect
may be a rugged appearance of element property boundaries. User
suggestions are welcome to overcome this deficiency. 

@section application Application

The method is used to output a property to a TSolid representing a 
subregion of a CSMP model. 

@section messages Messages 

The method will report an error and fail to write an output file, if it
cannot open the specified file with extension '.so' in text mode. 

If the target group is not defined (which it may due to a spelling error,
for instance), the method will also report this as an error. 
*/
template<uint32_t dim>
void  GoCadInterface<dim>::OutputVariableToTSolid( const Model<dim>& sgroup, 
                                                      const char* fname, 
                                                      const char* group_name, 
                                                      const char* var, 
                                                      long timestep ) const
 {
   const Region<dim>&  sg(sgroup.Region("Model"));
   sg.UpdateMemberIndexes();

   // 0. Finding the desired Region
   // ----------------------------
   const Region<dim>&  gref  = sgroup.Region( group_name );
    
   // 1. Property Access in the Model
   // ------------------------------------
   csmp::Index prop_key = sgroup.Database().StorageKey(var);

   if ( prop_key.type != SCALAR && prop_key.place == ELEMENT_INTEGRATION_POINT )
     {
        throw csmp::Exception( WARNING, "GoCadInterface<dim>::OutputVariableToTSolid", 
                        "With this method only scalar IntegrationPoint variables can be handled. Nothing is done..." ); 		       
        return;
     }

   // 3. Getting output file ready
   // ----------------------------
   ofstream  ofs;
   char      file[200], step[100];
   strcpy( file, fname ); 
   if ( timestep > 0 )
     {
        snprintf( step, sizeof(step), "%ld", timestep );
        strcat( file, step );
     }
   strcat( file, ".so" );
   
   ofs.open ( file, ios::out|ios::trunc );
   if ( !ofs )
     {
         throw csmp::Exception( FATAL_ERROR, "GoCadInterface<dim>::OutputVariableToTSolid", 
                                      "Output file could not be opened");
         return;
      }

   // 4. Writing data to file
   // -----------------------
   // 4.1 writing the file header
   // ---------------------------
   double  pmin, pmax;
   gref.MinMaxOf( var, pmin, pmax );
   char object_name[200];
   strcpy( object_name, group_name );
   if ( timestep > 0 ) strcat( object_name, step );
   GocadHeader  header( "TSolid", object_name, group_name );
   if      ( prop_key.type == SCALAR ) header.AddProperty( var, pmin, pmax );
   else if ( prop_key.type == VECTOR ) header.AddProperty( var, pmin, pmax, 99, 3 );
   else if ( prop_key.type == TENSOR ) header.AddProperty( var, pmin, pmax, 99, 9 );
   header.WriteToText( ofs );

   // 4.2 making a list of nodes that belong to the group
   // ---------------------------------------------------
   vector<size_t>    nodes(   gref.Nodes() );
   vector<size_t>    elmts(   gref.Cells() );

   size_t n(0U);
   for ( auto nit=gref.NodesBegin(); nit!=gref.NodesEnd(); nit++, n++ )    nodes[n] = (*nit)->Idx();
   for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); it++, n++ ) elmts[n] = (*it)->Idx();

   // 4.4 making a list of element properties if these are required
   // -------------------------------------------------------------
   //  All properties are output as NODE Properties.
   // Thus, if element properties are to be output one
   // must first create an intermediate array to store these 
   // and then output from this array.
   ScalarVariable              sc;
   VectorVariable<dim>          vc;
   TensorVariable<dim>          ts;
   Point<dim>                   xyz;
   vector<ScalarVariable >        sc_elmt_data;
   vector<VectorVariable<dim> >   vc_elmt_data;
   vector<TensorVariable<dim> >   ts_elmt_data;
   vector<size_t>                 node_ids;
   // node-id's of elements
   set<size_t>                    elmt_nds;
   
   if ( prop_key.place == ELEMENT ) 
     {
        if       ( prop_key.type == SCALAR ) { 
             sc_elmt_data.reserve( gref.Nodes() );
             for ( size_t i{0U}; i<gref.Nodes(); i++ ) sc_elmt_data.push_back( sc );
          }
        else if  ( prop_key.type == VECTOR ) {
             vc_elmt_data.reserve( gref.Nodes() );
             for ( size_t i{0U}; i<gref.Nodes(); i++ ) vc_elmt_data.push_back( vc );
          }
        else if  ( prop_key.type == TENSOR ) {
             ts_elmt_data.reserve(gref.Nodes() );
             for ( size_t i{0U}; i<gref.Nodes(); i++ ) ts_elmt_data.push_back( ts );
          }
        for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); it++ )
          {
              for ( auto j{0U}; j<(*it)->Nodes(); j++ ) {
                  elmt_nds.insert( (*it)->N(j)->Idx() );
                  switch ( prop_key.type )
                    {
                        case SCALAR:
                             (*it)->Read(prop_key, sc );
                             // triangles only
                             for ( auto i{0U}; i<(*it)->Nodes(); i++ ) sc_elmt_data[ (*it)->N(j)->Idx() ] = sc;
                          break;
                        case VECTOR:
                             (*it)->Read(prop_key, vc );
                             for ( auto i{0U}; i<(*it)->Nodes(); i++ ) vc_elmt_data[ (*it)->N(j)->Idx() ] = vc;
                          break;
                        case TENSOR:
                             (*it)->Read(prop_key, ts );
                             for ( auto i{0U}; i<(*it)->Nodes(); i++ ) ts_elmt_data[ (*it)->N(j)->Idx() ] = ts;
                        default:
                            cout <<"\nGoCadInterface<dim>::OutputVariableToTSolid: type of elmt. prop. not defined.\n";
                    }
             }
          }
     }

   // 5. Outputting data to GoCad text file
   // -------------------------------------
   switch( prop_key.place )
     {
        case NODE:
             for ( size_t i{0U}; i<nodes.size(); i++ )
               {
                  ofs <<"PVRTX "<< i+1 <<" ";
                  ofs << gref.N( nodes[i] )->x() <<" "<< gref.N( nodes[i] )->y() <<" "<< gref.N( nodes[i] )->z() <<" ";
                  switch( prop_key.type )
                    {  
                       case SCALAR:                   
                            gref.N( nodes[i] )->Read(prop_key, sc );
                            ofs << sc();
                         break;
                       case VECTOR:                   
                            gref.N( nodes[i] )->Read(prop_key, vc );
                            for ( auto j{0U}; j<dim; j++ )  ofs << vc[j] <<" ";
                         break;
                       case TENSOR:                   
                            gref.N( nodes[i] )->Read(prop_key, ts );
                            for ( auto j{0U}; j<dim; j++ )
                              for ( auto k=0; k<dim; k++ )
                                 ofs << ts(j,k) <<" ";
                       default:
                         throw out_of_range("GoCadInterface<dim>: variable type not handled");
                    }
                  if ( gref.N( nodes[i] )->AtBoundary() != NOT ) ofs <<" CNXYZ";
                  ofs << endl;
               } 
          break;
       case ELEMENT_INTEGRATION_POINT: {
             size_t counter{0};
             for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); it++, counter++ )
               {
                  for ( auto j{0U}; j<(*it)->IntegrationPoints(); j++ ) {
                        // getting constraint point coordinates
                        Point<dim>  pt((*it)->IntegrationPoint( j ));
                        (*it)->Read(prop_key, sc );
                        ofs <<"PVRTX "<< counter+1 <<" ";
                        ofs << pt[0] <<" "<< pt[1] <<" "<< pt[2] <<" ";
                        ofs << sc();
                        if ( atBoundary(*it) != NOT ) ofs <<" CNXYZ";
                        ofs << endl;
                    }
               } 
       }
           break;
       case ELEMENT: {
             size_t counter{0};
             for ( auto en_it=elmt_nds.begin(); en_it!=elmt_nds.end(); en_it++, counter++ )
               {
                  ofs <<"PVRTX "<< counter <<" ";
                  ofs << gref.N( (*en_it) )->x() <<" "<< gref.N( (*en_it) )->y() <<" "<< gref.N( (*en_it) )->z() <<" ";
                  switch( prop_key.type )
                    {  
                       case SCALAR:
                            ofs << sc_elmt_data[ (*en_it) ]() <<" ";
                         break;
                       case VECTOR:
                            for ( auto j{0U}; j<dim; j++ )
                               ofs << vc_elmt_data[ (*en_it) ][j] <<" ";
                         break;
                       case TENSOR:
                            for ( auto k=0; k<dim; k++ )
                              for ( auto j{0U}; j<dim; j++ )
                                 ofs << ts_elmt_data[ (*en_it) ](k,j) <<" ";
                       default:
                          throw out_of_range("GoCadInterface<dim>: variable type not handled");
                     }
                  if ( gref.N( (*en_it) )->AtBoundary() != NOT ) ofs <<" CNXYZ";
                  ofs << endl;
               }
       }
          break;
        default:
          cout <<"\nModel<dim>::OutputVariableToTSolid: ";
          cout <<"Property placement could not be parsed. Nothing was done."<< endl;
          return;
     }

    // 4.3 Writing plist (nodes that make up the tetrahedra)
    // -----------------------------------------------------
    map<size_t,size_t>  new_node_ids;

    for ( size_t i{0U}; i<nodes.size(); i++ ) new_node_ids[ nodes[i] ] = i;
    for ( size_t i{0U}; i<elmts.size(); i++ )
      {
      
  cout <<"\nGoCadInterface<dim>::OutputVariableToTSolid: this has not been properly debugged yet!"<< endl;
         ofs <<"TETRA ";
         for ( size_t j{0U}; j<4U; j++ ) 
           ofs << (*new_node_ids.find( gref.N(i)->Idx() )).second <<" ";
         ofs << endl;
      }
    ofs << "END"<< endl;    

    ofs.close();
    cout <<"\nGoCadInterface<dim>::OutputVariableToTSolid: file '"<< file <<"' written successfully."<< endl;

 } // end OutputVariableToTSolid








/**
 
Depending on whether the CSMP model is either two- or three-dimensional,
OutputRegionsToGoCadFiles() output the target variable from all Region
subregions inside the model to GoCad TSurf or GoCad TSolid files,
respectively.  

@section arguments Input Arguments 

The first argument gives the name of the output physical variable and 
the second argument depicts the current timestep in a transient CSMP
simulation. This timestep will be appended to the output filename before
the GoCad specific file extension ('.ts' for TSurf and '.so' for TSolid
files). The names of the output files will be composed of the Region names 
inside the CSP model and the name of the physical variable using a
hyphen as separator, for instance.  
 
granite-pressure1.ts
 

@section implementation Implementation

Method depends of the Region interfaces 'OutputVariableToTSurface()' and 
'OutputVariableToTSolid()'. 

@section application Application

OutputRegionsToGoCadFiles() presents a very convinient way to output 
selected subregions of a CSP model as GoCad input files for visualisation. 

@section messages Messages 

The method will report an error and exit without completing its task, if 
no Region objects were defined in the model. 
*/
template<uint32_t dim>
void GoCadInterface<dim>::OutputRegionsToGoCadFiles( const Model<dim>& sgroup, 
                                                       const char* property, 
                                                       long timestep ) const
 {
    const Region<dim>&  sg(sgroup.Region("Model"));
    sg.UpdateMemberIndexes();

    typename map<string,Region<dim> >::const_iterator  it;
    char        file_name[INFO_STRING];
    string  prop(property);
    for ( string::iterator t=prop.begin(); t!=prop.end(); t++ ) if ( *t == ' ' ) *t = '_';

    for ( it=sgroup.UniqueRegionsBegin(); it!=sgroup.UniqueRegionsEnd(); it++ )
      {
         strcpy( file_name, (*it).first.c_str() );
         strcat( file_name, "_" );
         strcat( file_name, prop.c_str() );
         if ( dim == 2 )
            OutputVariableToTSurface( sgroup, file_name, (*it).first.c_str(), property, timestep );
         else  
            OutputVariableToTSolid( sgroup, file_name, (*it).first.c_str(), property, timestep );  
      }

 } // end OutputRegionsToGoCadFiles




/**
 
Loops through plist making a map of node numbers. If this map has jumps 
in the numbering, these are detected when looping through it again. 
 */
template<uint32_t dim>
bool GoCadInterface<dim>::VerifyConsecutiveNodeNumbering( map<size_t,vector<int64_t> >& plist )
 const 
 {
    map<size_t,vector<int64_t> >::const_iterator  it;
    vector<int64_t>::const_iterator                  nit;
    set<int64_t>                                     node_numbers;
    set<int64_t>::const_iterator                     sit;
    size_t                                          counter(1);
    
    for ( it=plist.begin(); it!=plist.end(); it++ )
      for ( nit=(*it).second.begin(); nit!=(*it).second.end(); nit++ ) 
        node_numbers.insert( *nit );
          
    for ( sit=node_numbers.begin(); sit!=node_numbers.end(); sit++ )
      // if a gap in the numbering leads to a difference in incremented counter
      // and the ordered node numbers in the list, false is returned.
      if ( *sit != counter++ ) 
        {
           cout <<"\nGoCadInterface::VerifyConsecutiveNodeNumbering:"<< endl;
           cout <<"\n     Node numbering gap in 'plist': Current node: " << endl << endl;
           cout <<"     "<< *sit     <<" versus "<< counter <<", last correct node: ";
           cout << *(--sit) <<" versus "<< counter-2 << endl << endl;
           // to continue tracking
           counter = *(++sit);
           return false;
        }      
    
    cout <<"\nGoCadInterface::VerifyConsecutiveNodeNumbering: Yes, nodes are numbered consecutively."<< endl;
    return true;
  
 } // end VerifyConsecutiveNodeNumbering   






template class GoCadInterface<2U>;
template class GoCadInterface<3U>;

} // csmp


