#include "RhinoSurfaceReader.h"
#include "LinearTriangle3D.h"
#include "Box.h"
#include "VSet.h"

using namespace std;


namespace csmp {

SKM_RhinoSurfaceReader::SKM_RhinoSurfaceReader()
 {
 }
 
 
SKM_RhinoSurfaceReader::SKM_RhinoSurfaceReader( const char* raw_file )
 {
    InitializeFrom( raw_file, false );
 }
 
 
SKM_RhinoSurfaceReader::SKM_RhinoSurfaceReader( const SKM_RhinoSurfaceReader& r )
 {
    *this = r;
 }
 


SKM_RhinoSurfaceReader::~SKM_RhinoSurfaceReader()
 {
 }
 
 
SKM_RhinoSurfaceReader&  SKM_RhinoSurfaceReader::operator=( const SKM_RhinoSurfaceReader& r )
 {
    if ( &r != this ) objects = r.objects;
    return *this;
 }



bool SKM_RhinoSurfaceReader::InitializeFrom( const char* raw_file, bool erase_old )
 {
    // access Rhino ".raw" file
    char                   fname[200], intext[1000];
    const char* const      delims =" ,\t";
    char*                  token;
    mjl::Point3D           p1, p2, p3;
    list<mjl::Triangle3D>  triangles;
    double               x, y, z;
    int64_t                  triangle_counter(0);
    string                 oname;
    
    strcpy( fname, raw_file );
    strcat( fname, ".raw" );
    ifstream  ifs( fname );
    assert( ifs.is_open() );
    
    if ( erase_old )
      {
         cout <<"\nSKM_RhinoSurfaceReader::InitializeFrom: ";
         cout <<"Erasing old objects before reading new ones."<< endl;
         EraseObjects();
      }
    // reading Rhino surface object by object
    while ( !ifs.eof() )
      {
         // 0. each triangle's coordinates occupy a single line in the Rhino file 
         // reading line:  "spec-name  state\n"
         ifs.getline( intext, 200 );
         // exit condition
         if ( strlen( intext ) <= 1 ) break;
         
         // 1. If the line is a single alphanumeric expression,
         //    it must be the object name
         if ( isalpha(intext[0]) ) 
           {
              // if this name marks the first line after an object
              // this object is packed away.
              if ( triangle_counter != 0 )
                {
                   // the earlier object name is used as map key
                   objects[ oname ] = triangles;
                   triangles.erase( triangles.begin(), triangles.end() );
                   triangle_counter = 0;
                }
              // the new object name is parsed
              oname = strtok( intext, delims );
           }
         // 2. else, the line must hold the coordinates of a triangle
         //    these are 9 floating point numbers
         else 
           {
               // reading coordinates of Vertex 1
               token = strtok( intext, delims );
               x = atof( token );
               token = strtok( NULL, delims );
               y = atof( token );
               token = strtok( NULL, delims );
               z = atof( token );
               p1.Set( x, y, z );

                // reading coordinates of Vertex 1
               token = strtok( NULL, delims );
               x = atof( token );
               token = strtok( NULL, delims );
               y = atof( token );
               token = strtok( NULL, delims );
               z = atof( token );
               p2.Set( x, y, z );

               // reading coordinates of Vertex 1
               token = strtok( NULL, delims );
               x = atof( token );
               token = strtok( NULL, delims );
               y = atof( token );
               token = strtok( NULL, delims );
               z = atof( token );
               p3.Set( x, y, z );
               
               // storing triangle and counting triangles per object
               triangles.push_back( mjl::Triangle3D( p1, p2, p3, triangle_counter++ ) ); 
           }
       }

     // 3. Storing the last surface in the file as well
     if ( triangle_counter != 0 )
       {
          // the earlier object name is used as map key
          objects[ oname ] = triangles;
          triangles.erase( triangles.begin(), triangles.end() );
          triangle_counter = 0;
       }
    
    // 4. reporting the outcome of the reading
    if ( objects.empty() )
      cout <<"\nSKM_RhinoSurfaceReader::InitializeFrom: File: "<< fname <<" No objects were read.";
    else
      {
         cout <<"\nSKM_RhinoSurfaceReader::InitializeFrom: File: "<< fname <<" read successfully.";   
         cout <<"\n\tread "<< objects.size() <<" objects with following names and numbers of triangles:\n";
         for ( map<string,list<mjl::Triangle3D> >::const_iterator
               it=objects.begin(); it!=objects.end(); it++ )
           cout <<"\t\t'"<< (*it).first <<"':  "<< (*it).second.size() << endl;
         cout << endl; 
      }  
    // if data were read from file, an O.K. is returned   
    return !objects.empty();
 
 } // end 
 



void SKM_RhinoSurfaceReader::ObjectToPData( const string& obj_name,
                                            map<size_t,mjl::Point3D>&  points,
                                            map<size_t,vector<int64_t> >& plist,
                                            size_t poffset ) const
  {
     // finding the desired object in the map 
     map<string,list<mjl::Triangle3D> >::const_iterator it = objects.find(obj_name);  
     if ( it == objects.end() ) {
          cout <<"\nSKM_RhinoSurfaceReader::ObjectToPData: Object could not be found. ";
          cout <<"Nothing was done."<< endl;
          return;
      }
    cout <<"\nSKM_RhinoSurfaceReader::ObjectToPData: Converting object: "<< (*it).first << endl;  
    // zapping the supplied maps
    points.erase( points.begin(), points.end() );
    plist.erase( plist.begin(), plist.end() ); 
     
    // zooming through the triangle list, 
    // creating a unique set of nodes ordered by hashkeys created from the node coordinates
    // for this, the node coordinates are converted to strings 
    map<string,vector<double> >          pxyz;
    vector<double>                       coord(3);
    list<mjl::Triangle3D>::const_iterator  lit;
    string                                 key;
    char                                   num[20];
    int32_t                                  i, j;
     
    for ( lit=(*it).second.begin(); lit!=(*it).second.end(); lit++ )
      {
         // for all triangle verticies
         for ( i=0; i<3; i++ )
           {
              key.erase( key.begin(), key.end() );
              // for all the coordinate directions
              for ( j=0; j<3; j++ )
                {
                   coord[static_cast<uint32_t>(j)] = (*lit)[i][j];
                   // convert double coordinate value to string and add to hash key
                   sprintf( num, "%lf", (*lit)[i][j] );
                   key += num;
                }
              // adding new map entry
              pxyz[ key ] = coord;
           }
      }
      
    // creating a second numbered (0...n-1) map of vertex coordinates
    // and assigning the vertex points to the output map
    map<string,size_t>                              pxyz_ids;
    size_t                                          n;
    vector<int64_t>                                  ids(3);
    map<string,vector<double> >::const_iterator  pit;
    
    for ( n=poffset, pit=pxyz.begin(); pit!=pxyz.end(); pit++ )
      { 
         points[ n ] = mjl::Point3D( (*pit).second[0], (*pit).second[1], (*pit).second[2] );
         pxyz_ids[ (*pit).first ] = n++;
      }
   
   
    // creating the plist by searching the ID's for each vertex in the vertex list
    map<string,size_t>::const_iterator  cit;
    
    for ( n=poffset, lit=(*it).second.begin(); lit!=(*it).second.end(); lit++ )
      {
         // for all triangle verticies
         for ( i=0; i<3; i++ )
           {
              key.erase( key.begin(), key.end() );
              // for all the coordinate directions create key to search the map
              // --------------------------------------------------------------
              for ( j=0; j<3; j++ )
                {
                   coord[static_cast<uint32_t>(j)] = (*lit)[i][j];
                   // convert double coordinate value to string and add to hash key
                   sprintf( num, "%lf", (*lit)[i][j] );
                   key += num;
                }
              // search the ID map for the node number
              // -------------------------------------
              if ( (cit=pxyz_ids.find( key )) == pxyz_ids.end() )
                {
                   cout <<"\nSKM_RhinoSurfaceReader::ObjectToPData: ";
                   cout <<" Unable to find a vertex for the key: "<< key << endl;
                   cout <<"\nTerminating execution of ObjectToPData()." << endl;
                   return;
                }
              else ids[static_cast<uint32_t>(i)] = (*cit).second; 
           }
         // assigning entry to plist
         plist[ n++ ] = ids;
      }
     
  } // end ObjectToPData






void SKM_RhinoSurfaceReader::CreateNeighborPData( const map<size_t,vector<int64_t> >& plist,
                                                  map<size_t,vector<int64_t> >& pfverts,
                                                  vector<std::int8_t>& pbflags ) const
 {
     //  parent element id,  edge of p1 < p2
     map<pair<size_t,size_t>,int64_t> edge_map;
     
     // just in case
     if ( !pfverts.empty() ) pfverts.erase( pfverts.begin(), pfverts.end() );
     if ( !pbflags.empty() ) pbflags.erase( pbflags.begin(), pbflags.end() );
     
     // making duplicate edge list from plist where edges of opposite elements
     // have opposite node id numbers
     for ( auto pit=plist.begin(); pit!=plist.end(); pit++ ) {
          // edge 1 (01) counter-clockwise nodes
          edge_map.insert( make_pair( make_pair((*pit).second[0], (*pit).second[1]), static_cast<int32_t>((*pit).first) ) );
          // edge 2 (12)
          edge_map.insert( make_pair( make_pair((*pit).second[1], (*pit).second[2]), static_cast<int32_t>((*pit).first) ) );
          // edge 3 (20)
          edge_map.insert( make_pair( make_pair((*pit).second[2], (*pit).second[0]), static_cast<int32_t>((*pit).first) ) );
       }
 
     // finding the neighbor id's by using the opposite edge definition now to retrieve the neighbor edges
     // if there is no such neighbor, the edge is located on the model boundary and is flagged as IRREGULAR_OUTSIDE
     map<pair<size_t,size_t>,int64_t>::const_iterator  eit;
     vector<int64_t> nbors(3);

     for ( auto pit=plist.begin(); pit!=plist.end(); pit++ ) {
          // initializing neighbor ids
          // neighbor 1 (12)                              clockwise nodes 
          if ( (eit=edge_map.find( make_pair((*pit).second[2],(*pit).second[1]))) == edge_map.end() ) {
               nbors[0] = IRREGULAR_OUTSIDE;
               pbflags[ (*pit).second[2] ] = nbors[0];
               pbflags[ (*pit).second[1] ] = nbors[0];
            }
          else nbors[0] = (*eit).second;
          // neighbor 2 (20)
          if ( (eit=edge_map.find( make_pair((*pit).second[0],(*pit).second[2]))) == edge_map.end() ) {
               nbors[1] = IRREGULAR_OUTSIDE;
               pbflags[ (*pit).second[0] ] = nbors[1];
               pbflags[ (*pit).second[2] ] = nbors[1];
            }
          else nbors[1] = (*eit).second;
          // neighbor 3 (01)
          if ( (eit=edge_map.find( make_pair((*pit).second[1],(*pit).second[0]))) == edge_map.end() ) {
               nbors[2] = IRREGULAR_OUTSIDE;
               pbflags[ (*pit).second[1] ] = nbors[2];
               pbflags[ (*pit).second[0] ] = nbors[2];
            }
          else nbors[2] = (*eit).second;
          // adding new neighbor vector to map
          pfverts[ (*pit).first ] = nbors;
       }
 
 
 } // end CreateNeighborPData










size_t  SKM_RhinoSurfaceReader::SelectAndDescribeObjects( list<pair<string,string> >& selections,
                                                          bool erase_list_before ) const
 {
    if ( objects.empty() )
      {
         cout <<"\nSKM_RhinoSurfaceReader::SelectObjects: Error: No objects to select from."<< endl;
         return 0;
      }
    if ( erase_list_before ) selections.erase( selections.begin(), selections.end() );
 
    cout <<"\nSKM_RhinoSurfaceReader::SelectObjects: Please select objects from listing: "<< endl;
    // looping throug the objects and writing them to the output file
    map<string,list<mjl::Triangle3D> >::const_iterator  it;
    int answ(0);
    
    for ( it=objects.begin(); it!=objects.end(); it++ )
      {
         // selecting an object
         cout <<"\nObject: "<< (*it).first <<" (yes=1,no=0): "; 
         cin >> answ;
         
         // if an object is selected its geological type is queried
         if ( answ == 1 ) 
           {
              cout <<"\nGoCad 'GEOLOGICAL_TYPE', options: "<< endl;
              cout <<"  top (0)\n  intraformational (1)\n  fault (2)\n  unconformity (3)\n";
              cout <<"  intrusive(4)\n  topography (5)\n  boundary (6)\n  ghost (7)"<< endl; 
              cin >> answ;
              
              // storing the object
              switch ( answ )
                {
                   case 0: 
                     selections.push_back( pair<string,string>((*it).first, "top") );
                     break;
                   case 1: 
                     selections.push_back( pair<string,string>((*it).first, "intraformational") );
                     break;
                   case 2: 
                     selections.push_back( pair<string,string>((*it).first, "fault") );
                     break;
                   case 3: 
                     selections.push_back( pair<string,string>((*it).first, "unconformity") );
                     break;
                   case 4: 
                     selections.push_back( pair<string,string>((*it).first, "intrusive") );
                     break;
                   case 5: 
                     selections.push_back( pair<string,string>((*it).first, "topography") );
                     break;
                   case 6: 
                     selections.push_back( pair<string,string>((*it).first, "boundary") );
                     break;
                   case 7: 
                     selections.push_back( pair<string,string>((*it).first, "ghost") );
                }   
              answ = 0;
           }
      }  
     
    cout <<"\nThank you."<< endl;
 
    return selections.size();
 }




/**
 
Creates a map of unique points in the object, labeled (keys) by their point
number. A map is also created which contains entries of the node numbers
which make up each triangle.   
*/
bool SKM_RhinoSurfaceReader::PopObject( const char *obj_name,
                                        map<size_t,mjl::Point3D >&  points,
                                        map<size_t,vector<int64_t> >& plist,
                                        size_t poffset ) const
 {
     ObjectToPData( string(obj_name), points, plist, poffset );
     
     return true;
 }





/**
 
Puts a header into the GoCad TSurf file. This header specifies the 
defaults for the colors and properties with which the surface will
be drawn.  
*/
void SKM_RhinoSurfaceReader::WriteGocadHeader( const char* surf_name, 
                                               const char* GEOLOGICAL_TYPE, 
                                               ofstream& ofs ) const
 {
    assert( ofs.is_open() );
  // writing GoCad text file header
  // ------------------------------
  ofs <<"GOCAD TSurf 0.001"<< endl;
  ofs <<"HEADER {"<< endl;
  ofs <<"name:"<< surf_name << endl;
  ofs <<"*solid:true"<< endl;
  ofs <<"*painted:true"<< endl;
  // the material is used to discern surfaces by their names
  ofs <<"*painted*variable:material"<< endl;
  ofs <<"*mesh*color:gray30"<< endl;
  ofs <<"*cn:true"<< endl;
  ofs <<"*cn*size:0.2"<< endl;
  ofs <<"*cn*color:white"<< endl;
  ofs <<"*under_threshold:false"<< endl;
  ofs <<"*border_only:false"<< endl;
  ofs <<"*color_coded:false"<< endl;
  ofs <<"*normals:true"<< endl;
  ofs <<"*normals:true"<< endl;
  ofs <<"*solid*transparency:0"<< endl;
  ofs <<"*shrink_coef:1"<< endl;
  ofs <<"*solid*specular:gray70"<< endl;
  ofs <<"}"<< endl;
  
  // specifying the geological object type 
  ofs <<"GEOLOGICAL_TYPE "<< GEOLOGICAL_TYPE << endl;
  
  /* upon a time when properties will also be transmitted to GoCad
  
  ofs <<"PROPERTIES material"<< endl;
  ofs <<"NO_DATA_VALUES 1e-30"<< endl;
  ofs <<"PROPERTY_CLASSES material"<< endl;
  ofs <<"ESIZES 1"<< endl; // one-dimensional material variable
  
  // property class header (how property shall be displayed
  ofs <<"PROPERTY_CLASS_HEADER material {"<< endl;
  ofs <<"*low_clip:0.0"<< endl;
  ofs <<"*high_clip:"<< objects.size() << endl;
  ofs <<"*pclip:99"<< endl;
  ofs <<"*colormap:rainbow"<< endl;
  ofs <<"*colormap*transparency:true"<< endl;
  ofs <<"*colormap*low_clip_transparent:false"<< endl;
  ofs <<"*colormap*transparency_min:0"<< endl;
  ofs <<"*colormap*transparency_max:1"<< endl;
  ofs <<"*colormap*high_clip_transparent:false"<< endl;
  ofs <<"*colormap*nodata:false"<< endl;
  ofs <<"*cnp*symbol:point"<< endl;
  ofs <<"}"<< endl;
  */
  
 } // end WriteGocadHeader




// output object to VSet
void SKM_RhinoSurfaceReader::OutputObjectTo( const char* obj, VSet<3U>& vset ) const
 {
    if ( objects.empty() ) {
         cerr <<"\nSKM_RhinoSurfaceReader::OutputObjectTo: No objects are present.";
         return;
      }
   
     // looping throug the objects and writing them to the output file
    string  object(obj);
    map<string,list<mjl::Triangle3D> >::const_iterator  it=objects.find( object );
    if ( it == objects.end() ) {
         cerr <<"\nSKM_RhinoSurfaceReader::OutputObjectTo: Desired object could not be found."<< endl;
         return;
      }  

    // find object and create pxyz and plist arrays for the desired object
    map<size_t,mjl::Point3D>     points;
    map<size_t,vector<int64_t> >  plist;

    PopObject( object.c_str(), points, plist );
    
    // create 'pfverts' data
    vector<std::int8_t>         pbflags;
    map<size_t,vector<int64_t> > pfverts;

    CreateNeighborPData( plist, pfverts, pbflags );

    // storing data in 'vset'
    LinearTriangle3D  tri3;
    vset.Erase();
    vset.Resize( tri3.Nodes(),
                 tri3.Neighbors(),
                 tri3.ElementType(),
                 points.size(), plist.size() );
    
    // finite element type
    vset.SingleElementType( LINEAR_TRIANGLE3D );

    // px, py, pz
    for ( auto ptit=points.begin(); ptit!=points.end(); ptit++ )
      {
         vset.Px( (*ptit).first, (*ptit).second.x_ );
         vset.Py( (*ptit).first, (*ptit).second.y_ );
         vset.Pz( (*ptit).first, (*ptit).second.z_ );
      }

    // boundary flags and data
    vset.AddBFlags( pbflags.begin(), pbflags.end() );
    
    // plist, pfverts
    vset.AddPlist( plist.begin(), plist.end() );
    vset.AddPfverts( pfverts.begin(), pfverts.end() );

 } // end OutputObjectTo









/** Writes all objects from the Rhino output file into a single Gocad file.


@section application Application

Not necessarily recommended since GoCad often has difficulties in reading
individual surfaces.  
 */
void SKM_RhinoSurfaceReader::WriteObjectsToTSurf( const char* tsurf_file, 
                                                  const char* GEOLOGICAL_TYPE ) const
 {
   ofstream  ofs;
   char      file[200];
   strcpy( file, tsurf_file ); 
   strcat( file, ".ts" );

   ofs.open ( file, ios::out|ios::trunc );
   if ( !ofs )
     {
         cout <<"\nSKM_RhinoSurfaceReader::WriteObjectsToTSurf: "; 
         cout <<"Output file could not be opened"<< endl;
         return;
      }

   if ( objects.empty() )
     {
        cout <<"\nSKM_RhinoSurfaceReader::WriteObjectsToTSurf: No objects are present.";
        return;
     }
   
  WriteGocadHeader( tsurf_file, GEOLOGICAL_TYPE, ofs ); 
         
  // looping throug the objects and writing them to the output file
  map<string,list<mjl::Triangle3D> >::const_iterator  it;
    
  for ( it=objects.begin(); it!=objects.end(); it++ )
    {
       cout <<"\nSaving object: "<< (*it).first <<" to file."; 
       WriteObjectToTSurf( (*it).first.c_str(), ofs ); 
    }  
  // terminating object description
  ofs <<"END" << endl;
    
  ofs.close();
  cout <<"\nSKM_RhinoSurfaceReader::WriteObjectsToTSurf: "<< file <<" written successfully"<< endl;
}  
        




/**
 
Writes a series of files which contain one surface each corresponding to
the selected surfaces from the Rhino '.raw' dataset. The names of the 
objects are appended to the output file names which start with the
supplied file name. A hyphen is put between the filename and the 
object name.  

*/
void  SKM_RhinoSurfaceReader::WriteSelectedObjectsToTSurf( const char* tsurf_file ) const    
 {
   if ( objects.empty() )
     {
        cout <<"\nSKM_RhinoSurfaceReader::WriteSelectedObjectsToTSurf: No objects are present.";
        return;
     }
   
   // 1. selecting objects for output to file
   list<pair<string,string> >                  selections;
   list<pair<string,string> >::const_iterator  it;
   
   if ( SelectAndDescribeObjects( selections, false ) == 0 ) return;

   // 2. looping throug the selected objects and writing them to the output file
   //    preparing output file 
   ofstream  ofs;
   char      file[200];
    
  for ( it=selections.begin(); it!=selections.end(); it++ )
    {
       // setting the file name to be the input name plus the appended
       // object name
       strcpy( file, tsurf_file ); 
       strcat( file, "-" );
       strcat( file, (*it).first.c_str() );
       strcat( file, ".ts" );

       // opening the output file
       ofs.open ( file, ios::out|ios::trunc );
       if ( !ofs )
         {
             cout <<"\nSKM_RhinoSurfaceReader::WriteSelectedObjectsToTSurf: "; 
             cout <<"Output file could not be opened"<< endl;
             return;
         }

       // writing the file header  
       WriteGocadHeader( (*it).first.c_str(), (*it).second.c_str(), ofs ); 
         
       // writing object to file 
       cout <<"\nSaving object '"<< (*it).first <<"' in file: "<< file << endl; 
       WriteObjectToTSurf( (*it).first.c_str(), ofs ); 

       // terminating object description
       ofs <<"END" << endl;
       ofs.close();
    }  
    
  cout <<"\nSKM_RhinoSurfaceReader::WriteSelectedObjectsToTSurf: '"<< tsurf_file;
  cout <<"' file series written successfully"<< endl;
  
 } // end WriteSelectedObjectsToTSurf  
               
               
              
              
                     

/** The surface will be written to file without associating any properties
with it.  

*/
void SKM_RhinoSurfaceReader::WriteObjectToTSurf( const char* obj, ofstream& ofs ) const
 {
    assert( ofs.is_open() );
 
    map<size_t,mjl::Point3D>                     points;
    map<size_t,vector<int64_t> >                  plist;
    map<size_t,mjl::Point3D>::const_iterator     it;
    map<size_t,vector<int64_t> >::const_iterator  pit;

    // find object and create pxyz and plist arrays for the desired object
    PopObject( obj, points, plist, 0 );

    // write object to TSurf format into supplied file stream
    ofs <<"TFACE" << endl;
    
    // write vertex coordinates
    for ( it=points.begin(); it!=points.end(); it++ )
      {
         // identifier     node ID 1...n
         ofs <<"VRTX "<< (*it).first+1 <<" ";
         ofs << (*it).second.X() <<" "<< (*it).second.Y() <<" "<< (*it).second.Z() << endl;
      }

    // writing node id's per triangle as in plist
    for ( pit=plist.begin(); pit!=plist.end(); pit++ )
      {
         ofs <<"TRGL ";
         // again node ids must be augmented by 1 since Gocad counts 1...n
         for ( size_t i=0; i<3; i++ ) ofs << (*pit).second[i]+1 <<" ";
         ofs << endl;
      } 

    // terminating object description
    // ofs <<"END" << endl;

} // end WriteObjectToTSurf
     




size_t SKM_RhinoSurfaceReader::Objects() const
  { return objects.size(); }
  




void SKM_RhinoSurfaceReader::EraseObjects() 
  { objects.erase( objects.begin(), objects.end() ); }            




/**
 
Swaps the specified coordinate axis which are assumed to be labeled
0=x, 1=y, and 2=z. 

@section arguments Input Arguments 

The two integer arguments refer to the axis which shall be changed.  

@section application Application

To adapt surface data to another frame of reference. 

@section messages Messages 

If the axes do not exist (are below 0 or greater than 2) the method
will quit emitting a message. 
*/
void SKM_RhinoSurfaceReader::ExchangeCoordinateAxes( int axis_a, int axis_b )
 {
    if ( objects.empty() )
      {
         cout <<"\nSKM_RhinoSurfaceReader::ExchangeCoordinateAxes: No object data are present. ";
         cout <<"Nothing was done..."<< endl;
         return;
      }
    if ( axis_a < 0 || axis_a > 2 || axis_b < 0 || axis_b > 2 )
      {
         cout <<"\nSKM_RhinoSurfaceReader::ExchangeCoordinateAxes: Target axis does not exist: ";
         cout << axis_a <<" or "<< axis_b << endl;
         return;
      }

    // looping over the objects 
    mjl::Point3D  pt[3];
    double      swap;
    int64_t         pid;
    
    for ( map<string,list<mjl::Triangle3D> >::iterator
          oit=objects.begin(); oit!=objects.end(); oit++ )
      for ( list<mjl::Triangle3D>::iterator
            tit=(*oit).second.begin(); tit!=(*oit).second.end(); tit++ )
        {
           for ( int32_t i=0; i<3; i++ )
             {
                // changing point by point 
                // 1. getting the point 
                pt[i]         = (*tit)[i];
                // 2. changing the point coordinates
                swap          = pt[i][axis_a];
                pt[i](axis_a) = pt[i][axis_b];
                pt[i](axis_b) = swap;
             }
           // storing the reworked points re-organising the triangle if necessary
           pid = (*tit).id_;
           (*tit).Set( pt[0], pt[1], pt[2], pid );
       }

 } // end ExchangeCoordinateAxes


 
 
     
/**
 
Muliplies the triangle coordinates by the supplied factors and reorganises
the triangles if their sense of rotation changes. 

@section arguments Input Arguments 

The scale factors for the x, y, and z axis are entered. A factor of 
1. is to be entered if an axis shall not be scaled and a factor of -1.
if an axis shall be inverted. 

@section application Application

To scale an object for a new frame of reference.  

@section messages Messages 

Nothing is done if no objects are present. 
 */
void SKM_RhinoSurfaceReader::ScaleCoordinates( double x_fac, double y_fac, double z_fac )
 {
    if ( objects.empty() )
      {
         cout <<"\nSKM_RhinoSurfaceReader::ScaleCoordinates: No object data are present. ";
         cout <<"Nothing was done..."<< endl;
         return;
      }

    // looping over the objects 
    mjl::Point3D  pt[3];
    
    for ( map<string,list<mjl::Triangle3D> >::iterator oit = objects.begin(); oit!=objects.end(); oit++ )
      for ( list<mjl::Triangle3D>::iterator tit=(*oit).second.begin(); tit!=(*oit).second.end(); tit++ )
        {
           for ( int32_t i=0; i<3; i++ )
             {
                // changing point by point 
                // 1. getting the point 
                pt[i] = (*tit)[i];
                // 2. scaling the point coordinates
                pt[i](0) *= x_fac;
                pt[i](1) *= y_fac;
                pt[i](2) *= z_fac;
             }
           // storing the reworked points re-organising the triangle if necessary
           const long pid = (*tit).id_;
           (*tit).Set( pt[0], pt[1], pt[2], pid );
       }

 } // end TranslateCoordinates
 
 
 
 
     
/** Moves the triangle coordinates by the specified amounts in the x, y, and
z, directions. 

@section arguments Input Arguments 

The three method arguments specify the amount by which each point coordinate
is to be moved. A value of zero leaves the coordinates in place. 

@section application Application

To transform the stored surfaces to another space position. 

@section messages Messages 

If no objects are stored, Nothing is done.  
*/
void SKM_RhinoSurfaceReader::MoveCoordinates( double x_move, double y_move, double z_move )
 {
    if ( objects.empty() )
      {
         cout <<"\nSKM_RhinoSurfaceReader::MoveCoordinates: No object data are present. ";
         cout <<"Nothing was done..."<< endl;
         return;
      }

    // looping over the objects 
    map<string,list<mjl::Triangle3D> >::iterator  oit;  
    list<mjl::Triangle3D>::iterator                            tit;
    mjl::Point3D                                               pt[3];
    int32_t                                                     pid;
    
    for ( oit = objects.begin(); oit!=objects.end(); oit++ )
      for ( tit=(*oit).second.begin(); tit!=(*oit).second.end(); tit++ )
        {
           for ( int32_t i=0; i<3; i++ )
             {
                // changing point by point 
                // 1. getting the point 
                pt[i] = (*tit)[i];
                // 2. moving the point coordinates
                pt[i](0) += x_move;
                pt[i](1) += y_move;
                pt[i](2) += z_move;
             }
           // storing the reworked points re-organising the triangle if necessary
           pid = (*tit).id_;
           (*tit).Set( pt[0], pt[1], pt[2], pid );
       }

 } // end MoveCoordinates

} // end namespace csp




















