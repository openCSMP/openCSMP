// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "RhinoSurfaceReader.h"
#include "LinearTriangle3D.h"
#include "VSet.h"
#include "Box.h"

using namespace std;


namespace csmp {

RhinoSurfaceReader::RhinoSurfaceReader()
 {
 }
 
 
RhinoSurfaceReader::RhinoSurfaceReader( const char* raw_file )
 {
    InitializeFrom( raw_file, false );
 }
 
 
RhinoSurfaceReader::RhinoSurfaceReader( const RhinoSurfaceReader& r )
 {
    *this = r;
 }
 


RhinoSurfaceReader::~RhinoSurfaceReader()
 {
 }
 
 
RhinoSurfaceReader&  RhinoSurfaceReader::operator=( const RhinoSurfaceReader& r )
 {
    if ( &r != this ) objects = r.objects;
    return *this;
 }



#include <sstream>
#include <algorithm>

bool RhinoSurfaceReader::InitializeFrom(const char* raw_file, bool erase_old)
{
    std::string fname = std::string(raw_file) + ".raw";
    std::ifstream ifs(fname);
    
    if (!ifs.is_open()) {
        std::cerr << "Error: Could not open Rhino file: " << fname << std::endl;
        return false;
    }

    if (erase_old) {
        std::cout << "\nRhinoSurfaceReader::InitializeFrom: Erasing old objects." << std::endl;
        EraseObjects();
    }

    std::string line;
    std::string current_oname = "Default_Object";
    std::vector<LightWeightTriangle> triangles;
    int64_t triangle_counter = 0;

    while (std::getline(ifs, line)) {
        // Clean up leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        if (line.empty()) continue;

        // 1. If the line starts with a letter, it's a new object name
        if (std::isalpha(line[0])) {
            // Save the previous object before starting a new one
            if (!triangles.empty()) {
                objects[current_oname] = triangles;
                triangles.clear();
                triangle_counter = 0;
            }
            
            // Extract the first word as the name
            std::stringstream ss(line);
            ss >> current_oname;
        }
        // 2. Otherwise, parse the 9 coordinates (3 points)
        else {
            std::stringstream ss(line);
            double coords[9];
            bool success = true;

            for (int i = 0; i < 9; ++i) {
                if (!(ss >> coords[i])) {
                    success = false;
                    break;
                }
            }

            if (success) {
                Point<3> p1(coords[0], coords[1], coords[2]);
                Point<3> p2(coords[3], coords[4], coords[5]);
                Point<3> p3(coords[6], coords[7], coords[8]);

                // Replacing mjl::Triangle3D with your internal Triangle representation
                triangles.push_back(LightWeightTriangle(p1, p2, p3, triangle_counter++));
            }
        }
    }

    // 3. Store the final object
    if (!triangles.empty()) {
        objects[current_oname] = triangles;
    }

    // 4. Reporting
    if (objects.empty()) {
        std::cout << "\nRhinoSurfaceReader::InitializeFrom: " << fname << " - No objects found.";
    } else {
        std::cout << "\nRhinoSurfaceReader::InitializeFrom: " << fname << " read successfully.";
        std::cout << "\n\tRead " << objects.size() << " objects:\n";
        
        for (const auto& pair : objects) {
            std::cout << "\t\t'" << pair.first << "': " << pair.second.size() << " triangles\n";
        }
    }

    return !objects.empty();
} 




void RhinoSurfaceReader::ObjectToPData(const string& obj_name,
                                       map<size_t, Point<3>>& points,
                                       map<size_t, vector<size_t>>& plist,
                                       size_t poffset) const
{
    // 1. Find the object
    auto it = objects.find(obj_name);
    if (it == objects.end()) {
        cout << "\nRhinoSurfaceReader::ObjectToPData: Object '" << obj_name << "' not found." << endl;
        return;
    }

    cout << "\nRhinoSurfaceReader::ObjectToPData: Converting object: " << it->first << endl;

    points.clear();
    plist.clear();

    // Map to store unique coordinate strings to their assigned ID
    // Using a string key acts as a "poor man's" fuzzy-logic unifier for doubles
    map<string, size_t> pxyz_to_id;
    size_t next_point_id = poffset;
    size_t face_counter = poffset;

    // A helper to generate the coordinate key (replaces snprintf logic)
    auto generate_key = [](const Point<3>& p) {
        char buf[128];
        // Use %.6f or similar to ensure slight precision variations don't create duplicate points
        snprintf(buf, sizeof(buf), "%.6f,%.6f,%.6f", p[0], p[1], p[2]);
        return string(buf);
    };

    // 2. Single pass through the triangles
    for (const auto& tri : it->second) {
        vector<size_t> face_indices(3);

        for (uint32_t i = 0; i < 3; ++i) {
            // tri[i] returns the i-th Point<3> of the triangle
            string key = generate_key(tri[i]);

            auto existing = pxyz_to_id.find(key);
            if (existing == pxyz_to_id.end()) {
                // New unique point found
                size_t new_id = next_point_id++;
                pxyz_to_id[key] = new_id;
                points[new_id] = tri[i];
                face_indices[i] = new_id;
            } else {
                // Shared point found
                face_indices[i] = existing->second;
            }
        }
        
        // Store the face (connectivity)
        plist[face_counter++] = face_indices;
    }

    cout << "\tProcessed " << it->second.size() << " triangles into " 
         << points.size() << " unique vertices." << endl;
}





void RhinoSurfaceReader::CreateNeighborPData( const map<size_t,vector<size_t> >& plist,
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
               pbflags[ (*pit).second[2] ] = static_cast<int8_t>(nbors[0]);
               pbflags[ (*pit).second[1] ] = static_cast<int8_t>(nbors[0]);
            }
          else nbors[0] = (*eit).second;
          // neighbor 2 (20)
          if ( (eit=edge_map.find( make_pair((*pit).second[0],(*pit).second[2]))) == edge_map.end() ) {
               nbors[1] = IRREGULAR_OUTSIDE;
               pbflags[ (*pit).second[0] ] = static_cast<int8_t>(nbors[1]);
               pbflags[ (*pit).second[2] ] = static_cast<int8_t>(nbors[1]);
            }
          else nbors[1] = (*eit).second;
          // neighbor 3 (01)
          if ( (eit=edge_map.find( make_pair((*pit).second[1],(*pit).second[0]))) == edge_map.end() ) {
               nbors[2] = IRREGULAR_OUTSIDE;
               pbflags[ (*pit).second[1] ] = static_cast<int8_t>(nbors[2]);
               pbflags[ (*pit).second[0] ] = static_cast<int8_t>(nbors[2]);
            }
          else nbors[2] = (*eit).second;
          // adding new neighbor vector to map
          pfverts[ (*pit).first ] = nbors;
       }
 
 
 } // end CreateNeighborPData










size_t RhinoSurfaceReader::SelectAndDescribeObjects(list<pair<string, string>>& selections,
                                                    bool erase_list_before) const
{
    if (objects.empty()) {
        cout << "\nRhinoSurfaceReader::SelectObjects: Error: No objects to select from." << endl;
        return 0;
    }

    if (erase_list_before) {
        selections.clear();
    }

    // Define the geological types in an array for easy indexing
    static const vector<string> geo_types = {
        "top", "intraformational", "fault", "unconformity",
        "intrusive", "topography", "boundary", "ghost"
    };

    cout << "\nRhinoSurfaceReader::SelectObjects: Please select objects from listing: " << endl;

    for (const auto& [name, triangle_list] : objects) {
        int choice = 0;
        cout << "\nObject: " << name << " (yes=1, no=0): ";
        cin >> choice;

        if (choice == 1) {
            cout << "\nGoCad 'GEOLOGICAL_TYPE' options:" << endl;
            for (size_t i = 0; i < geo_types.size(); ++i) {
                cout << "  " << geo_types[i] << " (" << i << ")\n";
            }
            
            cout << "Selection: ";
            cin >> choice;

            // Validate input and store the pair
            if (choice >= 0 && static_cast<size_t>(choice) < geo_types.size()) {
                selections.push_back({name, geo_types[choice]});
            } else {
                cout << "Invalid selection. Skipping geological type assignment." << endl;
            }
        }
    }

    cout << "\nThank you." << endl;
    return selections.size();
}



/**
 
Creates a map of unique points in the object, labeled (keys) by their point
number. A map is also created which contains entries of the node numbers
which make up each triangle.   
*/
bool RhinoSurfaceReader::PopObject( const char *obj_name,
                                    map<size_t,Point<3> >&  points,
                                    map<size_t,vector<size_t> >& plist,
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
void RhinoSurfaceReader::WriteGocadHeader( const char* surf_name,
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




void RhinoSurfaceReader::OutputObjectTo(const char* obj, VSet<3U>& vset) const
{
    if (objects.empty()) {
        cerr << "\nRhinoSurfaceReader::OutputObjectTo: No objects are present.";
        return;
    }

    // 1. Find the object in the map
    string target_name(obj);
    auto it = objects.find(target_name);
    
    if (it == objects.end()) {
        cerr << "\nRhinoSurfaceReader::OutputObjectTo: Desired object '" << target_name << "' not found." << endl;
        return;
    }

    // 2. Generate unique points and connectivity (plist)
    // Note: I'm assuming PopObject was updated to your new Triangle type
    map<size_t, Point<3>> points;
    map<size_t, vector<size_t>> plist;
    PopObject(obj, points, plist);

    // 3. Topology calculation: Find neighbors and boundary flags
    vector<std::int8_t> pbflags;
    map<size_t, vector<int64_t>> pfverts;
    CreateNeighborPData(plist, pfverts, pbflags);

    // 4. Initialize VSet
    // Using a temporary to get element metadata (cleaner than hardcoding)
    LinearTriangle3D tri3;
    vset.Erase();
    vset.Resize(tri3.Nodes(),
                tri3.Neighbors(),
                tri3.ElementType(),
                points.size(), 
                plist.size());

    vset.SingleElementType(LINEAR_TRIANGLE3D);

    // 5. Transfer Point Coordinates
    for (const auto& [id, point] : points) {
        vset.Px(id, point[0]);
        vset.Py(id, point[1]);
        vset.Pz(id, point[2]);
    }

    // 6. Bulk transfer of boundary and connectivity data
    vset.AddBFlags(pbflags.begin(), pbflags.end());
    vset.AddPlist(plist.begin(), plist.end());
    vset.AddPfverts(pfverts.begin(), pfverts.end());

} // end OutputObjectTo








/** Writes all objects from the Rhino output file into a single Gocad file.


@section application Application

Not necessarily recommended since GoCad often has difficulties in reading
individual surfaces.  
 */
void RhinoSurfaceReader::WriteObjectsToTSurf(const char* tsurf_file,
                                             const char* GEOLOGICAL_TYPE) const
{
    // 1. Safety and File Path
    if (objects.empty()) {
        cout << "\nRhinoSurfaceReader::WriteObjectsToTSurf: No objects are present." << endl;
        return;
    }

    std::string filename = std::string(tsurf_file) + ".ts";
    ofstream ofs(filename, ios::out | ios::trunc);

    if (!ofs) {
        cout << "\nRhinoSurfaceReader::WriteObjectsToTSurf: ";
        cout << "Output file '" << filename << "' could not be opened" << endl;
        return;
    }

    // 2. Gocad Header
    // Note: Assuming WriteGocadHeader takes std::ostream& or ofstream&
    WriteGocadHeader(tsurf_file, GEOLOGICAL_TYPE, ofs);

    // 3. Loop through the objects using structured bindings
    // This replaces the complex map<string, list<mjl::Triangle3D>> iterator
    for (const auto& [name, triangles] : objects) {
        cout << "\nSaving object: " << name << " to file." << endl;
        
        // Assuming WriteObjectToTSurf handles the mjl-free triangles now
        WriteObjectToTSurf(name.c_str(), ofs);
    }

    // 4. Terminating GoCad TSurf format
    ofs << "END" << endl;

    ofs.close();
    cout << "\nRhinoSurfaceReader::WriteObjectsToTSurf: " << filename << " written successfully" << endl;
}
        




/**
 
Writes a series of files which contain one surface each corresponding to
the selected surfaces from the Rhino '.raw' dataset. The names of the 
objects are appended to the output file names which start with the
supplied file name. A hyphen is put between the filename and the 
object name.  

*/
void  RhinoSurfaceReader::WriteSelectedObjectsToTSurf( const char* tsurf_file ) const
 {
   if ( objects.empty() )
     {
        cout <<"\nRhinoSurfaceReader::WriteSelectedObjectsToTSurf: No objects are present.";
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
             cout <<"\nRhinoSurfaceReader::WriteSelectedObjectsToTSurf: ";
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
    
  cout <<"\nRhinoSurfaceReader::WriteSelectedObjectsToTSurf: '"<< tsurf_file;
  cout <<"' file series written successfully"<< endl;
  
 } // end WriteSelectedObjectsToTSurf  
               
               
              
              
                     

/** The surface will be written to file without associating any properties
with it.  

*/
void RhinoSurfaceReader::WriteObjectToTSurf( const char* obj, ofstream& ofs ) const
 {
    assert( ofs.is_open() );
 
    map<size_t,Point<3>>     points;
    map<size_t,vector<size_t> >  plist;

    // find object and create pxyz and plist arrays for the desired object
    PopObject( obj, points, plist, 0 );

    // write object to TSurf format into supplied file stream
    ofs <<"TFACE" << endl;
    
    // write vertex coordinates
    for ( auto it=points.begin(); it!=points.end(); it++ )
      {
         // identifier     node ID 1...n
         ofs <<"VRTX "<< (*it).first+1 <<" ";
         ofs << (*it).second[0] <<" "<< (*it).second[1] <<" "<< (*it).second[2] << endl;
      }

    // writing node id's per triangle as in plist
    for ( auto pit=plist.begin(); pit!=plist.end(); pit++ )
      {
         ofs <<"TRGL ";
         // again node ids must be augmented by 1 since Gocad counts 1...n
         for ( unsigned int i{0U}; i<3; i++ ) ofs << (*pit).second[i]+1 <<" ";
         ofs << endl;
      } 

    // terminating object description
    // ofs <<"END" << endl;

} // end WriteObjectToTSurf
     




size_t RhinoSurfaceReader::Objects() const
  { return objects.size(); }
  




void RhinoSurfaceReader::EraseObjects()
  { objects.erase( objects.begin(), objects.end() ); }            




void RhinoSurfaceReader::ExchangeCoordinateAxes(int axis_a, int axis_b)
{
    if (objects.empty()) return;
    if (axis_a < 0 || axis_a > 2 || axis_b < 0 || axis_b > 2) return;
    if (axis_a == axis_b) return;

    // A single swap of two axes is a reflection. 
    // Reflection flips the normal, so we MUST flip the winding order 
    // (swap nodes 1 and 2) to keep the normal pointing the same way.
    for (auto& [name, triangles] : objects) {
        for (auto& tri : triangles) {
            
            // 1. Swap the coordinates for all 3 vertices
            for ( uint32_t i = 0; i < 3; ++i) {
                Point<3> vertex = tri[i];
                std::swap(vertex[axis_a], vertex[axis_b]);
                tri.SetVertex(i, vertex); 
            }

            // 2. Correct the Winding Order
            // If the original was Node 0 -> 1 -> 2 (CCW)
            // After coordinate swap, it becomes CW.
            // Swapping Node 1 and Node 2 restores CCW.
            tri.ReverseWinding(); 
        }
    }
    
    cout << "Exchanged axes " << axis_a << " and " << axis_b << " and updated winding orders." << endl;
}
 
 
     
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
void RhinoSurfaceReader::ScaleCoordinates(double x_fac, double y_fac, double z_fac)
{
    if (objects.empty()) {
        cout << "\nRhinoSurfaceReader::ScaleCoordinates: No object data present. Nothing done..." << endl;
        return;
    }

    // Optimization: If all factors are 1.0, just exit
    if (x_fac == 1.0 && y_fac == 1.0 && z_fac == 1.0) return;

    // Use structured bindings (C++17) to iterate over name and triangle list
    for (auto& [name, triangles] : objects) {
        for (auto& tri : triangles) {
            Point<3> pts[3];
            
            for (uint32_t i = 0; i < 3; ++i) {
                // 1. Get the current vertex
                pts[i] = tri[i];

                // 2. Scale the coordinates
                pts[i][0] *= x_fac;
                pts[i][1] *= y_fac;
                pts[i][2] *= z_fac;
            }

            // 3. Update the triangle with the scaled points
            // Assuming tri.id_ is accessible or stored internally
            tri.Set(pts[0], pts[1], pts[2], tri.id_);
            
            // Note: If any scale factor is negative, it will flip the normal!
            // In a serious CAD/FEA context, you might check if (x_fac*y_fac*z_fac < 0)
            // and reverse the winding order if it does.
            if ((x_fac * y_fac * z_fac) < 0.0) {
                tri.ReverseWinding();
            }
        }
    }
    
    cout << "Scaled coordinates by factors: (" << x_fac << ", " << y_fac << ", " << z_fac << ")" << endl;
}
 
 
 
     
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
void RhinoSurfaceReader::MoveCoordinates( double x_move, double y_move, double z_move )
 {
    if ( objects.empty() )
      {
         cout <<"\nRhinoSurfaceReader::MoveCoordinates: No object data are present. ";
         cout <<"Nothing was done..."<< endl;
         return;
      }

    // looping over the objects 
    Point<3>  pt[3];
    
    for ( auto oit = objects.begin(); oit!=objects.end(); oit++ )
      for ( auto tit=(*oit).second.begin(); tit!=(*oit).second.end(); tit++ )
        {
           for ( int32_t i=0; i<3; i++ )
             {
                // changing point by point 
                // 1. getting the point 
                pt[i] = (*tit)[i];
                // 2. moving the point coordinates
                pt[i][0] += x_move;
                pt[i][1] += y_move;
                pt[i][2] += z_move;
             }
           // storing the reworked points re-organising the triangle if necessary
           auto pid = (*tit).id_;
           (*tit).Set( pt[0], pt[1], pt[2], pid );
       }

 } // end MoveCoordinates

} // end namespace csmp




















