#include "Tecplot_Interface.h"
#include "InterFace.h"
#include "Region.h"
#include "Model.h"
#include "CSMP_highLevelUtilities.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim>
void Tecplot_Interface<dim>
::OutputDataToTecplotFile( const Model<dim>&  sg,
                           const char*        file_name,
                           const char*        var_name,
                           long               timestep )
  {
    this->OutputDataToTecplotFile(sg, sg.Region("Model"), file_name, var_name, timestep );
  } // end OutputDataToTecplotFile


template<size_t dim>
void Tecplot_Interface<dim>
::OutputDataToTecplotFile( const Model<dim>&  sg,
                           const char*        group_name,
                           const char*        file_name,
                           const char*        var_name,
                           long               timestep )
  {
    this->OutputDataToTecplotFile(sg, sg.Region(group_name), file_name, var_name, timestep );
  }
    
template<size_t dim>
void Tecplot_Interface<dim>
::OutputGridToTecplotFile( const Model<dim>&  sg,
                           const char*        file_name,
                           long               timestep )
{
    this->OutputGridToTecplotFile(sg, sg.Region("Model"), file_name, timestep );
} // end OutputDataToTecplotFile


template<size_t dim>
void Tecplot_Interface<dim>
::OutputGridToTecplotFile( const Model<dim>&  sg,
                           const char*        group_name,
                           const char*        file_name,
                           long               timestep )
{
    this->OutputGridToTecplotFile(sg, sg.Region(group_name), file_name, timestep );
}
    
template<size_t dim>
void Tecplot_Interface<dim>
::OutputDataToTecplotFile( const Model<dim>&  sg,
                           const Region<dim>& gref,
                           const char*        file_name,
                           const char*        var_name,
                           long               timestep )
{
     csmp::Index  prop_key = sg.Database().StorageKey(var_name);
     string  variable(var_name);
     replaceWhiteSpaceBy( variable, '-' );

     if ( prop_key.place != NODE && prop_key.place != ELEMENT ) {
         cout << "\nTecplot_Interface<dim>::OutputDataToTecplotFile: Only NODE and ELEMENT data can be visualized. Nothing is done..." << endl;
         return;
       }

     if ( prop_key.type != SCALAR && prop_key.type != VECTOR ) {
         cout << "\nTecplot_Interface<dim>::OutputDataToTecplotFile: Only SCALAR and VECTOR variables can be visualized. Nothing is done..." << endl;
         return;
       }

     if ( prop_key.place == ELEMENT && prop_key.type != VECTOR ) {
         cout << "\nTecplot_Interface<dim>::OutputDataToTecplotFile: Only VECTOR data can be visualized for ELEMENT. Nothing is done..." << endl;
         return;
       }
       
     if ( sg.Mesh().HybridElementMesh() ) {
         cout << "\nTecplot_Interface<dim>::OutputDataToTecplotFile: Tecplot cannot visualize mixed element meshes Nothing is done..." << endl;
         return;
       }

     const Region<dim>&  super_group(sg.Region("Model"));

     if ( super_group.E(0)->FE_Type() != LINEAR_TRIANGLE &&
          super_group.E(0)->FE_Type() != LINEAR_TETRAHEDRON &&
          super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_TRIANGLE && 
          super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_TETRAHEDRON && 
          super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_QUADRILATERAL &&
          super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_HEXAHEDRON ) {
          cout << "\nTecplot_Interface<dim>::OutputDataToTecplotFile: This element type cannot be visualized. Nothing is done..." << endl;
          return;
       }

  
     // 0. creating list of element ID's of associated Region
     // --------------------------------------------------------------------
     vector<size_t>     elmt_ids;
     gref.MemberElementIndexes( elmt_ids );
                

     // 1. getting new node mapping and updating storage if geometry has changed
     // ------------------------------------------------------------------------
     if ( prop_key.place == NODE ) {
         NodeBasedTopology( sg, elmt_ids, plist, node_mapping );
         TransformPlist( sg, plist, transformed_plist );
       }
     else {
         PointBasedTopology( sg, elmt_ids, plist, node_mapping );
         transformed_plist = plist;
       }
     if ( prop_key.place == NODE ) NodeData( sg, prop_key, node_mapping, pxyz_data );
     else                          ElementPointData( sg, prop_key, node_mapping, pxyz_data );

     // 2. opening data output file in ascii format
     // -------------------------------------------
     string  region_name( gref.Name() );
     string  outfile( region_name.c_str() );
     outfile += "-";
     outfile += file_name;
     replaceWhiteSpaceBy( outfile, '-' );
     char  num[30];
     sprintf( num, "%ld", timestep );
     outfile += num;
     outfile += ".dat";
       
     ofstream ofs;
     ofs.open( outfile.c_str(), ios::out|ios::trunc );
     if ( !ofs )
       {
          cout <<"\nTecplot_Interface<dim>::OutputDataToTecplotFile: outfile cannot be opened."<< endl;
          return;
       }
      
     // 2.1 write file header
     ofs <<"TITLE = \""<< var_name <<" CSMP simulation of region " << region_name << " output at timestep "<< timestep <<"\""<< endl;
     if ( prop_key.place == NODE ) {
         if ( prop_key.type == SCALAR )
           ofs <<"VARIABLES = x,\t y,\t z,\t " << variable << ",\t node" << endl;
         else 
           ofs <<"VARIABLES = x,\t y,\t z,\t " << variable << "1,\t " << variable << "2,\t " << variable << "3,\t node" << endl;
         }
     else ofs <<"VARIABLES = x,\t y,\t z,\t " << variable << "1,\t " << variable << "2,\t " << variable << "3,\t element" << endl; 
     
     if ( prop_key.place == NODE ) {   
         if ( gref.E(0)->FE_Type() == LINEAR_TRIANGLE || 
              gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_TRIANGLE )
           ofs <<"ZONE T=\""<<region_name<<" mesh\" n="<< gref.Nodes() <<", e= "<< gref.Elements() <<", et = triangle";
         if ( gref.E(0)->FE_Type() == LINEAR_TETRAHEDRON ||
              gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_TETRAHEDRON  )
           ofs <<"ZONE T=\""<<region_name<<" mesh\" n="<< gref.Nodes() <<", e= "<< gref.Elements() <<", et = tetrahedron";
         if ( gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL )
           ofs <<"ZONE T=\""<<region_name<<" mesh\" n="<< gref.Nodes() <<", e= "<< gref.Elements() <<", et = quadrilateral";
         if ( gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_HEXAHEDRON )
           ofs <<"ZONE T=\""<<region_name<<" mesh\" n="<< gref.Nodes() <<", e= "<< gref.Elements() <<", et = brick";
         ofs <<", f = fepoint"<< endl;
       }

     ofs.setf( ios::scientific );
     
     // 3. writing node coordinates and values
     // --------------------------------------
     typename map<size_t,vector<double64> >::const_iterator  nit;
     size_t counter(1);
     for ( nit=pxyz_data.begin(); nit!=pxyz_data.end(); nit++ )
       {
          for ( size_t i=0; i<dim; i++ ) ofs << (*nit).second[i] <<" ";
          if ( dim == 2U ) ofs << 0.0 <<" ";
          // SCALAR
          if ( prop_key.type == SCALAR ) ofs << (*nit).second[3] <<" ";
          // VECTOR
          else {
              for ( size_t j=3; j<(*nit).second.size(); j++ )  ofs << (*nit).second[j] <<" ";
            }
          ofs << counter;
          ofs << endl;
          counter++;
       }
     ofs << endl;  
            
     // 4. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     if ( prop_key.place == NODE ) {
         typename map<size_t,vector<size_t> >::const_iterator  eit;
         typename vector<size_t>::const_iterator  it;

         for ( eit=transformed_plist.begin(); eit!=transformed_plist.end(); eit++ )
           {
              for ( it=(*eit).second.begin(); it!=(*eit).second.end(); it++ ) ofs << *it+1 <<" ";
              ofs << endl;
           }   
         ofs << endl;
       }
     
     ofs.close();
     cout <<"\nTecplot_Interface<"<< dim <<">::OutputDataToTecplotFile: file '";
     cout << outfile <<"' written successfully."<< endl;

  } // end OutputDataToTecplotFile

template<size_t dim>
void Tecplot_Interface<dim>
::OutputGridToTecplotFile( const Model<dim>&  sg,
                           const Region<dim>& gref,
                           const char*        file_name,
                           long               timestep )
{
    const Region<dim>&  super_group(sg.Region("Model"));
    
    if ( super_group.E(0)->FE_Type() != LINEAR_TRIANGLE &&
         super_group.E(0)->FE_Type() != LINEAR_TETRAHEDRON &&
         super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_TRIANGLE &&
         super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_TETRAHEDRON &&
         super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_QUADRILATERAL &&
         super_group.E(0)->FE_Type() != ISOPARAMETRIC_LINEAR_HEXAHEDRON ) {
        cout << "\nTecplot_Interface<dim>::OutputDataToTecplotFile: This element type cannot be visualized. Nothing is done..." << endl;
        return;
    }
    
    // 0. creating list of element ID's of associated Region
    // --------------------------------------------------------------------
    vector<size_t>     elmt_ids;
    gref.MemberElementIndexes( elmt_ids );
    
    // 1. getting new node mapping and updating storage if geometry has changed
    // ------------------------------------------------------------------------
    NodeBasedTopology( sg, elmt_ids, plist, node_mapping );
    //TransformPlist( sg, plist, transformed_plist );
    XyzData( sg, node_mapping, pxyz_data );
    
    // 2. opening data output file in ascii format
    // -------------------------------------------
    string  region_name( gref.Name() );
    string  outfile( region_name.c_str() );
    outfile += "-";
    outfile += file_name;
    replaceWhiteSpaceBy( outfile, '-' );
    char  num[30];
    sprintf( num, "%ld", timestep );
    outfile += num;
    outfile += ".dat";
    
    ofstream ofs;
    ofs.open( outfile.c_str(), ios::out|ios::trunc );
    if ( !ofs )
    {
        cout <<"\nTecplot_Interface<dim>::OutputDataToTecplotFile: outfile cannot be opened."<< endl;
        return;
    }
    
    // 2.1 write file header
    ofs <<"TITLE = \"CSMP mesh of region " << region_name << " output at timestep "<< timestep <<"\""<< endl;
    ofs <<"FILETYPE = GRID"<< endl;
    ofs <<"VARIABLES = \"x\", \"y\", \"z\"" << endl;
    
    ofs << "ZONE T = \""<<region_name<<"\", ";
    ofs << "N = " << gref.Nodes() << ", E = " << gref.Elements() << ", ";
    ofs << "DATAPACKING = POINT, ";
    if ( gref.E(0)->FE_Type() == LINEAR_TRIANGLE || gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_TRIANGLE )
        ofs << "ZONETYPE = FETRIANGLE,";
    else if ( gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL )
        ofs << "ZONETYPE = FEQUADRILATERAL,";
    else if ( gref.E(0)->FE_Type() == LINEAR_TETRAHEDRON || gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_TETRAHEDRON  )
        ofs << "ZONETYPE = FETETRAHEDRON,";
    else if ( gref.E(0)->FE_Type() == ISOPARAMETRIC_LINEAR_HEXAHEDRON )
        ofs << "ZONETYPE = FEBRICK";
    else if ( dim == 3 )
        ofs << "ZONETYPE = FEPOLYHEDRON,";
    else
        ofs << "ZONETYPE = FEPOLYGON";
    ofs << endl;
    
    ofs.setf( ios::scientific );
    
    
    // 3. writing node coordinates
    // --------------------------------------
    typename map<size_t,vector<double64> >::const_iterator  nit;
    typename map<size_t,vector<double64> >::const_iterator  nitEnd = pxyz_data.end();
    for ( nit=pxyz_data.begin(); nit != nitEnd; nit++ )
    {
        for ( size_t i=0; i < 3U; i++ ){
            ofs << (*nit).second[i];
        if( i != 2U )
            ofs <<" ";
        }
        ofs << endl;
    }
    ofs << endl;
    

    // 4. writing element nodes connection list
    // --------------------------------------
    typename map<size_t,vector<size_t> >::const_iterator  elit;
    typename map<size_t,vector<size_t> >::const_iterator  elitEnd = plist.end();
    typename vector<size_t>::const_iterator  elnit;
    typename vector<size_t>::const_iterator  elnitEnd;
    for ( elit=plist.begin(); elit != elitEnd; elit++ )
    {
        elnitEnd = (*elit).second.end();
        for ( elnit = (*elit).second.begin(); elnit != elnitEnd; )
        {
            ofs << (*elnit) + 1;
            if( ++elnit != elnitEnd )
                ofs <<" ";
        }
        ofs << endl;
    }
    ofs << endl;

    
    /*
    // 5. writing face neighbor connections list
    // -----------------------------------------------------
    typename map<size_t,vector<size_t> >::const_iterator  eit;
    typename vector<size_t>::const_iterator  it;
    for ( eit=transformed_plist.begin(); eit!=transformed_plist.end(); eit++ )
    {
        for ( it=(*eit).second.begin(); it!=(*eit).second.end(); it++ )
            ofs << (*it) + 1 <<" ";
        ofs << endl;
    }
    ofs << endl;
    
    ofs.close();
    cout <<"\nTecplot_Interface<"<< dim <<">::OutputDataToTecplotFile: file '";
    cout << outfile <<"' written successfully."<< endl;
    */
    
} // end OutputDataToTecplotFile



template<size_t dim>
void Tecplot_Interface<dim>::NodeBasedTopology( const Model<dim>& sg,
                                                const vector<size_t>& elmt_ids,
                                                map<size_t,vector<size_t> >& plist,
                                                map<size_t,size_t >& node_nums )
 {
    vector<size_t>  ids(3);
    vector<size_t>  pentry(3);
    
    plist.erase( plist.begin(), plist.end() );
    node_nums.erase( node_nums.begin(), node_nums.end() );
    
    const Region<dim>&  super_group(sg.Region("Model"));

    // 1. creating unique node number list, and plist by looping over the selected elements
    // ------------------------------------------------------------------------------------
    super_group.UpdateMemberIndexes();
    size_t  nodes(0U);
    
    for ( typename vector<size_t>::const_iterator
          lit=elmt_ids.begin(); lit!=elmt_ids.end(); lit++ )
      {
         const size_t eid( *lit );
         // getting node ids and renumbering them 0...n-1
         for ( size_t i=0U; i<super_group.E( eid )->Nodes(); i++ ) {
              pair<typename map<size_t,size_t>::iterator,bool>
                node_it = node_nums.insert( make_pair( super_group.E( eid )->N(i)->Idx(), nodes ) );
              if ( node_it.second == true ) nodes++;
           }
         
         // building the plist, minimizing the search by always using the smallest 
         // size of the map possible
         pentry.resize( super_group.E( eid )->Nodes() );
         for ( size_t i=0; i<super_group.E( eid )->Nodes(); i++ )
           {
              typename map<size_t,size_t>::const_iterator
                nit = node_nums.find( super_group.E( eid )->N(i)->Idx() );
              pentry[i] = (*nit).second;
           }
         
         // inserting new element id and empty vector into the plist
         pair<typename map<size_t,vector<size_t> >::iterator,bool>
           pit = plist.insert( make_pair( eid, vector<size_t>() ) );
         assert( pit.second == true );
         // initializing plist vector
         (*pit.first).second.reserve( super_group.E( eid )->Nodes() );
         for ( size_t i=0; i<super_group.E( eid )->Nodes(); i++ )
           (*pit.first).second.push_back( pentry[i] );
      }   
     
 } // end NodeBasedTopology



template<size_t dim>
void Tecplot_Interface<dim>::PointBasedTopology( const Model<dim>& sg,
                                                 const vector<size_t>& elmt_ids,
                                                 map<size_t,vector<size_t> >& plist,
                                                 map<size_t,size_t>& node_nums )
{
    size_t  nodes(0U);
    vector<size_t>  pentry(1);
    pair<size_t,vector<size_t> >  alpha_help;
    vector<size_t>                alpha_vec(1);

    const Region<dim>&  super_group(sg.Region("Model"));

    plist.erase( plist.begin(), plist.end() );
    node_nums.erase( node_nums.begin(), node_nums.end() );

    // 1. creating unique point number list, and plist by looping over the selected elements
    // ------------------------------------------------------------------------------------
    for ( typename vector<size_t>::const_iterator
         lit=elmt_ids.begin(); lit!=elmt_ids.end(); lit++ )
    {
        const size_t eid( *lit );
        // getting node ids and renumbering them including duplicates 0...elmts * npe's
        node_nums[ nodes ] = super_group.E( eid )->Idx();
        pentry[0] = nodes++;
        // inserting new element id and single-element vector into plist
        alpha_help.first  = eid;
        alpha_vec[0]      = pentry[0];
        alpha_help.second = alpha_vec;
        pair<typename map<size_t,vector<size_t> >::iterator,bool> pit = plist.insert( alpha_help );
        assert( pit.second == true );
    }

} // end PointBasedTopology




template<size_t dim>
void Tecplot_Interface<dim>::TransformPlist( const Model<dim>&  sg,
                                             const map<size_t,vector<size_t> >& plist,
                                             map<size_t,vector<size_t> >& tplist )
 {
    const Region<dim>&  super_group(sg.Region("Model"));
 
    map<size_t,vector<size_t> >::const_iterator  it;
    // triangle   tetrahedron
    vector<size_t> pentry(3), tentry(4);
    size_t         entries(1);
 
    tplist.erase( tplist.begin(), tplist.end() );
   
    for ( it=plist.begin(); it!=plist.end(); it++ )
      {
         switch ( super_group.E( (*it).first )->Nodes() )
           {
              case 2: // segments (vectors are just copied over)
                   tplist[ (*it).first ] = (*it).second;
                   break;
               
              case 3: // triangle (vectors are just copied over)
                   tplist[ (*it).first ] = (*it).second;
                break;
                
              case 4: // tetrahedron (vectors are just copied over)
                   tplist[ (*it).first ] = (*it).second;
                break;
              
              case 6: // quadratic triangle -> 4 triangles
                   // triangle 1
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[3];
                   pentry[2] = (*it).second[5];
                   tplist[ entries++ ] = pentry;
                   // triangle 2
                   pentry[0] = (*it).second[3];
                   pentry[1] = (*it).second[1];
                   pentry[2] = (*it).second[4];
                   tplist[ entries++ ] = pentry;
                   // triangle 3
                   pentry[0] = (*it).second[5];
                   pentry[1] = (*it).second[4];
                   pentry[2] = (*it).second[2];
                   tplist[ entries++ ] = pentry;
                   // triangle 4
                   pentry[0] = (*it).second[3];
                   pentry[1] = (*it).second[4];
                   pentry[2] = (*it).second[5];
                   tplist[ entries++ ] = pentry;
                break;

              case 7: // quadratic barycentric triangle -> 6 triangles
                   // triangle 1
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[3];
                   pentry[2] = (*it).second[6];
                   tplist[ entries++ ] = pentry;
                   // triangle 2
                   pentry[0] = (*it).second[3];
                   pentry[1] = (*it).second[1];
                   pentry[2] = (*it).second[6];
                   tplist[ entries++ ] = pentry;
                   // triangle 3
                   pentry[0] = (*it).second[6];
                   pentry[1] = (*it).second[1];
                   pentry[2] = (*it).second[4];
                   tplist[ entries++ ] = pentry;
                   // triangle 4
                   pentry[0] = (*it).second[6];
                   pentry[1] = (*it).second[4];
                   pentry[2] = (*it).second[2];
                   tplist[ entries++ ] = pentry;
                   // triangle 5
                   pentry[0] = (*it).second[5];
                   pentry[1] = (*it).second[6];
                   pentry[2] = (*it).second[2];
                   tplist[ entries++ ] = pentry;
                   // triangle 6
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[6];
                   pentry[2] = (*it).second[5];
                   tplist[ entries++ ] = pentry;
                break;

/**********************************************************************************************************
              case 8: // linear hexahedron (fix according to ANSYSs description: F1 to F6)
                   // quadrilateral 1 - Face 1 Bottom
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[1];
                   pentry[2] = (*it).second[2];
                   tplist[ entries++ ] = pentry;
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[2];
                   pentry[2] = (*it).second[3];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 2 - Face 2 Front
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[1];
                   pentry[2] = (*it).second[5];
                   tplist[ entries++ ] = pentry;
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[5];
                   pentry[2] = (*it).second[4];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 3 - Face 3 Right
                   pentry[0] = (*it).second[1];
                   pentry[1] = (*it).second[2];
                   pentry[2] = (*it).second[6];
                   tplist[ entries++ ] = pentry;
                   pentry[0] = (*it).second[1];
                   pentry[1] = (*it).second[6];
                   pentry[2] = (*it).second[5];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 4 - Face 4 Back
                   pentry[0] = (*it).second[3];
                   pentry[1] = (*it).second[2];
                   pentry[2] = (*it).second[6];
                   tplist[ entries++ ] = pentry;
                   pentry[0] = (*it).second[3];
                   pentry[1] = (*it).second[6];
                   pentry[2] = (*it).second[7];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 5 - Face Left
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[3];
                   pentry[2] = (*it).second[7];
                   tplist[ entries++ ] = pentry;
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[7];
                   pentry[2] = (*it).second[4];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 6 - Face Up
                   pentry[0] = (*it).second[4];
                   pentry[1] = (*it).second[5];
                   pentry[2] = (*it).second[6];
                   tplist[ entries++ ] = pentry;
                   pentry[0] = (*it).second[4];
                   pentry[1] = (*it).second[6];
                   pentry[2] = (*it).second[7];
                   tplist[ entries++ ] = pentry;
                break;
***********************************************************************************************************/

                   // Initial version crashes???
//**********************************************************************************************************
              case 8: // linear hexahedron
                   // quadrilateral 1
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[1];
                   pentry[2] = (*it).second[2];
                   pentry[3] = (*it).second[3];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 2
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[1];
                   pentry[2] = (*it).second[5];
                   pentry[3] = (*it).second[4];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 3
                   pentry[0] = (*it).second[1];
                   pentry[1] = (*it).second[2];
                   pentry[2] = (*it).second[6];
                   pentry[3] = (*it).second[5];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 4
                   pentry[0] = (*it).second[3];
                   pentry[1] = (*it).second[2];
                   pentry[2] = (*it).second[6];
                   pentry[3] = (*it).second[7];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 5
                   pentry[0] = (*it).second[0];
                   pentry[1] = (*it).second[3];
                   pentry[2] = (*it).second[7];
                   pentry[3] = (*it).second[4];
                   tplist[ entries++ ] = pentry;
                   // quadrilateral 6
                   pentry[0] = (*it).second[4];
                   pentry[1] = (*it).second[5];
                   pentry[2] = (*it).second[6];
                   pentry[3] = (*it).second[7];
                   tplist[ entries++ ] = pentry;
                   break;
//***********************************************************************************************************/
               
              case 10: // quadratic tetrahedron -> 13 tetrahedra
                   // 4 corner tetrahedra
                   // tetrahedron 1
                   tentry[0]=(*it).second[1];
                   tentry[1]=(*it).second[5];
                   tentry[2]=(*it).second[4];
                   tentry[3]=(*it).second[7];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 2
                   tentry[0]=(*it).second[2];
                   tentry[1]=(*it).second[6];
                   tentry[2]=(*it).second[5];
                   tentry[3]=(*it).second[8];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 3
                   tentry[0]=(*it).second[0];
                   tentry[1]=(*it).second[4];
                   tentry[2]=(*it).second[6];
                   tentry[3]=(*it).second[9];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 4
                   tentry[0]=(*it).second[7];
                   tentry[1]=(*it).second[8];
                   tentry[2]=(*it).second[9];
                   tentry[3]=(*it).second[3];
                   tplist[ entries++ ] = tentry;
                   // 4 tetrahedra which make up central octahedron
                   // tetrahedron 5
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[5];
                   tentry[2]=(*it).second[8];
                   tentry[3]=(*it).second[7];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 6
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[5];
                   tentry[2]=(*it).second[6];
                   tentry[3]=(*it).second[8];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 7
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[6];
                   tentry[2]=(*it).second[9];
                   tentry[3]=(*it).second[8];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 8
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[8];
                   tentry[2]=(*it).second[9];
                   tentry[3]=(*it).second[7];
                   tplist[ entries++ ] = tentry;
                break;
                
              case 11: // barycentric quadratic tetrahedron -> 13 tetrahedra
                   // basal 3 outer tetrahedra
                   // tetrahedron 1
                   tentry[0]=(*it).second[0];
                   tentry[1]=(*it).second[4];
                   tentry[2]=(*it).second[6];
                   tentry[3]=(*it).second[7];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 2
                   tentry[0]=(*it).second[1];
                   tentry[1]=(*it).second[5];
                   tentry[2]=(*it).second[4];
                   tentry[3]=(*it).second[8];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 3
                   tentry[0]=(*it).second[2];
                   tentry[1]=(*it).second[6];
                   tentry[2]=(*it).second[5];
                   tentry[3]=(*it).second[9];
                   tplist[ entries++ ] = tentry;
                   // top 3 tetrahedra
                   // tetrahedron 4
                   tentry[0]=(*it).second[7];
                   tentry[1]=(*it).second[8];
                   tentry[2]=(*it).second[10];
                   tentry[3]=(*it).second[3];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 5
                   tentry[0]=(*it).second[8];
                   tentry[1]=(*it).second[9];
                   tentry[2]=(*it).second[10];
                   tentry[3]=(*it).second[3];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 6
                   tentry[0]=(*it).second[7];
                   tentry[1]=(*it).second[10];
                   tentry[2]=(*it).second[9];
                   tentry[3]=(*it).second[3];
                   tplist[ entries++ ] = tentry;
                   // intermediate 3 side tetrahedra
                   // tetrahedron 7
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[7];
                   tentry[2]=(*it).second[8];
                   tentry[3]=(*it).second[10];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 8
                   tentry[0]=(*it).second[5];
                   tentry[1]=(*it).second[8];
                   tentry[2]=(*it).second[9];
                   tentry[3]=(*it).second[10];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 9
                   tentry[0]=(*it).second[6];
                   tentry[1]=(*it).second[7];
                   tentry[2]=(*it).second[10];
                   tentry[3]=(*it).second[9];
                   tplist[ entries++ ] = tentry;
                   // oblique 3 tetrahedra in the lower part of the parent tetrahedron
                   // tetrahedron 10
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[6];
                   tentry[2]=(*it).second[7];
                   tentry[3]=(*it).second[10];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 11
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[5];
                   tentry[2]=(*it).second[10];
                   tentry[3]=(*it).second[8];
                   tplist[ entries++ ] = tentry;
                   // tetrahedron 12
                   tentry[0]=(*it).second[6];
                   tentry[1]=(*it).second[5];
                   tentry[2]=(*it).second[9];
                   tentry[3]=(*it).second[10];
                   tplist[ entries++ ] = tentry;
                   // lower plane basal tetrahedron
                   // tetrahedron 13
                   tentry[0]=(*it).second[4];
                   tentry[1]=(*it).second[5];
                   tentry[2]=(*it).second[6];
                   tentry[3]=(*it).second[10];
                   tplist[ entries++ ] = tentry;
                break;

              default:
                   cout <<"\nElement ID: "<< super_group.E( (*it).first )->Idx();
                   cout <<" with "<< super_group.E( (*it).first )->Nodes() <<" nodes."<< endl;
                   throw csmp::Exception( CSMP_FATAL_ERROR, "Tecplot_Interface<dim>::TransformPlist", 
                                  "Unable to interpret how this element shall be broken in subelements");
           }
      }
 
 } // end TransformPlist

template<size_t dim>
void Tecplot_Interface<dim>::XyzData( const Model<dim>& sg,
                                     const map<size_t,size_t>& node_nums,
                                     map<size_t,vector<double64> >& pxyz_data )
{
    const Region<dim>&  super_group(sg.Region("Model"));

    pair<typename map<size_t,vector<double64> >::iterator,bool>  dit;
    typename map<size_t,size_t>::const_iterator                  nit;

    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );

    // extracting the node coordinates and property values from the Model
    // results are put into the first three elements of 'pxyz_data'
    // --------------------------------------------------------------------------
    for ( nit=node_nums.begin(); nit!=node_nums.end(); nit++ )
    {
        // inserting new element into the cordinate map node id
        dit = pxyz_data.insert( make_pair( (*nit).second, vector<double64>() ) );
        assert( dit.second == true );
        
        // node coordinates (reserves storage for three coordinates )
        (*dit.first).second.reserve(3);
        (*dit.first).second.push_back( super_group.N( (*nit).first )->x() );
        if ( dim == 1 ) (*dit.first).second.push_back( 0.0 );
        else (*dit.first).second.push_back( super_group.N( (*nit).first )->y() );
        if ( dim == 3 ) (*dit.first).second.push_back( super_group.N( (*nit).first )->z() );
        else            (*dit.first).second.push_back( 0. );
    }
    
} // end XyzData

template<size_t dim>
void Tecplot_Interface<dim>::NodeData( const Model<dim>& sg,
                                      const csmp::Index& prop_key,
                                      const map<size_t,size_t>& node_nums,
                                      map<size_t,vector<double64> >& pxyz_data )
{
    const Region<dim>&  super_group(sg.Region("Model"));

    if ( prop_key.place != NODE )
        throw csmp::Exception( CSMP_FATAL_ERROR, "Tecplot_Interface<dim>::NodeData",
                              "Method only applies to node properties" );

    pair<typename map<size_t,vector<double64> >::iterator,bool>  dit;
    typename map<size_t,size_t>::const_iterator                nit;
    ScalarVariable      sc;
    VectorVariable<dim>  vc;
    TensorVariable<dim>  ts;

    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );

        // extracting the node coordinates and property values from the Model
        // results are put into the first three elements of 'pxz_data'
        // --------------------------------------------------------------------------
    for ( nit=node_nums.begin(); nit!=node_nums.end(); nit++ )
        {
            // inserting new element into the cordinate map using new node coordinate number
        dit = pxyz_data.insert( make_pair((*nit).second, vector<double64>() ) );
        assert( dit.second == true );
        
            // node coordinates (reserves storage for three coordinates and a scalar data value)
        (*dit.first).second.reserve(4);
        (*dit.first).second.push_back( super_group.N( (*nit).first )->x() );
        (*dit.first).second.push_back( super_group.N( (*nit).first )->y() );
        if ( dim == 3 ) (*dit.first).second.push_back( super_group.N( (*nit).first )->z() );
        else            (*dit.first).second.push_back( 0. );
        
            // data values
        if ( prop_key.type == SCALAR )
            {
            super_group.N( (*nit).first )->Read( prop_key, sc );
            (*dit.first).second.push_back( sc() );
            }
        else if ( prop_key.type == VECTOR )
            {
                // storage for 3 coordinate values and 3 data values (in 2D 3rd place = 0.0)
            (*dit.first).second.reserve(6);
            super_group.N( (*nit).first )->Read( prop_key, vc );
                // variables have always 3 components since view screen is 3D
            for ( size_t j=0; j<2; j++ ) (*dit.first).second.push_back( vc[j] );
            if ( dim == 3 )            (*dit.first).second.push_back( vc[2] );
            else                       (*dit.first).second.push_back( 0.0 );
            }
        else if ( prop_key.type == TENSOR )
            {
                // storage for 3 coordinate values and 9 data values (in 2D 3rd row, column = 0.0)
            (*dit.first).second.reserve(12);
            super_group.N( (*nit).first )->Read( prop_key, ts );
                // variables have always 3 components since view screen is 3D
            if ( dim == 3 )
                for ( size_t k=0; k<3; k++ )
                    for ( size_t l=0; l<3; l++ ) (*dit.first).second.push_back( ts(k,l) );
            else {
                (*dit.first).second.push_back( ts(0,0) );
                (*dit.first).second.push_back( ts(0,1) );
                (*dit.first).second.push_back( 0. );
                (*dit.first).second.push_back( ts(1,0) );
                (*dit.first).second.push_back( ts(1,1) );
                (*dit.first).second.push_back( 0. );
                (*dit.first).second.push_back( 0. );
                (*dit.first).second.push_back( 0. );
                (*dit.first).second.push_back( 0. );
            }
            }
        }

} // end NodeData


template<size_t dim>
void Tecplot_Interface<dim>::ElementPointData( const Model<dim>& sg,
                                               const csmp::Index&     prop_key,
                                               const map<size_t,size_t>& node_nums,
                                               map<size_t,vector<double64> >& pxyz_data )
 {
    const Region<dim>&  super_group(sg.Region("Model"));
    
    if ( prop_key.place != ELEMENT )
      throw csmp::Exception( CSMP_FATAL_ERROR, "Tecplot_Interface<dim>::ElementPointData", 
                                   "Method only applies to element properties" );
    
    pair<typename map<size_t,vector<double64> >::iterator,bool>  dit;
    typename map<size_t,size_t>::const_iterator               nit;
    ScalarVariable      sc;
    VectorVariable<dim>  vc;
    TensorVariable<dim>  ts;
    VectorVariable<dim>  pt;
    
    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );
    
    // 1. extracting the element-center coordinates and property values from the Model
    //    results are put into the first three elements of 'pxz_data'
    // --------------------------------------------------------------------------        
    for ( nit=node_nums.begin(); nit!=node_nums.end(); nit++ )
      {
         // inserting new element into the cordinate map using new node coordinate number
         dit = pxyz_data.insert( make_pair( (*nit).first, vector<double64>() ) );
         assert( dit.second == true );
         
         // node coordinates (reserves storage for three coordinates and a scalar data value)
         (*dit.first).second.reserve(4);
         pt = super_group.E( (*nit).second )->BaryCenter();
         (*dit.first).second.push_back( pt[0] );
         (*dit.first).second.push_back( pt[1] );
         if ( dim == 3 ) (*dit.first).second.push_back( pt[2] );
         else                 (*dit.first).second.push_back( 0. );
      }
      
    // 2. Looping through the plist, associating the element properties with the nodal pxyz_data     
    // data values
    // -----------------------------------------------------------------------------------------
    typename map<size_t,vector<size_t> >::const_iterator  it;
    typename vector<size_t>::const_iterator                  pit;
    typename map<size_t,vector<double64> >::iterator        dit2;

    for ( it=plist.begin(); it!=plist.end(); it++ )
      {
         if ( prop_key.type == SCALAR )
           {
              super_group.E( (*it).first )->Read( prop_key, sc );
              // looping over the element's node 
              for ( pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   dit2 = pxyz_data.find( (*pit) );
                   (*dit2).second.push_back( sc() );
                }
           }
         else if ( prop_key.type == VECTOR )
           {
              super_group.E( (*it).first )->Read( prop_key, vc );
              for ( pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   dit2 = pxyz_data.find( (*pit) );
                  (*dit2).second.reserve(6);
                   if ( dim == 3U )
                     for ( size_t j=0; j<3; j++ ) (*dit2).second.push_back( vc[j] );
                   else {
                        (*dit2).second.push_back( vc[0] );
                        (*dit2).second.push_back( vc[1] );
                        (*dit2).second.push_back( 0.0 );
                     }
                }
             
           }
         else if ( prop_key.type == TENSOR )
           {
              super_group.E( (*it).first )->Read( prop_key, ts );
              for ( pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   dit2 = pxyz_data.find( (*pit) );
                  (*dit2).second.reserve(12);
                   if ( dim == 3U )
                     for ( size_t k=0; k<3; k++ )
                       for ( size_t l=0; l<3; l++ ) (*dit2).second.push_back( ts(k,l) );
                   else {
                      (*dit2).second.push_back( ts(0,0) );
                      (*dit2).second.push_back( ts(0,1) );
                      (*dit2).second.push_back( 0.0 );
                      (*dit2).second.push_back( ts(1,0) );
                      (*dit2).second.push_back( ts(1,1) );
                      (*dit2).second.push_back( 0.0 );
                      (*dit2).second.push_back( 0.0 );
                      (*dit2).second.push_back( 0.0 );
                      (*dit2).second.push_back( 0.0 );
                   }
                }
           }
      }
     
 } // end ElementPointData


template class Tecplot_Interface<1U>;
template class Tecplot_Interface<2U>;
template class Tecplot_Interface<3U>;
 
} // end namespace csp  
