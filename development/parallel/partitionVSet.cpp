#include "partitionVSet.h"
#include "ICEM_ElementSpecifications.h"
#include "FEM_Traits.h"

// break VSet into subsets on the basis of supplied ID's of the elements in the target 
// partitions

#define CSP_GLOBAL_DEBUGGING

using namespace std;

namespace csp {

/*F <H4>Function:</H4><CODE>
<!------------------------------------------------------------------------>
void createPfverts( VSet<csp_float,dim>&, 
                    bool isoparametric, 
                    const char* etypes_file );
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Establishes which elements form the neighbors of each elements. If there
are no neighbors the neighbor entry in the 'pfverts' vector is assigned
a negative value of REGION_BOUNDARY=-28. It is up to the user to 
apply more diagnostics here. - One could for instance establish which
elements correspond to the element but are located in a different part
of the model. <p>

<H4>Input Arguments:</H4><!----------------------------------------------->

The function adds the 'pfverts' information to the VSet supplied as
first argument. <p>

The second argument determines whether the element
types in the model shall be interpreted as isoparametric or not.
This is only of consequence if the face specifications for the 
isoparametric version of an element are different than for the 
analytically integrated element. <p>

The third method argument is the full name of the file from which 
the element specifications shall be read. This file must be present
in the directory of the executable. A basic version of it can be 
generated from the information in the header file where the prototype
of this function is specified. <p>

<H4>Output Arguments &amp; Return Value</H4><!---------------------------->

The method returns the modified VSet. <p>

<H4>Implementation:</H4><!------------------------------------------------>

Ordered sets are used as keys which describe the nodes which make
up the faces of the finite elements. This face list is then used 
to identify neighboring elements. <p>

<H4>Application:</H4><!--------------------------------------------------->

The function can be used to add missing 'pfverts' data to a finite 
element connectivity dataset stored in the supplied VSet. <p> 

<!------------------------------------------------------------------------>
tested: running draft1 SKM-20/5/2005 */
template<typename fT, stl_index dim>
void createPfverts( VSet<fT,dim>& vset, bool isoparametric, const char* etypes_file ) 
 {
    // 0. getting the finite element specifications from file
    list<FEM_Traits>  fems;
    ifstream  ifs(etypes_file);
    char str[256];
    ifs.getline( str, 256, '\n' );
    ifs.getline( str, 256, '\n' ); 
    FEM_Traits  file_trait;
    while ( file_trait.InitializeFromFileStream( ifs ) )
      fems.push_back( file_trait );
    ifs.close();

    // 1. making a vector of CSP finite-element types for fast lookup
    ICEM_ElementSpecifications  especs;
    vector<FEM_Traits>          csp_etypes( EXPERIMENTAL_ELEMENT+1U );
    int32                       etype;

    for ( list<FEM_Traits>::const_iterator it=fems.begin(); it!=fems.end(); it++ ) {
         csp_etypes[ etype=especs.CSP_TypeFromICEMType( (*it).etype.c_str(), isoparametric, dim ) ] = (*it);
         if ( !vset.MixedElementTypeMesh() and etype == vset.ElementType(1U) ) break; 
      }

    // 2. making an element correspondance list ordered by face keys
    typedef pair<pair<stl_index,stl_index>,pair<stl_index,stl_index> > etuple;
    map<set<stl_index>,etuple>  neighbor_elements;
    for ( stl_index i=1U; i<=vset.Elements(); i++ )
      {
 
         const int32 elmt = (vset.MixedElementTypeMesh()) ? vset.ElementType(i) : vset.ElementType(1U);
          stl_index   eface(0);

         // looping over the faces creating keys from their nodes
          for ( vector<vector<int> >::const_iterator 
                it=csp_etypes[elmt].FaceNodesBegin(); it!=csp_etypes[elmt].FaceNodesEnd(); it++, eface++ )
            {
               set<stl_index> face;
               // collecting the global face node ids into a set
               for ( vector<int>::const_iterator jt=(*it).begin(); jt!=(*it).end(); jt++ )
                 face.insert( vset.Plist( i, static_cast<stl_index>(*jt) ) );
               // if face can be inserted it is the first time and the parent element id
               // is stored as first value, else it is stored as second value
               pair<map<set<stl_index>,etuple>::iterator,bool>
               fit = neighbor_elements.insert( make_pair( face, make_pair( make_pair(i,eface), make_pair(0,0) ) ) );
               if ( !fit.second ) (*fit.first).second.second = make_pair(i,eface);
            }
      }

    // 3. pre-flagging all neighbors as internal boundaries
    cout <<"\ncreatePfverts: assigning newly computed neighbor element data..."<< endl;
    for ( deque<vector<int32> >::iterator it=vset.PfvertsBegin(); it!=vset.PfvertsEnd(); it++ )
      for ( vector<int32>::iterator pit=(*it).begin(); pit!=(*it).end(); pit++ )
        (*pit) = IRREGULAR_OUTSIDE;

    // 4. assigning neighbor information where possible
    for ( map<set<stl_index>,etuple>::const_iterator 
          it=neighbor_elements.begin(); it!=neighbor_elements.end(); it++ )
      // if the face is shared between two elements
      if ( (*it).second.first.first > 0 and (*it).second.second.first > 0 ) {
            // the neigbour information is conveyed to both of these
            // first element
            vset.Pfvert( (*it).second.first.first, (*it).second.first.second, (*it).second.second.first );
            // second element
            vset.Pfvert( (*it).second.second.first, (*it).second.second.second, (*it).second.first.first );
        }
        
} // end createPfverts




//cerr <<"\ncreatePfverts: newly created Pfverts data:"<< endl;
//int i(1);
//for ( deque<vector<int32> >::const_iterator it=vset.PfvertsBegin(); it!=vset.PfvertsEnd(); it++, i++ ) {
//     cerr <<"\nelmt.nbors "<< i <<": ";
//     for ( vector<int32>::const_iterator pit=(*it).begin(); pit!=(*it).end(); pit++ )
//       cerr << (*pit) <<" ";
//  }

// is instantiated automatically when the next method is called
//template void createPfverts( VSet<csp_float,2U>&, bool, const char* ); 
//template void createPfverts( VSet<csp_float,3U>&, bool, const char* ); 








/*F <H4>Function:</H4><CODE>
<!------------------------------------------------------------------------>
bool partitionVSet( VSet<fT,dim>&,
                    const vector<vector<stl_index> >&, 
                    vector<VSet<fT,dim> >&, 
       map<pair<stl_index,stl_index>,vector<pair<stl_index,stl_index> > >& ) 
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Uses the supplied element identification numbers (1..n) to divide the 
input VSet into subsets with a new node / element numbering.
Once the new vsets have been created, maps of pairs of local node
numbers (1..n) are initialized, establishing the correspondance of node
numbering across adjacent VSet partitions. These maps are stored in a new
map using the implicit ID's of the partitions as keys. Thus, for instance,
the map of nodes shared between regions 1 and 2 will have the unique access
key (1,2). <p>

<H4>Input Arguments:</H4><!----------------------------------------------->

The first function argument is a reference to the VSet which shall be 
partitioned. <p>

The second argument is a reference to a vector of vectors which must 
hold the global element ID numbers (1..n) of the groups of elements which shall
make the model partitions. <p>

<H4>Output Arguments &amp; Return Value</H4><!---------------------------->

The VSet partitions and node-correspondance maps are returned into the 
second and third method arguments, respectively. <p>

<H4>Implementation:</H4><!------------------------------------------------>

Partitions are implicitly numbered 0..n-1.<p>
Node in the node correspondance maps are numbered 1..n.<p>

<H4>Messages:</H4><!------------------------------------------------------>

The function checks the integrity of the supplied VSet before it starts
processing it. Any problems identified by these checks are reported. <p>

The element ID vectors are checked as well to make sure that they do not
contain element ID's that are outside the range of the input VSet and 
that they cover all elements which are in the supplied VSet. <p>

<!------------------------------------------------------------------------>
tested: running draft1 DC - 20-03-06 */
template<typename fT, stl_index dim>
bool partitionVSet( VSet<fT,dim>& vset, // input VSet
                    const vector<vector<stl_index> >& element_based_partitions, // target regions
		    std::vector<std::set<std::pair<stl_index,int> > >&  outerhalo,
                    std::vector<VSet<csp_float,dim> >&  vsets,	// to store sub-VSets
                    VSetConnectivity<fT,dim>& VSetConnectivity,
		    bool isoparametric ) // information about the model that shall be created
 {
    // --------------------------------------------------------------------
    // 1. initial checks and cleanup
    // --------------------------------------------------------------------
    cout <<"\npartitionVSet: checking the VSet before it is going to be partitioned..."<< endl;
    if ( vset.Vertices() == 0U  or vset.Elements() == 0 ) {
         cerr <<"\npartitionVSet: input VSet appears to be empty; returning."<< endl;
         return false;
      }
    if ( element_based_partitions.empty() ) {
         cerr <<"\npartitionVSet: input partition map is empty; returning."<< endl;
         return false;
      }
    else {
        uint32  counter(0U);
	    for ( vector<vector<stl_index> >::const_iterator
	          it=element_based_partitions.begin(); it!=element_based_partitions.end(); it++, counter++ )
	      if ( (*it).empty() ) {
	           cerr <<"\npartitionVSet: input partition "<< counter <<" (0..n-1) is empty; returning."<< endl;
	           return false;
	        }
     }
    // duplicates and completeness check
    stl_index       elmts(0U);
    set<stl_index>  unique_ids;
    bool            duplicates_found(false);
    for ( vector<vector<stl_index> >::const_iterator
          it=element_based_partitions.begin(); it!=element_based_partitions.end(); it++ )

//DC: comment out: Halo element should be accepted
/*
    for ( vector<stl_index>::const_iterator t=(*it).begin(); t!=(*it).end(); t++ ) {
           pair<set<stl_index>::iterator,bool> inserted = unique_ids.insert( *t );
           if ( !inserted.second ) {
                cerr <<"\npartitionVSet: duplicate element ID = "<< *t << " found int partition map."<< endl;
                duplicates_found = true;
             }
           elmts++;
        }
    if ( duplicates_found ) {
         cerr <<"\npartitionVSet: partitioning cannot be performed because of duplicate element IDs."<< endl;
         return false;
      }

    if ( elmts != vset.Elements() ) {
         cerr <<"\npartitionVSet: partitioning cannot be performed because ID collections in partition list are not complete."<< endl;
         return false;
      }
*/
    
    
    unique_ids.erase( unique_ids.begin(), unique_ids.end() );
    
    // scrutinizing the input VSet
    if ( !vset.CheckFixPData() ) {
          cerr <<"\npartitionVSet: input vset failed consistency checks; returning."<< endl;
          return false;
      }
    if ( vset.ConstraintPoints() != 0U ) {
          cerr <<"\npartitionVSet: function does not work for VSets with constraint points; returning."<< endl;
          return false;
      }


    // cleaning up the output vsets vector if necessary
    if ( !vsets.empty() ) {
         cerr <<"\npartitionVSet: input VSet vector was not empty; resetting it."<< endl;
         vsets.erase( vsets.begin(), vsets.end() );
      }
    vsets.resize( element_based_partitions.size() );
 
    
    
    // --------------------------------------------------------------------
    // 2. Initializing the new VSets
    // --------------------------------------------------------------------
    cout <<"\npartitionVSet: creating "<< element_based_partitions.size() <<" new VSets containing the partitions of the model..."<< endl;
    // case A: mono-element mesh (PElements record is empty)
    stl_index                           processed_vsets(0U);
    deque<stl_index>                    mixed_ele_plist, mixed_ele_pfvertlist;
    vector<map<vector<fT>,stl_index> >  partition_coordinates( element_based_partitions.size() );
    vector<fT>                          node_coordinate(dim);
    map<stl_index,stl_index>            node_id_map;
    std::vector<stl_index>              first_outerhalo;
    first_outerhalo.resize(element_based_partitions.size());

      for ( vector<vector<stl_index> >::const_iterator
            it=element_based_partitions.begin(); it!=element_based_partitions.end(); it++, processed_vsets++ )
        {
	   // put outerhalo in local vector
	   vector<stl_index> loc_outerhalo;
	   loc_outerhalo.resize(outerhalo[processed_vsets].size());
	   unsigned int j(0);
	   for (set<pair<stl_index,int> >::const_iterator i=outerhalo[processed_vsets].begin();
	        i!=outerhalo[processed_vsets].end();i++,j++)
             loc_outerhalo[j]=(*i).first;

           // 2.1 resizing the VSet and setting element types (O.K.)
           // --------------------------------------------------------------------
           cout <<"\n\tsetting up a new VSet for partition "<< processed_vsets + 1U <<" with ";
           // counting the nodes, plist and pfvert entries
           if ( !node_id_map.empty() ) node_id_map.erase( node_id_map.begin(), node_id_map.end() );
           stl_index  new_node_id(1U);
           for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++ )
	     {
               for ( vector<stl_index>::const_iterator i=vset.PlistBegin( *eid ); i!=vset.PlistEnd( *eid ); i++ )
                 {
		   vector<stl_index>::const_iterator f=find(loc_outerhalo.begin(),loc_outerhalo.end(),(*i));
		   if(f==loc_outerhalo.end()) // node is not an outerhalo node
                    {
		      // old and new node id's in plist are numbered from 1..n
                      pair<map<stl_index,stl_index>::const_iterator,bool> nit=node_id_map.insert( make_pair(*i,new_node_id) );
                      if ( nit.second ) new_node_id++;
		    }
                 }
               mixed_ele_plist.push_back( vset.PlistSize( *eid ) );
               mixed_ele_pfvertlist.push_back( vset.PfvertsSize( *eid ) );
             }

	   // Store 1st outerhalo node ID in vector and continu with outerhalo nodes
           first_outerhalo[processed_vsets]=new_node_id;

           for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++ )
             // second loop to add outerhalo nodes
             for ( vector<stl_index>::const_iterator i=vset.PlistBegin( *eid ); i!=vset.PlistEnd( *eid ); i++ )
              {
		vector<stl_index>::const_iterator f=find(loc_outerhalo.begin(),loc_outerhalo.end(),(*i));
		if(f!=loc_outerhalo.end()) // node is an outerhalo node
                 {
//cout<<"Outerhalo node, global id = "<<(*i)<<", local id = "<<new_node_id<<", x = "<<vset.Px(*i-1U)<<", y = "<<vset.Py(*i-1U)<< endl;
                    // old and new node id's in plist are numbered from 1..n
                    pair<map<stl_index,stl_index>::const_iterator,bool> nit=node_id_map.insert( make_pair(*i,new_node_id) );
                    if ( nit.second ) new_node_id++;
		 }
              }

cout<<"\n\nRank "<<processed_vsets<<", first outerhalo = "<<first_outerhalo[processed_vsets]<< endl;

	   vsets[processed_vsets].ResizeNodes( node_id_map.size() );
           vsets[processed_vsets].ResizePlist( mixed_ele_plist );
           vsets[processed_vsets].ResizePfverts( mixed_ele_pfvertlist );
           mixed_ele_plist.erase( mixed_ele_plist.begin(), mixed_ele_plist.end() );
           mixed_ele_pfvertlist.erase( mixed_ele_pfvertlist.begin(), mixed_ele_pfvertlist.end() );

           // conveying the information about element type
           if ( !vset.MixedElementTypeMesh() ) {
                 vector<int32>  elmt_types( 1U, vset.ElementType(1U) );
                 vsets[processed_vsets].AddElementTypes( elmt_types.begin(), elmt_types.end() );
             }
           else {
               vector<int32>  etypes;
               etypes.reserve( vsets[processed_vsets].Elements() );
               for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++ )
                 etypes.push_back( vset.ElementType( *eid ) );
               vsets[processed_vsets].AddElementTypes( etypes.begin(), etypes.end() );
               etypes.erase( etypes.begin(), etypes.end() );
            }
           cout << vsets[processed_vsets].Vertices() <<" nodes and "<< vsets[processed_vsets].Elements() <<" elements."<< endl;
           

           // 2.2 assigning 'plist' entries and node coordinates 'px,py,pz' (O.K.)
           // --------------------------------------------------------------------
           cout <<"\n\tassigning 'plist' entries and node coordinates..."<< endl;
           deque<vector<stl_index> >::iterator plt(vsets[processed_vsets].PlistBegin());
           for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++, plt++ ) {
	             // looping over the individual entries in the plist
	             for ( vector<stl_index>::iterator 
	                   pit1=vset.PlistBegin( *eid ), pit2=(*plt).begin(); pit1!=vset.PlistEnd( *eid ); pit1++, pit2++ )
	               {
	                   // assigning plist entries, but with a new node numbering (0..n-1)
	                   map<stl_index,stl_index>::const_iterator  key = node_id_map.find( *pit1 );
#ifdef CSP_GLOBAL_DEBUGGING	                 
	                   assert( key != node_id_map.end() );
#endif	        
                       // assigning new (local) node ids to the plist    
                       *pit2 = (*key).second; 
#ifdef CSP_GLOBAL_DEBUGGING	                 
                       assert( *pit2 > 0U and *pit2 <= vsets[processed_vsets].Vertices() );
#endif	            
	                   // assigning node coordinates (here coordinates are repeatedly overwritten=inefficient)
	                   // NB: Px,Py,Pz() functions take entries from 0..nodes-1
	                                    vsets[processed_vsets].Px( *pit2-1U, node_coordinate[0]=vset.Px( *pit1-1U ) );
	                   if ( dim != 1U ) vsets[processed_vsets].Py( *pit2-1U, node_coordinate[1]=vset.Py( *pit1-1U ) );
	                   if ( dim == 3U ) vsets[processed_vsets].Pz( *pit2-1U, node_coordinate[2]=vset.Pz( *pit1-1U ) );
	                   
	                   // remembering node coordinates in the partition
	                   partition_coordinates[processed_vsets].insert( make_pair( node_coordinate, *pit2 ) );
	               }
             }


           // 2.3 copying material properties (using unique node id's created above) (O.K.)
           // -----------------------------------------------------------------------------
           cout <<"\n\tcopying the material properties..."<< endl;
           // scalar variables
           if ( vset.ScalarPropertiesBegin() != vset.ScalarPropertiesEnd() )
             for ( typename map<string,FEM_Data<ScalarVariable<fT> > >::const_iterator
                   ct=vset.ScalarPropertiesBegin(); ct!=vset.ScalarPropertiesEnd(); ct++ ) 
               {
                   stl_index  n(0U);
                   cout <<"\n\tprocessing scalar variable: '"<< (*ct).first <<"'"<< endl;
                   bool add_data(true);
                   // reading the dataset
                   FEM_Data<ScalarVariable<fT> >  data1;
                   vset.Data( (*ct).first.c_str(), data1 );
                   // building the target dataset
                   const stl_index dsize = ( data1.Placement() == NODE ) ? vsets[processed_vsets].Vertices() 
                                                                         : vsets[processed_vsets].Elements();
                   // new dataset
                   FEM_Data<ScalarVariable<fT> >  data2( data1.Placement(), dsize );
                   // copying data selectively
                   if ( data1.Placement() == ELEMENT )
                     for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++ )
                       data2[n++] = data1[ *eid - 1U ];
                   else if ( data1.Placement() == NODE )
                     for ( map<stl_index,stl_index>::const_iterator t=node_id_map.begin(); t!=node_id_map.end(); t++ )
                        data2[(*t).second - 1U] = data1[ (*t).first - 1U ];
                   else {
                        cerr <<"\n\tplacement of scalar-type input data '"<< (*ct).first <<"' was not recognized."<< endl;
                        add_data = false;
                     }
                   if ( add_data ) vsets[processed_vsets].AddData( (*ct).first.c_str(), data2 );
               }
           // vector variables
           if ( vset.VectorPropertiesBegin() != vset.VectorPropertiesEnd() )
             for ( typename map<string,FEM_Data<VectorVariable<fT,dim> > >::const_iterator
                   ct=vset.VectorPropertiesBegin(); ct!=vset.VectorPropertiesEnd(); ct++ ) 
               {
                   stl_index  n(0U);
                   cout <<"\n\tprocessing vector variable: '"<< (*ct).first <<"'"<< endl;
                   bool add_data(true);
                   FEM_Data<VectorVariable<fT,dim> >  data1;
                   vset.Data( (*ct).first.c_str(), data1 );
                   const stl_index dsize = ( data1.Placement() == NODE ) ? vsets[processed_vsets].Vertices() 
                                                                         : vsets[processed_vsets].Elements();
                   FEM_Data<VectorVariable<fT,dim> >  data2( data1.Placement(), dsize );
                   if ( data1.Placement() == ELEMENT )
                     for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++ )
                       data2[n++] = data1[ *eid - 1U ];
                   else if ( data1.Placement() == NODE )
                     for ( map<stl_index,stl_index>::const_iterator t=node_id_map.begin(); t!=node_id_map.end(); t++ )
                       data2[(*t).second - 1U] = data1[ (*t).first -1U ];
                   else {
                        cerr <<"\n\tplacement of vector-type input data '"<< (*ct).first <<"' was not recognized."<< endl;
                        add_data = false;
                     }
                   if ( add_data ) vsets[processed_vsets].AddData( (*ct).first.c_str(), data2 );
               }
           // tensor variables
           if ( vset.TensorPropertiesBegin() != vset.TensorPropertiesEnd() )
             for ( typename map<string,FEM_Data<TensorVariable<fT,dim> > >::const_iterator
                   ct=vset.TensorPropertiesBegin(); ct!=vset.TensorPropertiesEnd(); ct++ ) 
               {
                   stl_index  n(0U);
                   cout <<"\n\tprocessing tensor variable: '"<< (*ct).first <<"'"<< endl;
                   bool add_data(true);
                   FEM_Data<TensorVariable<fT,dim> >  data1;
                   vset.Data( (*ct).first.c_str(), data1 );
                   const stl_index dsize = ( data1.Placement() == NODE ) ? vsets[processed_vsets].Vertices() 
                                                                         : vsets[processed_vsets].Elements();
                   FEM_Data<TensorVariable<fT,dim> >  data2( data1.Placement(), dsize );
                   if ( data1.Placement() == ELEMENT )
                     for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++ )
                       data2[n++] = data1[ *eid - 1U ];
                   else if ( data1.Placement() == NODE )
                     for ( map<stl_index,stl_index>::const_iterator t=node_id_map.begin(); t!=node_id_map.end(); t++ )
                       data2[(*t).second - 1U] = data1[ (*t).first -1U ];
                   else {
                        cerr <<"\n\tplacement of tensor-type input data '"<< (*ct).first <<"' was not recognized."<< endl;
                        add_data = false;
                     }
                   if ( add_data ) vsets[processed_vsets].AddData( (*ct).first.c_str(), data2 );
               }
             
             
           // 2.4 applying boundary flags 'pbflags' and 'pbvalues' values (O.K.)
           // ------------------------------------------------------------------
           std::vector<stl_index> OuterhaloAtBoundary;
           OuterhaloAtBoundary.resize(0);
           cout <<"\n\tcopying boundary flags and potential boundary values..."<< endl;
           typename map<stl_index,fT>::const_iterator vt=vset.BValuesBegin();
           for ( map<stl_index,int32>::const_iterator bt=vset.BFlagsBegin(); bt!=vset.BFlagsEnd(); bt++, vt++ )
             {
                // boundary flags have indices 1..n
                                 node_coordinate[0] = vset.Px( (*bt).first-1U ); 
                if ( dim != 1U ) node_coordinate[1] = vset.Py( (*bt).first-1U ); 
                if ( dim == 3U ) node_coordinate[2] = vset.Pz( (*bt).first-1U ); 
                
                // search corresponding local boundary node and make store it in the new VSet
                typename map<vector<fT>,stl_index>::const_iterator 
                sit = partition_coordinates[processed_vsets].find( node_coordinate );
                if ( sit != partition_coordinates[processed_vsets].end() )
                 {
                      vsets[processed_vsets].AddBFlag( (*sit).second, (*bt).second );
                      vsets[processed_vsets].AddBValue( (*sit).second, (*vt).second );
                      if((*sit).second>=first_outerhalo[processed_vsets])
                        OuterhaloAtBoundary.push_back((*sit).second);
                 }  
              }

          // All Outerhalo nodes need to be flagged IRREGULAR_OUTSIDE

          for ( std::deque<std::vector<stl_index> >::const_iterator el_iter=vsets[processed_vsets].PlistBegin(); el_iter!=vsets[processed_vsets].PlistEnd(); el_iter++ )
           for (stl_index j=0;j!=(*el_iter).size();j++)
              if(((*el_iter)[j])>=first_outerhalo[processed_vsets])
               {
                 // if outerhalo nodes is also a global region boundary
                 // then do nothing.
                 // otherwise assign IRREGULAR_OUTSIDE
                 std::vector<stl_index>::iterator f=find(OuterhaloAtBoundary.begin(),OuterhaloAtBoundary.end(),(*el_iter)[j]);
                 if(f==OuterhaloAtBoundary.end())
                  {
                    vsets[processed_vsets].AddBFlag( ((*el_iter)[j]), IRREGULAR_OUTSIDE );
                    vsets[processed_vsets].AddBValue( ((*el_iter)[j]), 1. );
                  }
               }
           // 2.5 creating and transferring 'pfverts' information 
           // --------------------------------------------------------------------
           cout <<"\n\tcopying and creating element neighbor information ('pfverts')..."<< endl;
           createPfverts( vsets[processed_vsets], isoparametric );            
           // adding the old pfverts information (pfverts must have been resized already)
           assert( vsets[processed_vsets].PfvertsBegin() != vsets[processed_vsets].PfvertsEnd() );
           stl_index  elmt_id(1U); 
           for ( vector<stl_index>::const_iterator eid=(*it).begin(); eid!=(*it).end(); eid++ ) {
	             // looping over the entries in the pfverts
	             vector<int32>::iterator pft2=vsets[processed_vsets].PfvertsBegin( elmt_id++ ); 
	             for ( vector<int32>::const_iterator pft1=vset.PfvertsBegin( *eid ); pft1!=vset.PfvertsEnd( *eid ); pft1++, pft2++ )
	               // if the old element neighbor was missing, the model boundary information is conveyed to new VSet
	               if ( *pft1 < 0 ) *pft2 = *pft1; 
             }

       } // end loop over partitions
     
     
     
    // --------------------------------------------------------------------
    // 3. Creating the node correspondance maps (node ids 1..n)
    // --------------------------------------------------------------------
    
    // DC: HERE INNER AND OUTER HALO NODES NEED TO BE RENUMBERED IN THE SAME WAY
    // SO ALSO MAKE THE SAME MAPS
    std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >  correspondance_map;

    cout <<"\npartitionVSet: creating node point correspondance maps for the shared boundaries of the new VSets..."<< endl;
    if ( !correspondance_map.empty() )  correspondance_map.erase( correspondance_map.begin(), correspondance_map.end() );
    //      node_id-p1,node_id_p2
    set<pair<stl_index,stl_index> >     shared_nodes;
    vector<pair<stl_index,stl_index> >  empty_vec;
    for ( stl_index i=0U; i<partition_coordinates.size(); i++ )
      for ( stl_index j=0U; j<partition_coordinates.size(); j++ ) 
        if ( i != j ) {
			 // finding the nodes which are shared among the partitions i and j
			 for ( typename map<vector<fT>,stl_index>::const_iterator 
			       pt=partition_coordinates[i].begin(); pt!=partition_coordinates[i].end(); pt++ ) {
			       typename map<vector<fT>,stl_index>::const_iterator foundit=partition_coordinates[j].find( (*pt).first );
			       if ( foundit != partition_coordinates[j].end() )
 		             shared_nodes.insert( make_pair( (*pt).second, (*foundit).second ) );
			   }
             
             // inserting the node-id pairs into the correspondence list                   
             if ( !shared_nodes.empty() ) {
                  //   key(partition#,partition#), shared boundary nodes(#-in part1, #-in part2)
                  pair<map<pair<stl_index,stl_index>,vector<pair<stl_index,stl_index> > >::iterator,bool>
                  cit = correspondance_map.insert( make_pair( make_pair(i,j), empty_vec ) );
                  assert( cit.second );
                  (*cit.first).second.reserve( shared_nodes.size() );
                  // now the correspondance map pairs are initialized
                  for ( set<pair<stl_index,stl_index> >::const_iterator 
                        sht=shared_nodes.begin(); sht!=shared_nodes.end(); sht++ )
                    (*cit.first).second.push_back( (*sht) );

                  shared_nodes.erase( shared_nodes.begin(), shared_nodes.end() );
               }
          }
    VSetConnectivity.Initialize(correspondance_map, first_outerhalo, outerhalo );
    VSetConnectivity.RemoveHaloToHaloConnections();

    cout <<"\npartitionVSet: Vset partitioned successfully."<< endl;
    return true;
    
 } // end partitionVSet
 
 
 
 

//cerr <<"\npartition coordinates["<< i <<"]: ";
//for ( typename map<vector<fT>,stl_index>::const_iterator 
//      pt=partition_coordinates[i].begin(); pt!=partition_coordinates[i].end(); pt++ ) 
//  cerr << (*pt).second <<" ("<< (*pt).first[0] <<","<< (*pt).first[1] <<") ";

template
bool partitionVSet( VSet<csp_float,2>&, 
                    const vector<vector<stl_index> >&, 
                    vector<set<pair<stl_index,int> > >&,
                    vector<VSet<csp_float,2> >&,		    
                    VSetConnectivity<csp_float,2>&,
		    bool );
template
bool partitionVSet( VSet<csp_float,3>&, 
                    const vector<vector<stl_index> >&, 
                    vector<set<pair<stl_index,int> > >&,
                    vector<VSet<csp_float,3> >&,		    
                    VSetConnectivity<csp_float,3>&,
		    bool );


} // end namespace csp
