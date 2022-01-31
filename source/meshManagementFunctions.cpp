QUESTIONS / TODO: 
  - A mesh patch needs a root element if no element neighbors in the patch are shared, however, a node may still be shared between two patches so that a root node traversal might find all elements. Not sure how to handle this? - I am tempted to only have root element pointers because the FE is smallest building block of a model.
 
DECISION: ONLY ELEMENT TRAVERSAL IS NEEDED


========================================================================================
    IDEAS
========================================================================================

DELETER

std::vector<int*> v;

template<typename T>
struct deleter : std::unary_function<const T*, void>
{
  void operator() (const T *ptr) const
  {
    delete ptr;
  }
};

// call deleter for each element , freeing them
std::for_each (v.begin (), v.end (), deleter<int> ());
v.clear ();























========================================================================================
      OLD CODE
========================================================================================

/**

Visit all elements without relying on their storage in a container.
The Idx numbering of elements and nodes is not altered by this method.

@return number of elements that were discovered.

@attention, this method always gets called when a region is first formed.

@attention nodes must have been assigned their parent elements for this method to work.

@author SKM 9/20/2008, CSMP Castasegna workshop, Switzerland.

*/
template<size_t dim>
size_t Region<dim>::AccumulateAll( const csmp::Node<dim>* root_node,
                                   bool reestablishNeighborConnectivity )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( root_node == nullptr )
       csmp_error.notice( FATAL_ERROR, "Region<dim>::AccumulateAll:", "Root node pointer is dangling!" );

     if ( root_node->Parent(0) == nullptr )
       csmp_error.notice( FATAL_ERROR, "Region<dim>::AccumulateAll:",
                         "Root node must have been assigned parent elements; else this method cannot operate." );

     if ( !this->elmt_vec_.empty() )
       csmp_error.notice( WARNING, "Region<dim>::AccumulateAll:",
                                   "Region is not empty; deleting all content." );
     this->elmt_vec_.clear();

     // 1. traversal of the existing mesh nodes to find all its elements
     set<csmp::Element<dim>*>       explored_elements;
     set<const csmp::Node<dim>*>    discovered_nodes;
     deque<const csmp::Node<dim>*>  current_nodes;
     // starting at the root element
     discovered_nodes.insert( root_node );
     current_nodes.push_back( root_node );

     // MESH TRAVERSAL
     while ( !current_nodes.empty() ) {
          const csmp::Node<dim>*  n_ptr(*current_nodes.begin());
          // for all parent elements of the current node
          for ( size_t i=0U; i<n_ptr->Parents(); i++ ) {
               // for all the nodes of each parent element
               for ( size_t j=0U; j<n_ptr->Parent(i)->Nodes(); j++ )
                 // if this node is not the one from which we started
                 if ( j != n_ptr->ParentNodeNumber(i) ) {
                      pair<typename set<const csmp::Node<dim>*>::iterator,bool>
                        new_node = discovered_nodes.insert( n_ptr->Parent(i)->N(j) );
                      if ( new_node.second ) current_nodes.push_back( n_ptr->Parent(i)->N(j) );
                   }
               // storing the explored element
               explored_elements.insert( n_ptr->Parent(i) );
            }
          // removing the node from the discovered (but not yet explored) deque
          current_nodes.pop_front();
       }

    // 2. assigning and trimming excess storage from the element pointer vector
    this->elmt_vec_.assign( explored_elements.begin(), explored_elements.end() );
    vector<csmp::Element<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

    // 3. creating pointers to the nodes of the identified elements
    set<csmp::Node<dim>*>  node_set;
    for ( typename vector<csmp::Element<dim>*>::const_iterator
          it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ )
      for ( typename vector<csmp::Node<dim>*>::const_iterator
            nit=(*it)->NodesBegin(); nit!=(*it)->NodesEnd(); nit++ )
        node_set.insert( (*nit) );

    // TODO: use emplace here?
    this->node_vec_.assign( node_set.begin(), node_set.end() );

    // 4. (re)connecting elements up to their neighbors
    //    TODO: this is a very time-consuming step; is there a speed-up?
    if ( reestablishNeighborConnectivity )
      this->EstablishNeighborConnectivity();

    // 5. identifying the boundaries
    this->IdentifyPerimeter();

    return this->elmt_vec_.size();

 } // end AccumulateAll








 
 
  
 // ------------------------------------
 // Junchul - method to find root Faces
 // ------------------------------------
      // traversal of existing mesh to find all its faces
      deque<set<Face<dim>*>>	explored_face_groups;
      set<Face<dim>*>			    discovered_faces;
      deque<Face<dim>*>		    current_faces;
      map<size_t, Face<dim>*> faces_map;
      for ( auto f : face_connector )
        faces_map[ f->Idx() ] = f;

      // starting at the first face		
      current_faces.push_back( faces_map.begin()->second );

      size_t group_idx = 0U;
      while ( !faces_map.empty() )
      {
        while ( !current_faces.empty() ) {
          const Face<dim>*  f_ptr( *current_faces.begin() );
          // for all neighbor faces of the current face
          for ( size_t i = 0U; i < f_ptr->Neighbors(); i++ ) {
              if ( f_ptr->Neighbor( i ) == nullptr ) continue;

              // if this neighbor is new one					
              auto new_face = discovered_faces.insert( f_ptr->Neighbor( i ) );
              if ( new_face.second ) {
                current_faces.push_back( f_ptr->Neighbor( i ) );
                faces_map.erase( f_ptr->Neighbor( i )->Idx() );
              }
            }
          // removing the face from the discovered (but not yet explored) deque
          current_faces.pop_front();
        }
        if ( discovered_faces.empty() ) {
            auto first_face = faces_map.begin()->second;
            discovered_faces.insert( first_face );
            explored_face_groups.push_back( discovered_faces );
            faces_map.erase( first_face->Idx() );
            discovered_faces.clear();
            group_idx++;
            continue;
          }
        explored_face_groups.push_back( discovered_faces );

        // if there is only one single face left, it creates a new group of faces.
        if ( faces_map.size() == 1 ) {
            discovered_faces.clear();
            auto first_face = faces_map.begin()->second;
            discovered_faces.insert( first_face );
            explored_face_groups.push_back( discovered_faces );
            break;
          }
        if ( !faces_map.empty() ) {
            auto first_face = faces_map.begin()->second;
            current_faces.push_back( first_face );
          }

        discovered_faces.clear();
        group_idx++;
      }

      for ( auto group : explored_face_groups )
        root_face_group_.push_back( (*group.begin()) );






// JUNCHUL CODE TO FIND ROOT INTERFACE OBJECTS

      // traversal of the existing mesh nodes to find all its interfaces
      deque<set<InterFace<dim>*>>	 explored_interface_groups;
      set<InterFace<dim>*>		     discovered_interfaces;
      deque<InterFace<dim>*>		   current_interfaces;
      map<size_t, InterFace<dim>*> interfaces_map;
      for ( auto f : interface_connector )
        interfaces_map[f->Idx()] = f;

      // starting at the first interface		
      current_interfaces.push_back( interfaces_map.begin()->second );

      size_t group_idx = 0U;
      while ( !interfaces_map.empty() )
      {
        while ( !current_interfaces.empty() ) {
          const InterFace<dim>*  n_ptr( *current_interfaces.begin() );
          // for all neighbor interfaces of the current interface
          for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
            if ( n_ptr->Neighbor( i ) == nullptr ) continue;

            // if this neighbor is a new one					
            auto new_interface = discovered_interfaces.insert( n_ptr->Neighbor( i ) );
            if ( new_interface.second ) {
              current_interfaces.push_back( n_ptr->Neighbor( i ) );
              interfaces_map.erase( n_ptr->Neighbor( i )->Idx() );
            }
          }
          // removing the interface from the discovered (but not yet explored) deque
          current_interfaces.pop_front();
        }
        if ( discovered_interfaces.size() == 0 ) {
          auto first_interface = interfaces_map.begin()->second;
          discovered_interfaces.insert( first_interface );
          explored_interface_groups.push_back( discovered_interfaces );
          interfaces_map.erase( first_interface->Idx() );
          discovered_interfaces.clear();
          group_idx++;
          continue;
        }
        explored_interface_groups.push_back( discovered_interfaces );

        // if there is only one single interface left, it creates a new group of interfaces.
        if ( interfaces_map.size() == 1 ) {
          discovered_interfaces.clear();
          auto first_interface = interfaces_map.begin()->second;
          discovered_interfaces.insert( first_interface );
          explored_interface_groups.push_back( discovered_interfaces );
          break;
        }
        if ( interfaces_map.size() > 0 ) {
          auto first_interface = interfaces_map.begin()->second;
          current_interfaces.push_back( first_interface );
        }

        discovered_interfaces.clear();
        group_idx++;
      }

      for ( auto group : explored_interface_groups )
        root_interface_group_.push_back( (*group.begin()) );
        
        
        
        
   if ( !node_connector.empty() ) {
    // assgining the root pointers
    // traversal of the existing mesh nodes to find all its nodes
    deque<set<Node<dim>*>>	explored_node_groups;
    set<Node<dim>*>			    discovered_nodes;
    deque<Node<dim>*>		    current_nodes;
    map<size_t, Node<dim>*> nodes_map;
    // TODO: why not use the node connector directly, it is already sorted, down to line 761
    for ( auto n : node_connector )
      nodes_map[n->Idx()] = n;

    // starting at the first node
    current_nodes.push_back( nodes_map.begin()->second );

    size_t group_idx = 0U;
    while ( !nodes_map.empty() )
    {
      while ( !current_nodes.empty() ) {
        auto n_ptr( *current_nodes.begin() );
        // for all neighbor nodes of the current node
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one	
          // TODO: is the purpose to flag a node that has already been visited?				
          auto new_node = discovered_nodes.insert( n_ptr->Neighbor( i ) );
          if ( new_node.second ) {
            current_nodes.push_back( n_ptr->Neighbor( i ) );
            nodes_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the node from the discovered (but not yet explored) deque
        current_nodes.pop_front();
      }
      if ( discovered_nodes.size() == 0 ) {
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        nodes_map.erase( first_node->Idx() );
        discovered_nodes.clear();
        group_idx++;
        continue;
      }
      // avoid the infinite iterative looping search
      bool infinite_looping( false );
      for ( size_t i = 0U; i < explored_node_groups.size(); i++)
        {
          if ( explored_node_groups[i] == discovered_nodes ) {
            infinite_looping = true;
            if ( csmp_error.Verbose() )
              cerr << "\nMeshManager<" << dim << ">::Initialize: node connections might be wrong..." << endl;
          }
        }

      if ( infinite_looping ) {
        discovered_nodes.clear();
        if ( nodes_map.size() > 0 ) {
          auto first_node = nodes_map.begin()->second;
          discovered_nodes.insert( first_node );
          nodes_map.erase( first_node->Idx() );
        }        
      }      
      
      explored_node_groups.push_back( discovered_nodes );      

      // if there is only one single node left, it creates a new group of nodes.
      if ( nodes_map.size() == 1 ) {
        discovered_nodes.clear();
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        if ( csmp_error.Verbose() )
          cerr << "\nMeshManager<" << dim << ">::Initialize: the node (" << first_node->Idx() << ") doesn't have any parents and neighbours..." << endl;
        break;
      }
      if ( nodes_map.size() > 0 ) {
        auto first_node = nodes_map.begin()->second;
        current_nodes.push_back( first_node );
      }

      discovered_nodes.clear();
      group_idx++;
    }

    for ( auto group : explored_node_groups )
      root_node_group_.push_back( (*group.begin()) );

  } // end nodes

  set<Element<dim>*>	explored_elmt_groups;
  if ( elmt_connector.size()>0 ) {
    // assgining the root pointers
    for ( auto root_node : root_node_group_ )
    {
      for ( size_t i = 0; i < root_node->Parents(); i++ ) {
        auto parent = root_node->Parent( i );
        if ( parent ) {
          explored_elmt_groups.insert( parent );
          break;
        }
      }
    }
  } // end elements

  for ( auto root_elmt : explored_elmt_groups )
    root_elmt_group_.push_back( root_elmt );

  if ( root_node_group_.size() != root_elmt_group_.size() ) {
    root_node_group_.clear();
    for ( auto root_elmt : root_elmt_group_ )
      root_node_group_.push_back( root_elmt->N( 0 ) );
  }
       
        
        
// MESH MODIFICATION JCK

/**
Removes specific elements according to the set of elmt_numbers.

@param  region_name The name of region which includes the elements which shall be removed.
@param  elmt_numbers the element numbers which shall be removed.
@return the number of the removed elements is returned.

@attention JCK code - watch out how this is called.

*/
template<size_t dim, template<size_t> class REGION_COMPLEX>
size_t RegionInterface<dim, REGION_COMPLEX>::RemoveElements( const char* region_name, const std::set<long>& elmt_numbers )
{
  // check whether region exists (should be a notice only, nothrow)
  if ( !ContainsRegion( region_name ) )
    throw csmp::Exception( WARNING,
                           "RegionsInterface<dim,REGION_COMPLEX>::RemoveElements",
                           "region did not exist: ",
                           region_name );


cerr<<"\n\n\n\n\n\nCalled RegionInterface<dim, REGION_COMPLEX>::RemoveElements: from JCK - watch for side effects.\n\n\n\n\n\n"<< endl;

// TODO: use symmetric difference to find new element vector for region
// TODO: ask mesh manager to remove range of elements

  // finding the region in the corresponding map
  typename map<string, csmp::Region<dim> >::iterator iterUniqueRegion( uniqueGroupMap_.find( string( region_name ) ) );



  // deleting the elements according to elmt_numbers
  auto& subdomain = iterUniqueRegion->second;
  REGION_COMPLEX<dim>* regionComplex( static_cast<REGION_COMPLEX<dim>* >(this) ) ;
  auto& meshMgr = regionComplex->Mesh();  
  auto& elementVector = subdomain.CellVector();
  const csmp::Index key_enr = regionComplex->Database().StorageKey( "element number" );

  long removed_elements( 0 );
  std::set<Node<dim>*> candidate_nodes;  
  for ( size_t i = 0U; i < elementVector.size(); i++ ) {
    // 1.1 Remove this element from its neighbour's connections
    Element<dim>* e = elementVector[i];
    if ( e->AtBoundary() != NOT && e->AtBoundary() != IRREGULAR ) continue;
    long idx = static_cast<long>(e->Read( key_enr ));
    if ( elmt_numbers.find( idx ) == elmt_numbers.end() ) continue;    
    
    // 1.2. find candidate nodes belonging to this element
    for ( size_t j = 0U; j < e->Nodes(); j++ ) {
      Node<dim>* n = e->N( j );
      candidate_nodes.insert( n );
    }    
    meshMgr.Erase( e );
    elementVector.erase( elementVector.begin() + i );
    elementVector.swap( elementVector );
    removed_elements++;
    i--;
  }
  for ( size_t i = 0U; i < elementVector.size(); i++ ) {
    Element<dim>* e = elementVector[i];
    // remove this element from its neighbour's connections
    auto& neighbourVector = e->NeighborElementVector();
    for ( size_t j = 0U; j < neighbourVector.size(); j++ ) {
      if ( neighbourVector[j] == NULL ) continue;
      auto& nnVector = neighbourVector[j]->NeighborElementVector();
      for ( size_t k = 0U; k < nnVector.size(); k++ ) {
        if ( nnVector[k] == NULL ) continue;
        auto& nn = nnVector[k];        
        if ( nn->Nodes() == 0 ) {
          nn = NULL;
        }
      }
    }
  }
  for(Node<dim>* n: candidate_nodes )
  {
    if ( n->Parents() == 0 )
      meshMgr.Erase( n );
  }
  subdomain.CreateNodePointerVector();
  subdomain.EstablishNeighborConnectivity();
  subdomain.IdentifyPerimeter();

  // remove those elements from the default region 'Model'
  const char* default_model = "Model";
  if ( ContainsRegion( default_model ) ) {
    // finding the region in the corresponding map
    typename std::map<std::string, csmp::Region<dim> >::iterator iterRegion( groupMap_.find( std::string( default_model ) ) );

    auto& defulat_domain = iterRegion->second;
    auto& elementVector = defulat_domain.ElementVector();

    for ( size_t i = 0U; i < elementVector.size(); i++ ) {
      // remove this element 
      Element<dim>* e = elementVector[i];
      if ( e->Nodes() == 0 ) {
        elementVector.erase( elementVector.begin() + i );
        elementVector.swap( elementVector );
        i--;
      }
    }
    for ( size_t i = 0U; i < elementVector.size(); i++ ) {
      Element<dim>* e = elementVector[i];
      // remove this element from its neighbour's connections
      auto& neighbourVector = e->NeighborElementVector();
      for ( size_t j = 0U; j < neighbourVector.size(); j++ ) {
        auto& ne = neighbourVector[j];
        if ( ne == NULL ) continue;
        if ( ne->Nodes() == 0 ) {
          ne = NULL;
        }
      }
    }
    defulat_domain.CreateNodePointerVector1();
    defulat_domain.EstablishNeighborConnectivity();
    defulat_domain.IdentifyPerimeter();
  }

  meshMgr.RebuildParentRelationships( subdomain.NodesBegin(), subdomain.NodesEnd() );
  
  return removed_elements;
  
} // end RemoveElements


// JUNCHUL'S CODE

/**
Deletes the corresponding face
*/
template<size_t dim>
void MeshManager<dim>::Erase( Face<dim>* face )
{
  // if the found face is the root face
  // change the root face into one of its neighbors and delete the found face.
  Face<dim>* new_root_face( nullptr );
  Face<dim>* root_face( nullptr );
  if ( root_face_group_.size() > 0 ) {
    for ( auto root : root_face_group_ )
      if ( face == root ) root_face = root;

    if ( root_face ) {
      for ( auto e : root_face->NeighborElementVector() ) {
        if ( root_face != e ) {
          new_root_face = e;
          break;
        }
      }
      if ( new_root_face == nullptr )
        cerr << "MeshManager::Erase(Face): this is the root face which cannot be deleted. \n";
      else
        root_face = new_root_face;
    }
  }

  // delete the face
  if ( face ) {
    // unassign the connections of its neighbors
    for ( auto n : face->NeighborElementVector() )
      if ( n != nullptr ) n->Unassign( face );

    delete face;
    face = nullptr;
    n_faces_--;
  }
  
} // end Erase( Face pointer )







/**
Deletes the corresponding interface
*/
template<size_t dim>
void MeshManager<dim>::Erase( InterFace<dim>* interface )
{
  // if the found interface is the root interface
  // change the root interface into one of its neighbors and delete the found interface.
  InterFace<dim>* new_root_interface( nullptr );
  InterFace<dim>* root_interface( nullptr );
  if ( root_interface_group_.size() > 0 ) {
    for ( auto root : root_interface_group_ )
      if ( interface == root ) root_interface = root;

    if ( root_interface ) {
      for ( auto e : root_interface->NeighborElementVector() ) {
        if ( root_interface != e ) {
          new_root_interface = e;
          break;
        }
      }
      if ( new_root_interface == nullptr )
        cerr << "MeshManager::Erase(InterFace): this is the root interface which cannot be deleted. \n";
      else
        root_interface = new_root_interface;
    }
  }

  // delete the interface
  if ( interface ) {
    // unassign the connections of its neighbors
    for ( auto n : interface->NeighborElementVector() )
      if ( n != NULL ) n->Unassign( interface );

    delete interface;
    interface = nullptr;
    n_interfaces_--;
  }
}


/**
    Update  root pointers of the mesh after modifying the mesh.
*/
template<size_t dim>
void MeshManager<dim>::Update( std::deque<Node<dim>*> nodes, std::deque<Element<dim>*> elmts )
{
  // forming root pointers of the nodes for contiguous subdomains
  if ( nodes.size()>0 ) {
    // finding the root pointers		
    deque<set<Node<dim>*>>	explored_node_groups;
    set<Node<dim>*>			discovered_nodes;
    deque<Node<dim>*>		current_nodes;
    map<size_t, Node<dim>*> nodes_map;
    for ( auto n : nodes )
      nodes_map[n->Idx()] = n;

    // starting at the first node
    current_nodes.push_back( nodes_map.begin()->second );

    size_t group_idx = 0U;
    while ( !nodes_map.empty() )
    {
      while ( !current_nodes.empty() ) {
        Node<dim>*  n_ptr( *current_nodes.begin() );
        // for all neighbor nodes of the current node
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_node = discovered_nodes.insert( n_ptr->Neighbor( i ) );
          if ( new_node.second ) {
            current_nodes.push_back( n_ptr->Neighbor( i ) );
            nodes_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the node from the discovered (but not yet explored) deque
        current_nodes.pop_front();
      }
      if ( discovered_nodes.size() == 0 ) {
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        nodes_map.erase( first_node->Idx() );
        discovered_nodes.clear();
        group_idx++;
        continue;
      }
      explored_node_groups.push_back( discovered_nodes );

      // if there is only one single node left, it creates a new group of nodes.
      if ( nodes_map.size() == 1 ) {
        discovered_nodes.clear();
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        break;
      }
      if ( nodes_map.size() > 0 ) {
        auto first_node = nodes_map.begin()->second;
        current_nodes.push_back( first_node );
      }

      discovered_nodes.clear();
      group_idx++;
    }

    //update the node root pointers
    std::deque<Node<dim>*> updated_root_node_group;
    for ( auto group : explored_node_groups )
    {
      deque<Node<dim>*> nodes_vec;
      nodes_vec.assign( group.begin(), group.end() );
      sort( nodes_vec.begin(), nodes_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_node_group.push_back( nodes_vec.front() );
      nodes_vec.clear();
    }

    root_node_group_.swap( updated_root_node_group );
  } // end nodes	

  if ( elmts.size()>0 ) {
    // assgining the root pointers
    set<Element<dim>*>	explored_elmt_groups;
    for ( auto root_node : root_node_group_ )
    {
      for ( size_t i = 0; i < root_node->Parents(); i++ ) {
        auto parent = root_node->Parent( i );
        if ( parent ) {
          explored_elmt_groups.insert( parent );
          break;
        }
      }
    }

    //update the element root pointers
    std::deque<Element<dim>*> updated_root_elmt_group;
    for ( auto root_elmt : explored_elmt_groups )
      updated_root_elmt_group.push_back( root_elmt );

    root_elmt_group_.swap( updated_root_elmt_group );

    if ( root_node_group_.size() != root_elmt_group_.size() ) {
      root_node_group_.clear();
      for ( auto root_elmt : root_elmt_group_ )
        root_node_group_.push_back( root_elmt->N( 0 ) );
    }
  } // end elements
}


/**
       seems to also create connectivity for the Face and InterFace objects
*/
template<size_t dim>
void MeshManager<dim>::Update()
{
  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( *this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root faces to find all its faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( *this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( *this, interfaces );
  sort( interfaces.begin(), interfaces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // forming root pointers of the nodes for contiguous subdomains
  if ( nodes.size()>0 ) {
    // finding the root pointers		
    deque<set<Node<dim>*>>	explored_node_groups;
    set<Node<dim>*>			discovered_nodes;
    deque<Node<dim>*>		current_nodes;
    map<size_t, Node<dim>*> nodes_map;
    for ( auto n : nodes )
      nodes_map[n->Idx()] = n;

    // starting at the first node
    current_nodes.push_back( nodes_map.begin()->second );

    size_t group_idx = 0U;
    while ( !nodes_map.empty() )
    {
      while ( !current_nodes.empty() ) {
        Node<dim>*  n_ptr( *current_nodes.begin() );
        // for all neighbor nodes of the current node
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_node = discovered_nodes.insert( n_ptr->Neighbor( i ) );
          if ( new_node.second ) {
            current_nodes.push_back( n_ptr->Neighbor( i ) );
            nodes_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the node from the discovered (but not yet explored) deque
        current_nodes.pop_front();
      }
      if ( discovered_nodes.size() == 0 ) {
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        nodes_map.erase( first_node->Idx() );
        discovered_nodes.clear();
        group_idx++;
        continue;
      }
      explored_node_groups.push_back( discovered_nodes );

      // if there is only one single node left, it creates a new group of nodes.
      if ( nodes_map.size() == 1 ) {
        discovered_nodes.clear();
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        break;
      }
      if ( nodes_map.size() > 0 ) {
        auto first_node = nodes_map.begin()->second;
        current_nodes.push_back( first_node );
      }

      discovered_nodes.clear();
      group_idx++;
    }

    //update the node root pointers
    std::deque<Node<dim>*> updated_root_node_group;
    for ( auto group : explored_node_groups )
    {
      deque<Node<dim>*> nodes_vec;
      nodes_vec.assign( group.begin(), group.end() );
      sort( nodes_vec.begin(), nodes_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_node_group.push_back( nodes_vec.front() );
      nodes_vec.clear();
    }

    root_node_group_.swap( updated_root_node_group );
  } // end nodes	

  if ( elmts.size()>0 ) {
    // assgining the root pointers
    set<Element<dim>*>	explored_elmt_groups;
    for ( auto root_node : root_node_group_ )
    {
      for ( size_t i = 0; i < root_node->Parents(); i++ ) {
        auto parent = root_node->Parent( i );
        if ( parent ) {
          explored_elmt_groups.insert( parent );
          break;
        }
      }
    }

    //update the element root pointers
    std::deque<Element<dim>*> updated_root_elmt_group;
    for ( auto root_elmt : explored_elmt_groups )
      updated_root_elmt_group.push_back( root_elmt );

    root_elmt_group_.swap( updated_root_elmt_group );

    if ( root_node_group_.size() != root_elmt_group_.size() ) {
      root_node_group_.clear();
      for ( auto root_elmt : root_elmt_group_ )
        root_node_group_.push_back( root_elmt->N( 0 ) );
    }
  } // end elements

    // forming root pointers of the faces for contiguous regions
  if ( faces.size() > 0 ) {
    // finding the root pointers
    // traversal of the existing mesh nodes to find all its faces
    deque<set<Face<dim>*>>		explored_face_groups;
    set<Face<dim>*>				discovered_faces;
    deque<Face<dim>*>			current_faces;
    map<size_t, Face<dim>*>		faces_map;
    for ( auto f : faces )
      faces_map[f->Idx()] = f;

    // starting at the first face		
    current_faces.push_back( faces.front() );

    size_t group_idx = 0U;
    while ( !faces_map.empty() )
    {
      while ( !current_faces.empty() ) {
        Face<dim>*  n_ptr( *current_faces.begin() );
        // for all neighbor faces of the current face
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_face = discovered_faces.insert( n_ptr->Neighbor( i ) );
          if ( new_face.second ) {
            current_faces.push_back( n_ptr->Neighbor( i ) );
            faces_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the face from the discovered (but not yet explored) deque
        current_faces.pop_front();
      }
      if ( discovered_faces.size() == 0 ) {
        auto first_face = faces_map.begin()->second;
        discovered_faces.insert( first_face );
        explored_face_groups.push_back( discovered_faces );
        faces_map.erase( first_face->Idx() );
        discovered_faces.clear();
        group_idx++;
        continue;
      }
      explored_face_groups.push_back( discovered_faces );

      // if there is only one single face left, it creates a new group of faces.
      if ( faces_map.size() == 1 ) {
        discovered_faces.clear();
        auto first_face = faces_map.begin()->second;
        discovered_faces.insert( first_face );
        explored_face_groups.push_back( discovered_faces );
        break;
      }
      if ( faces_map.size() > 0 ) {
        auto first_face = faces_map.begin()->second;
        current_faces.push_back( first_face );
      }

      discovered_faces.clear();
      group_idx++;
    }

    // update the face root pointers
    std::deque<Face<dim>*> updated_root_face_group;
    for ( auto group : explored_face_groups )
    {
      deque<Face<dim>*> faces_vec;
      faces_vec.assign( group.begin(), group.end() );
      sort( faces_vec.begin(), faces_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_face_group.push_back( faces_vec.front() );
      faces_vec.clear();
    }
    root_face_group_.swap( updated_root_face_group );
  } // end faces

    // forming root pointers of the interfaces for contiguous regions
  if ( interfaces.size() > 0 ) {
    // finding the root pointers
    // traversal of the existing mesh nodes to find all its interfaces
    deque<set<InterFace<dim>*>>		explored_interface_groups;
    set<InterFace<dim>*>			discovered_interfaces;
    deque<InterFace<dim>*>			current_interfaces;
    map<size_t, InterFace<dim>*>	interfaces_map;
    for ( auto f : interfaces )
      interfaces_map[f->Idx()] = f;

    // starting at the first interface		
    current_interfaces.push_back( interfaces.front() );

    size_t group_idx = 0U;
    while ( !interfaces_map.empty() )
    {
      while ( !current_interfaces.empty() ) {
        InterFace<dim>*  n_ptr( *current_interfaces.begin() );
        // for all neighbor interfaces of the current interface
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_interface = discovered_interfaces.insert( n_ptr->Neighbor( i ) );
          if ( new_interface.second ) {
            current_interfaces.push_back( n_ptr->Neighbor( i ) );
            interfaces_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the interface from the discovered (but not yet explored) deque
        current_interfaces.pop_front();
      }
      if ( discovered_interfaces.size() == 0 ) {
        auto first_interface = interfaces_map.begin()->second;
        discovered_interfaces.insert( first_interface );
        explored_interface_groups.push_back( discovered_interfaces );
        interfaces_map.erase( first_interface->Idx() );
        discovered_interfaces.clear();
        group_idx++;
        continue;
      }
      explored_interface_groups.push_back( discovered_interfaces );

      // if there is only one single interface left, it creates a new group of interfaces.
      if ( interfaces_map.size() == 1 ) {
        discovered_interfaces.clear();
        auto first_interface = interfaces_map.begin()->second;
        discovered_interfaces.insert( first_interface );
        explored_interface_groups.push_back( discovered_interfaces );
        break;
      }
      if ( interfaces_map.size() > 0 ) {
        auto first_interface = interfaces_map.begin()->second;
        current_interfaces.push_back( first_interface );
      }

      discovered_interfaces.clear();
      group_idx++;
    }

    // update the face root pointers
    std::deque<InterFace<dim>*> updated_root_face_group;
    for ( auto group : explored_interface_groups )
    {
      deque<InterFace<dim>*> interfaces_vec;
      interfaces_vec.assign( group.begin(), group.end() );
      sort( interfaces_vec.begin(), interfaces_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_face_group.push_back( interfaces_vec.front() );
      interfaces_vec.clear();
    }
    root_interface_group_.swap( updated_root_face_group );
  } // end interfaces
}


**
         Checks whether node already exists (in which case it returns it) or creates a new disconnected node and returns a pointer to it.

        @todo Node should be added in the context of creating a new element or element, face, interface that it is connected to
        @todo if there is a collocated new node this should be reflected in the Manifold
        
         @note checks 
*/
template<size_t dim> // TODO: Misleading name, is this function useful?
Node<dim>* MeshManager<dim>::AddIfUnique( Node<dim>& node )
{
  set<Element<dim>*>		discovered_elements;
  set<const Node<dim>*>	discovered_nodes;
  deque<const Node<dim>*>	current_nodes;
  for ( auto root_node : root_node_group_ ) {
    discovered_nodes.insert( root_node );
    current_nodes.push_back( root_node );
    while ( !current_nodes.empty() ) {
      const Node<dim>*  n_ptr( *current_nodes.begin() );
      for ( size_t i = 0U; i < n_ptr->Parents(); i++ ) {
        for ( size_t j = 0U; j < n_ptr->Parent( i )->Nodes(); j++ ) {
          if ( j != n_ptr->ParentNodeNumber( i ) ) {
            if ( *(n_ptr->Parent( i )->N( j )) == node )
              return n_ptr->Parent( i )->N( j );
            pair<typename set<const Node<dim>*>::const_iterator, bool>
              new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
            if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
          }
        }
        discovered_elements.insert( n_ptr->Parent( i ) );
      }
      current_nodes.pop_front();
    }
  }

  Node<dim>* new_node = new Node<dim>( node );
  *new_node = node;
  n_nodes_++;

  return new_node;
}


/**
Inserts the corresponding element
*/
template<size_t dim>
Element<dim>* MeshManager<dim>::AddIfUnique( Element<dim>& elmt )
{
  set<Element<dim>*>		discovered_elements;
  set<const Node<dim>*>	discovered_nodes;
  deque<const Node<dim>*>	current_nodes;
  for ( auto root_node : root_node_group_ ) {
    // starting at the root element
    discovered_nodes.insert( root_node );
    current_nodes.push_back( root_node );
    while ( !current_nodes.empty() ) {
      const Node<dim>*  n_ptr( *current_nodes.begin() );
      for ( size_t i = 0U; i < n_ptr->Parents(); i++ ) {
        for ( size_t j = 0U; j < n_ptr->Parent( i )->Nodes(); j++ ) {
          if ( j != n_ptr->ParentNodeNumber( i ) ) {
            pair<typename set<const Node<dim>*>::const_iterator, bool>
              new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
            if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
          }
        }
        if ( *(n_ptr->Parent( i )) == elmt ) return n_ptr->Parent( i );
        discovered_elements.insert( n_ptr->Parent( i ) );
      }
      current_nodes.pop_front();
    }
  }

  Element<dim>* new_elmt = new Element<dim>( elmt );
  *new_elmt = elmt;
  n_elmts_++;

  return new_elmt;
}


/**
Inserts the corresponding face
*/
template<size_t dim>
Face<dim>* MeshManager<dim>::AddIfUnique( Face<dim>& face )
{
  set<Face<dim>*>		discovered_faces;
  deque<Face<dim>*>	current_faces;

  Face<dim>* found_face = NULL;

  for ( auto root_face : root_face_group_ ) {
    discovered_faces.insert( root_face );
    current_faces.push_back( root_face );
    while ( !current_faces.empty() || found_face != NULL ) {
      Face<dim>*  n_ptr( *current_faces.begin() );
      for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
        if ( n_ptr->Neighbor( i ) == NULL ) continue;
        if ( *(n_ptr->Neighbor( i )) == face )
          return n_ptr->Neighbor( i );
        auto new_face = discovered_faces.insert( n_ptr->Neighbor( i ) );
        if ( new_face.second ) current_faces.push_back( n_ptr->Neighbor( i ) );
      }
      current_faces.pop_front();
    }
  }

  Face<dim>* new_face = new Face<dim>( face );
  *new_face = face;
  n_faces_++;
  return new_face;
}


/**
Inserts the corresponding interface
*/
template<size_t dim>
InterFace<dim>* MeshManager<dim>::AddIfUnique( InterFace<dim>& interface )
{
  set<InterFace<dim>*>	discovered_interfaces;
  deque<InterFace<dim>*>	current_interfaces;

  for ( auto root_interface : root_interface_group_ ) {
    discovered_interfaces.insert( root_interface );
    current_interfaces.push_back( root_interface );
    while ( !current_interfaces.empty() ) {
      InterFace<dim>*  n_ptr( *current_interfaces.begin() );
      for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
        if ( n_ptr->Neighbor( i ) == NULL ) continue;
        if ( *(n_ptr->Neighbor( i )) == interface ) return n_ptr->Neighbor( i );
        auto new_interface = discovered_interfaces.insert( n_ptr->Neighbor( i ) );
        if ( new_interface.second ) current_interfaces.push_back( n_ptr->Neighbor( i ) );
      }
      current_interfaces.pop_front();
    }
  }

  InterFace<dim>* new_interface = new InterFace<dim>( interface );
  *new_interface = interface;
  n_interfaces_++;
  return new_interface;
}






// =====================================================================================
SLOT MAP IDEA (Sean Middleditch)
// =====================================================================================
// Try it out!

Object:
  int index
  int version
  other data

SlotMap:
  Object objects[]
  int slots[]
  int freelist[]
  int count

  Get(id):
    index = indirection[id.index]
    if objects[index].version = id.version:
      return &objects[index]
    else:
      return null

  CreateObject():
    index = freelist.pop()

    objects[count].index = id
    objects[count].version += 1

    indirection[index] = count

    Object* object = &objects[count].object
    object.initialize()

    count += 1

    return object

  Remove(id):
    index = indirection[id.index]
    if objects[index].version = id.version:
      objects[index].version += 1
      objects[count - 1].version += 1

      swap(objects[index].data, objects[count - 1].data)
The indirection layer allows you to have a stable identifier (the index into the indirection layer, where entries do not move) for a resource that can move during compaction (the main object list).

The version tag allows you to store an ID to an object that might be deleted. For example, you have the id (10,1). The object with index 10 is deleted (say, your bullet hits an object and is destroyed). The object in that location of memory in the main object list then has its version number bumped, giving it (10,2). If you try to look up (10,1) again from a stale ID, the lookup returns that object through index 10, but can see that the version number has changed, so the ID is no longer valid.

This is the absolute fastest data structure you can have with a stable ID that allows objects to move in memory, which is important for data locality and cache coherence. This is faster than any implementation of a hash table possible; a hash table at the very least needs to calculate a hash (more instructions than a table lookup) and then has to follow the hash chain (either a linked list in the horrible case of std::unordered_map, or an open-addressed list in any not-stupid implementation of a hash table), and then has to do a value compare on each key (no more expensive, but possible less expensive, than the version tag check). A very good hash table (not the one in any implementation of the STL, as the STL mandates a hash table that optimizes for different use cases than you game about for a game object list) might save on one indirection, but wil be a larger bloated structure that you can't iterate over at maximum efficiency, and loses out in real-world performance for use cases like this one.

There are various improvements you can make to the base algorithm. Using something like a std::deque for the main object list, for instance; one extra layer of indirection, but allows objects to be inserted into a full list without invalidating any temporary pointers you've acquired from the slotmap.

You can also avoid storing the index inside the object, as the index can be calculated from the object's memory address (this - objects), and even better is only needed when removing the object in which case you already have the object's id (and hence index) as a parameter.

Apologies for the write-up; I don't feel it's the clearest description it could be. It's late and it's difficult to explain without spending more time than I have on code samples.
Share
Improve this answer
Follow
edited Aug 9 '12 at 7:10
answered Aug 9 '12 at 2:41

Sean Middleditch
41.1k33 gold badges8585 silver badges129129 bronze badges
1
You are trading off an extra deref and a high alloc/free cost (swap) every access for 'compact' storage. In my experience with video games, that's a bad trade :) YMMV of course. – Jeff Gates Aug 9 '12 at 8:32 
1
You don't actually do the dereference that often in real world scenarios. When you do, you can store the returned pointer locally, especially if you use the deque variant or know you won't be creating new objects while you have the pointer. Iterating over the collections is a very expensive and frequent operation, you need the stable id, you want memory compaction for volatile objects (like bullets, particles, etc), and the indirection is very efficient on modem hardware. This technique is used in more than a few very high performance commercial engines. :) – Sean Middleditch Aug 9 '12 at 20:37
1
In my experience: (1) Video games are judged on worst case performance, not average case performance. (2) You normally have 1 iteration over a collection per frame, thus compacting simply 'makes your worst case less frequent'. (3) You often have many allocs/frees in a single frame, high cost means you limit that capability. (4) You have unbounded derefs per frame (in games I've worked on, including Diablo 3, often the deref was the highest perf cost op after moderate optimization, >5% of server load ). I don't mean to dismiss other solutions, just pointing out my experiences and reasoning! – Jeff Gates Aug 10 '12 at 20:05
3
I love this data structure. I'm surprised it's not more well-known. It's simple and solves all the problems that have been making me bang my head for months. Thanks for sharing. – jbatez Nov 23 '13 at 3:17
2
Any newbie reading this should be very wary of this advice. This is a very misleading answer. "The answer is always to use an array or std::vector. Types like a linked list or a std::map are usually absolutely horrendous in games, and that definitely includes cases like collections of game objects." is greatly exaggerated. There is no "ALWAYS" answer, otherwise these other containers wouldn't have been created. To say maps/lists are "horrendous" is also hyperbole. There are A LOT of video games which use these. "Most Efficient" is not "Most Practical" and can be misread as a subjective "Best". – user50286 May 20 '15 at 7:04

  
        
