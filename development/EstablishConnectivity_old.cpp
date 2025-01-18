template<size_t dim>
template<template<size_t> class CELL>
void  MeshManager<dim>::BuildConnectivity( typename deque<CELL<dim>*>::iterator first,
                                           typename deque<CELL<dim>*>::iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( distance(first,last) == 0U ) {
         csmp_error.notice( WARNING, "MeshManager<dim>::BuildConnectivity:", "supplied cell vector is empty; nothing was done." );
         return;
      }
    cout << "\nMeshManager<"<< dim <<">::BuildConnectivity: Establishing CSMP FE neighbor connectivity...\n";
 
    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    cout << "\n\tBuilding a multimap of the faces of the cells...\n";
    //       key             face number,neighbor
    multimap<set<Node<dim>*>,pair<size_t,CELL<dim>*> >  volume_neighbor_keys,
                                                        surface_neighbor_keys,
                                                        line_neighbor_keys;
    vector<size_t>                 fnids;
    typename std::set<Node<dim>*>  key; // region, boundary and split boundary all use nodes

    while ( first != last ) {
          const size_t faces((*first)->Faces());
          for ( size_t face=0U; face<faces; ++face )
            {
               if ( (*first) == nullptr ) {
                    csmp_error.notice( ERROR, "MeshManager<dim>::BuildConnectivity:",
                                      "supplied input range contains null pointers; nothing was done." );
                    return;
                 }
               // creating face key of node pointers from indices of face nodes
               (*first)->FE()->NodesOfFace( face, fnids );
               const size_t nodes(fnids.size());
               for ( size_t j=0U; j<nodes; ++j )
                 key.insert( (*first)->N( fnids[j] ) );
                 
               // inserting newly generated keys into multimap
               if ( (*first)->IsVolumeElement() )
                 volume_neighbor_keys.insert( make_pair( key, make_pair( face, (*first) ) ) );
               else if ( (*first)->IsSurfaceElement() )
                 surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*first) ) ) );
               else // for all line elements
                 line_neighbor_keys.insert( make_pair( key, make_pair( face, (*first) ) ) );
               key.clear();
            }
          first++;
        }

    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    cout << "  (Re)building neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() )
      {

csmp_error.notice( ERROR, "MeshManager::BuildConnectivity:", "line element neighbor connectivity calculation most likely faulty.");

        CELL<dim>* e1Ptr(nullptr);
        CELL<dim>* e2Ptr(nullptr);

        cout << "\n\t\tline elements...";
        //                key                  n-face, neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,CELL<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                              it2(line_neighbor_keys.begin());
        it2++;

        while ( it2 != line_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first )
                {
                   assert( (*it1).second.second != nullptr );
                   assert( (*it2).second.second != nullptr );
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning the two cells face neighbors to one another
                   //                        face pointer  nbor face idx   neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // line elements
    
    // 2.2 surface elements
    // --------------------
    if ( dim >= 2U and !surface_neighbor_keys.empty() )
      {
        cout << "\n\t\tsurface elements...";
        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,CELL<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                                 it2(surface_neighbor_keys.begin());
        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first )
                {
                    assert( (*it1).second.second != nullptr );
                    assert( (*it2).second.second != nullptr );
                    CELL<dim>* e1Ptr((*it1).second.second);
                    CELL<dim>* e2Ptr((*it2).second.second);
                    if ( e1Ptr != e2Ptr ) {
                         ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                         ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                      }
                    else csmp_error.notice( WARNING, "MeshManager::BuildConnectivity:",
                                            "discovered potentially duplicate surface simplex.");
                    ++it1;
                    ++it2;
                }
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // surface elements
      
    // 2.3 volume elements
    // -------------------
    if ( dim == 3U and !volume_neighbor_keys.empty() ) {

        CELL<dim>* e1Ptr(nullptr);
        CELL<dim>* e2Ptr(nullptr);

        cout << "\n\t\tvolume elements...\n";
        //                key             n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,CELL<dim>*> >::iterator it1(volume_neighbor_keys.begin()),
                                                                              it2(volume_neighbor_keys.begin());
        it2++;

        while ( it2 != volume_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first )
                {
                   assert( (*it1).second.second != nullptr );
                   assert( (*it2).second.second != nullptr );
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr );
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == volume_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=3
    
 } // end RebuildConnectivity