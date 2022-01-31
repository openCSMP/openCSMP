
/// retrieves and returns the element ids of the first contiguous element patch that can be reached by mesh traversal from the starting element
template<size_t dim>
void floodFillViaIndexes( const Region<dim>&, size_t starting_idx,
                          std::set<size_t>& output_contiguous_subset );

/**
    Retrieves and returns the element ids of the first contiguous element patch
    that can be reached by mesh traversal from the starting element.
    
    @attention the elements in the region must have a unique numbering scheme.
*/
template<size_t dim>
void floodFillViaIndexes( const Region<dim>& gref, size_t starting_idx,
                          std::set<size_t>& elements_contiguous_subset )
 {
    assert( starting_idx < gref.Elements() );
    assert( gref.E(starting_idx) != NULL );
    // identifying the neighbors of the first element to be looked at
    deque<size_t>  neighbor_elements;
    const size_t  neighbors(gref.E(starting_idx)->Neighbors());
    for ( size_t i=0U; i<neighbors; i++ )
      if ( gref.E(starting_idx)->Neighbor(i) != NULL )
        neighbor_elements.push_back( gref.E(starting_idx)->Neighbor(i)->Idx() );
 
    // performing the floodfill, starting with an empty set
    if ( !elements_contiguous_subset.empty() )
      elements_contiguous_subset.clear();
   
    // insert the first element into the new subset
    elements_contiguous_subset.insert( starting_idx );
      
    // element set for subsequent passes
    deque<size_t>  new_neighbor_elements;
   
     while( !neighbor_elements.empty() )
       {
          // 1. loop over those neighbors that are not already part of the list
          for ( typename deque<size_t>::const_iterator
                idx=neighbor_elements.begin(); idx!=neighbor_elements.end(); ++idx )
            // if the element has not already been dealt with
            if ( elements_contiguous_subset.find( (*idx) ) == elements_contiguous_subset.end() ) {
                const size_t  neighbors(gref.E(*idx)->Neighbors());
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( size_t j=0U; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( gref.E(*idx)->Neighbor(j) != NULL ) {
                        assert( gref.E(*idx)->Neighbor(j)->Idx() < gref.Elements() );
                        new_neighbor_elements.push_back( gref.E(*idx)->Neighbor(j)->Idx() );
                        elements_contiguous_subset.insert( gref.E(*idx)->Neighbor(j)->Idx() );
                    }
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_elements = new_neighbor_elements;
            
          // 3. emptying neighbor set for the next loop
          new_neighbor_elements.clear();
      }
     
 } // end floodFillViaIndexes


template void floodFillViaIndexes( const Region<1U>&, size_t, set<size_t>& );
template void floodFillViaIndexes( const Region<2U>&, size_t, set<size_t>&  );
template void floodFillViaIndexes( const Region<3U>&, size_t, set<size_t>&  );





/**

Two versions, for Model and for Region computations are defined.
This method accumulates the basic 'Operands' of the finite-element equations
into the global solution matrix and the righthand vector. This is
done only if their condition flag (VARIABLE_FLAG) is equal to the essential
condition flag defined for this basic operand. (by default,
essential-condition flags of Operands are DIRICH(let)).

The method retrieves information from the mesh and the property managers.
If a group computation is carried out, the target group is accessed by
a pointer.

@section implementation Implementation

AssignEssentialConditions() does not decrease the size of the global
matrix 'G' or the righthand vector 'rh'. In stead, it zeros rows
in 'G' and then sets the G'ith element in question to one while the
Operand value is then placed into 'rh' (this technique is described for
fixed displacements (Dirichlet) in Segerlind, 1984, p. 417ff).

@section application Application

Because it simultaneously affects the solution matrix and the righthand
vector, AssignEssentialConditions() must be applied AFTER the accumulation
process has been completed.

When you assemble vector or tensor variables, you also have the option
of only assembling one of their components. To do this just flag the
components that you want to assemble as DIRICH.

@section messages Messages

The following errors and warnings will be reported:

@code
please call EstablishMatrixSetup() prior to this method.
@endcode

If EstablishMatrixSetup() has not been executed so far.

@code
No basic operands have been specified...
@endcode

If there are no dependent variables defined !!!

@code
Warning: So far no conditions are assigned to elements, faces, segments
@endcode

The dependent variable must be placed on the nodes:
*/
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions( const COMPUTATION_DOMAIN<dim>& gref )
 {
    if ( !setup_established_ )
      throw csmp::Exception( ERROR, "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
                      "please call EstablishMatrixSetup() prior to this method.");

    if ( basic_operands_.empty() )
      throw csmp::Exception( ERROR, "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
                             "No (basic) operands have been specified...");

    // ----------------------------------------------------
    // 2. if the PDE_IntegratorExperimental applies to a Region
    // ----------------------------------------------------
    //  This method accumulates basic 'Operands' in the finite-element equations
    //  into the global solution matrix and the righthand vector, if their condition
    //  flag is equal to the essential condition flag of the basic operand.
    //  (by default, the essential-condition flag of Operand's is DIRICH(let)).

    //  AssignEssentialConditions() must be applied AFTER accumulation process.
    //
    //  When you assemble vector or tensor variables, you also have the option
    //  of only assembling one of their components. To do this just set the
    //  components that you do not want to assemble to DBL_MAX.
    // -------------------------------------------------------------
     const size_t dim2(dim * dim);

     for ( operandsConstIterator
           it=test_operands_.begin(); it!=test_operands_.end(); it++ )
       {
         typename vector<csmp::Node<dim>*>::const_iterator  niter(gref.NodesBegin());
         csmp::Index prop_key = (*it).first.key;
         size_t      offset   = (*it).second;

         if ( prop_key.place != NODE )
           throw csmp::Exception( ERROR, "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
                                 "So far no conditions are assigned to elements, faces, segments");

          switch( prop_key.type )
           {
              case SCALAR:
                while ( niter != gref.NodesEnd() ) {
                      if ( (*niter)->Status( prop_key ) == DIRICH )
                        {
                           size_t  position = (*niter)->Idx() + offset;
                           G_.ZeroRow( position );
                           G_.Add( position, position, 1. );
                           rh_[ position ] = (*niter)->Read( prop_key );
                        }
                       niter++;
                    }
                break;
              case VECTOR: {
                  VectorVariable<dim>  vc;
                  while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, vc );
                        for ( size_t i=0U; i<dim; i++ )
                          if ( vc.Flag(i) == DIRICH ) {
                               size_t  position = (*niter)->Idx() * dim + i + offset;
                               G_.ZeroRow( position );
                               G_.Add( position, position, 1. );
                               rh_[ position ] = vc(i);
                            }
                        niter++;
                     }
                  }
                break;
              case TENSOR: {
                 TensorVariable<dim>  ts;
                 while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, ts );
                        for ( size_t i=0U; i<dim; i++ )
                          if ( ts.Flag(i) == DIRICH ) 
                            for ( size_t j=0U; j<dim; j++ )
                              {
                                size_t  position = (*niter)->Idx() * dim2 + i * dim + j + offset;
                                G_.ZeroRow( position );
                                G_.Add( position, position, 1. );
                                rh_[ position ] = ts(i,j);
                              }
                         ++niter;
                      }
                   }
                 break;
              case ARRAY: {
                  ArrayVariable  ar(prop_key.dataDepth);
                  while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, ar );
                        if ( ar.Flag() == DIRICH )
                        {
                            for ( size_t i=0U; i<prop_key.dataDepth; i++ )
                            {
                               size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                               G_.ZeroRow( position );
                               G_.Add( position, position, 1. );
                               rh_[ position ] = ar(i);
                            }
                        }
                        niter++;
                     }
                  }
                 break;
              case FLAGGEDARRAY: {
                  FlaggedArrayVariable  ar(prop_key.dataDepth);
                  while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, ar );
                        for ( size_t i=0U; i<prop_key.dataDepth; i++ )
                          if ( ar.Flag(i) == DIRICH ) {
                               size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                               G_.ZeroRow( position );
                               G_.Add( position, position, 1. );
                               rh_[ position ] = ar(i);
                            }
                        niter++;
                     }
                  }
                 break;
               default:
                 throw csmp::Exception( FATAL_ERROR,
                                       "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
                                       "Variable type not recognised by this method" );
        }
    } // end for

} // end AssignEssentialConditions



/**
     For a given vector of Dirichlet entries in the solution sparse matrix A and the righthand vector rhs
     this function performs a row / column elimination.
     The results are returned into SparseMatrix B.
     
     @param index_mapping integer mapping relating the row/column indices of the condensed matrix to those of the
     original matrix. This mapping is needed for the assignment of the solution results to the 
     Model.
*/
void eliminateDirichletConstraints( const map<size_t,double64>& Dirichlet_constraints,
                                    const csmp::SparseMatrix& A, vector<double64>& rhs,
                                    csmp::SparseMatrix& B, vector<long64>& index_mapping  )
{
    // 0. checking whether anything needs to be done
    if ( Dirichlet_constraints.empty() ) {
         B = A;
         index_mapping.resize( A.Rows() );
         iota( index_mapping.begin(), index_mapping.end(), 0U );
         return;
      }

    // 1. creating an old (full matrix) to new (eliminated matrix) index mapping
    const size_t dof(A.Rows());
    index_mapping.resize(dof);
    long64 new_index(0U);
  
    for ( size_t i=0U; i<dof; ++i ) {
         // if the index corresponds to a Dirichlet row, we can omitt it
         if ( Dirichlet_constraints.find(i) != Dirichlet_constraints.end() )
           index_mapping[i] = UNSPECIFIED;
         else {
             index_mapping[i] = new_index;
             new_index++;
          }
      }

    // 2. modifying the righthand vector taking into account the Dirichlet conditions
    for ( map<size_t,double64>::const_iterator
          it = Dirichlet_constraints.begin(); it != Dirichlet_constraints.end(); ++it )
      {
          // looping over rows, avoiding the zero elements
          const SparseMatrix::colsConstIterator rowEnd(A.RowEnd((*it).first));
          for ( SparseMatrix::colsConstIterator cit(A.RowBegin((*it).first)); cit!=rowEnd; ++cit )
            if ( (*cit).first != (*it).first )
              {
                 // dividing the remaining rhs elements by the column values from the eliminated rows
                 rhs[ (*cit).first ] = (*it).second * -A( (*cit).first, (*it).first );
              }
      }
  
     // 3. condensing the right-hand vector into non-Dirichlet elements only
     new_index = 0U;
     for ( size_t i=0U; i<dof; ++i )
      // if the index corresponds to a Dirichlet row, we can omitt it
      if ( index_mapping[i] != UNSPECIFIED ) {
            rhs[new_index] = rhs[i];
            new_index++;
         }
      rhs.resize( dof - Dirichlet_constraints.size() );
  
    // 4. generating the new condensed sparse matrix. Use the index_mapping to set matrix elements indices.
    if ( B.Entries() > 0U ) B.Erase();
    B.Resize( dof - Dirichlet_constraints.size() );
    size_t row(0U);
    // for all rows of the original matrix
    for ( size_t original_row=0U; original_row<dof; ++original_row ) {
        // if they are non zero
        if ( index_mapping[original_row] != UNSPECIFIED ) {
              const SparseMatrix::colsConstIterator rowEnd(A.RowEnd(original_row));
              for ( SparseMatrix::colsConstIterator cit(A.RowBegin(original_row)); cit!=rowEnd; ++cit )
                 // copying values if they are non zero
                if ( index_mapping[ (*cit).first ] != UNSPECIFIED )
                  B.Assign( row, index_mapping[ (*cit).first ], A( original_row, (*cit).first ) );
               row++;
          }
     }
  
} // end eliminateDirichletConstraints




    /// eliminates essential (Dirichlet) conditions, condensing the the solution matrix, rhs etc. to that of the remaining DOF
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::EliminateEssentialConditions( const COMPUTATION_DOMAIN<dim>& domain )
 {
    if ( !setup_established_ )
      throw csmp::Exception( ERROR, "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::EliminateEssentialConditions:",
                      "please call EstablishMatrixSetup() prior to this method.");

    if ( basic_operands_.empty() )
      throw csmp::Exception( ERROR, "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::EliminateEssentialConditions:",
                             "No (basic) operands have been specified...");

     // 1. establish the rows (and columns) to which Dirichlet conditions were applied and storing these in a map
     //    of index-value pairs
    map<size_t,double64>  Dirichlet_constraints;

    for ( operandsConstIterator
          it=test_operands_.begin(); it!=test_operands_.end(); it++ )
       {
         typename vector<csmp::Node<dim>*>::const_iterator  niter(domain.NodesBegin());
         const auto nodesEnd(domain.NodesEnd());
         csmp::Index prop_key = (*it).first.key;
         size_t      offset   = (*it).second;

         if ( prop_key.place != NODE )
           throw csmp::Exception( ERROR, "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::EliminateEssentialConditions:",
                                 "So far no conditions are assigned to elements, faces, segments");

          switch( prop_key.type )
           {
              case SCALAR:
                while ( niter != nodesEnd ) {
                      if ( (*niter)->Status( prop_key ) == DIRICH ) {
                           const size_t position = (*niter)->Idx() + offset;
                           Dirichlet_constraints.insert( make_pair( position, (*niter)->Read( prop_key ) ) );
                        }
                       niter++;
                    }
                break;
              case VECTOR: {
                  VectorVariable<dim>  vc;
                  while ( niter != nodesEnd ) {
                        (*niter)->Read( prop_key, vc );
                        for ( size_t i=0U; i<dim; i++ )
                          if ( vc.Flag(i) == DIRICH ) {
                               const size_t  position = (*niter)->Idx() * dim + i + offset;
                               Dirichlet_constraints.insert( make_pair( position, vc(i) ) );
                            }
                        niter++;
                     }
                  }
                break;
              case ARRAY: {
                  ArrayVariable  ar(prop_key.dataDepth);
                  while ( niter != nodesEnd ) {
                        (*niter)->Read( prop_key, ar );
                        if ( ar.Flag() == DIRICH ) {
                            for ( size_t i=0U; i<prop_key.dataDepth; i++ ) {
                                 const size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                                 Dirichlet_constraints.insert( make_pair( position, ar(i) ) );
                              }
                          }
                        niter++;
                     }
                  }
                 break;
              case FLAGGEDARRAY: {
                  FlaggedArrayVariable  ar(prop_key.dataDepth);
                  while ( niter != nodesEnd ) {
                        (*niter)->Read( prop_key, ar );
                        for ( size_t i=0U; i<prop_key.dataDepth; i++ )
                          if ( ar.Flag(i) == DIRICH ) {
                               const size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                               Dirichlet_constraints.insert( make_pair( position, ar(i) ) );
                            }
                        niter++;
                     }
                  }
                 break;
               default:
                 throw csmp::Exception( FATAL_ERROR,
                                       "PDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::EliminateEssentialConditions:",
                                       "Variable type not recognised by this method" );
        }
    } // end for
 
    cout <<"\n\nPDE_IntegratorExperimental<dim,COMPUTATION_DOMAIN>::EliminateEssentialConditions: eliminating "<< Dirichlet_constraints.size() <<" degrees of freedom.\n";
 
    // X. using the information about the Dirichlet constraints to do a row/ column elimination on the solution matrix
    csmp::SparseMatrix B;
    eliminateDirichletConstraints( Dirichlet_constraints, G_, rh_, B, Dirichlet_index_mapping_ );
    G_ = B;
   
 } // EliminateEssentialConditions
 
 
 





//
//  Clippings.cpp - code from Junchul Kim (2018) that was removed in latest CSMP revision
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 5/9/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//
/**

Explores all the nodes and the elements from the mesh by mesh traversal from the starting root nodes,
and returns the explored nodes and elements.

The method depends on correct neighbor information.

@param mesh  MeshMananger pointer
@param nodes  elmts: the method returns deques of pointers to the nodes and the elements

@section application Application

*/
template<size_t dim>
void exploreNodesAndElementsFromMesh( const MeshManager<dim>& mesh,
                                      deque<Node<dim>*>& nodes, deque<Element<dim>*>& elmts )
{
	// traversal of the existing mesh nodes to find all nodes and elements
	set<Element<dim>*>		explored_elements;
	set<Node<dim>*>	      discovered_nodes;
	deque<Node<dim>*>	    current_nodes;
	for (size_t g = 0U; g < mesh.NodeGroups(); g++) {
		auto root_node = mesh.RootNode(g);
		// starting at the root node
		discovered_nodes.insert(root_node);
		current_nodes.push_back(root_node);
		while (!current_nodes.empty()) {
			Node<dim>*  n_ptr(*current_nodes.begin());
			// for all parent elements of the current node
			for (size_t i = 0U; i < n_ptr->Parents(); i++) {
				// for all the nodes of each parent element
				for (size_t j = 0U; j < n_ptr->Parent(i)->Nodes(); j++)
					// if this node is not the one from which we started
					if (j != n_ptr->ParentNodeNumber(i)) {
						auto new_node = discovered_nodes.insert(n_ptr->Parent(i)->N(j));
						if (new_node.second) current_nodes.push_back(n_ptr->Parent(i)->N(j));
					}
				// storing the explored element
				explored_elements.insert(n_ptr->Parent(i));
			}
			// removing the node from the discovered (but not yet explored) deque
			current_nodes.pop_front();
		}
	}

	// assigning the explored nodes and elements from the node and element pointer vectors
	nodes.assign(discovered_nodes.begin(), discovered_nodes.end());
	elmts.assign(explored_elements.begin(), explored_elements.end());
}

template void exploreNodesAndElementsFromMesh( const MeshManager<1U>&, std::deque<Node<1U>*>&, std::deque<Element<1U>*>&);
template void exploreNodesAndElementsFromMesh( const MeshManager<2U>&, std::deque<Node<2U>*>&, std::deque<Element<2U>*>&);
template void exploreNodesAndElementsFromMesh( const MeshManager<3U>&, std::deque<Node<3U>*>&, std::deque<Element<3U>*>&);




/**

Explores all the faces from the mesh by mesh traversal from the starting root faces,
and returns the explored faces.

The method depends on correct neighbor information.

@param mesh  MeshMananger pointer
@param faces  the method returns deques of pointers to the faces

@section application Application

*/
template<size_t dim>
void exploreFacesFromMesh( const MeshManager<dim>& mesh, std::deque<Face<dim>*>& faces )
{
	// traversal of the existing mesh root faces to find all faces
	set<Face<dim>*>	discovered_faces;
	deque<Face<dim>*>	current_faces;

	for (size_t g = 0U; g < mesh.FaceGroups(); g++) {
		auto root_face = mesh.RootFace(g);
		// starting at the first face
		discovered_faces.insert(root_face);
		current_faces.push_back(root_face);
		while (!current_faces.empty()) {
			Face<dim>*  n_ptr(*current_faces.begin());
			// for all neighbor faces of the current face
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;

				// if this neighbor is new one
				auto new_face = discovered_faces.insert(n_ptr->Neighbor(i));
				if (new_face.second) current_faces.push_back(n_ptr->Neighbor(i));
			}
			// removing the face from the discovered (but not yet explored) deque
			current_faces.pop_front();
		}
	}
	
	// assigning the explored faces from the face pointer vectors
	faces.assign(discovered_faces.begin(), discovered_faces.end());
}

template void exploreFacesFromMesh( const MeshManager<1U>&, std::deque<Face<1U>*>&);
template void exploreFacesFromMesh( const MeshManager<2U>&, std::deque<Face<2U>*>&);
template void exploreFacesFromMesh( const MeshManager<3U>&, std::deque<Face<3U>*>&);




/**

Explores all the interfaces from the mesh by mesh traversal from the starting root interfaces,
and returns the explored interfaces.

The method depends on correct neighbor information.

@param mesh MeshMananger pointer
@param interfaces the method returns deques of pointers to the interfaces

@section application Application

*/
template<size_t dim>
void exploreInterFacesFromMesh( const MeshManager<dim>& mesh, std::deque<InterFace<dim>*>& interfaces )
{
	// traversal of the existing mesh root interfaces to find all interfaces
	set<InterFace<dim>*>	discovered_interfaces;
	deque<InterFace<dim>*>	current_interfaces;

	for (size_t g = 0U; g < mesh.InterFaceGroups(); g++) {
		auto root_interface = mesh.RootInterFace(g);
		// starting at the first interface
		discovered_interfaces.insert(root_interface);
		current_interfaces.push_back(root_interface);
		while (!current_interfaces.empty()) {
			InterFace<dim>*  n_ptr(*current_interfaces.begin());
			// for all neighbor interfaces of the current interface
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;

				// if this neighbor is new one
				auto new_interface = discovered_interfaces.insert(n_ptr->Neighbor(i));
				if (new_interface.second) current_interfaces.push_back(n_ptr->Neighbor(i));
			}
			// removing the interface from the discovered (but not yet explored) deque
			current_interfaces.pop_front();
		}
	}

	// assigning the explored interfaces from the interface pointer vectors
	interfaces.assign(discovered_interfaces.begin(), discovered_interfaces.end());
}

template void exploreInterFacesFromMesh( const MeshManager<1U>&, std::deque<InterFace<1U>*>& );
template void exploreInterFacesFromMesh( const MeshManager<2U>&, std::deque<InterFace<2U>*>& );
template void exploreInterFacesFromMesh( const MeshManager<3U>&, std::deque<InterFace<3U>*>& );


#endif /* Clippings */



