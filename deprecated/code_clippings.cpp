
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



/**
higherDimensionalNeighbors() - finds the higher-dim neighbor elements of
a dim-1 element embedded within the higher-dim mesh.

one inside element idx or two inside-outside element idxs are stored in the parameter 'in_out_elements'.

assumptions
- assumes that the nodes and elements in the entire model domain are numbered continuously

@author Eduardo Pezzuli
@date 2019

*/
template<uint32_t dim>
bool  higherDimensionalNeighbors( const Element<dim>& e, vector<Element<dim>*>& in_out_elements )
{
  // 1. looping over the parent elements of the nodes searching for the faces which are shared with the lower dimensional element
  // -----------------------------------------------------------------------------------------------------------------------------
  // making a set of element nodes to later identify faces by comparison
  set<size_t>  node_set, test_set;
  const uint32_t nodes( e.Nodes() );
  for ( auto i{0U}; i<nodes; ++i ) node_set.insert( e.N( i )->Idx() );
  map<Element<dim>*,uint32_t>  nbor_elmts;
  vector<uint32_t> fnids;
  for ( auto i{0U}; i<nodes; i++ ) {
    const auto parents( e.N( i )->Parents() );
    for ( auto j{0U}; j<parents; ++j ) {
      Element<dim>* eptr( e.N( i )->Parent( j ) );
      const auto faces( eptr->Faces() );
      for ( auto k{0U}; k<faces; ++k ) {
        eptr->FE()->NodesOfFace( k, fnids );
        auto fnodes( fnids.size() );
        for ( auto l = 0U; l<fnodes; ++l )
          test_set.insert( eptr->N( fnids[l] )->Idx() );
        // if the face is shared the element and its face are recorded
        if ( node_set == test_set ) {
          // storing a pointer to this element and its local face number
          // making sure that no duplicate is received
          nbor_elmts.insert( make_pair( eptr, k ) );
        }
        test_set.clear();
      }
    }
  }

  if( nbor_elmts.size() == 1U ){
    typename map<Element<dim>*,uint32_t>::const_iterator  nbit( nbor_elmts.begin() );
    assert( (*nbit).first != nullptr );
    in_out_elements.push_back( (*nbit).first );
  }
  else if ( nbor_elmts.size() == 2U ) {
    // 2. finding which of the neighbors is the inside one by projecting face normals onto lower dim element normal
    // -------------------------------------------------------------------------------------------------------------
    vector<double>  enrml, fnrml;
    e.UnitNormal( enrml );
    typename map<Element<dim>*,uint32_t>::const_iterator  nbit( nbor_elmts.begin() );
    bool inside_elmt_found( false );
    bool outside_elmt_found( false );
    pair<Element<dim>*, Element<dim>*> nbors;
    pair<uint32_t,uint32_t> faces;

    // first element
    // -------------
    assert( (*nbit).first != nullptr );
    faces.first = (*nbit).second;
    (*nbit).first->UnitNormalToFace( faces.first, fnrml );
    double dotproduct( 0. );
    for ( size_t k = 0U; k < dim; ++k )
      dotproduct += enrml[k] * fnrml[k];

    // if the projection is negative, the first element lies on the outside
    if ( dotproduct < 0. ) {
      nbors.second = (*nbit).first;
    }
    else {
      nbors.first = (*nbit).first;
    }
    nbit++;
    
    // second element
    // --------------
    assert( (*nbit).first != nullptr );
    faces.second = (*nbit).second;
    (*nbit).first->UnitNormalToFace( faces.second, fnrml );
    dotproduct = 0.;
    for ( auto k{0U}; k < dim; ++k )
      dotproduct += enrml[k] * fnrml[k];

    // if the projection is negative, the second element lies on the outside
    if ( dotproduct < 0. ) {
      // checking that we have no duplication here
      assert( outside_elmt_found == false );
      nbors.second = (*nbit).first;
    }
    else {
      assert( inside_elmt_found == false );
      nbors.first = (*nbit).first;
    }

    in_out_elements.push_back( nbors.first );
    in_out_elements.push_back( nbors.second );
    
  }
  else{
    return false;
  }

  return true;
} // end higherDimensionalNeighbors






/**
     higherDimensionalNeighbor() - finds a higher-dimensional element, one face of which
     matches  the input supplied lower-dimensional element.
     
     @return pointer to the higher-dimensional adjacent element or NULL when not found.
     
     @return face of the higher-dimensional element that matches the lower-dim element
     
     @return material ID of the higher-dimensional element; if it does not exist, NaN is returned.
     @note an undefined material key is accepted, but in this case NaN will be returned.
 
     The discovered element is considered to be located on the inside of the lower-dim element,
     hence its normal ought to be pointing away from it, else there is an orientation problem.
 
     @attention assumptions
     - assumes that the nodes and elements in the entire model domain are numbered continuously
     
     application
     - use this function for finding the higher-dimensional neighbor of a surface element that sits on the
       outside boundary of the model
 
     @test SKM 22/8/2018 - fixed a bug where element returned had lower spatial dimensional than supplied
     element.
 
*/
template<uint32_t dim>
const Element<dim>* const higherDimensionalNeighbor( const Element<dim>& e, const csmp::Index& mtrl_key,
                                                     uint32_t& local_face_number_of_e, double& material_ID  )
 {
     if constexpr ( dim == 3 ) assert( e.IsSurfaceElement() );
     if constexpr ( dim == 2 ) assert( e.IsLineElement() );
     assert( mtrl_key.type == SCALAR );
     assert( mtrl_key.place == ELEMENT || mtrl_key.place == UNDEFINED );

     // 1. looping over the parent elements of the nodes searching their faces for ones that are shared with the lower dimensional element
     // ----------------------------------------------------------------------------------------------------------------------------------
     // making a set of element nodes to later identify faces by comparison
     set<size_t>   node_set, test_set;
   
     const size_t  nodes(e.Nodes());
     for ( uint32_t i{0}; i<nodes; ++i ) node_set.insert(e.N(i)->Idx());
   
     const csmp::Element<dim>*  nbor_elmt(nullptr);
     vector<uint32_t>           fnids;
   
     // since the same element may be discovered by each of the face nodes
     // the loop is stopped after the first discovery
     for ( auto i{0U}; i<nodes; i++ ) {
          const auto parents(e.N(i)->Parents());
          for ( auto j{0U}; j<parents; ++j ) {
               const Element<dim>* const eptr(e.N(i)->Parent(j));
               // only if the element is not the same and also of a different type
               if ( eptr != &e and
                    eptr->FE_Type() != e.FE_Type() and
                    eptr->Nodes() >= e.Nodes() )
                 {
                   const auto faces(eptr->Faces());
                   for ( auto k{0U}; k<faces; ++k ) {
                         eptr->FE()->NodesOfFace( k, fnids );
                         size_t fnodes(fnids.size());
                         for ( size_t l{0U}; l<fnodes; ++l )
                           test_set.insert( eptr->N( fnids[l] )->Idx() );
                         // if the face is shared the element and its face are recorded
                         if ( node_set == test_set ) {
                              // storing the pointer to this element and its local face number
                              // making sure that no duplicates are received
                              nbor_elmt = eptr;
                              local_face_number_of_e = k;
                              break;
                           }
                         test_set.clear();
                     }
                 }
            }
       }
    assert( nbor_elmt != nullptr );

     if constexpr ( dim == 3 ) {
           if (e.IsVolumeElement())
             assert( nbor_elmt->IsVolumeElement() );
           else if (e.IsSurfaceElement())
             assert( nbor_elmt->IsVolumeElement() );
           else if (e.IsLineElement())
             assert( nbor_elmt->IsSurfaceElement() or  nbor_elmt->IsVolumeElement() );
       }
     if constexpr ( dim == 2 ) {
           if (e.IsSurfaceElement())
             assert( nbor_elmt->IsSurfaceElement() );
           else if (e.IsLineElement())
             assert( nbor_elmt->IsSurfaceElement() );
       }
   
    // 2. drawing the results
    // -------------------------------------------------------------------------------------------------------------
    material_ID = ( nbor_elmt != nullptr && mtrl_key.place != UNDEFINED )
                  ? nbor_elmt->Read(mtrl_key) : numeric_limits<double>::quiet_NaN();
   
    return move(nbor_elmt);
   
 } // end higherDimensionalNeighbor


// STUB
template<>
const Element<1U>* const higherDimensionalNeighbor( const Element<1U>& e, const csmp::Index&, uint32_t&, double& )
 {
    throw logic_error("higherDimensionalNeighbor(in BoundaryInterface: there should not be any boundaries in a 1D model.");
    return &e;
 }


template const Element<2U>* const higherDimensionalNeighbor( const Element<2U>&, const csmp::Index&, uint32_t&, double& );
template const Element<3U>* const higherDimensionalNeighbor( const Element<3U>&, const csmp::Index&, uint32_t&, double& );

// TESTING
/*
if ( dim == 2 && nbor_elmt->IsLineElement() ) {
     cerr <<"\nhigherDimensionalNeighbor: potentially found duplicate edge elements; candidates are:\n";
     e.Out();
     cerr <<"\n\nand:\n";
     nbor_elmt->Out();
  }
*/



/**
    Connects Element objects to their same-dimensional neighbors in as much as is possible.
    
    Where there are no neighbors the neighbor pointers will be nulled.
    
    @attention The assumption is made that all nodes in the model have a unique numbering.
    
    @author SKM 2012
    
        @todo deal with manifolds, disambiguating them on the basis of element orientation (only elements int the same plane or aligned elements should be neighbors)

        @TODO: no need to use set, use sort() then unique() on the result vector, searching will be much faster
*/
template<uint32_t dim>
void  establishNeighborConnectivity( vector<Element<dim>*>& simplexVector, bool unassign_neighbors_outside, bool verbose )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( simplexVector.empty() ) {
         csmp_error.notice( WARNING, "establishNeighborConnectivity( Element )", "supplied element vector is empty; nothing was done" );
         return;
      }

    cout << "\nestablishNeighborConnectivity( Element ): Establishing CSMP FE neighbor connectivity...\n";
 
    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    if (verbose) cout << "  Building element face list...\n";

    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<uint32_t,Element<dim>*> >  volume_neighbor_keys,
                                                             surface_neighbor_keys,
                                                             line_neighbor_keys;
    vector<uint32_t>               fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<Element<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( auto face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == nullptr ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( Element )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }
           // creating face key from idx's of face
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( auto j{0U}; j<fnids.size(); j++ )
               key.insert( (*it)->N( fnids[j] ) );
             
           // inserting newly generated keys into multimap
           if ( (*it)->FE()->IsVolumeElement() )
               volume_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else if ( (*it)->FE()->IsSurfaceElement() )
               surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else // for all line elements
               line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           key.clear();

           // unassign neighbors outside of the provided vector range
           if ( unassign_neighbors_outside )
           {
             // remove element from neighbor list of its neighbors
             const auto neighbors( (*it)->Neighbors() );
             for ( auto neighbor = 0; neighbor < neighbors; ++neighbor )
               if ( (*it)->Neighbor( neighbor ) != NULL ) {
                 const auto neighbor_neighbors( (*it)->Neighbor( neighbor )->Neighbors() );
                 for ( auto i = 0; i < neighbor_neighbors; i++ )
                   if ( (*it)->Neighbor( neighbor )->Neighbor( i ) != NULL )
                     if ( (*it)->Neighbor( neighbor )->Neighbor( i ) == (*it) ) {
                       (*it)->NeighborElementVector()[i] = NULL;
                     }
                 
                 (*it)->NeighborElementVector()[neighbor] = NULL;
               }
             (*it)->NeighborElementVector().clear();
           }
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)

    if (verbose) cout << "  Building element neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(nullptr);
        Element<dim>* e2Ptr(nullptr);

        if (verbose) cout << "\n\t\tline elements...";

        auto it1(line_neighbor_keys.begin()), it2(line_neighbor_keys.begin());

        it2++;

        while ( it2 != line_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
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
    if ( dim >= 2U and !surface_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(nullptr);
        Element<dim>* e2Ptr(nullptr);

        if (verbose) cout << "\n\t\tsurface elements...";

        auto it1(surface_neighbor_keys.begin()), it2(surface_neighbor_keys.begin());

        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                    // assigning eachothers faces
                    //                          face pointer                   nbor face idx        neighbor pointer
                    ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                    ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // surface elements
      
    // 2.3 volume elements
    // -------------------
    if ( dim == 3U and !volume_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(nullptr);
        Element<dim>* e2Ptr(nullptr);

        if (verbose) cout << "\n\t\tvolume elements...\n";

        auto it1(volume_neighbor_keys.begin()), it2(volume_neighbor_keys.begin());

        it2++;

        while ( it2 != volume_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ( (*it1).second.second != 0 and (*it2).second.second != 0 ) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == volume_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=3
    
 } // end establishNeighborConnectivity

// explicit instantiations
template void establishNeighborConnectivity<1U>( std::vector<csmp::Element<1U>*>&, bool, bool );
template void establishNeighborConnectivity<2U>( std::vector<csmp::Element<2U>*>&, bool, bool );
template void establishNeighborConnectivity<3U>( std::vector<csmp::Element<3U>*>&, bool, bool );




/** CONNECTIVITY BETWEEN INTERFACES
       
         Inside and outside must be considered.
*/
template<uint32_t dim>
void  establishNeighborConnectivity( std::vector<csmp::InterFace<dim>*>& simplexVector, bool unassign_neighbors_outside, bool verbose )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( simplexVector.empty() ) {
         csmp_error.notice( WARNING, "establishNeighborConnectivity( InterFace ):", "supplied element vector is empty; nothing was done." );
         return;
      }

    if (verbose) cout << "\nestablishNeighborConnectivity( InterFace ): Establishing CSMP FE neighbor connectivity...\n";

    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    if (verbose)  cout << "  Building element face list...\n";

    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<uint32_t,InterFace<dim>*> >  surface_neighbor_keys,
                                                               line_neighbor_keys;
    vector<uint32_t>               fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<InterFace<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( auto face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == NULL ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( InterFace )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }

           // creating face keys from idx's of interface for INSIDE & OUTSIDE
           (*it)->CurrentSide( INSIDE ); // just for node-vector
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( size_t j{0U}; j<fnids.size(); j++ )
             key.insert( (*it)->N( fnids[j] ) );

           // inserting newly generated key into multimap
           if ( (*it)->FE()->IsSurfaceElement() )
               surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else if ( (*it)->FE()->IsLineElement() )
               line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else{
               csmp_error.notice( ERROR, "establishNeighborConnectivity( InterFace )", "supplied element vector contains volumetric element! nothing was done" );
               return;
           }
           key.clear();

           // unassign neirghbors outside of the provided vector range
           if ( unassign_neighbors_outside )
           {
             // remove element from neighbor list of its neighbors
             const auto neighbors( (*it)->Neighbors() );
             for ( auto neighbor = 0; neighbor < neighbors; ++neighbor )
               if ( (*it)->Neighbor( neighbor ) != NULL ) {
                 const auto neighbor_neighbors( (*it)->Neighbor( neighbor )->Neighbors() );
                 for ( auto i = 0; i < neighbor_neighbors; i++ )
                   if ( (*it)->Neighbor( neighbor )->Neighbor( i ) != NULL )
                     if ( (*it)->Neighbor( neighbor )->Neighbor( i ) == (*it) ) {
                       (*it)->NeighborElementVector()[i] = NULL;
                     }

                 (*it)->NeighborElementVector()[neighbor] = NULL;
               }
             (*it)->NeighborElementVector().clear();
           }
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    if (verbose) cout << "  Building element neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() ) {

        InterFace<dim>* e1Ptr(NULL);
        InterFace<dim>* e2Ptr(NULL);

        if (verbose) cout << "\n\t\tline elements...";

        auto it1(line_neighbor_keys.begin()), it2(line_neighbor_keys.begin());

        it2++;

        while ( it2 != line_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                    //                        face pointer, neighbor pointer, side-of interface
                   (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                   (*it2).second.second->Assign( (*it2).second.first, e1Ptr );

                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=1

    // 2.2 surface elements
    // --------------------
    if ( dim >= 2U and !surface_neighbor_keys.empty() ) {

        InterFace<dim>* e1Ptr(NULL);
        InterFace<dim>* e2Ptr(NULL);

        if (verbose) cout << "\n\t\tsurface elements...";

        auto it1(surface_neighbor_keys.begin()), it2(surface_neighbor_keys.begin());

        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                    // assigning eachothers faces
                    //                        face pointer, neighbor pointer, side-of interface
                    (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                    (*it2).second.second->Assign( (*it2).second.first, e1Ptr );


                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=2


 } // end establishNeighborConnectivity

template void establishNeighborConnectivity<1U>( std::vector<csmp::InterFace<1U>*>&, bool, bool );
template void establishNeighborConnectivity<2U>( std::vector<csmp::InterFace<2U>*>&, bool, bool );
template void establishNeighborConnectivity<3U>( std::vector<csmp::InterFace<3U>*>&, bool, bool );





#endif /* Clippings */



