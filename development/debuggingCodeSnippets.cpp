//
//  debuggingCodeSnippets.cpp
//  Open CSMP++
//
// Code used to check things during the creation of OpenCSMP, grouped by functionality.
//
//  Created by Stephan Matthai on 25/4/2022.
//

#include "CSMP_definitions.h"

using namespace std;


// MESH MANAGER

#ifdef MESH_MANAGER_DEBUG
integrityCheck<dim,Element>( ElementsBegin(), ElementsEnd() );
if ( Faces() > 0 )
  integrityCheck<dim,Face>( FacesBegin(), FacesEnd() );
if ( InterFaces() > 0 ) {
      integrityCheck<dim,InterFace>( InterFacesBegin(), InterFacesEnd() );
     // add test for node manifolds
  }
#endif



    cout <<"\n\nprinting Face vector (face numbers & neighbors connected):\n";
    for ( const auto& fit : face_vector ) {
        cout <<"\n\tface "<< fit->Idx() <<": ";
        for ( auto i{0U}; i<fit->Neighbors(); i++ ) {
            if ( fit->Neighbor(i) == nullptr ) cout <<"none";
            else cout << fit->Neighbor(i)->Idx();
            cout <<", ";
          }
      }




// InterFace manifolds are assumed to be created before the creation of the InterFace by Duplicate(node)
   // 3. assigning node manifolds
   const size_t n_nodes{(*ifp).FE()->Nodes()};
   for ( auto i{0U}; i<n_nodes; ++i )
     // if the inside node is different from the outside node so that there needs to be a manifold
     if ( (*ifp).N(i,INSIDE) != (*ifp).N(i,OUTSIDE) ) {
         // 1. if both nodes are not yet manifolds
         if ( !(*ifp).N(i,INSIDE)->IsManifold() && !(*ifp).N(i,OUTSIDE)->IsManifold() ) {
              node_manifold_manager_->MergeManifolds( (*ifp).N(i,INSIDE)->Manifold(), (*ifp).N(i,OUTSIDE)->Manifold() );
              continue;
           }
         // 2. if both nodes are already manifolds, they are merged into single one
         if ( (*ifp).N(i,INSIDE)->IsManifold() && (*ifp).N(i,OUTSIDE)->IsManifold() ) {
              // if they are different from one-another, they are merged
              if ( (*ifp).N(i,INSIDE)->Manifold() != (*ifp).N(i,OUTSIDE)->Manifold() )
                node_manifold_manager_->MergeManifolds( (*ifp).N(i,INSIDE)->Manifold(),
                                                        (*ifp).N(i,OUTSIDE)->Manifold() );
              continue;
           }
         // 3. if the inside node is already a manifold and does not contain the second one
         //    because the second is not a manifold
         if ( (*ifp).N(i,INSIDE)->IsManifold() ) {
              // the outside node is added to it
              (*ifp).N(i,INSIDE)->Manifold()->Add( (*ifp).N(i,OUTSIDE), OUTSIDE );
           }
         // 4. if the outside node is already a manifold
         else if ( (*ifp).N(i,OUTSIDE)->IsManifold() ) {
              // the inside node is added to the outside nodes manifold
              (*ifp).N(i,OUTSIDE)->Manifold()->Add( (*ifp).N(i,INSIDE), INSIDE );
           }
       }





#ifdef MESH_MANAGER_DEBUG
integrityCheck<dim,Element>( ElementsBegin(), ElementsEnd() );
if ( Faces() > 0 )
  integrityCheck<dim,Face>( FacesBegin(), FacesEnd() );
if ( InterFaces() > 0 ) {
      integrityCheck<dim,InterFace>( InterFacesBegin(), InterFacesEnd() );
     // add test for node manifolds
  }
#endif




// getting the segment nodes from these face-nodes assuming that the first nodes are the corner nodes
switch( etype ) {
case ISOPARAMETRIC_LINEAR_TRIANGLE:
case ISOPARAMETRIC_QUADRATIC_TRIANGLE:
case LINEAR_TRIANGLE:
case LINEAR_TRIANGLE3D: {
      // segment 1, of the three segments given by the corner nodes of the triangular face
      auto sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N( fnids[1] ), it.first.first->N( fnids[2] ) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
      sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N(fnids[2]), it.first.first->N(fnids[0]) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
      sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N(fnids[0]), it.first.first->N(fnids[1]) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
    }
  break;
case ISOPARAMETRIC_LINEAR_QUADRILATERAL:
case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL:
case LINEAR_RECTANGLE: {
      auto sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N(0), it.first.first->N(1) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
      sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N(1), it.first.first->N(2) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
      sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N(2), it.first.first->N(3) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
      sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N(3), it.first.first->N(0) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
    }
  break;
case ISOPARAMETRIC_LINEAR_BAR:
case ISOPARAMETRIC_QUADRATIC_BAR:
case LINEAR_BAR: {
      // here the nodes at the endpoints are used
      auto sit = segment_nbors.insert( make_pair( set<Node<dim>*>{ it.first.first->N(0), it.first.first->N(1) }, set{iface_count} ) );
      if ( sit.second == false ) (*sit.first).second.insert( iface_count );
    }
  break;
default:
  csmp_error.notice( WARNING, "MeshManager<dim>::CreateInterFacesBetweenNodeSharingElements",
                     parseFiniteElementType(etype), "face type of element could not be parsed.");




/**
      Sets neighbor pointers of cells surrounding the domain to 'nullptr' if they were pointing to cells within the domain.
      
      The intention of this method is to avoid that these pointers wil accidentially be derefefenced causing crashes
      once the subdomain has been deleted.
      
      @todo: appears to have side effects.
      
      @todo IF POSSIBLE WE AVOID A DIRECT CONNECTION BETWEEN SUBDOMAIN AND MESHMANAGER
      
      SKM 9/2/2022
*/
template<uint32_t dim>
template<template<uint32_t> class CELL>
size_t MeshManager<dim>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<dim,CELL>& subdomain )
 {
    size_t n_detachments{0};
    const size_t n_cells{ subdomain.Cells() };
    for ( size_t i=subdomain.InteriorCells(); i < n_cells; ++i ) {
         const uint32_t n_perim_faces( subdomain.PerimeterFaces(i) );
         for ( auto j{0U}; j < n_perim_faces; ++j ) {
              auto p_face = subdomain.PerimeterFace( i, j );
              // detaching outside neighbor, if any
              if ( subdomain.E(i)->Neighbor(p_face) != nullptr ) {
                   const auto n_nbor_nbors{ subdomain.E(i)->Neighbor(p_face)->Neighbors() };
                   // loop over the outside neighbors faces until the one is found that matches the perimeter element
                   for ( auto k{0}; k<n_nbor_nbors; ++k )
                     if ( subdomain.E(i)->Neighbor(p_face)->Neighbor(k) == subdomain.E(i) ) {
                          // detach subdomain cell
                          subdomain.E(i)->Neighbor(p_face)->Neighbor(k)->Unassign( subdomain.E(i) );
                          n_detachments++;
                          break;
                       }
                }
           }
      }
      
   return n_detachments;
  
 } // end DetachOutsideNeighborsAlongPerimeter
  
template size_t MeshManager<1>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<1,Element>& );
template size_t MeshManager<2>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<2,Element>& );
template size_t MeshManager<3>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<3,Element>& );

template size_t MeshManager<1>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<1,Face>& );
template size_t MeshManager<2>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<2,Face>& );
template size_t MeshManager<3>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<3,Face>& );


 // TESTING
cerr <<"\nFace connectivity at boundary:\n";
while ( first2 != cellsEnd ) {
     cerr <<"\n\t"<< (*first2)->Idx() <<": ";
     for ( size_t face{0}; face < (*first2)->Faces(); ++face ) {
          if ( (*first2)->Neighbor(face) == nullptr ) cerr << face <<":NO ";
          else cerr << face <<":"<< (*first2)->Neighbor(face)->Idx() <<" ";
       }
     first2++;
  }


// DEBUGGING
cout <<"\n\n"<<"printing line element pairs (node, faces that share it):\n";
for ( auto it : cell_pairs ) {
     cout <<"\t\t"<< it.first->Idx() <<": ";
     for ( auto i : it.second ) cout << i.first->Idx() <<", ";
     cout << endl;
  }
  
  
  
  
// CLIPPING OF FACE NUMBERING FUNCTION FOR INTERFACE
    const size_t neighbors2x( f->Neighbors() * 2 );
    for ( size_t j{0U}; j<neighbors2x; ++j ) {
      InterFace<dim>* const ptr( f->Neighbor( j ) );
      if ( ptr != nullptr ) {
        // building search maps that we will use to find the shared interfaces
        // key=pointset   face iD
        map<set<Point<dim> >, pair<INTERFACE_SIDE, size_t> >   inner_elmt_faces, outer_elmt_faces;
        vector<uint32_t>  nids;
        // first element
        Element<dim>* e1 = f->InnerParent();
        for ( size_t face = 0U; face<e1->Faces(); ++face ) {
          e1->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j{0U}; j<nids.size(); ++j )
            face_key.insert( e1->N( nids[j] )->Coordinate() );
          outer_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
        }
        // second element
        Element<dim>* e2 = f->OuterParent();
        for ( size_t face = 0U; face<e2->Faces(); ++face ) {
          e2->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j{0U}; j<nids.size(); ++j )
            face_key.insert( e2->N( nids[j] )->Coordinate() );
          inner_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
        }

        // 2. finding the shared faces
        bool found( false );
        int64_t  inner_face_id( -1 ), outer_face_id( -1 );
        for ( auto& inner_face : inner_elmt_faces ) {
          for ( auto& outer_face : outer_elmt_faces ) {
            if ( inner_face.first == outer_face.first ) {
              inner_face_id = inner_face.second.second;
              outer_face_id = outer_face.second.second;
              found = true;
              break;
            }
          }
          if ( found ) break;
        }

        vset.Pfvert( eidx, j, static_cast<int32_t>(ptr->Idx()) );
      }
      else
        vset.Pfvert( eidx, j, REGION_BOUNDARY );
    }
// FACE NUMBERING



// CUT OUT BECAUSE BETTER DONE WITH HELP OF OCTREE
  /// only if there is not already a node at this location a new node is created, else a pointer to the nearest node is returned, no parent element  gets connected
  Node<dim>* const		 AddNodeAtUniqueLocation( const Point<dim>&, size_t nearby_node,
                                                const LocalVariables& node_variables,
                                                BOX_BOUNDARY = NOT,
                                                TOPOTYPE = MESH_VERTEX );

/**
    Inserts  a new Node at the desired point, but only if there is not already a Node there.
    
    @param nearby_node a node which the new node will be connected with.
    @param nvars the nodal variables this node will have to store
    @param bflag the boundary flag of the node (NOT if not at a boundary)
    @param topo the topological flag for the new node, i.e., whether it will be part of the models BREP
    @return if there is already a node at the point location, a pointer to that node is returned
    
    Search algorithm for a collocated node uses  "nearby" node as a starting point.
        
      Idea: start from nearby Node
        - loop over the neighbor nodes of the node ranking them in terms of their proximity from the target point
        - move to closest node and then repeat (remembering the shortest distance)
        - repeat until node is found while the distance decreases
        - if distance increases, the node does not exist and will be created
        - allow  to move across manifold member nodes in order to cross split boundaries
*/
template<uint32_t dim>
Node<dim>* const MeshManager<dim>::AddNodeAtUniqueLocation( const Point<dim>& pt,
                                                            size_t nearby_node,
                                                            const LocalVariables& nvars,
                                                            BOX_BOUNDARY bflag,
                                                            TOPOTYPE topo )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
   // if the node location needs to be compared with existing ndes
   if ( nearby_node >=nodes_.size() ) {
        csmp_error.Note( WARNING, "MeshManager<dim>::AddNodeAt",
                          "nearby Node not contained in Mesh:", to_string(nearby_node) );
        // using the last node
        nearby_node = nodes_.size() - 1U;
     }

   // searching the mesh tree for a node with the same location (using the provided point as a start location)
   Node<dim>*          nptr( &(*next(nodes_.begin(),static_cast<long>(nearby_node))) );
   double              new_distance(pt.DistanceTo(nptr->Coordinate())), old_distance(1e30);
   map<double,size_t>  distances;
   // estimating a tolerance on the basis of the distance of the point to the node and the first node
   const double tolerance = 1.0e-7 * (new_distance + pt.DistanceTo((*nodes_.begin()).Coordinate())) / 2.;
   while ( old_distance > new_distance )
     {
        // tree traversal: looping the neighbor nodes of the current node, finding the one that is the closest to the point
        const size_t n_nbors( nptr->Neighbors() );
        for ( uint32_t i{0U}; i<n_nbors; ++i )
          distances.insert( make_pair( pt.DistanceTo( nptr->Neighbor(i)->Coordinate() ), i ) );
        // since map defaults to less, its first entry is the node we want
        nptr = nptr->Neighbor( static_cast<uint32_t>((*distances.begin()).second) );
        old_distance = new_distance;
        new_distance = (*distances.begin()).first;
        assert( nptr != nullptr );
     }
   // TODO: deal with NodeManifolds - if IsManifold()...
   // if a node matching the point location was found, a pointer to it is returned
   if ( fabs(new_distance) < tolerance ) return nptr;
  
   // else a new node is created
   typename plf::colony<Node <dim>>::iterator
     nit = nodes_.emplace( Node<dim>( nodes_.size(), pt, nvars, bflag, topo ) );
     
   return &(*nit);
}


    // converting lower-dimensional regions to boundaries if they consist of two-dimensional elements located on the model boundary
    // - finding potential 2-dimensional regions
    {
      set<string> twoDimRegions;
      for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit ) {
           const auto spatialDim = (*rit).second.ElementSpatialDimensions();
           assert( spatialDim.first == 1U ); // making sure that all elements in the region have the same spatial dimension
           if ( spatialDim.second == 2U )    // highest spatial dimension of elements in region
             twoDimRegions.insert( (*rit).first );
        }
      // converting these vertical side boundary faces into a single boundary
      if ( !twoDimRegions.empty() ) {
           model.MergeRegions( twoDimRegions, "SIDES" );
           // setting all boundary nodes to IRREGULAR as is needed by CreateExternalBoundaryFrom()
           Region<3U>& sides = model.Region("SIDES");
           for ( auto& nit : sides.NodeVector() ) nit->AtBoundary( IRREGULAR );
           const bool check_topo_attributes_of_nodes{false};
           model.CreateExternalBoundaryFrom( "SIDES", check_topo_attributes_of_nodes );
           // the input region gets removed automatically
        }
    }

