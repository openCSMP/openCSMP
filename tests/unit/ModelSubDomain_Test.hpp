//
//  ModelSubDomain_Test.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#ifndef ModelSubDomain_Test_hpp
#define ModelSubDomain_Test_hpp

#include "Test.h"
#include "ModelSubDomain.h"

namespace csmp {

class ModelSubDomain_Test : public Test {
  public:
    virtual void run();

    /// compares node locations and connectivity
    template<size_t dim,template<size_t> class simplicial_complex>
    bool CompareModelSubdomains( const ModelSubDomain<dim,simplicial_complex>&,
                                 const ModelSubDomain<dim,simplicial_complex>&,
                                 bool verbose );
};

/** compares node locations and connectivity
    - uses point locations for node comparison

*/
template<size_t dim,template<size_t> class simplicial_complex>
bool ModelSubDomain_Test::CompareModelSubdomains( const ModelSubDomain<dim,simplicial_complex>& domain1,
                                                  const ModelSubDomain<dim,simplicial_complex>& domain2,
                                                  bool verbose )
 {
    // domain names
    if ( domain1.Name() != domain2.Name() ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: names mismatch: "<< domain1.Name() <<" vs. "<< domain2.Name() <<"\n";
         //return false;
      }
    _test( domain1.Name() == domain2.Name() );
   
    // 1. nodes
    // --------
    // number of nodes
    if ( domain1.Nodes() != domain2.Nodes() ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: node numbers do not match: "<< domain1.Nodes() <<" vs "<< domain2.Nodes() <<"\n";
         //return false;
      }
    _test( domain1.Nodes() == domain2.Nodes() );
   
    // interior nodes
    if ( domain1.InteriorNodes() != domain2.InteriorNodes() ) {
         if ( verbose ) {
              std::cerr <<"\ncompareModelSubdomains: number of interior nodes does not match: "<< domain1.InteriorNodes() <<" vs "<< domain2.InteriorNodes() <<"\n";
              std::cerr <<"\tinterior nodes domain 1 vs domain 2:\n";
              for ( auto nit=domain1.NodesBegin(); nit!=domain1.PerimeterNodesBegin(); ++nit )
                std::cerr << (*nit)->Idx() <<" ";
              std::cerr << std::endl;
              for ( auto nit=domain2.NodesBegin(); nit!=domain2.PerimeterNodesBegin(); ++nit )
                std::cerr << (*nit)->Idx() <<" ";
              std::cerr << std::endl;
           }
         //return false;
      }
    _test( domain1.InteriorNodes() == domain2.InteriorNodes() );
   
    // node boundary flags
    int node_flag_mismatches(0);
    auto nit2=domain2.NodesBegin();
    for ( auto nit1=domain1.NodesBegin(); nit1!=domain1.NodesEnd(); ++nit1, ++nit2 )
      if ( (*nit1)->AtBoundary() != (*nit2)->AtBoundary() ) node_flag_mismatches++;
   
    if ( node_flag_mismatches > 0 ) {
          std::cerr <<"\ncompareModelSubdomains: the BOX_BOUNDARY flags of the domains do not match.\n";
          //return false;
      }
    _test( node_flag_mismatches == 0 );
  
    // node coordinates (sorted ranges of points have to be created)
    // ranges for domain 1
    std::set<Point<dim> >  interior_points1;
    for ( auto nit=domain1.NodesBegin(); nit!=domain1.PerimeterNodesBegin(); ++nit )
      interior_points1.insert( (*nit)->Coordinate() );
    std::set<Point<dim> >  perimeter_points1;
    for ( auto nit=domain1.PerimeterNodesBegin(); nit!=domain1.NodesEnd(); ++nit )
      interior_points1.insert( (*nit)->Coordinate() );
    // ranges for domain 2
    std::set<Point<dim> >  interior_points2;
    for ( auto nit=domain2.NodesBegin(); nit!=domain2.PerimeterNodesBegin(); ++nit )
      interior_points2.insert( (*nit)->Coordinate() );
    std::set<Point<dim> >  perimeter_points2;
    for ( auto nit=domain2.PerimeterNodesBegin(); nit!=domain2.NodesEnd(); ++nit )
      interior_points2.insert( (*nit)->Coordinate() );
    // comparisons
    if ( interior_points1 != interior_points2 ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: locations of interior points do not match.\n";
         //return false;
      }
    _test( interior_points1 == interior_points2 );

    if ( perimeter_points1 != perimeter_points2 ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: locations of perimeter points do not match.\n";
         //return false;
      }
    _test( perimeter_points1 == perimeter_points2 );
   
    interior_points1.clear();
    interior_points2.clear();
    perimeter_points1.clear();
    perimeter_points2.clear();
   
    // if the points are all the same,
    // nodes and elements in both domains are renumbered to be able to compare connectivity
    domain1.RenumberNodes();
    domain2.RenumberNodes();
 
   
    // 2. elements, i.e. connectivity
    // ------------------------------
    // number of elements
    if ( domain1.Elements() != domain2.Elements() ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: number of elements does not match.\n";
         //return false;
      }
    _test( domain1.Elements() == domain2.Elements() );
   
    // number of interior elements
    if ( domain1.PerimeterElements() != domain2.PerimeterElements() ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: mismatch in number of perimeter elements.\n";
         //return false;
      }
    _test( domain1.PerimeterElements() == domain2.PerimeterElements() );
   
    // element connectivity of interior elements (not assuming that elements are in same order)
    std::set<std::vector<size_t> > plist_entries1;
    for ( auto it=domain1.ElementsBegin(); it!=domain1.PerimeterElementsBegin(); ++it ) {
         std::vector<size_t> nodes( (*it)->Nodes() );
         for ( size_t i=0U; i<(*it)->Nodes(); ++i ) {
              nodes[i] = (*it)->N(i)->Idx();
           }
         plist_entries1.insert( move(nodes) );
      }
    std::set<std::vector<size_t> > plist_entries2;
    for ( auto it=domain2.ElementsBegin(); it!=domain2.PerimeterElementsBegin(); ++it ) {
         std::vector<size_t> nodes( (*it)->Nodes() );
         for ( size_t i=0U; i<(*it)->Nodes(); ++i ) {
              nodes[i] = (*it)->N(i)->Idx();
           }
         plist_entries2.insert( move(nodes) );
      }
    // comparing plists
    bool plists_are_the_same(true);
    if ( plist_entries1 != plist_entries2 ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: (plist) mismatch in element-to-node connectivity.\n";
         plists_are_the_same = false;
         //return false;
      }
    _test( plists_are_the_same );

    plist_entries1.clear();
    plist_entries2.clear();
   
    // comparing the element neighbor connectivity
    domain1.RenumberElements();
    domain2.RenumberElements();

    std::set<std::vector<size_t> > pfverts_entries1;
    for ( auto it=domain1.ElementsBegin(); it!=domain1.PerimeterElementsBegin(); ++it ) {
         std::vector<size_t> nbors( (*it)->Neighbors(),0 );
         for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr ) {
                nbors[i] = (*it)->Neighbor(i)->Idx();
             }
         pfverts_entries1.insert( move(nbors) );
      }
    std::set<std::vector<size_t> > pfverts_entries2;
    for ( auto it=domain2.ElementsBegin(); it!=domain2.PerimeterElementsBegin(); ++it ) {
         std::vector<size_t> nbors( (*it)->Neighbors(),0 );
         for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr ) {
                nbors[i] = (*it)->Neighbor(i)->Idx();
             }
         pfverts_entries2.insert( move(nbors) );
      }
    // comparing neighbor connectivity lists
    bool pfverts_are_the_same(true);
    if ( pfverts_entries1 != pfverts_entries2 ) {
         if ( verbose )
           std::cerr <<"\ncompareModelSubdomains: (pfverts) mismatch in element to neighbor-element connectivity.\n";
         pfverts_are_the_same = false;
         //return false;
      }
    _test( pfverts_are_the_same );

    plist_entries1.clear();
    plist_entries2.clear();
   
    // all comparisons passed - domains have same nodes and elements.
    return true;
   
 } // end compareModelSubdomains


} // end csmp

#endif /* ModelSubDomain_Test_hpp */
