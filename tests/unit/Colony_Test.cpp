#include "Colony_Test.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "IsoparametricLinearTriangle.h"
#include "Element.h"
#include "Node.h"
#include "plf_colony.h"

using namespace std;

namespace csmp {

Colony_Test::Colony_Test()
{
}
	
  
  
Colony_Test::~Colony_Test()
{
}
	
  
  
/**   
    Testing:
 
      1. LengthInDirection() - quadrilateral 2D
      2. LengthInDirection() - quadrilateral 3D
 
      3. UnitNormalToFace() - actually testing FEPolicy functionality
 
*/
void Colony_Test::run()
{
  // verifying that pointers to objects in colony remaining after erase are not invalidated
  TestColonyWith_int();
  
  // is Element fit for storage in colony
  TestThatElementIsCopyableAssignableMovable();
  
  // testing standard operations with csmp::Element stored in colony
  TestColonyWith_Element();
}
	
  
  
bool Colony_Test::TestColonyWith_int()
{
   // COLONY: testing the use of a colony, with pointers (as opposed to iterators) to its elements
   // elements are deleted using the get_iterator() from pointer functionality
   plf::colony<int>  element_container{ 7, 3, 4, 9, 2, 1 };
   vector<int*>      int_ptrs;
   int_ptrs.reserve(6);
   for ( auto it=element_container.begin(); it!=element_container.end(); ++it )
     int_ptrs.push_back( &(*it) );
   
   cout <<"\nint colony:\n";
   for ( const auto& n : element_container ) cout <<" "<< n;
   cout << endl;

   cout <<"\ndereferenced int colony pointers:\n";
   for ( const auto& n : int_ptrs ) cout <<" "<< *n;
   cout << endl;

   // deleting specific elements
   element_container.erase( element_container.get_iterator( int_ptrs[2] ) );
   int_ptrs[2] = nullptr; // essential for pointer deletion from access vector to work
   element_container.erase( element_container.get_iterator( int_ptrs[4] ) );
   int_ptrs[4] = nullptr;
   // and the pointers pointing at them
   int_ptrs.erase( remove( int_ptrs.begin(), int_ptrs.end(), nullptr ), int_ptrs.end() );

   cout <<"\ndereferenced int colony pointers after erasure of elements 3 and 5:\n";
   for ( auto it=element_container.begin(); it!=element_container.end(); ++it ) {
        if ( (*it) != 4 && (*it) != 2 )
          // testing that only the targeted elements were erased
          _test( (*it) >= 1 && (*it) <= 9 );
        cout <<" "<< (*it);
     }
   cout << endl;

   cout <<"\ndereferenced colony pointers:\n"; // non-deleted elements are still valid
   for ( size_t n{0}; n<int_ptrs.size(); ++n ) {
         cout <<" "<< *int_ptrs[n];
         _test( *int_ptrs[n] == 7 || *int_ptrs[n] == 3 || *int_ptrs[n] == 9 || *int_ptrs[n] == 1 );
     }
   cout << endl;
   // int_ptrs.erase( next(int_ptrs.begin(),2) );
   // int_ptrs.erase( next(int_ptrs.begin(),4) );
   
   // adding and deleting a range of elements
   int_ptrs.reserve(25);
   for ( int i{10}; i<25; ++i ) {
        auto it = element_container.insert(i);
        int_ptrs.push_back( &(*it) );
     }
   cout <<"\nmain: colony grown by 15 elements:\n";
   int sum{0};
   auto int_it = element_container.begin();
   for ( size_t i{0ul}; i<element_container.size(); ++i ) {
        cout <<" "<< *int_ptrs[i] <<" prt->"<< (*int_it) <<" ";
        sum += (*int_it);
        int_it++;
     }
   cout << endl;
   
   // deleting the first 5 new elements via the pointers, then updating the pointer vector
   // (to find an element using next() and its original offset will stop working after the first insertion)
   for ( auto it=next(int_ptrs.begin(),4); it!=next(int_ptrs.begin(),9); ++it ) {
   // does not work because on cannot rely on the order of the elements in the container
   //  element_container.erase( element_container.get_iterator( (*next(int_ptrs.begin(),4)) ),
   //                           element_container.get_iterator( (*next(int_ptrs.begin(),9)) ) );
        element_container.erase( element_container.get_iterator( (*it) ) );
        (*it) = nullptr;
     }
   int_ptrs.erase( remove( int_ptrs.begin(), int_ptrs.end(), nullptr ), int_ptrs.end() );

   // test: pointers should still be valid even after insertions and erasures
   int checksum{0};
   cout <<"\nmain: colony after deleting 4 to 9:\n";
   for ( auto it=element_container.begin(); it!=element_container.end(); ++it ) {
        cout <<" "<< (*it);
        checksum += (*it);
     }
   cout << endl;
   _test( checksum == sum - 60 );

   cout <<"\nmain: colony access via pointers:\n";
   for ( auto it=int_ptrs.begin(); it!=int_ptrs.end(); ++it )
     cout <<" "<< *(*it);
   cout << endl;

   // more than doubling the size of the colony
   for ( int i{0}; i<100; ++i )
     element_container.insert( 99 );

   cout <<"\nmain: colony access via pointers after extra insertions into colony:\n";
   for ( auto it=int_ptrs.begin(); it!=int_ptrs.end(); ++it )
     cout <<" "<< *(*it);
   cout << endl;
   
   cout <<"\nmain: listing the elements of the colony using the next() iterator:\n";
   const size_t n_elements{ element_container.size() };
   for ( size_t i{0}; i<n_elements; ++i )
     cout <<" "<< (*next( element_container.begin(), i ));
   cout << endl;
   
   cout <<"\nmain: memory used by colony vs. stored data only: ";
   cout << element_container.memory() <<" vs "<< sizeof(int) * element_container.size() << endl;

   // demonstrating the small overhead for larger objects
   cout <<"\nmain: memory used by colony with one hundred 3.2 kbyte objects: ";
   struct Blob { double ary[500]; };
   plf::colony<Blob> blob_colony;
   blob_colony.insert( 500, Blob{} );
   cout << blob_colony.memory() <<" vs "<< sizeof(Blob) * blob_colony.size() << endl;
   
   return true;
   
} // end TestColonyWith_int




/**
  Check that Element meets the pre-requisites for storage in the colony.
*/
void Colony_Test::TestThatElementIsCopyableAssignableMovable()
 {
    // 1. creating one triangle and one quadrilateral
    const LocalVariables evars( 1, // scalarsVars,
                                2, // vectorVars,
                                1, // tensorVars,
                                0, // array_count,
                                0, // array_length,
                                0, // flag_array_count,
                                0, // flag_array_length,
                                9, // total_data_depth,
                                6 ); // total_flag_depth
                                
    const IntegrationPointVariables ivars;

  // Element<dim>( elmt_idx, fem, stencil, evars, cvars, mtrl_idx );
    FiniteVolumeStencil<2>  quad( "ISOPARAMETRIC_LINEAR_QUADRILATERAL" );
    FiniteVolumeStencil<2>  tria( "ISOPARAMETRIC_LINEAR_TRIANGLE" );
   
      // creating some elements
    IsoparametricLinearQuadrilateral fe_q(2U);
    IsoparametricLinearTriangle      fe_t;
    int32_t                          mtrl_idx{0};
    // this is the constructor used in MeshManager
    csmp::Element<2U>  e0( 0, &fe_q, &quad, evars, ivars, mtrl_idx ),
                       e1( 1, &fe_t, &tria, evars, ivars, mtrl_idx );
   
    csmp::Node<2U> n0, n1, n2, n3, n4;
    // node indexes
    n0.Idx(0); n1.Idx(1); n2.Idx(2); n3.Idx(3); n4.Idx(4);
    // and boundary flags
    n0.AtBoundary(CNR1); n1.AtBoundary(CNR2); n2.AtBoundary(CNR3); n3.AtBoundary(CNR4);
    n4.AtBoundary(TOP);

    // quad 1
    n0.x(0.); n0.y(0.);
    n1.x(2.); n1.y(0.);
    n2.x(2.); n2.y(2.);
    n3.x(0.); n3.y(2.);
    e0.Idx( 0 );
    e0.Assign( 0, &n0 );
    e0.Assign( 1, &n1 );
    e0.Assign( 2, &n2 );
    e0.Assign( 3, &n3 );
    
    // tria above
    n4.x(1.); n4.y(3.);
    e1.Idx( 1 );
    e1.Assign( 0, &n2 );
    e1.Assign( 1, &n4 );
    e1.Assign( 2, &n3 );
    
    // neighbors (only the non-null ones need to be assigned)
    e0.Assign( 3, &e1 );
    e1.Assign( 2, &e0 );
    
    // 2. testing normal copying, assignment and move operation
    //---------------------------------------------------------
    // element 0
    csmp::Element<2U>  e0_copy( e0 );
    _test( e0_copy == e0 );
    csmp::Element<2U>  e0_assigned = e0;
    _test( e0_assigned == e0 );
    csmp::Element<2U>  e0_moved( std::move(e0_copy) );
    _test( e0_moved == e0 );
    // comparison with triangle
    _test( !(e0 == e1) );
    
    // element 1
    csmp::Element<2U>  e1_copy( e1 );
    _test( e1_copy == e1 );
    csmp::Element<2U>  e1_assigned = e0;
    _test( e1_assigned == e0 );

    // 3. testing same operations done via the colony
    //-----------------------------------------------
    // checking that element gets correctly inserted into colony
    /* Note from colony documentation:
       In general T must meet the requirements of Erasable, CopyAssignable and CopyConstructible.
       However, if emplace is utilized to insert elements into the colony, and no functions
       which involve copying or moving are utilized, T is only required to meet the requirements of Erasable.
       If move-insert is utilized instead of emplace, T must also meet the requirements of MoveConstructible.
    */
    
    // copying and assignment
    plf::colony< Element<2> >  elmt_colony;
    auto it_e0 = elmt_colony.insert( e0_assigned );
    auto it_e1 = elmt_colony.insert( e1_copy );
    _test( (*it_e0) == e0 );
    _test( (*it_e1) == e1 );
    auto it_e2 = elmt_colony.insert( e0_moved );
    elmt_colony.erase( it_e0 );
    elmt_colony.erase( it_e2 );
    _test( (*it_e1) == e1 );

    // move semantics
    plf::colony< Element<2> >  elmt_colony2;
    it_e0 = elmt_colony.emplace( e0_assigned );
    it_e1 = elmt_colony.emplace( e1_copy );
    _test( (*it_e0) == e0 );
    _test( (*it_e1) == e1 );
    it_e2 = elmt_colony.emplace( e0_moved );
    elmt_colony.erase( it_e0 );
    elmt_colony.erase( it_e2 );
    _test( (*it_e1) == e1 );

 } // end TestThatElementIsCopyableAssignableMovable




bool Colony_Test::TestColonyWith_Element()
 {
  // making sure that everything is like in a simulation
  const LocalVariables evars( 1, // scalarsVars,
                              2, // vectorVars,
                              1, // tensorVars,
                              0, // array_count,
                              0, // array_length,
                              0, // flag_array_count,
                              0, // flag_array_length,
                              9, // total_data_depth,
                              6 ); // total_flag_depth
                              
  const IntegrationPointVariables ivars;

// Element<dim>( elmt_idx, fem, stencil, evars, cvars, mtrl_idx );
  FiniteVolumeStencil<2>  quad( "ISOPARAMETRIC_LINEAR_QUADRILATERAL" );
  FiniteVolumeStencil<2>  tria( "ISOPARAMETRIC_LINEAR_TRIANGLE" );
 
    // creating some elements
  IsoparametricLinearQuadrilateral fe_q(2U);
  IsoparametricLinearTriangle      fe_t;
  int32_t                          mtrl_idx{0};
  // this is the constructor used in MeshManager
  csmp::Element<2U>  e0( 0, &fe_q, &quad, evars, ivars, mtrl_idx ),
                     e1( 1, &fe_q, &quad, evars, ivars, mtrl_idx ),
                     e2( 2, &fe_t, &tria, evars, ivars, mtrl_idx ),
                     e3( 3, &fe_t, &tria, evars, ivars, mtrl_idx ),
                     e4( 4, &fe_t, &tria, evars, ivars, mtrl_idx );
 
  csmp::Node<2U> n0, n1, n2, n3, n4, n5, n6, n7;
  // quad 1
  n0.x(0.); n0.y(0.);
  n1.x(2.); n1.y(0.);
  n2.x(2.); n2.y(2.);
  n3.x(0.); n3.y(2.);
  e0.Idx( 0 );
  e0.Assign( 0, &n0 );
  e0.Assign( 1, &n1 );
  e0.Assign( 2, &n2 );
  e0.Assign( 3, &n3 );
  // adjacent quad
  n4.x(7.5); n4.y(0.);
  n5.x(7.5); n5.y(2.);
  e1.Idx( 1 );
  e1.Assign( 0, &n1 );
  e1.Assign( 1, &n4 );
  e1.Assign( 2, &n5 );
  e1.Assign( 3, &n2 );
  // tria above
  n6.x(7.5); n6.y(3.);
  n7.x(0.5); n7.y(3.);
  e2.Idx( 2 );
  e2.Assign( 0, &n5 );
  e2.Assign( 1, &n6 );
  e2.Assign( 2, &n7 );
  // tria left of the previous one
  e3.Idx( 3 );
  e3.Assign( 0, &n2 );
  e3.Assign( 1, &n5 );
  e3.Assign( 2, &n7 );
  // tria top left
  e4.Idx( 4 );
  e4.Assign( 0, &n2 );
  e4.Assign( 1, &n7 );
  e4.Assign( 2, &n3 );
  
  // node indexes
  n0.Idx(0); n1.Idx(1); n2.Idx(2); n3.Idx(3); n4.Idx(4); n5.Idx(5); n6.Idx(6); n7.Idx(7);
  // and boundary flags
  n0.AtBoundary(CNR1); n1.AtBoundary(BOTTOM); n2.AtBoundary(NOT); n3.AtBoundary(LEFT);
  n4.AtBoundary(CNR2); n5.AtBoundary(RIGHT); n6.AtBoundary(CNR3); n7.AtBoundary(CNR4);
  
  // neighbors (only the non-null ones need to be assigned)
  e0.Assign( 1, &e1 ); e0.Assign( 2, &e4 );
  e1.Assign( 2, &e3 ); e1.Assign( 3, &e0 );
  e2.Assign( 1, &e3 );
  e3.Assign( 0, &e2 ); e3.Assign( 1, &e4 ); e3.Assign( 2, &e1 );
  e4.Assign( 1, &e0 ); e4.Assign( 2, &e3 );


  // putting copies of these into a colony
  // -------------------------------------
  plf::colony< Element<2> >  elmt_colony;
  // -------------------------------------
  vector <Element<2>*>  eptr( 5 );
  eptr[0] = &( *elmt_colony.insert( e0 ) );
  eptr[1] = &( *elmt_colony.insert( e1 ) );
  eptr[2] = &( *elmt_colony.insert( e2 ) );
  eptr[3] = &( *elmt_colony.insert( e3 ) );
  eptr[4] = &( *elmt_colony.insert( e4 ) );
  if ( verbose_ ) {
       cout <<"\nTestColonyWith_Element2: initial elements.\n";
       for ( const auto& i : elmt_colony ) i.Out();
    }
  
  // copy construction of some more colonies
  plf::colony< Element<2> >  elmt_colony2( elmt_colony ), elmt_colony3( elmt_colony );
  
  // erasure of single element (2) from colony 2
  elmt_colony2.erase( next(elmt_colony2.begin(),2) );
  _test( elmt_colony2.size() == 4 );
    
  // erasing the rest (note return iterator)
  for ( auto it = elmt_colony2.begin(); it != elmt_colony2.end(); /* ++it */ ) {
       it = elmt_colony2.erase(it);
    }
  _test( elmt_colony2.empty() );
  if ( verbose_ ) cout <<"\nelmt_colony2 size: "<< elmt_colony2.size() << endl;
  
  // erasing range of elements from colony 3 using iterators
  // (after this erasure the order of the elements in the colony is not necessarily the original one anymore)
  elmt_colony3.erase( next(elmt_colony3.begin(),2), next(elmt_colony3.begin(),4) );
  _test( elmt_colony3.size() == 3 );
  if ( verbose_ ) cout <<"\nelmt_colony3 size: "<< elmt_colony3.size() << endl;
  
  // using get_iterator with points to delete elements 0,3,5 from colony 1
  elmt_colony.erase( elmt_colony.get_iterator( eptr[0] ) );
  elmt_colony.erase( elmt_colony.get_iterator( eptr[2] ) );
  elmt_colony.erase( elmt_colony.get_iterator( eptr[4] ) );
  _test( elmt_colony.size() == 2 );
  
  // adding those elements back to 'elmt_colony'
  elmt_colony.insert( e0 );
  elmt_colony.insert( e3 );
  elmt_colony.insert( e4 );
  
  // checking the content of the colony with the range operator
  for ( auto i : elmt_colony ) {
       _test( i.FE() != nullptr );
       _test( i.FE_Type() != UNKNOWN );
       _test( i.Idx() < 5 );
    }
  // checking the content of the colony with iterators
  for ( auto it=elmt_colony.begin(); it!=elmt_colony.end(); ++it ) {
       _test( (*it).FE() != nullptr );
       _test( (*it).FE_Type() != UNKNOWN );
       _test( (*it).Idx() < 5 );
    }
  
  if ( verbose_ ) {
      cout <<"\nTestColonyWith_Element2: rebuilt colony with all original elements.\n";
      for ( auto i : elmt_colony ) i.Out();
    }
  
  return true;

} // end TestColonyWith_Element
  
  
  
  
} //end namespace csmp





