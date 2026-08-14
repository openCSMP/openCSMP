#include "Element_Test.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "IsoparametricLinearTriangle.h"
#include "Element.h"
#include "PropertyDatabase.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

Element_Test::Element_Test()
{
}
  
  
/**   
    Testing:
 
      1. LengthInDirection() - quadrilateral 2D
      2. LengthInDirection() - quadrilateral 3D
 
      3. UnitNormalToFace() - actually testing FEPolicy functionality
 
*/
void Element_Test::run()
  {
    MoveSemanticsTest();

    ElementLengthTest2D();
    ElementLengthTest3D();
    
    VariableAccessAndIterators();
  }




void Element_Test::MoveSemanticsTest()
 {
  // 0. creating a mini-element patch for testing
  // --------------------------------------------
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
  // this is the constructor used by the MeshManager
  csmp::Element<2U>  e0( 0, &fe_q, &quad, evars, ivars, mtrl_idx ),
                     e1( 1, &fe_q, &quad, evars, ivars, mtrl_idx ),
                     e2( 2, &fe_t, &tria, evars, ivars, mtrl_idx ),
                     e3( 3, &fe_t, &tria, evars, ivars, mtrl_idx ),
                     e4( 4, &fe_t, &tria, evars, ivars, mtrl_idx );
                     
  // writing some non-NaN data into the variable storage
  Element<2>::Data edata;
  edata.flags = { PLAIN, ANY, INIT_GUESS, INIT_COND, FIELD_DATA, PERIODIC, DIRICH, NEUMANN, ROBIN };
  edata.data  = { 0., 1., 2., 3., 4., 5., 6., 7., 8. };
  e0.LVS( edata );
  e1.LVS( edata );
  e2.LVS( edata );
  e3.LVS( edata );
  e4.LVS( edata );
 
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
  
  
  // 1. testing the copy and move assigments and constructors
  // --------------------------------------------------------

  // copy constructor
  // ----------------
  if ( verbose_ ) cout <<"\n\n"<<"Element_Test::MoveSemanticsTest: "<< endl;
  if ( verbose_ ) cout <<"\n\t"<<"testing copy construction:"<< endl;
  csmp::Element<2U> e1_copy( e1 );
  // inbuilt assignment operator
  _test( e1_copy == e1 );
  CompareElements( e1_copy, e1 );
    
  // assignment operator
  // -------------------
  if ( verbose_ ) cout <<"\n\t"<<"testing assignment:"<< endl;
  csmp::Element<2U> e1_assigned( &fe_t );
  e1_assigned = e1;
  _test( e1_assigned == e1 );
  CompareElements( e1_assigned, e1 );
    
  // move constructor
  // ----------------
  if ( verbose_ ) cout <<"\n\t"<<"testing move constructor:"<< endl;
  csmp::Element<2U> elmt_to_move( e1 );
  // move construction triggered by 'vector::emplace_back'
  vector<csmp::Element<2U>> moved_elmt;
  moved_elmt.emplace_back( std::move(elmt_to_move) ); 
  _test( moved_elmt[0] == e1 );
  CompareElements( moved_elmt[0], e1 );

  // move assignment operator
  // ------------------------
  if ( verbose_ ) cout <<"\n\t"<<"testing move assignment:"<< endl;
  csmp::Element<2U> e1_moved2( &fe_t );
  e1_moved2 = std::move( e1_copy );
  _test( e1_moved2 == e1 );
  CompareElements( e1_moved2, e1 );

  if ( verbose_ ) cout <<"\n\t"<<"testing moved variable storage:"<< endl;
  if ( verbose_ ) moved_elmt[0].OutLVS();
  
  // forced move assignment to get a completely initialised element
  // --------------------------------------------------------------
  if ( verbose_ ) cout <<"\n\t"<<"testing forced move assignment:"<< endl;
  csmp::Element<2U> e0_move_assigned( &fe_q, &quad );
  e0_move_assigned = std::move( e0 );
  _test( e0_move_assigned.Idx() == 0 );
  _test( e0_move_assigned.N(0) == &n0 );
  _test( e0_move_assigned.N(1) == &n1 );
  _test( e0_move_assigned.N(2) == &n2 );
  _test( e0_move_assigned.N(3) == &n3 );
  _test( e0_move_assigned.N(0)->AtBoundary() == CNR1 );
  _test( e0_move_assigned.N(1)->AtBoundary() == BOTTOM );
  _test( e0_move_assigned.N(2)->AtBoundary() == NOT );
  _test( e0_move_assigned.N(3)->AtBoundary() == LEFT );
  _test( e0_move_assigned.Neighbor(0) == nullptr );
  _test( e0_move_assigned.Neighbor(1) == &e1 );
  _test( e0_move_assigned.Neighbor(2) == &e4 );
  _test( e0_move_assigned.Neighbor(3) == nullptr );

 } // end MoveSemanticsTest



template<uint32_t dim>
void Element_Test::CompareElements( const Element<dim>& a, const Element<dim>& b )
 {
   // testing the policies
   _test( a.FE() == b.FE() );
   _test( a.FV() == b.FV() );
   if ( a.FE() && b.FE() )
     _test( a.FE_Type() == b.FE_Type() );
   
   _test( a.Idx()         == b.Idx() );
   _test( a.Material_ID() == b.Material_ID() );
   _test( a.Region_ID()   == b.Region_ID() );

   // local variable storage
//   _test( dynamic_cast<typename Element<dim>::LocalVariableStorage>(a) == b.LVS() ); // dynamic_cast<typename Element<dim>::LocalVariableStorage>(b) );
   auto adata = a.LVS();
   auto bdata = b.LVS();
   _test( adata.flags == bdata.flags );
   _test( adata.data  == bdata.data );

  // connected nodes
  _test( a.Nodes() == b.Nodes() );
  if ( a.Nodes() == b.Nodes() )
    for ( uint32_t i{0u}; i<a.Nodes(); ++i ) {
         _test( a.N(i) == b.N(i) );
         // making sure that there are nodes present
         if ( a.N(i) && b.N(i) )
           _test( a.N(i)->AtBoundary() == b.N(i)->AtBoundary() );
      }
      
  // neighbors
  _test( a.Neighbors() == b.Neighbors() );
  _test( a.ConnectedNeighbors() == b.ConnectedNeighbors() );
  if ( a.Neighbors() == b.Neighbors() )
     for ( uint32_t i{0u}; i<a.Neighbors(); ++i ) {
          _test( a.Neighbor(i) == b.Neighbor(i) );
          if ( (a.Neighbor(i) && b.Neighbor(i)) &&
               (a.Neighbor(i)->FE() && b.Neighbor(i)->FE()) )
            _test( a.Neighbor(i)->FE_Type() == b.Neighbor(i)->FE_Type() );
          // uses element member function for comparison
          _test( a.Neighbor(i) == b.Neighbor(i) );
       }
 }





/**
    Create a fully-constructed element, including variable storage and ascertaining
    that variables cannot be written to when a constant iterator is used.
*/
void Element_Test::VariableAccessAndIterators()
 {
    PropertyDatabase<2> dbase("CSMP-variables.txt");
    
    // making sure that everything is like in a simulation
    const LocalVariables evars = dbase.LocalVariablesAt( NODE );
    /*
                              ( 1, // scalarsVars,
                                2, // vectorVars,
                                1, // tensorVars,
                                0, // array_count,
                                0, // array_length,
                                0, // flag_array_count,
                                0, // flag_array_length,
                                9, // total_data_depth,
                                6 ); // total_flag_depth
    */
    const IntegrationPointVariables ivars = dbase.IntegrationPointVariablesAt( ELEMENT );

  // Element<dim>( elmt_idx, fem, stencil, evars, cvars, mtrl_idx );
    FiniteVolumeStencil<2>  quad_fv( "ISOPARAMETRIC_LINEAR_QUADRILATERAL" );
   
      // creating some elements
    IsoparametricLinearQuadrilateral fe_q(2U);
    int32_t                          mtrl_idx{0};
    // this is the constructor used by the MeshManager
    csmp::Element<2U>  quad( 0, &fe_q, &quad_fv, evars, ivars, mtrl_idx );
   
    // Nodes - constructor: Node( size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT );
    csmp::Node<2U> n0( 0, Point<2>(0.,0.), evars, CNR1, EXTERIOR_POINT ),
                   n1( 1, Point<2>(1.,0.), evars, CNR2, EXTERIOR_POINT ),
                   n2( 2, Point<2>(1.,1.), evars, CNR3, EXTERIOR_POINT ),
                   n3( 3, Point<2>(0.,1.), evars, CNR4, EXTERIOR_POINT );
    
    // assigning the nodes
    quad.Assign( 0, &n0 );
    quad.Assign( 1, &n1 );
    quad.Assign( 2, &n2 );
    quad.Assign( 3, &n3 );
    
    // assigning some property values
    const csmp::Index nvar_key = dbase.StorageKey("nodal variable");
    n0.Store( nvar_key, makeScalar(ANY,1.) );
    n1.Store( nvar_key, makeScalar(ANY,2.) );
    n2.Store( nvar_key, makeScalar(ANY,3.) );
    n3.Store( nvar_key, makeScalar(ANY,4.) );
    
    const csmp::Index evar_key = dbase.StorageKey("element variable");
    quad.Store( evar_key, makeScalar(PLAIN,1.0e-12) );
    
    // testing write access to node variable
    const csmp::Element<2U>&  quad_ref = quad;
    
    for ( auto nit=quad_ref.NodesBegin(); nit!=quad_ref.NodesEnd(); nit++  ) {
          double var = (*nit)->Read( nvar_key );
          _test( var >= 1. && var <= 4. );
          (*nit)->Store( nvar_key, makeScalar(ANY,3.) );
      }
      
    // testing a copy of the element
    csmp::Element<2U>  quad1( quad );
    _test( approximatelyEqual( quad1.Read(evar_key),1.0e-12) );
    _test( approximatelyEqual( quad1.Volume(),1.0) );
    
    // and an assignment
    csmp::Element<2U>  quad2 = quad1;
    _test( approximatelyEqual( quad2.Read(evar_key),1.0e-12) );
    _test( approximatelyEqual( quad2.Volume(),1.0) );

    // and move semantics
    csmp::Element<2U>  quad3( Element<2U>( 0, &fe_q, &quad_fv, evars, ivars, mtrl_idx ) );
    _test( isnan(quad3.Read(evar_key)) );
    _test( quad3.FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL );

    cout <<"\nElement_Test:VariableAccessAndIterators: done"<< endl;
    
} // end VariableAccessAndIterators



  
void Element_Test::ElementLengthTest2D()
{
  IsoparametricLinearQuadrilateral fe(2U);
  csmp::Element<2U>  e( &fe );
 
  csmp::Node<2U> n1, n2, n3, n4;
  n1.x(8.55096); n1.y(8.59504);
  n2.x(12.3416); n2.y(0.396694);
  n3.x(5.37741); n3.y(-4.71625);
  n4.x(5.90634); n4.y(4.27548);
    
  e.Idx( 1 );
  e.Assign( 0, &n1 );
  e.Assign( 1, &n2 );
  e.Assign( 2, &n3 );
  e.Assign( 3, &n4 );

  VectorVariable<2U> direction(PLAIN, PLAIN);
  double fLength(0.);
  
  //test 1
  //direction
  direction(0)=4.31956-12.3416;
  direction(1)=8.8595+4.45179;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 9.763, fTolerance);
  
  //test 2
  //direction
  direction(0)=12.2643-8.36609;
  direction(1)=7.48967+5.60591;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 13.663, fTolerance);
  
  //test 3a
  //direction
  direction(0)=-1.85124-19.0413;
  direction(1)=7.09642-7.09642;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);
  
  //test 3b
  //direction
  direction(0)=4-10;
  direction(1)=-3+3;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);

  //test 4
  //direction
  direction(0)=0-0;
  direction(1)=3-0;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 13.311, fTolerance);
  
} // end ElementLengthTest2D
  
  
  
  
void Element_Test::ElementLengthTest3D()
{
  IsoparametricLinearQuadrilateral fe(3U);
  csmp::Element<3U> e( &fe );

  Node<3U> n1, n2, n3, n4;
  n1.x(5.90634); n1.y(5       ); n1.z(7.);
  n2.x(10     ); n2.y(8.59504 ); n2.z(3.);
  n3.x(12.3416); n3.y(0.396694); n3.z(0.);
  n4.x(5.37741); n4.y(-4.71625); n4.z(0.);
 
  e.Idx( 1 );
  e.Assign( 0, &n1 );
  e.Assign( 1, &n2 );
  e.Assign( 2, &n3 );
  e.Assign( 3, &n4 );
  
  VectorVariable<3U> direction(PLAIN, PLAIN);
  double fLength(0.);
  
  //test 1
  //direction
  direction(0)=12+12;
  direction(1)=-18-18;
  direction(2)=-7-7;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\n3d Length: " << fLength << std::endl;
  _equal(fLength, 9.568, fTolerance);
  
   //test 1
  //direction
  direction(0)=20;
  direction(1)=0;
  direction(2)=0;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\n3d Length: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);
}
  
  
  





   
  
  
} //end namespace csmp





