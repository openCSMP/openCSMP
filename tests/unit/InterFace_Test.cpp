#include "InterFace_Test.h"
#include "Element.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearLineElement.h"
#include "IsoparametricQuadraticTriangle.h"
#include "IsoparametricQuadraticLineElement.h"
#include "IsoparametricLinearTetrahedron.h"
#include "Exception.h"


using namespace std;

namespace csmp {

InterFace_Test::InterFace_Test()
{
}

InterFace_Test::~InterFace_Test(){

}


/**   
    Testing:
 
     InterFace creation and operations
 
*/

void InterFace_Test::run()
{
std::cout << "InterFace Test is Being Run " << std::endl;

    Assign_tests();
    Geometry_tests();

    NodeCoordinateMatrix_linear_test();

    Assign_Geometry_quadratic_test();
    Assign_Geometry_quadratic2_test();

    Assign_tests_linear_3D();


}
	
void InterFace_Test::Assign_tests(){

    std::cout << "Running Test: Assign_Tests" << std::endl;

    IsoparametricLinearTriangle fe(2U);
    csmp::Element<2U>  e1( &fe ), e2(&fe), e_nb_in1(&fe), e_nb_in2(&fe), e_nb_ou1(&fe), e_nb_ou2(&fe);
    csmp::Node<2U> n1, n2, n3, n4, n5, n6;

    // Coordinates needed as some methods rely on co-location! So make sure that nodes on (2,4) and (3,6) are co-located
    n1.x(1.0); n1.y(1.0);
    n2.x(2.0); n2.y(1.0); //topologically collocated with n4
    n3.x(2.0); n3.y(2.0); //topologically collocated with n6
    n4.x(2.0); n4.y(1.0); //topologically collocated with n2
    n5.x(3.0); n5.y(1.0);
    n6.x(2.0); n6.y(2.0); //topologically collocated with n3
    ConstructInterFaceIngredients(e1, e2, n1, n2, n3, n4, n5, n6, e_nb_in1, e_nb_in2, e_nb_ou1, e_nb_ou2);

    //Knowledge of face ID
    uint32_t faceID_of_e1 = 0;
    uint32_t faceID_of_e2 = 1;

    //Interface Construction
    // ---------------------
    // simplemost
    IsoparametricLinearLineElement face_FE;
    FiniteVolumeStencil<2>         line_stencil("ISOPARAMETRIC_LINEAR_BAR");
    const LocalVariables           lvars; // empty
    IntegrationPointVariables      ivars; // empty
    InterFace<2> if_obj0( &face_FE, &line_stencil, lvars, ivars ), if_obj1( &face_FE, &line_stencil, lvars, ivars ),
                 if_obj2( &face_FE, &line_stencil, lvars, ivars ), if_obj3( &face_FE, &line_stencil, lvars, ivars );


    ///Beginning Use and Tests
    /// -----------------------------------------------------------------
    ///Assign 0
    //Only assign elements
    if_obj0.Assign( &e1, faceID_of_e1, &e2, faceID_of_e2 );
    _test(e1 == *(if_obj0.InnerParent()) );
    _test(e2 == *(if_obj0.OuterParent()) );
    //assigns nodes manually
    if_obj0.Assign(0, &n2, INSIDE);
    if_obj0.Assign(1, &n3, INSIDE);
    if_obj0.Assign(0, &n6, OUTSIDE);
    if_obj0.Assign(1, &n4, OUTSIDE);
    //test nodes assignment
    _test( n2 ==  *(if_obj0.InnerParent()->N( if_obj0.ParentNodeNumber(size_t(0), INSIDE) ) ) );
    _test( n3 ==  *(if_obj0.InnerParent()->N( if_obj0.ParentNodeNumber(size_t(1), INSIDE) ) ) );
    _test( n6 ==  *(if_obj0.OuterParent()->N( if_obj0.ParentNodeNumber(size_t(2), OUTSIDE) ) ) );
    _test( n4 ==  *(if_obj0.OuterParent()->N( if_obj0.ParentNodeNumber(size_t(3), OUTSIDE) ) ) );


    ///Assign 1
    if_obj1.Assign(&e1, faceID_of_e1, INSIDE);
    if_obj1.Assign(&e2, faceID_of_e2, OUTSIDE);
    _test(e1 == *(if_obj1.InnerParent()) );
    _test(e2 == *(if_obj1.OuterParent()) );
    //Manually assign must then be used, this has already been tested.

    ///Assign 2
    //Assigns Element and nodes
    if_obj2.AssignElementsAndNodes( &e1, &e2 );
    //Test Element sides
    _test( e1 == *(if_obj2.InnerParent()) );
    _test( e2 == *(if_obj2.OuterParent()) );
    //test node assigned
    _test( n2 ==  *(if_obj2.InnerParent()->N( if_obj2.ParentNodeNumber(size_t(0), INSIDE) ) ) );
    _test( n3 ==  *(if_obj2.InnerParent()->N( if_obj2.ParentNodeNumber(size_t(1), INSIDE) ) ) );
    _test( n6 ==  *(if_obj2.OuterParent()->N( if_obj2.ParentNodeNumber(size_t(2), OUTSIDE) ) ) );
    _test( n4 ==  *(if_obj2.OuterParent()->N( if_obj2.ParentNodeNumber(size_t(3), OUTSIDE) ) ) );


    ///Assign 3
    if_obj3.Assign( &e1, faceID_of_e1, &e2, faceID_of_e2 );
    // connect the nodes
    if_obj3.InitialiseNodeVector();
    //Testing correct element sides
    _test(e1 == *(if_obj3.InnerParent()) );
    _test(e2 == *(if_obj3.OuterParent()) );
    //Testing Node numbering
    _test( n2 ==  *(if_obj3.InnerParent()->N( if_obj3.ParentNodeNumber(size_t(0), INSIDE) ) ) );
    _test( n3 ==  *(if_obj3.InnerParent()->N( if_obj3.ParentNodeNumber(size_t(1), INSIDE) ) ) );
    _test( n6 ==  *(if_obj3.OuterParent()->N( if_obj3.ParentNodeNumber(size_t(2), OUTSIDE) ) ) );
    _test( n4 ==  *(if_obj3.OuterParent()->N( if_obj3.ParentNodeNumber(size_t(3), OUTSIDE) ) ) );

    ///Neighbour Testing
    //Neighbour assign functionality
    if_obj0.Assign(0, &if_obj1);
    if_obj0.Assign(1, &if_obj2);
    _test( if_obj0.Neighbors() == 2 );
    _test( if_obj0.ConnectedNeighbors() == 2);
    _test( if_obj1 == *if_obj0.Neighbor(0) );
    _test( if_obj2 == *if_obj0.Neighbor(1) );
    //testing face before disconnected neighbors
    _test(if_obj0.Faces() == 2);
    //Unassagning neighbours
    if_obj0.Unassign( &if_obj1 );
    _test( if_obj0.Neighbors() == 2);
    _test( if_obj0.ConnectedNeighbors() == 1);
    _test( if_obj0.Neighbor(1) == &if_obj2);
    _test( if_obj0.Neighbor(0) == nullptr );
    //Testing Face
    _test(if_obj0.Faces() == 2);

    //Testing Node call
    _test( &n2 == if_obj0.N(0,INSIDE));
    _test( &n6 == if_obj0.N(0,OUTSIDE));


    //Testing matching nodes
    for (uint32_t n{0U}; n < face_FE.Nodes() ; ++n){
      _test( if_obj0.MatchingN(n,INSIDE)->Coordinate() == if_obj0.MatchingN(n,OUTSIDE)->Coordinate());
      _test( if_obj3.MatchingN(n,INSIDE)->Coordinate() == if_obj3.MatchingN(n,OUTSIDE)->Coordinate());
      _test( if_obj2.MatchingN(n,INSIDE)->Coordinate() == if_obj2.MatchingN(n,OUTSIDE)->Coordinate());
    }


    std::cout << "Finished Running: Assign tests" << std::endl;
}





void InterFace_Test::Geometry_tests(){

    std::cout << "Running Tests: Geometry tests" << std::endl;
    //TODO : Change ordering of local nodes so that node number 3 is actually ordered as 0 in element
    IsoparametricLinearTriangle fe(2U);
    csmp::Element<2U>  e1( &fe ), e2(&fe), e_nb_in1(&fe), e_nb_in2(&fe), e_nb_ou1(&fe), e_nb_ou2(&fe);

    csmp::Node<2U> n1, n2, n3, n4, n5, n6;
    n1.x(1.0); n1.y(1.0);
    n2.x(2.0); n2.y(1.0); //topologically collocated with n4
    n3.x(2.0); n3.y(3.0); //topologically collocated with n6
    n4.x(2.0); n4.y(1.0); //topologically collocated with n2
    n5.x(3.0); n5.y(1.0);
    n6.x(2.0); n6.y(3.0); //topologically collocated with n3

    ConstructInterFaceIngredients(e1, e2, n1, n2, n3, n4, n5, n6, e_nb_in1, e_nb_in2, e_nb_ou1, e_nb_ou2);
    //Knowledge of face ID -- a result of this particular constrution
    uint32_t faceID_of_e1 = 0;
    uint32_t faceID_of_e2 = 1;

    //Construc InterFace objects
    IsoparametricLinearLineElement face_FE;
    FiniteVolumeStencil<2>         line_stencil("ISOPARAMETRIC_LINEAR_BAR");
    const LocalVariables           lvars; // empty
    IntegrationPointVariables      ivars; // empty

    InterFace<2>  if_obj0( &face_FE, &line_stencil, lvars, ivars ),
                  if_obj1( &face_FE, &line_stencil, lvars, ivars ),
                  if_obj2( &face_FE, &line_stencil, lvars, ivars );

    //Assign 0
    if_obj0.AssignElementsAndNodes( &e1, &e2 );

    //Assign 1
    if_obj1.Assign(&e1, faceID_of_e1, INSIDE);
    if_obj1.Assign(&e2, faceID_of_e2, OUTSIDE);
    //assigns nodes manually
    if_obj1.Assign(0, &n2, INSIDE);
    if_obj1.Assign(1, &n3, INSIDE);
    if_obj1.Assign(0, &n6, OUTSIDE);
    if_obj1.Assign(1, &n4, OUTSIDE);
    //Assign 2
    if_obj2.Assign( &e1, faceID_of_e1, &e2, faceID_of_e2 );
    if_obj2.InitialiseNodeVector();



    ///Beginning Geometry Tests
    /// -------------------------------------------------------------------------------------------------

    //Area test - done for all automatic assign methods
    //ONLY 2D Tested!
    _test(if_obj0.Area(INSIDE)  == 2);
    _test(if_obj1.Area(INSIDE)  == 2);
    _test(if_obj2.Area(INSIDE)  == 2);
    _test(if_obj0.Area(OUTSIDE) == 2);
    _test(if_obj1.Area(OUTSIDE) == 2);
    _test(if_obj2.Area(OUTSIDE) == 2);
    
/* SKM: disrupts testing due to the exceptions thrown    
#ifdef NDEBUG
    try {
        if_obj0.Area(MIDDLE);
        _test(false);
    } catch (csmp::Exception throw_area) {
        _test(true);
    }
    try {
        if_obj1.Area(MIDDLE);
        _test(false);
    } catch (csmp::Exception throw_area) {
        _test(true);
    }
    try {
        if_obj2.Area(MIDDLE);
        _test(false);
    } catch (csmp::Exception throw_area) {
        _test(true);
    }
#endif
*/

    //Unit Normal Test
    Point<2U> nrml, n_inside{1.,0.}, n_outside{-1.,0.}, n_middle;
    n_middle = n_inside;

    //starting test
    nrml = if_obj0.UnitNormal(INSIDE);
    _test(nrml == n_inside);
    nrml = if_obj1.UnitNormal(INSIDE);
    _test(nrml == n_inside);
    nrml = if_obj2.UnitNormal(INSIDE);
    _test(nrml == n_inside);
    nrml = if_obj0.UnitNormal(OUTSIDE);
    _test(nrml == n_outside);
    nrml = if_obj1.UnitNormal(OUTSIDE);
    _test(nrml == n_outside);
    nrml = if_obj2.UnitNormal(OUTSIDE);
    _test(nrml == n_outside);

    // Testing that unit normal throws exception when it has no base element
#ifdef NDEBUG //TODO: Remove later so tests in debug still appear
    try {
        if_obj0.UnitNormal(nrml, MIDDLE);
        _test( false);                                        //test failed
    } catch ( csmp::Exception normal_thrown ) { _test(true);} //test passed
    try {
        if_obj1.UnitNormal(nrml, MIDDLE);
        _test( false );
    } catch (csmp::Exception normal_throw) { _test(true); }
    try {
        if_obj2.UnitNormal(nrml, MIDDLE);
        _test(false);
    } catch (csmp::Exception normal_throw) { _test(true); }
#endif
    //Testing normal without middle element
    /* should throw exception
    Point<2> nrml_in_out = if_obj0.UnitNormal();
    _test( n_middle == VectorVariable<2>(nrml_in_out));
    nrml_in_out = if_obj1.UnitNormal();
    _test( n_middle == VectorVariable<2>(nrml_in_out));
    nrml_in_out = if_obj2.UnitNormal();
    _test( n_middle == VectorVariable<2>(nrml_in_out));
   */


    //Test Barycenter of interface
    _test( if_obj0.BaryCenter() == Point<2>(std::vector<double>(2.,2.) ));
    _test( if_obj1.BaryCenter() == Point<2>(std::vector<double>(2.,2.) ));
    _test( if_obj2.BaryCenter() == Point<2>(std::vector<double>(2.,2.) ));

    //Testing length in direction
    VectorVariable<2> vec;
    vec.Component(0, 0.0);
    vec.Component(1, 1.0);
    _test( if_obj0.LengthInDirection(vec) == 2.0);
    _test( if_obj1.LengthInDirection(vec) == 2.0);
    _test( if_obj2.LengthInDirection(vec) == 2.0);

    vec.Component(0, 1.0);
    vec.Component(1, 1.0);
    _equal( sqrt(2.0), if_obj0.LengthInDirection(vec), 1.0E-4);
    _equal( sqrt(2.0), if_obj1.LengthInDirection(vec), 1.0E-4);
    _equal( sqrt(2.0), if_obj2.LengthInDirection(vec), 1.0E-4);

    VectorVariable<2> direction( nrml );
    _test( if_obj0.LengthInDirection(direction) == 0.0);
    _test( if_obj1.LengthInDirection(direction) == 0.0);
    _test( if_obj2.LengthInDirection(direction) == 0.0);


    //Testing Coordinate matrix member Element, Face, InterFace
    DenseMatrix<DM_MIN> XY_inside(2, 2);
    XY_inside.AssignRow(0, n2.Coordinate());
    XY_inside.AssignRow(1, n3.Coordinate());

    DenseMatrix<DM_MIN> XY_outside(2,2);
    XY_outside.AssignRow(0, n6.Coordinate());
    XY_outside.AssignRow(1, n4.Coordinate());

    DenseMatrix<DM_MIN> XY(2,2);
    if_obj0.CoordinateMatrix();
    _test( XY_inside == if_obj0.FE()->XY );
    if ( verbose_ ) {
        XY_inside.Out();
        if_obj0.FE()->XY.Out();
      }
    if_obj0.CurrentSide( OUTSIDE );
    if_obj0.Idx( 734 ); // to prompt update
    if_obj0.CoordinateMatrix();
    _test( XY_outside == if_obj0.FE()->XY );
    if ( verbose_ ) {
        XY_outside.Out();
        if_obj0.FE()->XY.Out();
      }
}



void InterFace_Test::NodeCoordinateMatrix_linear_test( ) {


  std::cout << "Running Test: NodeCoordinateMatrix_test" << std::endl;

  IsoparametricLinearTriangle fe(2U);
  csmp::Element<2U>  e1( &fe ), e2(&fe), e_nb_in1(&fe), e_nb_in2(&fe), e_nb_ou1(&fe), e_nb_ou2(&fe);
  csmp::Node<2U> n1, n2, n3, n4, n5, n6;

  // Coordinates needed as some methods rely on co-location! So make sure that nodes on (2,4) and (3,6) are co-located
  n1.x(1.0); n1.y(1.0);
  n2.x(2.0); n2.y(1.0); //topologically collocated with n4
  n3.x(3.0); n3.y(3.0); //topologically collocated with n6
  n4.x(2.0); n4.y(1.0); //topologically collocated with n2
  n5.x(3.0); n5.y(1.0);
  n6.x(3.0); n6.y(3.0); //topologically collocated with n3
  ConstructInterFaceIngredients(e1, e2, n1, n2, n3, n4, n5, n6, e_nb_in1, e_nb_in2, e_nb_ou1, e_nb_ou2);

  //Knowledge of face ID
  uint32_t faceID_of_e1 = 0;
  uint32_t faceID_of_e2 = 1;

  //Interface Construction
  IsoparametricLinearLineElement face_FE;
  FiniteVolumeStencil<2>         line_stencil("ISOPARAMETRIC_LINEAR_BAR");
  const LocalVariables           lvars; // empty
  IntegrationPointVariables      ivars; // empty

  InterFace<2>  if_obj0( &face_FE, &line_stencil, lvars, ivars ),
                if_obj1( &face_FE, &line_stencil, lvars, ivars ),
                if_obj2( &face_FE, &line_stencil, lvars, ivars );

  //Assign 0
  if_obj0.AssignElementsAndNodes( &e1, &e2 );

  //Assign 1
  if_obj1.Assign(&e1, faceID_of_e1, INSIDE);
  if_obj1.Assign(&e2, faceID_of_e2, OUTSIDE);
  //assigns nodes manually
  if_obj1.Assign(0, &n2, INSIDE);
  if_obj1.Assign(1, &n3, INSIDE);
  if_obj1.Assign(0, &n6, OUTSIDE);
  if_obj1.Assign(1, &n4, OUTSIDE);
  //Assign 2
  if_obj2.Assign( &e1, faceID_of_e1, &e2, faceID_of_e2 );
  if_obj2.InitialiseNodeVector();



  ///Beginning Geometry Tests
  /// -------------------------------------------------------------------------------------------------

  //Area test - done for all automatic assign methods
  //ONLY 2D Tested!
  _test(if_obj0.Area(INSIDE)  == std::sqrt(5));
  _test(if_obj1.Area(INSIDE)  == std::sqrt(5));
  _test(if_obj2.Area(INSIDE)  == std::sqrt(5));
  _test(if_obj0.Area(OUTSIDE) == std::sqrt(5));
  _test(if_obj1.Area(OUTSIDE) == std::sqrt(5));
  _test(if_obj2.Area(OUTSIDE) == std::sqrt(5));


  //Unit Normal Test
  Point<2U> nrml, n_inside{2.,-1.}, n_outside{-2.,1.}, n_middle;
  n_inside.NormalizeLengthTo(1.); n_outside.NormalizeLengthTo();
  n_middle = n_inside;

  //starting test
  nrml = if_obj0.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj1.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj2.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj0.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);
  nrml = if_obj1.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);
  nrml = if_obj2.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);


  ///Beginning NodeCoordinateMatrix_Tests and Finite ELement shape functions
  /// -------------------------------------------------------------------------------------------------


  double l0 = std::sqrt(0.75*0.75 + 1.5*1.5);
  double shape_function_N0 =  l0 / std::sqrt(5);
  double shape_function_N1 =  1.0-shape_function_N0;
  vector<double> N{0,0}, xyz{2.25,1.5};


  ///Starting test to check that node coordinate matrix is consistent with node-numbering of the face, and therefore uses
  /// the correct finite element interpolation functions N = [N1,N2] with each node

  //Inside
  if_obj0.CurrentSide(INSIDE);
  if_obj0.N_AtGlobalPoint(N, xyz );
  _test( shape_function_N0 == N[0] );
  _test( shape_function_N1 == N[1] );


  N[0]=0; N[1]=0;
  face_FE.CurrentID( face_FE.CurrentID() + 1 );
  if_obj1.CurrentSide(INSIDE);
  if_obj1.N_AtGlobalPoint(N, xyz );
  _test( shape_function_N0 == N[0] );
  _test( shape_function_N1 == N[1] );

  N[0]=0; N[1]=0;
  face_FE.CurrentID( face_FE.CurrentID() + 1 );
  if_obj2.CurrentSide(INSIDE);
  if_obj2.N_AtGlobalPoint(N, xyz );
  _test( shape_function_N0 == N[0] );
  _test( shape_function_N1 == N[1] );

  std::cout << "Shape function N0 " << N[0] << " = " << shape_function_N0 << std::endl;
  std::cout << "Shape function N1 " << N[1] << " = " << shape_function_N1 << std::endl;

  //Outside - NOW THE SHAPE FUNCTIONS HAVE SWITCHED SINCE NODE 0 is now at the top
  N[0]=0; N[1]=0;
  face_FE.CurrentID( face_FE.CurrentID() + 1 );
  if_obj0.CurrentSide(OUTSIDE);
  if_obj0.N_AtGlobalPoint(N, xyz );
  _test( shape_function_N1 == N[0] );
  _test( shape_function_N0 == N[1] );

  N[0]=0; N[1]=0;
  face_FE.CurrentID( face_FE.CurrentID() + 1 );
  if_obj1.CurrentSide(OUTSIDE);
  if_obj1.N_AtGlobalPoint(N, xyz );
  _test( shape_function_N1 == N[0] );
  _test( shape_function_N0 == N[1] );

  N[0]=0; N[1]=0;
  face_FE.CurrentID( face_FE.CurrentID() + 1 );
  if_obj2.CurrentSide(OUTSIDE);
  if_obj2.N_AtGlobalPoint(N, xyz );
  _test( shape_function_N1 == N[0] );
  _test( shape_function_N0 == N[1] );

  std::cout << "Shape function N1 " << N[0] << " = " << shape_function_N1 << std::endl;
  std::cout << "Shape function N0 " << N[1] << " = " << shape_function_N0 << std::endl;


}




void InterFace_Test::Assign_Geometry_quadratic_test(){


  std::cout << "Running Test: NodeCoordinateMatrix_quadratic_test" << std::endl;

  IsoparametricQuadraticTriangle fe(2U);
  csmp::Element<2U>  e1( &fe ), e2(&fe), e_nb_in1(&fe), e_nb_in2(&fe), e_nb_ou1(&fe), e_nb_ou2(&fe);
  csmp::Node<2U> n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12;

  // Coordinates needed as some methods rely on co-location! So make sure that nodes on (2,4) and (3,6) are co-located
  n1.x(1.0); n1.y(1.0);
  n2.x(2.0); n2.y(1.0); //topologically collocated with n4
  n3.x(3.0); n3.y(3.0); //topologically collocated with n6
  n4.x(2.0); n4.y(1.0); //topologically collocated with n2
  n5.x(3.0); n5.y(1.0);
  n6.x(3.0); n6.y(3.0); //topologically collocated with n3

  n7.x(1.5); n7.y(1.0);
  n8.x(2.5); n8.y(2.0); //collocated with n12
  n9.x(2.0); n9.y(2.0);

  n10.x(2.5); n10.y(1.0);
  n11.x(3.0); n11.y(2.0);
  n12.x(2.5); n12.y(2.0); //collocated with n8


  ConstructInterFaceIngredients_quadratic(e1, e2, n1, n2, n3, n4, n5, n6, e_nb_in1, e_nb_in2, e_nb_ou1, e_nb_ou2, n7, n8, n9, n10, n11, n12);

  //Knowledge of face ID
  uint32_t faceID_of_e1 = 0;
  uint32_t faceID_of_e2 = 1;

  //Interface Construction
  IsoparametricQuadraticLineElement face_FE;
  FiniteVolumeStencil<2>         line_stencil("ISOPARAMETRIC_QUADRATIC_BAR");
  const LocalVariables           lvars; // empty
  IntegrationPointVariables      ivars; // empty

  InterFace<2>  if_obj0( &face_FE, nullptr, lvars, ivars ),
                if_obj1( &face_FE, nullptr, lvars, ivars ),
                if_obj2( &face_FE, nullptr, lvars, ivars );

  //Assign 0
  if_obj0.AssignElementsAndNodes( &e1, &e2 );

  //Assign 1
  if_obj1.Assign(&e1, faceID_of_e1, INSIDE);
  if_obj1.Assign(&e2, faceID_of_e2, OUTSIDE);
  //assigns nodes manually
  if_obj1.Assign(0, &n2, INSIDE);
  if_obj1.Assign(1, &n3, INSIDE);
  if_obj1.Assign(0, &n6, OUTSIDE);
  if_obj1.Assign(1, &n4, OUTSIDE);
  //quadratic
  if_obj1.Assign(2, &n8, INSIDE);
  if_obj1.Assign(2, &n12,OUTSIDE);

  //Assign 2
  if_obj2.Assign( &e1, faceID_of_e1, &e2, faceID_of_e2 );
  if_obj2.InitialiseNodeVector();


  //Testing Assignment
  _test( &n2 == if_obj0.N(0,INSIDE));
  _test( &n3 == if_obj0.N(1,INSIDE));
  _test( &n8 == if_obj0.N(2,INSIDE));

  _test( &n6  == if_obj0.N(0,OUTSIDE));
  _test( &n4  == if_obj0.N(1,OUTSIDE));
  _test( &n12 == if_obj0.N(2,OUTSIDE));


  //Testing Assignment
  _test( &n2 == if_obj1.N(0,INSIDE));
  _test( &n3 == if_obj1.N(1,INSIDE));
  _test( &n8 == if_obj1.N(2,INSIDE));

  _test( &n6  == if_obj1.N(0,OUTSIDE));
  _test( &n4  == if_obj1.N(1,OUTSIDE));
  _test( &n12 == if_obj1.N(2,OUTSIDE));


  //Testing Assignment
  _test( &n2 == if_obj2.N(0,INSIDE));
  _test( &n3 == if_obj2.N(1,INSIDE));
  _test( &n8 == if_obj2.N(2,INSIDE));

  _test( &n6  == if_obj2.N(0,OUTSIDE));
  _test( &n4  == if_obj2.N(1,OUTSIDE));
  _test( &n12 == if_obj2.N(2,OUTSIDE));


  ///Beginning Geometry Tests
  /// -------------------------------------------------------------------------------------------------

  //Area test - done for all automatic assign methods
  //ONLY 2D Tested!
  _equal(if_obj0.Area(INSIDE)  , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj1.Area(INSIDE)  , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj2.Area(INSIDE)  , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj0.Area(OUTSIDE) , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj1.Area(OUTSIDE) , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj2.Area(OUTSIDE) , std::sqrt(5), numeric_limits<double>::epsilon() );


  //Unit Normal Test
  Point<2U> nrml, n_inside{2.,-1.}, n_outside{-2.,1.}, n_middle;
  n_inside.NormalizeLengthTo(1.); n_outside.NormalizeLengthTo();
  n_middle = n_inside;

  //starting test
  nrml = if_obj0.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj1.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj2.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj0.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);
  nrml = if_obj1.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);
  nrml = if_obj2.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);


  //Testing matching nodes
  for (uint32_t n{0U}; n < face_FE.Nodes() ; ++n){
    _test( if_obj0.MatchingN(n,INSIDE)->Coordinate() == if_obj0.MatchingN(n,OUTSIDE)->Coordinate());
    _test( if_obj1.MatchingN(n,INSIDE)->Coordinate() == if_obj1.MatchingN(n,OUTSIDE)->Coordinate());
    _test( if_obj2.MatchingN(n,INSIDE)->Coordinate() == if_obj2.MatchingN(n,OUTSIDE)->Coordinate());
  }


}








//Same as before but with different element node numbering
void InterFace_Test::Assign_Geometry_quadratic2_test(){



  std::cout << "Running Test: Assign_Geometry_quadratic2_test" << std::endl;

  IsoparametricQuadraticTriangle fe(2U);
  csmp::Element<2U>  e1( &fe ), e2(&fe), e_nb_in1(&fe), e_nb_in2(&fe), e_nb_ou1(&fe), e_nb_ou2(&fe);
  csmp::Node<2U> n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12;

  // Coordinates needed as some methods rely on co-location! So make sure that nodes on (2,4) and (3,6) are co-located
  n1.x(1.0); n1.y(1.0);
  n2.x(2.0); n2.y(1.0); //topologically collocated with n4
  n3.x(3.0); n3.y(3.0); //topologically collocated with n6
  n4.x(2.0); n4.y(1.0); //topologically collocated with n2
  n5.x(3.0); n5.y(1.0);
  n6.x(3.0); n6.y(3.0); //topologically collocated with n3

  n7.x(1.5); n7.y(1.0);
  n8.x(2.5); n8.y(2.0); //collocated with n12
  n9.x(2.0); n9.y(2.0);

  n10.x(2.5); n10.y(1.0);
  n11.x(3.0); n11.y(2.0);
  n12.x(2.5); n12.y(2.0); //collocated with n8


  ConstructInterFaceIngredients_quadratic(e1, e2, n1, n2, n3, n5, n6, n4, e_nb_in1, e_nb_in2, e_nb_ou1, e_nb_ou2, n7, n8, n9, n11, n12, n10);

  //Knowledge of face ID
  uint32_t faceID_of_e1 = 0;
  uint32_t faceID_of_e2 = 0;

  //Interface Construction
  IsoparametricQuadraticLineElement face_FE;
  FiniteVolumeStencil<2>         line_stencil("ISOPARAMETRIC_QUADRATIC_BAR");
  const LocalVariables           lvars; // empty
  IntegrationPointVariables      ivars; // empty

  InterFace<2>  if_obj0( &face_FE, nullptr, lvars, ivars ),
                if_obj1( &face_FE, nullptr, lvars, ivars ),
                if_obj2( &face_FE, nullptr, lvars, ivars );

  //Assign 0
  if_obj0.AssignElementsAndNodes( &e1, &e2 );

  //Assign 1
  if_obj1.Assign(&e1, faceID_of_e1, INSIDE);
  if_obj1.Assign(&e2, faceID_of_e2, OUTSIDE);
  //assigns nodes manually
  if_obj1.Assign(0, &n2, INSIDE);
  if_obj1.Assign(1, &n3, INSIDE);
  if_obj1.Assign(0, &n6, OUTSIDE);
  if_obj1.Assign(1, &n4, OUTSIDE);
  //quadratic
  if_obj1.Assign(2, &n8, INSIDE);
  if_obj1.Assign(2, &n12,OUTSIDE);

  //Assign 2
  if_obj2.Assign( &e1, faceID_of_e1, &e2, faceID_of_e2 );
  if_obj2.InitialiseNodeVector();


  //Testing Assignment
  _test( &n2 == if_obj0.N(0,INSIDE));
  _test( &n3 == if_obj0.N(1,INSIDE));
  _test( &n8 == if_obj0.N(2,INSIDE));

  _test( &n6  == if_obj0.N(0,OUTSIDE));
  _test( &n4  == if_obj0.N(1,OUTSIDE));
  _test( &n12 == if_obj0.N(2,OUTSIDE));


  //Testing Assignment
  _test( &n2 == if_obj1.N(0,INSIDE));
  _test( &n3 == if_obj1.N(1,INSIDE));
  _test( &n8 == if_obj1.N(2,INSIDE));

  _test( &n6  == if_obj1.N(0,OUTSIDE));
  _test( &n4  == if_obj1.N(1,OUTSIDE));
  _test( &n12 == if_obj1.N(2,OUTSIDE));


  //Testing Assignment
  _test( &n2 == if_obj2.N(0,INSIDE));
  _test( &n3 == if_obj2.N(1,INSIDE));
  _test( &n8 == if_obj2.N(2,INSIDE));

  _test( &n6  == if_obj2.N(0,OUTSIDE));
  _test( &n4  == if_obj2.N(1,OUTSIDE));
  _test( &n12 == if_obj2.N(2,OUTSIDE));


  ///Beginning Geometry Tests
  /// -------------------------------------------------------------------------------------------------

  //Area test - done for all automatic assign methods
  //ONLY 2D Tested!
  _equal(if_obj0.Area(INSIDE)  , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj1.Area(INSIDE)  , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj2.Area(INSIDE)  , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj0.Area(OUTSIDE) , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj1.Area(OUTSIDE) , std::sqrt(5), numeric_limits<double>::epsilon() );
  _equal(if_obj2.Area(OUTSIDE) , std::sqrt(5), numeric_limits<double>::epsilon() );


  //Unit Normal Test
  Point<2U> nrml, n_inside{2.,-1.}, n_outside{-2.,1.}, n_middle;
  n_inside.NormalizeLengthTo(1.); n_outside.NormalizeLengthTo();
  n_middle = n_inside;

  //starting test
  nrml = if_obj0.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj1.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj2.UnitNormal(INSIDE);
  _test(nrml == n_inside);
  nrml = if_obj0.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);
  nrml = if_obj1.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);
  nrml = if_obj2.UnitNormal(OUTSIDE);
  _test(nrml == n_outside);



  //Testing matching nodes
  for (uint32_t n{0U}; n < face_FE.Nodes() ; ++n){
    _test( if_obj0.MatchingN(n,INSIDE)->Coordinate() == if_obj0.MatchingN(n,OUTSIDE)->Coordinate());
    _test( if_obj1.MatchingN(n,INSIDE)->Coordinate() == if_obj1.MatchingN(n,OUTSIDE)->Coordinate());
    _test( if_obj2.MatchingN(n,INSIDE)->Coordinate() == if_obj2.MatchingN(n,OUTSIDE)->Coordinate());
  }

}








//For 3d Triangular elements
void InterFace_Test::Assign_tests_linear_3D(){

  std::cout << "Running test: Linear 3d tetrahedra Assign" << std::endl;

  //Note: Ansys conventions on tetrahedral elements have the ordering of face 3 to have a clockwise ordering of nodes.
  // Afterwhich if the nodes of the face are called, they are given in a counter clockwise ordering , assumign that the nodes are inserted in this clockwise way.

  IsoparametricLinearTetrahedron  fe;
  csmp::Element<3U>  e1( &fe ), e2(&fe), e_nb_in1(&fe), e_nb_in2(&fe), e_nb_ou1(&fe), e_nb_ou2(&fe);
  csmp::Node<3U> n0, n1, n2, n3, n4, n5, n6, n7;

  n0.x(1.0); n0.y(1.0); n0.z( 0.0);       //topologically collocated with n6
  n1.x(2.0); n1.y(2.0); n1.z( 0.0);       //topologically collocated with n4
  n2.x(3.0); n2.y(1.0); n2.z( 0.0);       //topologically collocated with n5
  n3.x(2.0); n3.y(1.0); n3.z(-1.0);
  n4.x(2.0); n4.y(2.0); n4.z( 0.0);
  n5.x(1.0); n5.y(1.0); n5.z( 0.0);
  n6.x(3.0); n6.y(1.0); n6.z( 0.0);
  n7.x(2.0); n7.y(1.0); n7.z( 1.0);

  ConstructInterFaceIngredients_3D(e1, e2, n0, n1, n2, n3, n4, n5, n6, n7, e_nb_in1, e_nb_in2, e_nb_ou1, e_nb_ou2);

  //Knowledge of face ID
  size_t faceID_of_e1 = 3;
  size_t faceID_of_e2 = 3;


  //Testing unit normals of linear tetrahedron face
  VectorVariable<3> nrml;
  e1.UnitNormalToFace(3, nrml);
  _equal(  0.0, nrml[0], 0.0001);
  _equal(  0.0, nrml[1], 0.0001);
  _equal(  1.0, nrml[2], 0.0001);
  e2.UnitNormalToFace(3, nrml);
  _equal(  0.0, nrml[0], 0.0001);
  _equal(  0.0, nrml[1], 0.0001);
  _equal( -1.0, nrml[2], 0.0001);


  //Interface Construction
  IsoparametricLinearTriangle face_FE(3U);
  const LocalVariables           lvars; // empty
  IntegrationPointVariables      ivars; // empty


  InterFace<3> if_obj0( &face_FE, nullptr, lvars, ivars),
               if_obj1( &face_FE, nullptr, lvars, ivars),
               if_obj2( &face_FE, nullptr, lvars, ivars);



  //Assign 0
  if_obj0.AssignElementsAndNodes(&e1, &e2);

  //Assign 1
  if_obj1.Assign(&e1, faceID_of_e1, INSIDE);
  if_obj1.Assign(&e2, faceID_of_e2, OUTSIDE);
  //nodes manually done
  if_obj1.Assign(0, &n0, INSIDE);
  if_obj1.Assign(1, &n2, INSIDE);
  if_obj1.Assign(2, &n1, INSIDE);
  if_obj1.Assign(0, &n4, OUTSIDE);
  if_obj1.Assign(1, &n6, OUTSIDE);
  if_obj1.Assign(2, &n5, OUTSIDE);

  //assign 2
  if_obj2.Assign( &e1, faceID_of_e1, &e2, faceID_of_e2 );
  if_obj2.InitialiseNodeVector();


  ///Beginning Use and Tests for 3D
  /// -----------------------------------------------------------------

  ///Assign 0
  //test nodes manual assignment  and parent node number calling functionality
  _test(e1 == *(if_obj0.InnerParent()) );
  _test(e2 == *(if_obj0.OuterParent()) );

  _test( n0 == *(if_obj0.N(0,INSIDE)) );
  _test( n2 == *(if_obj0.N(1,INSIDE)) );
  _test( n1 == *(if_obj0.N(2,INSIDE)) );
  _test( n4 == *(if_obj0.N(0,OUTSIDE)) ); //NOTE: CURRENT CONSTRUCTION DOES NOT GUARANTEE THAT N(i,OUTSIDE) = N(I,INSIDE)
  _test( n6 == *(if_obj0.N(1,OUTSIDE)) );
  _test( n5 == *(if_obj0.N(2,OUTSIDE)) );

  //-- THIS DOES NOT GUARANTEE MATCHING NODE NUMBER -
  _test( n0 ==  *(if_obj0.InnerParent()->N( if_obj0.ParentNodeNumber(size_t(0), INSIDE) ) ) );
  _test( n2 ==  *(if_obj0.InnerParent()->N( if_obj0.ParentNodeNumber(size_t(1), INSIDE) ) ) );
  _test( n1 ==  *(if_obj0.InnerParent()->N( if_obj0.ParentNodeNumber(size_t(2), INSIDE) ) ) );
  _test( n4 ==  *(if_obj0.OuterParent()->N( if_obj0.ParentNodeNumber(size_t(3), OUTSIDE) ) ) );
  _test( n6 ==  *(if_obj0.OuterParent()->N( if_obj0.ParentNodeNumber(size_t(4), OUTSIDE) ) ) );
  _test( n5 ==  *(if_obj0.OuterParent()->N( if_obj0.ParentNodeNumber(size_t(5), OUTSIDE) ) ) );


  ///ASSIGN 1
  _test(e1 == *(if_obj1.InnerParent()) );
  _test(e2 == *(if_obj1.OuterParent()) );

  _test( n0 == *(if_obj1.N(0,INSIDE)) );
  _test( n2 == *(if_obj1.N(1,INSIDE)) );
  _test( n1 == *(if_obj1.N(2,INSIDE)) );
  _test( n4 == *(if_obj1.N(0,OUTSIDE)) ); //NOTE: CURRENT CONSTRUCTION DOES NOT GUARANTEE THAT N(i,OUTSIDE) = N(I,INSIDE)
  _test( n6 == *(if_obj1.N(1,OUTSIDE)) );
  _test( n5 == *(if_obj1.N(2,OUTSIDE)) );

  //test node assigned -- nodes should match the other side.
  _test( n0 ==  *(if_obj1.InnerParent()->N( if_obj1.ParentNodeNumber(size_t(0), INSIDE) ) ) );
  _test( n2 ==  *(if_obj1.InnerParent()->N( if_obj1.ParentNodeNumber(size_t(1), INSIDE) ) ) );
  _test( n1 ==  *(if_obj1.InnerParent()->N( if_obj1.ParentNodeNumber(size_t(2), INSIDE) ) ) );
  _test( n4 ==  *(if_obj1.OuterParent()->N( if_obj1.ParentNodeNumber(size_t(3), OUTSIDE) ) ) );
  _test( n6 ==  *(if_obj1.OuterParent()->N( if_obj1.ParentNodeNumber(size_t(4), OUTSIDE) ) ) );
  _test( n5 ==  *(if_obj1.OuterParent()->N( if_obj1.ParentNodeNumber(size_t(5), OUTSIDE) ) ) );


  ///ASSIGN 2
  _test(e1 == *(if_obj2.InnerParent()) );
  _test(e2 == *(if_obj2.OuterParent()) );

  _test( n0 == *(if_obj2.N(0,INSIDE)) );
  _test( n2 == *(if_obj2.N(1,INSIDE)) );
  _test( n1 == *(if_obj2.N(2,INSIDE)) );
  _test( n4 == *(if_obj2.N(0,OUTSIDE)) ); //NOTE: CURRENT CONSTRUCTION DOES NOT GUARANTEE THAT N(i,OUTSIDE) = N(I,INSIDE)
  _test( n6 == *(if_obj2.N(1,OUTSIDE)) );
  _test( n5 == *(if_obj2.N(2,OUTSIDE)) );

  //test node assigned -- nodes should match the other side.
  _test( n0 ==  *(if_obj2.InnerParent()->N( if_obj2.ParentNodeNumber(size_t(0), INSIDE) ) ) );
  _test( n2 ==  *(if_obj2.InnerParent()->N( if_obj2.ParentNodeNumber(size_t(1), INSIDE) ) ) );
  _test( n1 ==  *(if_obj2.InnerParent()->N( if_obj2.ParentNodeNumber(size_t(2), INSIDE) ) ) );
  _test( n4 ==  *(if_obj2.OuterParent()->N( if_obj2.ParentNodeNumber(size_t(3), OUTSIDE) ) ) );
  _test( n6 ==  *(if_obj2.OuterParent()->N( if_obj2.ParentNodeNumber(size_t(4), OUTSIDE) ) ) );
  _test( n5 ==  *(if_obj2.OuterParent()->N( if_obj2.ParentNodeNumber(size_t(5), OUTSIDE) ) ) );


  Point<3U> n_inside(  std::vector<double>{0., 0., 1.0}),
            n_outside( std::vector<double>{0., 0., -1.0});

  //starting test
  Point<3U> normal;
  normal = if_obj0.UnitNormal(INSIDE);
  _test(normal == n_inside);
  normal = if_obj1.UnitNormal(INSIDE);
  _test(normal == n_inside);
  normal = if_obj2.UnitNormal(INSIDE);
  _test(normal == n_inside);
  normal = if_obj0.UnitNormal(OUTSIDE);
  _test(normal == n_outside);
  normal = if_obj1.UnitNormal(OUTSIDE);
  _test(normal == n_outside);
  normal = if_obj2.UnitNormal(OUTSIDE);
  _test(normal == n_outside);


  //Testing matching nodes
  for (uint32_t n{0U}; n < face_FE.Nodes() ; ++n){
    _test( if_obj0.MatchingN(n,INSIDE)->Coordinate() == if_obj0.MatchingN(n,OUTSIDE)->Coordinate());
    _test( if_obj1.MatchingN(n,INSIDE)->Coordinate() == if_obj1.MatchingN(n,OUTSIDE)->Coordinate());
    _test( if_obj2.MatchingN(n,INSIDE)->Coordinate() == if_obj2.MatchingN(n,OUTSIDE)->Coordinate());
  }


  /*
  //Now testing that nodes are indeed matched
  //testing inside nodes
  _test( n0 == *(if_obj2.N(0,INSIDE)) );
  _test( n2 == *(if_obj2.N(1,INSIDE)) );
  _test( n1 == *(if_obj2.N(2,INSIDE)) );
  //testing outside nodes match with inside nodes (only testing linear functionality
  _test( if_obj2.N(0,INSIDE)->Coordinate() == if_obj2.N(0,OUTSIDE)->Coordinate() );
  _test( if_obj2.N(1,INSIDE)->Coordinate() == if_obj2.N(1,OUTSIDE)->Coordinate() );
  _test( if_obj2.N(2,INSIDE)->Coordinate() == if_obj2.N(2,OUTSIDE)->Coordinate() );
*/

}

















//Nodes in1, in2, out0, out1 are on the interface object.
void InterFace_Test::ConstructInterFaceIngredients(Element<2> &e_in, Element<2> &e_out,
                                                   Node<2> &n_in0,  Node<2> &n_in1,  Node<2> &n_in2,
                                                   Node<2> &n_out0, Node<2> &n_out1, Node<2> &n_out2,
                                                   Element<2> &e_nb_in1, Element<2> &e_nb_in2,
                                                   Element<2> &e_nb_ou1, Element<2> &e_nb_ou2)
{
    n_in0.ResizeParentStorage(3); //outer node
    n_in1.ResizeParentStorage(2); //inner nodes...
    n_in2.ResizeParentStorage(2); //...
    n_out0.ResizeParentStorage(2);
    n_out1.ResizeParentStorage(3);
    n_out2.ResizeParentStorage(2); //outer node

    //Assigning Actual Parents To nodes - parent node index to pnode with element e1 in parent elem pointers
    n_in0.Assign(0,&e_in);
    n_in1.Assign(1,&e_in);
    n_in2.Assign(2,&e_in);
    n_out0.Assign(0,&e_out);
    n_out1.Assign(1,&e_out);
    n_out2.Assign(2,&e_out);

    ///Calibrating Neighbours
    //Assigning Neighbour elements to Nodes.
    //Inside
    n_in0.Assign(1, &e_nb_in1);
    n_in1.Assign(0, &e_nb_in1);
    n_in2.Assign(0, &e_nb_in2);
    //Outside
    n_out0.Assign(0, &e_nb_ou1);
    n_out1.Assign(0, &e_nb_ou2);
    n_out2.Assign(1, &e_nb_ou2);

    //Assigning Nodes to Element in InterFace
    //Inside
    e_in.Idx( 1 );
    e_in.Assign( 0, &n_in0 );
    e_in.Assign( 1, &n_in1 );
    e_in.Assign( 2, &n_in2 );
    //Outside
    e_out.Idx(2);
    e_out.Assign( 0, &n_out0 );
    e_out.Assign( 1, &n_out1 );
    e_out.Assign( 2, &n_out2 );

    //Assigning Neighbour element nodes
    //inside
    e_nb_in1.Idx(3);
    e_nb_in1.Assign(0, &n_in1);
    e_nb_in1.Assign(1, &n_in0);
    e_nb_in2.Idx(4);
    e_nb_in2.Assign(0, &n_in2);
    //outside
    e_nb_ou1.Idx(5);
    e_nb_ou1.Assign(0, &n_out0);
    e_nb_ou2.Idx(6);
    e_nb_ou2.Assign(0, &n_out1);
    e_nb_ou2.Assign(1, &n_out2);
    //Note, not all nodes of IsoparametricLinearTriangle are defined for neighbour elements
}

//Nodes in1, in2, out0, out1 are on the interface object.
void InterFace_Test::ConstructInterFaceIngredients_quadratic(Element<2> &e_in, Element<2> &e_out,
                                                   Node<2> &n_in0,  Node<2> &n_in1,  Node<2> &n_in2,
                                                   Node<2> &n_out0, Node<2> &n_out1, Node<2> &n_out2,
                                                   Element<2> &e_nb_in1, Element<2> &e_nb_in2,
                                                   Element<2> &e_nb_ou1, Element<2> &e_nb_ou2,
                                                   Node<2> &n_in3,  Node<2> &n_in4,  Node<2> &n_in5,
                                                   Node<2> &n_out3, Node<2> &n_out4, Node<2> &n_out5)
{
    n_in0.ResizeParentStorage(3); //outer node
    n_in1.ResizeParentStorage(2); //inner nodes...
    n_in2.ResizeParentStorage(2); //...
    n_out0.ResizeParentStorage(2);
    n_out1.ResizeParentStorage(3);
    n_out2.ResizeParentStorage(2); //outer node

    n_in3.ResizeParentStorage(2); //outer node
    n_in4.ResizeParentStorage(1); //inner nodes...
    n_in5.ResizeParentStorage(2); //...
    n_out3.ResizeParentStorage(2);
    n_out4.ResizeParentStorage(2);
    n_out5.ResizeParentStorage(1); //outer node

    //Assigning Actual Parents To nodes - parent node index to pnode with element e1 in parent elem pointers
    n_in0.Assign(0,&e_in);
    n_in1.Assign(1,&e_in);
    n_in2.Assign(2,&e_in);
    n_out0.Assign(0,&e_out);
    n_out1.Assign(1,&e_out);
    n_out2.Assign(2,&e_out);

    n_in3.Assign(3,&e_in);
    n_in4.Assign(4,&e_in);
    n_in5.Assign(5,&e_in);
    n_out3.Assign(3,&e_out);
    n_out4.Assign(4,&e_out);
    n_out5.Assign(5,&e_out);

    ///Calibrating Neighbours
    //Assigning Neighbour elements to Nodes.
    //Inside
    n_in0.Assign(1, &e_nb_in1);
    n_in1.Assign(0, &e_nb_in1);
    n_in2.Assign(0, &e_nb_in2);
    //Outside
    n_out0.Assign(0, &e_nb_ou1);
    n_out1.Assign(0, &e_nb_ou2);
    n_out2.Assign(1, &e_nb_ou2);


    //Assigning Nodes to Element in InterFace
    //Inside
    e_in.Idx( 1 );
    e_in.Assign( 0, &n_in0 );
    e_in.Assign( 1, &n_in1 );
    e_in.Assign( 2, &n_in2 );
    //Outside
    e_out.Idx(2);
    e_out.Assign( 0, &n_out0 );
    e_out.Assign( 1, &n_out1 );
    e_out.Assign( 2, &n_out2 );

    //midside
    //inside
    e_in.Assign( 3, &n_in3 );
    e_in.Assign( 4, &n_in4 );
    e_in.Assign( 5, &n_in5 );
    //Outside
    e_out.Assign( 3, &n_out3 );
    e_out.Assign( 4, &n_out4 );
    e_out.Assign( 5, &n_out5 );

    //Assigning Neighbour element nodes
    //inside
    e_nb_in1.Idx(3);
    e_nb_in1.Assign(0, &n_in1);
    e_nb_in1.Assign(1, &n_in0);
    e_nb_in2.Idx(4);
    e_nb_in2.Assign(0, &n_in2);
    //outside
    e_nb_ou1.Idx(5);
    e_nb_ou1.Assign(0, &n_out0);
    e_nb_ou2.Idx(6);
    e_nb_ou2.Assign(0, &n_out1);
    e_nb_ou2.Assign(1, &n_out2);
    //Note, not all nodes of IsoparametricLinearTriangle are defined for neighbour elements
}












void InterFace_Test::ConstructInterFaceIngredients_3D(Element<3> &e_in, Element<3> &e_out,
                                                      Node<3> &n_in0, Node<3> &n_in1, Node<3> &n_in2, Node<3> &n_in3,
                                                      Node<3> &n_out0, Node<3> &n_out1, Node<3> &n_out2, Node<3> &n_out3,
                                                      Element<3> &e_nb_in1, Element<3> &e_nb_in2,
                                                      Element<3> &e_nb_ou1, Element<3> &e_nb_ou2)
{
  //setting node storage size
  n_in0.ResizeParentStorage(2); //outer node
  n_in1.ResizeParentStorage(2); //inner nodes...
  n_in2.ResizeParentStorage(2); //...
  n_in3.ResizeParentStorage(2);
  n_out0.ResizeParentStorage(2);//inner nodes
  n_out1.ResizeParentStorage(2);//...
  n_out2.ResizeParentStorage(2);
  n_out3.ResizeParentStorage(2); //outer node

  //Assigning Actual Parents To nodes - parent node index to pnode with element e1 in parent elem pointers
  n_in0.Assign(0,&e_in);
  n_in1.Assign(1,&e_in);
  n_in2.Assign(2,&e_in);
  n_in3.Assign(3,&e_in);
  n_out0.Assign(0,&e_out);
  n_out1.Assign(1,&e_out);
  n_out2.Assign(2,&e_out);
  n_out3.Assign(3,&e_out);

  //Assigning Neighbour elements to Nodes.
  //Inside
  n_in0.Assign(0, &e_nb_in1);
  n_in3.Assign(1, &e_nb_in1);
  n_in2.Assign(2, &e_nb_in1);
  n_in1.Assign(0, &e_nb_in2);
  n_in2.Assign(1, &e_nb_in2);
  n_in3.Assign(2, &e_nb_in2);
  //Outside
  n_out0.Assign(0, &e_nb_ou1);
  n_out1.Assign(1, &e_nb_ou1);
  n_out3.Assign(2, &e_nb_ou1);
  n_out0.Assign(0, &e_nb_ou2);
  n_out3.Assign(1, &e_nb_ou2);
  n_out2.Assign(2, &e_nb_ou2);

  //Assigning Nodes to Element in InterFace
  //Inside
  e_in.Idx( 1 );
  e_in.Assign( 0, &n_in0 );
  e_in.Assign( 1, &n_in1 );
  e_in.Assign( 2, &n_in2 );
  e_in.Assign( 3, &n_in3 );
  //Outside
  e_out.Idx(2);
  e_out.Assign( 0, &n_out0 );
  e_out.Assign( 1, &n_out1 );
  e_out.Assign( 2, &n_out2 );
  e_out.Assign( 3, &n_out3 );

  //Assigning Neighbour element nodes
  //inside
  e_nb_in1.Idx(3);
  e_nb_in1.Assign(0, &n_in0);
  e_nb_in1.Assign(1, &n_in3);
  e_nb_in1.Assign(2, &n_in2);
  e_nb_in2.Idx(4);
  e_nb_in2.Assign(0, &n_in1);
  e_nb_in2.Assign(1, &n_in2);
  e_nb_in2.Assign(2, &n_in3);

  //outside
  e_nb_ou1.Idx(5);
  e_nb_ou1.Assign(0, &n_out0);
  e_nb_ou1.Assign(1, &n_out1);
  e_nb_ou1.Assign(2, &n_out3);

  e_nb_ou2.Idx(6);
  e_nb_ou2.Assign(0, &n_out1);
  e_nb_ou2.Assign(1, &n_out3);
  e_nb_ou2.Assign(2, &n_out2);

}



















} //end namespace csmp





