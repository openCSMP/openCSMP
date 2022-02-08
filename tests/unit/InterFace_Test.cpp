#include "InterFace_Test.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearLineElement.h"
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
    size_t faceID_of_e1 = 0;
    size_t faceID_of_e2 = 1;

    //Interface Construction
    // ---------------------
    // simplemost
    IsoparametricLinearLineElement face_FE;
    FiniteVolumeStencil<2>         line_stencil("ISOPARAMETRIC_LINEAR_BAR");
    InterFace<2> if_obj0( &face_FE, &line_stencil ), if_obj1( &face_FE, &line_stencil ),
                 if_obj2( &face_FE, &line_stencil ), if_obj3( &face_FE, &line_stencil );


    ///Beginning Use and Tests
    /// -----------------------------------------------------------------
    ///Assign 0
    //Only assigns element
    if_obj0.Assign(&e1, &e2,  false);
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
    if_obj2.Assign(&e1, &e2, true) ;
    //Test Element sides
    _test( e1 == *(if_obj2.InnerParent()) );
    _test( e2 == *(if_obj2.OuterParent()) );
    //test node assigned
    _test( n2 ==  *(if_obj2.InnerParent()->N( if_obj2.ParentNodeNumber(size_t(0), INSIDE) ) ) );
    _test( n3 ==  *(if_obj2.InnerParent()->N( if_obj2.ParentNodeNumber(size_t(1), INSIDE) ) ) );
    _test( n6 ==  *(if_obj2.OuterParent()->N( if_obj2.ParentNodeNumber(size_t(2), OUTSIDE) ) ) );
    _test( n4 ==  *(if_obj2.OuterParent()->N( if_obj2.ParentNodeNumber(size_t(3), OUTSIDE) ) ) );


    ///Assign 3
    if_obj3.Assign(&e1, faceID_of_e1, &e2, faceID_of_e2, true);
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
    n6.x(2.0); n6.y(3.0);//topologically collocated with n3

    ConstructInterFaceIngredients(e1, e2, n1, n2, n3, n4, n5, n6, e_nb_in1, e_nb_in2, e_nb_ou1, e_nb_ou2);
    //Knowledge of face ID -- a result of this particular constrution
    size_t faceID_of_e1 = 0;
    size_t faceID_of_e2 = 1;

    //Construc InterFace objects
    IsoparametricLinearLineElement face_FE;
    FiniteVolumeStencil<2>         line_stencil("ISOPARAMETRIC_LINEAR_BAR");
    InterFace<2> if_obj0( &face_FE, &line_stencil ), if_obj1( &face_FE, &line_stencil ), if_obj2( &face_FE, &line_stencil );

    //Assign 0
    if_obj0.Assign(&e1, &e2, true);
    //Assign 1
    if_obj1.Assign(&e1, faceID_of_e1, INSIDE);
    if_obj1.Assign(&e2, faceID_of_e2, OUTSIDE);
    //assigns nodes manually
    if_obj1.Assign(0, &n2, INSIDE);
    if_obj1.Assign(1, &n3, INSIDE);
    if_obj1.Assign(0, &n6, OUTSIDE);
    if_obj1.Assign(1, &n4, OUTSIDE);
    //Assign 2
    if_obj2.Assign(&e1, faceID_of_e1, &e2, faceID_of_e2, true);


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
    VectorVariable<2> nrml;
    VectorVariable<2> n_inside, n_outside, n_middle;
    n_inside.Component(0,  1.0); //True unit normal to
    n_inside.Component(1,  0.0);
    n_outside.Component(0,-1.0);
    n_outside.Component(1, 0.0);
    n_middle = n_inside;

    //starting test
    if_obj0.UnitNormal(nrml, INSIDE);
    _test(nrml == n_inside);
    if_obj1.UnitNormal(nrml, INSIDE);
    _test(nrml == n_inside);
    if_obj2.UnitNormal(nrml, INSIDE);
    _test(nrml == n_inside);
    if_obj0.UnitNormal(nrml, OUTSIDE);
    _test(nrml == n_outside);
    if_obj1.UnitNormal(nrml, OUTSIDE);
    _test(nrml == n_outside);
    if_obj2.UnitNormal(nrml, OUTSIDE);
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

    _test( if_obj0.LengthInDirection(nrml) == 0.0);
    _test( if_obj1.LengthInDirection(nrml) == 0.0);
    _test( if_obj2.LengthInDirection(nrml) == 0.0);


    //Testing Coordinate matrix -- WHICH ONE???
    DenseMatrix<DM_MIN> XY_inside(2, 2);
    XY_inside.AssignRow(0, n2.Coordinate());
    XY_inside.AssignRow(1, n3.Coordinate());

    DenseMatrix<DM_MIN> XY_outside(2,2);
    XY_outside.AssignRow(0, n6.Coordinate());
    XY_outside.AssignRow(1, n4.Coordinate());

    DenseMatrix<DM_MIN> XY(2,2);
    if_obj0.NodeCoordinateMatrix( XY , INSIDE);
    //_test( XY_inside == XY);
    if_obj0.NodeCoordinateMatrix( XY, OUTSIDE);

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



} //end namespace csmp





