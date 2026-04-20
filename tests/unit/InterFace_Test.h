#ifndef CSMP_INTERFACE_TEST_H
#define CSMP_INTERFACE_TEST_H

#include "Test.h"
#include "InterFace.h"
#include "DenseMatrix.h"
#include "Element.h"

namespace csmp
{

/**
   Tests the construction of InterFace via its constructors and reconstructors
       as well as the functionality of the interface.
       
       @todo only 2D so far; expand to 3D. 
*/
class InterFace_Test : public Test
{

  public:
    InterFace_Test();

    ~InterFace_Test();
  
    virtual void run() override;
  
    void Assign_tests();
    void Geometry_tests();

    void NodeCoordinateMatrix_linear_test();
    void Assign_Geometry_quadratic_test();
    void Assign_Geometry_quadratic2_test();

    void Assign_tests_linear_3D();





    ///Data for test

    void ConstructInterFaceIngredients(  Element<2> &e_in, Element<2> &e_out,
                                         Node<2> &n_in0,  Node<2> &n_in1,  Node<2> &n_in2,
                                         Node<2> &n_out0, Node<2> &n_out1, Node<2> &n_out2,
                                         Element<2> &e_nb_in1, Element<2> &e_nb_in2,
                                         Element<2> &e_nb_ou1, Element<2> &e_nb_ou2 );

    void ConstructInterFaceIngredients_quadratic(  Element<2> &e_in, Element<2> &e_out,
                                         Node<2> &n_in0,  Node<2> &n_in1,  Node<2> &n_in2,
                                         Node<2> &n_out0, Node<2> &n_out1, Node<2> &n_out2,
                                         Element<2> &e_nb_in1, Element<2> &e_nb_in2,
                                         Element<2> &e_nb_ou1, Element<2> &e_nb_ou2,
                                         Node<2> &n_in3,  Node<2> &n_in4,  Node<2> &n_in5,
                                         Node<2> &n_out3, Node<2> &n_out4, Node<2> &n_out5);



    void ConstructInterFaceIngredients_3D( Element<3> &e_in, Element<3> &e_out,
                                               Node<3> &n_in0,  Node<3> &n_in1,  Node<3> &n_in2, Node<3> &n_in4,
                                               Node<3> &n_out0, Node<3> &n_out1, Node<3> &n_out2, Node<3> &n_out4,
                                               Element<3> &e_nb_in1, Element<3> &e_nb_in2,
                                               Element<3> &e_nb_ou1, Element<3> &e_nb_ou2 );

  private:
    const static bool verbose_ = true;
};

} //end csmp

#endif
