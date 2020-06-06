#ifndef CSMP_INTERFACE_TEST_H
#define CSMP_INTERFACE_TEST_H

#include "Test.h"
#include "InterFace.h"
#include "DenseMatrix.h"

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
    explicit InterFace_Test( bool verbose=false );

    ~InterFace_Test();
  
    virtual void run() override;
  
    void Assign_tests();
    void Geometry_tests();

    void ConstructInterFaceIngredients(  Element<2> &e_in, Element<2> &e_out,
                                         Node<2> &n_in0,  Node<2> &n_in1,  Node<2> &n_in2,
                                         Node<2> &n_out0, Node<2> &n_out1, Node<2> &n_out2,
                                         Element<2> &e_nb_in1, Element<2> &e_nb_in2,
                                         Element<2> &e_nb_ou1, Element<2> &e_nb_ou2 );

  private:
    const bool verbose_;
};

} //end csmp

#endif
