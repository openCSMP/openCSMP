#include "DenseMatrix_Example.h"

#include "DenseMatrix.h"

using namespace std;

namespace csmp{

void DenseMatrix_Example::Specifications()
{
  SetTitle( "Dense Matrix Operations" );
  SetDifficulty( 1 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "SKM" );
  AddDescription( "basic operations with csmp::DenseMatrix" );
  AddDescription( "source in: DenseMatrix_Example.cpp" );
}


/** *****************************************************************************************
    This program shows how to do simple matrix computations using the
    DenseMatrix class. The Following calculation is done:

 				Ax + b = c
 
    Where:
        A = [N x N] matrix
        x,b & c	= [N x 1] vectors

        N is dimension = 3
 
   **************************************************************************************** */
void DenseMatrix_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
    //ostream &cout = *GetStream();

    unsigned int N = 3;
    DenseMatrix<3>	A(N,N,0);
    DenseMatrix<3>  x(N,1,0),
                    b(N,1,0),
                    c(N,1,0);

    // Fill the matrix A and vectors x and b with (arbitrary) values
    A.FillCol( 0, 2. );	//		2.	0.	3.
    A.FillCol( 1, 3. );	//	A = 2.	0.	3.
    A.ZeroCol( 2 );		//		2.	0.	3.
    x.FillRow( 0, 1.);
    x.FillRow( 1, 3.);
    x.FillRow( 2, 2.);
    b.FillRow( 0, 6.);
    b.FillRow( 1, 3.);
    b.FillRow( 2, 5.);

    cout <<"\nExample: input matrix A:";
    A.Out();

    cout <<"\nExample: input vectors x and b: ";
    x.Out();
    b.Out();

    // Computation
    c  = A * x;
    c += b;

    // Output
    for ( size_t i=0; i<N; i++ )
      cout << "c = " << c(i,0) << endl;

    // or
    cout <<"\nExample: Result c = A * x + b";
    c.Out();

} // end Run

} // csmp
