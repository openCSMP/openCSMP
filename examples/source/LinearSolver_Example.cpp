#include "LinearSolver_Example.h"
#include "GaussJordan_Solver.h"

using namespace std;

namespace csmp {

void LinearSolver_Example::Specifications()
{
  SetTitle( "Linear System solution" );
  SetDifficulty( 1 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "Roman M.N." );
  AddDescription( "Solving simple linear equation: A*x=b " );
  AddDescription( "source in: LinearSolver_Example.cpp" );
}

void LinearSolver_Example
::SetupMatrix( SparseMatrix& A)
{
    // define Types
    const int BlockSize = 1;
    // make a block compressed row matrix with five point stencil
    const int BW2=31, N=BW2*BW2;
    const double factor( static_cast<double>(4+(BlockSize-1)) );
    A.Resize(N);
    for (size_t r = 0; r < N; r+=BlockSize )
    {
        // diagonal element
        for (size_t i=0; i<BlockSize; ++i)
            for (size_t j=0; j<BlockSize; ++j)
                if ( i == j )
                    A.Assign(r+i,r+j, factor );
                else
                    A.Assign(r+i,r+j,-1.0);
        //off-diagonal elements
        auto row=r/BW2;
        auto col=r%BW2;
        if (col-1>=0)
            for (size_t i=0; i<BlockSize; i++)
                A.Assign(r+i,r+i-1,-1.0);
        if (col+1<BW2)
            for (size_t i=0; i<BlockSize; i++)
                A.Assign(r+i,r+i+1,-1.0);
        if (row-1>=0)
            for (size_t i=0; i<BlockSize; i++)
                A.Assign(r+i,r+i-BW2,-1.0);
        if (row+1<BW2)
            for (size_t i=0; i<BlockSize; i++)
                A.Assign(r+i,r+i+BW2,-1.0);
    }

    std::vector<double> x(N,0.0);
    x[0]=1.0;
    x[N-1]=2.0;

    //cout<<"\nmatrix:\n";
    //A.Out();
    cout<<"\nsolution:\n";
    for (size_t i=0; i<N; i++)
        cout<< "x[ "<< i <<"] = "<< x[i] <<"\n";
    cout<<"\n";
}

void LinearSolver_Example
::SetupVectors( SparseMatrix&  A,
                std::vector<double>& b,
                std::vector<double>& x)
{
    // set up system:
    const int BW2=31, N=BW2*BW2;

    // prescribe known solution
    x.resize(N,0.0);
    x[0]=1.0;
    x[N-1]=2.0;

    // set right hand side accordingly
    b.resize(N,0.0);
    A.MultiplyWith(x,b);

    // initial guess
    fill( x.begin(),x.end(), 1.0 ) ;
    for (size_t i=0; i<N; i++)
      x[i] = static_cast<double>(i)*0.1;

    cout<<"\nrhs:\n";
    for (size_t i=0; i<N; i++)
        cout<< "b[ "<< i <<"] = "<< b[i] <<"\n";
    cout<<"\n";
}

#if defined(CSMP_WITH_DUNE_ISTL)

void LinearSolver_Example
::SetupMatrix( Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& A)
{
    // define Types
    const int BlockSize = 1;
    typedef Dune::FieldVector<double,BlockSize> VB;
    typedef Dune::BlockVector<VB> Vector;
    typedef Dune::FieldMatrix<double,BlockSize,BlockSize> MB;
    typedef Dune::BCRSMatrix<MB>  Matrix;

    // build little blocks
    MB D=0;
    for (int i=0; i<BlockSize; i++)
      for (int j=0; j<BlockSize; j++)
        if (i==j)
            D[i][j] = 4+(BlockSize-1);
        else
            D[i][j] = -1;

    MB E=0;
    for (int i=0; i<BlockSize; i++)
      E[i][i] = -1;

    // make a block compressed row matrix with five point stencil
    const int BW2=31, N=BW2*BW2;
    A.setSize( N,N,5*N );
    A.setBuildMode(Dune::BCRSMatrix<MB>::row_wise);
    for (Matrix::CreateIterator i=A.createbegin(); i!=A.createend(); ++i)
    {
      int row=i.index()/BW2;
      int col=i.index()%BW2;
      i.insert(i.index());
      if (col-1>=0)  i.insert(i.index()-1);
      if (col+1<BW2) i.insert(i.index()+1);
      if (row-1>=0)  i.insert(i.index()-BW2);
      if (row+1<BW2) i.insert(i.index()+BW2);
    }
    for (Matrix::RowIterator i=A.begin(); i!=A.end(); ++i)
      for (Matrix::ColIterator j=(*i).begin(); j!=(*i).end(); ++j)
        if ( i.index()==j.index() )
          (*j) = D;
        else
          (*j) = E;

    // Output Matrix content
    //printmatrix(std::cout,A,"system matrix","row",10,2);

    Vector x(N);
    x=0;
    x[0]=1;
    x[N-1]=2;

    // Output Solution content
    printvector(std::cout,x,"solution","x");

}

void LinearSolver_Example
::SetupVectors( const Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& A,
                Dune::BlockVector<Dune::FieldVector<double,1> >& b,
                Dune::BlockVector<Dune::FieldVector<double,1> >& x)
{
    // define Types
    const int BlockSize = 1;
    typedef Dune::FieldVector<double,BlockSize> VB;
    typedef Dune::BlockVector<VB> Vector;

    // set up system:
    const int BW2=31, N=BW2*BW2;

    // prescribe known solution
    x.resize(N);
    x=0;
    x[0]=1;
    x[N-1]=2;

    // set right hand side accordingly
    b.resize(N);
    b=0;
    A.umv(x,b);

    // initial guess
    x=1;
    for (int i=0; i<N; i++)
      x[i] = i*0.1;

    // Output Solution content
    printvector(std::cout,b,"rhs","b");
}

#endif


void LinearSolver_Example::Run()
{
    // setup system
    SparseMatrix A;
    std::vector<double> b,x;
    const int nunknowns(1);
    SetupMatrix(A);
    SetupVectors(A,b,x);

    #if !defined(CSMP_WITH_DUNE_ISTL)

    GaussJordan_Solver solver;

    solver.Solve(A,b,x,nunknowns);

    #elif defined(CSMP_WITH_DUNE_ISTL)

    int choice;
    cout<<"\nThe available solvers for the system A*x=b are:\n";
    cout<<"01) Jacobi Solver\n";
    cout<<"02) ILU0 Solver\n";
    cout<<"03) ILUn Solver\n";
    cout<<"04) MINRES Solver with SSOR Preconditioner\n";
    cout<<"05) CG Solver with Jacobian Preconditioner\n";
    cout<<"06) CG Solver with SSOR Preconditioner\n";
    cout<<"07) CG Solver with ILU0 Preconditioner\n";
    cout<<"08) CG Solver with ILUn Preconditioner\n";
    cout<<"09) BiCGSTAB Solver with Jacobian Preconditioner\n";
    cout<<"10) BiCGSTAB Solver with SSOR Preconditioner\n";
    cout<<"11) BiCGSTAB Solver with ILU0 Preconditioner\n";
    cout<<"12) BiCGSTAB Solver with ILUn Preconditioner\n";
    cout<<"13) Iterative Solver with AMG Preconditioner smoothed by SSOR\n";
    cout<<"14) Iterative Solver with AMG Preconditioner smoothed by SOR\n";
    cout<<"15) CG Solver with AMG Preconditioner smoothed by SSOR\n";
    cout<<"16) BiCGSTAB Solver with AMG Preconditioner smoothed by SSOR\n";
    cout<<"17) BiCGSTAB Solver with AMG Preconditioner smoothed by SOR\n";
    cout<<"Please choose the Solver you would like to examine:";
    cin >> choice;

    // define Types
    const int BlockSize (1);
    typedef Dune::FieldVector<double,BlockSize> VB;
    typedef Dune::FieldMatrix<double,BlockSize,BlockSize> MB;
    typedef Dune::BlockVector<VB> Vector;
    typedef Dune::BCRSMatrix<MB>  Matrix;

    // setup system
    //Matrix A;
    //Vector b,x;
    //SetupMatrix(A);
    //SetupVectors(A,b,x);

    DuneISTL_Settings  sets;
    sets.Set_reduction( 1.0e-15 );
    sets.Set_usesuperlu(false);
    if( choice == -2 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_ExplicitDiagonal> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    if( choice == -1 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_Base<Dune::SeqJac,Dune::LoopSolver> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 1 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_LOOP_Jac> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 2 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_ILU0<Dune::LoopSolver> > solver;
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 3 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_ILUn<Dune::LoopSolver> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 4 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_MINRES_SSOR> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 5 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_CG_Jac> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 6 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_CG_SSOR> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 7 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_CG_ILU0> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 8 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_CG_ILUn> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 9 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_BCGS_Jac> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 10 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_BCGS_SSOR> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 11 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_BCGS_ILU0> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 12 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_BCGS_ILUn> solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 13 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_LS_AMG_SSOR<3U,Matrix,Vector> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 14 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_LS_AMG_SOR<3U,Matrix,Vector> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 15 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_CG_AMG_SSOR<3U,Matrix,Vector> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 16 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_BCGS_AMG_SSOR<3U,Matrix,Vector> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else if( choice == 17 )
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_BCGS_AMG_SOR<3U,Matrix,Vector> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    else
    {
        DuneISTL_Solver<Matrix,Vector,Vector,DuneISTL_SEQ_AMG<3U,Matrix,Vector,Dune::SeqSOR,Dune::BiCGSTABSolver> > solver(&sets);
        solver.Solve(A,b,x,nunknowns);
    }
    //printvector(std::cout,x,"numerical solution","r");
    #endif

    cout<<"\nnumerical solution:\n";
    const size_t N = b.size();
    for (size_t i=0; i<N; i++)
        cout<< "r[ "<< i <<"] = "<< x[i]<< "\n";
    cout<<"\n";
    
} // end Run

} // csmp
