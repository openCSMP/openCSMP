#include "Face_Test.h"
#include "Face.h"
#include "Element.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearTetrahedron.h"
#include "FiniteVolumeStencil.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "variableOperations.h"
#include "compareFloats.h"

using namespace std;

namespace csmp
{

  template<uint32_t dim>
  bool faceUnitNormalPointsOutward( const Face<dim>& f )
    {
      VectorVariable<dim> faceUN( PLAIN, 0. ), faceToInner( PLAIN, 0. );
      f.UnitNormal(faceUN);
      const Point<dim> bcFace( f.BaryCenter() ), bcInner( f.Parent(INSIDE)->BaryCenter() );
      for( size_t d(0); d < dim; ++d )
        faceToInner(d) = bcInner[d] - bcFace[d];
      if( dotProduct(faceToInner,faceUN) > 0. )
        return false;
      return true;
    }




void Face_Test::run()
  {
    // 0. some basics
    // --------------
    // creating a plain face
  	FiniteElement* feptr= new IsoparametricLinearTriangle( 4 );
	  FiniteVolumeStencil<DIM>* fvptr = new FiniteVolumeStencil<DIM>( "ISOPARAMETRIC_LINEAR_TRIANGLE" );

    // creating a face with some scalar variables
    const size_t face_number{1};
    Face<DIM> face( face_number, feptr, fvptr,
                    LocalVariables( 1,0,0,0,0,0,0,1,1 ),
                    IntegrationPointVariables() );

    // testing read/write of variables
    Index scalarKey( SCALAR, FACE, 0 );
    const double scalarValue( 1.1 ); const double tolerance( 1.0E-5 );
    ScalarVariable scalarVariable( PLAIN, scalarValue );
    face.Store( scalarKey, scalarVariable );
    _equal( face.Read( scalarKey ), scalarValue, tolerance );

    // connecting the face to its nodes
    // creating 4 nodes
    Node<DIM> n0, n1, n2, n3, n4;
    n0.x( 1.5 ); n0.y( 0. ); n0.z(  1. );
    n1.x( 1.5 ); n1.y( 0. ); n1.z( -0.5 );
    n2.x( 0.0 ); n2.y( 0. ); n2.z(  0. );
    n3.x( 2.0 ); n3.y( 1. ); n3.z(  0. );
    n4.x( 4.0 ); n4.y( 0. ); n4.z(  0. );
    n0.Idx(0); n1.Idx(1); n2.Idx(2); n3.Idx(3); n4.Idx(4);
    face.Assign( 0, &n0 );
    face.Assign( 1, &n1 );
    face.Assign( 2, &n3 );

    // creating 2 linear tetrahedra (1 on each side of the face, 0 on inside)
    IsoparametricLinearTetrahedron tetra( 4 );
    FiniteVolumeStencil<DIM> tet_stencil( "ISOPARAMETRIC_LINEAR_TETRAHEDRON" );
    Element<DIM>  tet1( &tetra, &tet_stencil ), tet2( &tetra, &tet_stencil );
    tet1.Idx(0);
    tet1.Assign( 0, &n0 );
    tet1.Assign( 1, &n1 );
    tet1.Assign( 2, &n2 );
    tet1.Assign( 3, &n3 );
    _test( tet1.Volume() > 0.1 );
    tet2.Idx(2);
    tet2.Assign( 0, &n0 );
    tet2.Assign( 1, &n4 );
    tet2.Assign( 2, &n1 );
    tet2.Assign( 3, &n3 );
    _test( tet2.Volume() > 0.1 );
    // neighbors
    tet1.Assign( 2, &tet2 );
    tet2.Assign( 1, &tet1 );
    // face higher-dimensional neighbors (Assign function)
    face.Assign( &tet1, &tet2 );
    _test( face.InnerParentFaceID() == 2 );
    _test( face.OuterParentFaceID() == 1 );
    
    // testing face normal( test requires higher-dimensional neighbors)
    _test( FaceUnitNormalPointsOutward( face ) );
    

    // 1. creating more Face objects from the same data, but with different constructors
    // ---------------------------------------------------------------------------------
    // are they the same?                     lin.intpol, isoparametric
    FiniteElementManager             fem_manager( DIM, 1, true );
    FiniteVolumeStencilManager<DIM>  fvm_manager( fem_manager );
    
    // constructor with auto-detection of shared face
    Face<DIM> face2( fem_manager, fvm_manager,
                     &tet1, &tet2,
                     LocalVariables( 1,0,0,0,0,0,0,1,1 ),
                     IntegrationPointVariables() );
    _test( face2 == face );
    
    // construction from a lower-dimensional element
    Element<DIM>  triangle( feptr, fvptr,
                            LocalVariables( 1,0,0,0,0,0,0,1,1 ),
                            IntegrationPointVariables() );
    triangle.Assign( 0, &n0 );
    triangle.Assign( 1, &n1 );
    triangle.Assign( 2, &n3 );
    // from lower-dim element Face constructor
    const size_t inner_parent_face_id{2}, outer_parent_face_id{1};
    Face<DIM> face3( triangle, &tet1, &tet2, inner_parent_face_id, outer_parent_face_id,
                     LocalVariables( 1,0,0,0,0,0,0,1,1 ), IntegrationPointVariables() );
    _test( face3 == face );
    
    // test construction of an Face on face 0 of element tet1
    const size_t boundary_face{0};
    Face<DIM> face4( tet1, fem_manager.E( ISOPARAMETRIC_LINEAR_TRIANGLE ),
                     fvm_manager, boundary_face,
                     LocalVariables( 1,0,0,0,0,0,0,1,1 ), IntegrationPointVariables() );
    _test( face4.N(0) == &n1 );
    _test( face4.N(1) == &n2 );
    _test( face4.N(2) == &n3 );
    _test( face4.InnerParent() == &tet1 );
    
    // TODO: test this edge constructor
    /*
    template<uint32_t dim>
    Face<dim>::Face( csmp::FiniteElement* FE_type_of_boundary_face,
                     const FiniteVolumeStencilManager<dim>& fvm_manager,
                     Element<dim>* const parent_of_face1,
                     Element<dim>* const parent_of_face2,
                     size_t parent_elmt1_segm_id,
                     size_t parent_elmt2_segm_id,
                     const std::vector<Node<dim>*>&  edge_nodes,
                     const LocalVariables& ep,
                     const IntegrationPointVariables& ip )
    */

    // 2. copy construction, assignment  and move operations
    // ---------------------------------------------------------------------------------
    Face<DIM> face5( face4 );
    _test( face5 == face4 );

    Face<DIM> face6 = face3;
    _test( face6 == face3 );
    
    // move operations
    // copy
    Face<DIM> face7( move(face5) );
    _test( face7 == face4 );
    // assignment
    face5 = move(face6);
    _test( face5 == face3 );
    
    
    // 3. Member functions of Face
    // ---------------------------------------------------------------------------------
    // assigning neighbor faces
    face3.Assign( 0, &face );
    face3.Assign( 1, &face2 );
    face3.Assign( 2, &face2 );
    _test( face3.ConnectedNeighbors() == 3 );
    // removing a specific neighbor face that exists multiple times in the Face neighbor array.
    face3.Unassign( &face2 );
    face3.Unassign( &face2 );
    _test( face3.Neighbor(1) == nullptr );
    _test( face3.Neighbor(2) == nullptr );
    
    // center of gravity
    Point<DIM>     bctr = face3.BaryCenter();
    vector<double> N;
    face3.N_AtGlobalPoint( N, bctr.Coordinates() );
    _test( N[0] > 0. ); // is inside
    _test( N[1] > 0. ); // is inside
    _test( N[2] > 0. ); // is inside
    
    double length = face.LengthInDirection( VectorVariable<DIM>{ANY,ANY,ANY,0.,0.,1.} );
    _test( approximatelyEqual(length,1.5) );

    // cleanup
	  delete feptr;
	  delete fvptr;

    return;
  }



bool Face_Test::FaceUnitNormalPointsOutward( const Face<3>& f )
  {
    return faceUnitNormalPointsOutward(f);
  }


} // csmp
