#include "Face_Test.h"
#include "Model.h"
#include "IsoparametricLinearTetrahedron.h"
#include "FiniteVolumeStencil.h"
#include "variableOperations.h"

namespace csmp
{

  template<size_t dim>
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
    // creating a plain face
	FiniteElement* feptr= new IsoparametricLinearTetrahedron( 4 );
	FiniteVolumeStencil<DIM>* fvptr = new FiniteVolumeStencil<DIM>( "ISOPARAMETRIC_LINEAR_TETRAHEDRON" );

    size_t face_number{1};
    Face<DIM> face( face_number, feptr, fvptr,
                    LocalVariables( 1,0,0,0,0,0,0,1,1 ),
                    IntegrationPointVariables() );

    // read/write
    Index scalarKey( SCALAR, FACE, 0 );
    const double scalarValue( 1.1 ); const double tolerance( 1.0E-5 );
    ScalarVariable scalarVariable( PLAIN, scalarValue );
    face.Store( scalarKey, scalarVariable );
    _equal( face.Read( scalarKey ), scalarValue, tolerance );

    // creating 5 nodes which will set up 2 linear tetrahedrons
    Node<DIM> n1,n2,n3,n4,n5;
    n1.x( 2. ); n1.y( 7. ); n1.z( 10. );                     

    _test( true );

	delete feptr;
	delete fvptr;

    return;
  }


  bool Face_Test::FaceUnitNormalPointsOutward( const Face<3>& f )
    {
      return faceUnitNormalPointsOutward(f);
    }


} // csmp
