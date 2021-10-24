#include "CsmpIntro1_TestCase.h"

#include "ANSYS_Model3D.h"
#include "VTU_Interface.h"
#include "Region.h"
#include "Boundary.h"

#include "Visitor.h"


using namespace std;

namespace csmp {

template<size_t dim>
class VolumeOutput_Visitor : public Visitor<dim> 
  {
  public:
    VolumeOutput_Visitor( Model<dim>& model );
    virtual ~VolumeOutput_Visitor() {}

    virtual void Visit(Element<dim>* element);

  private:
    Model<dim>& model_;
  };



  template<size_t dim>
  VolumeOutput_Visitor<dim>::VolumeOutput_Visitor( Model<dim>& model)
    : Visitor<dim>( MODEL, ELEMENT ), model_( model )
    {
    }


  template<size_t dim>
  void VolumeOutput_Visitor<dim>::Visit( Element<dim>* element )
    {
      Index volumeKey( model_.Database().StorageKey( "conductivity" ) );

      element->Store( volumeKey, makeScalar( PLAIN, element->Volume() ) );
    } // Visit()




void CsmpIntro1_TestCase::run()
  {
    // this is a 3D problem
    enum{DIM=3};

    // csmp variables
    // [ 1. ]
    // [ PLAIN ]
    const ScalarVariable oneScalar( PLAIN, 1. );
    // [ 1. ]
    // [ 1. ]
    // [ 1. ]
    // [ PLAIN ]
    // [ PLAIN ]
    // [ PLAIN ]
    const VectorVariable<DIM> oneVector( PLAIN, 1. );
    // [ 1.  1.  1. ]
    // [ 1.  1.  1. ]
    // [ 1.  1.  1. ]
    // [ PLAIN PLAIN PLAIN ]
    // [ PLAIN PLAIN PLAIN ]
    // [ PLAIN PLAIN PLAIN ]
    const TensorVariable<DIM> oneTensor( PLAIN, 1. );


    // establishing a box shaped csmp::Model providing icem mesh & regions file and variables file
    const string modelName( "BoxHalfs3D" );
    const string varFileName( "RSP-variables.txt" );
    ANSYS_Model3D model( modelName.data(), varFileName.data() );


    // testing for proper region instantiation
    _test( model.ContainsRegion("Model") );
    _test( model.ContainsRegion("MATRIX_LEFT") );
    _test( model.ContainsRegion("MATRIX_RIGHT") );
    _test( !model.ContainsRegion("MATRIX_right") );

    // testing for proper box boundary instantiation
    _test( model.ContainsBoundary("LEFT") );
    _test( model.ContainsBoundary("RIGHT") );
    _test( model.ContainsBoundary("TOP") );
    _test( model.ContainsBoundary("BOTTOM") );
    _test( model.ContainsBoundary("FRONT") );
    _test( model.ContainsBoundary("BACK") );
    _test( !model.ContainsBoundary("back") );

    // referencing regions
    Region<DIM>& modelRegion( model.Region("Model") );
    Region<DIM>& matrixLeft( model.Region("MATRIX_LEFT") );
    Region<DIM>& matrixRight( model.Region("MATRIX_RIGHT") );


    // referencing boundaries
    Boundary<DIM>& left( model.Boundary("LEFT") );
    Boundary<DIM>& right( model.Boundary("RIGHT") );
    Boundary<DIM>& top( model.Boundary("TOP") );
    Boundary<DIM>& front( model.Boundary("FRONT") );
    Boundary<DIM>& back( model.Boundary("BACK") );
    Boundary<DIM>& bottom( model.Boundary("BOTTOM") );


    // writing to a model (dispatches to its nodes)
    model.InputPropertyValue( "saturation oil", oneScalar );


    // VTU output interface
    VTU_Interface<DIM> vtu( model );
    if ( verbose_ )
       vtu.OutputDataToVTU( "OilSaturation1", "saturation oil", "Model", static_cast<int>(0) );
    
    // writing to a region & outputting
    matrixLeft.InputPropertyValue( "saturation oil", makeScalar( DIRICH, 20. ) );

    // new output
    if ( verbose_ ) {
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", "Model", static_cast<int>(0) );

        // outputting subdomains
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", left,        static_cast<int>(0) );
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", right,       static_cast<int>(0) );
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", bottom,      static_cast<int>(0) );
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", front,       static_cast<int>(0) );
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", back,        static_cast<int>(0) );
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", top,         static_cast<int>(0) );
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", matrixRight, static_cast<int>(0) );
        vtu.OutputDataToVTU( "OilSaturation2", "saturation oil", matrixLeft,  static_cast<int>(0) );
      }
    
    VolumeOutput_Visitor<DIM> volumeVisitor( model );
    model.Accept( volumeVisitor );

    if ( verbose_ ) {
        vtu.OutputDataToVTU( "VolumeVariable", "conductivity", matrixRight, static_cast<int>(0) );
        vtu.OutputDataToVTU( "VolumeVariable", "conductivity", matrixLeft,  static_cast<int>(0) );
      }
  }


} // csmp
