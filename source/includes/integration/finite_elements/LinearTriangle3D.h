#ifndef LINEAR_TRIANGLE_3D_H
#define LINEAR_TRIANGLE_3D_H

#include "FiniteElement.h"

namespace csmp {

/**
 
@brief Straight-sided, analytically integrated triangular finite element 
suspended in 3D space; linear element interpolation functions.

@date 2000
@author S.K. Matthaei
@author Stephen G. Roberts
@author S. Geiger

 
@section motivation Motivation

To carry out computations using linear triangular and tetrahedral 
elements with exact integration in which surfaces inside the tetrahedral
mesh are represented by triangles oriented freely in space. This triangle
has to be put into the FiniteElementManager in
order to be accessible for computations.  
 

@section collaboration Collaboration

The LinearTriangle3D collaborates with the 'Element' class within 
a bridge pattern proving polymorphic behaviour of the 'Element'.
 

@section implementation Implementation

See CSMP User's guide for a full listing of shape functions and their
derivatives.  
 
 
@section examples Application Examples

The shape function values and their derivatives can be output to VTK
by the public method 'OutputToVTK()'. This result can be visualized
with the following VTK tcl script.  

@code
catch {load vtktcl}
source vtkInt.tcl

vtkUnstructuredGridReader reader
    reader SetFileName    "etest1.vtk"
    reader SetScalarsName "sum_N"
    reader SetVectorsName "derivative_N"

vtkDataSetMapper gridMapper
    gridMapper SetInput [reader GetOutput]
    gridMapper ScalarVisibilityOff
    
vtkActor gridActor
    gridActor SetMapper gridMapper
    [gridActor GetProperty] SetColor 0.1 0.1 0.1
    [gridActor GetProperty] SetRepresentationToWireframe

vtkDataSetMapper modelMapper
    modelMapper SetInput [reader GetOutput]
    modelMapper ScalarVisibilityOn
    reader Update; #force update for scalar range
    eval modelMapper SetScalarRange [[reader GetOutput] GetScalarRange]
    # red=high, blue=low
    [modelMapper GetLookupTable] SetHueRange .667 0.
    
vtkActor modelActor
    modelActor SetMapper modelMapper
    [modelActor GetProperty] SetColor .2 .2 .2
    [modelActor GetProperty] SetOpacity 0.7

 vtkHedgeHog hhog1
    hhog1 SetInput [reader GetOutput]
    hhog1 SetScaleFactor 2.0 

 vtkPolyDataMapper uvMapper
    uvMapper SetInput [hhog1 GetOutput]
 vtkActor uvActor
    uvActor SetMapper uvMapper
    [uvActor GetProperty] SetColor 0. 0. 0.

vtkAxes axes
    axes SetOrigin 0. 0. 0.
    axes SetScaleFactor 10.
vtkPolyDataMapper axesMapper
    axesMapper SetInput [axes GetOutput]
    axesMapper ScalarVisibilityOff
    
vtkActor axesActor
    axesActor SetMapper axesMapper
    [axesActor GetProperty] SetDiffuseColor 0.1 0.1 0.1
    
vtkVectorText textX
    textX SetText "X"
vtkPolyDataMapper textMapperX
    textMapperX SetInput [textX GetOutput]
    
vtkFollower textActorX
    textActorX SetMapper textMapperX
    textActorX SetScale .5 .5 .5
    textActorX AddPosition 11. 0. 0.
    [textActorX GetProperty] SetDiffuseColor 0.1 0.1 0.1
   
vtkVectorText textY
    textY SetText "Y"
vtkPolyDataMapper textMapperY
    textMapperY SetInput [textY GetOutput]
    
vtkFollower textActorY
    textActorY SetMapper textMapperY
    textActorY SetScale .5 .5 .5 
    textActorY AddPosition 0. 11. 0
    [textActorY GetProperty] SetDiffuseColor 0.1 0.1 0.1
    
vtkVectorText textZ
    textZ SetText "Z"
vtkPolyDataMapper textMapperZ
    textMapperZ SetInput [textZ GetOutput]
    
vtkFollower textActorZ
    textActorZ SetMapper textMapperZ
    textActorZ SetScale .5 .5 .5 
    textActorZ AddPosition 0. 0. 11.
    [textActorZ GetProperty] SetDiffuseColor 0.1 0.1 0.1

vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

ren1 AddActor modelActor
ren1 AddActor gridActor
ren1 AddActor uvActor

ren1 AddActor axesActor
ren1 AddActor textActorX
ren1 AddActor textActorY
ren1 AddActor textActorZ

ren1 SetBackground 1 1 1
renWin SetSize 900 700

iren Initialize
iren SetUserMethod {wm deiconify .vtkInteract}

renWin SetFileName "element.ppm"
renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .
@endcode
 
*/
class LinearTriangle3D : public FiniteElement {
  public:
    LinearTriangle3D();
    ~LinearTriangle3D();

    virtual double       Volume();
    virtual double       AspectRatio();
    virtual double       InnerRadius();
    virtual void           EdgeLengths( std::vector<double>& vec );
    virtual void           NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void           NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual void           CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual size_t         CornerNodes() const  { return 3; }
    virtual void           CornerNodes( std::vector<size_t>& ids ) const;
    virtual void           UnitNormal( std::vector<double>& unrml ) const;
    virtual void           UnitNormalToFace( size_t face, std::vector<double>& unrml ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return LINEAR_BAR; };

    /// element interpolation functions at the point xzy
    virtual void           N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void           N_AtBaryCenter( std::vector<double>& N );

    /// first derivatives of element interpolation functions
    virtual void           dN( DenseMatrix<DM_MIN>& );
    /// note that the derivative is constant so that the point location will not be considerd
    virtual   double     dN_At( DenseMatrix<DM_MIN>&, const std::vector<double>& xyz );
    virtual   double     dN_AtBarycenter( DenseMatrix<DM_MIN>& );

    /// integral over interpolation functions
    virtual void           IntegralNN( DenseMatrix<DM_MIN>& M );
    
    virtual void           OutputNodeDataToVTK( const char* file_name,
                                                const char* var_name,
                                                DenseMatrix<DM_MIN>& DATA ) const;

    void                   OutputToVTK( const char* file_name );
};


} // end namespace csmp

#endif




