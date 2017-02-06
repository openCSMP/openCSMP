#ifndef CSMP_FINITE_ELEMENT_H
#define CSMP_FINITE_ELEMENT_H

#include "DenseMatrix.h"

namespace csmp {
/**
@file FiniteElement.h

@author S.K. Matthai
@date 1996

@todo  (3) Make FiniteElements singletons
@todo  (3) Manage elements in a more generic way 

@addtogroup CSMPglobalEnums
@{
*/


//                  CSMP TYPE                                             ANSYS - ICEMTYPE
//                  =========                                             ================

enum CSMP_FEM_TYPE{ UNKNOWN,
                    LINEAR_BAR,    										                  // BAR_2            = 2,
                    QUADRATIC_BAR, 										                  // BAR_3            = 3,
                    CUBIC_BAR, 										                      // BAR_4  
                    LINEAR_TRIANGLE,   
                    LINEAR_TRIANGLE3D, 									                // TRI_3            = 8,
                    BARYCENTRIC_LINEAR_TRIANGLE, 			    	           	// TRI_3_X          = 9,
                    QUADRATIC_TRIANGLE, 								                // 2D & 3D TRI_6    = 10,
                    BARYCENTRIC_QUADRATIC_TRIANGLE, 					          // TRI_6_X          = 11,
                    CUBIC_TRIANGLE,
                    LINEAR_TETRAHEDRON, 								                // TETRA_4          = 4,
                    QUADRATIC_TETRAHEDRON, 								              // TETRA_10         = 5,
                    BARYCENTRIC_QUADRATIC_TETRAHEDRON,                  // PYRA_5           = 18,
                    CUBIC_TETRAHEDRON,
                    ISOPARAMETRIC_LINEAR_BAR,    						            // BAR_2            = 2,
                    ISOPARAMETRIC_QUADRATIC_BAR, 						            // BAR_3            = 3,
                    ISOPARAMETRIC_CUBIC_BAR, 
                    ISOPARAMETRIC_LINEAR_TRIANGLE,  					          // TRI_3            = 8,
                    ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE, 		      // TRI_3_X          = 9,
                    ISOPARAMETRIC_QUADRATIC_TRIANGLE, 					        // 2D & 3D TRI_6    = 10,
                    ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE, 		  // TRI_6_X          = 11,
                    ISOPARAMETRIC_CUBIC_TRIANGLE,
                    ISOPARAMETRIC_LINEAR_TETRAHEDRON, 					        // TETRA_4          = 4,
                    ISOPARAMETRIC_QUADRATIC_TETRAHEDRON, 			          // TETRA_10         = 5,
                    ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON,    // TETRA_11
                    ISOPARAMETRIC_CUBIC_TETRAHEDRON,
                    ISOPARAMETRIC_LINEAR_PYRAMID,						            // PYRA_5           = 18,
                    ISOPARAMETRIC_QUADRATIC_PYRAMID13, 					        // PYRA_13          = 24,
                    ISOPARAMETRIC_QUADRATIC_PYRAMID14,     				      // PYRA_14          = 22,
                    ISOPARAMETRIC_CUBIC_PYRAMID,
                    ISOPARAMETRIC_LINEAR_PRISM,  						            // PENTA_6          = 12,
                    ISOPARAMETRIC_QUADRATIC_PRISM15,         			      // PENTA_15         = 13,
                    ISOPARAMETRIC_QUADRATIC_PRISM18,           			    // PENTA_18         = 21,
                    ISOPARAMETRIC_CUBIC_PRISM,
                    ISOPARAMETRIC_LINEAR_QUADRILATERAL,                 // QUAD_4           = 14,
                    ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL,     // QUAD_4_X         = 15,
                    ISOPARAMETRIC_QUADRATIC_QUADRILATERAL,              // QUAD_8           = 16,
                    ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL,  // QUAD_8_X         = 17,
                    ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9,  			      // QUAD_9           = 19,
                    ISOPARAMETRIC_CUBIC_QUADRILATERAL,
                    ISOPARAMETRIC_LINEAR_HEXAHEDRON,              		  // HEXA_8           = 6,
                    ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20,       		    // HEXA_20          = 7,
                    ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27, 				      // HEXA_27          = XX,
                    ISOPARAMETRIC_CUBIC_HEXAHEDRON,
                    ZERO_DIMENSIONAL_FACE,                              // the face of a line element 
                    POINT_ELEMENT,
                    POLYGONAL_ELEMENT,
                    POLYHEDRAL_ELEMENT,
                    EXPERIMENTAL_ELEMENT };

enum ELEMENT_DIMENSION { LINE=1, SURFACE=2, VOLUME=3, HYPER_DIMENSIONAL=4 };

/**
@}
*/

bool isTriangularElement( CSMP_FEM_TYPE );
bool isQuadrilateralElement( CSMP_FEM_TYPE );
bool isLineElement( CSMP_FEM_TYPE );

ELEMENT_DIMENSION  parseFiniteElementDimension( CSMP_FEM_TYPE );
CSMP_FEM_TYPE      parseFiniteElementTypeEnum( int32 csp_etype );
CSMP_FEM_TYPE      parseFiniteElementType( const std::string& etype );
const char*        parseFiniteElementType( int32 etype );


/**

@brief CSMP's base class of all finite elements whose functionality 
is provided through polymorphism.

@author S.K. Matthai
@author Stephen G. Roberts
@date 1996
 
@section motivatio Motivation

Base class for CSMP's finite element library.  
 
@section design Design Intent

Provide a unified interface for all finite element calculations
carried out by CSMP.  

Store the matrices necessary for finite element computations so that 
these do not have to be constructed/destructed by every method that 
uses them.  

@section collaborations Collaborations
 
Forms a bridge pattern with the Element class which provides an
interface to it. Used by any FEM calculations facilitated by the 
Algorithm class. 

@todo make this class a Singleton and avoid pointer connections with elements

*/
class FiniteElement {
  public:
    FiniteElement( CSMP_FEM_TYPE csp_fem_type, 
                   bool isoparametric, bool uses_local_coordinates, 
                   size_t order_of_shape_functions );

    virtual ~FiniteElement() {}

    /// assigns integer value used to avoid repeating the same operation
    void               CurrentID( size_t id );
  
    /// returns integer value for comparison to avoid repeating the same operation
    size_t             CurrentID() const;
  
    /// returns ID value set to FEM during construction
    static size_t      InitialID();
  
    /// return minimum spatial dimension in which this element can exist
    size_t             Dim() const;
  
    size_t             Nodes() const;
  
    /// returns the number of edges of the FE, i.e. the number of connections between nodes
    size_t             Segments() const;
  
    /// returns the number of sides the finite element has
    size_t             Faces() const;
  
    /// neighbor elements that may be connected to this FE in its characteristic space
    size_t             Neighbors() const;
  
    /// returns how many nodes make up a particular face
    size_t             NodesPerFace( size_t ) const;
  
    /// @todo (1) SKM: WHAT IS THIS
    size_t             IntegrationPointNeighbors() const;
  
    /// returns number of quadrature points used by current integration scheme
    size_t             IntegrationPoints() const;
  
    /// returns order of interpolation scheme: linear=1, quadratic=2, cubic=3
    size_t             Interpolation() const;
  
    /// for numerically integrated elements, reports whether the order of the interpolation functions is the same as that of the shape functions
    bool               Isoparametric() const;
  
    /// reports whether element is defined in a local (r,s,t) coordinate system; analytically integrated elements are not
    bool               UsesLocalCoordinates() const;
  
    /// returns the interpolation order of the shape functions as opposed to the interpolation (basis) functions
    size_t             OrderOfShapeFunctions() const;
  
    bool               IsLineElement() const;
    bool               IsSurfaceElement() const;
    bool               IsVolumeElement() const;
  
    CSMP_FEM_TYPE      ElementType() const;
  
    /// computes and returns the volume of the element which is an area for a surface- and a length for a line element
    virtual double64   Volume();
  
    /// returns the ratio between the maximum and minimum spatial extent of the element
    virtual double64   AspectRatio();
  
    /// returns the radius of the largest sphere that could be inscribed in the element
    virtual double64   InnerRadius();
  
    /// reports the lengths of all edges (connections between nodes) in the currrent element
    virtual void       EdgeLengths( std::vector<double64>& vec );
  
    /// returns 0..nodes-1 (local) node numbers of the supplied boundary nodes, i.e., reports the actual ordering
    virtual void       ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes, 
                                                   std::vector<size_t>& fnids );
                                                   
    // all of the following give 0...nodes-1 (local) node numbers
    /// reports the egde nodes (0..nodes-1) which can be more than 2 in higher-order elements
    virtual void       NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
  
    /// reports the face nodes (0..nodes-1) of the given face
    virtual void       NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
  
    /// reports only the corner nodes (0..nodes-1) of the element
    virtual void       CornerNodes( std::vector<size_t>& ids ) const;

    /// reports only the mid-edge nodes (0..nodes-1) of the element
    virtual void       MidSideNodes( std::vector<size_t>& ids ) const;

    /// reports the nodes (0..nodes-1) of the element which are not shared with any other element
    virtual void       InteriorNodes( std::vector<size_t>& ids ) const;
  
    virtual size_t     CornerNodes() const;
    virtual size_t     MidSideNodes() const;
    virtual size_t     InteriorNodes() const;
  
    /// reports the element type that would match that of its corresponding face
    virtual CSMP_FEM_TYPE ElementTypeOfFace( size_t face ) const;
  
    /// reports the element type that would match that of its corresponding segment
    virtual CSMP_FEM_TYPE ElementTypeOfSegment( size_t segment ) const;

    /// returns the unit normal where it makes sense (line and surface elements); expects XY matrix to be initialised with node coordinates
    virtual void       UnitNormal( std::vector<double64>& unrml ) const;

    // TODO: SKM: add more efficient unit-normal calculation that returns point 
    /// returns the unit normal to element face (for planar elements, normal lies in that plane; expects XY matrix to be initialised with node coordinates
    virtual void       UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;
  
    /// computes unit normal to face using a parametric to physical space transformation available in elements with a local coordinate framework
    virtual void       UnitNormalAtFaceBarycenter( size_t face, std::vector<double64>& nrml );
  
    /// accessor for coordinate matrix which is associated with the element
    double64           XYZ( size_t i, size_t j ) const;

    /// mutator for coordinate matrix which is associated with the element
    void               XYZ( size_t i, size_t j, double64 val );

    /// for element types defined using a local coordinate system, returns the node coordinates into the argument matrix
    virtual void       ReferenceCoordinates( DenseMatrix<DM_MIN>& matCoords ) const;

    /// using the interpolation functions, the supplied suite of scalar variable values is extrapolated from the Gauss points to the nodes
    virtual   void      ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double64>& IVAR,
                                                                    std::vector<double64>&       NVAR ) const;

    /// returns the global coordinates of the element integration point
    virtual   void      IntegrationPoint( size_t i, std::vector<double64>& xyz ) const;
  
    /// returns the integration weight of the desired quadrature point
    virtual   double64  WeightAtIntegrationPoint( size_t i ) const;
  
    /// reports the local node numbers in counter clockwise order
    virtual   void      CounterClockwiseNodes( std::vector<size_t>& ids ) const;

    /// reports the values of the interpolation functions at the point given in global coordinates; @attention slow for numerically integrated elements
    virtual   void      N( std::vector<double64>& N, const std::vector<double64>& xyz );
  
    /// reports the values of the interpolation functions at the given quadrature point
    virtual   void      N_AtIntegrationPoint( size_t ip, std::vector<double64>& N );
  
    /// reports the values of the interpolation functions at the center of gravity of the element
    virtual   void      N_AtBaryCenter( std::vector<double64>& N );
    
    /// first derivatives of interpolation functions for elements where this is a single constant value
    virtual   void      dN( DenseMatrix<DM_MIN>& );
  
    /// returns first derivative of interpolaton functions at global point; in isoparametric elements, the determinant of the Jacobian is returned as well
    virtual   double64  dN_At( DenseMatrix<DM_MIN>&, const std::vector<double64>& xyz );

    /// returns first derivative of interpolaton functions at quadrature point; determinant of the Jacobian is returned as well
    virtual   double64  dN_AtIntegrationPoint( DenseMatrix<DM_MIN>&, size_t gauss_point );

    /// returns first derivative of interpolaton functions at node; determinant of the Jacobian is returned as well
    virtual   double64  dN_AtNode( DenseMatrix<DM_MIN>&, size_t node );

    /// returns first derivative of interpolaton functions at center of gravity; determinant of the Jacobian is returned as well
    virtual   double64  dN_AtBarycenter( DenseMatrix<DM_MIN>& );
    
    // for isoparametric elements Jacobian coordinate transformations (at a point within the element)
    // are returned into the protected matrices JAC and JINV
  
    /// initialises JAC = DN*XY (intpl. derivative * coordinate matrix) parametric-to-physical space transformation matrix for given quadrature point
    virtual   void      JacobianAtIntegrationPoint( size_t ip );

    /// initialises JAC = DN*XY (intpl. derivative * coordinate matrix) parametric-to-physical space transformation matrix for given quadrature point defined in parametric space (r,s,t)
    virtual   void      JacobianAt( const std::vector<double64>& rst );
  
    /// initialises JAC = DN*XY (intpl. derivative * coordinate matrix) parametric-to-physical space transformation matrix for 1D element using the supplied intpol.f.derivatives
    virtual   void      Jacobian( const std::vector<double64>& dnr ); // 1D

    /// initialises JAC = DN*XY (intpl. derivative * coordinate matrix) parametric-to-physical space transformation matrix for 2D element using the supplied intpol.f.derivatives
    virtual   void      Jacobian( const std::vector<double64>& dnr, const std::vector<double64>& dns ); // 2D

    /// initialises JAC = DN*XY (intpl. derivative * coordinate matrix) parametric-to-physical space transformation matrix for 3D element using the supplied intpol.f.derivatives
    virtual   void      Jacobian( const std::vector<double64>& dnr, const std::vector<double64>& dns,
                                  const std::vector<double64>& dnt ); // 3D
  
    /// initializes the member matrix JINV
    virtual   double64  JacobianInverse();
  
    /// returns the determinant of the member Jacobian matrix JAC that must have been initialised before
    virtual   double64  JacobianDeterminant();
	
    /// initialises interpolation function product matrix N^T x N for analytically integrated elements
    virtual   void      IntegralNN( DenseMatrix<DM_MIN>& );
  
    /// initialises interpolation function derivatives product matrix DN^T x DN for element where derivatives are constant
    virtual   void      IntegraldNdN( DenseMatrix<DM_MIN>& );
    
    // local interpolation functions in elements that use a local coordinate system (r,s,t)

    /// 1D element interpolation functions N(r)
    virtual void  Nr(  double64 r, std::vector<double64>& NRST ) const;
    /// 1D element interpolation functions derivatives N(r)/dr
    virtual void  dNr( double64 r, std::vector<double64>& DNR ) const;

    /// 2D element interpolation functions N(r,s)
    virtual void  Nrs( double64 r, double64 s, std::vector<double64>& NRST ) const;
    /// 2D element interpolation functions derivatives N(r,s)/dr
    virtual void  dNr( double64 r, double64 s, std::vector<double64>& DNR ) const;
    /// 2D element interpolation functions derivatives N(r,s)/ds
    virtual void  dNs( double64 r, double64 s, std::vector<double64>& DNS ) const;

    /// 3D element interpolation functions N(r,s,t)
    virtual void  Nrst( double64 r, double64 s, double64 t, std::vector<double64>& NRST ) const;
    /// 2D element interpolation functions derivatives N(r,s,t)/dr
    virtual void  dNr( double64 r,  double64 s, double64 t, std::vector<double64>& DNR ) const;
    /// 2D element interpolation functions derivatives N(r,s,t)/ds
    virtual void  dNs( double64 r,  double64 s, double64 t, std::vector<double64>& DNS ) const;
    /// 2D element interpolation functions derivatives N(r,s,t)/dt
    virtual void  dNt( double64 r,  double64 s, double64 t, std::vector<double64>& DNT ) const;
    
    /// prints finite element properties to screen
    virtual void  Out(std::ostream& os) const;
  
    /// prints the current finite element into a visualisation toolkit (VTK) textfile for visual examination
    virtual void  OutputNodeDataToVTK( const char* file_name, const char* var_name, 
                                       DenseMatrix<DM_MIN>& DATA ) const;

  protected:
    FiniteElement( const FiniteElement& e );
    FiniteElement& operator=( const FiniteElement& e );
    void           Isoparametric( bool isoparam );
    void           UsesLocalCoordinates( bool uses );
    void           LineElement();
    void           SurfaceElement();
    void           VolumeElement();
    void           ElementType( CSMP_FEM_TYPE etype );

    size_t dim,       /**< spatial dimension of element */
           itp,       /**< degree of interpolation */
           npf,       /**< nodes per face */
           npe,       /**< nodes per element  */
           spe,       /**< segments per element  */
           fpe,       /**< faces per element  */
           epe,       /**< neighbors of element */
           nne,       /**< typical number of elements that share each node */
           cne,       /**< typical number of elements that share each integration point */
           gpe;       /**< Gauss points per element for numerical integration */

    DenseMatrix<DM_MIN>   M;         /**< test-function coefficients */
                     
  public:
    DenseMatrix<DM_MIN>   XY;        /**< node coordinates */
    DenseMatrix<DM_MIN>   JAC;       /**< Jacobian matrix */
    DenseMatrix<DM_MIN>   JINV;      /**< the Jacobians inverse */

    mutable std::vector<double64>  NRST, ///< convenience storage for interpolation function values
                                   DNR,  ///< convenience storage for interpolation function derivatives
                                   DNS,  ///< convenience storage for interpolation function derivatives
                                   DNT;  ///< convenience storage for interpolation function derivatives
    mutable std::vector<size_t>    IDX;  ///< convenience storage for integer vectors

  private:
    FiniteElement();
    bool    isoparametric,               ///< true if interpolation order of basis- and shape functions is the same
            uses_local_coordinates;      ///< for elements defined in a local (r,s,t) coordinate system
  
    size_t  order_of_shape_functions,    ///< order of element shape functions
            object_id;                   ///< number used to track operations done on a particular element
    
    enum { LINE, SURFACE, VOLUME } element_category; ///< element classifier
    enum CSMP_FEM_TYPE             csp_fem_type;     ///< elements supported by CSMP, see definition above
  
    /// stub for all virtual functions, providing feedback to users if member functions are not defined for particular element type
    void InstructUser( const char* method ) const;
};


// inlined methods

inline bool isTriangularElement( CSMP_FEM_TYPE etype )
 {
    // most likely case first for performance
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return true;
    if ( etype == LINEAR_TRIANGLE ) return true;
    if ( etype == LINEAR_TRIANGLE3D ) return true;
    if ( etype == QUADRATIC_TRIANGLE ) return true;
    if ( etype == BARYCENTRIC_LINEAR_TRIANGLE ) return true;
    if ( etype == BARYCENTRIC_QUADRATIC_TRIANGLE ) return true;
    if ( etype == CUBIC_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE ) return true;
    return false;
 }
 
inline bool isQuadrilateralElement( CSMP_FEM_TYPE etype )
 {
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return true;
    return false;
 }
 
inline bool isLineElement( CSMP_FEM_TYPE etype )
 {
    if ( etype == ISOPARAMETRIC_LINEAR_BAR ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_BAR ) return true;
    if ( etype == LINEAR_BAR ) return true;
    if ( etype == QUADRATIC_BAR ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_BAR ) return true;
    if ( etype == CUBIC_BAR ) return true;
    return false;
 }


inline void    FiniteElement::CurrentID( size_t id ) { object_id = id; }
inline size_t  FiniteElement::CurrentID() const      { return object_id; }

/// initalizing to a value that makes sure that ID does not equal initial element idx
inline size_t  FiniteElement::InitialID() { return UINT_MAX; }


inline bool  FiniteElement::Isoparametric() const        { return isoparametric; }

inline bool  FiniteElement::UsesLocalCoordinates() const { return uses_local_coordinates; }

inline bool  FiniteElement::IsLineElement()    const     
 { if ( element_category == LINE )   return true; return false; }
 
inline bool  FiniteElement::IsSurfaceElement() const     
 { if ( element_category == SURFACE ) return true; return false; }
 
inline bool  FiniteElement::IsVolumeElement()  const     
 { if ( element_category == VOLUME )  return true; return false; }

inline void  FiniteElement::ElementType( CSMP_FEM_TYPE etype ) { csp_fem_type = etype; }

inline CSMP_FEM_TYPE  FiniteElement::ElementType() const       { return csp_fem_type; }

inline  size_t  FiniteElement::OrderOfShapeFunctions() const { return order_of_shape_functions; }

inline  size_t  FiniteElement::Interpolation()    const { return itp; }
inline  size_t  FiniteElement::Dim()              const { return dim; }
inline  size_t  FiniteElement::Nodes()            const { return npe; }
inline  size_t  FiniteElement::Segments()         const { return spe; }
inline  size_t  FiniteElement::Faces()            const { return fpe; }
inline  size_t  FiniteElement::Neighbors()        const { return epe; }
inline  size_t  FiniteElement::NodesPerFace( size_t ) const { return npf; }
inline  size_t  FiniteElement::IntegrationPointNeighbors() const { return cne; }
inline  size_t  FiniteElement::IntegrationPoints() const { return gpe; }

inline  double64  FiniteElement::XYZ( size_t i, size_t j ) const { return XY(i,j); }
inline  void    FiniteElement::XYZ( size_t i, size_t j, double64 val ) { XY(i,j) = val; }

} // csmp

#endif
