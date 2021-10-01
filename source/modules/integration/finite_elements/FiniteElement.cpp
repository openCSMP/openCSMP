#include "FiniteElement.h"
#include "CSMP_definitions.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

/**

This default constructor set the element specifications to a linear
triangular element which uses global coordinates and has linear 
shape functions.  

@section application Application 

This constructor is used to construct the FiniteElement baseclass, thus
it is lazy in that it leaves it to the FiniteElement subclass to initialize
many of its protected variable. This is done, since the FiniteElement is
quasi an abstract class. Only the subclass knows the critical values.  
*/
FiniteElement::FiniteElement() 
 : XY(3,2), M(3,3), JAC(2,2), JINV(2,2), 
   isoparametric(false), 
   uses_local_coordinates(false),
   order_of_shape_functions(1),
   element_category(SURFACE),
   csp_fem_type(UNKNOWN),
   object_id( InitialID() )  
 {
 }


FiniteElement::FiniteElement( CSMP_FEM_TYPE csp_fem_type,
                              bool isoparametric, 
                              bool uses_local_coordinates,  
                              size_t order_of_shape_functions )
 : XY(3,2), M(3,3), JAC(2,2), JINV(2,2), 
   isoparametric(isoparametric), 
   uses_local_coordinates(uses_local_coordinates),
   order_of_shape_functions(order_of_shape_functions),
   element_category(SURFACE),
   csp_fem_type(csp_fem_type),
   object_id( InitialID() )
 {
 }


FiniteElement::FiniteElement( const FiniteElement& e ) 
  { 
     *this = e; 
  }



FiniteElement&  FiniteElement::operator=( const FiniteElement& e )
 {
    if ( this != &e ) 
      {
        dim = e.dim;
        itp = e.itp;
        npe = e.npe;
        spe = e.spe;
        fpe = e.fpe;
        epe = e.epe;
        nne = e.nne;
        cne = e.cne;
        gpe = e.gpe;
        isoparametric            = e.isoparametric;
        uses_local_coordinates   = e.uses_local_coordinates;
        order_of_shape_functions = e.order_of_shape_functions;
        element_category         = e.element_category;
        csp_fem_type             = e.csp_fem_type;
        object_id                = e.object_id;
        XY  = e.XY;
        M   = e.M;
      }
    return *this;
 }

void  FiniteElement::Isoparametric( bool isoparam ) { isoparametric = isoparam; }
void  FiniteElement::UsesLocalCoordinates( bool uses ) { uses_local_coordinates = uses; }

// inlined methods

bool isTriangularElement( CSMP_FEM_TYPE etype )
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
 
bool isQuadrilateralElement( CSMP_FEM_TYPE etype )
 {
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return true;
	  if ( etype == LINEAR_RECTANGLE ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return true;
    return false;
 }
 
bool isLineElement( CSMP_FEM_TYPE etype )
 {
    if ( etype == ISOPARAMETRIC_LINEAR_BAR ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_BAR ) return true;
    if ( etype == LINEAR_BAR ) return true;
    if ( etype == QUADRATIC_BAR ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_BAR ) return true;
    if ( etype == CUBIC_BAR ) return true;
    return false;
 }



bool isTriangular( CSMP_FEM_TYPE etype )
 {
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ||
         etype == LINEAR_TRIANGLE || 
         etype == LINEAR_TRIANGLE3D || 
         etype == BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == QUADRATIC_TRIANGLE ||
         etype == BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == CUBIC_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_CUBIC_TRIANGLE  )
      return true;

    return false;
 }


bool isQuadrilateral( CSMP_FEM_TYPE etype )
 {
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL || 
         etype == LINEAR_QUADRILATERAL || 
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ||
         etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL )
      return true;

    return false;
 }


bool isHexahedral( CSMP_FEM_TYPE etype )
 {
    if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON || 
         etype == LINEAR_QUADRILATERAL || 
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ||
         etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ||
         etype == LINEAR_CUBOID )
      return true;

    return false;
 } 


void    FiniteElement::CurrentID( size_t id ) { object_id = id; }
size_t  FiniteElement::CurrentID() const      { return object_id; }

/// initalizing to a value that makes sure that ID does not equal initial element idx
size_t  FiniteElement::InitialID() { return UINT_MAX; }


bool  FiniteElement::Isoparametric() const        { return isoparametric; }

bool  FiniteElement::UsesLocalCoordinates() const { return uses_local_coordinates; }

bool  FiniteElement::IsLineElement()    const     
 { if ( element_category == LINE )   return true; return false; }
 
bool  FiniteElement::IsSurfaceElement() const     
 { if ( element_category == SURFACE ) return true; return false; }
 
bool  FiniteElement::IsVolumeElement()  const     
 { if ( element_category == VOLUME )  return true; return false; }
 
bool FiniteElement::IsSimplex() const
 {
    if ( csp_fem_type == ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return true;
    if ( csp_fem_type == ISOPARAMETRIC_LINEAR_TRIANGLE ) return true;
    if ( csp_fem_type == ISOPARAMETRIC_LINEAR_BAR ) return true;
    if ( csp_fem_type == LINEAR_TETRAHEDRON ) return true;
    if ( csp_fem_type == LINEAR_TRIANGLE3D ) return true;
    if ( csp_fem_type == LINEAR_TRIANGLE ) return true;
    if ( csp_fem_type == LINEAR_BAR ) return true;
    return false;
 }

void  FiniteElement::ElementType( CSMP_FEM_TYPE etype ) { csp_fem_type = etype; }

CSMP_FEM_TYPE  FiniteElement::ElementType() const       { return csp_fem_type; }

 size_t  FiniteElement::OrderOfShapeFunctions() const { return order_of_shape_functions; }

 size_t  FiniteElement::Interpolation()    const { return itp; }
 size_t  FiniteElement::Dim()              const { return dim; }
 size_t  FiniteElement::Nodes()            const { return npe; }
 size_t  FiniteElement::Segments()         const { return spe; }
 size_t  FiniteElement::Faces()            const { return fpe; }
 size_t  FiniteElement::Neighbors()        const { return epe; }
 size_t  FiniteElement::NodesPerFace( size_t ) const { return npf; }
 size_t  FiniteElement::IntegrationPointNeighbors() const { return cne; }
 size_t  FiniteElement::IntegrationPoints() const { return gpe; }

 double64  FiniteElement::XYZ( size_t i, size_t j ) const { return XY(i,j); }
 void    FiniteElement::XYZ( size_t i, size_t j, double64 val ) { XY(i,j) = val; }



void FiniteElement::InstructUser( const char* method ) const
 {
     cout <<"\nFiniteElement::InstructUser: called by: "<< method << endl;
     cout <<"Info: The method you have called is a ";
     cout <<"virtual function of the base class FiniteElement. It does nothing ";
     cout <<"except for printing this message."<< endl;
     cout <<"The reason why you get this message is that the specific desired method: '";
     cout << method <<"' is not defined for the element type you are using.";
     cout <<"To implement the desired method, go to the FiniteElement subclass ";
     cout <<" corresponding to your element and define an implementation of the method."<< endl;
     cout <<"Values of supplied arguments follow: "<< endl << endl; 
 }


void FiniteElement::LineElement()     { element_category = LINE; }
void FiniteElement::SurfaceElement()  { element_category = SURFACE; }
void FiniteElement::VolumeElement()   { element_category = VOLUME; }


double64 FiniteElement::AspectRatio()
  {
     InstructUser("FiniteElement::AspectRatio()");
     cout <<"\ncalled by object: "<< object_id << endl;
     return 0.0;
  }
  
double64 FiniteElement::InnerRadius()
  {
     InstructUser("FiniteElement::InnerRadius()");
     cout <<"\ncalled by object: "<< object_id << endl;
     return 0.0;
  }
  
void FiniteElement::EdgeLengths( vector<double64>& vec )
  {
     InstructUser("FiniteElement::EdgeLengths(vector<double64>)");
     out(vec);
     cout <<"\ncalled by object: "<< object_id << endl;
  }



void FiniteElement::NodesOfSegment( size_t sid, vector<size_t>& snids ) const
  {
     cout <<"\nFiniteElement::NodesOfSegment: Returns the local node ID numbers of ";
     cout <<"the nodes which constitute the element segment with the entered (local) ID number. "<< endl;
     cout <<"\ncalled by object: "<< object_id <<" for face "<< sid << endl;
     out(snids);
     throw invalid_argument("FiniteElement::NodesOfSegment");
  }


void FiniteElement::NodesOfFace( size_t fid, vector<size_t>& fnids ) const
  {
     cout <<"\nFiniteElement::NodesOfFace: Returns the local node ID numbers of ";
     cout <<"the nodes which constitute the element face with the entered ID number. ";
     cout <<"In triangular and tetrahedral elements the faces lie opposite of ";
     cout <<"the nodes with the same ID." << endl;
     cout <<"\ncalled by object: "<< object_id <<" for face "<< fid << endl;
     out(fnids);
     throw invalid_argument("FiniteElement::NodesOfFace");
  }


  
    
void FiniteElement::IntegraldNdN( DenseMatrix<DM_MIN>& M )
  {
     InstructUser("FiniteElement::IntegraldNdN( DenseMatrix<DM_MIN>& M )");
     M.Out();
     cout <<"\ncalled by object: "<< object_id << endl;
     throw invalid_argument("FiniteElement::IntegraldNdN");
  }
  

  
// other information

void FiniteElement::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                vector<size_t>& fnids )
  {
     InstructUser("FiniteElement::ConsecutiveNodesAtBoundary");
     cout <<"\nnodes at boundary: "<< endl;
     out(bnodes);
     out(fnids);
     cout <<"\ncalled by object: "<< object_id << endl;
     throw invalid_argument("FiniteElement::ConsecutiveNodesAtBoundary");
  }


void FiniteElement::CornerNodes( vector<size_t>& ids ) const
 {
    InstructUser("FiniteElement::CornerNodes");
    out(ids);
    throw invalid_argument("FiniteElement::CornerNodes");
 }
 
 
void FiniteElement::MidSideNodes( vector<size_t>& ids ) const
 {
    InstructUser("FiniteElement::MidSideNodes");
    out(ids);
    throw invalid_argument("FiniteElement::MidSideNodes");
 }

void FiniteElement::InteriorNodes( vector<size_t>& ids ) const
 {
    InstructUser("FiniteElement::InteriorNodes");
    out(ids);
    throw invalid_argument("FiniteElement::InteriorNodes");
 }

size_t FiniteElement::CornerNodes() const
 {
    InstructUser("FiniteElement::CornerNodes");
    throw invalid_argument("FiniteElement::CornerNodes");
    return 0;
 }
 
 
size_t FiniteElement::MidSideNodes() const
 {
    InstructUser("FiniteElement::MidSideNodes");
    throw invalid_argument("FiniteElement::MidSideNodes");
    return 0;
 }

size_t FiniteElement::InteriorNodes() const
 {
    InstructUser("FiniteElement::InteriorNodes");
    throw invalid_argument("FiniteElement::InteriorNodes");
    return 0;
 }


CSMP_FEM_TYPE FiniteElement::ElementTypeOfFace( size_t ) const
  {
     /// @todo (2-P) Remove typid by name fct
     InstructUser("FiniteElement::ElementTypeOfFace");
     cout <<"\nDerived Element: "<< typeid(this).name() << endl;
     cout <<"\nThis method must be defined for the FE element implementation which you are using"<< endl;
     throw invalid_argument("FiniteElement::ElementTypeOfFace");
     
     return UNKNOWN;
  }


CSMP_FEM_TYPE FiniteElement::ElementTypeOfSegment( size_t ) const
  {
     /// @todo (2-P) Remove typid by name fct
     InstructUser("FiniteElement::ElementTypeOfSegment");
     cout <<"\nDerived Element: "<< typeid(this).name() << endl;
     cout <<"\nThis method must be defined for the FE element implementation which you are using"<< endl;
     throw invalid_argument("FiniteElement::ElementTypeOfSegment");
     
     return UNKNOWN;
  }



void  FiniteElement::UnitNormal( vector<double64>& ) const
  {
     InstructUser("FiniteElement::UnitNormal");
     cout <<"\nThis method is not defined for the FE element implementation which you are using"<< endl;
     cout <<"\nPerhaps, the element type is not planar such that a normal is not an applicable concept"<< endl;
     throw invalid_argument("FiniteElement::UnitNormal(vector<double>)");
  }


/**
    Virtual function stub, overload method in desired subclass of FiniteElement
    
    method should return the unit normal to element face (for planar elements, this normal will lie in that plane
    
    method expects that the XY matrix is initialised with node coordinates
*/
void  FiniteElement::UnitNormalToFace( size_t, std::vector<double64>& ) const
 {
     InstructUser("FiniteElement::UnitNormalToFace");
     cout <<"\nThis method is not defined for the FE element implementation which you are using"<< endl;
     cout <<"\nyou need to overload (=define) the method in the desired FiniteElement subclass"<< endl;
     throw invalid_argument("FiniteElement::UnitNormalToFace(size_t,vector<double>)");
 }



/**
     This method is intended for numerically integrated element types with higher order 
     shape functions. As these elements can get deformed in physical space, it is important 
     to specify the precise location of the unit normal.
     
     @note this method cannot be constant because it modifies matrices that are data members
     of the FiniteElement class.
*/
void  FiniteElement::UnitNormalAtFaceBarycenter( size_t, std::vector<double64>& )
 {
     InstructUser("FiniteElement::UnitNormalAtFaceBarycenter");
     cout <<"\nThis method is not defined for the FE element implementation which you are using"<< endl;
     cout <<"\nyou need to overload (=define) the method in the desired FiniteElement subclass"<< endl;
     throw invalid_argument("FiniteElement::UnitNormalAtFaceBarycenter(size_t,vector<double>)");
 }



double64 FiniteElement::WeightAtIntegrationPoint( size_t ) const
  {
     throw invalid_argument("FiniteElement::WeightAtIntegrationPoint");
     return 0.5 * (1. / gpe); 
  }



void  FiniteElement::IntegrationPoint( size_t i, std::vector<double64>& xyz ) const
 { 
    cout <<"\nFiniteElement::IntegrationPoint: called for integration point "<< i;
    cout <<" and output data vector 'xyz' of size "<< xyz.size() << endl;
    cout <<"\nThis method which returns the position of the integration point ";
    cout <<" in global coordinates must be overloaeded for the element of interest."<< endl;
    throw invalid_argument("FiniteElement::IntegrationPoint");
 }




void  FiniteElement::CounterClockwiseNodes( vector<size_t>& ids ) const
 {
     InstructUser("FiniteElement::CounterClockwiseNodes");
     out(ids);
     cout <<"\nThis method must be defined by the FE element which you are using"<< endl;
     throw invalid_argument("FiniteElement::CounterClockwiseNodes");
 }


/** Extrapolation member( gauss points to nodes )

Expects input vector to be of size nvars * gauss points and in format
[ var1 @ gp1, var2 @ gp1, var3 @ gp1,
  var1 @ gp2, var2 @ gp2, var3 @ gp2 for all gauss points, with nvar == 3 ]

  */
void  FiniteElement::ExtrapolateIntegrationPointVariableToNodes( size_t nvars, 
                                                                 const vector<double64>&, 
                                                                 vector<double64>& ) const
 {
     InstructUser("FiniteElement::ExtrapolateIntegrationPointVariableToNodes");
     cout <<"\nThis method must be defined for the FE element which you are using"<< endl;
     cout <<"\nSupplied data for N-variables = "<< nvars << endl;
     throw invalid_argument("FiniteElement::ExtrapolateIntegrationPointVariableToNodes");
 } 




double64 FiniteElement::Volume()
 {
    InstructUser("FiniteElement::Volume");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    throw invalid_argument("FiniteElement::Volume");
    return 0.0;
 } 


void   FiniteElement::N( vector<double64>& N, const vector<double64>& xyz )
 {
    InstructUser("FiniteElement::N");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (N vector<double64>, xyz vector<double64>):"<< endl;
    out( N );
    for ( size_t i=0; i<xyz.size(); i++ ) cout << xyz[i] <<" ";
    cout << endl;
    throw invalid_argument("FiniteElement::N");
 } 


void   FiniteElement::N_AtIntegrationPoint( size_t ip, vector<double64>& N )
 {
    InstructUser("FiniteElement::N_AtIntegrationPoint");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (ip, N vector<double64>):"<< ip << endl;
    out( N );
    throw invalid_argument("FiniteElement::N_AtIntegrationPoint");
 } 


void   FiniteElement::N_AtBaryCenter( vector<double64>& N )
 {
    InstructUser("FiniteElement::N_AtBaryCenter");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (ip, N vector<double64>):"<< endl;
    out( N );
    throw invalid_argument("FiniteElement::N_AtBaryCenter");
 }

void FiniteElement::Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K)
{
	InstructUser("FiniteElement::Integral_dNT_K_dN");
	cout << "\nThis method is not defined for the FE element type which you are using" << endl;
	cout << "\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M):" << endl;
	M.Out();
	throw invalid_argument("FiniteElement::dN");
}


void   FiniteElement::dN( DenseMatrix<DM_MIN>& M ) // coefficients 
 {
    InstructUser("FiniteElement::dN");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M):"<< endl;
    M.Out();
    throw invalid_argument("FiniteElement::dN");
 } 
    
    
double64 FiniteElement::dN_At( DenseMatrix<DM_MIN>& M, const vector<double64>& xyz  ) // derivatives at 'xy'
 {
    InstructUser("FiniteElement::dN_At");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M, vector<double64> xyz):"<< endl;
    M.Out();
    for ( size_t i=0; i<xyz.size(); i++ ) cout << xyz[i] <<" ";
    cout << endl;
    throw invalid_argument("FiniteElement::dN_At");
    return 0.0;
 } 


double64 FiniteElement::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point ) 
 {
    InstructUser("FiniteElement::dN_AtIntegrationPoint");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (gp, dof, DenseMatrix<DM_MIN> M):"<< endl;
    cout << gauss_point <<", "<< endl;
    M.Out();
    throw invalid_argument("FiniteElement::dN_AtIntegrationPoint");
    return 0.0;
 } 


double64 FiniteElement::dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node ) 
 {
    InstructUser("FiniteElement::dN_AtNode");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (node, dof, DenseMatrix<DM_MIN> M): ";
    cout << node <<", "<< endl;
    M.Out();
    throw invalid_argument("FiniteElement::dN_AtNode");
    return 0.0;
 } 


double64 FiniteElement::dN_AtBarycenter( DenseMatrix<DM_MIN>& M ) 
 {
    InstructUser("FiniteElement::dN_AtBarycenter");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M): " << endl;
    M.Out();
    throw invalid_argument("FiniteElement::dN_AtBarycenter");
    return 0.0;
 } 
    

void   FiniteElement::IntegralNN( DenseMatrix<DM_MIN>& M )
 {
    InstructUser("FiniteElement::IntegralNN");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M):"<< endl;
    M.Out();
    throw invalid_argument("FiniteElement::IntegralNN");
 } 



void   FiniteElement::OutputNodeDataToVTK( const char* file_name, const char* var_name, 
                                           DenseMatrix<DM_MIN>& DATA ) const
 {
    InstructUser("FiniteElement::OutputNodeDataToVTK");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (filename, varname, DATA): "<< endl;
    cout <<"'"<< file_name <<"', '"<< var_name <<"'"<< endl;
    DATA.Out();
    throw invalid_argument("FiniteElement::OutputNodeDataToVTK");
 } 



void   FiniteElement::JacobianAtIntegrationPoint( size_t ip )
 {
    InstructUser("FiniteElement::JacobianAtIntegrationPoint");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"Just write a new method which calculates the Jacobian matrix for the element in question.";
    cout <<"\nMethod arguments for single-element mesh output, input (integration point): "<< endl;
    cout <<"'"<< ip <<"', '"<< endl;
    throw invalid_argument("FiniteElement::JacobianAtIntegrationPoint");
 }




    // 1D
void   FiniteElement::Nr( double64 r, vector<double64>& NR ) const
 {
    InstructUser("FiniteElement::Nr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinate r: "<< r << endl;
    out( NR );
    cout << endl;
    throw invalid_argument("FiniteElement::Nr");
 }



    // 1D
void   FiniteElement::Nr( double64 r, double64* NR ) const
 {
    InstructUser("FiniteElement::Nr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinate r: "<< r << endl;
    cout << endl;
    throw invalid_argument("FiniteElement::Nr");
 }


void   FiniteElement::dNr( double64 r, vector<double64>& DNR ) const
 {
    InstructUser("FiniteElement::dNr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinate r: "<< r << endl;
    out( DNR );
    cout << endl;
    throw invalid_argument("FiniteElement::dNr");
 }


    // 2D
void   FiniteElement::Nrs( double64 r, double64 s, vector<double64>& NRS ) const
 {
    InstructUser("FiniteElement::Nrs");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s: "<< r <<" "<< s << endl;
    out( NRS );
    cout << endl;
    throw invalid_argument("FiniteElement::Nrs");
 }

    // 2D
void   FiniteElement::Nrs( double64 r, double64 s, double64* NRS ) const
 {
    InstructUser("FiniteElement::Nrs");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s: "<< r <<" "<< s << endl;
    cout << endl;
    throw invalid_argument("FiniteElement::Nrs");
 }


void   FiniteElement::dNr( double64 r, double64 s, vector<double64>& DNR ) const
 {
    InstructUser("FiniteElement::dNr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s: "<< r <<" "<< s << endl;
    out( DNR );
    cout << endl;
    throw invalid_argument("FiniteElement::dNr");
 }


void   FiniteElement::dNs( double64 r, double64 s, vector<double64>& DNS ) const
 {
    InstructUser("FiniteElement::dNs");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s: "<< r <<" "<< s << endl;
    out( DNS );
    cout << endl;
    throw invalid_argument("FiniteElement::dNs");
 }


    // 3D
void   FiniteElement::Nrst( double64 r, double64 s, double64 t, vector<double64>& NRST ) const
 {
    InstructUser("FiniteElement::Nrst");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    out( NRST );
    cout << endl;
    throw invalid_argument("FiniteElement::Nrst");
 }


    // 3D
void   FiniteElement::Nrst( double64 r, double64 s, double64 t, double64* NRST ) const
 {
    InstructUser("FiniteElement::Nrst");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    cout << endl;
    throw invalid_argument("FiniteElement::Nrst");
 }


void   FiniteElement::dNr( double64 r,  double64 s, double64 t, vector<double64>& DNR ) const
 {
    InstructUser("FiniteElement::dNr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    out( DNR );
    cout << endl;
    throw invalid_argument("FiniteElement::dNr");
 }


void   FiniteElement::dNs( double64 r,  double64 s, double64 t, vector<double64>& DNS ) const
 {
    InstructUser("FiniteElement::dNs");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    out( DNS );
    cout << endl;
    throw invalid_argument("FiniteElement::dNs");
 }


void   FiniteElement::dNt( double64 r,  double64 s, double64 t, vector<double64>& DNT ) const
 {
    InstructUser("FiniteElement::dNt");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    out( DNT );
    cout << endl;
    throw invalid_argument("FiniteElement::dNt");
 }


void FiniteElement::JacobianAt( const std::vector<double64>& rst )
{
    InstructUser("FiniteElement::JacobianAt");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< rst[0] <<" "<< rst[1] <<" "<< rst[2] << endl;
    out( DNT );
    cout << endl;
    throw invalid_argument("FiniteElement::JacobianAt");
}


/**

Returns the reference coordinates of the element. 

DenseMatrix<DM_MIN> & matCoords; -> Holds the reference coordinates of the finite element. 
*/
void FiniteElement::ReferenceCoordinates(DenseMatrix<DM_MIN> & /*matCoords*/) const
{
    InstructUser("FiniteElement::ReferenceCoordinates");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    out( DNT );
    cout << endl;
    throw invalid_argument("FiniteElement::ReferenceCoordinates");
}
	

/**

Compute the Jacobian matrix for the local coordinate point 'rst' in the 
finite element. 

@param dnr The local shape function derivatives at the point of interest. These
are provided as pointers to inbuilt arrays of type 'double64'. 
Whether a third argument is supplied or not determines whether a 2D
or a 3D Jacobian is output.  

@return The Jacobian matrix is returned into the protected matrix JAC.
*/
void  FiniteElement::Jacobian( const vector<double64>& dnr ) // 1D
 {
    JAC.Resize(dim,dim);
    JAC(0,0) = static_cast<double64>(0.0);
    for ( size_t j=0; j<npe; j++ ) JAC(0,0) += dnr[j] * XY(j,0);
      
 } // end Jacobian (1D)




void  FiniteElement::Jacobian( const vector<double64>& dnr, const vector<double64>& dns ) // 2D
 {
    JAC.Resize(dim,dim);

    for ( size_t i=0; i<dim; i++ ) {
        JAC(0,i) = JAC(1,i) = static_cast<double64>(0.0);
        for ( size_t j=0; j<npe; j++ )
          {
             JAC(0,i) += dnr[j] * XY(j,i);
             JAC(1,i) += dns[j] * XY(j,i);
          }
      }
      
 } // end Jacobian (2D)


void FiniteElement::Jacobian( const vector<double64>& dnr, const vector<double64>& dns, const vector<double64>& dnt )
 {
    JAC.Resize(dim,dim);

    for ( size_t i=0; i<dim; i++ ) {
        JAC(0,i) = JAC(1,i) = JAC(2,i) = static_cast<double64>(0.);
        for ( size_t j=0; j<npe; j++ )
          {
             JAC(0,i) += dnr[j] * XY(j,i);
             JAC(1,i) += dns[j] * XY(j,i);
             JAC(2,i) += dnt[j] * XY(j,i);
          }
      }
      
 } // end Jacobian (3D)


double64 FiniteElement::JacobianDeterminant()
{
    if ( dim == 3U ) 
		 return JAC(0,0) * ( JAC(1,1) * JAC(2,2) - JAC(1,2) * JAC(2,1) ) -
                JAC(0,1) * ( JAC(1,0) * JAC(2,2) - JAC(1,2) * JAC(2,0) ) +
                JAC(0,2) * ( JAC(1,0) * JAC(2,1) - JAC(1,1) * JAC(2,0) );
    else if ( dim == 2U ) 
         // compute determinant  
         return JAC(0,0)*JAC(1,1) - JAC(1,0)*JAC(0,1);
    //dim == 1
    return JAC(0,0);
}
	



/** 
    Computes the determinant J of the Jacobian matrix computed for this element/integration point
    
    @attention if the calculation fails, an ERROR is reported and some details about the Jacobian
    are printed to screen.
 
*/
double64  FiniteElement::JacobianInverse()
 {
    JINV.Resize(dim,dim);
    if ( dim == 1U ) {
         JINV(0,0) = JAC(0,0);
         return JAC(0,0);
      }
    
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( dim == 2U ) {
         // compute determinant  
         double64 detJ = JAC(0,0)*JAC(1,1) - JAC(1,0)*JAC(0,1);
    
         // inversion of J
         const double64 dum = JAC(0,0) / detJ;
         JINV(0,0)  =  JAC(1,1) / detJ;
         JINV(0,1)  = -JAC(0,1) / detJ;
         JINV(1,0)  = -JAC(1,0) / detJ;
         JINV(1,1)  =  dum;
         
         if ( detJ <= 0. ) {
              cerr <<"\n\nFiniteElement::JacobianInverse(2D): element "<< CurrentID() <<": erroneous determinant of 2D Jacobian matrix: ";
              cerr << detJ << endl;
              cerr <<"\ncaused by element of type: "<< parseFiniteElementType(csp_fem_type) << endl;
              for ( size_t i=0; i<Nodes(); i++ )
                {
                  cerr<<" Node( "<<i<<" ): ";
                  for ( size_t j=0; j<XY.Cols(); j++ )
                      cerr << XY(i,j) <<" ";
                  cerr<<endl;
                }
              csmp_error.notice( WARNING, "FiniteElement::JacobianInverse:",
                                "the value of the Jacobian is negative; check node-numbering.");

//              return fabs(detJ);
           }
      
         return detJ;
      }
      
     // 3D case
     const double64 detJ = JAC(0,0) * ( JAC(1,1) * JAC(2,2) - JAC(1,2) * JAC(2,1) ) -
                           JAC(0,1) * ( JAC(1,0) * JAC(2,2) - JAC(1,2) * JAC(2,0) ) +
                           JAC(0,2) * ( JAC(1,0) * JAC(2,1) - JAC(1,1) * JAC(2,0) );        

     if ( detJ <= 0. ) {
          std::cerr <<"\n\nFiniteElement::JacobianInverse(3D): element "<< CurrentID() <<": erroneous determinant of 3D Jacobian matrix: ";
          std::cerr << detJ << std::endl;
          cerr <<"\ncaused by element of type: "<< parseFiniteElementType(csp_fem_type) << endl;
          for ( size_t i=0; i<Nodes(); i++ )
            {
              cerr <<" Node( "<<i<<" ): ";
              for ( size_t j=0; j<XY.Cols(); j++ )
                  cerr << XY(i,j) <<" ";
              cerr <<endl;
            }
          csmp_error.notice( WARNING, "FiniteElement::JacobianInverse:",
                            "the value of the Jacobian is negative; check node-numbering.");
       
//          detJ = fabs(detJ);
       }

     double64 det1(1.0 / detJ);
     
     JINV(0,0) = ( JAC(1,1) * JAC(2,2) - JAC(1,2) * JAC(2,1) ) *  det1;
     JINV(1,0) = ( JAC(1,0) * JAC(2,2) - JAC(1,2) * JAC(2,0) ) * -det1;
     JINV(2,0) = ( JAC(1,0) * JAC(2,1) - JAC(1,1) * JAC(2,0) ) *  det1;
     JINV(0,1) = ( JAC(0,1) * JAC(2,2) - JAC(0,2) * JAC(2,1) ) * -det1;
     JINV(1,1) = ( JAC(0,0) * JAC(2,2) - JAC(0,2) * JAC(2,0) ) *  det1;
     JINV(2,1) = ( JAC(0,0) * JAC(2,1) - JAC(0,1) * JAC(2,0) ) * -det1;
     JINV(0,2) = ( JAC(0,1) * JAC(1,2) - JAC(0,2) * JAC(1,1) ) *  det1;
     JINV(1,2) = ( JAC(0,0) * JAC(1,2) - JAC(0,2) * JAC(1,0) ) * -det1;
     JINV(2,2) = ( JAC(0,0) * JAC(1,1) - JAC(0,1) * JAC(1,0) ) *  det1;

     return detJ;
     
 } // end



void  FiniteElement::Out() const
 {
    cout <<"\nFiniteElement::Out: "<< endl;
    if ( isoparametric ) cout <<"\t isoparametric ";
    if      ( element_category == LINE )    cout <<"LINE element ";
    else if ( element_category == SURFACE ) cout <<"SURFACE element ";
    else if ( element_category == VOLUME )  cout <<"VOLUME element ";
    cout <<"of ANSYS-CSMP type: "<< csp_fem_type;
    if ( uses_local_coordinates ) cout <<" and using a local coordinate system. ";
    cout << endl;
    cout <<"\n\tspatial dimension of element:        "<< dim;
    cout <<"\n\tdegree of interpolation functions:   "<< itp; 
    cout <<"\n\tnodes per element face (if applic.): "<< npf; 
    cout <<"\n\tnodes:                               "<< npe;
    cout <<"\n\tsegments per element:                "<< spe; 
    cout <<"\n\tfaces per element:                   "<< fpe; 
    cout <<"\n\tneighbors of element:                "<< epe; 
    cout <<"\n\tintegration points per element:      "<< gpe; 
    cout <<"\n\tapproximate elements that share each node:             "<< nne;
    cout <<"\n\tapproximate elements that share each constraint point: "<< cne << endl;
    
    cout <<"\nCurrent node coordinate matrix: ";  
    XY.Out();
    if ( !isoparametric ) {
         cout <<"\nCurrent test function coefficient matrix: ";  
         M.Out();
      }         
    else {
         cout <<"\nCurrent Jacobian matrix: ";  
         JAC.Out();     
         cout <<"\nCurrent inverse of Jacobian matrix: ";  
         JINV.Out();     
      }

 } // end Out
 
 
 
 
ELEMENT_DIMENSION  parseFiniteElementDimension( CSMP_FEM_TYPE etype )
  {
     if ( etype == ISOPARAMETRIC_LINEAR_BAR ) return LINE;
     if ( etype == ISOPARAMETRIC_QUADRATIC_BAR ) return LINE; 						
     if ( etype == ISOPARAMETRIC_CUBIC_BAR ) return LINE; 
     if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE ) return SURFACE;  					
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ) return SURFACE;
     if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return SURFACE; 					 
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) return SURFACE; 		
     if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE ) return SURFACE;
     if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return VOLUME; 					
     if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return VOLUME; 			    
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return VOLUME;   
     if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ) return VOLUME;
     if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ) return VOLUME;						 
     if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return VOLUME; 					
     if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return VOLUME;     				
     if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID ) return VOLUME;
     if ( etype == ISOPARAMETRIC_LINEAR_PRISM ) return VOLUME;  						
     if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ) return VOLUME;         			
     if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ) return VOLUME;           			
     if ( etype == ISOPARAMETRIC_CUBIC_PRISM ) return VOLUME;
     if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return SURFACE;                  
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return SURFACE;    
     if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return SURFACE;            
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ) return SURFACE;  
     if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return SURFACE;  		
     if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) return SURFACE;
     if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ) return VOLUME;              		
     if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return VOLUME;       		
     if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return VOLUME; 				
     if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) return VOLUME;
     if ( etype == LINEAR_CUBOID) return VOLUME;
     if ( etype == LINEAR_RECTANGLE) return SURFACE;
     if ( etype == LINEAR_BAR ) return LINE;
     if ( etype == QUADRATIC_BAR ) return LINE;
     if ( etype == CUBIC_BAR ) return LINE;
     if ( etype == LINEAR_TRIANGLE ) return SURFACE;
     if ( etype == LINEAR_TRIANGLE3D ) return SURFACE;
     if ( etype == BARYCENTRIC_LINEAR_TRIANGLE ) return SURFACE;
     if ( etype == QUADRATIC_TRIANGLE ) return SURFACE;
     if ( etype == BARYCENTRIC_QUADRATIC_TRIANGLE ) return SURFACE;
     if ( etype == CUBIC_TRIANGLE ) return SURFACE;
     if ( etype == LINEAR_TETRAHEDRON ) return VOLUME;
     if ( etype == QUADRATIC_TETRAHEDRON ) return VOLUME;
     if ( etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return VOLUME;
     if ( etype == CUBIC_TETRAHEDRON ) return VOLUME;

     cerr <<"\nparseFiniteElementDimension: Could not identify dimension of element type: "<< etype << endl;
     return static_cast<ELEMENT_DIMENSION>(0);

 } // end parseFiniteElementDimension

 
 
 
 
/// converts enum names
CSMP_FEM_TYPE  parseFiniteElementTypeEnum( int8_t etype )
  {
     if ( etype == UNKNOWN ) return UNKNOWN;
	 if (etype == LINEAR_RECTANGLE) return LINEAR_RECTANGLE;
	 if (etype == LINEAR_CUBOID) return LINEAR_CUBOID;
     if ( etype == LINEAR_BAR ) return LINEAR_BAR;    										
     if ( etype == QUADRATIC_BAR ) return QUADRATIC_BAR;
     if ( etype == CUBIC_BAR ) return CUBIC_BAR;										   
     if ( etype == LINEAR_TRIANGLE ) return LINEAR_TRIANGLE;   
     if ( etype == LINEAR_TRIANGLE3D ) return LINEAR_TRIANGLE3D; 									
     if ( etype == BARYCENTRIC_LINEAR_TRIANGLE ) return BARYCENTRIC_LINEAR_TRIANGLE; 			    	 
     if ( etype == QUADRATIC_TRIANGLE ) return QUADRATIC_TRIANGLE; 								
     if ( etype == BARYCENTRIC_QUADRATIC_TRIANGLE ) return BARYCENTRIC_QUADRATIC_TRIANGLE; 					
     if ( etype == CUBIC_TRIANGLE ) return CUBIC_TRIANGLE;
     if ( etype == LINEAR_TETRAHEDRON ) return LINEAR_TETRAHEDRON; 								
     if ( etype == QUADRATIC_TETRAHEDRON ) return QUADRATIC_TETRAHEDRON; 								
     if ( etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return BARYCENTRIC_QUADRATIC_TETRAHEDRON;                  
     if ( etype == CUBIC_TETRAHEDRON ) return CUBIC_TETRAHEDRON;
     if ( etype == ISOPARAMETRIC_LINEAR_BAR ) return ISOPARAMETRIC_LINEAR_BAR;    						
     if ( etype == ISOPARAMETRIC_QUADRATIC_BAR ) return ISOPARAMETRIC_QUADRATIC_BAR; 						
     if ( etype == ISOPARAMETRIC_CUBIC_BAR ) return ISOPARAMETRIC_CUBIC_BAR; 
     if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE ) return ISOPARAMETRIC_LINEAR_TRIANGLE;  					
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE;
     if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return ISOPARAMETRIC_QUADRATIC_TRIANGLE; 					 
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE; 		
     if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE ) return ISOPARAMETRIC_CUBIC_TRIANGLE;
     if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return ISOPARAMETRIC_LINEAR_TETRAHEDRON; 					
     if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return ISOPARAMETRIC_QUADRATIC_TETRAHEDRON; 			    
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON;   
     if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ) return ISOPARAMETRIC_CUBIC_TETRAHEDRON;
     if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ) return ISOPARAMETRIC_LINEAR_PYRAMID;						 
     if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return ISOPARAMETRIC_QUADRATIC_PYRAMID13; 					
     if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return ISOPARAMETRIC_QUADRATIC_PYRAMID14;     				
     if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID ) return ISOPARAMETRIC_CUBIC_PYRAMID;
     if ( etype == ISOPARAMETRIC_LINEAR_PRISM ) return ISOPARAMETRIC_LINEAR_PRISM;  						
     if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ) return ISOPARAMETRIC_QUADRATIC_PRISM15;         			
     if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ) return ISOPARAMETRIC_QUADRATIC_PRISM18;           			
     if ( etype == ISOPARAMETRIC_CUBIC_PRISM ) return ISOPARAMETRIC_CUBIC_PRISM;
     if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return ISOPARAMETRIC_LINEAR_QUADRILATERAL;                  
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL;    
     if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL;            
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL;  
     if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;  		
     if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) return ISOPARAMETRIC_CUBIC_QUADRILATERAL;
     if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ) return ISOPARAMETRIC_LINEAR_HEXAHEDRON;              		
     if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20;       		
     if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27; 				
     if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) return ISOPARAMETRIC_CUBIC_HEXAHEDRON;
     if ( etype == EXPERIMENTAL_ELEMENT ) return EXPERIMENTAL_ELEMENT; 

     cout <<"\nparseFiniteElementType: Could not identify element type: "<< etype << endl;
     return UNKNOWN;
 
 } // end parseFiniteElementType


 /// converts enum names into text
CSMP_FEM_TYPE  parseFiniteElementType( const std::string& etype )
  {
     if ( etype == "UNKNOWN" ) return UNKNOWN;
	 if ( etype == "LINEAR_RECTANGLE") return LINEAR_RECTANGLE;
	 if ( etype == "LINEAR_CUBOID") return LINEAR_CUBOID;
     if ( etype == "LINEAR_BAR" ) return LINEAR_BAR;
     if ( etype == "QUADRATIC_BAR" ) return QUADRATIC_BAR;
     if ( etype == "CUBIC_BAR" ) return CUBIC_BAR;
     if ( etype == "LINEAR_TRIANGLE" ) return LINEAR_TRIANGLE;
     if ( etype == "LINEAR_TRIANGLE3D" ) return LINEAR_TRIANGLE3D;
     if ( etype == "BARYCENTRIC_LINEAR_TRIANGLE" ) return BARYCENTRIC_LINEAR_TRIANGLE;
     if ( etype == "QUADRATIC_TRIANGLE" ) return QUADRATIC_TRIANGLE;
     if ( etype == "BARYCENTRIC_QUADRATIC_TRIANGLE" ) return BARYCENTRIC_QUADRATIC_TRIANGLE;
     if ( etype == "CUBIC_TRIANGLE" ) return CUBIC_TRIANGLE;
     if ( etype == "LINEAR_TETRAHEDRON" ) return LINEAR_TETRAHEDRON;
     if ( etype == "QUADRATIC_TETRAHEDRON" ) return QUADRATIC_TETRAHEDRON;
     if ( etype == "BARYCENTRIC_QUADRATIC_TETRAHEDRON" ) return BARYCENTRIC_QUADRATIC_TETRAHEDRON;
     if ( etype == "CUBIC_TETRAHEDRON" ) return CUBIC_TETRAHEDRON;
     if ( etype == "ISOPARAMETRIC_LINEAR_BAR" ) return ISOPARAMETRIC_LINEAR_BAR;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_BAR" ) return ISOPARAMETRIC_QUADRATIC_BAR;
     if ( etype == "ISOPARAMETRIC_CUBIC_BAR" ) return ISOPARAMETRIC_CUBIC_BAR;
     if ( etype == "ISOPARAMETRIC_LINEAR_TRIANGLE" ) return ISOPARAMETRIC_LINEAR_TRIANGLE;
     if ( etype == "ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE" ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_TRIANGLE" ) return ISOPARAMETRIC_QUADRATIC_TRIANGLE;
     if ( etype == "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE" ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE;
     if ( etype == "ISOPARAMETRIC_CUBIC_TRIANGLE" ) return ISOPARAMETRIC_CUBIC_TRIANGLE;
     if ( etype == "ISOPARAMETRIC_LINEAR_TETRAHEDRON" ) return ISOPARAMETRIC_LINEAR_TETRAHEDRON;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_TETRAHEDRON" ) return ISOPARAMETRIC_QUADRATIC_TETRAHEDRON;
     if ( etype == "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON" ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON;
     if ( etype == "ISOPARAMETRIC_CUBIC_TETRAHEDRON" ) return ISOPARAMETRIC_CUBIC_TETRAHEDRON;
     if ( etype == "ISOPARAMETRIC_LINEAR_PYRAMID" ) return ISOPARAMETRIC_LINEAR_PYRAMID;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_PYRAMID13" ) return ISOPARAMETRIC_QUADRATIC_PYRAMID13;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_PYRAMID14" ) return ISOPARAMETRIC_QUADRATIC_PYRAMID14;
     if ( etype == "ISOPARAMETRIC_CUBIC_PYRAMID" ) return ISOPARAMETRIC_CUBIC_PYRAMID;
     if ( etype == "ISOPARAMETRIC_LINEAR_PRISM" ) return ISOPARAMETRIC_LINEAR_PRISM;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_PRISM15" ) return ISOPARAMETRIC_QUADRATIC_PRISM15;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_PRISM18" ) return ISOPARAMETRIC_QUADRATIC_PRISM18;
     if ( etype == "ISOPARAMETRIC_CUBIC_PRISM" ) return ISOPARAMETRIC_CUBIC_PRISM;
     if ( etype == "ISOPARAMETRIC_LINEAR_QUADRILATERAL" ) return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
     if ( etype == "ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL" ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL" ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL;
     if ( etype == "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL" ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9" ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;
     if ( etype == "ISOPARAMETRIC_CUBIC_QUADRILATERAL" ) return ISOPARAMETRIC_CUBIC_QUADRILATERAL;
     if ( etype == "ISOPARAMETRIC_LINEAR_HEXAHEDRON" ) return ISOPARAMETRIC_LINEAR_HEXAHEDRON;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20" ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20;
     if ( etype == "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27" ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27;
     if ( etype == "ISOPARAMETRIC_CUBIC_HEXAHEDRON" ) return ISOPARAMETRIC_CUBIC_HEXAHEDRON;
     if ( etype == "EXPERIMENTAL_ELEMENT" ) return EXPERIMENTAL_ELEMENT;

     cout <<"\nparseFiniteElementType: Could not identify element type: "<< etype << endl;
     return UNKNOWN;

 } // end parseFiniteElementType

 
 
 
 
 
 /// converts enum names into text
 const char* parseFiniteElementType( int8_t etype )
  {
     if ( etype == UNKNOWN ) return "UNKNOWN";
	   if ( etype == LINEAR_RECTANGLE) return "LINEAR_RECTANGLE";
	   if ( etype == LINEAR_CUBOID) return "LINEAR_CUBOID";
     if ( etype == LINEAR_BAR ) return "LINEAR_BAR";    										                                    // BAR_2     =2,
     if ( etype == QUADRATIC_BAR ) return "QUADRATIC_BAR";
     if ( etype == CUBIC_BAR ) return "CUBIC_BAR";										                                            // BAR_4  
     if ( etype == LINEAR_TRIANGLE ) return "LINEAR_TRIANGLE";   
     if ( etype == LINEAR_TRIANGLE3D ) return "LINEAR_TRIANGLE3D"; 									                                // TRI_3     =8, 
     if ( etype == BARYCENTRIC_LINEAR_TRIANGLE ) return "BARYCENTRIC_LINEAR_TRIANGLE"; 			    		                        // TRI_3_X   =9, 
     if ( etype == QUADRATIC_TRIANGLE ) return "QUADRATIC_TRIANGLE"; 								                                //  2D & 3D TRI_6    =10, 
     if ( etype == BARYCENTRIC_QUADRATIC_TRIANGLE ) return "BARYCENTRIC_QUADRATIC_TRIANGLE"; 					                    // TRI_6_X  =11,
     if ( etype == CUBIC_TRIANGLE ) return "CUBIC_TRIANGLE";
     if ( etype == LINEAR_TETRAHEDRON ) return "LINEAR_TETRAHEDRON"; 								                                // TETRA_4   =4, 
     if ( etype == QUADRATIC_TETRAHEDRON ) return "QUADRATIC_TETRAHEDRON"; 								                            // TETRA_10  =5,
     if ( etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return "BARYCENTRIC_QUADRATIC_TETRAHEDRON";                                  // PYRA_5   =18,
     if ( etype == CUBIC_TETRAHEDRON ) return "CUBIC_TETRAHEDRON";
     if ( etype == ISOPARAMETRIC_LINEAR_BAR ) return "ISOPARAMETRIC_LINEAR_BAR";    						                        // BAR_2     =2,
     if ( etype == ISOPARAMETRIC_QUADRATIC_BAR ) return "ISOPARAMETRIC_QUADRATIC_BAR"; 						                        // BAR_3     =3,
     if ( etype == ISOPARAMETRIC_CUBIC_BAR ) return "ISOPARAMETRIC_CUBIC_BAR"; 
     if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE ) return "ISOPARAMETRIC_LINEAR_TRIANGLE";  					                    // TRI_3     =8, 
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ) return "ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE"; 	         	// TRI_3_X   =9,
     if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return "ISOPARAMETRIC_QUADRATIC_TRIANGLE"; 					                //  2D & 3D TRI_6    =10, 
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) return "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE"; 		    // TRI_6_X  =11,
     if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE ) return "ISOPARAMETRIC_CUBIC_TRIANGLE";
     if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return "ISOPARAMETRIC_LINEAR_TETRAHEDRON"; 					                // TETRA_4   =4, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return "ISOPARAMETRIC_QUADRATIC_TETRAHEDRON"; 			                    // TETRA_10  =5,
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON";      // TETRA_11
     if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ) return "ISOPARAMETRIC_CUBIC_TETRAHEDRON";
     if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ) return "ISOPARAMETRIC_LINEAR_PYRAMID";						                    // PYRA_5   =18,  
     if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return "ISOPARAMETRIC_QUADRATIC_PYRAMID13"; 					                // PYRA_13  =24, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return "ISOPARAMETRIC_QUADRATIC_PYRAMID14";     				                // PYRA_14  =22, 
     if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID ) return "ISOPARAMETRIC_CUBIC_PYRAMID";
     if ( etype == ISOPARAMETRIC_LINEAR_PRISM ) return "ISOPARAMETRIC_LINEAR_PRISM";  						                        // PENTA_6  =12, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ) return "ISOPARAMETRIC_QUADRATIC_PRISM15";         			                    // PENTA_15 =13, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ) return "ISOPARAMETRIC_QUADRATIC_PRISM18";           			                // PENTA_18 =21,
     if ( etype == ISOPARAMETRIC_CUBIC_PRISM ) return "ISOPARAMETRIC_CUBIC_PRISM";
     if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return "ISOPARAMETRIC_LINEAR_QUADRILATERAL";                                // QUAD_4   =14, 
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return "ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL";        // QUAD_4_X =15, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL";                          // QUAD_8   =16, 
     if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ) return "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL";  // QUAD_8_X =17, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9";  			            // QUAD_9   =19,
     if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) return "ISOPARAMETRIC_CUBIC_QUADRILATERAL";
     if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ) return "ISOPARAMETRIC_LINEAR_HEXAHEDRON";              		                // HEXA_8    =6, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20";       		                // HEXA_20   =7, 
     if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27"; 				            // HEXA_27   =5,
     if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) return "ISOPARAMETRIC_CUBIC_HEXAHEDRON";
     if ( etype == EXPERIMENTAL_ELEMENT ) return "EXPERIMENTAL_ELEMENT"; 

     cout <<"\nparseFiniteElementType: Could not identify element type: "<< etype << endl;
     return "UNKNOWN";
 
 } // end parseFiniteElementType


 FV_FACET_TYPE parseFacetType( const std::string& ftype )
 {
     if (ftype == "POINT_FACET") return POINT_FACET;
     if (ftype == "UNIT_LINEAR_FACET") return UNIT_LINEAR_FACET;
     if (ftype == "TRIANGULAR_FACET") return TRIANGULAR_FACET;
     if (ftype == "QUADRILATERAL_FACET") return QUADRILATERAL_FACET;
     cout << "\nparseFacetType: Could not identify facet type " << ftype << '\n';
     return static_cast<FV_FACET_TYPE>(0);
 }

 const char* parseFacetType( int8_t ftype )
 {
     switch (ftype)
     {
         case POINT_FACET: return "POINT_FACET";
         case UNIT_LINEAR_FACET: return "UNIT_LINEAR_FACET";
         case TRIANGULAR_FACET: return "TRIANGULAR_FACET";
         case QUADRILATERAL_FACET: return "QUADRILATERAL_FACET";
         default:
             cout << "\nparseFacetType: Could not identify facet type " << ftype << '\n';
             return "(unknown)";
     }
 } // end parseFacetType


 
 } // end namespace csp
 
 
 
 
