#include "FiniteElement.h"
#include "CSMP_definitions.h"
#include "ErrorHandler.h"

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
   object_id_( InitialID() )
 {
   cout <<"\n"<<"FiniteElement::FiniteElement: default constructor called."<< endl;
   Out();
 }


FiniteElement::FiniteElement( CSMP_FEM_TYPE csp_fem_type,
                              bool isoparametric, 
                              bool uses_local_coordinates,  
                              uint32_t order_of_shape_functions )
 : XY(3,2), M(3,3), JAC(2,2), JINV(2,2), 
   isoparametric(isoparametric), 
   uses_local_coordinates(uses_local_coordinates),
   order_of_shape_functions(order_of_shape_functions),
   element_category(SURFACE),
   csp_fem_type(csp_fem_type),
   object_id_( InitialID() )
 {
 }


void  FiniteElement::Isoparametric( bool isoparam ) { isoparametric = isoparam; }
void  FiniteElement::UsesLocalCoordinates( bool uses ) { uses_local_coordinates = uses; }


// methods



void    FiniteElement::CurrentID( size_t id ) noexcept { object_id_ = id; }
size_t  FiniteElement::CurrentID() const noexcept { return object_id_; }

/// initalizing to a value that makes sure that ID does not equal initial element idx
size_t  FiniteElement::InitialID() noexcept { return std::numeric_limits<uint32_t>::max(); }

bool  FiniteElement::Isoparametric() const noexcept { return isoparametric; }

bool  FiniteElement::UsesLocalCoordinates() const noexcept { return uses_local_coordinates; }

bool  FiniteElement::IsLine()    const noexcept
 { if ( element_category == LINE )   return true; return false; }
 
bool  FiniteElement::IsSurface() const noexcept
 { if ( element_category == SURFACE ) return true; return false; }
 
bool  FiniteElement::IsVolume()  const noexcept    
 { if ( element_category == VOLUME )  return true; return false; }
 

void  FiniteElement::ElementType( CSMP_FEM_TYPE etype ) noexcept { csp_fem_type = etype; }

CSMP_FEM_TYPE  FiniteElement::ElementType() const noexcept { return csp_fem_type; }

 uint32_t  FiniteElement::OrderOfShapeFunctions() const noexcept { return order_of_shape_functions; }

 uint32_t  FiniteElement::Interpolation()    const noexcept { return itp; }
 uint32_t  FiniteElement::Dim()              const noexcept { return dim; }
 uint32_t  FiniteElement::Nodes()            const noexcept { return npe; }
 uint32_t  FiniteElement::Segments()         const noexcept { return spe; }
 uint32_t  FiniteElement::Faces()            const noexcept { return fpe; }
 uint32_t  FiniteElement::Neighbors()        const noexcept { return epe; }
 uint32_t  FiniteElement::NodesPerFace( uint32_t ) const noexcept { return npf; }
// uint32_t  FiniteElement::IntegrationPointNeighbors() const { return cne; }
 uint32_t  FiniteElement::IntegrationPoints() const noexcept { return gpe; }

 double  FiniteElement::XYZ( uint32_t i, uint32_t j ) const { return XY(i,j); }
 void    FiniteElement::XYZ( uint32_t i, uint32_t j, double val ) { XY(i,j) = val; }



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


double FiniteElement::AspectRatio()
  {
     InstructUser("FiniteElement::AspectRatio()");
     cout <<"\ncalled by object: "<< object_id_ << endl;
     return 0.0;
  }
  
double FiniteElement::InnerRadius()
  {
     InstructUser("FiniteElement::InnerRadius()");
     cout <<"\ncalled by object: "<< object_id_ << endl;
     return 0.0;
  }
  
void FiniteElement::EdgeLengths( vector<double>& vec )
  {
     InstructUser("FiniteElement::EdgeLengths(vector<double>)");
     copy( vec.begin(), vec.end(), ostream_iterator<double>(cout, " ") );
     cout <<"\ncalled by object: "<< object_id_ << endl;
  }



void FiniteElement::NodesOfSegment( uint32_t sid, vector<uint32_t>& snids ) const
  {
     cout <<"\nFiniteElement::NodesOfSegment: Returns the local node ID numbers of ";
     cout <<"the nodes which constitute the element segment with the entered (local) ID number. "<< endl;
     cout <<"\ncalled by object: "<< object_id_ <<" for face "<< sid << endl;
     copy( snids.begin(), snids.end(), ostream_iterator<uint32_t>(cout, " ") );
     throw invalid_argument("FiniteElement::NodesOfSegment");
  }


vector<uint32_t>  FiniteElement::NodesOfFace( uint32_t fid ) const
  {
     cout <<"\nFiniteElement::NodesOfFace: Returns the local node ID numbers of ";
     cout <<"the nodes which constitute the element face with the entered ID number. ";
     cout <<"In triangular and tetrahedral elements the faces lie opposite of ";
     cout <<"the nodes with the same ID." << endl;
     cout <<"\ncalled by object: "<< object_id_ <<" for face "<< fid << endl;
     throw invalid_argument("FiniteElement::NodesOfFace");
  }


std::vector<uint32_t> FiniteElement::CornerNodesOfFace( uint32_t face_id ) const
 {
     cout <<"\nFiniteElement::CornerNodesOfFace: Returns the local node ID numbers of ";
     cout <<"the corner nodes of the element face with the entered ID number. ";
     cout <<"In triangular and tetrahedral elements the faces lie opposite of ";
     cout <<"the nodes with the same ID." << endl;
     cout <<"\ncalled by object: "<< object_id_ <<" for face "<< face_id << endl;
     throw invalid_argument("FiniteElement::CornerNodesOfFace");
     return vector<uint32_t>{};
 }


vector<uint32_t>  FiniteElement::NodesConnectedTo( uint32_t node_id ) const
  {
     cout <<"\nFiniteElement::NodesConnectedTo: Returns the local node ID numbers of ";
     cout <<"the corner nodes of the element face with the entered ID number. ";
     cout <<"In triangular and tetrahedral elements the faces lie opposite of ";
     cout <<"the nodes with the same ID." << endl;
     cout <<"\ncalled by object: "<< object_id_ <<" for node "<< node_id << endl;
     throw invalid_argument("FiniteElement::NodesConnectedTo");
     return vector<uint32_t>{};
  }

    
void FiniteElement::IntegraldNdN( DenseMatrix<DM_MIN>& DM )
  {
     InstructUser("FiniteElement::IntegraldNdN( DenseMatrix<DM_MIN>& M )");
     DM.Out();
     cout <<"\ncalled by object: "<< object_id_ << endl;
     throw invalid_argument("FiniteElement::IntegraldNdN");
  }
  

  
// other information

void FiniteElement::CornerNodes( vector<uint32_t>& ids ) const
 {
    InstructUser("FiniteElement::CornerNodes");
    copy( ids.begin(), ids.end(), ostream_iterator<uint32_t>(cout, " ") );
    throw invalid_argument("FiniteElement::CornerNodes");
 }
 
 
void FiniteElement::MidSideNodes( vector<uint32_t>& ids ) const
 {
    InstructUser("FiniteElement::MidSideNodes");
    copy( ids.begin(), ids.end(), ostream_iterator<uint32_t>(cout, " ") );
    throw invalid_argument("FiniteElement::MidSideNodes");
 }

void FiniteElement::InteriorNodes( vector<uint32_t>& ids ) const
 {
    InstructUser("FiniteElement::InteriorNodes");
    copy( ids.begin(), ids.end(), ostream_iterator<uint32_t>(cout, " ") );
    throw invalid_argument("FiniteElement::InteriorNodes");
 }

uint32_t FiniteElement::CornerNodes() const
 {
    InstructUser("FiniteElement::CornerNodes");
    throw invalid_argument("FiniteElement::CornerNodes");
    return 0;
 }
 
 
uint32_t FiniteElement::MidSideNodes() const
 {
    InstructUser("FiniteElement::MidSideNodes");
    throw invalid_argument("FiniteElement::MidSideNodes");
    return 0;
 }

uint32_t FiniteElement::InteriorNodes() const
 {
    InstructUser("FiniteElement::InteriorNodes");
    throw invalid_argument("FiniteElement::InteriorNodes");
    return 0;
 }


CSMP_FEM_TYPE FiniteElement::ElementTypeOfFace( uint32_t ) const
  {
     /// @todo (2-P) Remove typid by name fct
     InstructUser("FiniteElement::ElementTypeOfFace");
     cout <<"\nDerived Element: "<< typeid(this).name() << endl;
     cout <<"\nThis method must be defined for the FE element implementation which you are using"<< endl;
     throw invalid_argument("FiniteElement::ElementTypeOfFace");
     
     return UNKNOWN;
  }


CSMP_FEM_TYPE FiniteElement::ElementTypeOfSegment( uint32_t ) const
  {
     /// @todo (2-P) Remove typid by name fct
     InstructUser("FiniteElement::ElementTypeOfSegment");
     cout <<"\nDerived Element: "<< typeid(this).name() << endl;
     cout <<"\nThis method must be defined for the FE element implementation which you are using"<< endl;
     throw invalid_argument("FiniteElement::ElementTypeOfSegment");
     
     return UNKNOWN;
  }



vector<double>  FiniteElement::UnitNormal() const
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
void  FiniteElement::UnitNormalToFace( uint32_t, std::vector<double>& ) const
 {
     InstructUser("FiniteElement::UnitNormalToFace");
     cout <<"\nThis method is not defined for the FE element implementation which you are using"<< endl;
     cout <<"\nyou need to overload (=define) the method in the desired FiniteElement subclass"<< endl;
     throw invalid_argument("FiniteElement::UnitNormalToFace( uint32_t,vector<double>)");
 }



/**
     This method is intended for numerically integrated element types with higher order 
     shape functions. As these elements can get deformed in physical space, it is important 
     to specify the precise location of the unit normal.
     
     @note this method cannot be constant because it modifies matrices that are data members
     of the FiniteElement class.
*/
void  FiniteElement::UnitNormalAtFaceBarycenter( uint32_t, std::vector<double>& )
 {
     InstructUser("FiniteElement::UnitNormalAtFaceBarycenter");
     cout <<"\nThis method is not defined for the FE element implementation which you are using"<< endl;
     cout <<"\nyou need to overload (=define) the method in the desired FiniteElement subclass"<< endl;
     throw invalid_argument("FiniteElement::UnitNormalAtFaceBarycenter( uint32_t,vector<double>)");
 }



double FiniteElement::WeightAtIntegrationPoint( uint32_t ) const
  {
     throw invalid_argument("FiniteElement::WeightAtIntegrationPoint");
     return 0.5 * (1. / gpe); 
  }



void  FiniteElement::IntegrationPoint( uint32_t i, std::vector<double>& xyz ) const
 { 
    cout <<"\nFiniteElement::IntegrationPoint: called for integration point "<< i;
    cout <<" and output data vector 'xyz' of size "<< xyz.size() << endl;
    cout <<"\nThis method which returns the position of the integration point ";
    cout <<" in global coordinates must be overloaeded for the element of interest."<< endl;
    throw invalid_argument("FiniteElement::IntegrationPoint");
 }




/** Extrapolation member( gauss points to nodes )

Expects input vector to be of size nvars * gauss points and in format
[ var1 @ gp1, var2 @ gp1, var3 @ gp1,
  var1 @ gp2, var2 @ gp2, var3 @ gp2 for all gauss points, with nvar == 3 ]

  */
void  FiniteElement::ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                                 const vector<double>&, 
                                                                 vector<double>& ) const
 {
     InstructUser("FiniteElement::ExtrapolateIntegrationPointVariableToNodes");
     cout <<"\nThis method must be defined for the FE element which you are using"<< endl;
     cout <<"\nSupplied data for N-variables = "<< nvars << endl;
     throw invalid_argument("FiniteElement::ExtrapolateIntegrationPointVariableToNodes");
 } 




double FiniteElement::Volume()
 {
    InstructUser("FiniteElement::Volume");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    throw invalid_argument("FiniteElement::Volume");
    return 0.0;
 } 


void   FiniteElement::N( vector<double>& N, const vector<double>& xyz )
 {
    InstructUser("FiniteElement::N");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (N vector<double>, xyz vector<double>):"<< endl;
    copy( N.begin(), N.end(), ostream_iterator<double>(cout, " ") );
    for ( auto i{0U}; i<xyz.size(); i++ ) cout << xyz[i] <<" ";
    cout << endl;
    throw invalid_argument("FiniteElement::N");
 } 


void   FiniteElement::N_AtIntegrationPoint( uint32_t ip, vector<double>& N )
 {
    InstructUser("FiniteElement::N_AtIntegrationPoint");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (ip, N vector<double>):"<< ip << endl;
    copy( N.begin(), N.end(), ostream_iterator<double>(cout, " ") );
    throw invalid_argument("FiniteElement::N_AtIntegrationPoint");
 } 


void   FiniteElement::N_AtBaryCenter( vector<double>& N )
 {
    InstructUser("FiniteElement::N_AtBaryCenter");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (ip, N vector<double>):"<< endl;
    copy( N.begin(), N.end(), ostream_iterator<double>(cout, " ") );
    throw invalid_argument("FiniteElement::N_AtBaryCenter");
 }

void FiniteElement::Integral_dNT_K_dN(DenseMatrix<DM_MIN>& DM, DenseMatrix<DM_MIN>& K)
{
	InstructUser("FiniteElement::Integral_dNT_K_dN");
	cout << "\nThis method is not defined for the FE element type which you are using" << endl;
	cout << "\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M):" << endl;
	DM.Out();
	throw invalid_argument("FiniteElement::dN");
}


void   FiniteElement::dN( DenseMatrix<DM_MIN>& DM ) // coefficients
 {
    InstructUser("FiniteElement::dN");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M):"<< endl;
    DM.Out();
    throw invalid_argument("FiniteElement::dN");
 } 
    
    
double FiniteElement::dN_At( DenseMatrix<DM_MIN>& DM, const vector<double>& xyz  ) // derivatives at 'xy'
 {
    InstructUser("FiniteElement::dN_At");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M, vector<double> xyz):"<< endl;
    DM.Out();
    for ( auto i{0U}; i<xyz.size(); i++ ) cout << xyz[i] <<" ";
    cout << endl;
    throw invalid_argument("FiniteElement::dN_At");
    return 0.0;
 } 


double FiniteElement::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& DM, uint32_t gauss_point )
 {
    InstructUser("FiniteElement::dN_AtIntegrationPoint");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (gp, dof, DenseMatrix<DM_MIN> M):"<< endl;
    cout << gauss_point <<", "<< endl;
    DM.Out();
    throw invalid_argument("FiniteElement::dN_AtIntegrationPoint");
    return 0.0;
 } 


double FiniteElement::dN_AtNode( DenseMatrix<DM_MIN>& DM, uint32_t node )
 {
    InstructUser("FiniteElement::dN_AtNode");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (node, dof, DenseMatrix<DM_MIN> M): ";
    cout << node <<", "<< endl;
    DM.Out();
    throw invalid_argument("FiniteElement::dN_AtNode");
    return 0.0;
 } 


double FiniteElement::dN_AtBarycenter( DenseMatrix<DM_MIN>& DM )
 {
    InstructUser("FiniteElement::dN_AtBarycenter");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M): " << endl;
    DM.Out();
    throw invalid_argument("FiniteElement::dN_AtBarycenter");
    return 0.0;
 } 
    

void   FiniteElement::IntegralNN( DenseMatrix<DM_MIN>& DM )
 {
    InstructUser("FiniteElement::IntegralNN");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"\nMethod arguments for single-element mesh output, input (DenseMatrix<DM_MIN> M):"<< endl;
    DM.Out();
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



void   FiniteElement::JacobianAtIntegrationPoint( uint32_t ip )
 {
    InstructUser("FiniteElement::JacobianAtIntegrationPoint");
    cout <<"\nThis method is not defined for the FE element type which you are using"<< endl;
    cout <<"Just write a new method which calculates the Jacobian matrix for the element in question.";
    cout <<"\nMethod arguments for single-element mesh output, input (integration point): "<< endl;
    cout <<"'"<< ip <<"', '"<< endl;
    throw invalid_argument("FiniteElement::JacobianAtIntegrationPoint");
 }




    // 1D
void   FiniteElement::Nr( double r, vector<double>& NR ) const
 {
    InstructUser("FiniteElement::Nr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinate r: "<< r << endl;
    copy( NR.begin(), NR.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::Nr");
 }



void   FiniteElement::dNr( double r, vector<double>& vDNR ) const
 {
    InstructUser("FiniteElement::dNr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinate r: "<< r << endl;
    copy( vDNR.begin(), vDNR.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::dNr");
 }


    // 2D
void   FiniteElement::Nrs( double r, double s, vector<double>& NRS ) const
 {
    InstructUser("FiniteElement::Nrs");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s: "<< r <<" "<< s << endl;
    copy( NRS.begin(), NRS.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::Nrs");
 }

void   FiniteElement::dNr( double r, double s, vector<double>& vDNR ) const
 {
    InstructUser("FiniteElement::dNr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s: "<< r <<" "<< s << endl;
    copy( vDNR.begin(), vDNR.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::dNr");
 }


void   FiniteElement::dNs( double r, double s, vector<double>& vDNS ) const
 {
    InstructUser("FiniteElement::dNs");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s: "<< r <<" "<< s << endl;
    copy( vDNS.begin(), vDNS.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::dNs");
 }


    // 3D
void   FiniteElement::Nrst( double r, double s, double t, vector<double>& vNRST ) const
 {
    InstructUser("FiniteElement::Nrst");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    copy( vNRST.begin(), vNRST.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::Nrst");
 }




void   FiniteElement::dNr( double r,  double s, double t, vector<double>& vDNR ) const
 {
    InstructUser("FiniteElement::dNr");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    copy( vDNR.begin(), vDNR.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::dNr");
 }


void   FiniteElement::dNs( double r,  double s, double t, vector<double>& vDNS ) const
 {
    InstructUser("FiniteElement::dNs");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    copy( vDNS.begin(), vDNS.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::dNs");
 }


void   FiniteElement::dNt( double r,  double s, double t, vector<double>& vDNT ) const
 {
    InstructUser("FiniteElement::dNt");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< r <<" "<< s <<" "<< t << endl;
    copy( vDNT.begin(), vDNT.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::dNt");
 }


void FiniteElement::JacobianAt( const std::vector<double>& rst )
{
    InstructUser("FiniteElement::JacobianAt");
    cout <<"\nThis method is not defined for the FE element type which you are using."<< endl;
    if ( !uses_local_coordinates || !isoparametric )
      cout <<"\tThis function uses local coordinates which do not exist for the current element type."<< endl;
    cout <<"\nMethod arguments: "<< endl;
    cout <<"\nLocal coordinates r,s,t: "<< rst[0] <<" "<< rst[1] <<" "<< rst[2] << endl;
    copy( DNT.begin(), DNT.end(), ostream_iterator<double>(cout, " ") );
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
    copy( DNT.begin(), DNT.end(), ostream_iterator<double>(cout, " ") );
    cout << endl;
    throw invalid_argument("FiniteElement::ReferenceCoordinates");
}
	

/**

Compute the Jacobian matrix for the local coordinate point 'rst' in the 
finite element. 

@param dnr The local shape function derivatives at the point of interest. These
are provided as pointers to inbuilt arrays of type 'double'. 
Whether a third argument is supplied or not determines whether a 2D
or a 3D Jacobian is output.  

The Jacobian matrix is returned into the protected matrix JAC.
*/
void  FiniteElement::Jacobian( const vector<double>& dnr ) // 1D
 {
    JAC.Resize(dim,dim);
    JAC(0,0) = static_cast<double>(0.0);
    for ( auto j{0U}; j<npe; j++ ) JAC(0,0) += dnr[j] * XY(j,0);
      
 } // end Jacobian (1D)




void  FiniteElement::Jacobian( const vector<double>& dnr, const vector<double>& dns ) // 2D
 {
    JAC.Resize(dim,dim);

    for ( uint32_t i{0U}; i<dim; i++ ) {
        JAC(0,i) = JAC(1,i) = static_cast<double>(0.0);
        for ( uint32_t j{0U}; j<npe; j++ )
          {
             JAC(0,i) += dnr[j] * XY(j,i);
             JAC(1,i) += dns[j] * XY(j,i);
          }
      }
      
 } // end Jacobian (2D)


void FiniteElement::Jacobian( const vector<double>& dnr, const vector<double>& dns, const vector<double>& dnt )
 {
    JAC.Resize(dim,dim);

    for ( uint32_t i{0U}; i<dim; i++ ) {
        JAC(0,i) = JAC(1,i) = JAC(2,i) = static_cast<double>(0.);
        for ( uint32_t j{0U}; j<npe; j++ )
          {
             JAC(0,i) += dnr[j] * XY(j,i);
             JAC(1,i) += dns[j] * XY(j,i);
             JAC(2,i) += dnt[j] * XY(j,i);
          }
      }
      
 } // end Jacobian (3D)


double FiniteElement::JacobianDeterminant()
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
double  FiniteElement::JacobianInverse()
 {
    JINV.Resize(dim,dim);
    if ( dim == 1U ) {
         JINV(0,0) = JAC(0,0);
         return JAC(0,0);
      }
    
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( dim == 2U ) {
         // compute determinant  
         const double detJ = JAC(0,0)*JAC(1,1) - JAC(1,0)*JAC(0,1);
    
         // inversion of J
         const double dum = JAC(0,0) / detJ;
         JINV(0,0)  =  JAC(1,1) / detJ;
         JINV(0,1)  = -JAC(0,1) / detJ;
         JINV(1,0)  = -JAC(1,0) / detJ;
         JINV(1,1)  =  dum;
         
         if ( detJ <= 0. ) {
              cerr <<"\n\nFiniteElement::JacobianInverse(2D): element "<< CurrentID() <<": erroneous determinant of 2D Jacobian matrix: ";
              cerr << std::defaultfloat << detJ << endl;
              cerr <<"\ncaused by element of type: "<< parseFiniteElementType(csp_fem_type) << endl;
              for ( uint32_t i{0U}; i<Nodes(); i++ )
                {
                  cerr<<" Node("<<i<<"): "<< std::scientific;
                  for ( uint32_t j{0U}; j<XY.Cols(); j++ )
                      cerr << XY(i,j) <<" ";
                  cerr << std::defaultfloat << endl;
                }
              cerr <<"(are the nodes perhaps numbered clockwise?), node coordinate matrix:";
              string file_name{ parseFiniteElementType( ElementType() ) };
              // to print a scalar, the data matrix only needs 1 row
              DenseMatrix<DM_MIN> DATA( 1, Nodes() );
              // values increase linearly from first to last node (so that node numbering direction can be seen)
              for ( auto j{0u}; j<Nodes(); ++j ) DATA(0,j) = j;
              OutputNodeDataToVTK( file_name.c_str(), "error_code", DATA );
              csmp_error.Note( WARNING, "FiniteElement::JacobianInverse:",
                                "the value of the Jacobian is negative; check node-numbering.");
           }
      
         return detJ;
      }
      
     // 3D case
     const double detJ = JAC(0,0) * ( JAC(1,1) * JAC(2,2) - JAC(1,2) * JAC(2,1) ) -
                         JAC(0,1) * ( JAC(1,0) * JAC(2,2) - JAC(1,2) * JAC(2,0) ) +
                         JAC(0,2) * ( JAC(1,0) * JAC(2,1) - JAC(1,1) * JAC(2,0) );        

     if ( detJ <= 0. ) {
          std::cerr <<"\n\nFiniteElement::JacobianInverse(3D): element "<< CurrentID() <<": erroneous determinant of 3D Jacobian matrix: ";
          std::cerr << std::defaultfloat << detJ << std::endl;
          cerr <<"\ncaused by element of type: "<< parseFiniteElementType(csp_fem_type) << endl;
          for ( auto i{0U}; i<Nodes(); i++ )
            {
              cerr <<" Node("<<i<<"): "<< std::scientific;
              for ( auto j{0U}; j<XY.Cols(); j++ )
                  cerr << XY(i,j) <<" ";
              cerr << std::defaultfloat << endl;
            }
          csmp_error.Note( WARNING, "FiniteElement::JacobianInverse:",
                            "the value of the Jacobian is negative; check node-numbering.");
       
//          detJ = fabs(detJ);
       }

     double det1(1.0 / detJ);
     
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
 
 
 
/**
      Ordered in sequence of most common queries to speed up
*/
CELL_SHAPE  parseFiniteElementDimension( CSMP_FEM_TYPE etype ) noexcept
  {
    switch (etype)
    {
        // LINE elements
        case ISOPARAMETRIC_LINEAR_BAR:
        case ISOPARAMETRIC_QUADRATIC_BAR:
        case ISOPARAMETRIC_CUBIC_BAR:
        case LINEAR_BAR:
        case QUADRATIC_BAR:
        case CUBIC_BAR:
            return LINE;

        // SURFACE elements
        case ISOPARAMETRIC_LINEAR_TRIANGLE:
        case ISOPARAMETRIC_QUADRATIC_TRIANGLE:
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE:
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE:
        case ISOPARAMETRIC_CUBIC_TRIANGLE:
        case ISOPARAMETRIC_LINEAR_QUADRILATERAL:
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL:
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL:
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL:
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9:
        case ISOPARAMETRIC_CUBIC_QUADRILATERAL:
        case LINEAR_TRIANGLE:
        case LINEAR_TRIANGLE3D:
        case BARYCENTRIC_LINEAR_TRIANGLE:
        case QUADRATIC_TRIANGLE:
        case BARYCENTRIC_QUADRATIC_TRIANGLE:
        case CUBIC_TRIANGLE:
        case LINEAR_RECTANGLE:
        case POLYGONAL_ELEMENT:
            return SURFACE;

        // VOLUME elements
        case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
        case LINEAR_TETRAHEDRON:
        case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON:
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON:
        case ISOPARAMETRIC_CUBIC_TETRAHEDRON:
        case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20:
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27:
        case ISOPARAMETRIC_CUBIC_HEXAHEDRON:
        case LINEAR_CUBOID:
        case ISOPARAMETRIC_LINEAR_PRISM:
        case ISOPARAMETRIC_QUADRATIC_PRISM15:
        case ISOPARAMETRIC_QUADRATIC_PRISM18:
        case ISOPARAMETRIC_CUBIC_PRISM:
        case ISOPARAMETRIC_LINEAR_PYRAMID:
        case ISOPARAMETRIC_QUADRATIC_PYRAMID13:
        case ISOPARAMETRIC_QUADRATIC_PYRAMID14:
        case ISOPARAMETRIC_CUBIC_PYRAMID:
        case POLYHEDRAL_ELEMENT:
        case QUADRATIC_TETRAHEDRON:
        case BARYCENTRIC_QUADRATIC_TETRAHEDRON:
        case CUBIC_TETRAHEDRON:
            return VOLUME;

        // POINT elements
        case ZERO_DIMENSIONAL_FACE:
        case POINT_ELEMENT:
            return POINT;

        default:
            cerr << "\nparseFiniteElementDimension: Could not identify dimension of element type: " << etype << endl;
            return static_cast<CELL_SHAPE>(UNSPECIFIED);
    }

 } // end parseFiniteElementDimension

 
 
 
 
/// converts enum names
CSMP_FEM_TYPE  parseFiniteElementTypeEnum( int8_t etype ) noexcept
  {
    switch (etype)
    {
        // LINE elements
        case LINEAR_BAR: 
        case QUADRATIC_BAR:
        case CUBIC_BAR:
        case ISOPARAMETRIC_LINEAR_BAR:
        case ISOPARAMETRIC_QUADRATIC_BAR:
        case ISOPARAMETRIC_CUBIC_BAR:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // TRIANGLE elements
        case LINEAR_TRIANGLE:
        case LINEAR_TRIANGLE3D:
        case BARYCENTRIC_LINEAR_TRIANGLE:
        case QUADRATIC_TRIANGLE:
        case BARYCENTRIC_QUADRATIC_TRIANGLE:
        case CUBIC_TRIANGLE:
        case ISOPARAMETRIC_LINEAR_TRIANGLE:
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE:
        case ISOPARAMETRIC_QUADRATIC_TRIANGLE:
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE:
        case ISOPARAMETRIC_CUBIC_TRIANGLE:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // TETRAHEDRON elements
        case LINEAR_TETRAHEDRON:
        case QUADRATIC_TETRAHEDRON:
        case BARYCENTRIC_QUADRATIC_TETRAHEDRON:
        case CUBIC_TETRAHEDRON:
        case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
        case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON:
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON:
        case ISOPARAMETRIC_CUBIC_TETRAHEDRON:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // PYRAMID elements
        case ISOPARAMETRIC_LINEAR_PYRAMID:
        case ISOPARAMETRIC_QUADRATIC_PYRAMID13:
        case ISOPARAMETRIC_QUADRATIC_PYRAMID14:
        case ISOPARAMETRIC_CUBIC_PYRAMID:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // PRISM elements
        case ISOPARAMETRIC_LINEAR_PRISM:
        case ISOPARAMETRIC_QUADRATIC_PRISM15:
        case ISOPARAMETRIC_QUADRATIC_PRISM18:
        case ISOPARAMETRIC_CUBIC_PRISM:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // QUADRILATERAL elements
        case LINEAR_RECTANGLE:
        case ISOPARAMETRIC_LINEAR_QUADRILATERAL:
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL:
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL:
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL:
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9:
        case ISOPARAMETRIC_CUBIC_QUADRILATERAL:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // HEXAHEDRON elements
        case LINEAR_CUBOID:
        case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20:
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27:
        case ISOPARAMETRIC_CUBIC_HEXAHEDRON:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // POINT elements
        case ZERO_DIMENSIONAL_FACE:
        case POINT_ELEMENT:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // POLYGON / POLYHEDRAL
        case POLYGONAL_ELEMENT:
        case POLYHEDRAL_ELEMENT:
            return static_cast<CSMP_FEM_TYPE>(etype);

        // experimental / other
        case EXPERIMENTAL_ELEMENT:
            return static_cast<CSMP_FEM_TYPE>(etype);

        default:
            cerr << "\nparseFiniteElementTypeEnum: Could not identify element type: " 
                 << static_cast<int>(etype) << endl;
            return UNKNOWN;
    }
 
 } // end parseFiniteElementType







 /// converts enum names into text
CSMP_FEM_TYPE  parseFiniteElementType( const std::string& etype ) noexcept
  {
    // static map so it is initialized only once
    static const std::unordered_map<std::string, CSMP_FEM_TYPE> typeMap {
        {"UNKNOWN", UNKNOWN},
        {"LINEAR_RECTANGLE", LINEAR_RECTANGLE},
        {"LINEAR_CUBOID", LINEAR_CUBOID},
        {"LINEAR_BAR", LINEAR_BAR},
        {"QUADRATIC_BAR", QUADRATIC_BAR},
        {"CUBIC_BAR", CUBIC_BAR},
        {"LINEAR_TRIANGLE", LINEAR_TRIANGLE},
        {"LINEAR_TRIANGLE3D", LINEAR_TRIANGLE3D},
        {"BARYCENTRIC_LINEAR_TRIANGLE", BARYCENTRIC_LINEAR_TRIANGLE},
        {"QUADRATIC_TRIANGLE", QUADRATIC_TRIANGLE},
        {"BARYCENTRIC_QUADRATIC_TRIANGLE", BARYCENTRIC_QUADRATIC_TRIANGLE},
        {"CUBIC_TRIANGLE", CUBIC_TRIANGLE},
        {"LINEAR_TETRAHEDRON", LINEAR_TETRAHEDRON},
        {"QUADRATIC_TETRAHEDRON", QUADRATIC_TETRAHEDRON},
        {"BARYCENTRIC_QUADRATIC_TETRAHEDRON", BARYCENTRIC_QUADRATIC_TETRAHEDRON},
        {"CUBIC_TETRAHEDRON", CUBIC_TETRAHEDRON},
        {"ISOPARAMETRIC_LINEAR_BAR", ISOPARAMETRIC_LINEAR_BAR},
        {"ISOPARAMETRIC_QUADRATIC_BAR", ISOPARAMETRIC_QUADRATIC_BAR},
        {"ISOPARAMETRIC_CUBIC_BAR", ISOPARAMETRIC_CUBIC_BAR},
        {"ISOPARAMETRIC_LINEAR_TRIANGLE", ISOPARAMETRIC_LINEAR_TRIANGLE},
        {"ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE", ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE},
        {"ISOPARAMETRIC_QUADRATIC_TRIANGLE", ISOPARAMETRIC_QUADRATIC_TRIANGLE},
        {"ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE", ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE},
        {"ISOPARAMETRIC_CUBIC_TRIANGLE", ISOPARAMETRIC_CUBIC_TRIANGLE},
        {"ISOPARAMETRIC_LINEAR_TETRAHEDRON", ISOPARAMETRIC_LINEAR_TETRAHEDRON},
        {"ISOPARAMETRIC_QUADRATIC_TETRAHEDRON", ISOPARAMETRIC_QUADRATIC_TETRAHEDRON},
        {"ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON", ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON},
        {"ISOPARAMETRIC_CUBIC_TETRAHEDRON", ISOPARAMETRIC_CUBIC_TETRAHEDRON},
        {"ISOPARAMETRIC_LINEAR_PYRAMID", ISOPARAMETRIC_LINEAR_PYRAMID},
        {"ISOPARAMETRIC_QUADRATIC_PYRAMID13", ISOPARAMETRIC_QUADRATIC_PYRAMID13},
        {"ISOPARAMETRIC_QUADRATIC_PYRAMID14", ISOPARAMETRIC_QUADRATIC_PYRAMID14},
        {"ISOPARAMETRIC_CUBIC_PYRAMID", ISOPARAMETRIC_CUBIC_PYRAMID},
        {"ISOPARAMETRIC_LINEAR_PRISM", ISOPARAMETRIC_LINEAR_PRISM},
        {"ISOPARAMETRIC_QUADRATIC_PRISM15", ISOPARAMETRIC_QUADRATIC_PRISM15},
        {"ISOPARAMETRIC_QUADRATIC_PRISM18", ISOPARAMETRIC_QUADRATIC_PRISM18},
        {"ISOPARAMETRIC_CUBIC_PRISM", ISOPARAMETRIC_CUBIC_PRISM},
        {"ISOPARAMETRIC_LINEAR_QUADRILATERAL", ISOPARAMETRIC_LINEAR_QUADRILATERAL},
        {"ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL", ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL},
        {"ISOPARAMETRIC_QUADRATIC_QUADRILATERAL", ISOPARAMETRIC_QUADRATIC_QUADRILATERAL},
        {"ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL", ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL},
        {"ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9", ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9},
        {"ISOPARAMETRIC_CUBIC_QUADRILATERAL", ISOPARAMETRIC_CUBIC_QUADRILATERAL},
        {"ISOPARAMETRIC_LINEAR_HEXAHEDRON", ISOPARAMETRIC_LINEAR_HEXAHEDRON},
        {"ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20", ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20},
        {"ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27", ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27},
        {"ISOPARAMETRIC_CUBIC_HEXAHEDRON", ISOPARAMETRIC_CUBIC_HEXAHEDRON},
        {"ZERO_DIMENSIONAL_FACE", ZERO_DIMENSIONAL_FACE},
        {"POINT_ELEMENT", POINT_ELEMENT},
        {"POLYGONAL_ELEMENT", POLYGONAL_ELEMENT},
        {"POLYHEDRAL_ELEMENT", POLYHEDRAL_ELEMENT},
        {"EXPERIMENTAL_ELEMENT", EXPERIMENTAL_ELEMENT}
    };

    auto it = typeMap.find(etype);
    if (it != typeMap.end()) return it->second;

    std::cerr << "\nparseFiniteElementType: Could not identify element type: " << etype << std::endl;
    return UNKNOWN;

 } // end parseFiniteElementType

 
 
 
 
 
 /// converts enum names into text
 const char* parseFiniteElementType( int8_t etype ) noexcept
  {
    switch (etype)
    {
        case UNKNOWN: return "UNKNOWN";
        case LINEAR_RECTANGLE: return "LINEAR_RECTANGLE";
        case LINEAR_CUBOID: return "LINEAR_CUBOID";

        case LINEAR_BAR: return "LINEAR_BAR";
        case QUADRATIC_BAR: return "QUADRATIC_BAR";
        case CUBIC_BAR: return "CUBIC_BAR";

        case LINEAR_TRIANGLE: return "LINEAR_TRIANGLE";
        case LINEAR_TRIANGLE3D: return "LINEAR_TRIANGLE3D";
        case BARYCENTRIC_LINEAR_TRIANGLE: return "BARYCENTRIC_LINEAR_TRIANGLE";
        case QUADRATIC_TRIANGLE: return "QUADRATIC_TRIANGLE";
        case BARYCENTRIC_QUADRATIC_TRIANGLE: return "BARYCENTRIC_QUADRATIC_TRIANGLE";
        case CUBIC_TRIANGLE: return "CUBIC_TRIANGLE";

        case LINEAR_TETRAHEDRON: return "LINEAR_TETRAHEDRON";
        case QUADRATIC_TETRAHEDRON: return "QUADRATIC_TETRAHEDRON";
        case BARYCENTRIC_QUADRATIC_TETRAHEDRON: return "BARYCENTRIC_QUADRATIC_TETRAHEDRON";
        case CUBIC_TETRAHEDRON: return "CUBIC_TETRAHEDRON";

        case ISOPARAMETRIC_LINEAR_BAR: return "ISOPARAMETRIC_LINEAR_BAR";
        case ISOPARAMETRIC_QUADRATIC_BAR: return "ISOPARAMETRIC_QUADRATIC_BAR";
        case ISOPARAMETRIC_CUBIC_BAR: return "ISOPARAMETRIC_CUBIC_BAR";

        case ISOPARAMETRIC_LINEAR_TRIANGLE: return "ISOPARAMETRIC_LINEAR_TRIANGLE";
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE: return "ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE";
        case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return "ISOPARAMETRIC_QUADRATIC_TRIANGLE";
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE: return "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE";
        case ISOPARAMETRIC_CUBIC_TRIANGLE: return "ISOPARAMETRIC_CUBIC_TRIANGLE";

        case ISOPARAMETRIC_LINEAR_TETRAHEDRON: return "ISOPARAMETRIC_LINEAR_TETRAHEDRON";
        case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: return "ISOPARAMETRIC_QUADRATIC_TETRAHEDRON";
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON: return "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON";
        case ISOPARAMETRIC_CUBIC_TETRAHEDRON: return "ISOPARAMETRIC_CUBIC_TETRAHEDRON";

        case ISOPARAMETRIC_LINEAR_PYRAMID: return "ISOPARAMETRIC_LINEAR_PYRAMID";
        case ISOPARAMETRIC_QUADRATIC_PYRAMID13: return "ISOPARAMETRIC_QUADRATIC_PYRAMID13";
        case ISOPARAMETRIC_QUADRATIC_PYRAMID14: return "ISOPARAMETRIC_QUADRATIC_PYRAMID14";
        case ISOPARAMETRIC_CUBIC_PYRAMID: return "ISOPARAMETRIC_CUBIC_PYRAMID";

        case ISOPARAMETRIC_LINEAR_PRISM: return "ISOPARAMETRIC_LINEAR_PRISM";
        case ISOPARAMETRIC_QUADRATIC_PRISM15: return "ISOPARAMETRIC_QUADRATIC_PRISM15";
        case ISOPARAMETRIC_QUADRATIC_PRISM18: return "ISOPARAMETRIC_QUADRATIC_PRISM18";
        case ISOPARAMETRIC_CUBIC_PRISM: return "ISOPARAMETRIC_CUBIC_PRISM";

        case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return "ISOPARAMETRIC_LINEAR_QUADRILATERAL";
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL: return "ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL";
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL";
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL: return "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL";
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9: return "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9";
        case ISOPARAMETRIC_CUBIC_QUADRILATERAL: return "ISOPARAMETRIC_CUBIC_QUADRILATERAL";

        case ISOPARAMETRIC_LINEAR_HEXAHEDRON: return "ISOPARAMETRIC_LINEAR_HEXAHEDRON";
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20: return "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20";
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: return "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27";
        case ISOPARAMETRIC_CUBIC_HEXAHEDRON: return "ISOPARAMETRIC_CUBIC_HEXAHEDRON";

        case ZERO_DIMENSIONAL_FACE: return "ZERO_DIMENSIONAL_FACE";
        case POINT_ELEMENT: return "POINT_ELEMENT";
        case POLYGONAL_ELEMENT: return "POLYGONAL_ELEMENT";
        case POLYHEDRAL_ELEMENT: return "POLYHEDRAL_ELEMENT";

        case EXPERIMENTAL_ELEMENT: return "EXPERIMENTAL_ELEMENT";

        default:
            std::cerr << "\nparseFiniteElementType: Could not identify element type: " 
                      << static_cast<int>(etype) << std::endl;
            return "UNKNOWN";
    }
 
 } // end parseFiniteElementType



 /// converts enum names into text
 const char* parseAbbreviated_FE_Type( int8_t etype ) noexcept
  {
    switch (etype)
    {
        case UNKNOWN: return "?";

        case LINEAR_RECTANGLE: return "RECT";
        case LINEAR_CUBOID: return "CUBOID";

        case LINEAR_BAR: return "BAR";
        case QUADRATIC_BAR: return "BAR^2";
        case CUBIC_BAR: return "BAR^3";

        case LINEAR_TRIANGLE: return "TRIA";
        case LINEAR_TRIANGLE3D: return "TRIA3D";
        case BARYCENTRIC_LINEAR_TRIANGLE: return "BTRIA";
        case QUADRATIC_TRIANGLE: return "TRI^2";
        case BARYCENTRIC_QUADRATIC_TRIANGLE: return "BTRI^2";
        case CUBIC_TRIANGLE: return "TRI^3";

        case LINEAR_TETRAHEDRON: return "TET";
        case QUADRATIC_TETRAHEDRON: return "TET^2";
        case BARYCENTRIC_QUADRATIC_TETRAHEDRON: return "BTET^2";
        case CUBIC_TETRAHEDRON: return "TET^3";

        case ISOPARAMETRIC_LINEAR_BAR: return "IBAR";
        case ISOPARAMETRIC_QUADRATIC_BAR: return "IBAR^2";
        case ISOPARAMETRIC_CUBIC_BAR: return "IBAR^3";

        case ISOPARAMETRIC_LINEAR_TRIANGLE: return "ITRIA";
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE: return "IBTRIA";
        case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return "ITRIA^2";
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE: return "IBTRIA^2";
        case ISOPARAMETRIC_CUBIC_TRIANGLE: return "ITRIA^3";

        case ISOPARAMETRIC_LINEAR_TETRAHEDRON: return "ITET";
        case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: return "ITET^2";
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON: return "IBTET^2";
        case ISOPARAMETRIC_CUBIC_TETRAHEDRON: return "ITET^3";

        case ISOPARAMETRIC_LINEAR_PYRAMID: return "IPYR";
        case ISOPARAMETRIC_QUADRATIC_PYRAMID13: return "IPYR13";
        case ISOPARAMETRIC_QUADRATIC_PYRAMID14: return "IPYR14";
        case ISOPARAMETRIC_CUBIC_PYRAMID: return "IPYR^3";

        case ISOPARAMETRIC_LINEAR_PRISM: return "IPRISM";
        case ISOPARAMETRIC_QUADRATIC_PRISM15: return "IPRISM15";
        case ISOPARAMETRIC_QUADRATIC_PRISM18: return "IPRISM18";
        case ISOPARAMETRIC_CUBIC_PRISM: return "IPRISM^3";

        case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return "IQUAD";
        case ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL: return "IBQUAD";
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return "IQUAD^2";
        case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL: return "IBQUAD^2";
        case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9: return "IQUAD9^2";
        case ISOPARAMETRIC_CUBIC_QUADRILATERAL: return "IQUAD^3";

        case ISOPARAMETRIC_LINEAR_HEXAHEDRON: return "IHEX";
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20: return "IHEX20";
        case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: return "IHEX27";
        case ISOPARAMETRIC_CUBIC_HEXAHEDRON: return "IHEX^3";

        case ZERO_DIMENSIONAL_FACE: return "0D_FACE";
        case POINT_ELEMENT: return "POINT";
        case POLYGONAL_ELEMENT: return "POLYGON";
        case POLYHEDRAL_ELEMENT: return "POLYHED";

        case EXPERIMENTAL_ELEMENT: return "EXPERIMENTAL";

        default:
            std::cerr << "\nparseAbbreviated_FE_Type: Could not identify element type: "
                      << static_cast<int>(etype) << std::endl;
            return "?";
    }
 
 } // end parseAbbreviated_FE_Type



// hashing function for strings
constexpr unsigned int hash(const char* str, int h = 0) {
    return !str[h] ? 5381 : (hash(str, h + 1) * 33) ^ static_cast<unsigned char>(str[h]);
}

FV_FACET_TYPE parseFacetType(const std::string& ftype) noexcept
{
    switch(hash(ftype.c_str()))
    {
        case hash("POINT_FACET"): return POINT_FACET;
        case hash("UNIT_LINEAR_FACET"): return UNIT_LINEAR_FACET;
        case hash("TRIANGULAR_FACET"): return TRIANGULAR_FACET;
        case hash("QUADRILATERAL_FACET"): return QUADRILATERAL_FACET;
        default:
            std::cerr << "\nparseFacetType: Could not identify facet type " << ftype << '\n';
            return static_cast<FV_FACET_TYPE>(0);
    }
}


const char* parseFacetType( int8_t ftype ) noexcept
 {
     switch (ftype)
     {
         case POINT_FACET: return "POINT_FACET";
         case UNIT_LINEAR_FACET: return "UNIT_LINEAR_FACET";
         case TRIANGULAR_FACET: return "TRIANGULAR_FACET";
         case QUADRILATERAL_FACET: return "QUADRILATERAL_FACET";
         default:
             cerr << "\nparseFacetType: Could not identify facet type " << ftype << '\n';
             return "(unknown)";
     }
 } // end parseFacetType


/*
bool isTriangularElement( CSMP_FEM_TYPE );
bool isQuadrilateralElement( CSMP_FEM_TYPE );
bool isLineElement( CSMP_FEM_TYPE );
bool isTriangular( CSMP_FEM_TYPE );
bool isQuadrilateral( CSMP_FEM_TYPE );
bool isHexahedral( CSMP_FEM_TYPE );
*/
CSMP_FEM_TYPE finiteElementTypeOfSharedFace(CSMP_FEM_TYPE etype1, CSMP_FEM_TYPE etype2, bool isoparametric) noexcept
{
    if (!isoparametric) {
        cerr << "\nfiniteElementTypeOfSharedFace: only isoparametric elements are handled so far.";
        return UNKNOWN;
    }

    // Handle line elements first
    if (isLineElement(etype1) || isLineElement(etype2)) return ZERO_DIMENSIONAL_FACE;

    // Determine the shared face type
    switch (etype1) {
        case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
            switch (etype2) {
                case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
                case ISOPARAMETRIC_LINEAR_PRISM:
                case ISOPARAMETRIC_LINEAR_PYRAMID:
                    return ISOPARAMETRIC_LINEAR_TRIANGLE;
                default:
                    if (isTriangular(etype2) || isQuadrilateral(etype2))
                        return ISOPARAMETRIC_LINEAR_BAR;
                    break;
            }
            return ISOPARAMETRIC_LINEAR_TRIANGLE; // fallback for any unmatched case
        case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
            switch (etype2) {
                case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
                case ISOPARAMETRIC_LINEAR_PRISM:
                case ISOPARAMETRIC_LINEAR_PYRAMID:
                    return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
                default:
                    if (isTriangular(etype2) || isQuadrilateral(etype2))
                        return ISOPARAMETRIC_LINEAR_BAR;
                    break;
            }
            return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
        case ISOPARAMETRIC_LINEAR_PRISM:
            switch (etype2) {
                case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
                    return ISOPARAMETRIC_LINEAR_TRIANGLE;
                case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
                    return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
                default:
                    if (isTriangular(etype2) || isQuadrilateral(etype2))
                        return ISOPARAMETRIC_LINEAR_BAR;
                    break;
            }
            break;
        case ISOPARAMETRIC_LINEAR_PYRAMID:
            switch (etype2) {
                case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
                    return ISOPARAMETRIC_LINEAR_TRIANGLE;
                case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
                    return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
                default:
                    if (isTriangular(etype2) || isQuadrilateral(etype2))
                        return ISOPARAMETRIC_LINEAR_BAR;
                    break;
            }
            break;
        default:
            if (isTriangular(etype1) || isQuadrilateral(etype1) ||
                isTriangular(etype2) || isQuadrilateral(etype2))
                return ISOPARAMETRIC_LINEAR_BAR;
            break;
    }

    return UNKNOWN;
}



/**
 * @brief Determines (for a Simplex element mesh consisting of straight line, triangle, and tetrahedral elements) the corresponding quadratic CSMP_FEM_TYPE for any given linear type.
 * @param isoparametric The VData mesh to inspect.
 * @return The quadratic counterpart if found, otherwise UNKNOWN.
 */
CSMP_FEM_TYPE getQuadraticType( bool isoparametric, int8_t etype ) noexcept
  {
    CSMP_FEM_TYPE qetype = UNKNOWN;

    if (isoparametric) {
        switch (etype) {
            case ISOPARAMETRIC_LINEAR_BAR:         qetype = ISOPARAMETRIC_QUADRATIC_BAR;          break;
            case ISOPARAMETRIC_LINEAR_TRIANGLE:    qetype = ISOPARAMETRIC_QUADRATIC_TRIANGLE;     break;
            case ISOPARAMETRIC_LINEAR_TETRAHEDRON: qetype = ISOPARAMETRIC_QUADRATIC_TETRAHEDRON;  break;
            case ISOPARAMETRIC_LINEAR_HEXAHEDRON:  qetype = ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20; break;
            case ISOPARAMETRIC_LINEAR_PRISM:       qetype = ISOPARAMETRIC_QUADRATIC_PRISM15;      break;
            case ISOPARAMETRIC_LINEAR_PYRAMID:     qetype = ISOPARAMETRIC_QUADRATIC_PYRAMID13;    break;
            default:
                std::cerr << "\nGetQuadraticType: Unsupported Isoparametric element '" 
                          << parseFiniteElementType(etype) << "' probably not a simplex element." << std::endl;
        }
    } else {
        // Analytically integrated elements
        switch (etype) {
            case LINEAR_BAR:                       qetype = QUADRATIC_BAR;           break;
            case LINEAR_TRIANGLE:
            case LINEAR_TRIANGLE3D:                qetype = QUADRATIC_TRIANGLE;      break;
            case LINEAR_TETRAHEDRON:               qetype = QUADRATIC_TETRAHEDRON;   break;
            default:
                std::cerr << "\nGetQuadraticType: Unsupported Analytic element '" 
                          << parseFiniteElementType(etype) << "' probably not a simplex element." << std::endl;
        }
    }

    return qetype;
}

 
 } // end namespace csmp
 
 
 
 
