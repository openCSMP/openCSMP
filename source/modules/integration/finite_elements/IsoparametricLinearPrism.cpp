#include "IsoparametricLinearPrism.h"
#include "Exception.h"
#include "TriangularFacet.h"
#include "QuadrilateralFacet.h"

using namespace std;

namespace csmp {

IsoparametricLinearPrism::IsoparametricLinearPrism( size_t integrationPoints )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_PRISM, true, true, 1U ),
      NXYZ(6,3),
      IP(integrationPoints,3)
 {
    //AAM, 12.07.02
    dim = 3;
    itp = 1;
    // npf= 3/4?
    npf = 4;
    npe = 6;
    fpe = 5;
    spe = 9;
    epe = 5;
    // nne - typical number of finite elements, sharing each node
    // 8?
    nne = 8;
    cne = 0;
    // Integration points
    gpe = integrationPoints;

    M.Resize(npe,npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_LINEAR_PRISM);

    // base class matrices
    XY.Resize(npe,dim);
    JAC.Resize(dim,dim);
    JINV.Resize(dim,dim);

    // base class vectors
    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    /* Element nymbering scheme ...

          ^ t
          |
          |
        4 o

          |          o 6
          |
        5
    o     |
          |
        1 o -------- o 3  --> s
         /
        /
       /
      /
     /
    o 2
    r

    Note that this numbering scheme is similar to used by ANSYSCFD
    */

    // initializing local node coordinates
    // local node coordinates are defined as r==ksi, s==nu, t==mu.
    // similar to the base of the tetrahedron

    NXYZ(0,0) = 0.0, NXYZ(0,1) = 0.0, NXYZ(0,2) = -1.0;
    NXYZ(1,0) = 1.0, NXYZ(1,1) = 0.0, NXYZ(1,2) = -1.0;
    NXYZ(2,0) = 0.0, NXYZ(2,1) = 1.0, NXYZ(2,2) = -1.0;
    NXYZ(3,0) = 0.0, NXYZ(3,1) = 0.0, NXYZ(3,2) =  1.0;
    NXYZ(4,0) = 1.0, NXYZ(4,1) = 0.0, NXYZ(4,2) =  1.0;
    NXYZ(5,0) = 0.0, NXYZ(5,1) = 1.0, NXYZ(5,2) =  1.0;

    W.resize( gpe );

    //Numeric integration 2nd order, number of point m=6
    if(integrationPoints==1)
     {
         IP(0,0)=1.0/3.0; IP(0,1)=1.0/3.0; IP(0,2)=0.0;
         W[0]=1.0;
     }
     else if (integrationPoints==6 )
         // Generate 6 integration points for full numeric integration + weights
         GenerateIntegrationPoints(IP,W);
     else
     {
         cout<<"IsoparametricLinearPrism::IsoparametricLinearHexahedron: Number of integration points should be 1 or 6"<<endl;
         throw std::range_error("***ERROR: IsoparametricLinearPrism::IsoparametricLinearHexahedron: N of integration points should be 1 or 6");
     }
 }


IsoparametricLinearPrism::~IsoparametricLinearPrism()
 {
 }




/**

  void IsoparametricLinearPrism::Nrst( double64 r, double64 s, vector<double64>& N ) const



Description:

Returns the value of the element interpolation functions at the point 'rst'
in local coordinates.

@section arguments Input Arguments

The floating point coordinates 'r', 's' and 't' of the point at which the
interpolation shall be carried out.

@param N The fourth method argument is the vector into which the values of the
n=nodes interpolation functions at the point 'rst' will be returned.

@section implementation Implementation

See in source code header file.

@section application Application

Method is used to compute property values at the integration points of
the element.
*/
void IsoparametricLinearPrism::Nrst(
                    double64 r,
                    double64 s,
                    double64 t,
                    std::vector<double64>& N ) const
{
   N.resize(npe);

   const double64 tPlus   = 1.0+t;
   const double64 tMinus  = 1.0-t;
   const double64 L1      = 1.0-r-s;
   const double64 L2      = r;
   const double64 L3      = s;

   //cout<<" IsoparametricLinearPrism::Nrst L1="<<L1<<", L2="<<L2<<", L3="<<L3<<" Sum="<<L1+L2+L3<<endl;

   N[0] = 0.5*L1*tMinus;
   N[1] = 0.5*L2*tMinus;
   N[2] = 0.5*L3*tMinus;
   N[3] = 0.5*L1*tPlus;
   N[4] = 0.5*L2*tPlus;
   N[5] = 0.5*L3*tPlus;

 }

/**

This method and the complementary method dNs() compute the shape function
derivates with respect to the local coordinate axis 'r', and 's'.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

@param DNR The shape function derivatives are returned into the second method
argument which is a floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.
*/
void IsoparametricLinearPrism::dNr (
                double64,
                double64,
                double64 t,
                std::vector<double64>& DNR ) const
{
   double64 tPlus   = 1.0+t;
   double64 tMinus  = 1.0-t;

   DNR.resize(npe);

   DNR[0] = -0.5*tMinus;
   DNR[1] =  0.5*tMinus;
   DNR[2] =  0.0;
   DNR[3] = -0.5*tPlus;
   DNR[4] =  0.5*tPlus;
   DNR[5] =  0.0;

}

void IsoparametricLinearPrism::dNs(
                double64,
                double64,
                double64 t,
                std::vector<double64>& DNS ) const
{
   double64 tPlus   = 1.0+t;
   double64 tMinus  = 1.0-t;

   DNS.resize(npe);

   DNS[0] = -0.5*tMinus;
   DNS[1] =  0.0;
   DNS[2] =  0.5*tMinus;
   DNS[3] = -0.5*tPlus;
   DNS[4] =  0.0;
   DNS[5] =  0.5*tPlus;

}


void IsoparametricLinearPrism::dNt(
                double64 r,
                double64 s,
                double64,
                std::vector<double64>& DNT ) const
{

   DNT.resize(npe);

   const double64 rMinus  =   1.0-r;

   DNT[0] = -0.5*(rMinus-s);
   DNT[1] = -0.5*r;
   DNT[2] = -0.5*s;
   DNT[3] =  0.5*(rMinus-s);
   DNT[4] =  0.5*r;
   DNT[5] =  0.5*s;

}



//tested AAM  ok1
void
 IsoparametricLinearPrism::JacobianAtIntegrationPoint( size_t gauss_point )
 {
    assert( gauss_point < gpe );
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }




void
IsoparametricLinearPrism::N_AtBaryCenter( std::vector<double64>& N )
 {
    N.resize(npe);

    double64 OneThird=1.0/3.0;
    Nrst(OneThird,OneThird,0.0,N);

 }


/**

Computes the interpolation function values at the specified integration
point 'ip'.

@param N The interpolation function values are returned into the third argument.
*/
void
IsoparametricLinearPrism::N_AtIntegrationPoint( size_t ip, std::vector<double64>& N )
 {
    assert( ip < gpe );
     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), N );
 } // end N_AtIntegrationPoint







/**

Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local corner node ID numbers
for the element.

*/
void
IsoparametricLinearPrism::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
    ids[5] = 5;
}




/** Returns the counter-clockwise local node numbering for the element.

@param ids Returns the counter-clockwise local node numbering for the element.
*/
void
IsoparametricLinearPrism::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
    ids[5] = 5;
}


void
IsoparametricLinearPrism::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    snids.resize(2);
     if ( segm_id == 0 ) {
         snids[0] = 0;
         snids[1] = 1;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 1;
         snids[1] = 2;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 2;
         snids[1] = 0;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 0;
         snids[1] = 3;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 1;
         snids[1] = 4;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 2;
         snids[1] = 5;
      }
    else if ( segm_id == 6 ) {
         snids[0] = 3;
         snids[1] = 4;
      }
    else if ( segm_id == 7 ) {
         snids[0] = 4;
         snids[1] = 5;
      }
    else if ( segm_id == 8 ) {
         snids[0] = 5;
         snids[1] = 3;
      }
    else
    std::cout <<"\nIsoparametricLinearPrism::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment




/**
     Starting with the smallest local node number, the face nodes are returned
     in counter clockwise order from the outside looking in.
*/
void IsoparametricLinearPrism::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {
    // bottom face
    if( face_id == 0  )
      {
         fnids.resize(3);
         fnids[0] = 0;
         fnids[1] = 2;
         fnids[2] = 1;
          }
    // first side face (in counter-clockwise order)
    else if ( face_id == 1 )
      {
         fnids.resize(4);
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 4;
         fnids[3] = 3;
      }
    // second side face
    else if ( face_id == 2 )
      {
         fnids.resize(4);
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 5;
         fnids[3] = 4;
      }
    // third side face (in the back)
    else if ( face_id == 3 )
      {
         fnids.resize(4);
         fnids[0] = 0;
         fnids[1] = 3;
         fnids[2] = 5;
         fnids[3] = 2;
      }
    // top face
    else if ( face_id == 4 )
      {
         fnids.resize(3);
         fnids[0] = 3;
         fnids[1] = 4;
         fnids[2] = 5;
      }
    else
    std::cerr <<"\nIsoparametricLinearPrism::NodesOfFace: Invalid Face ID requested: "<< face_id << std::endl;
 }


CSMP_FEM_TYPE  IsoparametricLinearPrism::ElementTypeOfFace( size_t face )  const
 {
    assert( face < fpe );
    if ( face == 0U || face == 4U )
        return ISOPARAMETRIC_LINEAR_TRIANGLE;

    return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
 }

double64 IsoparametricLinearPrism::WeightAtIntegrationPoint( size_t i ) const { return W[i]; }




void
IsoparametricLinearPrism::GenerateIntegrationPoints(
                                                    DenseMatrix<DM_MIN> &Ip,
                                                    std::vector<double64>& We)
{
    const size_t numberOfTriaIntegrationPoints=3; // order 2
    const size_t numberOfDimensionsInPlane=2;
    DenseMatrix<DM_MIN> IPTRIA(numberOfTriaIntegrationPoints,numberOfDimensionsInPlane);
    const double64 constA=0.577350269189626, constT1=0.66666666666667, constT2=0.16666666666667;

    IPTRIA(0,0)=constT1; IPTRIA(0,1)=constT2;
    IPTRIA(1,0)=constT2; IPTRIA(1,1)=constT1;
    IPTRIA(2,0)=constT2; IPTRIA(2,1)=constT2;

    vector<double64> WTRIA(numberOfTriaIntegrationPoints);
    // Weights for the base integration
    WTRIA[0]=constT2;
    WTRIA[1]=constT2;
    WTRIA[2]=constT2;

    const size_t numberOfLineIntegrationPoints=2;
    vector<double64> IPLINE(numberOfLineIntegrationPoints),WLINE(numberOfLineIntegrationPoints);

    IPLINE[0]=-constA;  IPLINE[1]=constA;
    WLINE [0]= 1.0;     WLINE [1]=1.0;

    size_t i=0;
    for(size_t j=0;j<numberOfTriaIntegrationPoints;j++)
    {
    for(size_t k=0;k<numberOfLineIntegrationPoints;k++)
        {
            Ip(i,0)=IPTRIA(j,0);
            Ip(i,1)=IPTRIA(j,1);
            Ip(i,2)=IPLINE[k];
            We[i++]=WTRIA[j]*WLINE[k];
        }
    }
}

double64
IsoparametricLinearPrism::AreaOfBase(
                            std::vector<double64> &V1XYZ,
                            std::vector<double64> &V2XYZ,
                            std::vector<double64> &V3XYZ
                            )
{

    double64 Det1=V1XYZ[1]*(V2XYZ[2]-V3XYZ[2])-V2XYZ[1]*(V1XYZ[2]-V3XYZ[2])+V3XYZ[1]*(V1XYZ[2]-V2XYZ[2]);
    double64 Det2=V1XYZ[2]*(V2XYZ[0]-V3XYZ[0])-V2XYZ[2]*(V1XYZ[0]-V3XYZ[0])+V3XYZ[2]*(V1XYZ[0]-V2XYZ[0]);
    double64 Det3=V1XYZ[0]*(V2XYZ[1]-V3XYZ[1])-V2XYZ[0]*(V1XYZ[1]-V3XYZ[1])+V3XYZ[0]*(V1XYZ[1]-V2XYZ[1]);


    return 0.5*sqrt	(Det1*Det1+Det2*Det2+Det3*Det3);

}



/**

Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN6 The interpolation-function derivative matrix is returned into the
second method argument.

*/
void
IsoparametricLinearPrism::dN( DenseMatrix<DM_MIN>& DN6 )
{
    DN6.Resize(dim,npe);
    M.Resize(dim,1);

     // Jacobian transformation to global coordinate system
     for ( size_t i=0; i<npe; i++ )
       {
          // here the global coordinates come in
          dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
          dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
          dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

          // compute Jacobian matrix, its determinant and inversex
          Jacobian( DNR, DNS, DNT );
          JacobianInverse();
          ////////// Debug Printout////////////////////////////////
          //  cout<<" At vertice #: "<<i<<endl;
          //	cout<<" Jacobian Matrix: "<<endl;
          //	JAC.Out(cout);
          //	cout<<" Jacobian Inverse Matrix: "<<endl;
          //JINV.Out(cout);
          /////////////////////////////////////////////////////////
          //for ( size_t j=0; j<dim; j++ ) TEMP(j,0) = DN8(j,i);
          M(0,0)=DNR[i]; M(1,0)=DNS[i]; M(2,0)=DNT[i];

          // 3x3 * 3x1 = 3x1 gives the global DN entries
          JINV *= M;
          DN6(0,i) = JINV(0,0);
          DN6(1,i) = JINV(1,0);
          DN6(2,i) = JINV(2,0);
          JINV.Resize(dim,dim);
       }
} // end dN




/**

Computes derivatives of shape functions at a point XY given in global coordinates
within the element (also given in global coordinates).

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at the point XY are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument.
*/
double64
IsoparametricLinearPrism::dN( DenseMatrix<DM_MIN>& DN2,
                              const vector<double64>& xyz  )
  {

    vector<double64> rst(dim);

    PhysicalToParametric(rst, xyz);

    // here find derivatives
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    DN2.Resize(dim,dim);
    DN2  = JINV;

    DenseMatrix<DM_MIN> DN(dim,npe);

    dN(DN);

    DN2*=DN;

    /////////////////////////////// Debug printout ///////////////////////////////////////////////
    cout<<" IsoparametricLinearPrism::dN  For given xyz=("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<"), rst=("<<
    rst[0]<<","<<rst[1]<<","<<rst[2]<<")"<<endl;
    cout<<" IsoparametricLinearPrism::dN  Matrix DN2: "<<endl;
    DN2.Out(cout);
    getchar();
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

    return detJ;

 }

/** Outputs the shape functions at the point 'xyz' which must lie
within the Hexahedron.
*/
void
IsoparametricLinearPrism::N(
                    std::vector<double64>& N,
                    const std::vector<double64>& xyz
                    )
{
    vector<double64> rst(dim);

    PhysicalToParametric(rst, xyz);

    Nrst(rst[0],rst[1],rst[2], N);
}




/**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J at the desired node
point.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the
point i (i=0...i=gauss points-1).

@section arguments Input Arguments

The first argument is a reference to the parent Element object,
the second argument is the result matrix, the third argument indicates
the gauss point at which the derivative matrix is computed and the
fourth argument specifies the degrees of freedom per node, for which
the B matrix shall be transformed.

@return The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.
*/
double64
IsoparametricLinearPrism::dN_AtNode( DenseMatrix<DM_MIN>& B, size_t nd )
 {
    assert( nd < npe );
    dNr( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR );
    dNs( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS );
    dNt( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 6 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];
    B(0,4) = DNR[4], B(1,4) = DNS[4], B(2,4) = DNT[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5], B(2,5) = DNT[5];

    B = JINV * B;

    return detJ;
 }



/// Projection function from rst->xyz
void
IsoparametricLinearPrism::ParametricToPhysical(std::vector<double64> &rst, std::vector<double64>& xyz)
{
    Nrst(rst[0],rst[1],rst[2], NRST );

    for(size_t i=0; i<dim; i++) xyz[i]=0.0;

    for(size_t i=0; i<npe; i++)
    {
        xyz[0]+=XY(i,0)*NRST[i];
        xyz[1]+=XY(i,1)*NRST[i];
        xyz[2]+=XY(i,2)*NRST[i];
    }
}



/** Projection function from xyz->rst
 Starts from the 2/3 2/3 0 point in  the parametric space
 */
void
IsoparametricLinearPrism::PhysicalToParametric(
                                       std::vector<double64>& rSt,
                                       const std::vector<double64>& xyz
                                       )
{

    vector<double64> outxyz(dim);
    vector<double64> rstHatK(dim);

    std::vector<double64> distanceFromGivenPointLinf(dim,0.0);
    double64 distanceFromGivenPointL2;

    // Find largest and smallest segments in order to define precision
    vector<double64> vec(spe);
    EdgeLengths( vec );
    double64 seg_max(vec[0]), seg_min(vec[0]);
    for ( size_t i=1; i<spe; i++ )
    {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
    }

    const double64 geometricTolerance = 0.005*seg_min;

    // First guess as BaryCenter
    rstHatK[0] = 1.0/3.0;
    rstHatK[1] = 1.0/3.0;
    rstHatK[2] = 0.0;

    ParametricToPhysical( rstHatK, outxyz);

    distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
    distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
    distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );

    distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                     distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                     distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

    if( distanceFromGivenPointL2 > geometricTolerance )
    //if( (distanceFromGivenPointLinf[0] > geometricTolerance) ||
    //    (distanceFromGivenPointLinf[1] > geometricTolerance) ||
    //    (distanceFromGivenPointLinf[2] > geometricTolerance)  )
    {

        vector<double64> rstHatK_PlusOne(dim);
        double64 minDistanceFromGivenPoint;

        const double64 constantMu               = 1.0;
        const size_t numberOfFirstIterrations   = 5;
        const size_t maxNumberOfIterrations     = 20;
        const size_t numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double64 incrementR = 1.0/(numberOfFirstIterrations-1);
        const double64 incrementS = 1.0/(numberOfFirstIterrations-1);
        const double64 incrementT = 2.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = 0.0;
        for(size_t i=0;i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = 0.0;
            for(size_t j=0;j<numberOfFirstIterrations;j++)
            {
                rstHatK_PlusOne[2] = -1.0;
                for(size_t k=0;k<numberOfFirstIterrations;k++)
                {
                    ParametricToPhysical( rstHatK_PlusOne, outxyz);

                    distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
                    distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
                    distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );


                    distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                                     distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                                     distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

                    if( minDistanceFromGivenPoint > distanceFromGivenPointL2 )
                    {
                        for(size_t i=0; i<dim; i++)
                            rstHatK[i] = rstHatK_PlusOne[i];

                        minDistanceFromGivenPoint = distanceFromGivenPointL2;
                    }
                    rstHatK_PlusOne[2] += incrementT;
                }
                rstHatK_PlusOne[1] += incrementS;
            }
            rstHatK_PlusOne[0] += incrementR;
        }

        ParametricToPhysical( rstHatK, outxyz);

        distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
        distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
        distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );

        distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                         distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                         distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

        size_t iteration = 1;

        /// Newton-Raphson iterations
        while ( ( distanceFromGivenPointL2 > geometricTolerance ) && ( iteration < maxNumberOfIterrations ) )
        //while ( ( ( distanceFromGivenPointLinf[0] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[1] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[2] > geometricTolerance )    )
        //      &&  ( iteration < maxNumberOfIterrations )                )
        {
            // Out of range check
            if( (rstHatK[0]< 0.) || (rstHatK[0]>1.) ||
                (rstHatK[1]< 0.) || (rstHatK[1]>1.) ||
                (rstHatK[2]<-1.) || (rstHatK[2]>1.)  )
            {
                rstHatK[2] = 0.0;
                rstHatK[2] = 0.0;
                rstHatK[2] = 2.0;

                break;
            }

            // Inverse Jacobian calculations
            if( iteration < numberOfIterationsWhenJacobiIsNotConstant )
            {
                dNr( rstHatK[0], rstHatK[1], rstHatK[2], DNR );
                dNs( rstHatK[0], rstHatK[1], rstHatK[2], DNS );
                dNt( rstHatK[0], rstHatK[1], rstHatK[2], DNT );
                Jacobian( DNR, DNS, DNT );
                // Check whether Jacobian is positive ( might be not true for the point outside the element )
                //const double64 detJ = JacobianDeterminant();
                //if( detJ > 0.0 )
                JacobianInverse();
            }

            // Newton iteration
            rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) +JINV(2,0)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) +JINV(2,1)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[2] = rstHatK[2] - constantMu*(JINV(0,2)*(outxyz[0]-xyz[0]) + JINV(1,2)*(outxyz[1]-xyz[1]) +JINV(2,2)*(outxyz[2]-xyz[2]));

            for(size_t i=0; i<dim; i++)
                rstHatK[i] = rstHatK_PlusOne[i];

            ParametricToPhysical( rstHatK, outxyz);
            distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
            distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
            distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );
            distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                             distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                             distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

            iteration++;

        }

        JAC.Zero();

        if( iteration == maxNumberOfIterrations )
        {
            cout<<" IsoparametricLinearPrism::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            cout<<" IsoparametricLinearPrism::PhysicalToParametric: Real Point:\t"
                <<"x = "<<xyz[0]<<" ;\t"
                <<"y = "<<xyz[1]<<" ;\t"
                <<"z = "<<xyz[2]<<"\n";
            cout<<" IsoparametricLinearPrism::PhysicalToParametric: Found Point:\t"
               <<"x = "<<outxyz[0]<<" ;\t"
               <<"y = "<<outxyz[1]<<" ;\t"
               <<"z = "<<outxyz[2]<<"\n";
            cout<<" IsoparametricLinearPrism::PhysicalToParametric: Distance:\t"
               <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
               <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
               <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            cout<<" IsoparametricLinearPrism::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<" ;\t"
               <<"t = "<<rstHatK[2]<<"\n";

            std::vector<double64> N(npe,0.0);
            Nrst(rstHatK[0], rstHatK[1], rstHatK[2], N );
            cout<<" IsoparametricLinearPrism::PhysicalToParametric: Shape functions:\t"
                <<"N[0] = "<<N[0]<<" ;\t"
                <<"N[1] = "<<N[1]<<" ;\t"
                <<"N[2] = "<<N[2]<<" ;\t"
                <<"N[3] = "<<N[3]<<" ;\t"
                <<"N[4] = "<<N[4]<<" ;\t"
                <<"N[5] = "<<N[5]<<"\n";

            csmp::Exception( CSMP_WARNING, "IsoparametricLinearPrism::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for(size_t i=0; i<dim; i++)
        rSt[i] = rstHatK[i];
}




/**
     Uses the TriangularFacet and QuadrilateralFacet classes to compute the normals to the prism faces.
     Based of the face nodes as specified in NodesOfFace()

     @author SKM 15/2/2016
     
     @test OK
*/
void  IsoparametricLinearPrism::UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const
 {
     assert( face < Faces() );
     unrml.resize(3);
   
     // triangular face (bottom)
     if ( face == 0 ) { // - counter-clockwise nodes (ouside looking in):
          Point<3> nrml = normalOfTriangle( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                            Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                            Point<3>(XY(1,0),XY(1,1),XY(1,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     // quadrilateral face (front left)
     if ( face == 1 ) { //
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                               Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                               Point<3>(XY(4,0),XY(4,1),XY(4,2)),
                                               Point<3>(XY(3,0),XY(3,1),XY(3,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     // quadrilateral face (front right)
     if ( face == 2 ) { //
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                               Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                               Point<3>(XY(5,0),XY(5,1),XY(5,2)),
                                               Point<3>(XY(4,0),XY(4,1),XY(4,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     // quadrilateral face (back)
     if ( face == 3 ) { //
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                               Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                               Point<3>(XY(5,0),XY(5,1),XY(5,2)),
                                               Point<3>(XY(2,0),XY(2,1),XY(2,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }
   
     // triangular face (top)
     if ( face == 4 ) { //
          Point<3> nrml = normalOfTriangle( Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                            Point<3>(XY(4,0),XY(4,1),XY(4,2)),
                                            Point<3>(XY(5,0),XY(5,1),XY(5,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }
   
 } // end UnitNormalToFace (prism)







 /**

Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
\*/
double64  IsoparametricLinearPrism::AspectRatio()
{
   vector<double64> vec(spe);

   EdgeLengths( vec );

   double64 seg_max(vec[0]), seg_min(vec[0]);

   // find largest segment
   for ( size_t i=1; i<spe; i++ ) {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
     }

   return seg_max / seg_min;
}

/**

@param len The lengths of the 12 segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as simple Euclidian distance between corner
vertices of the edge.
*/
void
IsoparametricLinearPrism::EdgeLengths( std::vector<double64>& len )
{
    double64 sum;
    len.resize(spe);

    // segment 1
    sum     = (XY(1,0)-XY(0,0)) * (XY(1,0)-XY(0,0));
    sum    += (XY(1,1)-XY(0,1)) * (XY(1,1)-XY(0,1));
    sum    += (XY(1,2)-XY(0,2)) * (XY(1,2)-XY(0,2));
    len[0]  = sqrt(sum);
    // segment 2
    sum     = (XY(2,0)-XY(1,0)) * (XY(2,0)-XY(1,0));
    sum    += (XY(2,1)-XY(1,1)) * (XY(2,1)-XY(1,1));
    sum    += (XY(2,2)-XY(1,2)) * (XY(2,2)-XY(1,2));
    len[1]  = sqrt(sum);
    // segment 3
    sum     = (XY(0,0)-XY(2,0)) * (XY(0,0)-XY(2,0));
    sum    += (XY(0,1)-XY(2,1)) * (XY(0,1)-XY(2,1));
    sum    += (XY(0,2)-XY(2,2)) * (XY(0,2)-XY(2,2));
    len[2]  = sqrt(sum);
    // segment 4
    sum     = (XY(4,0)-XY(3,0)) * (XY(4,0)-XY(3,0));
    sum    += (XY(4,1)-XY(3,1)) * (XY(4,1)-XY(3,1));
    sum    += (XY(4,2)-XY(3,2)) * (XY(4,2)-XY(3,2));
    len[3]  = sqrt(sum);
    // segment 5
    sum     = (XY(5,0)-XY(4,0)) * (XY(5,0)-XY(4,0));
    sum    += (XY(5,1)-XY(4,1)) * (XY(5,1)-XY(4,1));
    sum    += (XY(5,2)-XY(4,2)) * (XY(5,2)-XY(4,2));
    len[4]  = sqrt(sum);
    // segment 6
    sum     = (XY(3,0)-XY(5,0)) * (XY(3,0)-XY(5,0));
    sum    += (XY(3,1)-XY(5,1)) * (XY(3,1)-XY(5,1));
    sum    += (XY(3,2)-XY(5,2)) * (XY(3,2)-XY(5,2));
    len[5]  = sqrt(sum);
    // segment 7
    sum     = (XY(3,0)-XY(0,0)) * (XY(3,0)-XY(0,0));
    sum    += (XY(3,1)-XY(0,1)) * (XY(3,1)-XY(0,1));
    sum    += (XY(3,2)-XY(0,2)) * (XY(3,2)-XY(0,2));
    len[6]  = sqrt(sum);
    // segment 8
    sum     = (XY(4,0)-XY(1,0)) * (XY(4,0)-XY(1,0));
    sum    += (XY(4,1)-XY(1,1)) * (XY(4,1)-XY(1,1));
    sum    += (XY(4,2)-XY(1,2)) * (XY(4,2)-XY(1,2));
    len[7] = sqrt(sum);
    // segment 9
    sum     = (XY(5,0)-XY(2,0)) * (XY(5,0)-XY(2,0));
    sum    += (XY(5,1)-XY(2,1)) * (XY(5,1)-XY(2,1));
    sum    += (XY(5,2)-XY(2,2)) * (XY(5,2)-XY(2,2));
    len[8] = sqrt(sum);
 } // end EdgeLengths



/**

Computes Integral N dA = sum_1...n Wi Ji Ni. This formulation also
gives a meaningful area if the element boundaries are curved. For a
description of numerical integration  of isoparametric quadratic
triangular elements, refer to Cook et al. (1989) 3rd Ed., p. 183.

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The area (m2) of the finite element.
*/
double64
IsoparametricLinearPrism::Volume()
{
    double64   area; // determinant
    size_t  i;

    // numerical integration:
    // looping over the 4 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( area=0.0, i=0; i<gpe; i++ )
      {
         dNr( IP(i,0), IP(i,1), IP(i,2), DNR );
         dNs( IP(i,0), IP(i,1), IP(i,2), DNS );
         dNt( IP(i,0), IP(i,1), IP(i,2), DNT );

         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix
         Jacobian( DNR, DNS, DNT );
         area += JacobianInverse() * W[i];
      }


    //cout<<" Volume of Prism="<<area<<endl;

    return area;
}



double64
IsoparametricLinearPrism::VolumeOfRegularPrism()
{
     vector<double64> V1XYZ(dim), V2XYZ(dim), V3XYZ(dim);
     for(size_t i=0; i<dim; i++)
        {
            V1XYZ[i]=XY(0,i);V2XYZ[i]=XY(1,i);V3XYZ[i]=XY(2,i);
        }
     double64  areaBottom=AreaOfBase(V1XYZ,V2XYZ,V3XYZ);

    for(size_t i=0; i<dim; i++)
        {
            V1XYZ[i]=XY(3,i);V2XYZ[i]=XY(4,i);V3XYZ[i]=XY(5,i);
        }

    double64  areaTop=AreaOfBase(V1XYZ,V2XYZ,V3XYZ);
    //cout<<"IsoparametricLinearPrism::Volume Base area="<<areaBottom<<"Top Area="<<areaTop<<endl;

    if(areaTop!=areaBottom)
        cout<<" IsoparametricLinearPrism::Volume ***WARNING: Irregular prism: => Volume aproximate ..."<<endl;

    vector<double64> len(9);
    EdgeLengths(len);

    if(len[6]!=len[7] && len[7]!=len[8] && len[6]!=len[8])
        cout<<" IsoparametricLinearPrism::Volume ***WARNING: Irregular prism => Volume is aproximate ..."<<endl;

    double64 avergeHeight=(len[6]+len[7]+len[8])/3.0;

    return areaBottom*avergeHeight;
}



 /**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the gauss
point i (i=0...i=gauss points-1).

@section arguments Input Arguments

The first argument is a reference to the parent Element object,
the second argument is the result matrix, the third argument indicates
the gauss point at which the derivative matrix is computed and the
fourth argument specifies the degrees of freedom per node, for which
the B matrix shall be transformed.

@param B The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.

*/
double64
IsoparametricLinearPrism::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, size_t gauss_point )
 {
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
    // cout<<" IsoparametricLinearPrism::dN_AtIntegrationPoint: No integration points used ..."<<endl;
    assert( gauss_point < gpe );
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 8 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];
    B(0,4) = DNR[4], B(1,4) = DNS[4], B(2,4) = DNT[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5], B(2,5) = DNT[5];

    B = JINV * B;

    return detJ;
 }



double64
IsoparametricLinearPrism::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    //double64 TwoThirds=2.0*sqrt(0.5)/3.0;

    double64 OneThird=1.0/3.0;

    dNr( OneThird,OneThird, 0.0, DNR );
    dNs( OneThird,OneThird, 0.0, DNS );
    dNt( OneThird,OneThird, 0.0, DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];
    B(0,4) = DNR[4], B(1,4) = DNS[4], B(2,4) = DNT[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5], B(2,5) = DNT[5];

    B = JINV * B;

    return detJ;
 }


 /**

This method attempts to calculate the inner radius of the curved-sided
quad by talking 1/6 of the perimeter and dividing this measure by
the area of the quad.

The procedure is empirically based giving a good match if the quadss
are relatively even-sided and have straight edges.

@section arguments Input Arguments

The parent element is queried for its node coordinates.

*/
double64
IsoparametricLinearPrism::InnerRadius()
{
   vector<double64> segms(spe);
   double64         sum(0.0);

   EdgeLengths( segms );
   for ( size_t i=0; i<spe; i++ ) sum += segms[i];

   //cout<<" IsoparametricLinearPrism::InnerRadius: WARNING: function not tetsted"<<endl;

   if(AspectRatio()>4.0)
   cout<<" IsoparametricLinearPrism::InnerRadius: ***WARNING: function not applicable for CURRENT HAR element"<<endl;


   return Volume() / (sum/6.);
}




/** Returns 0 local node ids of the nodes located at the midsides of
the element.

@param ids Returns an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

*/
void
IsoparametricLinearPrism::MidSideNodes(std::vector<size_t>& ids) const
 {
    cout<<" IsoparametricLinearPrism::MidSideNodes WARNING: MidSideNodes not present "<<endl;
    ids[0]=0;
 }




/**

Method returns a vector with a size of 3, containing the consecutively
ordered local node numbers of the side of the element which lies at the
indicated model boundary. If the element is not on the model boundary an
error is reported and the vector is initialized to unspecified.

@section arguments Input Arguments

The parent Element, the target boundary of the the Model of elements,
and a Meschach vector which will hold the the local indices of the identified
nodes (0...5).

@param fnids The resulting local node id's are returned into the third method argument.

@section implementation Implementation

While only the element knows which nodes are located at the mode boundary,
the FiniteElement knows in which order these appear.

@section application Application

To assign Neumann boundary conditions with a PDE operator for surface
integrals.
*/
void
IsoparametricLinearPrism::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                      vector<size_t>& fnids )
 {
    fnids.resize(bnodes.size());

     if ( bnodes.size() != 4 && bnodes.size() != 3 )
       throw csmp::Exception( CSMP_ERROR, "IsoparametricLinearPrism::ConsecutiveNodesAtBoundary",
               "Cannot resolve node sequence for element boundary",
               "Probably because element lies at two boundaries simultaneously" );


 } // end ConsecutiveNodesAtBoundary




/**

As the linear character of the element do not expect high accuracy, one
integration point is enough and used in most of the applications (FE, V1
by Zienkewitch).

Because in the local coordinate system, the pyramid is always straight
sided, and the derivatives of the interpolation functions are bilinear,
variable values computed at integration points can simply be extrapolated
to the nodes.  This involves the steps:

1. Compute the interpolation function coefficients for the pyramid which
   is defined by the one integration point.

2. For each node point compute the values of the interpolation functions
   and use these to extrapolate the values of the variables at the nodes.
*/
void
IsoparametricLinearPrism::ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                      const vector<double64>& IVAR,
                                                                      vector<double64>& NVAR )
const
{
  if(gpe!=1) throw csmp::Exception( CSMP_ERROR,	 "IsoparametricLinearPrism::PhysicalToParametric",
                                          "Number of Integration points should be 1");

 assert( gpe == 1);
 assert( IVAR.size() >= (gpe*nvars) );

   NVAR.resize( npe * nvars );

  // Define nodal values as bi-linear variation of the integration points values
  // See Zienkewitch, pp. 351, for example
    for ( size_t i=0; i<npe; i++ )
        {
            for ( size_t k=0; k<nvars; k++ )
            {
                NVAR[i*nvars + k] = IVAR[k];
            }
        }
} // end ExtrapolateIntegrationPointVariableToNodes (vectors)






/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be uptodate
void  IsoparametricLinearPrism::IntegrationPoint( size_t ip,
                                                  vector<double64>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(3U);
    xyz[0]=xyz[1]=xyz[2]=0.;
     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

    for( size_t i=0U; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint


void
IsoparametricLinearPrism::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);

    matCoords = NXYZ;
}

void IsoparametricLinearPrism::JacobianAt( const std::vector<double64>& rst )
{
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR,
              DNS,
              DNT );
}

void
IsoparametricLinearPrism::OutputNodeDataToVTK( const char* file_name,
                                       const char* var_name,
                                       DenseMatrix<DM_MIN>& DATA ) const
  {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     sprintf( elmt, "%lu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cout <<"\nIsoparametricLinearPrism::OutputNodeDataToVTK ";
           cout <<"Output file could not be opened."<< endl;
           return;
       }

     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;

     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << npe <<" float"<< endl;
     for ( size_t i=0; i<npe; i++ )
       {
          for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point)) -corect
     // Cells
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< 7 << endl;

     ofs << 6 <<" 0 1 2 3 4 5" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 13 << endl; // VTK_WEDGE
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     // Unfortunately the data can only be output as nodal variables
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );

     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" float"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1x9
           for ( size_t i=0; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 10
          for ( size_t i=0; i<DATA.Cols(); i++ ) {
               for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricLinearPrism::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

 } // end OutputNodeDataToVTK



} // end namespace csmp



