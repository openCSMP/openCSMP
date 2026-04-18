#include "IsoparametricLinearPyramid.h"
#include "Exception.h"
#include "triangularFacet.h"
#include "QuadrilateralFacet.h"

using namespace std;

namespace csmp {

/// @todo (3) Check streched PYRS? - PhysicalToParametric fails???
/// @note Barycenter of pyramid?

IsoparametricLinearPyramid::IsoparametricLinearPyramid( uint32_t integrationPoints )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_PYRAMID, true, true, 1U ),
    NXYZ(5,3),
    IP(integrationPoints,3)
 {
    //AAM, 07.02
    dim = 3;
    itp = 1;
    // 3/4?
    npf = 3;
    npe = 5;
    fpe = 5;
    spe = 8;
    epe = 5;
    // nne - typical number of finite elements, sharing each node
    nne = 6;
    cne = 0;
    // Integration points ???// 0!
    gpe = integrationPoints;

    M.Resize(npe,npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_LINEAR_PYRAMID);

      // base class matrices
    XY.Resize(npe,dim);
    JAC.Resize(dim,dim);
    JINV.Resize(dim,dim);

    // base class vectors
    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    // initializing local node coordinates, r, s, t
    NXYZ(0,0) =-1.0; NXYZ(0,1) =-1.0; NXYZ(0,2) =  0.0;
    NXYZ(1,0) = 1.0; NXYZ(1,1) =-1.0; NXYZ(1,2) =  0.0;
    NXYZ(2,0) = 1.0; NXYZ(2,1) = 1.0; NXYZ(2,2) =  0.0;
    NXYZ(3,0) =-1.0; NXYZ(3,1) = 1.0; NXYZ(3,2) =  0.0;
    NXYZ(4,0) = 0.0; NXYZ(4,1) = 0.0; NXYZ(4,2) =  1.0;

    W.resize( gpe );

    if(integrationPoints==1)
    {
        IP(0,0)=0.0; IP(0,1)=0.0; IP(0,2)=0.25;
        W[0]= 4.0/3.0; // volume is 1.3333333
    }
    else if ( integrationPoints==5 ) {
        IP(0,0)= 0.0; IP(0,1)= 0.0; IP(0,2)= 0.585410196624968;
        IP(1,0)=-0.5; IP(1,1)=-0.5; IP(1,2)= 0.138196601125011;
        IP(2,0)= 0.5; IP(2,1)=-0.5; IP(2,2)= 0.138196601125011;
        IP(3,0)= 0.5; IP(3,1)= 0.5; IP(3,2)= 0.138196601125011;
        IP(4,0)=-0.5; IP(4,1)= 0.5; IP(4,2)= 0.138196601125011;
        W[0] = W[1] = W[2] = W[3] = W[4] =   0.266666666666667;
    }
    else if ( integrationPoints==8 )
        // Generate 8 integration points for full numeric integration + weights
        GenerateIntegrationPoints(IP,W);
    else
    {
        cerr<<"IsoparametricLinearPyramid::IsoparametricLinearHexahedron: Number of integration points should be 1 or 8"<<endl;
        throw std::range_error("***ERROR: IsoparametricLinearPyramid::IsoparametricLinearHexahedron: N of integration points should be 1 or 8");
    }

}





/** Returns the value of the element interpolation functions at the point 'rst'
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
/// Roman, 2013. Corrected(optimised) approximation according to Bedrosian paper
void IsoparametricLinearPyramid::Nrst( double r,
                                       double s,
                                       double t,
                                       std::vector<double>& N ) const
{

   /// Zganski et. al.,1996. "A new family of finite elements: The pyramidal elements"
   /// G.Bedrosian, 1992. "Shape functions for three-dimentional finite element analysis"

    N.resize(npe);

   const double tMinus(1. - t);

   double fRationalTerm( 0. );
   if( t != 1. )
       fRationalTerm = r * s / tMinus;

   N[0] = 0.25*( tMinus + fRationalTerm - r - s );
   N[1] = 0.25*( tMinus - fRationalTerm + r - s );
   N[2] = 0.25*( tMinus + fRationalTerm + r + s );
   N[3] = 0.25*( tMinus - fRationalTerm - r + s );
   N[4] = t;

//   // alternative approximation
//   const double rPlus(1.+r);
//   const double sPlus(1.+s);
//   const double rMinus(1.-r);
//   const double sMinus(1.-s);
//   const double tMinus(1.-t);
//   N[0] = 0.25*rMinus*sMinus*tMinus;
//   N[1] = 0.25*rPlus *sMinus*tMinus;
//   N[2] = 0.25*rPlus *sPlus *tMinus;
//   N[3] = 0.25*rMinus*sPlus *tMinus;
//   N[4] = t;
}

void IsoparametricLinearPyramid::Nrst( double r,
                                       double s,
                                       double t,
                                       double* N ) const
{

   /// Zganski et. al.,1996. "A new family of finite elements: The pyramidal elements"
   /// G.Bedrosian, 1992. "Shape functions for three-dimentional finite element analysis"

   const double tMinus(1.-t);

   double fRationalTerm( 0. );
   if( t != 1. )
       fRationalTerm = r*s/tMinus;

   N[0] = 0.25*( tMinus + fRationalTerm - r - s );
   N[1] = 0.25*( tMinus - fRationalTerm + r - s );
   N[2] = 0.25*( tMinus + fRationalTerm + r + s );
   N[3] = 0.25*( tMinus - fRationalTerm - r + s );
   N[4] = t;

//   // alternative approximation
//   const double rPlus(1.+r);
//   const double sPlus(1.+s);
//   const double rMinus(1.-r);
//   const double sMinus(1.-s);
//   const double tMinus(1.-t);
//   N[0] = 0.25*rMinus*sMinus*tMinus;
//   N[1] = 0.25*rPlus *sMinus*tMinus;
//   N[2] = 0.25*rPlus *sPlus *tMinus;
//   N[3] = 0.25*rMinus*sPlus *tMinus;
//   N[4] = t;
}


/**

This method and the complementary method dNs() compute the shape function
derivates with respect to the local coordinate axis 'r', and 's'.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

@param DNR is the return value, i.e. a vector with the local-coordinate
system interpolation function derivates in the local r coordinate direction.

The shape function derivatives are returned into the second method
argument which is a floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.
*/
void IsoparametricLinearPyramid::dNr ( double,
                                       double s,
                                       double t,
                                       vector<double>& DNR ) const
{
   DNR.resize(npe);

   if ( t != 1. ) {
        // AAM
        DNR[0] = 0.25*(-1.+t+s)/(1.-t);
        DNR[1] =-0.25*(-1.+t+s)/(1.-t);
        DNR[2] = 0.25*( 1.-t+s)/(1.-t);
        DNR[3] =-0.25*( 1.-t+s)/(1.-t);
        DNR[4] = 0.;
        /* text book
        DNR[0] = -0.25 * (1. - s) / (1. - t);
        DNR[1] =  0.25 * (1. - s) / (1. - t);
        DNR[2] =  0.25 * (1. + s) / (1. - t);
        DNR[3] = -0.25 * (1. + s) / (1. - t);
        DNR[4] =  0.;
        */
     }
   else {//t==1
        DNR[0] =-0.25;
        DNR[1] = 0.25;
        DNR[2] = 0.25;
        DNR[3] =-0.25;
        DNR[4] = 0.;
        /*
        DNR[0] = 0.;
        DNR[1] = 0.;
        DNR[2] = 0.;
        DNR[3] = 0.;
        DNR[4] = 0.;
        */
     }

//   // alternative approximation
//   const double sPlus(1.0+s);
//   const double sMinus(1.0-s);
//   const double tMinus(1.-t);

//   DNR[0] = -0.25*sMinus*tMinus;
//   DNR[1] =  0.25*sMinus*tMinus;
//   DNR[2] =  0.25*sPlus *tMinus;
//   DNR[3] = -0.25*sPlus *tMinus;
//   DNR[4] =  0.0;

}

void IsoparametricLinearPyramid::dNs(
                double r,
                double,
                double t,
                vector<double>& DNS ) const
{
   DNS.resize(npe);

   if ( t != 1.0 ) {
        DNS[0] = 0.25*(-1.+t+r)/(1.-t);
        DNS[1] =-0.25*( 1.-t+r)/(1.-t);
        DNS[2] = 0.25*( 1.-t+r)/(1.-t);
        DNS[3] =-0.25*(-1.+t+r)/(1.-t);
        DNS[4] = 0.;
        /*
        DNS[0] = -0.25 * (1. - r) / (1. - t);
        DNS[1] = -0.25 * (1. + r) / (1. - t);
        DNS[2] =  0.25 * (1. + r) / (1. - t);
        DNS[3] =  0.25 * (1. - r) / (1. - t);
        DNS[4] =  0.;
        */
   }
   else {//t==1
        DNS[0] =-0.25;
        DNS[1] =-0.25;
        DNS[2] = 0.25;
        DNS[3] = 0.25;
        DNS[4] = 0.0;
        /*
        DNS[0] = 0.;
        DNS[1] = 0.;
        DNS[2] = 0.;
        DNS[3] = 0.;
        DNS[4] = 0.;
        */
   }

//   // alternative approximation
//   const double rPlus(1.+r);
//   const double rMinus(1.-r);
//   const double tMinus(1.-t);

//   DNS[0] = -0.25*rMinus*tMinus;
//   DNS[1] = -0.25*rPlus *tMinus;
//   DNS[2] =  0.25*rPlus *tMinus;
//   DNS[3] =  0.25*rMinus*tMinus;
//   DNS[4] =  0.0;

}

void IsoparametricLinearPyramid::dNt(
                double r,
                double s,
                double t,
                vector<double>& DNT ) const
{
   DNT.resize(npe);

   if ( t != 1. ) {
       const double RS( r*s );
       const double TMinus2( (1.0-t)*(1.0-t) );
       const double RST(RS/TMinus2);
       DNT[0] = -0.25+0.25*RST;
       DNT[1] = -0.25-0.25*RST;
       DNT[2] = -0.25+0.25*RST;
       DNT[3] = -0.25-0.25*RST;
       DNT[4] =  1.;
       /*
       DNT[0] = -0.25 * (1. - r) / (1. - s);
       DNT[1] = -0.25 * (1. + r) / (1. - s);
       DNT[2] = -0.25 * (1. + r) / (1. + s);
       DNT[3] = -0.25 * (1. - r) / (1. + s);
       DNT[4] =  1.;
       */
   }
   else { //t==1
       DNT[0] =-0.25;
       DNT[1] =-0.25;
       DNT[2] =-0.25;
       DNT[3] =-0.25;
       DNT[4] = 1.0;
       /*
       DNT[0] = 0.;
       DNT[1] = 0.;
       DNT[2] = 0.;
       DNT[3] = 0.;
       DNT[4] = 1.;
       */
    }

//   // alternative approximation
//   const double rPlus(1.+r);
//   const double sPlus(1.+s);
//   const double rMinus(1.-r);
//   const double sMinus(1.-s);

//   DNT[0] = -0.25*rMinus*sMinus;
//   DNT[1] = -0.25*rPlus *sMinus;
//   DNT[2] = -0.25*rPlus *sPlus;
//   DNT[3] = -0.25*rMinus*sPlus;
//   DNT[4] =  1.0;

}



/** Returns local node ids of the 4 nodes located at the corners of
the Linear Pyramid element.

@param ids Returns an unsigned integer vector with the 5 local corner node ID numbers
for the element.

*/
void
IsoparametricLinearPyramid::CornerNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(5);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
}



/** Returns the counter-clockwise local node numbering for the element.

@param ids returns the counter-clockwise local node numbering for the element.

*/
void
IsoparametricLinearPyramid::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
 }



void
IsoparametricLinearPyramid::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
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
         snids[1] = 3;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 3;
         snids[1] = 0;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 0;
         snids[1] = 4;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 1;
         snids[1] = 4;
      }
    else if ( segm_id == 6 ) {
         snids[0] = 2;
         snids[1] = 4;
      }
    else if ( segm_id == 7 ) {
         snids[0] = 3;
         snids[1] = 4;
      }
    else
    std::cout <<"\nIsoparametricLinearPyramid::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment



/**
     Returns the local node numbers of the faces in counter clockwise order,
     from the outside looking in. All number lists start with the smallest node number.
     
     @test SKM OK
*/
/*
void IsoparametricLinearPyramid::NodesOfFace( uint32_t face_id, std::vector<uint32_t>& fnids ) const
 {

    if( face_id == 4 )
      {
         fnids.resize(4);
         fnids[0] = 0;
         fnids[1] = 3;
         fnids[2] = 2;
         fnids[3] = 1;
      }
    else if ( face_id == 0 )
      {
         fnids.resize(3);
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 4;

      }
    else if ( face_id == 1 )
      {
         fnids.resize(3);
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 4;

      }
    else if ( face_id == 2 )
      {
         fnids.resize(3);
         fnids[0] = 2;
         fnids[1] = 3;
         fnids[2] = 4;

      }
    else if ( face_id == 3 )
      {
         fnids.resize(3);
         fnids[0] = 0;
         fnids[1] = 4;
         fnids[2] = 3;
      }
    else
    std::cout <<"\nIsoparametricLinearPyramid::NodesOfFace: Invalid Face ID requested: "<< face_id << std::endl;
 }
*/


uint32_t IsoparametricLinearPyramid::NodesPerFace( uint32_t face_id ) const noexcept
 {
		switch (face_id) {
        case 0: return 3U;
        case 1: return 3U;
        case 2: return 3U;
        case 3: return 3U;
        case 4: return 4U;
      }
    cerr <<"\nIsoparametricLinearPyramid::NodesPerFace: face "<< face_id <<" does not exist.";
    return npf;
 }



vector<uint32_t>  IsoparametricLinearPyramid::NodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{0,1,4};
        case 1: return vector<uint32_t>{1,2,4};
        case 2: return vector<uint32_t>{2,3,4};
        case 3: return vector<uint32_t>{0,4,3};
        case 4: return vector<uint32_t>{0,3,2,1};
      }
    cerr <<"\nIsoparametricLinearPyramid::NodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }


vector<uint32_t>  IsoparametricLinearPyramid::CornerNodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{0,1,4};
        case 1: return vector<uint32_t>{1,2,4};
        case 2: return vector<uint32_t>{2,3,4};
        case 3: return vector<uint32_t>{0,4,3};
        case 4: return vector<uint32_t>{0,3,2,1};
      }
    cerr <<"\nIsoparametricLinearPyramid::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }


vector<uint32_t>  IsoparametricLinearPyramid::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        // local corner node numbers are returned in ascending order
        case 0: return vector<uint32_t>{1,3,4};
        case 1: return vector<uint32_t>{0,2,4};
        case 2: return vector<uint32_t>{1,3,4};
        case 3: return vector<uint32_t>{0,2,4};
        case 4: return vector<uint32_t>{0,1,2,3};
        default:
          cerr <<"\nIsoparametricLinearPyramid::NodesConnectedTo: node "<< node_id <<" does not exist.";
      }
    return vector<uint32_t>{};
  }


CSMP_FEM_TYPE
IsoparametricLinearPyramid::ElementTypeOfFace( uint32_t face )  const
 {
    assert( face < fpe );
    if ( face == 4U ) return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
    return ISOPARAMETRIC_LINEAR_TRIANGLE;
 }

double IsoparametricLinearPyramid::WeightAtIntegrationPoint( uint32_t i ) const { return W[i]; }



/// Method to generate 8 integration points, needed for hex-like full bi-linear interpolation
/// @note Normally for linear pyramid number of integration poins as required by accuracy is 1
void
IsoparametricLinearPyramid::GenerateIntegrationPoints(DenseMatrix<DM_MIN> &Ip, std::vector<double>& We)
{
     const uint32_t numberOfQuadIntegrationPoints=4;
     const uint32_t numberOfDimensionsInPlane=2;
     DenseMatrix<DM_MIN> IPQUAD(numberOfQuadIntegrationPoints,numberOfDimensionsInPlane);
     const double constA=0.577350269189626;

     IPQUAD(0,0)=-constA; IPQUAD(0,1)=-constA;
     IPQUAD(1,0)= constA; IPQUAD(1,1)=-constA;
     IPQUAD(2,0)= constA; IPQUAD(2,1)= constA;
     IPQUAD(3,0)=-constA; IPQUAD(3,1)= constA;

     vector<double> WQUAD(numberOfQuadIntegrationPoints);

     WQUAD[0]=1.0;
     WQUAD[1]=1.0;
     WQUAD[2]=1.0;
     WQUAD[3]=1.0;

     const uint32_t numberOfPointsInTdimension=2;
     vector<double> T(numberOfPointsInTdimension);
     T[0]=0.455848155988775; T[1]=0.877485177344559;
     vector<double> b(numberOfPointsInTdimension);
     b[0]=0.100785882079825; b[1]=0.232547451253508;
     uint32_t i{0U};
     for( uint32_t j{0U};j<numberOfQuadIntegrationPoints;j++)
        {
        for( uint32_t k{0U};k<numberOfPointsInTdimension;k++)
            {
                Ip(i,0)=T[k]*IPQUAD(j,0);
                Ip(i,1)=T[k]*IPQUAD(j,1);
                Ip(i,2)=1-T[k];
                We[i++]=WQUAD[j]*b[k];
            }
        }
}




/**

Volume of terahedron - analytic composition of tet volumes
For degenerate configurations?

 Uses My coordinate system:
         z|
          |
          |
          |
          *
         / \
       /     \
     /         \
   x            y      csp   csp   csp
 to transfer to CSP: x==Z, y==X, z==Y
                     0=>2  1=>0  2==1

@warning Coordinate axes mapped, but work for cube at COFC?? */
double
IsoparametricLinearPyramid::VolumeOfTetra(
                            uint32_t verticeIndex1,
                            uint32_t verticeIndex2,
                            uint32_t verticeIndex3,
                            uint32_t verticeIndex4
                            )
{

    uint32_t i=verticeIndex1,j=verticeIndex2,k=verticeIndex3,l=verticeIndex4;

    double d1Times1=	 XY(j,2)*XY(k,0)*XY(l,1)-XY(j,2)*XY(k,1)*XY(l,0)-XY(j,0)*XY(k,2)*XY(l,1)+
                         XY(j,0)*XY(k,1)*XY(l,2)+XY(j,1)*XY(k,2)*XY(l,0)-XY(j,1)*XY(k,0)*XY(l,2);

    if(d1Times1<0.0) cout<<" IsoparametricLinearPyramid::VolTetra d1 <0.0 "<<endl;
    //if(d1Times1<0.0)	throw std::range_error(" LinearHexahedron::VolTetra d1 <0.0 ");

    double d2TimesXYi0=	-XY(i,2)*XY(k,0)*XY(l,1)+XY(i,2)*XY(k,1)*XY(l,0)+XY(i,2)*XY(j,0)*XY(l,1)-
                         XY(i,2)*XY(j,0)*XY(k,1)-XY(i,2)*XY(j,1)*XY(l,0)+XY(i,2)*XY(j,1)*XY(k,0);

    if(d2TimesXYi0<0.0)	cout<<" IsoparametricLinearPyramid::VolTetra d2TimesXYi0 <0.0"<<endl;
    //if(d2TimesXYi0<0.0) throw std::range_error(" LinearHexahedron::VolTetra d2TimesXYi0 <0.0");

    double d3TimesXYi1= -XY(i,1)*XY(k,2)*XY(l,0)+XY(i,1)*XY(k,0)*XY(l,2)+XY(i,1)*XY(j,2)*XY(l,0)-
                         XY(i,1)*XY(j,2)*XY(k,0)-XY(i,1)*XY(j,0)*XY(l,2)+XY(i,1)*XY(j,0)*XY(k,2);

    if(d3TimesXYi1<0.0) cout<<" IsoparametricLinearPyramid::VolTetra D3TimesXYi1 <0.0"<<endl;
    // throw std::range_error(" LinearHexahedron::VolTetra D3TimesXYi1 <0.0");

    double d4TimesXYi2= -XY(i,1)*XY(k,2)*XY(l,0)+XY(i,1)*XY(k,0)*XY(l,2)+XY(i,1)*XY(j,2)*XY(l,0)-
                         XY(i,1)*XY(j,2)*XY(k,0)-XY(i,1)*XY(j,0)*XY(l,2)+XY(i,1)*XY(j,0)*XY(k,2);

    if(d4TimesXYi2<0.0) cout<< " IsoparametricLinearPyramid::VolTetra d4TimesXYi2 <0.0"<<endl;
    //if(d4TimesXYi2<0.0) throw std::range_error(" LinearHexahedron::VolTetra d4TimesXYi2 <0.0");

    double volume =fabs(d1Times1) + fabs(d2TimesXYi0) + fabs(d3TimesXYi1) + fabs(d4TimesXYi2);

    return volume/6.0;
}



/// @warning Not applicable for irregular pyramid
double
IsoparametricLinearPyramid::VolumeOfPyramid()
{
    //double volume= VolumeOfTetra(0,1,2,4)+VolumeOfTetra(0,2,3,4);
    double volume= VolumeOfTetra(0,1,2,4)+VolumeOfTetra(3,0,2,4);
    cout<<" IsoparametricLinearPyramid::VolumeOfPyramid() V1="<<VolumeOfTetra(0,1,2,4)<<", V2="<<VolumeOfTetra(3,0,2,4)<<endl;
  return volume;
}




/**

Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN5 The interpolation-function derivative matrix is returned into the
second method argument.

*/
void
IsoparametricLinearPyramid::dN( DenseMatrix<DM_MIN>& DN5 )
{
    DN5.Resize(dim,npe);
    M.Resize(dim,1);

     // Jacobian transformation to global coordinate system
     for ( uint32_t i{0U}; i<npe; i++ )
       {
          // here the global coordinates come in
          dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
          dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
          dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

          // compute Jacobian matrix, its determinant and inversex
          Jacobian( DNR, DNS, DNT );
          JacobianInverse();

          ////////// Debug Printout////////////////////////////////
          //cout<<" At vertice #: "<<i<<endl;
          //cout<<" Jacobian Matrix: "<<endl;
          //JAC.Out();
          //cout<<" Jacobian Inverse Matrix: "<<endl;
          //JINV.Out();
          //getchar();
          /////////////////////////////////////////////////////////
          //for ( uint32_t j{0U}; j<dim; j++ ) TEMP(j,0) = DN8(j,i);
          M(0,0)=DNR[i]; M(1,0)=DNS[i]; M(2,0)=DNT[i];

          //DN5(0,i)=DNR[i]; DN5(1,i)=DNS[i]; DN5(2,i)=DNT[i];

          // 3x3 * 3x1 = 3x1 gives the global DN entries
          JINV *= M;
          DN5(0,i) = JINV(0,0);
          DN5(1,i) = JINV(1,0);
          DN5(2,i) = JINV(2,0);
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
double
IsoparametricLinearPyramid::dN( DenseMatrix<DM_MIN>& DN2,
                                const vector<double>& xyz  )
  {

    vector<double> rst(dim);

    PhysicalToParametric(rst, xyz);

    // here find derivatives
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    DN2.Resize(dim,dim);
    DN2  = JINV;

    DenseMatrix<DM_MIN> DN5(dim,npe);
    dN(DN5);
    DN2*=DN5;

    /////////////////////////////// Debug printout ///////////////////////////////////////////////
    //cout<<" IsoparametricLinearPyramid::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
    //rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
    //cout<<" IsoparametricLinearPyramid::dN  Matrix DN2: "<<endl;
    //DN2.Out();
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

    return detJ;

 }


/** Outputs the shape functions at the point 'xyz' which must lie
within the Hexahedron.
*/
void
IsoparametricLinearPyramid::N(
                    std::vector<double>& N,
                    const std::vector<double>& xyz
                    )
{
    vector<double> rst(dim);

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

tested AAM ok*/
double
IsoparametricLinearPyramid::dN_AtNode( DenseMatrix<DM_MIN>& B, uint32_t nd )
 {
    assert( nd < npe );
    dNr( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR );
    dNs( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS );
    dNt( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    // compose matrix DN = 3 x 5 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0]; B(2,0) = DNT[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1]; B(2,1) = DNT[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2]; B(2,2) = DNT[2];
    B(0,3) = DNR[3]; B(1,3) = DNS[3]; B(2,3) = DNT[3];
    B(0,4) = DNR[4]; B(1,4) = DNS[4]; B(2,4) = DNT[4];

    B = JINV * B;

    return detJ;
 }

//********************************************************************************
// Projection function from rst->xyz
//
//********************************************************************************
void
IsoparametricLinearPyramid::ParametricToPhysical(std::vector<double> &rst, std::vector<double>& xyz)
{
    Nrst(rst[0],rst[1],rst[2], NRST );

    for( uint32_t i{0U}; i<dim; i++)
        xyz[i]=0.0;

    for( uint32_t i{0U}; i<npe; i++)
    {
        xyz[0]+=XY(i,0)*NRST[i];
        xyz[1]+=XY(i,1)*NRST[i];
        xyz[2]+=XY(i,2)*NRST[i];
    }
}







//********************************************************************************
// Projection function from xyz->rst, based on Newton-Raphson iteration
//
// Speed of convergence is controlled by the parameter constantMu==1.0-2.0
//********************************************************************************
void
IsoparametricLinearPyramid::PhysicalToParametric(
                                       std::vector<double>& rSt,
                                       const std::vector<double>& xyz
                                    )
{

    vector<double> outxyz(dim);
    vector<double> rstHatK(dim);

    std::vector<double> distanceFromGivenPointLinf(dim,0.0);
    double distanceFromGivenPointL2;

    // Find largest and smallest segments in order to define precision
    vector<double> vec(spe);
    EdgeLengths( vec );
    double seg_max(vec[0]), seg_min(vec[0]);
    for ( uint32_t i=1; i<spe; i++ )
    {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
    }

    const double geometricTolerance = 0.005*seg_min;

    // First guess as BaryCenter
    rstHatK[0] = 0.0;
    rstHatK[1] = 0.0;
    rstHatK[2] = 0.25;

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

        vector<double> rstHatK_PlusOne(dim);
        double minDistanceFromGivenPoint;

        const double constantMu               = 1.0;
        const uint32_t numberOfFirstIterrations   = 5;
        const uint32_t maxNumberOfIterrations     = 20;
        const uint32_t numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double incrementR = 2.0/(numberOfFirstIterrations-1);
        const double incrementS = 2.0/(numberOfFirstIterrations-1);
        const double incrementT = 1.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = -1.0;
        for( uint32_t i{0U};i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = -1.0;
            for( uint32_t j{0U};j<numberOfFirstIterrations;j++)
            {
                rstHatK_PlusOne[2] = 0.0;
                for( uint32_t k{0U};k<numberOfFirstIterrations;k++)
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
                        for( uint32_t l=0; l<dim; l++)
                            rstHatK[l] = rstHatK_PlusOne[l];

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

        uint32_t iteration = 1;

        /// Newton-Raphson iterations
        while ( ( distanceFromGivenPointL2 > geometricTolerance ) && ( iteration < maxNumberOfIterrations ) )
        //while ( ( ( distanceFromGivenPointLinf[0] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[1] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[2] > geometricTolerance )    )
        //      &&  ( iteration < maxNumberOfIterrations )                )
        {
            // Out of range check
            if( (rstHatK[2]<0.) || (rstHatK[2]>1.)||
                (rstHatK[0]<rstHatK[2]-1.0) || (rstHatK[0]>1.0-rstHatK[2]) ||
                (rstHatK[1]<rstHatK[2]-1.0) || (rstHatK[1]>1.0-rstHatK[2])  )
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
                //const double detJ = JacobianDeterminant();
                //if( detJ > 0.0 )
                JacobianInverse();
            }

            // Newton iteration
            rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) +JINV(2,0)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) +JINV(2,1)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[2] = rstHatK[2] - constantMu*(JINV(0,2)*(outxyz[0]-xyz[0]) + JINV(1,2)*(outxyz[1]-xyz[1]) +JINV(2,2)*(outxyz[2]-xyz[2]));

            for( uint32_t i{0U}; i<dim; i++)
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
            cout<<" IsoparametricLinearPyramid::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            cout<<" IsoparametricLinearPyramid::PhysicalToParametric: Real Point:\t"
                <<"x = "<<xyz[0]<<" ;\t"
                <<"y = "<<xyz[1]<<" ;\t"
                <<"z = "<<xyz[2]<<"\n";
            cout<<" IsoparametricLinearPyramid::PhysicalToParametric: Found Point:\t"
               <<"x = "<<outxyz[0]<<" ;\t"
               <<"y = "<<outxyz[1]<<" ;\t"
               <<"z = "<<outxyz[2]<<"\n";
            cout<<" IsoparametricLinearPyramid::PhysicalToParametric: Distance:\t"
               <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
               <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
               <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            cout<<" IsoparametricLinearPyramid::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<" ;\t"
               <<"t = "<<rstHatK[2]<<"\n";

            std::vector<double> N(npe,0.0);
            Nrst(rstHatK[0], rstHatK[1], rstHatK[2], N );
            cout<<" IsoparametricLinearPyramid::PhysicalToParametric: Shape functions:\t"
                <<"N[0] = "<<N[0]<<" ;\t"
                <<"N[1] = "<<N[1]<<" ;\t"
                <<"N[2] = "<<N[2]<<" ;\t"
                <<"N[3] = "<<N[3]<<" ;\t"
                <<"N[4] = "<<N[4]<<"\n";

            csmp::Exception( WARNING, "IsoparametricLinearPyramid::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for( uint32_t i{0U}; i<dim; i++)
      rSt[i] = rstHatK[i];
}





/**
     Uses the TriangularFacet and QuadrlateralFacet classes to compute the normals to the prism faces.
     Based of the face nodes as specified in NodesOfFace()

     @author SKM 15/2/2016
     
     @test OK
*/
void  IsoparametricLinearPyramid::UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const
 {
     assert( face < Faces() );
     unrml.resize(3);
   
     // quadrilateral facet at the basis of pyramid
     if ( face == 4 ) { // 
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                               Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                               Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                               Point<3>(XY(1,0),XY(1,1),XY(1,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     // triangular facet 1
     if ( face == 0 ) { //
          Point<3> nrml = normalOfTriangle( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                            Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                            Point<3>(XY(4,0),XY(4,1),XY(4,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     // triangular facet 2
     if ( face == 1 ) { //
          Point<3> nrml = normalOfTriangle( Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                            Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                            Point<3>(XY(4,0),XY(4,1),XY(4,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     // triangular facet 3
     if ( face == 2 ) { //
          Point<3> nrml = normalOfTriangle( Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                            Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                            Point<3>(XY(4,0),XY(4,1),XY(4,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }
   
     // triangular facet 4
     if ( face == 3 ) { //
          Point<3> nrml = normalOfTriangle( Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                            Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                            Point<3>(XY(4,0),XY(4,1),XY(4,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }
   
 } // end UnitNormalToFace (prism)









/** Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double  IsoparametricLinearPyramid::AspectRatio()
{
   vector<double> vec(spe);

   EdgeLengths( vec );

   double seg_max(vec[0]), seg_min(vec[0]);

   // find largest segment
   for ( uint32_t i=1; i<spe; i++ )
   {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
   }

   return seg_max / seg_min;
}





/** Segment length is calculated directly as the length of linear segments.

@param len The lengths of the 12 segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as simple Euclidian distance between corner
vertices of the edge.

*/
void
IsoparametricLinearPyramid::EdgeLengths( std::vector<double>& len )
{
    double sum;
    len.resize(spe);

    // segment 1 1-2
    sum     = (XY(1,0)-XY(0,0)) * (XY(1,0)-XY(0,0));
    sum    += (XY(1,1)-XY(0,1)) * (XY(1,1)-XY(0,1));
    sum    += (XY(1,2)-XY(0,2)) * (XY(1,2)-XY(0,2));
    len[0]  = sqrt(sum);
    // segment 2 3-2
    sum     = (XY(2,0)-XY(1,0)) * (XY(2,0)-XY(1,0));
    sum    += (XY(2,1)-XY(1,1)) * (XY(2,1)-XY(1,1));
    sum    += (XY(2,2)-XY(1,2)) * (XY(2,2)-XY(1,2));
    len[1]  = sqrt(sum);
    // segment 3 4-3
    sum     = (XY(3,0)-XY(2,0)) * (XY(3,0)-XY(2,0));
    sum    += (XY(3,1)-XY(2,1)) * (XY(3,1)-XY(2,1));
    sum    += (XY(3,2)-XY(2,2)) * (XY(3,2)-XY(2,2));
    len[2]  = sqrt(sum);
    // segment 4 0-4
    sum     = (XY(0,0)-XY(3,0)) * (XY(0,0)-XY(3,0));
    sum    += (XY(0,1)-XY(3,1)) * (XY(0,1)-XY(3,1));
    sum    += (XY(0,2)-XY(3,2)) * (XY(0,2)-XY(3,2));
    len[3]  = sqrt(sum);
    // segment  5-1
    sum     = (XY(4,0)-XY(0,0)) * (XY(4,0)-XY(0,0));
    sum    += (XY(4,1)-XY(0,1)) * (XY(4,1)-XY(0,1));
    sum    += (XY(4,2)-XY(0,2)) * (XY(4,2)-XY(0,2));
    len[4]  = sqrt(sum);
    // segment 6 5-2
    sum     = (XY(4,0)-XY(1,0)) * (XY(4,0)-XY(1,0));
    sum    += (XY(4,1)-XY(1,1)) * (XY(4,1)-XY(1,1));
    sum    += (XY(4,2)-XY(1,2)) * (XY(4,2)-XY(1,2));
    len[5]  = sqrt(sum);
    // segment 7
    sum     = (XY(4,0)-XY(2,0)) * (XY(4,0)-XY(2,0));
    sum    += (XY(4,1)-XY(2,1)) * (XY(4,1)-XY(2,1));
    sum    += (XY(4,2)-XY(2,2)) * (XY(4,2)-XY(2,2));
    len[6]  = sqrt(sum);
    // segment 8
    sum     = (XY(4,0)-XY(3,0)) * (XY(4,0)-XY(3,0));
    sum    += (XY(4,1)-XY(3,1)) * (XY(4,1)-XY(3,1));
    sum    += (XY(4,2)-XY(3,2)) * (XY(4,2)-XY(3,2));
    len[7] = sqrt(sum);

 } // end EdgeLengths




double
IsoparametricLinearPyramid::Volume()
 {
    double  area(0.); // determinant

    // numerical integration:
    // looping over the 8 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( auto i{0U}; i<gpe; ++i )
      {
         dNr( IP(i,0), IP(i,1), IP(i,2), DNR );
         dNs( IP(i,0), IP(i,1), IP(i,2), DNS );
         dNt( IP(i,0), IP(i,1), IP(i,2), DNT );

         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix
         Jacobian( DNR, DNS, DNT );
         area += JacobianInverse() * W[i];
      }

    return area;
 }


/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

@param N The interpolation function values are returned into the third argument.

*/
inline void
IsoparametricLinearPyramid::N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N )
 {
    assert( ip < gpe );

    // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), N );

 } // end N_AtIntegrationPoint


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
double IsoparametricLinearPyramid::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, uint32_t gauss_point )
 {
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
    assert( gauss_point < gpe );

    //if(W[0]==0.0 && W[1]==0.0 && W[2]==0.0 && W[3]==0.0) Volume();

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    // compose matrix DN = 3 x 8 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0]; B(2,0) = DNT[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1]; B(2,1) = DNT[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2]; B(2,2) = DNT[2];
    B(0,3) = DNR[3]; B(1,3) = DNS[3]; B(2,3) = DNT[3];
    B(0,4) = DNR[4]; B(1,4) = DNS[4]; B(2,4) = DNT[4];

    B = JINV * B;

    return detJ;
 }


 void
 IsoparametricLinearPyramid::JacobianAtIntegrationPoint( uint32_t gauss_point )
 {

    assert(gauss_point<gpe);

    //if(W[0]==0.0 && W[1]==0.0 && W[2]==0.0 && W[3]==0.0) Volume();

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }



void
IsoparametricLinearPyramid::N_AtBaryCenter( std::vector<double>& N )
 {
    N.resize(npe);
    Nrst(0.0,0.0,0.25,N);
 }
 
 
 

 /**

This method attempts to calculate the inner radius of the curved-sided
quad by talking 1/2 of the perimeter and dividing this measure by
the area of the quad.

The procedure is empirically based giving a good match if the quadss
are relatively even-sided and have straight edges. Requires further check
for the case of the hexahedron, as it is directly borrowed from the
Tetrahedron

@section arguments Input Arguments

The parent element is queried for its node coordinates.

*/
double
IsoparametricLinearPyramid::InnerRadius()
{
   static  vector<double> segms(spe);
   double                 sum(0.0);

   EdgeLengths( segms );
   for ( uint32_t i{0U}; i<spe; i++ ) sum += segms[i];

   if(AspectRatio()>4.0)
      cerr<<"\nIsoparametricLinearPyramid::InnerRadius: WARNING: function not applicable for this high element aspect ratio.\n"<<endl;

   return Volume() / (sum/6.);
}




/** Returns 0 local node ids of the nodes located at the midsides of
the element.

@param ids Returns an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

*/
void
IsoparametricLinearPyramid::MidSideNodes(std::vector<uint32_t>& ids) const
 {
    cerr<<"\nIsoparametricLinearPyramid::MidSideNodes WARNING: linear element has no MidSideNodes.\n"<<endl;
    ids[0]=0;
 }
 
 
 

/**

Method returns a vector with a size of 3 or 4, containing the consecutively
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
/*
void
IsoparametricLinearPyramid::ConsecutiveNodesAtBoundary( const vector<uint32_t>& bnodes,
                                                        vector<uint32_t>& fnids )
 {
    fnids.resize(bnodes.size());

     if ( bnodes.size() !=4 || bnodes.size() !=3  )
       throw csmp::Exception( ERROR, "IsoparametricLinearPyramid::ConsecutiveNodesAtBoundary",
               "Cannot resolve node sequence for element boundary",
               "Probably because element lies at two boundaries simultaneously" );

 } // end ConsecutiveNodesAtBoundary
*/





double
IsoparametricLinearPyramid::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    // Barycenter of the regular square-sided pyramid is located on the symmetry 1/4 from base
    dNr( 0.0, 0.0, 0.25, DNR );
    dNs( 0.0, 0.0, 0.25, DNS );
    dNt( 0.0, 0.0, 0.25, DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    // compose matrix DN = 3 x 5 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0]; B(1,0) = DNS[0]; B(2,0) = DNT[0];
    B(0,1) = DNR[1]; B(1,1) = DNS[1]; B(2,1) = DNT[1];
    B(0,2) = DNR[2]; B(1,2) = DNS[2]; B(2,2) = DNT[2];
    B(0,3) = DNR[3]; B(1,3) = DNS[3]; B(2,3) = DNT[3];
    B(0,4) = DNR[4]; B(1,4) = DNS[4]; B(2,4) = DNT[4];

    B = JINV * B;

    return detJ;
 }

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
IsoparametricLinearPyramid::ExtrapolateIntegrationPointVariableToNodes(
                                                        uint32_t nvars,
                                                        const vector<double>& IVAR,
                                                        vector<double>& NVAR
                                                                )
const
{
   assert( IVAR.size() >= (gpe*nvars) );
   NVAR.resize( npe * nvars );
  // Define nodal values as bi-linear variation of the integration points values
  // See Zienkewitch, pp. 351, for example
    for ( uint32_t i{0U}; i<npe; i++ )
      for ( uint32_t k{0U}; k<nvars; k++ )
        NVAR[i*nvars + k] = IVAR[k];

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)





// integration point location transformed into global coordinates
// NB: the matrix XYZ must be uptodate
void  IsoparametricLinearPyramid::IntegrationPoint( uint32_t ip,
                                                    vector<double>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(3U);
    xyz[0]=xyz[1]=xyz[2]=0.;

     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

    for( auto i{0U}; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint


void
IsoparametricLinearPyramid::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);

    matCoords = NXYZ;
}


void IsoparametricLinearPyramid::JacobianAt( const std::vector<double>& rst )
{
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR,
              DNS,
              DNT );
}




void IsoparametricLinearPyramid::OutputNodeDataToVTK( const char* file_name,
                                                      const char* var_name,
                                                      DenseMatrix<DM_MIN>& DATA ) const
  {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     snprintf( elmt, sizeof(elmt), "%zu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cout <<"\nIsoparametricLinearPyramid::OutputNodeDataToVTK ";
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
     for ( uint32_t i{0U}; i<npe; i++ )
       {
          for ( uint32_t j{0U}; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point))
     // Cells
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< 6 << endl;
     //
     // pyramid
     ofs << 5 <<" 0 1 2 3 4" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 14 << endl; // VTK_PYRAMID
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
           for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 10
          for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) {
               for ( uint32_t j{0U}; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricLinearPyramid::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

 } // end OutputNodeDataToVTK




// REVIEW
// ========================================================================================
/*  SKM 17/10/2024

    Accuracy of pyramid integration is very low and it is unclear whether this is due to
    issues with the placement of the integration points (8 by default, but locations differ from those
    specified in UG) or whether there is an error with the implementation of the interpolation functions
    and - or their derivatives.
    
    The UG4 code also implements the linear pyramid, using the following 8 quadrature points
    and interpolation functions defined in 'quadrature.cpp' and 'lagrange.h (line 1976ff)' respectively.
    How these shape functions are implented is declared in 'local_shape_function_set.h':
    
    Here are these quadrature points:
    
*/


} // end namespace csmp
