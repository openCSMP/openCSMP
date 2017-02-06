#include "IsoparametricLinearHexahedron.h"
#include "QuadrilateralFacet.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"

using namespace std;

namespace csmp {

IsoparametricLinearHexahedron::IsoparametricLinearHexahedron(size_t integrationPoints)
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_HEXAHEDRON, true, true, 1U ),
    NXYZ(8,3),
    IP(integrationPoints,3)
 {
    assert( integrationPoints == 1 or
            integrationPoints == 4 or
            integrationPoints == 8 );
   
    //AAM, 23.05.02
    dim = 3;
    itp = 1;
    npf = 4;
    npe = 8;
    fpe = 6;
    spe = 12;
    epe = 6;
    // nne - typical number of finite elements, sharing each node
    nne = 6;
    cne = 0;
    // Integration points
    gpe = integrationPoints;

    M.Resize(npe,npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_LINEAR_HEXAHEDRON);

    // base class matrices
    XY.Resize(npe,dim);
    JAC.Resize(dim,dim);
    JINV.Resize(dim,dim);

    // base class vectors
    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    // initializing local node coordinates
    // local node coordinates are defined as r==ksi, s==nu, t==mu.
    NXYZ(0,0) =-1.0;
    NXYZ(1,0) = 1.0;
    NXYZ(2,0) = 1.0;
    NXYZ(3,0) =-1.0;
    NXYZ(4,0) =-1.0;
    NXYZ(5,0) = 1.0;
    NXYZ(6,0) = 1.0;
    NXYZ(7,0) =-1.0;
    // nodal y-coordinates
    NXYZ(0,1) =-1.0;
    NXYZ(1,1) =-1.0;
    NXYZ(2,1) = 1.0;
    NXYZ(3,1) = 1.0;
    NXYZ(4,1) =-1.0;
    NXYZ(5,1) =-1.0;
    NXYZ(6,1) = 1.0;
    NXYZ(7,1) = 1.0;
    // nodal z-coordinates
    NXYZ(0,2) = -1.0;
    NXYZ(1,2) = -1.0;
    NXYZ(2,2) = -1.0;
    NXYZ(3,2) = -1.0;
    NXYZ(4,2) =  1.0;
    NXYZ(5,2) =  1.0;
    NXYZ(6,2) =  1.0;
    NXYZ(7,2) =  1.0;


    W.resize(gpe);
    if( integrationPoints == 1 )
    {
        IP(0,0) = 0.0,	IP(0,1) = 0.0, 	IP(0,2) = 0.0;
        W[0] = 8.0;
    }
    else if( integrationPoints == 4 )
    {
        // 4-point integration scheme
        // Numeric integration 2nd order, number of point m=4
        // Weights of the interpolation functions are 2, see Lo, p. 68
        W[0] = 2.0;
        W[1] = 2.0;
        W[2] = 2.0;
        W[3] = 2.0;

        // Location of the integration points in r-s-t coordinates see Lo, p. 68
        // note that L1 = 1 - r - s - t
        double64 oneDivSqrtThree=1./sqrt(3.0);
        double64 sqrtTwoThirds=sqrt(2.0/3.0);

        IP(0,0) = 0.0,              IP(0,1) = -sqrtTwoThirds, 	IP(0,2) = -oneDivSqrtThree;
        IP(1,0) = 0.0,              IP(1,1) = sqrtTwoThirds,	IP(1,2) = -oneDivSqrtThree;
        IP(2,0) = -sqrtTwoThirds,   IP(2,1) = 0.0, 				IP(2,2) =  oneDivSqrtThree;
        IP(3,0) =  sqrtTwoThirds,   IP(3,1) = 0.0, 				IP(3,2) =  oneDivSqrtThree;
    }
    else if( integrationPoints == 8 )
    {
        // full 8-node integration scheme, see Akin, 1982, p.100
        for(size_t i=0; i<integrationPoints;i++)
            W[i] = 1.0;

        double64 a1 = 0.577350269189626;

        IP(0,0) =-a1,	IP(0,1) =-a1, 	IP(0,2) =-a1;
        IP(1,0) = a1,	IP(1,1) =-a1, 	IP(1,2) =-a1;
        IP(2,0) = a1,	IP(2,1) = a1, 	IP(2,2) =-a1;
        IP(3,0) =-a1,	IP(3,1) = a1, 	IP(3,2) =-a1;
        IP(4,0) =-a1,	IP(4,1) =-a1, 	IP(4,2) = a1;
        IP(5,0) = a1,	IP(5,1) =-a1, 	IP(5,2) = a1;
        IP(6,0) = a1,	IP(6,1) = a1, 	IP(6,2) = a1;
        IP(7,0) =-a1,	IP(7,1) = a1, 	IP(7,2) = a1;
    }
    else
    {
        cout<<"IsoparametricLinearHexahedron::IsoparametricLinearHexahedron: Number of integration points should be 1,4 or 8"<<endl;
        throw std::range_error
        ("***ERROR: IsoparametricLinearHexahedron::IsoparametricLinearHexahedron: N of integration points should be 1,4 or 8");
    }
}




IsoparametricLinearHexahedron::~IsoparametricLinearHexahedron()
 {
 ;
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
 Coordinate axes mapped, but work for cube at COFC??*/
double64
IsoparametricLinearHexahedron::VolumeOfTetra(
                            size_t verticeIndex1,
                            size_t verticeIndex2,
                            size_t verticeIndex3,
                            size_t verticeIndex4 )
{

    size_t i=verticeIndex1,j=verticeIndex2,k=verticeIndex3,l=verticeIndex4;

    double d1Times1=	 XY(j,2)*XY(k,0)*XY(l,1)-XY(j,2)*XY(k,1)*XY(l,0)-XY(j,0)*XY(k,2)*XY(l,1)+
                         XY(j,0)*XY(k,1)*XY(l,2)+XY(j,1)*XY(k,2)*XY(l,0)-XY(j,1)*XY(k,0)*XY(l,2);

    if(d1Times1<0.0) cout<<" IsoparametricLinearHexahedron::VolTetra d1 <0.0 "<<endl;
    //if(d1Times1<0.0)	throw std::range_error(" IsoparametricLinearHexahedron::VolTetra d1 <0.0 ");

    double d2TimesXYi0=	-XY(i,2)*XY(k,0)*XY(l,1)+XY(i,2)*XY(k,1)*XY(l,0)+XY(i,2)*XY(j,0)*XY(l,1)-
                         XY(i,2)*XY(j,0)*XY(k,1)-XY(i,2)*XY(j,1)*XY(l,0)+XY(i,2)*XY(j,1)*XY(k,0);

    if(d2TimesXYi0<0.0)	cout<<" IsoparametricLinearHexahedron::VolTetra d2TimesXYi0 <0.0"<<endl;
    //if(d2TimesXYi0<0.0) throw std::range_error(" IsoparametricLinearHexahedron::VolTetra d2TimesXYi0 <0.0");

    double d3TimesXYi1= -XY(i,1)*XY(k,2)*XY(l,0)+XY(i,1)*XY(k,0)*XY(l,2)+XY(i,1)*XY(j,2)*XY(l,0)-
                         XY(i,1)*XY(j,2)*XY(k,0)-XY(i,1)*XY(j,0)*XY(l,2)+XY(i,1)*XY(j,0)*XY(k,2);

    if(d3TimesXYi1<0.0) cout<<" IsoparametricLinearHexahedron::VolTetra D3TimesXYi1 <0.0"<<endl;
    // throw std::range_error(" IsoparametricLinearHexahedron::VolTetra D3TimesXYi1 <0.0");

    double d4TimesXYi2= -XY(i,1)*XY(k,2)*XY(l,0)+XY(i,1)*XY(k,0)*XY(l,2)+XY(i,1)*XY(j,2)*XY(l,0)-
                         XY(i,1)*XY(j,2)*XY(k,0)-XY(i,1)*XY(j,0)*XY(l,2)+XY(i,1)*XY(j,0)*XY(k,2);

    if(d4TimesXYi2<0.0) cout<< " IsoparametricLinearHexahedron::VolTetra d4TimesXYi2 <0.0"<<endl;
    //if(d4TimesXYi2<0.0) throw std::range_error(" IsoparametricLinearHexahedron::VolTetra d4TimesXYi2 <0.0");

    double64 volume =fabs(d1Times1) + fabs(d2TimesXYi0) + fabs(d3TimesXYi1) + fabs(d4TimesXYi2);

    return volume/6.0;
}





/**

Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

@param ids  an unsigned integer vector with the 4 local corner node ID numbers
for the element.

*/
void
IsoparametricLinearHexahedron::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(8);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
    ids[5] = 5;
    ids[6] = 6;
    ids[7] = 7;
 }



/**

Returns local node numbers.

@param ids  an unsigned integer vector with the 4 local corner node ID numbers
for the element.

*/
void
IsoparametricLinearHexahedron::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
    ids[5] = 5;
    ids[6] = 6;
    ids[7] = 7;
 }



void
IsoparametricLinearHexahedron::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    snids.resize(2);
    /*if ( segm_id == 0 ) {
         snids[0] = 1;
         snids[1] = 0;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 2;
         snids[1] = 1;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 3;
         snids[1] = 2;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 0;
         snids[1] = 3;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 5;
         snids[1] = 4;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 6;
         snids[1] = 5;
      }
    else if ( segm_id == 6 ) {
         snids[0] = 7;
         snids[1] = 6;
      }
    else if ( segm_id == 7 ) {
         snids[0] = 4;
         snids[1] = 7;
      }
    else if ( segm_id == 8 ) {
         snids[0] = 4;
         snids[1] = 0;
      }
    else if ( segm_id == 9 ) {
         snids[0] = 5;
         snids[1] = 1;
      }
    else if ( segm_id == 10 ) {
         snids[0] = 6;
         snids[1] = 2;
      }
    else if ( segm_id == 11 ) {
         snids[0] = 7;
         snids[1] = 3;
      }*/
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
         snids[1] = 5;
      }
    else if ( segm_id == 6 ) {
         snids[0] = 2;
         snids[1] = 6;
      }
    else if ( segm_id == 7 ) {
         snids[0] = 3;
         snids[1] = 7;
      }
    else if ( segm_id == 8 ) {
         snids[0] = 4;
         snids[1] = 5;
      }
    else if ( segm_id == 9 ) {
         snids[0] = 5;
         snids[1] = 6;
      }
    else if ( segm_id == 10 ) {
         snids[0] = 6;
         snids[1] = 7;
      }
    else if ( segm_id == 11 ) {
         snids[0] = 7;
         snids[1] = 4;
      }
 }


/** 
     Face numbering in counter-clockwise order from the outside looking in
     
     @test SKM 2/3/2016

*/
void
IsoparametricLinearHexahedron::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {
    fnids.resize(4);

    if      ( face_id == 0 )
      {
         fnids[0] = 0;
         fnids[1] = 3;
         fnids[2] = 2;
         fnids[3] = 1;
      }
    else if ( face_id == 1 )
      {
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 5;
         fnids[3] = 4;
      }
    else if ( face_id == 2 )
      {
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 6;
         fnids[3] = 5;
      }
    else if ( face_id == 3 )
      {
         fnids[0] = 2;
         fnids[1] = 3;
         fnids[2] = 7;
         fnids[3] = 6;
      }
    else if ( face_id == 4 )
      {
         fnids[0] = 0;
         fnids[1] = 4;
         fnids[2] = 7;
         fnids[3] = 3;
      }
    else if ( face_id == 5 )
      {
         fnids[0] = 4;
         fnids[1] = 5;
         fnids[2] = 6;
         fnids[3] = 7;
      }
    else
    std::cerr <<"\nIsoparametricLinearHexahedron::NodesOfFace: Invalid Face ID requested: "<< face_id << std::endl;
 }


CSMP_FEM_TYPE  IsoparametricLinearHexahedron::ElementTypeOfFace( size_t )  const
 {
    return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
 }



double64 IsoparametricLinearHexahedron::WeightAtIntegrationPoint( size_t i ) const { return W[i]; }





/**

Function gives analytic volume of Hex, based on summ of six consistuting elements of the Tets
Indexing of ANSYS is correct?
*/
double64
IsoparametricLinearHexahedron::VolumeOfHexa()
{
      double64 volume= VolumeOfTetra(0,1,3,4)+VolumeOfTetra(4,1,3,5)+VolumeOfTetra(4,5,3,7)+
                        VolumeOfTetra(1,2,3,6)+VolumeOfTetra(3,1,6,5)+VolumeOfTetra(5,6,3,7);

    cout<<" Vertices:"<<endl;

    for(size_t i=0;i<8;i++)
        cout<<" Hex N_"<<i<<", ("<< XY(i,0)<<","<<XY(i,1)<<","<<XY(i,2)<<");"<<endl;

    cout<<" Volumes : "<<VolumeOfTetra(0,1,3,4)<<" "<<VolumeOfTetra(4,1,3,5)<<" "<<VolumeOfTetra(4,5,3,7)
    <<" "<<VolumeOfTetra(1,2,3,6)<<" "<<VolumeOfTetra(3,1,6,5)<<" "<<VolumeOfTetra(5,6,3,7)<<endl;

    return volume;
}



void
IsoparametricLinearHexahedron::OutputElementToRhino( const char* file_name,
                                                     const bool NodesOn,
                                                     const bool IP_PointsVolOn,
                                                     const bool IP_PointsFaceOn	)
{
    char  outfile[NAME_STRING], elmt[30];
    strcpy( outfile, file_name );
    sprintf( elmt, "%lu", CurrentID() );
    strcat( outfile, elmt );
    strcat( outfile, ".txt" );

     // 0. opening data output file in ascii format
     ofstream rhinos;
     const long precisionOfOutput=20;
     rhinos.precision(precisionOfOutput);
     rhinos.open( outfile, ios::out|ios::trunc );
     if ( !rhinos )
       {
           cout <<"\nIsoparametricLinearHexahedron::OutputElementToRhino ";
           cout <<"Output file could not be opened."<< endl;
           return;
       }

     rhinos<<"# Volumetric FV sector  Vol="<<Volume()<<endl;
     // 1. writing the faces and points as well of the element
     // -----------------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     vector<size_t> fnids(npf);
     for(size_t i=0;i<fpe;i++)
     {

     NodesOfFace(i,fnids);
     rhinos<<"SrfPt ";
        for(size_t j=0;j<npf;j++){
            rhinos<<"(";
            for(size_t k=0;k<dim;k++)
                if(k!=dim-1) rhinos<<COORD(fnids[j],k)<<",";
                    else
                     rhinos<<COORD(fnids[j],k)<<") ";
                }
        rhinos<<" "<<endl;

        if(NodesOn) {
            for(size_t j=0;j<npf;j++){
                rhinos<<"point (";
                for(size_t k=0;k<dim;k++)
                if(k!=dim-1) rhinos<<COORD(fnids[j],k)<<",";
                    else
                     rhinos<<COORD(fnids[j],k)<<") ";
                }
            }
        rhinos<<endl;
     }

 // 2. output of the volume integration points, based on flag IP_PointsVolOn
 // -----------------------------------
 if(IP_PointsVolOn) {
 DenseMatrix<DM_MIN>  IPPHYS;
 IntegrationPointsFromParToPhys(IPPHYS);
 for(size_t i=0;i<gpe;i++)
    {

    rhinos<<"point (";

    for(size_t k=0;k<dim;k++)
        if(k!=dim-1) rhinos<<IPPHYS(i,k)<<",";
                    else
                     rhinos<<IPPHYS(i,k)<<") ";
    }
 }
 // 3. Set up layer for the element
 // -----------------------------------
 rhinos << endl<<"SelAll Layer "<<endl;

 // 4. Output of the faces integration points, based on bool IP_PointsFaceOn
 // -----------------------------------
 if( IP_PointsFaceOn ) {
   vector<size_t>  nIDs(npf);
   rhinos<<"# Element's faces integration points: "<<endl;
     for(size_t i=0;i<fpe;i++)
       {
        IsoparametricLinearQuadrilateral isoLinQuad(3);
        Element<3U>  element1( &isoLinQuad );
        vector<Node<3U> >  nodes(npf);
        NodesOfFace(i,fnids);
        for(size_t j=0;j<npf;j++) {
                 nodes[j].x(COORD(fnids[j],0));
                 nodes[j].y(COORD(fnids[j],1));
                 nodes[j].z(COORD(fnids[j],2));
            }
        nIDs[0]=1;
        nIDs[1]=2;
        nIDs[2]=3;
        nIDs[3]=4;

        DenseMatrix<DM_MIN>  IPPHYSQ;
        isoLinQuad.IntegrationPointsFromParToPhys(IPPHYSQ);
        const size_t numberIPofQuad=4;
        rhinos<<"# Face "<<i+1<<", area="<<element1.Volume()<<endl;
            for(size_t i=0;i<numberIPofQuad;i++)
            {
            rhinos<<"point (";
            for(size_t k=0;k<dim;k++)
            if(k!=dim-1) rhinos<<IPPHYSQ(i,k)<<",";
                        else
                         rhinos<<IPPHYSQ(i,k)<<") ";
            }
           rhinos<<endl;
      }
   } // if IP_PointsFaceOn
   rhinos << endl<<"SelAll Layer "<<endl;
}// finished OutputElementToRhino



/**

Computes the value of the element interpolation functions at the point 'rst'
in local coordinates.

@section arguments Input Arguments:

The floating point coordinates 'r', 's' and 't' of the point at which the
interpolation shall be carried out.

@param N the fourth method argument is the vector into which the values of the
n=nodes interpolation functions at the point 'rst' will be returned.

@section implementation Implementation

See in source code header file.

@section application Application

Method is used to compute property values at the integration points of
the element.
*/
void IsoparametricLinearHexahedron::Nrst(
                    double64 r,
                    double64 s,
                    double64 t,
                    std::vector<double64>& N ) const
{
   const double64 rPlus(1.0+r);
   const double64 sPlus(1.0+s);
   const double64 tPlus(1.0+t);
   const double64 rMinus(1.0-r);
   const double64 sMinus(1.0-s);
   const double64 tMinus(1.0-t);

   N.resize(npe);
   N[0] = 0.125*rMinus*sMinus*tMinus;
   N[1] = 0.125*rPlus*sMinus*tMinus;
   N[2] = 0.125*rPlus*sPlus*tMinus;
   N[3] = 0.125*rMinus*sPlus*tMinus;

   N[4] = 0.125*rMinus*sMinus*tPlus;
   N[5] = 0.125*rPlus*sMinus*tPlus;
   N[6] = 0.125*rPlus*sPlus*tPlus;
   N[7] = 0.125*rMinus*sPlus*tPlus;
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

@code
   //        ri             si               ti
   DNR[0] = -1./8. * (1. + -1. * s) * (1. + -1. * t);
   DNR[1] =  1./8. * (1. + -1. * s) * (1. + -1. * t);
   DNR[2] =  1./8. * (1. +  1. * s) * (1. + -1. * t);
   DNR[3] = -1./8. * (1. +  1. * s) * (1. + -1. * t);
   DNR[4] = -1./8. * (1. + -1. * s) * (1. +  1. * t);
   DNR[5] =  1./8. * (1. + -1. * s) * (1. +  1. * t);
   DNR[6] =  1./8. * (1. +  1. * s) * (1. +  1. * t);
   DNR[7] = -1./8. * (1. +  1. * s) * (1. +  1. * t);
@endcode

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.

 */
void IsoparametricLinearHexahedron::dNr(
                double64 /*r*/,
                double64 s,
                double64 t,
                std::vector<double64>& DNR ) const
{
   const double64 sPlus(1.0+s);
   const double64 tPlus(1.0+t);
   const double64 sMinus(1.0-s);
   const double64 tMinus(1.0-t);

   DNR.resize(npe);
   DNR[0] = -0.125*sMinus*tMinus;
   DNR[1] =  0.125*sMinus*tMinus;
   DNR[2] =  0.125*sPlus*tMinus;
   DNR[3] = -0.125*sPlus*tMinus;
   DNR[4] = -0.125*sMinus*tPlus;
   DNR[5] =  0.125*sMinus*tPlus;
   DNR[6] =  0.125*sPlus*tPlus;
   DNR[7] = -0.125*sPlus*tPlus;
}


/** full form
   //        si             ri               ti
   DNS[0] = -1./8. * (1. + -1. * r) * (1. + -1. * t);
   DNS[1] = -1./8. * (1. +  1. * r) * (1. + -1. * t);
   DNS[2] =  1./8. * (1. +  1. * r) * (1. + -1. * t);
   DNS[3] =  1./8. * (1. + -1. * r) * (1. + -1. * t);
   DNS[4] = -1./8. * (1. + -1. * r) * (1. +  1. * t);
   DNS[5] = -1./8. * (1. +  1. * r) * (1. +  1. * t);
   DNS[6] =  1./8. * (1. +  1. * r) * (1. +  1. * t);
   DNS[7] =  1./8. * (1. + -1. * r) * (1. +  1. * t);
*/
// tested: OK AAM
void IsoparametricLinearHexahedron::dNs(
                double64 r,
                double64 /*s*/,
                double64 t,
                std::vector<double64>& DNS ) const
{
   const double64 rPlus(1.0+r);
   const double64 tPlus(1.0+t);
   const double64 rMinus(1.0-r);
   const double64 tMinus(1.0-t);

   DNS.resize(npe);
   DNS[0] = -0.125*rMinus*tMinus;
   DNS[1] = -0.125*rPlus*tMinus;
   DNS[2] =  0.125*rPlus*tMinus;
   DNS[3] =  0.125*rMinus*tMinus;
   DNS[4] = -0.125*rMinus*tPlus;
   DNS[5] = -0.125*rPlus*tPlus;
   DNS[6] =  0.125*rPlus*tPlus;
   DNS[7] =  0.125*rMinus*tPlus;
}


/** full form
   //        ti             ri               si
   DNT[0] = -1./8. * (1. + -1. * r) * (1. + -1. * s);
   DNT[1] = -1./8. * (1. +  1. * r) * (1. + -1. * s);
   DNT[2] = -1./8. * (1. +  1. * r) * (1. +  1. * s);
   DNT[3] = -1./8. * (1. + -1. * r) * (1. +  1. * s);
   DNT[4] =  1./8. * (1. + -1. * r) * (1. + -1. * s);
   DNT[5] =  1./8. * (1. +  1. * r) * (1. + -1. * s);
   DNT[6] =  1./8. * (1. +  1. * r) * (1. +  1. * s);
   DNT[7] =  1./8. * (1. + -1. * r) * (1. +  1. * s);
*/
void IsoparametricLinearHexahedron::dNt(
                double64 r,
                double64 s,
                double64 /*t*/,
                std::vector<double64>& DNT ) const
{
   const double64 rPlus(1.0+r);
   const double64 sPlus(1.0+s);
   const double64 rMinus(1.0-r);
   const double64 sMinus(1.0-s);

   DNT.resize(npe);
   DNT[0] = -0.125*rMinus*sMinus;
   DNT[1] = -0.125*rPlus*sMinus;
   DNT[2] = -0.125*rPlus*sPlus;
   DNT[3] = -0.125*rMinus*sPlus;
   DNT[4] =  0.125*rMinus*sMinus;
   DNT[5] =  0.125*rPlus*sMinus;
   DNT[6] =  0.125*rPlus*sPlus;
   DNT[7] =  0.125*rMinus*sPlus;
}



/** Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN8 the interpolation-function derivative matrix is returned into the
second method argument.

 */
void
IsoparametricLinearHexahedron::dN( DenseMatrix<DM_MIN>& DN8 )
{
    DN8.Resize(dim,npe);
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
          //cout<<" At vertice #: "<<i<<endl;
          //cout<<" Jacobian Matrix: "<<endl;
          //JAC.Out(cout);
          //cout<<" Jacobian Inverse Matrix: "<<endl;
          //JINV.Out(cout);
          /////////////////////////////////////////////////////////
          M(0,0)=DNR[i]; M(1,0)=DNS[i]; M(2,0)=DNT[i];

          // 3x3 * 3x1 = 3x1 gives the global DN entries
          JINV *= M;
          DN8(0,i) = JINV(0,0);
          DN8(1,i) = JINV(1,0);
          DN8(2,i) = JINV(2,0);

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
IsoparametricLinearHexahedron::dN_At( DenseMatrix<DM_MIN>& DN2,
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

    DenseMatrix<DM_MIN> DN5(dim,npe);
    dN(DN5);
    DN2*=DN5;
    /////////////////////////////// Debug printout ///////////////////////////////////////////////
    cout<<" IsoparametricLinearHexahedron::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
    rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
    cout<<" IsoparametricLinearHexahedron::dN  Matrix DN2: "<<endl;
    DN2.Out(cout);
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

    return detJ;

 }


/**

Outputs the shape functions at the point 'xyz' which must lie
within the Hexahedron.
*/
void
IsoparametricLinearHexahedron::N(
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
IsoparametricLinearHexahedron::dN_AtNode( DenseMatrix<DM_MIN>& B, size_t nd )
 {
    assert( nd < npe );
    dNr( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR );
    dNs( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS );
    dNt( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT );

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
    B(0,6) = DNR[6], B(1,6) = DNS[6], B(2,6) = DNT[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7], B(2,7) = DNT[7];

    B = JINV * B;

    return detJ;
 }

//********************************************************************************
// Projection function from rst->xyz
//
//********************************************************************************
void IsoparametricLinearHexahedron::ParametricToPhysical( vector<double64>& rst,
                                                          vector<double64>& xyz )
{
    Nrst(rst[0],rst[1],rst[2], NRST );

    for(size_t i=0; i<dim; i++)
        xyz[i]=0.0;

    for(size_t i=0; i<npe; i++)
    {
        xyz[0]+=XY(i,0)*NRST[i];
        xyz[1]+=XY(i,1)*NRST[i];
        xyz[2]+=XY(i,2)*NRST[i];

    }
}



void IsoparametricLinearHexahedron::IntegrationPoint( size_t ip,
                                                      vector<double64>& xyz ) const
{
   assert( ip < gpe );
   Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

   xyz.resize( XY.Cols() );
   xyz[0] = xyz[1] = xyz[2] = 0.;

   for ( size_t i=0; i<npe; i++) {
       xyz[0] += XY(i,0) * NRST[i];
       xyz[1] += XY(i,1) * NRST[i];
       xyz[2] += XY(i,2) * NRST[i];
    }
}



void
IsoparametricLinearHexahedron::ReferenceCubeParametricToPhysical( vector<double64>& rst,
                                                                  vector<double64>& xyz )
{

vector<double64> N(npe);

Nrst(rst[0],rst[1],rst[2], N );

for(size_t i=0; i<dim; i++) xyz[i]=0.0;

for(size_t i=0; i<npe; i++)
    {
    xyz[0]+=XY(i,0)*N[i];
    xyz[1]+=XY(i,1)*N[i];
    xyz[2]+=XY(i,2)*N[i];
    }
    cout<<" ("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")"<<endl;

}


//********************************************************************************
// Projection function from xyz->rst
// In first part over regular grid looks for closest point to the given one;
// Then the point is used as initial guess for a Newton-Raphson iteration with
// constant Jacobian matrix
//********************************************************************************

void
IsoparametricLinearHexahedron::PhysicalToParametric(
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
    rstHatK[0] = 0.0;
    rstHatK[1] = 0.0;
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

        const double64 incrementR = 2.0/(numberOfFirstIterrations-1);
        const double64 incrementS = 2.0/(numberOfFirstIterrations-1);
        const double64 incrementT = 2.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = -1.0;
        for(size_t i=0;i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = -1.0;
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
            if( (rstHatK[0]<-1.) || (rstHatK[0]>1.) ||
                (rstHatK[1]<-1.) || (rstHatK[1]>1.) ||
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
            cout<<" IsoparametricLinearHexahedron::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            cout<<" IsoparametricLinearHexahedron::PhysicalToParametric: Real Point:\t"
                <<"x = "<<xyz[0]<<" ;\t"
                <<"y = "<<xyz[1]<<" ;\t"
                <<"z = "<<xyz[2]<<"\n";
            cout<<" IsoparametricLinearHexahedron::PhysicalToParametric: Found Point:\t"
               <<"x = "<<outxyz[0]<<" ;\t"
               <<"y = "<<outxyz[1]<<" ;\t"
               <<"z = "<<outxyz[2]<<"\n";
            cout<<" IsoparametricLinearHexahedron::PhysicalToParametric: Distance:\t"
               <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
               <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
               <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            cout<<" IsoparametricLinearHexahedron::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<" ;\t"
               <<"t = "<<rstHatK[2]<<"\n";

            std::vector<double64> N(npe,0.0);
            Nrst(rstHatK[0], rstHatK[1], rstHatK[2], N );
            cout<<" IsoparametricLinearHexahedron::PhysicalToParametric: Shape functions:\t"
                <<"N[0] = "<<N[0]<<" ;\t"
                <<"N[1] = "<<N[1]<<" ;\t"
                <<"N[2] = "<<N[2]<<" ;\t"
                <<"N[3] = "<<N[3]<<" ;\t"
                <<"N[4] = "<<N[4]<<" ;\t"
                <<"N[5] = "<<N[5]<<" ;\t"
                <<"N[6] = "<<N[6]<<" ;\t"
                <<"N[7] = "<<N[7]<<"\n";

            csmp::Exception( WARNING, "IsoparametricLinearHexahedron::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for(size_t i=0; i<dim; i++)
        rSt[i] = rstHatK[i];

}




/**
     Uses the QuadrilateralFacet to compute the normals to its faces.
     Based of the face nodes as specified in NodesOfFace().

     @author SKM 15/2/2016
     
     @test OK
*/
void  IsoparametricLinearHexahedron::UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const
 {
     assert( face < Faces() );
     unrml.resize(3);
   
     if ( face == 0 ) { // bottom - counter-clockwise nodes (ouside looking in):
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                               Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                               Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                               Point<3>(XY(1,0),XY(1,1),XY(1,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     if ( face == 1 ) { // front
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                               Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                               Point<3>(XY(5,0),XY(5,1),XY(5,2)),
                                               Point<3>(XY(4,0),XY(4,1),XY(4,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     if ( face == 2 ) { // right
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                               Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                               Point<3>(XY(6,0),XY(6,1),XY(6,2)),
                                               Point<3>(XY(5,0),XY(5,1),XY(5,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     if ( face == 3 ) { // back
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                               Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                               Point<3>(XY(7,0),XY(7,1),XY(7,2)),
                                               Point<3>(XY(6,0),XY(6,1),XY(6,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }
       
     if ( face == 4 ) { // left
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                               Point<3>(XY(4,0),XY(4,1),XY(4,2)),
                                               Point<3>(XY(7,0),XY(7,1),XY(7,2)),
                                               Point<3>(XY(3,0),XY(3,1),XY(3,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
          return;
       }

     if ( face == 5 ) { // top
          Point<3> nrml = normalAtFacetCenter( Point<3>(XY(4,0),XY(4,1),XY(4,2)),
                                               Point<3>(XY(5,0),XY(5,1),XY(5,2)),
                                               Point<3>(XY(6,0),XY(6,1),XY(6,2)),
                                               Point<3>(XY(7,0),XY(7,1),XY(7,2)) );
          unrml[0] = nrml[0];
          unrml[1] = nrml[1];
          unrml[2] = nrml[2];
       }
   
 } // end UnitNormalToFace (hexahedron)





 /** Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double64  IsoparametricLinearHexahedron::AspectRatio()
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

//********************************************************************************
// Projection function from rst->xyz
//
//********************************************************************************
/** Segment length is calculated directly as the length of linear segments.


@param len The lengths of the 12 segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as simple Euclidian distance between corner
vertices of the edge.
 */
void
IsoparametricLinearHexahedron::EdgeLengths( std::vector<double64>& len )
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
    sum     = (XY(3,0)-XY(2,0)) * (XY(3,0)-XY(2,0));
    sum    += (XY(3,1)-XY(2,1)) * (XY(3,1)-XY(2,1));
    sum    += (XY(3,2)-XY(2,2)) * (XY(3,2)-XY(2,2));
    len[2]  = sqrt(sum);
    // segment 4
    sum     = (XY(0,0)-XY(3,0)) * (XY(0,0)-XY(3,0));
    sum    += (XY(0,1)-XY(3,1)) * (XY(0,1)-XY(3,1));
    sum    += (XY(0,2)-XY(3,2)) * (XY(0,2)-XY(3,2));
    len[3]  = sqrt(sum);
    // segment 5
    sum     = (XY(5,0)-XY(4,0)) * (XY(5,0)-XY(4,0));
    sum    += (XY(5,1)-XY(4,1)) * (XY(5,1)-XY(4,1));
    sum    += (XY(5,2)-XY(4,2)) * (XY(5,2)-XY(4,2));
    len[4]  = sqrt(sum);
    // segment 6
    sum     = (XY(6,0)-XY(5,0)) * (XY(6,0)-XY(5,0));
    sum    += (XY(6,1)-XY(5,1)) * (XY(6,1)-XY(5,1));
    sum    += (XY(6,2)-XY(5,2)) * (XY(6,2)-XY(5,2));
    len[5]  = sqrt(sum);
    // segment 7
    sum     = (XY(7,0)-XY(6,0)) * (XY(7,0)-XY(6,0));
    sum    += (XY(7,1)-XY(6,1)) * (XY(7,1)-XY(6,1));
    sum    += (XY(7,2)-XY(6,2)) * (XY(7,2)-XY(6,2));
    len[6]  = sqrt(sum);
    // segment 8
    sum     = (XY(4,0)-XY(7,0)) * (XY(4,0)-XY(7,0));
    sum    += (XY(4,1)-XY(7,1)) * (XY(4,1)-XY(7,1));
    sum    += (XY(4,2)-XY(7,2)) * (XY(4,2)-XY(7,2));
    len[7] = sqrt(sum);
    // segment 9
    sum     = (XY(4,0)-XY(0,0)) * (XY(4,0)-XY(0,0));
    sum    += (XY(4,1)-XY(0,1)) * (XY(4,1)-XY(0,1));
    sum    += (XY(4,2)-XY(0,2)) * (XY(4,2)-XY(0,2));
    len[8] = sqrt(sum);
    // segment 10
    sum     = (XY(5,0)-XY(1,0)) * (XY(5,0)-XY(1,0));
    sum    += (XY(5,1)-XY(1,1)) * (XY(5,1)-XY(1,1));
    sum    += (XY(5,2)-XY(1,2)) * (XY(5,2)-XY(1,2));
    len[9] = sqrt(sum);
    // segment 11
    sum     = (XY(6,0)-XY(2,0)) * (XY(6,0)-XY(2,0));
    sum    += (XY(6,1)-XY(2,1)) * (XY(6,1)-XY(2,1));
    sum    += (XY(6,2)-XY(2,2)) * (XY(6,2)-XY(2,2));
    len[10] = sqrt(sum);
    // segment 12
    sum     = (XY(7,0)-XY(3,0)) * (XY(7,0)-XY(3,0));
    sum    += (XY(7,1)-XY(3,1)) * (XY(7,1)-XY(3,1));
    sum    += (XY(7,2)-XY(3,2)) * (XY(7,2)-XY(3,2));
    len[11] = sqrt(sum);
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
IsoparametricLinearHexahedron::Volume()
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

    //area/6?
    //cout<<" Volume of Hex="<<area<<endl;

    return area;
}


/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

@param N the interpolation function values are returned into the third argument.

*/
inline void
IsoparametricLinearHexahedron::N_AtIntegrationPoint( size_t ip, std::vector<double64>& N )

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

@return The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.
*/
double64
IsoparametricLinearHexahedron::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, size_t gauss_point )
 {
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
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
    B(0,6) = DNR[6], B(1,6) = DNS[6], B(2,6) = DNT[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7], B(2,7) = DNT[7];

    B = JINV * B;

    return detJ;
 }

 void
 IsoparametricLinearHexahedron::JacobianAtIntegrationPoint( size_t gauss_point )
 {
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }

// tested: (values come from the regular interpolation functions)
void
IsoparametricLinearHexahedron::N_AtBaryCenter( std::vector<double64>& N )
 {
    N.resize(npe);
    N[0] = 0.125;
    N[1] = 0.125;
    N[2] = 0.125;
    N[3] = 0.125;
    N[4] = 0.125;
    N[5] = 0.125;
    N[6] = 0.125;
    N[7] = 0.125;
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
IsoparametricLinearHexahedron::InnerRadius()
{
   vector<double64>  segms(spe);
   double64          sum(0.0);

   EdgeLengths( segms );
   for ( size_t i=0; i<spe; i++ ) sum += segms[i];

   // Function will produce unrelible value for high aspect ratio elements
   if(AspectRatio()>4.0)
   cout<<" IsoparametricLinearHexahedron::InnerRadius: ***WARNING: function not applicable for CURRENT HAR element"<<endl;

    // Originally sum/2.
   return Volume() / (sum/6.);
}

/** Returns 0 local node ids of the nodes located at the midsides of
the element.

@param ids an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

*/
void
IsoparametricLinearHexahedron::MidSideNodes(std::vector<size_t>& ids) const
 {
    cout<<" IsoparametricLinearHexahedron::MidSideNodes WARNING: MidSideNodes not present "<<endl;
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

@param fnids the resulting local node id's are returned into the third method argument.

@section implementation Implementation

While only the element knows which nodes are located at the mode boundary,
the FiniteElement knows in which order these appear.

@section application Application

To assign Neumann boundary conditions with a PDE operator for surface
integrals.
*/
void
IsoparametricLinearHexahedron::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                           vector<size_t>& fnids )
 {
     fnids.resize(bnodes.size());

     if ( bnodes.size() != 4 )
       throw csmp::Exception( ERROR, "IsoparametricLinearHexahedron::ConsecutiveNodesAtBoundary",
               "Cannot resolve node sequence for element boundary",
               "Probably because element lies at two boundaries simultaneously" );


 } // end ConsecutiveNodesAtBoundary


double64
IsoparametricLinearHexahedron::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    dNr( 0.0, 0.0, 0.0, DNR );
    dNs( 0.0, 0.0, 0.0, DNS );
    dNt( 0.0, 0.0, 0.0, DNT );

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
    B(0,6) = DNR[6], B(1,6) = DNS[6], B(2,6) = DNT[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7], B(2,7) = DNT[7];

    B = JINV * B;

    return detJ;
 }


 void
 IsoparametricLinearHexahedron::IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS)
 {
    IPPHYS.Resize(gpe,dim);
    vector<double64> outxyz(dim);
    vector<double64> rst(dim);

    for(size_t i=0;i<gpe;i++)
        {
        rst[0]=IP(i,0);rst[1]=IP(i,1);rst[2]=IP(i,2);
        ParametricToPhysical( rst, outxyz);
        for(size_t j=0;j<dim;j++) IPPHYS(i,j)=outxyz[j];
        }
   }


  void
  IsoparametricLinearHexahedron::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
  {
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);
    matCoords = NXYZ;
  }

/**

Because, in the local coordinate system, the quadrilateral is always straight
sided, and the derivatives of the interpolation functions are bilinear,
variable values computed at integration points can simply be extrapolated
to the nodes.  This involves the steps:

1. Compute the interpolation function coefficients for the quad which
   is defined by the four integration points.

2. For each node point compute the values of the interpolation functions
   and use these to extrapolate the values of the variables at the nodes.
*/
void
IsoparametricLinearHexahedron::ExtrapolateIntegrationPointVariableToNodes(
                                                        size_t nvars,
                                                        const vector<double64>& IVAR,
                                                        vector<double64>& NVAR
                                                                )
const
{

 assert( gpe == 8 || gpe==1);
 assert( IVAR.size() >= (gpe*nvars) );

   NVAR.resize( npe * nvars );

/*
  const size_t    dim1(2);
  // 0. Building a finite element of specific shape
  // ----------------------------------------------
  LinearQuadrilateral     linear_quadrilateral(dim1);

  // nodes
  Node<dim1>  node1, node2, node3, node4;
  node1.Idx( 1 );
  node1.x(rr[0] );node1.y(ss[0]);
  node2.Idx( 2 );
  node2.x( rr[1]);node2.y( ss[1] );
  node3.Idx( 3 );
  node3.x( rr[2] );node3.y( ss[2] );
  node4.Idx( 4 );
  node4.x( rr[3] );node4.y( ss[3] );

   // element 1
  Element<dim1>  element1( &linear_quadrilateral );
  element1.Idx( 1 );
  element1.ConnectToNodes ( 0, node1 );
  element1.ConnectToNodes ( 1, node2 );
  element1.ConnectToNodes ( 2, node3 );
  element1.ConnectToNodes ( 3, node4 );

*/
  // Define nodal values as bi-linear variation of the integration points values
  // See Zienkewitch, pp. 351, for example

  if(gpe==1)
  {
    for ( size_t i=0; i<npe; i++ )
      for ( size_t k=0; k<nvars; k++ )
        NVAR[i*nvars + k] = IVAR[k];
  }
  else
  if(gpe==8)
  {
    DenseMatrix<DM_MIN> MATRIX_A(gpe,gpe), TEMP_IP(gpe,1), TEMP_N(npe,1);
    // coefficients of the matrix A
    double64 a=0.25*(5+3.0*sqrt(3.0)), b=-0.25*(sqrt(3.0)+1), c=0.25*(sqrt(3.0)-1.0), d=0.25*(5-3.*sqrt(3.));

      MATRIX_A(0,0)=MATRIX_A(1,1)=MATRIX_A(2,2)=MATRIX_A(3,3)
     =MATRIX_A(4,4)=MATRIX_A(5,5)=MATRIX_A(6,6)=MATRIX_A(7,7) =a;

      MATRIX_A(0,1)=MATRIX_A(0,3)=MATRIX_A(0,4)=MATRIX_A(1,0)
     =MATRIX_A(1,2)=MATRIX_A(1,3)=MATRIX_A(1,5)=MATRIX_A(2,1)
     =MATRIX_A(2,3)=MATRIX_A(2,6)=MATRIX_A(3,0)=MATRIX_A(3,2)
     =MATRIX_A(3,7)=MATRIX_A(4,0)=MATRIX_A(4,5)=MATRIX_A(4,7)
     =MATRIX_A(5,1)=MATRIX_A(5,4)=MATRIX_A(5,6)
     =MATRIX_A(6,2)=MATRIX_A(6,5)=MATRIX_A(6,7)
     =MATRIX_A(7,3)=MATRIX_A(7,4)=MATRIX_A(7,6)=b;

      MATRIX_A(0,2)=MATRIX_A(0,5)=MATRIX_A(0,7)
     =MATRIX_A(1,3)=MATRIX_A(1,4)=MATRIX_A(1,6)
     =MATRIX_A(2,0)=MATRIX_A(2,5)=MATRIX_A(2,7)
     =MATRIX_A(3,1)=MATRIX_A(3,4)=MATRIX_A(3,6)
     =MATRIX_A(4,1)=MATRIX_A(4,3)=MATRIX_A(4,6)
     =MATRIX_A(5,0)=MATRIX_A(5,2)=MATRIX_A(5,7)
     =MATRIX_A(6,1)=MATRIX_A(6,3)=MATRIX_A(6,4)
     =MATRIX_A(7,0)=MATRIX_A(7,2)=MATRIX_A(7,5)=c;

      MATRIX_A(0,6)=MATRIX_A(1,7)=MATRIX_A(2,4)=MATRIX_A(3,5)
     =MATRIX_A(4,2)=MATRIX_A(5,3)=MATRIX_A(6,0)=MATRIX_A(7,1)=d;


  for ( size_t i=0; i<nvars; i++ )
    {
        for ( size_t k=0; k<gpe; k++ ) TEMP_IP(k,0)=IVAR[k*nvars +i];

            TEMP_N=MATRIX_A*TEMP_IP;
            //TEMP_N.Out(cout);
            for ( size_t j=0; j<npe; j++ )
            {
                NVAR[j*nvars+i]=TEMP_N(j,0);
            }
    }
  }

 //cout<<" LinearQuadrilateral::ExtrapolateIntegrationPointVariableToNodes Not yet done... "<<endl;

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)





void IsoparametricLinearHexahedron::JacobianAt( const std::vector<double64>& rst )
{
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR,
              DNS,
              DNT );
}

////
// Check the format of the file for VTK representation
////

void
IsoparametricLinearHexahedron::OutputNodeDataToVTK( const char* file_name,
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
           cout <<"\nIsoparametricLinearHexahedron::OutputNodeDataToVTK ";
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
     ofs <<"CELLS "<< 1 <<" "<< 9 << endl;
     // (1+4)X8
     // the 4 corner hexahedra
     ofs << 8 <<" 0 1 2 3 4 5 6 7" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 12 << endl; // VTK_HEX
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
     cout <<"\nIsoparametricLinearHexahedron::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

 } // end OutputNodeDataToVTK



} // end namespace csmp
