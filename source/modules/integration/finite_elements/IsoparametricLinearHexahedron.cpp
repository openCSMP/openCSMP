// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IsoparametricLinearHexahedron.h"
#include "QuadrilateralFacet.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp {

IsoparametricLinearHexahedron::IsoparametricLinearHexahedron( uint32_t integrationPoints )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_LINEAR_HEXAHEDRON, true, true, 1U ),
    NXYZ(8,3), V_(28), Vx_(28), Vy_(28), Vz_(28),
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
        IP(0,0) = 0.0;	IP(0,1) = 0.0; 	IP(0,2) = 0.0;
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
        double oneDivSqrtThree=1./sqrt(3.0);
        double sqrtTwoThirds=sqrt(2.0/3.0);

        IP(0,0) = 0.0;              IP(0,1) = -sqrtTwoThirds; 	IP(0,2) = -oneDivSqrtThree;
        IP(1,0) = 0.0;              IP(1,1) = sqrtTwoThirds;	IP(1,2) = -oneDivSqrtThree;
        IP(2,0) = -sqrtTwoThirds;   IP(2,1) = 0.0; 				IP(2,2) =  oneDivSqrtThree;
        IP(3,0) =  sqrtTwoThirds;   IP(3,1) = 0.0; 				IP(3,2) =  oneDivSqrtThree;
    }
    else if( integrationPoints == 8 )
    {
        // full 8-node integration scheme, see Akin, 1982, p.100
        for(uint32_t i{0U}; i<integrationPoints;i++)
            W[i] = 1.0;

        const double a1 = 0.577350269189626;

        IP(0,0) =-a1;	IP(0,1) =-a1;	IP(0,2) =-a1;
        IP(1,0) = a1;	IP(1,1) =-a1;	IP(1,2) =-a1;
        IP(2,0) = a1;	IP(2,1) = a1;	IP(2,2) =-a1;
        IP(3,0) =-a1;	IP(3,1) = a1;	IP(3,2) =-a1;
        IP(4,0) =-a1;	IP(4,1) =-a1;	IP(4,2) = a1;
        IP(5,0) = a1;	IP(5,1) =-a1;	IP(5,2) = a1;
        IP(6,0) = a1;	IP(6,1) = a1;	IP(6,2) = a1;
        IP(7,0) =-a1;	IP(7,1) = a1;	IP(7,2) = a1;
    }
    else
    {
        cout<<"IsoparametricLinearHexahedron::IsoparametricLinearHexahedron: Number of integration points should be 1,4 or 8"<<endl;
        throw std::range_error
        ("***ERROR: IsoparametricLinearHexahedron::IsoparametricLinearHexahedron: N of integration points should be 1,4 or 8");
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
 Coordinate axes mapped, but work for cube at COFC??*/
double
IsoparametricLinearHexahedron::VolumeOfTetra(
                            uint32_t verticeIndex1,
                            uint32_t verticeIndex2,
                            uint32_t verticeIndex3,
                            uint32_t verticeIndex4 )
{

    const uint32_t i=verticeIndex1,j=verticeIndex2,k=verticeIndex3,l=verticeIndex4;

    const double d1Times1=	 XY(j,2)*XY(k,0)*XY(l,1)-XY(j,2)*XY(k,1)*XY(l,0)-XY(j,0)*XY(k,2)*XY(l,1)+
                         XY(j,0)*XY(k,1)*XY(l,2)+XY(j,1)*XY(k,2)*XY(l,0)-XY(j,1)*XY(k,0)*XY(l,2);

    if(d1Times1<0.0) cout<<" IsoparametricLinearHexahedron::VolTetra d1 <0.0 "<<endl;
    //if(d1Times1<0.0)	throw std::range_error(" IsoparametricLinearHexahedron::VolTetra d1 <0.0 ");

    const double d2TimesXYi0=	-XY(i,2)*XY(k,0)*XY(l,1)+XY(i,2)*XY(k,1)*XY(l,0)+XY(i,2)*XY(j,0)*XY(l,1)-
                         XY(i,2)*XY(j,0)*XY(k,1)-XY(i,2)*XY(j,1)*XY(l,0)+XY(i,2)*XY(j,1)*XY(k,0);

    if(d2TimesXYi0<0.0)	cout<<" IsoparametricLinearHexahedron::VolTetra d2TimesXYi0 <0.0"<<endl;
    //if(d2TimesXYi0<0.0) throw std::range_error(" IsoparametricLinearHexahedron::VolTetra d2TimesXYi0 <0.0");

    const double d3TimesXYi1= -XY(i,1)*XY(k,2)*XY(l,0)+XY(i,1)*XY(k,0)*XY(l,2)+XY(i,1)*XY(j,2)*XY(l,0)-
                         XY(i,1)*XY(j,2)*XY(k,0)-XY(i,1)*XY(j,0)*XY(l,2)+XY(i,1)*XY(j,0)*XY(k,2);

    if(d3TimesXYi1<0.0) cout<<" IsoparametricLinearHexahedron::VolTetra D3TimesXYi1 <0.0"<<endl;
    // throw std::range_error(" IsoparametricLinearHexahedron::VolTetra D3TimesXYi1 <0.0");

    const double d4TimesXYi2= -XY(i,1)*XY(k,2)*XY(l,0)+XY(i,1)*XY(k,0)*XY(l,2)+XY(i,1)*XY(j,2)*XY(l,0)-
                         XY(i,1)*XY(j,2)*XY(k,0)-XY(i,1)*XY(j,0)*XY(l,2)+XY(i,1)*XY(j,0)*XY(k,2);

    if(d4TimesXYi2<0.0) cout<< " IsoparametricLinearHexahedron::VolTetra d4TimesXYi2 <0.0"<<endl;
    //if(d4TimesXYi2<0.0) throw std::range_error(" IsoparametricLinearHexahedron::VolTetra d4TimesXYi2 <0.0");

    const double volume =fabs(d1Times1) + fabs(d2TimesXYi0) + fabs(d3TimesXYi1) + fabs(d4TimesXYi2);

    return volume/6.0;
}





/**

Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

@param ids  an unsigned integer vector with the 4 local corner node ID numbers
for the element.

*/
void
IsoparametricLinearHexahedron::CornerNodes( std::vector<uint32_t>& ids ) const
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
IsoparametricLinearHexahedron::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
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
IsoparametricLinearHexahedron::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
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
// SKM refactored 9/8/22, returns empty vector if face_id is out of range

vector<uint32_t>  IsoparametricLinearHexahedron::NodesOfFace( uint32_t face_id ) const
 {
    assert( face_id < Faces() );
    switch (face_id) {
        case 0: return vector<uint32_t>{ 0, 3, 2, 1 };
        case 1: return vector<uint32_t>{ 0, 1, 5, 4 };
        case 2: return vector<uint32_t>{ 1, 2, 6, 5 };
        case 3: return vector<uint32_t>{ 2, 3, 7, 6 };
        case 4: return vector<uint32_t>{ 0, 4, 7, 3 };
        case 5: return vector<uint32_t>{ 4, 5, 6, 7 };
      }
    cerr <<"\nIsoparametricLinearHexahedron::NodesOfFace: Invalid Face ID requested: "<< face_id << endl;
    return vector<uint32_t>{};
    
 } // end NodesOfFace






vector<uint32_t>  IsoparametricLinearHexahedron::CornerNodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{0,3,2,1};
        case 1: return vector<uint32_t>{0,1,5,4};
        case 2: return vector<uint32_t>{1,2,6,5};
        case 3: return vector<uint32_t>{2,3,7,6};
        case 4: return vector<uint32_t>{0,4,7,3};
        case 5: return vector<uint32_t>{4,5,6,7};
      }
    cerr <<"\nIsoparametricLinearHexahedron::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }



/// returns the local  numbers of the nodes at the other end of the sgment that the argument node is on
vector<uint32_t>  IsoparametricLinearHexahedron::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        // local corner node numbers are returned in ascending order
        case 0: return vector<uint32_t>{1,3,4};
        case 1: return vector<uint32_t>{0,2,5};
        case 2: return vector<uint32_t>{1,3,6};
        case 3: return vector<uint32_t>{0,2,7};
        case 4: return vector<uint32_t>{0,5,7};
        case 5: return vector<uint32_t>{1,4,6};
        case 6: return vector<uint32_t>{2,5,7};
        case 7: return vector<uint32_t>{3,4,6};
        default:
          cerr <<"\nIsoparametricLinearHexahedron::NodesConnectedTo: node "<< node_id <<" does not exist.";
      }
    return vector<uint32_t>{};
  }




CSMP_FEM_TYPE  IsoparametricLinearHexahedron::ElementTypeOfFace( uint32_t )  const
 {
    return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
 }



double IsoparametricLinearHexahedron::WeightAtIntegrationPoint( uint32_t i ) const { return W[i]; }





/**

Function gives analytic volume of Hex, based on summ of six consistuting elements of the Tets
Indexing of ANSYS is correct?
*/
double IsoparametricLinearHexahedron::VolumeOfHexa()
{
      double volume= VolumeOfTetra(0,1,3,4)+VolumeOfTetra(4,1,3,5)+VolumeOfTetra(4,5,3,7)+
                     VolumeOfTetra(1,2,3,6)+VolumeOfTetra(3,1,6,5)+VolumeOfTetra(5,6,3,7);

    cout<<" Vertices:"<<endl;

    for( uint32_t i{0U};i<8;i++)
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
    snprintf( elmt, sizeof(elmt), "%zu", CurrentID() );
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
     for( uint32_t i{0U};i<fpe;i++)
     {
        auto fnids = NodesOfFace(i);
        rhinos<<"SrfPt ";
        for( uint32_t j{0U};j<npf;j++){
            rhinos<<"(";
            for( uint32_t k{0U};k<dim;k++)
                if(k!=dim-1) rhinos<<COORD(fnids[j],k)<<",";
                    else
                     rhinos<<COORD(fnids[j],k)<<") ";
                }
        rhinos<<" "<<endl;

        if(NodesOn) {
            for( uint32_t j{0U};j<npf;j++){
                rhinos<<"point (";
                for( uint32_t k{0U};k<dim;k++)
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
 for( uint32_t i{0U};i<gpe;i++)
    {

    rhinos<<"point (";

    for( uint32_t k{0U};k<dim;k++)
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
   vector<uint32_t>  nIDs(npf);
   rhinos<<"# Element's faces integration points: "<<endl;
     for(auto i{0U};i<fpe;i++)
       {
        IsoparametricLinearQuadrilateral isoLinQuad(3);
        Element<3U>  element1( &isoLinQuad );
        vector<Node<3U> >  nodes(npf);
        auto fnids = NodesOfFace(i);
        for( uint32_t j{0U};j<npf;j++) {
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
        const uint32_t numberIPofQuad=4;
        rhinos<<"# Face "<<i+1<<", area="<<element1.Volume()<<endl;
            for( uint32_t n{0u};n<numberIPofQuad;n++)
            {
            rhinos<<"point (";
            for(uint32_t k{0u};k<dim;k++)
            if(k!=dim-1) rhinos<<IPPHYSQ(n,k)<<",";
                        else
                         rhinos<<IPPHYSQ(n,k)<<") ";
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
                    double r,
                    double s,
                    double t,
                    std::vector<double>& N ) const
{
   const double rPlus(1.0+r);
   const double sPlus(1.0+s);
   const double tPlus(1.0+t);
   const double rMinus(1.0-r);
   const double sMinus(1.0-s);
   const double tMinus(1.0-t);

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

void IsoparametricLinearHexahedron::Nrst(
                    double r,
                    double s,
                    double t,
                    double* N ) const
{
   const double rPlus(1.0+r);
   const double sPlus(1.0+s);
   const double tPlus(1.0+t);
   const double rMinus(1.0-r);
   const double sMinus(1.0-s);
   const double tMinus(1.0-t);

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
                double /*r*/,
                double s,
                double t,
                std::vector<double>& DNR ) const
{
   const double sPlus(1.0+s);
   const double tPlus(1.0+t);
   const double sMinus(1.0-s);
   const double tMinus(1.0-t);

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
                double r,
                double /*s*/,
                double t,
                std::vector<double>& DNS ) const
{
   const double rPlus(1.0+r);
   const double tPlus(1.0+t);
   const double rMinus(1.0-r);
   const double tMinus(1.0-t);

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
                double r,
                double s,
                double /*t*/,
                std::vector<double>& DNT ) const
{
   const double rPlus(1.0+r);
   const double sPlus(1.0+s);
   const double rMinus(1.0-r);
   const double sMinus(1.0-s);

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
double
IsoparametricLinearHexahedron::dN_At( DenseMatrix<DM_MIN>& DN2,
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
    cout<<" IsoparametricLinearHexahedron::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
    rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
    cout<<" IsoparametricLinearHexahedron::dN  Matrix DN2: "<<endl;
    DN2.Out();
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

    return detJ;

 }


/**

Outputs the shape functions at the point 'xyz' which must lie
within the Hexahedron.
*/
void
IsoparametricLinearHexahedron::N(
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
*/
double
IsoparametricLinearHexahedron::dN_AtNode( DenseMatrix<DM_MIN>& BMAT, uint32_t nd )
 {
    assert( nd < npe );
    dNr( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR );
    dNs( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS );
    dNt( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    BMAT.Resize(dim,npe);
    BMAT(0,0) = DNR[0]; BMAT(1,0) = DNS[0]; BMAT(2,0) = DNT[0];
    BMAT(0,1) = DNR[1]; BMAT(1,1) = DNS[1]; BMAT(2,1) = DNT[1];
    BMAT(0,2) = DNR[2]; BMAT(1,2) = DNS[2]; BMAT(2,2) = DNT[2];
    BMAT(0,3) = DNR[3]; BMAT(1,3) = DNS[3]; BMAT(2,3) = DNT[3];
    BMAT(0,4) = DNR[4]; BMAT(1,4) = DNS[4]; BMAT(2,4) = DNT[4];
    BMAT(0,5) = DNR[5]; BMAT(1,5) = DNS[5]; BMAT(2,5) = DNT[5];
    BMAT(0,6) = DNR[6]; BMAT(1,6) = DNS[6]; BMAT(2,6) = DNT[6];
    BMAT(0,7) = DNR[7]; BMAT(1,7) = DNS[7]; BMAT(2,7) = DNT[7];

    BMAT = JINV * BMAT;

    return detJ;
 }

//********************************************************************************
// Projection function from rst->xyz
//
//********************************************************************************
void IsoparametricLinearHexahedron::ParametricToPhysical( vector<double>& rst,
                                                          vector<double>& xyz )
{
    Nrst(rst[0],rst[1],rst[2], NRST );

    for(uint32_t i{0U}; i<dim; i++)
        xyz[i]=0.0;

    for(uint32_t i{0U}; i<npe; i++)
    {
        xyz[0]+=XY(i,0)*NRST[i];
        xyz[1]+=XY(i,1)*NRST[i];
        xyz[2]+=XY(i,2)*NRST[i];

    }
}



void IsoparametricLinearHexahedron::IntegrationPoint( uint32_t ip,
                                                      vector<double>& xyz ) const
{
   assert( ip < gpe );
   Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

   xyz.resize( XY.Cols() );
   xyz[0] = xyz[1] = xyz[2] = 0.;

   for ( uint32_t i{0U}; i<npe; i++) {
       xyz[0] += XY(i,0) * NRST[i];
       xyz[1] += XY(i,1) * NRST[i];
       xyz[2] += XY(i,2) * NRST[i];
    }
}



void
IsoparametricLinearHexahedron::ReferenceCubeParametricToPhysical( vector<double>& rst,
                                                                  vector<double>& xyz )
{

vector<double> N(npe);

Nrst(rst[0],rst[1],rst[2], N );

for(uint32_t i{0U}; i<dim; i++) xyz[i]=0.0;

for(uint32_t i{0U}; i<npe; i++)
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
    for ( uint32_t i{1u}; i<spe; i++ )
    {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
    }

    const double geometricTolerance = 0.005*seg_min;

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

        vector<double> rstHatK_PlusOne(dim);
        double minDistanceFromGivenPoint;

        const double constantMu               = 1.0;
        const auto numberOfFirstIterrations   = 5;
        const auto maxNumberOfIterrations     = 20;
        const auto numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double incrementR = 2.0/(numberOfFirstIterrations-1);
        const double incrementS = 2.0/(numberOfFirstIterrations-1);
        const double incrementT = 2.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = -1.0;
        for(auto i{0U};i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = -1.0;
            for(auto j{0U};j<numberOfFirstIterrations;j++)
            {
                rstHatK_PlusOne[2] = -1.0;
                for(auto k=0;k<numberOfFirstIterrations;k++)
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
                        for(uint32_t l=0u; l<dim; l++)
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

        auto iteration = 1;

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
                //const double detJ = JacobianDeterminant();
                //if( detJ > 0.0 )
                JacobianInverse();
            }

            // Newton iteration
            rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) +JINV(2,0)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) +JINV(2,1)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[2] = rstHatK[2] - constantMu*(JINV(0,2)*(outxyz[0]-xyz[0]) + JINV(1,2)*(outxyz[1]-xyz[1]) +JINV(2,2)*(outxyz[2]-xyz[2]));

            for(auto i{0U}; i<dim; i++)
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

            std::vector<double> N(npe,0.0);
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

    for(auto i{0U}; i<dim; i++)
        rSt[i] = rstHatK[i];

}




/**
     Uses the QuadrilateralFacet to compute the normals to its faces.
     Based of the face nodes as specified in NodesOfFace().

     @author SKM 15/2/2016
     
     @test OK
*/
void  IsoparametricLinearHexahedron::UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const
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
double  IsoparametricLinearHexahedron::AspectRatio()
{
   vector<double> vec(spe);

   EdgeLengths( vec );

   double seg_max(vec[0]), seg_min(vec[0]);

   // find largest segment
   for ( uint32_t i{1u}; i<spe; i++ ) {
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
IsoparametricLinearHexahedron::EdgeLengths( std::vector<double>& len )
{
    double sum;
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
double
IsoparametricLinearHexahedron::Volume()
{
    double   area{0.}; // determinant

    // numerical integration:
    // looping over the 4 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( auto i{0U}; i<gpe; i++ )
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
IsoparametricLinearHexahedron::N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N )

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
double
IsoparametricLinearHexahedron::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& BMAT, uint32_t gauss_point )
 {
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
    assert( gauss_point < gpe );
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    // compose matrix DN = 3 x 8 in global coordinates
    // by multiplication of JINV with local DN
    BMAT.Resize(dim,npe);
    BMAT(0,0) = DNR[0]; BMAT(1,0) = DNS[0]; BMAT(2,0) = DNT[0];
    BMAT(0,1) = DNR[1]; BMAT(1,1) = DNS[1]; BMAT(2,1) = DNT[1];
    BMAT(0,2) = DNR[2]; BMAT(1,2) = DNS[2]; BMAT(2,2) = DNT[2];
    BMAT(0,3) = DNR[3]; BMAT(1,3) = DNS[3]; BMAT(2,3) = DNT[3];
    BMAT(0,4) = DNR[4]; BMAT(1,4) = DNS[4]; BMAT(2,4) = DNT[4];
    BMAT(0,5) = DNR[5]; BMAT(1,5) = DNS[5]; BMAT(2,5) = DNT[5];
    BMAT(0,6) = DNR[6]; BMAT(1,6) = DNS[6]; BMAT(2,6) = DNT[6];
    BMAT(0,7) = DNR[7]; BMAT(1,7) = DNS[7]; BMAT(2,7) = DNT[7];

    BMAT = JINV * BMAT;

    return detJ;
 }

 void
 IsoparametricLinearHexahedron::JacobianAtIntegrationPoint( uint32_t gauss_point )
 {
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }

// tested: (values come from the regular interpolation functions)
void
IsoparametricLinearHexahedron::N_AtBaryCenter( std::vector<double>& N )
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
double
IsoparametricLinearHexahedron::InnerRadius()
{
   vector<double>  segms(spe);
   double          sum(0.0);

   EdgeLengths( segms );
   for ( auto i{0U}; i<spe; i++ ) sum += segms[i];

   // Function will produce unreliable value for high aspect ratio elements
   if (AspectRatio()>4.0)
     cerr<<"\nIsoparametricLinearHexahedron::InnerRadius: WARNING: function not applicable for this high element aspect ratio.\n"<<endl;

    // Originally sum/2.
   return Volume() / (sum/6.);
}

/** Returns 0 local node ids of the nodes located at the midsides of
the element.

@param ids an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

*/
void
IsoparametricLinearHexahedron::MidSideNodes(std::vector<uint32_t>& ids) const
 {
    cerr<<"\nIsoparametricLinearHexahedron::MidSideNodes: WARNING: MidSideNodes not present.\n"<<endl;
    ids[0]=0;
 }

 

double
IsoparametricLinearHexahedron::dN_AtBarycenter( DenseMatrix<DM_MIN>& BMAT )
 {
    dNr( 0.0, 0.0, 0.0, DNR );
    dNs( 0.0, 0.0, 0.0, DNS );
    dNt( 0.0, 0.0, 0.0, DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    BMAT.Resize(dim,npe);
    BMAT(0,0) = DNR[0]; BMAT(1,0) = DNS[0]; BMAT(2,0) = DNT[0];
    BMAT(0,1) = DNR[1]; BMAT(1,1) = DNS[1]; BMAT(2,1) = DNT[1];
    BMAT(0,2) = DNR[2]; BMAT(1,2) = DNS[2]; BMAT(2,2) = DNT[2];
    BMAT(0,3) = DNR[3]; BMAT(1,3) = DNS[3]; BMAT(2,3) = DNT[3];
    BMAT(0,4) = DNR[4]; BMAT(1,4) = DNS[4]; BMAT(2,4) = DNT[4];
    BMAT(0,5) = DNR[5]; BMAT(1,5) = DNS[5]; BMAT(2,5) = DNT[5];
    BMAT(0,6) = DNR[6]; BMAT(1,6) = DNS[6]; BMAT(2,6) = DNT[6];
    BMAT(0,7) = DNR[7]; BMAT(1,7) = DNS[7]; BMAT(2,7) = DNT[7];

    BMAT = JINV * BMAT;

    return detJ;
 }


 void
 IsoparametricLinearHexahedron::IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS)
 {
    IPPHYS.Resize(gpe,dim);
    vector<double> outxyz(dim);
    vector<double> rst(dim);

    for(auto i{0U};i<gpe;i++)
        {
        rst[0]=IP(i,0);rst[1]=IP(i,1);rst[2]=IP(i,2);
        ParametricToPhysical( rst, outxyz);
        for(auto j{0U};j<dim;j++) IPPHYS(i,j)=outxyz[j];
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
                                                        uint32_t nvars,
                                                        const vector<double>& IVAR,
                                                        vector<double>& NVAR
                                                                )
const
{

 assert( gpe == 8 || gpe==1);
 assert( IVAR.size() >= (gpe*nvars) );

   NVAR.resize( npe * nvars );

/*
  const auto    dim1(2);
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
    for ( uint32_t i{0U}; i<npe; i++ )
      for ( uint32_t k=0u; k<nvars; k++ )
        NVAR[i*nvars + k] = IVAR[k];
  }
  else
  if(gpe==8)
  {
    DenseMatrix<DM_MIN> MATRIX_A(gpe,gpe), TEMP_IP(gpe,1), TEMP_N(npe,1);
    // coefficients of the matrix A
    double a=0.25*(5+3.0*sqrt(3.0)), b=-0.25*(sqrt(3.0)+1), c=0.25*(sqrt(3.0)-1.0), d=0.25*(5-3.*sqrt(3.));

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


  for ( uint32_t i{0U}; i<nvars; i++ )
    {
        for ( uint32_t k=0; k<gpe; k++ ) TEMP_IP(k,0)=IVAR[k*nvars +i];

            TEMP_N=MATRIX_A*TEMP_IP;
            //TEMP_N.Out();
            for ( auto j{0U}; j<npe; j++ )
            {
                NVAR[j*nvars+i]=TEMP_N(j,0);
            }
    }
  }

 //cout<<" LinearQuadrilateral::ExtrapolateIntegrationPointVariableToNodes Not yet done... "<<endl;

} // end ExtrapolateIntegrationPointVariableToNodes (vectors)





void IsoparametricLinearHexahedron::JacobianAt( const std::vector<double>& rst )
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
     // NB: Coordinate matrix must already be initialised
     outputDataToVTK( file_name, var_name, VTK_HEXAHEDRON, CurrentID(), XY, DATA );

 } // end OutputNodeDataToVTK




// The analytical integration is correct for aligned or rotated cuboids
void IsoparametricLinearHexahedron::Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K)
{
	double vol36 = 36.*Volume(), 
		     dxy = K(2, 2) * (XY(6, 0) - XY(0, 0))*(XY(6, 0) - XY(0, 0)) * (XY(6, 1) - XY(0, 1))*(XY(6, 1) - XY(0, 1)) / (vol36),
		     dxz = K(1, 1) * (XY(6, 0) - XY(0, 0))*(XY(6, 0) - XY(0, 0)) * (XY(6, 2) - XY(0, 2))*(XY(6, 2) - XY(0, 2)) / (vol36),
		     dyz = K(0, 0) * (XY(6, 1) - XY(0, 1))*(XY(6, 1) - XY(0, 1)) * (XY(6, 2) - XY(0, 2))*(XY(6, 2) - XY(0, 2)) / (vol36), 
		    dxyz = 4.*(dyz + dxz + dxy),
             xpy = dxy + dyz,   
		     xpz = dxy + dxz, 
		     ypz = dyz + dxz, 
		   xpypz = dxy + dyz + dxz;

	V_.resize(36);
	V_ = {dxyz, 2*xpz - 4*dyz, -2*ypz + dxy, 2*xpy - 4*dxz, -4*dxy + 2*ypz, -2*xpy + dxz, -xpypz, -2*xpz + dyz,
		  dxyz, 2*xpy - 4*dxz,  dxy - 2*ypz,  -2*xpz + dxz, -4*dxy + 2*ypz,  dyz - 2*xpz, -xpypz,
		  dxyz, 2*xpz - 4*dyz,       -xpypz,   dyz - 2*xpz, -4*dxy + 2*ypz, -2*xpy + dxz,
		  dxyz,   dyz - 2*xpz,       -xpypz,  -2*xpy + dxz, -4*dxy + 2*ypz,
		  dxyz, 2*xpz - 4*dyz, -2*ypz + dxy, 2*xpy - 4*dxz,
		  dxyz, 2*xpy - 4*dxz, -2*ypz + dxy,
		  dxyz, 2*xpz - 4*dyz,
	      dxyz};

	uint32_t k(0), j;
	for (uint32_t i = 0u; i < 8; ++i) { // 8 = npe
		for (j = 0; j < i; ++j) M(i, j) = M(j, i);
		for (j = i; j < 8; ++j) M(i, j) = V_[k++];
	}
	/*
	V_ = {         dxyz, 2*xpz - 4*dyz,  -2*ypz + dxy, 2*xpy - 4*dxz, -4*dxy + 2*ypz, -2*xpy + dxz, -xpypz, -2*xpz + dyz,
          2*xpz - 4*dyz,          dxyz, 2*xpy - 4*dxz,  dxy - 2*ypz,  -2*xpz + dxz, -4*dxy + 2*ypz,  dyz - 2*xpz, -xpypz,
		   -2*ypz + dxy, 2*xpy - 4*dxz,          dxyz, 2*xpz - 4*dyz,        -xpypz,   dyz - 2*xpz, -4*dxy + 2*ypz, -2*xpy + dxz,
		2 * xpy - 4 * dxz,  dxy - 2 * ypz, 2 * xpz - 4 * dyz, dxyz,   dyz - 2*xpz,       -xpypz,  -2*xpy + dxz, -4*dxy + 2*ypz,
		-4 * dxy + 2 * ypz,   -2 * xpz + dxz,        -xpypz,   dyz - 2 * xpz, dxyz, 2*xpz - 4*dyz, -2*ypz + dxy, 2*xpy - 4*dxz,
		-2 * xpy + dxz,  -4 * dxy + 2 * ypz,   dyz - 2 * xpz,       -xpypz, 2 * xpz - 4 * dyz,  dxyz, 2*xpy - 4*dxz, -2*ypz + dxy,
		-xpypz,   dyz - 2 * xpz,-4 * dxy + 2 * ypz,   -2 * xpy + dxz,  -2 * ypz + dxy, 2 * xpy - 4 * dxz, dxyz, 2*xpz - 4*dyz,
		-2 * xpz + dyz,  -xpypz, -2 * xpy + dxz, -4 * dxy + 2 * ypz,  2 * xpy - 4 * dxz,  -2 * ypz + dxy, 2 * xpz - 4 * dyz, dxyz};
	*/
}



} // end namespace csmp
