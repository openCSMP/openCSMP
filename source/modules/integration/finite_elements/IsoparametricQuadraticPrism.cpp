#include "IsoparametricQuadraticPrism.h"
#include <fstream>
#include <cstring>
#include "Exception.h"
#include <climits>

using namespace std;

namespace csmp {

IsoparametricQuadraticPrism::IsoparametricQuadraticPrism()
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_QUADRATIC_PRISM18, true, true, 2U ),
      DN(3,18),
      NXYZ(18,3),
      IP(6,3)
 {
    //AAM, 15.08.02
    dim = 3;
    itp = 2;
    npf = 9;
    npe = 18;
    fpe = 5;
    spe = 9;
    epe = 5;
    // nne - typical number of finite elements, sharing each node
    nne = 6;
    cne = 0;
    // Integration points
    gpe = 6;

    // What is M?
    M.Resize(npe,npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_QUADRATIC_PRISM18);

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
    // r-s-t/x-y-z
    NXYZ(0,0) = 0.0, NXYZ(0,1) = 0.0, NXYZ(0,2) = -1.0;
    NXYZ(1,0) = 1.0, NXYZ(1,1) = 0.0, NXYZ(1,2) = -1.0;
    NXYZ(2,0) = 0.0, NXYZ(2,1) = 1.0, NXYZ(2,2) = -1.0;
    NXYZ(3,0) = 0.0, NXYZ(3,1) = 0.0, NXYZ(3,2) =  1.0;
    NXYZ(4,0) = 1.0, NXYZ(4,1) = 0.0, NXYZ(4,2) =  1.0;
    NXYZ(5,0) = 0.0, NXYZ(5,1) = 1.0, NXYZ(5,2) =  1.0;

    NXYZ(6,0) = 0.5;NXYZ(6,1) = 0.0;NXYZ(6,2) = -1.0; //7
    NXYZ(7,0) = 0.5;NXYZ(7,1) = 0.5;NXYZ(7,2) = -1.0; //8
    NXYZ(8,0) = 0.0;NXYZ(8,1) = 0.5;NXYZ(8,2) = -1.0; //9

    NXYZ(9,0) = 0.0;NXYZ(9,1) = 0.0;NXYZ(9,2) =  0.0; //10
    NXYZ(10,0)= 1.0;NXYZ(10,1)= 0.0;NXYZ(10,2) = 0.0; //11
    NXYZ(11,0)= 0.0;NXYZ(11,1)= 1.0;NXYZ(11,2) = 0.0; //12

    NXYZ(12,0)= 0.5;NXYZ(12,1)= 0.0;NXYZ(12,2) = 1.0; //13
    NXYZ(13,0)= 0.5;NXYZ(13,1)= 0.5;NXYZ(13,2) = 1.0; //14
    NXYZ(14,0)= 0.0;NXYZ(14,1)= 0.5;NXYZ(14,2) = 1.0; //15

    NXYZ(15,0)= 0.5;NXYZ(15,1)= 0.0;NXYZ(15,2) = 0.0; //16
    NXYZ(16,0)= 0.5;NXYZ(16,1)= 0.5;NXYZ(16,2) = 0.0; //17
    NXYZ(17,0)= 0.0;NXYZ(17,1)= 0.5;NXYZ(17,2) = 0.0; //18

    W.resize( gpe );
    // Numeric integration 2nd order, number of point m=6
    GenerateIntegrationPoints(IP,W);
 }



IsoparametricQuadraticPrism::~IsoparametricQuadraticPrism()
 {
 }








/**

Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local corner node ID numbers
for the element.

*/
void
IsoparametricQuadraticPrism::CornerNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(6);
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
IsoparametricQuadraticPrism::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0]  = 0;
    ids[1]  = 6;
    ids[2]  = 1;
    ids[3]  = 7;
    ids[4]  = 2;
    ids[5]  = 8;
    ids[6]  = 9;
    ids[7]  =15;
    ids[8]  =10;
    ids[9]  =16;
    ids[10] =11;
    ids[11] =17;
    ids[12] =3;
    ids[13] =12;
    ids[14] =4;
    ids[15] =13;
    ids[16] =5;
    ids[17] =14;
 }



void
IsoparametricQuadraticPrism::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
 {
    snids.resize(3);
    if ( segm_id == 0 ) {
         snids[0] = 0;
         snids[1] = 1;
         snids[2] = 6;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 2;
         snids[1] = 1;
         snids[2] = 7;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 2;
         snids[1] = 0;
         snids[2] = 8;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 4;
         snids[1] = 3;
         snids[2] = 12;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 5;
         snids[1] = 4;
         snids[2] = 13;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 5;
         snids[1] = 3;
         snids[2] = 14;
      }
    else if ( segm_id == 6 ) {
         snids[0] = 3;
         snids[1] = 0;
         snids[2] = 9;
      }
    else if ( segm_id == 7 ) {
         snids[0] = 4;
         snids[1] = 1;
         snids[2] = 10;
      }
    else if ( segm_id == 8 ) {
         snids[0] = 5;
         snids[1] = 2;
         snids[2] = 11;
      }
    else
    std::cerr <<"\nIsoparametricQuadraticPrism::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment




/**

For this element, the faces are numbered such that the lower left closest is 1 ->
*/
void
IsoparametricQuadraticPrism::NodesOfFace( uint32_t face_id, std::vector<uint32_t>& fnids ) const
 {
    fnids.resize(9);

    // triangular face at base
    if ( face_id == 0 )
      {
        fnids.resize(6);
         // corner nodes
         fnids[0] = 0;
         fnids[1] = 2;
         fnids[2] = 1;
         // midside nodes
         fnids[3] = 8;
         fnids[4] = 7;
         fnids[5] = 6;
     }
    else if ( face_id == 1 )
      {
      fnids.resize(8);
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 4;
         fnids[3] = 3;
         //
         fnids[4] = 6;
         fnids[5] = 10;
         fnids[6] = 12;
         fnids[7] = 9;
      }
    else if ( face_id == 2 )
      {
      fnids.resize(8);
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 5;
         fnids[3] = 4;
         //
         fnids[4] = 7;
         fnids[5] = 11;
         fnids[6] = 13;
         fnids[7] = 10;
      }
    else if ( face_id == 3 )
      {
      fnids.resize(8);
         fnids[0] = 0;
         fnids[1] = 3;
         fnids[2] = 5;
         fnids[3] = 2;
         //
         fnids[4] = 9;
         fnids[5] = 14;
         fnids[6] = 11;
         fnids[7] = 8;
      }
    // triangular face at top
    else if ( face_id == 4 )
      {
      fnids.resize(6);
         fnids[0] = 3;
         fnids[1] = 4;
         fnids[2] = 5;
         fnids[3] = 12;
         fnids[4] = 13;
         fnids[5] = 14;
      }
    else
    std::cerr <<"\nIsoparametricQuadraticPrism::NodesOfFace: Invalid Face ID requested: "<< face_id <<std::endl;
 }



vector<uint32_t>  IsoparametricQuadraticPrism::CornerNodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{0,2,1};
        case 1: return vector<uint32_t>{0,1,4,3};
        case 2: return vector<uint32_t>{1,2,5,4};
        case 3: return vector<uint32_t>{0,3,5,2};
        case 4: return vector<uint32_t>{3,4,5};
      }
    cerr <<"\nIsoparametricQuadraticPrism::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }



vector<uint32_t>  IsoparametricQuadraticPrism::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        // local corner node numbers are returned in ascending order
        case 0: return vector<uint32_t>{1,2,3};
        case 1: return vector<uint32_t>{0,2,4};
        case 2: return vector<uint32_t>{0,1,5};
        case 3: return vector<uint32_t>{0,4,5};
        case 4: return vector<uint32_t>{1,3,5};
        case 5: return vector<uint32_t>{2,3,4};
        default:
          cerr <<"\nIsoparametricQuadraticPrism::NodesConnectedTo: node "<< node_id <<" does not exist.";
      }
    return vector<uint32_t>{};
  }



CSMP_FEM_TYPE
IsoparametricQuadraticPrism::ElementTypeOfFace( uint32_t )  const
 {
    return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;
 }

double IsoparametricQuadraticPrism::WeightAtIntegrationPoint( uint32_t i ) const { return W[i]; }



void
IsoparametricQuadraticPrism::GenerateIntegrationPoints(
                                                    DenseMatrix<DM_MIN> &Ip,
                                                    std::vector<double>& We)
{
    const uint32_t numberOfTriaIntegrationPoints=3; // order 2
    const uint32_t numberOfDimensionsInPlane=2;
    DenseMatrix<DM_MIN> IPTRIA(numberOfTriaIntegrationPoints,numberOfDimensionsInPlane);
    const double constA=0.577350269189626, constT1=0.66666666666667, constT2=0.16666666666667;

    IPTRIA(0,0)=constT1; IPTRIA(0,1)=constT2;
    IPTRIA(1,0)=constT2; IPTRIA(1,1)=constT1;
    IPTRIA(2,0)=constT2; IPTRIA(2,1)=constT2;

    vector<double> WTRIA(numberOfTriaIntegrationPoints);
    // Weights for the base integration
    WTRIA[0]=constT2;
    WTRIA[1]=constT2;
    WTRIA[2]=constT2;

    const uint32_t numberOfLineIntegrationPoints=2;
    vector<double> IPLINE(numberOfLineIntegrationPoints),WLINE(numberOfLineIntegrationPoints);

    IPLINE[0]=-constA;  IPLINE[1]=constA;
    WLINE [0]= 1.0;     WLINE [1]=1.0;

    uint32_t i=0;
    for(uint32_t j=0;j<numberOfTriaIntegrationPoints;j++)
    {
    for(uint32_t k=0;k<numberOfLineIntegrationPoints;k++)
        {
            Ip(i,0)=IPTRIA(j,0);
            Ip(i,1)=IPTRIA(j,1);
            Ip(i,2)=IPLINE[k];
            We[i++]=WTRIA[j]*WLINE[k];
        }
    }
}

void
IsoparametricQuadraticPrism::OutputNodeDataToVTK( const char* file_name,
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
           cerr <<"\nIsoparametricQuadraticPrism::OutputNodeDataToVTK ";
           cerr <<"Output file could not be opened."<< endl;
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
     for ( uint32_t i=0; i<npe; i++ )
       {
          for ( uint32_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point)) -corect
     // Cells
     // -----------------------------------------------------
     ofs <<"CELLS "<< 8 <<" "<< 56 << endl;
     // (1+4)X8
     // the 4 corner quads +4 midside nodes in form of prism
     ofs << 6 <<" 0 6 8 9 15 17" << endl;
     ofs << 6 <<" 6 1 7 15 10 16" << endl;
     ofs << 6 <<" 6 7 8 15 16 17" << endl;
     ofs << 6 <<" 8 7 2 17 16 11" << endl;

     ofs << 6 <<" 9 15 17 3 12 14" << endl;
     ofs << 6 <<" 15 10 16 12 4 13" << endl;
     ofs << 6 <<" 15 16 17 12 13 14" << endl;
     ofs << 6 <<" 17 16 11 14 13 5" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 8 << endl;
     ofs << 13 << endl; // VTK_WEDGE
     ofs << 13 << endl; // VTK_WEDGE
     ofs << 13 << endl; // VTK_WEDGE
     ofs << 13 << endl; // VTK_WEDGE
     ofs << 13 << endl; // VTK_WEDGE
     ofs << 13 << endl; // VTK_WEDGE
     ofs << 13 << endl; // VTK_WEDGE
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
           for ( uint32_t i=0; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 10
          for ( uint32_t i=0; i<DATA.Cols(); i++ ) {
               for ( uint32_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricQuadraticPrism::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

 } // end OutputNodeDataToVTK






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
void IsoparametricQuadraticPrism::Nrst(
                    double r, // L2
                    double s, // L3
                    double t,
                    std::vector<double>& N ) const
{
   N.resize(npe);
   const double L1=1.-r-s;
   N[0] =  0.5*L1*(2*L1-1.)*t*(t-1.);
   N[1] =  0.5*r*(2*r-1.)*t*(t-1.);
   N[2] =  0.5*s*(2*s-1.)*t*(t-1.);
   N[3] =  0.5*L1*(2*L1-1.)*t*(t+1.);
   N[4] =  0.5*r*(2*r-1.)*t*(t+1.);
   N[5] =  0.5*s*(2*s-1.)*t*(t+1.);
   N[6] =  2*L1*r*t*(t-1.);
   N[7] =  2*r*s*t*(t-1.);

   N[8] =  2*s*L1*t*(t-1.);
   N[9] =  L1*(2*L1-1.)*(1.-t*t);
   N[10]=  r*(2*r-1.)*(1.-t*t);
   N[11]=  s*(2*s-1.)*(1.-t*t);
   N[12]=  2*L1*r*t*(t+1.);
   N[13]=  2*r*s*t*(t+1.);
   N[14]=  2*s*L1*t*(t+1.);

   N[15]=  4*L1*r*(1.-t*t);
   N[16]=  4*r*s*(1.-t*t);
   N[17]=  4*s*L1*(1.-t*t);
}

void IsoparametricQuadraticPrism::Nrst(
                    double r, // L2
                    double s, // L3
                    double t,
                    double* N ) const
{
   const double L1=1.-r-s;
   N[0] =  0.5*L1*(2*L1-1.)*t*(t-1.);
   N[1] =  0.5*r*(2*r-1.)*t*(t-1.);
   N[2] =  0.5*s*(2*s-1.)*t*(t-1.);
   N[3] =  0.5*L1*(2*L1-1.)*t*(t+1.);
   N[4] =  0.5*r*(2*r-1.)*t*(t+1.);
   N[5] =  0.5*s*(2*s-1.)*t*(t+1.);
   N[6] =  2*L1*r*t*(t-1.);
   N[7] =  2*r*s*t*(t-1.);

   N[8] =  2*s*L1*t*(t-1.);
   N[9] =  L1*(2*L1-1.)*(1.-t*t);
   N[10]=  r*(2*r-1.)*(1.-t*t);
   N[11]=  s*(2*s-1.)*(1.-t*t);
   N[12]=  2*L1*r*t*(t+1.);
   N[13]=  2*r*s*t*(t+1.);
   N[14]=  2*s*L1*t*(t+1.);

   N[15]=  4*L1*r*(1.-t*t);
   N[16]=  4*r*s*(1.-t*t);
   N[17]=  4*s*L1*(1.-t*t);
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
void IsoparametricQuadraticPrism::dNr (
                double r,
                double s,
                double t,
                std::vector<double>& DNR ) const
{
   const double L1=1-r-s;
   DNR.resize(npe);
   DNR[0] =-0.5*(1.-2.*r-2.*s)*t*(t-1.)-1.0*L1*t*(t-1.);
   DNR[1] = 0.5*(2.*r-1.)*t*(t-1.)+1.0*r*t*(t-1.);
   DNR[2] = 0.0;

   DNR[3] = -0.5*(1.-2.*r-2.*s)*t*(t+1.)-1.0*L1*t*(t+1.);
   DNR[4] = 0.5*(2.*r-1.)*t*(t+1.)+1.0*r*t*(t+1.);
   DNR[5] = 0.0;

   DNR[6] = -2.0*r*t*(t-1.)+2.0*L1*t*(t-1.);
   DNR[7] = 2.0*s*t*(t-1.);

   DNR[8] = -2.0*s*t*(t-1.);
   DNR[9] = -(1.0-2.0*r-2.0*s)*(1.-t*t)-2.0*L1*(1.-t*t);
   DNR[10]= (2.0*r-1.)*(1.-t*t)+2.0*r*(1.-t*t);

   DNR[11]= 0.0;

   DNR[12]=-2.0*r*t*(t+1.)+2.0*L1*t*(t+1.);
   DNR[13]= 2.0*s*t*(t+1.);
   DNR[14]=  -2.0*s*t*(t+1.);
   DNR[15]=-4.0*r*(1.-t*t)+4.0*L1*(1.-t*t);

   DNR[16]=  4.0*s*(1.-t*t);
   DNR[17]=  -4.0*s*(1.-t*t);
}


void IsoparametricQuadraticPrism::dNs(
                double r,
                double s,
                double t,
                std::vector<double>& DNS ) const
{
   const double L1=1.-r-s;
   DNS.resize(npe);
   DNS[0] =-.5*(1.-2.*r-2.*s)*t*(t-1.)-1.0*L1*t*(t-1.);
   DNS[1] =0.0;
   DNS[2] =.5*(2.*s-1.)*t*(t-1.)+1.0*s*t*(t-1.);

   DNS[3] =-.5*(1.-2.*r-2.*s)*t*(t+1.)-1.0*L1*t*(t+1.);
   DNS[4] =0.0;
   DNS[5] =.5*(2.*s-1.)*t*(t+1.)+1.0*s*t*(t+1.);

   DNS[6] = -2.0*r*t*(t-1.);
   DNS[7] =2.0*r*t*(t-1.);

   DNS[8] = 2.0*L1*t*(t-1.)-2.0*s*t*(t-1.);
   DNS[9] =-(1.0-2.0*r-2.0*s)*(1.-t*t)-2.0*L1*(1.-t*t);
   DNS[10]= 0.0;
   DNS[11] =(2.0*s-1.)*(1.-t*t)+2.0*s*(1.-t*t);

   DNS[12] = -2.0*r*t*(t+1.);
   DNS[13] =2.0*r*t*(t+1.);
   DNS[14] = -2.0*s*t*(t+1.)+2.0*L1*t*(t+1.);

   DNS[15] = -4.0*r*(1.-t*t);
   DNS[16] =  4.0*r*(1.-t*t);
   DNS[17] =  4.0*L1*(1.-t*t)-4.0*s*(1.-t*t);

}


void IsoparametricQuadraticPrism::dNt(
                double r,
                double s,
                double t,
                std::vector<double>& DNT ) const
{

   const double L1=1.-r-s;
   DNT.resize(npe);
   DNT[0] =.5*L1*(1.-2.*r-2.*s)*(t-1.)+.5*L1*(1.-2.*r-2.*s)*t;
   DNT[1] =.5*r*(2.*r-1.)*(t-1.)+.5*r*(2.*r-1.)*t;
   DNT[2] =.5*s*(2.*s-1.)*(t-1.)+.5*s*(2.*s-1.)*t;

   DNT[3] =.5*L1*(1.-2.*r-2.*s)*(t+1.)+.5*L1*(1.-2.*r-2.*s)*t;
   DNT[4] =.5*r*(2.*r-1.)*(t+1.)+.5*r*(2.*r-1.)*t;
   DNT[5] =.5*s*(2.*s-1.)*(t+1.)+.5*s*(2.*s-1.)*t;

   DNT[6] =2.0*L1*r*(t-1.)+2.0*L1*r*t;
   DNT[7] =2.0*r*s*(t-1.)+2.0*r*s*t;

   DNT[8] =2.0*s*L1*(t-1.)+2.0*s*L1*t;
   DNT[9] =-2*L1*(1.0-2.0*r-2.0*s)*t;
   DNT[10]=-2*r*(2.0*r-1.)*t;

   DNT[11]= -2*s*(2.0*s-1.)*t;

   DNT[12] =  2.0*L1*r*(t+1.)+2.0*L1*r*t;
   DNT[13] =  2.0*r*s*(t+1.)+2.0*r*s*t;
   DNT[14] =  2.0*L1*s*(t+1.)+2.0*s*L1*t;
   DNT[15] =  -8.0*L1*r*t;

   DNT[16] = -8.0*r*s*t;
   DNT[17] = -8.0*s*L1*t;
}



/** Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN18 The interpolation-function derivative matrix is returned into the
second method argument.

*/
void
IsoparametricQuadraticPrism::dN( DenseMatrix<DM_MIN>& DN18 )
{
    DN18.Resize(dim,npe);
    M.Resize(dim,1);

    //DenseMatrix<DM_MIN> DERIVS(dim,npe);


     // Jacobian transformation to global coordinate system
     for ( uint32_t i=0; i<npe; i++ )
       {
          // here the global coordinates come in
          dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
          dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
          dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

         //cout<<" Node="<<i+1<<" ("<<XY(i,0)<<","<<XY(i,1)<<","<<XY(i,2)<<")"<<endl;
         //cout<<" N-point point: ("<<NXYZ(i,0)<<","<<NXYZ(i,1)<<","<< NXYZ(i,2)<<")"<<endl;
         //cout<<" Ders : ("<<DNR[i]<<","<<DNS[i]<<","<<DNT[i]<<")"<<endl;

         //DERIVS(0,i)=DNR[i];DERIVS(1,i)=DNS[i];DERIVS(2,i)=DNT[i];

          // compute Jacobian matrix, its determinant and inverse
        Jacobian( DNR, DNS, DNT );

        //cout<<" Jacobian Matrix: ";
        //JAC.Out();

        JacobianInverse();
        ///////// Debug Printout////////////////////////////////
        //cout<<" Jacobian Inverse Matrix: ";
        //JINV.Out();
        //getchar();
        /////////////////////////////////////////////////////////
        M(0,0)=DNR[i]; M(1,0)=DNS[i]; M(2,0)=DNT[i];
        // 3x3 * 3x1 = 3x1 gives the global DN entries
        JINV *= M;
        DN18(0,i) = JINV(0,0);
        DN18(1,i) = JINV(1,0);
        DN18(2,i) = JINV(2,0);
        JINV.Resize(dim,dim);
       }

       //OutputNodeDataToVTK( "etest1", "derivatives",DERIVS);

} // end dN

/** Computes derivatives of shape functions at a point XY given in global coordinates
within the element (also given in global coordinates).

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at the point XY are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument.

*/
double
IsoparametricQuadraticPrism::dN_At( DenseMatrix<DM_MIN>& DN2,
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

    DenseMatrix<DM_MIN> DN(dim,npe);

    dN(DN);

    DN2*=DN;

    /////////////////////////////// Debug printout ///////////////////////////////////////////////
    //cout<<" IsoparametricQuadraticPrism::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
    //rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
    //cout<<" IsoparametricQuadraticPrism::dN  Matrix DN2: "<<endl;
    //DN2.Out();
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

    return detJ;

 }

/** Outputs the shape functions at the point 'xyz' which must lie
within the Hexahedron.
 */
void
IsoparametricQuadraticPrism::N(
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
IsoparametricQuadraticPrism::dN_AtNode( DenseMatrix<DM_MIN>& B, uint32_t nd )
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
    B.Resize(dim,npe);
    for(uint32_t inode=0;inode<npe; inode++)
     {B(0,inode) = DNR[inode], B(1,inode) = DNS[inode], B(2,inode) = DNT[inode]; }

    B = JINV * B;

    return detJ;
 }

/** Projection function from rst->xyz
*/
void
IsoparametricQuadraticPrism::ParametricToPhysical(std::vector<double> &rst, std::vector<double>& xyz)
{
Nrst(rst[0],rst[1],rst[2], NRST );

for(uint32_t i=0; i<dim; i++) xyz[i]=0.0;

for(uint32_t i=0; i<npe; i++)
    {
    xyz[0]+=XY(i,0)*NRST[i];
    xyz[1]+=XY(i,1)*NRST[i];
    xyz[2]+=XY(i,2)*NRST[i];
    }
}

/** Projection function from xyz->rst
*/
void
IsoparametricQuadraticPrism::PhysicalToParametric(
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

        vector<double> rstHatK_PlusOne(dim);
        double minDistanceFromGivenPoint;

        const double constantMu               = 1.0;
        const uint32_t numberOfFirstIterrations   = 5;
        const uint32_t maxNumberOfIterrations     = 20;
        const uint32_t numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double incrementR = 1.0/(numberOfFirstIterrations-1);
        const double incrementS = 1.0/(numberOfFirstIterrations-1);
        const double incrementT = 2.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = 0.0;
        for(uint32_t i=0;i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = 0.0;
            for(uint32_t j=0;j<numberOfFirstIterrations;j++)
            {
                rstHatK_PlusOne[2] = -1.0;
                for(uint32_t k=0;k<numberOfFirstIterrations;k++)
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
                        for(uint32_t i=0; i<dim; i++)
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

        uint32_t iteration = 1;

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
                //const double detJ = JacobianDeterminant();
                //if( detJ > 0.0 )
                JacobianInverse();
            }

            // Newton iteration
            rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) +JINV(2,0)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) +JINV(2,1)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[2] = rstHatK[2] - constantMu*(JINV(0,2)*(outxyz[0]-xyz[0]) + JINV(1,2)*(outxyz[1]-xyz[1]) +JINV(2,2)*(outxyz[2]-xyz[2]));

            for(uint32_t i=0; i<dim; i++)
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
            cout<<" IsoparametricQuadraticPrism::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            cout<<" IsoparametricQuadraticPrism::PhysicalToParametric: Real Point:\t"
                <<"x = "<<xyz[0]<<" ;\t"
                <<"y = "<<xyz[1]<<" ;\t"
                <<"z = "<<xyz[2]<<"\n";
            cout<<" IsoparametricQuadraticPrism::PhysicalToParametric: Found Point:\t"
               <<"x = "<<outxyz[0]<<" ;\t"
               <<"y = "<<outxyz[1]<<" ;\t"
               <<"z = "<<outxyz[2]<<"\n";
            cout<<" IsoparametricQuadraticPrism::PhysicalToParametric: Distance:\t"
               <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
               <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
               <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            cout<<" IsoparametricQuadraticPrism::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<" ;\t"
               <<"t = "<<rstHatK[2]<<"\n";

            std::vector<double> N(npe,0.0);
            Nrst(rstHatK[0], rstHatK[1], rstHatK[2], N );
            cout<<" IsoparametricQuadraticPrism::PhysicalToParametric: Shape functions:\t"
                <<"N[0] = "<<N[0]<<" ;\t"
                <<"N[1] = "<<N[1]<<" ;\t"
                <<"N[2] = "<<N[2]<<" ;\t"
                <<"N[3] = "<<N[3]<<" ;\t"
                <<"N[4] = "<<N[4]<<" ;\t"
                <<"N[5] = "<<N[5]<<" ;\t"
                <<"N[6] = "<<N[6]<<" ;\t"
                <<"N[7] = "<<N[7]<<" ;\t"
                <<"N[8] = "<<N[8]<<" ;\t"
                <<"N[9] = "<<N[9]<<" ;\t"
                <<"N[10] = "<<N[10]<<" ;\t"
                <<"N[11] = "<<N[11]<<" ;\t"
                <<"N[12] = "<<N[12]<<" ;\t"
                <<"N[13] = "<<N[13]<<" ;\t"
                <<"N[14] = "<<N[14]<<" ;\t"
                <<"N[15] = "<<N[15]<<" ;\t"
                <<"N[16] = "<<N[16]<<" ;\t"
                <<"N[17] = "<<N[17]<<"\n";

            csmp::Exception( WARNING, "IsoparametricQuadraticPrism::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for(uint32_t i=0; i<dim; i++)
        rSt[i] = rstHatK[i];

}


 /** Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double  IsoparametricQuadraticPrism::AspectRatio()
{
   vector<double> vec(spe);

   EdgeLengths( vec );

   double seg_max(vec[0]), seg_min(vec[0]);

   // find largest segment
   for ( uint32_t i=1; i<spe; i++ ) {
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
IsoparametricQuadraticPrism::EdgeLengths( std::vector<double>& len )
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
    sum     = (XY(2,0)-XY(0,0)) * (XY(2,0)-XY(0,0));
    sum    += (XY(2,1)-XY(0,1)) * (XY(2,1)-XY(0,1));
    sum    += (XY(2,2)-XY(0,2)) * (XY(2,2)-XY(0,2));
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
    sum     = (XY(5,0)-XY(3,0)) * (XY(5,0)-XY(3,0));
    sum    += (XY(5,1)-XY(3,1)) * (XY(5,1)-XY(3,1));
    sum    += (XY(5,2)-XY(3,2)) * (XY(5,2)-XY(3,2));
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

/*

Computes Integral N dA = sum_1...n Wi Ji Ni. This formulation also
gives a meaningful area if the element boundaries are curved. For a
description of numerical integration  of isoparametric quadratic
triangular elements, refer to Cook et al. (1989) 3rd Ed., p. 183.

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The area (m2) of the finite element.
*/
/*
double
IsoparametricQuadraticPrism::VolumeThroughHex()
{
    // determinant
    //uint32_t  i=0;

  const uint32_t    dim1(3);

  QuadraticHexahedron     quadratic_hexahedron;

  Node<dim1>  n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19,
  n20, n21, n22, n23, n24, n25, n26;

  n0.Idx( 1 ),	n0.x(XYZ(0,0)),	n0.y(XYZ(0,1)),	n0.z(XYZ(0,2));
  n1.Idx( 2 ),	n1.x(XYZ(1,0)), n1.y(XYZ(1,1)),	n1.z(XYZ(1,2));
  n2.Idx( 3 ),	n2.x(XYZ(2,0)),	n2.y(XYZ(2,1)),	n2.z(XYZ(2,2));
  n3.Idx( 4 ),	n3.x(XYZ(2,0)),	n3.y(XYZ(2,1)),	n3.z(XYZ(2,2));

  n4.Idx( 5 ),	n4.x(XYZ(3,0)),	n4.y(XYZ(3,1)),	n4.z(XYZ(3,2));
  n5.Idx( 6 ),	n5.x(XYZ(4,0)),	n5.y(XYZ(4,1)),	n5.z(XYZ(4,2));
  n6.Idx( 7 ),	n6.x(XYZ(5,0)),	n6.y(XYZ(5,1)),	n6.z(XYZ(5,2));
  n7.Idx( 8 ),	n7.x(XYZ(5,0)),	n7.y(XYZ(5,1)),	n7.z(XYZ(5,2));

  n8.Idx( 9 ),	n8.x(XYZ(6,0)),	n8.y(XYZ(6,1)),	n8.z(XYZ(6,2));
  n9.Idx( 10 ),	n9.x(XYZ(7,0)),	n9.y(XYZ(7,1)),	n9.z(XYZ(7,2));
  n10.Idx( 11 ),	n10.x(XYZ(2,0)),n10.y(XYZ(2,1)),n10.z(XYZ(2,2));
  n11.Idx( 12 ),	n11.x(XYZ(8,0)),n11.y(XYZ(8,1)),n11.z(XYZ(8,2));
  n12.Idx( 13 ),	n12.x(XYZ(9,0)),n12.y(XYZ(9,1)),n12.z(XYZ(9,2));
  n13.Idx( 14 ),	n13.x(XYZ(10,0)),n13.y(XYZ(10,1)),n13.z(XYZ(10,2));
  n14.Idx( 15 ),	n14.x(XYZ(11,0)),n14.y(XYZ(11,1)),n14.z(XYZ(11,2));
  n15.Idx( 16 ),	n15.x(XYZ(11,0)),n15.y(XYZ(11,1)),n15.z(XYZ(11,2));
  n16.Idx( 17 ),	n16.x(XYZ(12,0)),n16.y(XYZ(12,1)),n16.z(XYZ(12,2));
  n17.Idx( 18 ),	n17.x(XYZ(13,0)),n17.y(XYZ(13,1)),n17.z(XYZ(13,2));
  n18.Idx( 19 ),	n18.x(XYZ(5,0)),n18.y(XYZ(5,1)),n18.z(XYZ(5,2));
  n19.Idx( 20 ),	n19.x(XYZ(14,0)),n19.y(XYZ(14,1)),n19.z(XYZ(14,2));

  n20.Idx( 21 ),	n20.x((XYZ(0,0)+XYZ(1,0)+XYZ(2,0))/3.0), n20.y((XYZ(0,1)+XYZ(1,1)+XYZ(2,1))/3.0), n20.z((XYZ(0,2)+XYZ(1,2)+XYZ(2,2))/3.0);
  n21.Idx( 22 ),	n21.x(XYZ(15,0)), 		n21.y(XYZ(15,1)),	n21.z(XYZ(15,2));
  n22.Idx( 23 ),	n22.x(XYZ(16,0)), 		n22.y(XYZ(16,1)),	n22.z(XYZ(16,2));
  n23.Idx( 24 ),	n23.x(XYZ(11,0)),	 	n23.y((XYZ(11,1)),	n23.z((XYZ(11,2));
  n24.Idx( 25 ),	n24.x(XYZ(17,0)),  		n24.y((XYZ(17,1)),	n24.z((XYZ(17,2));

  n25.Idx( 26 ),	n25.x((XYZ(3,0)+XYZ(4,0)+XYZ(5,0))/3.0)), n25.y((XYZ(3,1)+XYZ(4,1)+XYZ(5,1))/3.0)),
  n25.z((XYZ(3,2)+XYZ(4,2)+XYZ(5,2))/3.0));

  n26.Idx( 27 ),	n26.x((XYZ(9,0)+XYZ(10,0)+XYZ(11,0))/3.0), n26.y((XYZ(9,1)+XYZ(10,1)+XYZ(11,1))/3.0),
  n26.z((XYZ(9,2)+XYZ(10,2)+XYZ(11,2))/3.0);

   // element 1
  Element<dim1>  element1( & quadratic_hexahedron);

  element1.Idx( 1 );
  element1.ConnectToNodes ( 0, n0 );
  element1.ConnectToNodes ( 1, n1 );
  element1.ConnectToNodes ( 2, n2 );
  element1.ConnectToNodes ( 3, n3 );
  element1.ConnectToNodes ( 4, n4 );
  element1.ConnectToNodes ( 5, n5 );
  element1.ConnectToNodes ( 6, n6 );
  element1.ConnectToNodes ( 7, n7 );
  element1.ConnectToNodes ( 8, n8 );
  element1.ConnectToNodes ( 9, n9 );
  element1.ConnectToNodes ( 10, n10 );
  element1.ConnectToNodes ( 11, n11 );
  element1.ConnectToNodes ( 12, n12 );
  element1.ConnectToNodes ( 13, n13 );
  element1.ConnectToNodes ( 14, n14 );
  element1.ConnectToNodes ( 15, n15 );
  element1.ConnectToNodes ( 16, n16 );
  element1.ConnectToNodes ( 17, n17 );
  element1.ConnectToNodes ( 18, n18 );
  element1.ConnectToNodes ( 19, n19 );
  element1.ConnectToNodes ( 20, n20 );
  element1.ConnectToNodes ( 21, n21 );
  element1.ConnectToNodes ( 22, n22 );
  element1.ConnectToNodes ( 23, n23 );
  element1.ConnectToNodes ( 24, n24 );
  element1.ConnectToNodes ( 25, n25 );
  element1.ConnectToNodes ( 26, n26 );

  double   areaQ=element1.Volume();
  cout<<" Quadratic Volume of Prism="<<areaQ<<endl;

  Node<dim1>  n0, n1, n2, n3, n4, n5, n6, n7, n8;
  IsoparametricLinearHexahedron     linear_hexahedron(8);
  n0.Idx( 1 ),	n0.x(XYZ(0,0)),	n0.y(XYZ(0,1)),	n0.z(XYZ(0,2));
  n1.Idx( 2 ),	n1.x(XYZ(1,0)), n1.y(XYZ(1,1)),	n1.z(XYZ(1,2));
  n2.Idx( 3 ),	n2.x(XYZ(2,0)),	n2.y(XYZ(2,1)),	n2.z(XYZ(2,2));
  n3.Idx( 4 ),	n3.x(XYZ(2,0)),	n3.y(XYZ(2,1)),	n3.z(XYZ(2,2));

  n4.Idx( 5 ),	n4.x(XYZ(3,0)),	n4.y(XYZ(3,1)),	n4.z(XYZ(3,2));
  n5.Idx( 6 ),	n5.x(XYZ(4,0)),	n5.y(XYZ(4,1)),	n5.z(XYZ(4,2));
  n6.Idx( 7 ),	n6.x(XYZ(5,0)),	n6.y(XYZ(5,1)),	n6.z(XYZ(5,2));
  n7.Idx( 8 ),	n7.x(XYZ(5,0)),	n7.y(XYZ(5,1)),	n7.z(XYZ(5,2));
   // element 2
  Element<dim1>  element2( & linear_hexahedron);

  element2.Idx( 2 );
  element2.ConnectToNodes ( 0, n0 );
  element2.ConnectToNodes ( 1, n1 );
  element2.ConnectToNodes ( 2, n2 );
  element2.ConnectToNodes ( 3, n3 );
  element2.ConnectToNodes ( 4, n4 );
  element2.ConnectToNodes ( 5, n5 );
  element2.ConnectToNodes ( 6, n6 );
  element2.ConnectToNodes ( 7, n7 );

  double   areaL=element2.Volume();

  DenseMatrix<DM_MIN>  IPPHYS(gpe,dim);
  linear_hexahedron.IntegrationPointsFromParToPhys(IPPHYS);

  //IP=IPPHYS;
  //for(uint32_t i=0; i<gpe; i++) W[i]=2.0;

 return areaL;
}
*/


/**

Computes Integral N dA = sum_1...n Wi Ji Ni. This formulation also
gives a meaningful area if the element boundaries are curved.
For a  description of numerical integration  of isoparametric quadratic
triangular elements, refer to Cook et al. (1989) 3rd Ed., p. 183.

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The area (m2) of the finite element.
*/
double
IsoparametricQuadraticPrism::Volume()
{
    double   area{0.}; // determinant

    // numerical integration:
    // looping over the 4 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( auto i=0; i<gpe; i++ )
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
    //cout<<" IsoparametricQuadraticPrism::Volume() Volume of Prism="<<area<<", Versus through Hex="<<VolumeThroughHex()<<endl;

    return area;
}


/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

@param N The interpolation function values are returned into the third argument.

*/
inline void
IsoparametricQuadraticPrism::N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N )
 {
    assert( ip < gpe );
     // local interpolation function values
    //if(W[0]==0.0 && W[1]==0.0 && W[2]==0.0 && W[3]==0.0) Volume();
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
the second method argument (matrix d of Akin, p420). The method also
returns the determinantof the Jacobian matrix since it is often needed
in integration procedures.
*/
double
IsoparametricQuadraticPrism::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, uint32_t gauss_point )
 {
    //
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

    // compose matrix DN = 3 x 27 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);

    // Forming maTRIX delta Akin, p.420
    for(uint32_t inode=0;inode<npe; inode++)
         {B(0,inode) = DNR[inode], B(1,inode) = DNS[inode], B(2,inode) = DNT[inode]; }

    B = JINV * B;

    return detJ;
 }


 void
 IsoparametricQuadraticPrism::JacobianAtIntegrationPoint( uint32_t gauss_point )
 {

    assert(gauss_point<gpe);

    //if(W[0]==0.0 && W[1]==0.0 && W[2]==0.0 && W[3]==0.0) Volume();

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }

/**

Method returns a vector with a size of 8, containing the consecutively
ordered local node numbers of the side of the element which lies at the
indicated model boundary. If the element is not on the model boundary an
error is reported and the vector is initialized to unspecified.

@section arguments Input Arguments

The parent Element, the target boundary of the the Model of elements,
and a  vector which will hold the the local indices of the identified
nodes (0...8).

@param fnids The resulting local node id's are returned into the third method argument.

@section implementation Implementation

While only the element knows which nodes are located at the mode boundary,
the FiniteElement knows in which order these appear.

@section application Application

To assign Neumann boundary conditions with a PDE operator for surface
integrals.
*/
void
IsoparametricQuadraticPrism::ConsecutiveNodesAtBoundary( const vector<uint32_t>& bnodes,
                                                         vector<uint32_t>& fnids )
 {
    fnids.resize(bnodes.size());

     if ( bnodes.size() != 9 || bnodes.size() != 6  )
       throw csmp::Exception( ERROR, "IsoparametricQuadraticPrism::ConsecutiveNodesAtBoundary",
                      "Cannot resolve node sequence for element boundary",
                      "Probably because element lies at two boundaries simultaneously" );

 } // end ConsecutiveNodesAtBoundary


/** Returns 4 local node ids of the nodes 4-9 located at the midsides of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

 */
void
IsoparametricQuadraticPrism::MidSideNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(9);
    ids[0] = 6;
    ids[1] = 7;
    ids[2] = 8;
    ids[3] = 9;
    ids[4] = 10;
    ids[5] = 11;
    ids[6] = 12;
    ids[7] = 13;
    ids[8] = 14;
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
IsoparametricQuadraticPrism::InnerRadius()
{
   vector<double> segms(spe);
   double                 sum(0.0);

   EdgeLengths( segms );
   for ( uint32_t i=0; i<spe; i++ ) sum += segms[i];

   if(AspectRatio()>4.0)
     cerr<<"\nIsoparametricQuadraticPrism:::InnerRadius: WARNING: function not applicable for this high element aspect ratio.\n"<<endl;

   return Volume() / (sum/6.);
}



/// @todo (3) Check barycenter
double
IsoparametricQuadraticPrism::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {

    dNr( 0.0, 0.0, 0.0, DNR );
    dNs( 0.0, 0.0, 0.0, DNS );
    dNt( 0.0, 0.0, 0.0, DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);

    for(uint32_t i=0; i<npe; i++)
    {
    B(0,i) = DNR[i], B(1,i) = DNS[i], B(2,i) = DNT[i];
    }

    B = JINV * B;

    return detJ;
 }


/// @todo (3) Check barycenter
void
IsoparametricQuadraticPrism::N_AtBaryCenter( std::vector<double>& N )
 {
    N.resize(npe);

    Nrst(0.0,0.0,0.0, N);
 }

 /**

Method uses the four integration points of the element to define linear
interpolation functions (just like in the linear tetrahedron). With these
functions variable values specified for the integration points are
linearily extrapolated to the elements nodes. Since the extrapolation
functions are linear and the distance between the integration points and
the nodes (in the local coordinate system which is used) is small, the
extrapolation is exact to numerical precision.

@section arguments Input Arguments

The first argument specifies how many variable values are supplied by
the IVAR vector for each integration point. This allows, for instance,
to extrapolate vector or tensor variables with this method which are
put into IVAR sequentially. Thus, IVAR has the dimensions 4 x vars-per-
integration point. The third method argument NVAR will hold the
extrapolated values for each node. Thus, its size should be
10 x vars-per-integration point.

@param NVAR The extrapolated values for each node of the element are returned into
the third method argument NVAR.

@section application Application

Sometimes, property values calculated at integration points are more
exact than if one calculates them at the nodes. In this case the method
offers the possibility to linearily extrapolate the exact values to the
nodes.
*/
void
IsoparametricQuadraticPrism::ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                                         const vector<double>& IVAR,
                                                                         vector<double>&       NVAR )
 const
{
   static bool          first_call(true);
   uint32_t            i,j,k;

   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   if ( IVAR.size() != (gpe*nvars) )
     throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticPrism::ExtrapolateIntegrationPointVariableToNodes",
                                  "Input vector must have 'nvars' x 6 entries");

   NVAR.resize( npe * nvars );

   // 1. Compute the local interpolation function coefficients for the hex which
   //    is defined by the four integration points.
   if ( first_call ) {
       // checking starting conditions
       if ( gpe != 4 )
       throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticPrism::ExtrapolateIntegrationPointVariableToNodes",
         "This method expects six integration points on which extrapolation functions will be based on" );

        first_call = false;
     }

    vector<double> Nj(npe), sum1(npe);
// test function coefficients
    for(k=0;k<nvars;k++)
        {
        for(j=0;j<npe;j++)
            {
            sum1[j]=0.;
            for ( i=0; i<gpe; i++ )
                {
                    //N_AtIntegrationPoint(i,Nj);
                    sum1[j]+=Nj[j]*IVAR[k*nvars+i];
                }
                sum1[j]/=static_cast<double>(gpe);
            NVAR[k*npe+j]=sum1[j];
          }
        }
  } // end ExtrapolateIntegrationPointVariableToNodes (STL vectors)



/// @todo Check  at all
void
IsoparametricQuadraticPrism::IntegralNN( DenseMatrix<DM_MIN>& IntNN )
{
    uint32_t  i; // index;
      IntNN.Resize(gpe,gpe);
    vector<double> Ni(npe);
      static DenseMatrix<DM_MIN>  BEHAT(gpe,gpe);
      static DenseMatrix<DM_MIN>  GDER(dim,gpe); //global der at IP

    // numerical integration:
    // looping over the 4 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for (i=0; i<gpe; i++ )
      {

        N_AtIntegrationPoint(i,Ni);
        //global derivatevs of interpolation functions at IP
        // Step 5. Calculate first order global derivateves of the interpolation functions
        dN_AtIntegrationPoint(GDER,i);
        // Step 6. For properties PROP(x,i) with local variation use nodal values and interpolation functions to
        //         calculate (interpolate) the properties to integration points;
        // Here integration points coincide with nodes: 0-24, 1-22, 2-21, 3-23, 4-20, 5-25
/* SKM fix: commented unused code
        if(i==0) index=24;
        else
            if(i==1) index=22;
        else
            if(i==2) index=21;
        else
            if(i==3) index=23;
        else
            if(i==4) index=20;
        else
            if(i==5) index=25;
*/
        // Step 7. Execute matrix operations, defining the matrix integrand. That involves the sum of
        //         products of element properties PROP(1,i), Ni and GDER.
        // Interpolate nodal values (properties) to IPs - function here

        //for(uint32_t j=0;j<dim;j++)
        //	{
        //		BEHAT(1,i)+=PROP(index,1)*Ni[i]*GDER(j,i);
        //	}

        //????????????????????????????????
        // Laplace equation on the 27_node hexa ===> should be NumIntegral...object ?
        //for(uint32_t j=0;j<gpe;j++)
        //	{
        //	for(uint32_t k=0;k<gpe;i++)
        //		BEHAT(k,j)+=ValOfJacobian * W[i]*(GDER(1,k)*GDER(1,j)+GDER(2,k)*GDER(2,j)+GDER(0,k)*GDER(0,j));
        // 	}
        // should be here - IP values (properties)*Ni*GDER

        // Step 8. Multiply resulting matrix by W[i] and ValOfJacobianInverse and add to the previous
        //         contributions to the element matrix;
        //	IntNN(1,i)+= BEHAT(1,i)* ValOfJacobian * W[i];
        // Extrapolate to the nodes?
         IntNN=BEHAT;
      }

 } // end method





/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be uptodate
void  IsoparametricQuadraticPrism::IntegrationPoint( uint32_t ip,
                                                     vector<double>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(3U); xyz[0]=xyz[1]=xyz[2]=0.;
     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

    for( auto i{0}; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint


//AP 2006
void
IsoparametricQuadraticPrism::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);

    matCoords = NXYZ;
}


} // namespace csmp

