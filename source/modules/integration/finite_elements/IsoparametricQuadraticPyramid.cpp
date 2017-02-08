#include "IsoparametricQuadraticPyramid.h"
#include "Exception.h"
#include <fstream>
#include <cstring>
#include <climits>

using namespace std;

namespace csmp {


IsoparametricQuadraticPyramid::IsoparametricQuadraticPyramid()
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_QUADRATIC_PYRAMID13, true, true, 2U ),
      DN(3,14),
      NXYZ(14,3),
      IP(8,3)
 {
    //AAM, 15.08.02
   dim = 3;
   itp = 2;
   npf = 9;
   npe = 14;
   fpe = 5;
   spe = 8;
   epe = 5;
   // nne - typical number of finite elements, sharing each node
   nne = 6;
   cne = 0;
   // Integration points - using linear 4-point reference hex, therefore at curved bnds will be errors ...
   //
   gpe = 8;

   M.Resize(npe,npe);

   UsesLocalCoordinates(true);
   Isoparametric(true);
   VolumeElement();
   ElementType(ISOPARAMETRIC_QUADRATIC_PYRAMID13);

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
   NXYZ(0,0) =-1.0;NXYZ(0,1) =-1.0;NXYZ(0,2) =  0.0; //1
   NXYZ(1,0) = 1.0;NXYZ(1,1) =-1.0;NXYZ(1,2) =  0.0; //2
   NXYZ(2,0) = 1.0;NXYZ(2,1) = 1.0;NXYZ(2,2) =  0.0; //3
   NXYZ(3,0) =-1.0;NXYZ(3,1) = 1.0;NXYZ(3,2) =  0.0; //4
   NXYZ(4,0) = 0.0;NXYZ(4,1) = 0.0;NXYZ(4,2) =  1.0; //5
   NXYZ(5,0) = 0.0;NXYZ(5,1) =-1.0;NXYZ(5,2) =  0.0; //6
   NXYZ(6,0) = 1.0;NXYZ(6,1) = 0.0;NXYZ(6,2) =  0.0; //7
   NXYZ(7,0) = 0.0;NXYZ(7,1) = 1.0;NXYZ(7,2) =  0.0; //8
   NXYZ(8,0) =-1.0;NXYZ(8,1) = 0.0;NXYZ(8,2) =  0.0; //9
   NXYZ(9,0) =-0.5;NXYZ(9,1) =-0.5;NXYZ(9,2) =  0.5; //10
   NXYZ(10,0)= 0.5;NXYZ(10,1)=-0.5;NXYZ(10,2) = 0.5; //11
   NXYZ(11,0)= 0.5;NXYZ(11,1)= 0.5;NXYZ(11,2) = 0.5; //12
   NXYZ(12,0)=-0.5;NXYZ(12,1)= 0.5;NXYZ(12,2) = 0.5; //13
   NXYZ(13,0)= 0.0;NXYZ(13,1)= 0.0;NXYZ(13,2) = 0.0; //14

   W.resize( gpe );
   // Generate 8 integration points for full numeric integration + weights
   GenerateIntegrationPoints(IP,W);

 }



IsoparametricQuadraticPyramid::~IsoparametricQuadraticPyramid()
 {
 }


/** Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local corner node ID numbers
for the element.

*/
void
IsoparametricQuadraticPyramid::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(5);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
    ids[4] = 4;
 }

/** Returns the counter-clockwise local node numbering for the element.

@param ids Returns the counter-clockwise local node numbering for the element.

*/
void
IsoparametricQuadraticPyramid::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0]  = 0;
    ids[1]  = 5;
    ids[2]  = 1;
    ids[3]  = 6;
    ids[4]  = 2;
    ids[5]  = 7;
    ids[6]  = 3;
    ids[7]  = 8;
    ids[8]  = 9;
    ids[9]  =10;
    ids[10] =11;
    ids[11] =12;
    ids[12] =13;
    ids[13] =14;
 }




void
IsoparametricQuadraticPyramid::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    std::cout <<"\nIsoparametricQuadraticPyramid::NodesOfSegment: Not tested yet !"<< std::endl;
    snids.resize(3);
    if ( segm_id == 0 ) {
         snids[0] = 0;
         snids[1] = 1;
         snids[2] = 5;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 1;
         snids[1] = 2;
         snids[2] = 6;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 2;
         snids[1] = 3;
         snids[2] = 7;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 3;
         snids[1] = 0;
         snids[2] = 8;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 0;
         snids[1] = 4;
         snids[2] = 9;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 1;
         snids[1] = 4;
         snids[2] = 10;
      }
    else if ( segm_id == 6 ) {
         snids[0] = 2;
         snids[1] = 4;
         snids[2] = 11;
      }
    else if ( segm_id == 7 ) {
         snids[0] = 3;
         snids[1] = 4;
         snids[2] = 12;
      }
    else
    std::cout <<"\nIsoparametricQuadraticPyramid::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment




/**

For this element, the faces are numbered such that the lower left closest is 1 ->4
*/
/// @todo ANSYS convention
void
IsoparametricQuadraticPyramid::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {

    fnids.resize(9);
    if      ( face_id == 0 )
      {
        fnids.resize(9);
        fnids[0] = 0;
        fnids[1] = 8;
        fnids[2] = 1;
        fnids[3] = 9;
        fnids[4] = 2;
        fnids[5] = 10;
        fnids[6] = 3;
        fnids[7] = 11;
        fnids[8] = 20;
      }
    else if ( face_id == 1 )
      {
         fnids[0] = 0;
         fnids[1] = 8;
         fnids[2] = 1;
         fnids[3] = 13;
         fnids[4] = 5;
         fnids[5] = 16;
         fnids[6] = 4;
         fnids[7] = 12;
         fnids[8] = 21;
      }
    else if ( face_id == 2 )
      {
         fnids[0] = 1;
         fnids[1] = 9;
         fnids[2] = 2;
         fnids[3] = 14;
         fnids[4] = 6;
         fnids[5] = 17;
         fnids[6] = 5;
         fnids[7] = 13;
         fnids[8] = 22;
      }
    else if ( face_id == 3 )
      {
         fnids[0] = 2;
         fnids[1] = 14;
         fnids[2] = 6;
         fnids[3] = 18;
         fnids[4] = 7;
         fnids[5] = 15;
         fnids[6] = 3;
         fnids[7] = 10;
         fnids[8] = 23;
      }
    else if ( face_id == 4 )
      {
         fnids[0] = 0;
         fnids[1] = 11;
         fnids[2] = 3;
         fnids[3] = 15;
         fnids[4] = 7;
         fnids[5] = 19;
         fnids[6] = 4;
         fnids[7] = 12;
         fnids[8] = 24;
      }
    else if ( face_id == 5 )
      {
         fnids[0] = 4;
         fnids[1] = 16;
         fnids[2] = 5;
         fnids[3] = 17;
         fnids[4] = 6;
         fnids[5] = 18;
         fnids[6] = 7;
         fnids[7] = 19;
         fnids[8] = 25;
      }
    else
    std::cout <<"\nIsoparametricQuadraticPyramid::NodesOfFace: Invalid Face ID requested: "<< face_id << std::endl;
 }

/// @todo (3) Check ANSYS face numbering convention
CSMP_FEM_TYPE
IsoparametricQuadraticPyramid::ElementTypeOfFace( size_t face )  const
 {
    assert( face < fpe );
    if ( face == 0U ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;
    return ISOPARAMETRIC_QUADRATIC_TRIANGLE;
 }


double64 IsoparametricQuadraticPyramid::WeightAtIntegrationPoint( size_t i ) const { return W[i]; }


void
IsoparametricQuadraticPyramid::GenerateIntegrationPoints( DenseMatrix<DM_MIN> &Ip,
                                                          std::vector<double64>& We)
{
    const size_t numberOfQuadIntegrationPoints(4);
    const size_t numberOfDimensionsInPlane(2);
    DenseMatrix<DM_MIN> IPQUAD(numberOfQuadIntegrationPoints,numberOfDimensionsInPlane);
    const double64 constA(0.577350269189626);

    IPQUAD(0,0)=-constA; IPQUAD(0,1)=-constA;
    IPQUAD(1,0)= constA; IPQUAD(1,1)=-constA;
    IPQUAD(2,0)= constA; IPQUAD(2,1)= constA;
    IPQUAD(3,0)=-constA; IPQUAD(3,1)= constA;

    vector<double64> WQUAD(numberOfQuadIntegrationPoints);

    // Weights for the base integration
    WQUAD[0]=1.0;
    WQUAD[1]=1.0;
    WQUAD[2]=1.0;
    WQUAD[3]=1.0;

    const size_t numberOfPointsInTdimension=2;
    vector<double64> T(numberOfPointsInTdimension);
    T[0]=0.455848155988775; T[1]=0.877485177344559;
    vector<double64> b(numberOfPointsInTdimension);
    b[0]=0.100785882079825; b[1]=0.232547451253508;
    size_t i=0;
    for(size_t j=0;j<numberOfQuadIntegrationPoints;j++)
    {
    for(size_t k=0;k<numberOfPointsInTdimension;k++)
        {
            Ip(i,0)=T[k]*IPQUAD(j,0);
            Ip(i,1)=T[k]*IPQUAD(j,1);
            Ip(i,2)=1-T[k];
            We[i++]=WQUAD[j]*b[k];
        }
    }
}


void
IsoparametricQuadraticPyramid::OutputNodeDataToVTK( const char* file_name,
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
           cout <<"\nIsoparametricQuadraticPyramid::OutputNodeDataToVTK ";
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
     ofs <<"CELLS "<< 5 <<" "<< 36 << endl;
     // (1+4)X8
     // the 4 corner quads +4 midside nodes in form of polyline
     ofs << 8 <<" 0 5 1 6 2 7 3 8" << endl;
     ofs << 6 <<" 0 5 1 10 4 9"<< endl;
     ofs << 6 <<" 1 6 2 11 4 10" << endl;
     ofs << 6 <<" 2 11 4 12 3 7" << endl;
     ofs << 6 <<" 0 9 4 12 3 8" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 5 << endl;
     ofs << 4 << endl; // VTK_POLYLINE
     ofs << 4 << endl; // VTK_POLYLINE
     ofs << 4 << endl; // VTK_POLYLINE
     ofs << 4 << endl; // VTK_POLYLINE
     ofs << 4 << endl; // VTK_POLYLINE
     ofs << 4 << endl; // VTK_POLYLINE
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
     cout <<"\nIsoparametricQuadraticPyramid::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

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

Method can be used to compute property values at the integration points of
the element.
*/
void IsoparametricQuadraticPyramid::Nrst(
                    double64 r,
                    double64 s,
                    double64 t,
                    std::vector<double64>& N ) const
{
  N.resize(npe);

  N[13]= (1-r*r)*(1-s*s)*(1-t)*(0.5-r*r);

  // These shape functions are the corrected version of the ones that appear in the paper:
  // "A new family of finite elements: The pyramidal elements" by Zganski et al 1996
  // The isoparametric quadratic pyramid is a second order complete element. AP2006

  // Bedrosian paper - corrected, AAM, 20/08/04
  N[0] = 0.25*(-r-s-1)*((1-r)*(1-s)-t+r*s*t/(1-t)) + 0.25*N[13];
  N[1] = 0.25*( r-s-1)*((1+r)*(1-s)-t-r*s*t/(1-t)) + 0.25*N[13];
  N[2] = 0.25*( r+s-1)*((1+r)*(1+s)-t+r*s*t/(1-t)) + 0.25*N[13];
  N[3] = 0.25*(-r+s-1)*((1-r)*(1+s)-t-r*s*t/(1-t)) + 0.25*N[13];
  N[4] = t*(2*t-1);

  N[5] = (1+r-t)*(1-r-t)*(1-s-t)/(2*(1-t)) - 0.5*N[13];
  N[6] = (1+s-t)*(1-s-t)*(1+r-t)/(2*(1-t)) - 0.5*N[13];
  N[7] = (1+r-t)*(1-r-t)*(1+s-t)/(2*(1-t)) - 0.5*N[13];
  N[8] = (1+s-t)*(1-s-t)*(1-r-t)/(2*(1-t)) - 0.5*N[13];

  N[9] = t*(1-r-t)*(1-s-t)/(1-t);
  N[10]= t*(1+r-t)*(1-s-t)/(1-t);
  N[11]= t*(1+r-t)*(1+s-t)/(1-t);
  N[12]= t*(1-r-t)*(1+s-t)/(1-t);

}

/**

This method and the complementary methods dNr(), dNt() compute the shape function
derivates with respect to the local coordinate axis 'r','s' and 't'.

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
void IsoparametricQuadraticPyramid::dNr (
                double64 r,
                double64 s,
                double64 t,
                std::vector<double64>& DNR ) const
{
 DNR.resize(npe);

 if(t==1)
 {
   DNR[0] =  0.25;
   DNR[1] =  -0.25;
   DNR[2] =  -0.25;
   DNR[3] =  0.25;
   DNR[4] =  0.;
   DNR[5] =  0.;
   DNR[6] =  0.;
   DNR[7] =  0.;
   DNR[8] =  0.;
   DNR[9] =  -1.;
   DNR[10]=  1.;
   DNR[11]=  1.;
   DNR[12]=  -1.;
   DNR[13]=  0.;
 }
 else
 {
   DNR[13]= -2*r*(1-s*s)*(1-t)*(.5-r*r)-2*(1-r*r)*(1-s*s)*(1-t)*r;

   DNR[0] =- 0.25*(1.-r)*(1.-s)+0.25*t-0.25*r*s*t/(1.-t)+0.25*(-r-s-1.)*(-1.+s+s*t/(1.-t))   + 0.25*DNR[13];
   DNR[1] =  0.25*(r+1.)*(1.-s)-0.25*t-0.25*r*s*t/(1.-t)+0.25*(r-s-1.)*(1.-s-s*t/(1.-t)) 		 + 0.25*DNR[13];
   DNR[2] =  0.25*(r+1.)*(s+1.)-0.25*t+0.25*r*s*t/(1.-t)+0.25*(r+s-1.)*(s+1.+s*t/(1.-t)) 		 + 0.25*DNR[13];
   DNR[3] =- 0.25*(1.-r)*(s+1.)+0.25*t+0.25*r*s*t/(1.-t)+0.25*(-r+s-1.)*(-s-1.-s*t/(1.-t)) 	 + 0.25*DNR[13];
   DNR[4] =  0.;
   DNR[5] = (1.-r-t)*(1.-s-t)/(2.-2.*t)-(1.+r-t)*(1.-s-t)/(2.-2.*t)    - 0.5*DNR[13];
   DNR[6] = (1.+s-t)*(1.-s-t)/(2.-2.*t) 															 - 0.5*DNR[13];
   DNR[7] = (1.-r-t)*(1.+s-t)/(2.-2.*t)-(1.+s-t)*(1.+r-t)/(2.-2.*t)  	 - 0.5*DNR[13];
   DNR[8] =-(1.+s-t)*(1.-s-t)/(2.-2.*t) 															 - 0.5*DNR[13];
   DNR[9] =- t*(1.-s-t)/(1.-t);
   DNR[10]=  t*(1.-s-t)/(1.-t);
   DNR[11]=  t*(1.+s-t)/(1.-t);
   DNR[12]= -t*(1.+s-t)/(1.-t);
 }
}


/**

This method and the complementary methods dNr(), dNt() compute the shape function
derivates with respect to the local coordinate axis 'r','s' and 't'.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

@param DNS The shape function derivatives are returned into the second method
argument which is a floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.
*/

void IsoparametricQuadraticPyramid::dNs(
                double64 r,
                double64 s,
                double64 t,
                std::vector<double64>& DNS ) const
{
 DNS.resize(npe);

 if(t==1) //then r=0 and s=0
 {
   DNS[0] =  0.25;
   DNS[1] =  0.25;
   DNS[2] = -0.25;
   DNS[3] = -0.25;
   DNS[4] =  0.;
   DNS[5] =  0.;
   DNS[6] =  0.;
   DNS[7] =  0.;
   DNS[8] =  0.;
   DNS[9] =  -1.;
   DNS[10]=  -1.;
   DNS[11]=  1.;
   DNS[12]=  1.;
   DNS[13]=  0.;
 }
 else
 {
   DNS[13]= -2*(1-r*r)*s*(1-t)*(.5-r*r);

   DNS[0] = -0.25*(1.-r)*(1.-s)+0.25*t-0.25*r*s*t/(1.-t)+0.25*(-r-s-1.)*(-1.+r+r*t/(1.-t))    + 0.25*DNS[13];
   DNS[1] = -0.25*(r+1.)*(1.-s)+0.25*t+0.25*r*s*t/(1.-t)+0.25*(r-s-1.)*(-r-1.-r*t/(1.-t)) 	  + 0.25*DNS[13];
   DNS[2] =  0.25*(r+1.)*(s+1.)-0.25*t+0.25*r*s*t/(1.-t)+0.25*(r+s-1.)*(r+1.+r*t/(1.-t)) 		  + 0.25*DNS[13];
   DNS[3] =  0.25*(1.-r)*(s+1.)-0.25*t-0.25*r*s*t/(1.-t)+0.25*(-r+s-1.)*(1.-r-r*t/(1.-t)) 	  + 0.25*DNS[13];
   DNS[4] =  0.;
   DNS[5] =-(1.+r-t)*(1.-r-t)/(2.-2.*t)                                - 0.5*DNS[13];
   DNS[6] = (1.+r-t)*(1.-s-t)/(2.-2.*t)-(1.+s-t)*(1.+r-t)/(2.-2.*t) 	 - 0.5*DNS[13];
   DNS[7] = (1.+r-t)*(1.-r-t)/(2.-2.*t) 															 - 0.5*DNS[13];
   DNS[8] = (1.-r-t)*(1.-s-t)/(2.-2.*t)-(1.-r-t)*(1.+s-t)/(2.-2.*t) 	 - 0.5*DNS[13];
   DNS[9] = -t*(1.-r-t)/(1.-t);
   DNS[10]= -t*(1.+r-t)/(1.-t);
   DNS[11]=  t*(1.+r-t)/(1.-t);
   DNS[12]=  t*(1.-r-t)/(1.-t);

 }
}


/**

This method and the complementary methods dNr(), dNt() compute the shape function
derivates with respect to the local coordinate axis 'r','s' and 't'.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

@param DNT The shape function derivatives are returned into the second method
argument which is a floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.
*/

void IsoparametricQuadraticPyramid::dNt(
                double64 r,
                double64 s,
                double64 t,
                std::vector<double64>& DNT ) const
{
 DNT.resize(npe);
 if(t==1)//then r=0 and s=0
 {
  DNT[0] = 0.25*(-r-s-1.)*(-1.) - 0.25*0.5;
  DNT[1] = 0.25*(r-s-1.)*(-1.)  - 0.25*0.5;
  DNT[2] = 0.25*(r+s-1.)*(-1.)  - 0.25*0.5;
  DNT[3] = 0.25*(-r+s-1.)*(-1.) - 0.25*0.5;
  DNT[4] = 4.0-1.;
  DNT[5] = 0.5*0.5;
  DNT[6] = 0.5*0.5;
  DNT[7] = 0.5*0.5;
  DNT[8] = 0.5*0.5;
  DNT[9] = 1-2.;
  DNT[10]= 1-2.;
  DNT[11]= 1-2.;
  DNT[12]= 1-2.;
  DNT[13]= -0.5;
 }
 else
 {
  DNT[13]= -(1-r*r)*(1-s*s)*(.5-r*r);

  DNT[0] = 0.25*(-r-s-1.)*(-1.+r*s/(1.-t)+r*s*t/pow(1.-t,2)) + 0.25*DNT[13];
  DNT[1] = 0.25*(r-s-1.)*(-1.-r*s/(1.-t)-r*s*t/pow(1.-t,2))  + 0.25*DNT[13];
  DNT[2] = 0.25*(r+s-1.)*(-1.+r*s/(1.-t)+r*s*t/pow(1.-t,2))  + 0.25*DNT[13];
  DNT[3] = 0.25*(-r+s-1.)*(-1.-r*s/(1.-t)-r*s*t/pow(1.-t,2)) + 0.25*DNT[13];
  DNT[4] = 4.0*t-1.;
  DNT[5] = (-(1.-r-t)*(1.-s-t)-(1.+r-t)*(1.-s-t)-(1.+r-t)*(1.-r-t))/(2.-2.*t)+2.*(1.+r-t)*(1.-r-t)*(1.-s-t)/pow(2.-2.*t,2) - 0.5*DNT[13];
  DNT[6] = (-(1.+r-t)*(1.-s-t)-(1.+s-t)*(1.+r-t)-(1.+s-t)*(1.-s-t))/(2.-2.*t)+2.*(1.+s-t)*(1.-s-t)*(1.+r-t)/pow(2.-2.*t,2) - 0.5*DNT[13];
  DNT[7] = (-(1.-r-t)*(1.+s-t)-(1.+s-t)*(1.+r-t)-(1.+r-t)*(1.-r-t))/(2.-2.*t)+2.*(1.+r-t)*(1.-r-t)*(1.+s-t)/pow(2.-2.*t,2) - 0.5*DNT[13];
  DNT[8] = (-(1.-r-t)*(1.-s-t)-(1.-r-t)*(1.+s-t)-(1.+s-t)*(1.-s-t))/(2.-2.*t)+2.*(1.+s-t)*(1.-s-t)*(1.-r-t)/pow(2.-2.*t,2) - 0.5*DNT[13];
  DNT[9] = ((1.-r-t)*(1.-s-t)-t*(1.-s-t)-t*(1.-r-t))/(1.-t)+t*(1.-r-t)*(1.-s-t)/pow(1.-t,2);
  DNT[10]= ((1.+r-t)*(1.-s-t)-t*(1.-s-t)-t*(1.+r-t))/(1.-t)+t*(1.+r-t)*(1.-s-t)/pow(1.-t,2);
  DNT[11]= ((1.+r-t)*(1.+s-t)-t*(1.+s-t)-t*(1.+r-t))/(1.-t)+t*(1.+r-t)*(1.+s-t)/pow(1.-t,2);
  DNT[12]= ((1.-r-t)*(1.+s-t)-t*(1.+s-t)-t*(1.-r-t))/(1.-t)+t*(1.-r-t)*(1.+s-t)/pow(1.-t,2);

 }
} // dNt





/** Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN14 The interpolation-function derivative matrix is returned into the
second method argument.
*/
void
IsoparametricQuadraticPyramid::dN( DenseMatrix<DM_MIN>& DN14 )
{
    DN14.Resize(dim,npe);
    M.Resize(dim,1);
    bool ldbug(true);

     // Jacobian transformation to global coordinate system
     for ( size_t i=0; i<npe; i++ )
      {
          // here the global coordinates come in
          dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
          dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
          dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

          if(ldbug){
             cout<< " IsoparametricQuadraticPyramid::dN() :"<<endl;
             cout<<" Node="<<i+1<<" ("<<XY(i,0)<<","<<XY(i,1)<<","<<XY(i,2)<<")"<<endl;
             cout<<" N-point point: ("<<NXYZ(i,0)<<","<<NXYZ(i,1)<<","<< NXYZ(i,2)<<")"<<endl;
             cout<<" Ders : ("<<DNR[i]<<","<<DNS[i]<<","<<DNT[i]<<")"<<endl;
             //if(i==12) {DNR[i]=-DNR[i]; DNS[i]=-DNS[i]; DNT[i]=DNT[i];}
             //if(i==13) {DNR[i]=DNR[i]; DNS[i]=-DNS[i]; DNT[i]=DNT[i];}
             getchar();
           }
          // compute Jacobian matrix, its determinant and inverse
          Jacobian( DNR, DNS, DNT );
        if(ldbug){
          cout<<" At node  #: "<<i+1<<endl<<endl;
          cout<<" Jacobian Matrix: "<<endl;
          JAC.Out(cout);
        }
          JacobianInverse();
          ////////// Debug Printout////////////////////////////////
          if(ldbug){
            cout<<" Jacobian Inverse Matrix: "<<endl;
            JINV.Out(cout);
          }
          /////////////////////////////////////////////////////////
          M(0,0)=DNR[i]; M(1,0)=DNS[i]; M(2,0)=DNT[i];
          // 3x3 * 3x1 = 3x1 gives the global DN entries
          JINV *= M;
          DN14(0,i) = JINV(0,0);
          DN14(1,i) = JINV(1,0);
          DN14(2,i) = JINV(2,0);
          JINV.Resize(dim,dim);
      }
} // end dN



/**

Computes derivatives of shape functions at a point XYZ given in global coordinates
within the element (also given in global coordinates).

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at the point XY are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument.

*/
double64
IsoparametricQuadraticPyramid::dN_At( DenseMatrix<DM_MIN>& DN2,
                                                      const vector<double64>& xyz  )
  {
    bool ldebug(false);
    vector<double64> rst(dim);
    DN2.Resize(dim,dim);

    PhysicalToParametric(rst, xyz);

    // here find derivatives
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    DN2  = JINV;

    DenseMatrix<DM_MIN> DN5(dim,npe);
    dN(DN5);
    DN2*=DN5;

    /////////////////////////////// Debug printout ///////////////////////////////////////////////
    if(ldebug) {
      cout<<" IsoparametricQuadraticPyramid::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
      rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
      cout<<" IsoparametricQuadraticPyramid::dN  Matrix DN2: "<<endl;
      DN2.Out(cout);
    }
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

    return detJ;

 } // dN(XYZ)



/** Outputs the shape functions at the point 'xyz' which must lie
inside the Pyramid.

@section arguments Input Arguments

Vector of the point's physical coordinates.

@param N Vector of the shape functions at the nodes of a pyramid.

*/
void
IsoparametricQuadraticPyramid::N( std::vector<double64>& N,
                                  const std::vector<double64>& xyz )
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
IsoparametricQuadraticPyramid::dN_AtNode( DenseMatrix<DM_MIN>& B, size_t nd )
 {
    assert( nd < npe );
    dNr( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR );
    dNs( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS );
    dNt( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 14 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    for(size_t inode=0;inode<npe; inode++)
      {
        B(0,inode) = DNR[inode], B(1,inode) = DNS[inode], B(2,inode) = DNT[inode];
      }

    B = JINV * B;

    return detJ;
 }

/** Method implements projection from parametric to physical space, using shape
functions.

@section arguments Input Arguments

Vector of point's parametric coordinates rst.

@param xyz Vector of point's physical coordinates xyz.

@section implementation Implementation

Method uses standard implementation of shape functions.
*/


void
IsoparametricQuadraticPyramid::ParametricToPhysical(std::vector<double64> &rst, std::vector<double64>& xyz)
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

/**

Method implements projection from physical to parametric space, using Newton-Raphson
iterative method.

@section arguments Input Arguments

Vector of the point's physical coordinates xyz.

@param rSt Vector of the point's parametric coordinates rst.

@section implementation Implementation

Method uses standard implementation of the Newton-Raphson iterative method.
*/
void
IsoparametricQuadraticPyramid::PhysicalToParametric(
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

        vector<double64> rstHatK_PlusOne(dim);
        double64 minDistanceFromGivenPoint;

        const double64 constantMu               = 1.0;
        const size_t numberOfFirstIterrations   = 5;
        const size_t maxNumberOfIterrations     = 20;
        const size_t numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double64 incrementR = 2.0/(numberOfFirstIterrations-1);
        const double64 incrementS = 2.0/(numberOfFirstIterrations-1);
        const double64 incrementT = 1.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = -1.0;
        for(size_t i=0;i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = -1.0;
            for(size_t j=0;j<numberOfFirstIterrations;j++)
            {
                rstHatK_PlusOne[2] = 0.0;
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
            cout<<" IsoparametricQuadraticPyramid::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            cout<<" IsoparametricQuadraticPyramid::PhysicalToParametric: Real Point:\t"
                <<"x = "<<xyz[0]<<" ;\t"
                <<"y = "<<xyz[1]<<" ;\t"
                <<"z = "<<xyz[2]<<"\n";
            cout<<" IsoparametricQuadraticPyramid::PhysicalToParametric: Found Point:\t"
               <<"x = "<<outxyz[0]<<" ;\t"
               <<"y = "<<outxyz[1]<<" ;\t"
               <<"z = "<<outxyz[2]<<"\n";
            cout<<" IsoparametricQuadraticPyramid::PhysicalToParametric: Distance:\t"
               <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
               <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
               <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            cout<<" IsoparametricQuadraticPyramid::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<" ;\t"
               <<"t = "<<rstHatK[2]<<"\n";

            std::vector<double64> N(npe,0.0);
            Nrst(rstHatK[0], rstHatK[1], rstHatK[2], N );
            cout<<" IsoparametricQuadraticPyramid::PhysicalToParametric: Shape functions:\t"
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
                <<"N[13] = "<<N[13]<<"\n";

            csmp::Exception( CSMP_WARNING, "IsoparametricQuadraticPyramid::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for(size_t i=0; i<dim; i++)
        rSt[i] = rstHatK[i];

}


 /**

Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double64  IsoparametricQuadraticPyramid::AspectRatio()
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




/** Segment length is calculated directly as the length of linear segments.

@param len The lengths of the 12 segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as simple Euclidian distance between corner
vertices of the edge.
*/
void
IsoparametricQuadraticPyramid::EdgeLengths( std::vector<double64>& len )
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
    sum     = (XY(4,0)-XY(0,0)) * (XY(4,0)-XY(0,0));
    sum    += (XY(4,1)-XY(0,1)) * (XY(4,1)-XY(0,1));
    sum    += (XY(4,2)-XY(0,2)) * (XY(4,2)-XY(0,2));
    len[4]  = sqrt(sum);
    // segment 6
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


// Standard numeric integration over 8 points
// Tested AAM
/** Computes the volume of a given element.

@return Method returns volume of pyramid element.

*/
double64
IsoparametricQuadraticPyramid::Volume()
{
    bool ldebug(false);
    double64   area; // determinant
    size_t  i;

    // numerical integration:
    // looping over the 8 Gauss points calculating determinant
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

    if(ldebug) cout<<" IsoparametricQuadraticPyramid::Volume() Volume="<<area<<endl;
    return area;

}

/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

The index of the integration point.

@param N The interpolation function values are returned into the second argument.

*/
void IsoparametricQuadraticPyramid::N_AtIntegrationPoint( size_t ip, std::vector<double64>& N )
 {
    assert( ip < gpe );

    N.resize(npe);

    if(W[0]==0.0 && W[1]==0.0 && W[2]==0.0 && W[3]==0.0) Volume();
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
the second method argument (matrix d of Akin, p420). The method also
returns the determinantof the Jacobian matrix since it is often needed
in integration procedures.
*/
double64
IsoparametricQuadraticPyramid::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, size_t gauss_point )
 {
    //
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
    assert( gauss_point < gpe );

    if(W[0]==0.0 && W[1]==0.0 && W[2]==0.0 && W[3]==0.0) Volume();

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 14 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);

      // Forming maTRIX delta Akin, p.420
      for(size_t inode=0;inode<npe; inode++)
         {
            B(0,inode) = DNR[inode], B(1,inode) = DNS[inode], B(2,inode) = DNT[inode];
         }

    B = JINV * B;

    return detJ;
 }


 void
 IsoparametricQuadraticPyramid::JacobianAtIntegrationPoint( size_t gauss_point )
 {

    assert(gauss_point<gpe);

    if(W[0]==0.0 && W[1]==0.0 && W[2]==0.0 && W[3]==0.0) Volume();

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }

/**

Method returns a vector with a size of 3, containing the consecutively
ordered local node numbers of the side of the element which lies at the
indicated model boundary. If the element is not on the model boundary an
error is reported and the vector is initialized to unspecified.

@section arguments Input Arguments

The parent Element, the target boundary of the the Model of elements,
and a vector which will hold the the local indices of the identified
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
IsoparametricQuadraticPyramid::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                           vector<size_t>& fnids )
 {
     fnids.resize(bnodes.size());

     if ( bnodes.size() != 8 || bnodes.size() != 6  )
       throw csmp::Exception( CSMP_ERROR, "IsoparametricQuadraticPyramid::ConsecutiveNodesAtBoundary",
               "Cannot resolve node sequence for element boundary",
               "Probably because element lies at two boundaries simultaneously" );


 } // end ConsecutiveNodesAtBoundary


/**

Returns 4 local node ids of the nodes 4-9 located at the midsides of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

*/
void
IsoparametricQuadraticPyramid::MidSideNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(8);
    ids[0] = 5;
    ids[1] = 6;
    ids[2] = 7;
    ids[3] = 8;
    ids[4] = 9;
    ids[5] = 10;
    ids[6] = 11;
    ids[7] = 12;
 }


/** Computes derivateves of interpolation functions at barycenter of element.

@return The interpolation function values are returned into the first argument.
Method returns determinanat of the Jacobian transformation matrix from
parametric to physical space.
*/
double64
IsoparametricQuadraticPyramid::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    double64 oneFourth(0.25);

    dNr( 0.0, 0.0, oneFourth, DNR );
    dNs( 0.0, 0.0, oneFourth, DNS );
    dNt( 0.0, 0.0, oneFourth, DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);

    for(size_t i=0; i<npe; i++)
      B(0,i) = DNR[i], B(1,i) = DNS[i], B(2,i) = DNT[i];

    B = JINV * B;
   return detJ;
 }


/** Computes interpolation functions at barycenter of element.

@param N The interpolation function values are returned into first argument.

*/
void
IsoparametricQuadraticPyramid::N_AtBaryCenter( std::vector<double64>& N )
 {
    N.resize(npe);

    double64 oneFourth(0.25);

    Nrst(0.0,0.0,oneFourth, N);
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
double64
IsoparametricQuadraticPyramid::InnerRadius()
{
   vector<double64> segms(spe);
   double64                 sum(0.0);

   EdgeLengths( segms );
   for ( size_t i=0; i<spe; i++ ) sum += segms[i];

   if(AspectRatio()>4.0)
   cout<<" IsoparametricQuadraticPyramid::InnerRadius: ***WARNING: function not applicable for CURRENT HAR element"<<endl;

   return Volume() / (sum/6.);
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
IsoparametricQuadraticPyramid::ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                        const vector<double64>& IVAR,
                                                                        vector<double64>&       NVAR )
 const
{
   static bool first_call(true);
   size_t      i,j,k;

   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   if ( IVAR.size() != (gpe*nvars) )
     throw csmp::Exception( CSMP_FATAL_ERROR, "IsoparametricQuadraticPyramid::ExtrapolateIntegrationPointVariableToNodes",
                                  "Input vector must have 'nvars' x 6 entries");

   NVAR.resize( npe * nvars );

   // 1. Compute the local interpolation function coefficients for the hex which
   //    is defined by the four integration points.
   if ( first_call ) {
       // checking starting conditions
       if ( gpe != 6 )
       throw csmp::Exception( CSMP_FATAL_ERROR, "IsoparametricQuadraticPyramid::ExtrapolateIntegrationPointVariableToNodes",
         "This method expects six integration points on which extrapolation functions will be based on" );

        first_call = false;
     }

    vector<double64> Nj(npe), sum1(npe);

  //test function coefficients
    for(k=0;k<nvars;k++)
        {
        for(j=0;j<npe;j++)
            {
            sum1[j]=0.0;
            for ( i=0; i<gpe; i++ )
                {
                    //N_AtIntegrationPoint(i,Nj);
                    sum1[j]+=Nj[j]*IVAR[k*nvars+i];
                }
                sum1[j]/=static_cast<double64>(gpe);
            NVAR[k*npe+j]=sum1[j];
          }
        }
  } // end ExtrapolateIntegrationPointVariableToNodes (STL vectors)




/** integration point location transformed into global coordinates

@note the matrix XYZ must be uptodate

*/
void  IsoparametricQuadraticPyramid::IntegrationPoint( size_t ip,
                                                       vector<double64>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(3U); xyz[0]=xyz[1]=xyz[2]=0.;
     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

    for( size_t i=0U; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint

// AP 2006
void IsoparametricQuadraticPyramid::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords = NXYZ;
}



} // end namespace csmp
