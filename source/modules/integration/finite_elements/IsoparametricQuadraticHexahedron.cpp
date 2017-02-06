#include "IsoparametricQuadraticHexahedron.h"
#include <fstream>
#include <cstring>
#include "Exception.h"
#include <climits>
#include "SparseMatrix.h"

using namespace std;


namespace csmp {


/**

@warning there seems to be problem with derivative at nodes 13 and 14, 13 - singular JACINV, 14 - orientation

*/
IsoparametricQuadraticHexahedron::IsoparametricQuadraticHexahedron( size_t integrationPoints )
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27, true, true, 2U ),
      DN(3,27),
      NXYZ(27,3),
      IP(integrationPoints,3),
      JETAL(3,3)
 {
    assert( integrationPoints == 4 or
            integrationPoints == 6 or
            integrationPoints == 8 );
   
    //AAM, 17.07.02
    dim = 3;
    itp = 2;
    npf = 9;
    npe = 27;
    fpe = 6;
    spe = 12;
    epe = 6;
    // nne - typical number of finite elements, sharing each node
    nne = 6;
    cne = 0;
    // Integration points
    gpe = integrationPoints;


    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27);

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
    // r-s-t/x-y-z triple-checked
    NXYZ(0,0) =-1.0;NXYZ(0,1) =-1.0;NXYZ(0,2) = -1.0; //1
    NXYZ(1,0) = 1.0;NXYZ(1,1) =-1.0;NXYZ(1,2) = -1.0; //2
    NXYZ(2,0) = 1.0;NXYZ(2,1) = 1.0;NXYZ(2,2) = -1.0; //3
    NXYZ(3,0) =-1.0;NXYZ(3,1) = 1.0;NXYZ(3,2) = -1.0; //4
    NXYZ(4,0) =-1.0;NXYZ(4,1) =-1.0;NXYZ(4,2) =  1.0; //5
    NXYZ(5,0) = 1.0;NXYZ(5,1) =-1.0;NXYZ(5,2) =  1.0; //6
    NXYZ(6,0) = 1.0;NXYZ(6,1) = 1.0;NXYZ(6,2) =  1.0; //7
    NXYZ(7,0) =-1.0;NXYZ(7,1) = 1.0;NXYZ(7,2) =  1.0; //8
    NXYZ(8,0) = 0.0;NXYZ(8,1) =-1.0;NXYZ(8,2) = -1.0; //9
    NXYZ(9,0) = 1.0;NXYZ(9,1) = 0.0;NXYZ(9,2) = -1.0; //10
    NXYZ(10,0)= 0.0;NXYZ(10,1)= 1.0;NXYZ(10,2) =-1.0; //11
    NXYZ(11,0)=-1.0;NXYZ(11,1)= 0.0;NXYZ(11,2) =-1.0; //12
    NXYZ(12,0)=-1.0;NXYZ(12,1)=-1.0;NXYZ(12,2) = 0.0; //13
    NXYZ(13,0)= 1.0;NXYZ(13,1)=-1.0;NXYZ(13,2) = 0.0; //14
    NXYZ(14,0)= 1.0;NXYZ(14,1)= 1.0;NXYZ(14,2) = 0.0; //15
    NXYZ(15,0)=-1.0;NXYZ(15,1)= 1.0;NXYZ(15,2) = 0.0; //16
    NXYZ(16,0)= 0.0;NXYZ(16,1)=-1.0;NXYZ(16,2) = 1.0; //17
    NXYZ(17,0)= 1.0;NXYZ(17,1)= 0.0;NXYZ(17,2) = 1.0; //18
    NXYZ(18,0)= 0.0;NXYZ(18,1)= 1.0;NXYZ(18,2) = 1.0; //19
    NXYZ(19,0)=-1.0;NXYZ(19,1)= 0.0;NXYZ(19,2) = 1.0; //20
    NXYZ(20,0)= 0.0;NXYZ(20,1)= 0.0;NXYZ(20,2) =-1.0; //21
    NXYZ(21,0)= 0.0;NXYZ(21,1)=-1.0;NXYZ(21,2) = 0.0; //22
    NXYZ(22,0)= 1.0;NXYZ(22,1)= 0.0;NXYZ(22,2) = 0.0; //23
    NXYZ(23,0)= 0.0;NXYZ(23,1)= 1.0;NXYZ(23,2) = 0.0; //24
    NXYZ(24,0)=-1.0;NXYZ(24,1)= 0.0;NXYZ(24,2) = 0.0; //25
    NXYZ(25,0)= 0.0;NXYZ(25,1)= 0.0;NXYZ(25,2) = 1.0; //26
    NXYZ(26,0)= 0.0;NXYZ(26,1)= 0.0;NXYZ(26,2) = 0.0; //27

    // For refernce etalon of Jacobian transformation matrix
    JETAL(0,0)=0.5;JETAL(0,1)=0.0;JETAL(0,2)=0.0;
    JETAL(1,0)=0.0;JETAL(1,1)=0.0;JETAL(1,2)=-0.5;
    JETAL(2,0)=0.0;JETAL(2,1)=0.5;JETAL(2,2)=0.0;

    W.resize(gpe);

    if( integrationPoints == 4 )
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
        const double64 oneDivSqrtThree=1./sqrt(3.0);
        const double64 sqrtTwoThirds=sqrt(2.0/3.0);

        IP(0,0) = 0.0, 			IP(0,1) = -sqrtTwoThirds, 	IP(0,2) = -oneDivSqrtThree;
        IP(1,0) = 0.0, 			IP(1,1) = sqrtTwoThirds,	IP(1,2) = -oneDivSqrtThree;
        IP(2,0) = -sqrtTwoThirds,IP(2,1) = 0.0, 				IP(2,2) =  oneDivSqrtThree;
        IP(3,0) =  sqrtTwoThirds,IP(3,1) = 0.0, 				IP(3,2) =  oneDivSqrtThree;
    }
    else if( integrationPoints == 6 )
    {
        //Numeric integration 3rd order, number of point m=6
        // Weights of the interpolation functions are 4/3, see Lo, p. 68

        double64 fourThirds=4.0/3.0;

        W[0] = fourThirds;
        W[1] = fourThirds;
        W[2] = fourThirds;
        W[3] = fourThirds;
        W[4] = fourThirds;
        W[5] = fourThirds;

        // Location of the integration points in r-s-t coordinates see Lo, p. 68

        IP(0,0) = -1.0,         IP(0,1) = 0.0, 				IP(0,2) = 0.0; // Node 25
        IP(1,0) =  1.0, 		IP(1,1) = 0.0,				IP(1,2) = 0.0; // Node 23
        IP(2,0) =  0.0,			IP(2,1) =-1.0, 				IP(2,2) = 0.0; // Node 22
        IP(3,0) =  0.0,			IP(3,1) = 1.0, 				IP(3,2) = 0.0; // Node 24
        IP(4,0) =  0.0,			IP(4,1) = 0.0, 				IP(4,2) =-1.0; // Node 21
        IP(5,0) =  0.0,			IP(5,1) = 0.0, 				IP(5,2) = 1.0; // Node 26

        /*
        IP(0,0) = 1.0/sqrt(6.0),		IP(0,1) = 1.0/sqrt(2.0),	IP(0,2) = -1.0/sqrt(3.0);
        IP(1,0) = 1.0/sqrt(6.0),		IP(1,1) =-1.0/sqrt(2.0),	IP(1,2) = -1.0/sqrt(3.0);
        IP(2,0) =-1.0/sqrt(6.0),		IP(2,1) = 1.0/sqrt(2.0),	IP(2,2) =  1.0/sqrt(3.0);
        IP(3,0) =-1.0/sqrt(6.0),		IP(3,1) =-1.0/sqrt(2.0),	IP(3,2) =  1.0/sqrt(3.0);
        IP(4,0) =-sqrt(2.0/3.0),		IP(4,1) = 0.0, 				IP(4,2) = -1.0/sqrt(3.0);
        IP(5,0) = sqrt(2.0/3.0),		IP(5,1) = 0.0, 				IP(5,2) =  1.0/sqrt(3.0);

        PROP(0,0)=PROP(1,0)=PROP(2,0)=PROP(3,0)=PROP(4,0)=PROP(5,0)=10.0;
        */
    }
    else if( integrationPoints == 8 )
    {
        // full 8-node Gauss integration scheme, see Akin, 1982, p.100
        // for extrapolation it should be used only!
        for(size_t i=0; i<integrationPoints;i++)
            W[i]=1.0;

        const double64 a1=0.577350269189626;

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
         cerr<<"***ERROR: IsoparametricQuadraticHexahedron::IsoparametricLinearHexahedron: N of integration points should be 4 or 8"<<endl;
         throw std::range_error
         ("IsoparametricQuadraticHexahedron::IsoparametricLinearHexahedron: Number of integration points should be 4 or 8" );
    }

 }



IsoparametricQuadraticHexahedron::~IsoparametricQuadraticHexahedron()
 {
 }


void
 IsoparametricQuadraticHexahedron::JacobianAtIntegrationPoint( size_t gauss_point )
 {
    assert(gauss_point<gpe);

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }









/** Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local corner node ID numbers
for the element
*/
void
IsoparametricQuadraticHexahedron::CornerNodes( std::vector<size_t>& ids ) const
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



/** Returns the counter-clockwise local node numbering for the element.

@param ids Returns the counter-clockwise local node numbering for the element.

*/
void
IsoparametricQuadraticHexahedron::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0]  = 0;
    ids[1]  = 8;
    ids[2]  = 1;
    ids[3]  = 9;
    ids[4]  = 2;
    ids[5]  =10;
    ids[6]  = 3;
    ids[7]  =11;
    ids[8]  =12;
    ids[9]  =13;
    ids[10] =14;
    ids[11] =15;
    ids[12] =4;
    ids[13] =16;
    ids[14] =5;
    ids[15] =17;
    ids[16] =6;
    ids[17] =18;
    ids[18] =7;
    ids[19] =19;
    ids[20] =20;
    ids[21] =21;
    ids[22] =22;
    ids[23] =23;
    ids[24] =24;
    ids[25] =26;
    ids[26] =25;
 }




void
IsoparametricQuadraticHexahedron::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    std::cout <<"\nIsoparametricQuadraticHexahedron::NodesOfSegment: This function has not been checked yet !"<< std::endl;
    snids.resize(3);
    if ( segm_id == 0 ) {
         snids[0] = 1;
         snids[1] = 0;
         snids[2] = 8;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 2;
         snids[1] = 1;
         snids[2] = 8;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 3;
         snids[1] = 2;
         snids[2] = 10;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 0;
         snids[1] = 3;
         snids[2] = 11;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 5;
         snids[1] = 4;
         snids[2] = 12;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 6;
         snids[1] = 5;
         snids[2] = 13;
      }
    else if ( segm_id == 6 ) {
         snids[0] = 7;
         snids[1] = 6;
         snids[2] = 14;
      }
    else if ( segm_id == 7 ) {
         snids[0] = 4;
         snids[1] = 7;
         snids[2] = 15;
      }
    else if ( segm_id == 8 ) {
         snids[0] = 4;
         snids[1] = 0;
         snids[2] = 16;
      }
    else if ( segm_id == 9 ) {
         snids[0] = 5;
         snids[1] = 1;
         snids[2] = 17;
      }
    else if ( segm_id == 10 ) {
         snids[0] = 6;
         snids[1] = 2;
         snids[2] = 18;
      }
    else if ( segm_id == 11 ) {
         snids[0] = 7;
         snids[1] = 3;
         snids[2] = 19;
      }

 } // end NodesOfSegment





/** For this element, the faces are numbered such that the lower left closest is 1 ->4
 */
void
IsoparametricQuadraticHexahedron::NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const
 {
    fnids.resize(9);

    if      ( face_id == 0 )
      {
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
    std::cerr <<"\nIsoparametricQuadraticHexahedron::NodesOfFace: Invalid Face ID requested: "<< face_id <<std::endl;
 }


CSMP_FEM_TYPE
IsoparametricQuadraticHexahedron::ElementTypeOfFace( size_t )  const
 {
    return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;
 }

double64 IsoparametricQuadraticHexahedron::WeightAtIntegrationPoint( size_t i ) const { return W[i]; }


void
IsoparametricQuadraticHexahedron::OutputNodeDataToVTK( const char* file_name,
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
           cerr <<"\nIsoparametricQuadraticHexahedron::OutputNodeDataToVTK ";
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
     for ( size_t i=0; i<npe; i++ )
       {
          for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point)) -corect
     // Cells
     // -----------------------------------------------------
     //ofs <<"CELLS "<< 6 <<" "<< 60 << endl;
     ofs <<"CELLS "<< 8 <<" "<< 72 << endl;
     // (1+4)X8
     // the 4 corner quads +4 midside nodes in form of polyline
     // as application of the
     /*
     ofs << 9 <<" 0 8 1 9 2 10 3 11 0" << endl;
     ofs << 9 <<" 0 8 1 13 5 16 4 12 0" << endl;
     ofs << 9 <<" 1 9 2 14 6 17 5 13 1" << endl;
     ofs << 9 <<" 2 14 6 18 7 15 3 10 2" << endl;
     ofs << 9 <<" 3 15 7 19 4 12 0 11 3" << endl;
     ofs << 9 <<" 4 16 5 17 6 18 7 19 4" << endl;
     ofs << endl;
     */
     ofs << 8 <<" 0 8 20 11 12 21 26 24" << endl;
     ofs << 8 <<" 8 1 9 20 21 13 22 26"  << endl;
     ofs << 8 <<" 20 9 2 10 26 22 14 23" << endl;
     ofs << 8 <<" 11 20 10 3 24 26 23 15" << endl;
     ofs << 8 <<" 12 21 26 24 4 16 25 19" << endl;
     ofs << 8 <<" 21 13 22 26 16 5 17 25" << endl;
     ofs << 8 <<" 26 22 14 23 25 17 6 18" << endl;
     ofs << 8 <<" 24 26 23 15 19 25 18 7" << endl;
     ofs << endl;


     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 8 << endl;
     ofs << 12 << endl; // VTK_HEX
     ofs << 12 << endl; // VTK_HEX
     ofs << 12 << endl; // VTK_HEX
     ofs << 12 << endl; // VTK_HEX
     ofs << 12 << endl; // VTK_HEX
     ofs << 12 << endl; // VTK_HEX
     ofs << 12 << endl; // VTK_HEX
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
     cout <<"\nIsoparametricQuadraticHexahedron::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

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
void IsoparametricQuadraticHexahedron::Nrst(
                    double64 r,
                    double64 s,
                    double64 t,
                    std::vector<double64>& N ) const
{
   const double64 rPlus=1.0+r;
   const double64 sPlus=1.0+s;
   const double64 tPlus=1.0+t;
   const double64 rMinus=1.0-r;
   const double64 sMinus=1.0-s;
   const double64 tMinus=1.0-t;

   N.resize(npe);
    // corner nodes
   N[0] = -0.125*r*s*t*rMinus*sMinus*tMinus;
   N[1] =  0.125*r*s*t*rPlus *sMinus*tMinus;
   N[2] = -0.125*r*s*t*rPlus *sPlus *tMinus;
   N[3] =  0.125*r*s*t*rMinus*sPlus *tMinus;
   N[4] =  0.125*r*s*t*rMinus*sMinus*tPlus;
   N[5] = -0.125*r*s*t*rPlus *sMinus*tPlus;
   N[6] =  0.125*r*s*t*rPlus *sPlus *tPlus;
   N[7] = -0.125*r*s*t*rMinus*sPlus *tPlus;

   N[8] =   0.25*s*t*(1-r*r)*sMinus *tMinus;
   N[9] =  -0.25*r*t*rPlus  *(1-s*s)*tMinus;
   N[10] = -0.25*s*t*(1-r*r)*sPlus  *tMinus;
   N[11] =  0.25*r*t*rMinus *(1-s*s)*tMinus;

    // Sign error? should be -
   N[12] =0.25*r*s*rMinus*sMinus*(1-t*t);
   // Should be +
   N[13] =-0.25*r*s*rPlus *sMinus*(1-t*t);

   //N[12] =-0.25*r*s*rPlus *sMinus*(1-t*t);
   //N[13] = 0.25*r*s*rMinus*sMinus*(1-t*t);

   N[14] = 0.25*r*s*rPlus *sPlus *(1-t*t);
   N[15] =-0.25*r*s*rMinus*sPlus *(1-t*t);

   N[16] =-0.25*s*t*(1-r*r)*sMinus *tPlus;
   N[17] = 0.25*r*t*rPlus  *(1-s*s)*tPlus;
   N[18] = 0.25*s*t*(1-r*r)*sPlus  *tPlus;
   N[19] =-0.25*r*t*rMinus *(1-s*s)*tPlus;

   N[20] =-0.5*t*(1-r*r)*(1-s*s)*tMinus;
   N[21] =-0.5*s*(1-r*r)*sMinus *(1-t*t);
   N[22] = 0.5*r*rPlus  *(1-s*s)*(1-t*t);

   N[23] = 0.5*s*(1-r*r)*sPlus  *(1-t*t);
   N[24] =-0.5*r*rMinus *(1-s*s)*(1-t*t);
   N[25] = 0.5*t*(1-r*r)*(1-s*s)*tPlus;
   N[26] = (1-r*r)*(1-s*s)*(1-t*t);

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
void IsoparametricQuadraticHexahedron::dNr (  double64 r,
                                              double64 s,
                                              double64 t,
                                              std::vector<double64>& DNR ) const
{
   DNR.resize(npe);

    DNR[0] = -0.125*s*t*(1.0 - r)*(1.0 - s)*(1.0 - t) + 0.125*r*s*t*(1.0 - s)*(1.0 - t);

    DNR[1] = 0.125*s*t*(1.0 + r)*(1.0 - s)*(1.0 - t) + 0.125*r*s*t*(1.0 - s)*(1.0 - t);

    DNR[2] = -0.125*s*t*(1.0 + r)*(1.0 + s)*(1.0 - t) - 0.125*r*s*t*(1.0 + s)*(1.0 - t);

    DNR[3] = 0.125*s*t*(1.0 - r)*(1.0 + s)*(1.0 - t) - 0.125*r*s*t*(1.0 + s)*(1.0 - t);

    DNR[4] = 0.125*s*t*(1.0 - r)*(1.0 - s)*(1.0 + t) - 0.125*r*s*t*(1.0 - s)*(1.0 + t);

    DNR[5] = -0.125*s*t*(1.0 + r)*(1.0 - s)*(1.0 + t) - 0.125*r*s*t*(1.0 - s)*(1.0 + t);

    DNR[6] = 0.125*s*t*(1.0 + r)*(1.0 + s)*(1.0 + t) + 0.125*r*s*t*(1.0 + s)*(1.0 + t);

    DNR[7] = -0.125*s*t*(1.0 - r)*(1.0 + s)*(1.0 + t) + 0.125*r*s*t*(1.0 + s)*(1.0 + t);

    DNR[8] = -0.50*r*s*t*(1.0 - s)*(1.0 - t);

    DNR[9] = -0.25*t*(1.0 + r)*(1 - s*s)*(1.0 - t) - 0.25*r*t*(1 - s*s)*(1.0 - t);

    DNR[10] = 0.50*r*s*t*(1.0 + s)*(1.0 - t);

    DNR[11] = 0.25*t*(1.0 - r)*(1 - s*s)*(1.0 - t) - 0.25*r*t*(1 - s*s)*(1.0 - t);

    DNR[12] = 0.25*s*(1.0 - r)*(1.0 - s)*(1 - t*t) - 0.25*r*s*(1.0 - s)*(1 - t*t);

    DNR[13] = -0.25*s*(1.0 + r)*(1.0 - s)*(1 - t*t) - 0.25*r*s*(1.0 - s)*(1 - t*t);

    DNR[14] = 0.25*s*(1.0 + r)*(1.0 + s)*(1 - t*t) + 0.25*r*s*(1.0 + s)*(1 - t*t);

    DNR[15] = -0.25*s*(1.0 - r)*(1.0 + s)*(1 - t*t) + 0.25*r*s*(1.0 + s)*(1 - t*t);

    DNR[16] = 0.50*r*s*t*(1.0 - s)*(1.0 + t);

    DNR[17] = 0.25*t*(1.0 + r)*(1 - s*s)*(1.0 + t) + 0.25*r*t*(1 - s*s)*(1.0 + t);

    DNR[18] = -0.50*r*s*t*(1.0 + s)*(1.0 + t);

    DNR[19] = -0.25*t*(1.0 - r)*(1 - s*s)*(1.0 + t) + 0.25*r*t*(1 - s*s)*(1.0 + t);

    DNR[20] = 1.0*r*t*(1 - s*s)*(1.0 - t);

    DNR[21] = 1.0*r*s*(1.0 - s)*(1 - t*t);

    DNR[22] = 0.5*(1.0 + r)*(1 - s*s)*(1 - t*t)+0.5*r*(1 - s*s)*(1 - t*t);

    DNR[23] = -1.0*r*s*(1.0 + s)*(1 - t*t);

    DNR[24] = -0.5*(1.0 - r)*(1 - s*s)*(1 - t*t) + 0.5*r*(1 - s*s)*(1 - t*t);

    DNR[25] = -1.0*r*t*(1 - s*s)*(1.0 + t);

    DNR[26] = -2.*r*(1 - s*s)*(1 - t*t);
}




void IsoparametricQuadraticHexahedron::dNs( double64 r,
                                            double64 s,
                                            double64 t,
                                            std::vector<double64>& DNS ) const
{
   const double64 rPlus=1.0+r;
   const double64 sPlus=1.0+s;
   const double64 tPlus=1.0+t;
   const double64 rMinus=1.0-r;
   const double64 sMinus=1.0-s;
   const double64 tMinus=1.0-t;

   DNS.resize(npe);
   DNS[0] =-.125*r*t*(1-r)*(1-s)*(1-t)+.125*r*s*t*(1-r)*(1-t);
   DNS[1] = .125*r*t*(1+r)*(1-s)*(1-t)-.125*r*s*t*(1+ r)*(1-t);
   DNS[2] =-.125*r*t*(1+r)*(1+s)*(1-t)-.125*r*s*t*(1+r)*(1-t);
   DNS[3] = .125*r*t*(1-r)*(1+s)*(1-t)+.125*r*s*t*(1-r)*(1-t);

   DNS[4] = .125*r*t*(1-r)*(1-s)*(1+t)-.125*r*s*t*(1-r)*(1+t);
   DNS[5] =-.125*r*t*(1+r)*(1-s)*(1+t)+.125*r*s*t*(1+r)*(1+t);
   DNS[6] = .125*r*t*(1+r)*(1+s)*(1+t)+.125*r*s*t*(1+r)*(1+t);
   DNS[7] =-.125*r*t*(1-r)*(1+s)*(1+t)-.125*r*s*t*(1-r)*(1+t);

   DNS[8]  = .25*t*(1-r*r)*(1-s)*(1-t)-.25*s*t*(1-r*r)*(1-t);
   DNS[9]  = .50*r*t*(1+r)*s*(1-t);
   DNS[10] =-.25*t*(1-r*r)*(1+s)*(1-t)-.25*s*t*(1-r*r)*(1-t);
   DNS[11] =-.50*r*s*t*(1-r)*(1-t);

   DNS[12] = .25*r*(1-r)*(1-s)*(1-t*t)-.25*r*s*(1-r)*(1-t*t);// - or
   DNS[13] =-.25*r*(1+r)*(1-s)*(1-t*t)+.25*r*s*(1+r)*(1-t*t);

   DNS[14] = .25*r*(1+r)*(1+s)*(1-t*t)+.25*r*s*(1+r)*(1-t*t);
   DNS[15] =-.25*r*(1-r)*(1+s)*(1-t*t)-.25*r*s*(1-r)*(1-t*t);

   DNS[16] =-0.25*(1-r*r)*tPlus *(t*sMinus-s*t);
   DNS[17] =-0.5*r*s*t*rPlus*tPlus;
   DNS[18] = 0.25*(1-r*r)*tPlus *(t*sPlus +s*t);
   DNS[19] = 0.5*r*s*t*rMinus*tPlus;

   DNS[20] = t*s*(1-r*r)*tMinus;
   DNS[21] =-0.5*(1-r*r)*(1-t*t)*(1-2*s);
   DNS[22] = -s*r*(1-t*t)*rPlus;
   DNS[23] = 0.5*(1-r*r)*(1-t*t)*(1+2*s);

   DNS[24] = s*r*rMinus*(1-t*t);
   DNS[25] =-s*t*(1-r*r)*tPlus;
   DNS[26] =-2.0*s*(1-r*r)*(1-t*t);

}


void IsoparametricQuadraticHexahedron::dNt( double64 r,
                                            double64 s,
                                            double64 t,
                                            std::vector<double64>& DNT ) const
{
   const double64 rPlus=1.0+r;
   const double64 sPlus=1.0+s;
   const double64 tPlus=1.0+t;
   const double64 rMinus=1.0-r;
   const double64 sMinus=1.0-s;

   DNT.resize(npe);
   DNT[0] =-.125*r*s*(1-r)*(1-s)*(1-t)+.125*r*s*t*(1-r)*(1-s);
   DNT[1] = .125*r*s*(1+r)*(1-s)*(1-t)-.125*r*s*t*(1+r)*(1-s);
   DNT[2] =-.125*r*s*(1+r)*(1+s)*(1-t)+.125*r*s*t*(1+r)*(1+s);
   DNT[3] = .125*r*s*(1-r)*(1+s)*(1-t)-.125*r*s*t*(1-r)*(1+s);

   DNT[4] = .125*r*s*(1-r)*(1-s)*(1+t)+.125*r*s*t*(1-r)*(1-s);
   DNT[5] =-.125*r*s*(1+r)*(1-s)*(1+t)-.125*r*s*t*(1+r)*(1-s);
   DNT[6] = .125*r*s*(1+r)*(1+s)*(1+t)+.125*r*s*t*(1+r)*(1+s);
   DNT[7] =-.125*r*s*(1-r)*(1+s)*(1+t)-.125*r*s*t*(1-r)*(1+s);

   DNT[8] = .25*s*(1-r*r)*(1-s)*(1-t)-.25*s*t*(1-r*r)*(1-s);
   DNT[9] =-.25*r*(1+r)*(1-s*s)*(1-t)+.25*r*t*(1+r)*(1-s*s);
   DNT[10]=-.25*s*(1-r*r)*(1+s)*(1-t)+.25*s*t*(1-r*r)*(1+s);
   DNT[11]= .25*r*(1-r)*(1-s*s)*(1-t)-.25*r*t*(1-r)*(1-s*s);

   DNT[12] = -.50*r*s*t*(1-r)*(1-s);
   DNT[13] =  .50*r*s*t*(1+r)*(1-s);

   DNT[14] = -.50*r*s*t*(1+r)*(1+s);
   DNT[15] =  .50*r*s*t*(1-r)*(1+s);

   DNT[16] = -0.25*(1-r*r)*sMinus*(s*tPlus+s*t);
   DNT[17] =  0.25*(1-s*s)*rPlus *(r*tPlus+r*t);
   DNT[18] =  0.25*(1-r*r)*sPlus *(s*tPlus+s*t);
   DNT[19] = -0.25*(1-s*s)*rMinus*(r*tPlus+r*t);

   DNT[20] = -0.5*(1-r*r)*(1-s*s)*(1-2*t);
   DNT[21] =  s*t*(1-r*r)*sMinus;
   DNT[22] = -r*t*(1-s*s)*rPlus;
   DNT[23] = -s*t*(1-r*r)*sPlus;

   DNT[24] = r*t*rMinus*(1-s*s);
   DNT[25] = 0.5*(1-r*r)*(1-s*s)*(1+2*t);
   DNT[26] =-2.0*t*(1-r*r)*(1-s*s);
}




/**

Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the triangle and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@param DN27 The interpolation-function derivative matrix is returned into the
second method argument.

*/
void
IsoparametricQuadraticHexahedron::dN( DenseMatrix<DM_MIN>& DN27 )
{
    DN27.Resize(dim,npe);
    M.Resize(dim,1);

    //DenseMatrix<DM_MIN> DERIVS(dim,npe);

     // Jacobian transformation to global coordinate system
     for ( size_t i=0; i<npe; i++ )
       {

       dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
       dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
       dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

       ////////// Debug Printout////////////////////////////////
       cout<<" Node="<<i+1<<" ("<<XY(i,0)<<","<<XY(i,1)<<","<<XY(i,2)<<")"<<endl;
       cout<<" N-point point: ("<<NXYZ(i,0)<<","<<NXYZ(i,1)<<","<< NXYZ(i,2)<<")"<<endl;
       cout<<" Ders : ("<<DNR[i]<<","<<DNS[i]<<","<<DNT[i]<<")"<<endl;
       ////////// Debug Printout////////////////////////////////

       ////////// Debug Printout////////////////////////////////
       //getchar();
       //DERIVS(0,i)=DNR[i];DERIVS(1,i)=DNS[i];DERIVS(2,i)=DNT[i];
       //if(i==12) getchar();//{DNR[i]=-DNR[i]; DNS[i]=-DNS[i]; DNT[i]=DNT[i];}
       //if(i==13) {DNR[i]=DNR[i]; DNS[i]=-DNS[i]; DNT[i]=DNT[i];}
       //getchar();
       ////////// Debug Printout////////////////////////////////
        JAC.Zero();

          // compute Jacobian matrix, its determinant and inverse
        if(i!=12 && i!=13)
        {
        ////////// Debug Printout////////////////////////////////
        Jacobian( DNR, DNS, DNT );

        cout<<" Jacobian Matrix: "<<endl;
        JAC.Out(cout);
        /////////////////////////////////////////////////////////

        JacobianInverse();

        ////////// Debug Printout////////////////////////////////
        cout<<" Jacobian Inverse Matrix: "<<endl;
        JINV.Out(cout);
        /////////////////////////////////////////////////////////
        }
        else
        if(i==12 || i==13) {
        Jacobian( DNR, DNS, DNT );
        //JAC=JETAL;

        ////////// Debug Printout////////////////////////////////
        cout<<" Jacobian Matrix: "<<endl;
        JAC.Out(cout);
        /////////////////////////////////////////////////////////
        getchar();
         JacobianInverse();
        ////////// Debug Printout////////////////////////////////
        cout<<" Jacobian Inverse Matrix: "<<endl;
        JINV.Out(cout);
        /////////////////////////////////////////////////////////
          }

          M(0,0)=DNR[i]; M(1,0)=DNS[i]; M(2,0)=DNT[i];
          // 3x3 * 3x1 = 3x1 gives the global DN entries
          JINV *= M;
          DN27(0,i) = JINV(0,0);
          DN27(1,i) = JINV(1,0);
          DN27(2,i) = JINV(2,0);
          JINV.Resize(dim,dim);
          JINV.Zero();
     }
     //OutputNodeDataToVTK( "etest1", "derivatives",DERIVS);

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
IsoparametricQuadraticHexahedron::dN( DenseMatrix<DM_MIN>& DN2,
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

    /////////////////////////////// Debug printout ///////////////////////////////////////////////
    cout<<" IsoparametricQuadraticHexahedron::dN  For given xyz=("<<xyz[1]<<","<<xyz[2]<<","<<xyz[3]<<"), rst=("<<
    rst[1]<<","<<rst[2]<<","<<rst[3]<<")"<<endl;
    cout<<" IsoparametricQuadraticHexahedron::dN  Matrix DN2: "<<endl;
    DN2.Out(cout);
    /////////////////////////////// Debug printout ///////////////////////////////////////////////

    return detJ;

 }

/** Outputs the shape functions at the point 'xyz' which must lie
within the Hexahedron.
*/
void
IsoparametricQuadraticHexahedron::N(
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
IsoparametricQuadraticHexahedron::dN_AtNode( DenseMatrix<DM_MIN>& B, size_t nd )
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
    for(size_t inode=0;inode<npe; inode++)
     {B(0,inode) = DNR[inode], B(1,inode) = DNS[inode], B(2,inode) = DNT[inode]; }

    B = JINV * B;

    return detJ;
 }

/** Projection function from rst->xyz
*/
void
IsoparametricQuadraticHexahedron::ParametricToPhysical(std::vector<double64> &rst, std::vector<double64>& xyz)
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
 */
void
IsoparametricQuadraticHexahedron::PhysicalToParametric(
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
            cout<<" IsoparametricQuadraticHexahedron::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            cout<<" IsoparametricQuadraticHexahedron::PhysicalToParametric: Real Point:\t"
                <<"x = "<<xyz[0]<<" ;\t"
                <<"y = "<<xyz[1]<<" ;\t"
                <<"z = "<<xyz[2]<<"\n";
            cout<<" IsoparametricQuadraticHexahedron::PhysicalToParametric: Found Point:\t"
               <<"x = "<<outxyz[0]<<" ;\t"
               <<"y = "<<outxyz[1]<<" ;\t"
               <<"z = "<<outxyz[2]<<"\n";
            cout<<" IsoparametricQuadraticHexahedron::PhysicalToParametric: Distance:\t"
               <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
               <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
               <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            cout<<" IsoparametricQuadraticHexahedron::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<" ;\t"
               <<"t = "<<rstHatK[2]<<"\n";

            std::vector<double64> N(npe,0.0);
            Nrst(rstHatK[0], rstHatK[1], rstHatK[2], N );
            cout<<" IsoparametricQuadraticHexahedron::PhysicalToParametric: Shape functions:\t"
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
                <<"N[17] = "<<N[17]<<" ;\t"
                <<"N[18] = "<<N[18]<<" ;\t"
                <<"N[19] = "<<N[19]<<" ;\t"
                <<"N[20] = "<<N[20]<<" ;\t"
                <<"N[21] = "<<N[21]<<" ;\t"
                <<"N[22] = "<<N[22]<<" ;\t"
                <<"N[23] = "<<N[23]<<" ;\t"
                <<"N[24] = "<<N[24]<<" ;\t"
                <<"N[25] = "<<N[25]<<" ;\t"
                <<"N[26] = "<<N[26]<<"\n";

            csmp::Exception( WARNING, "IsoparametricQuadraticHexahedron::PhysicalToParametric",
                          "Newton-Raphson iteration not converged");
        }

    }

    for(size_t i=0; i<dim; i++)
        rSt[i] = rstHatK[i];

}


 /** Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double64  IsoparametricQuadraticHexahedron::AspectRatio()
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
IsoparametricQuadraticHexahedron::EdgeLengths( std::vector<double64>& len )
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
IsoparametricQuadraticHexahedron::Volume()
{
    double64   area; // determinant
    size_t  i;

    // numerical integration:
    // looping over the 6 Gauss points calculating determinant
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

    //cout<<" Volume of QuadraticHex="<<area<<endl;
    return area;
}





/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

The interpolation function values are returned into the third argument.

*/
inline void
IsoparametricQuadraticHexahedron::N_AtIntegrationPoint( size_t ip, std::vector<double64>& N )
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
the second method argument (matrix d of Akin, p420). The method also
returns the determinantof the Jacobian matrix since it is often needed
in integration procedures.

*/
double64
IsoparametricQuadraticHexahedron::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, size_t gauss_point )
 {
    //
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
    assert( gauss_point < gpe );
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 27 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);

    // Forming maTRIX delta Akin, p.420
    for(size_t inode=0;inode<npe; inode++)
         {B(0,inode) = DNR[inode], B(1,inode) = DNS[inode], B(2,inode) = DNT[inode]; }

    B = JINV * B;

    return detJ;
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
IsoparametricQuadraticHexahedron::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                              vector<size_t>& fnids )
 {
    fnids.resize(bnodes.size());

     if ( bnodes.size() != 8 )
       throw csmp::Exception( ERROR, "IsoparametricQuadraticHexahedron::ConsecutiveNodesAtBoundary",
               "Cannot resolve node sequence for element boundary",
               "Probably because element lies at two boundaries simultaneously" );


 } // end ConsecutiveNodesAtBoundary






/** Returns 4 local node ids of the nodes 4-9 located at the midsides of
the quadratic tetrahedral element.

@param ids Returns an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

*/
void
IsoparametricQuadraticHexahedron::MidSideNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(12);
    ids[0] = 8;
    ids[1] = 9;
    ids[2] = 10;
    ids[3] = 11;
    ids[4] = 12;
    ids[5] = 13;
    ids[6] = 14;
    ids[7] = 15;
    ids[8] = 16;
    ids[9] = 17;
    ids[10] = 18;
    ids[11] = 19;
 }


double64
IsoparametricQuadraticHexahedron::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
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

    for(size_t i=0; i<npe; i++)
    {
    B(0,i) = DNR[i], B(1,i) = DNS[i], B(2,i) = DNT[i];
    }

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
IsoparametricQuadraticHexahedron::InnerRadius()
{
   vector<double64> segms(spe);
   double64                 sum(0.0);

   EdgeLengths( segms );
   for ( size_t i=0; i<spe; i++ ) sum += segms[i];

   // Function will produce unrelible value for high aspect ratio elements
   if(AspectRatio()>4.0)
   cout<<" IsoparametricLinearHexahedron::InnerRadius: ***WARNING: function not applicable for CURRENT HAR element"<<endl;


   return Volume() / (sum/6.);
}

void
IsoparametricQuadraticHexahedron::N_AtBaryCenter( std::vector<double64>& N )
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
IsoparametricQuadraticHexahedron::ExtrapolateIntegrationPointVariableToNodes(
                                                               size_t nvars,
                                                               const vector<double64>& IVAR,
                                                               vector<double64>&       NVAR )
 const
{
   static bool  first_call(true);


   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   if ( IVAR.size() != (gpe*nvars) )
     throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticHexahedron::ExtrapolateIntegrationPointVariableToNodes",
                                  "Input vector must have 'nvars' x 6 entries");

   NVAR.resize( npe * nvars );

   // 1. Compute the local interpolation function coefficients for the hex which
   //    is defined by the four integration points.
   if ( first_call ) {
       // checking starting conditions
       if ( gpe != 8 )
       throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticHexahedron::ExtrapolateIntegrationPointVariableToNodes",
         "This method expects 8 integration points on which extrapolation functions will be based on" );

        first_call = false;
     }

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

  const size_t cornerNodes=8;

  for ( size_t i=0; i<nvars; i++ )
    {
        for ( size_t k=0; k<gpe; k++ ) TEMP_IP(k,0)=IVAR[k*nvars +i];

            TEMP_N=MATRIX_A*TEMP_IP;
            //TEMP_N.Out(cout);
            for ( size_t j=0; j<cornerNodes; j++ )
            {
                NVAR[j*nvars+i]=TEMP_N(j,0);
            }
            //As well as midside nodes - base F1
            NVAR[8*nvars+i]=0.5*(TEMP_N(0,0)+TEMP_N(1,0));
            NVAR[9*nvars+i]=0.5*(TEMP_N(1,0)+TEMP_N(2,0));
            NVAR[10*nvars+i]=0.5*(TEMP_N(2,0)+TEMP_N(3,0));
            NVAR[11*nvars+i]=0.5*(TEMP_N(3,0)+TEMP_N(0,0));
            // middle cross section
            NVAR[12*nvars+i]=0.5*(TEMP_N(0,0)+TEMP_N(4,0));
            NVAR[13*nvars+i]=0.5*(TEMP_N(1,0)+TEMP_N(5,0));
            NVAR[14*nvars+i]=0.5*(TEMP_N(2,0)+TEMP_N(6,0));
            NVAR[15*nvars+i]=0.5*(TEMP_N(3,0)+TEMP_N(7,0));
            // upper F6
            NVAR[16*nvars+i]=0.5*(TEMP_N(4,0)+TEMP_N(5,0));
            NVAR[17*nvars+i]=0.5*(TEMP_N(5,0)+TEMP_N(6,0));
            NVAR[18*nvars+i]=0.5*(TEMP_N(6,0)+TEMP_N(7,0));
            NVAR[19*nvars+i]=0.5*(TEMP_N(7,0)+TEMP_N(4,0));
            // Center nodes of faces
            NVAR[20*nvars+i]=0.25*(TEMP_N(0,0)+TEMP_N(1,0)+TEMP_N(2,0)+TEMP_N(3,0));
            NVAR[21*nvars+i]=0.25*(TEMP_N(0,0)+TEMP_N(1,0)+TEMP_N(5,0)+TEMP_N(4,0));
            NVAR[22*nvars+i]=0.25*(TEMP_N(1,0)+TEMP_N(2,0)+TEMP_N(6,0)+TEMP_N(5,0));
            NVAR[23*nvars+i]=0.25*(TEMP_N(3,0)+TEMP_N(2,0)+TEMP_N(6,0)+TEMP_N(7,0));
            NVAR[24*nvars+i]=0.25*(TEMP_N(0,0)+TEMP_N(3,0)+TEMP_N(7,0)+TEMP_N(4,0));
            NVAR[25*nvars+i]=0.25*(TEMP_N(4,0)+TEMP_N(5,0)+TEMP_N(6,0)+TEMP_N(7,0));

            NVAR[26*nvars+i]=0.125*(TEMP_N(0,0)+TEMP_N(1,0)+TEMP_N(2,0)+TEMP_N(3,0)+(TEMP_N(4,0)+TEMP_N(5,0)+TEMP_N(6,0)+TEMP_N(7,0)));
    }
  } // end ExtrapolateIntegrationPointVariableToNodes (STL vectors)

/*
// tested AAM ok1
void
IsoparametricQuadraticHexahedron::IntegralNN( DenseMatrix<DM_MIN>& IntNN )
{
    double64   area=0.0;
    size_t  i, index;
    IntNN.Resize(gpe,gpe);
    vector<double64> Ni(npe);
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
        double64 ValOfJacobian=dN_AtIntegrationPoint(GDER,i);
        // Step 6. For properties PROP(x,i) with local variation use nodal values and interpolation functions to
        //         calculate (interpolate) the properties to integration points;
        // Here integration points coincide with nodes: 0-24, 1-22, 2-21, 3-23, 4-20, 5-25
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

        // Step 7. Execute matrix operations, defining the matrix integrand. That involves the sum of
        //         products of element properties PROP(1,i), Ni and GDER.
        // Interpolate nodal values (properties) to IPs - function here

        //for(size_t j=0;j<dim;j++)
        //	{
        //		BEHAT(1,i)+=PROP(index,1)*Ni[i]*GDER(j,i);
        //	}

        //????????????????????????????????
        // Laplace equation on the 27_node hexa ===> should be NumIntegral...object ?
        //for(size_t j=0;j<gpe;j++)
        //	{
        //	for(size_t k=0;k<gpe;i++)
        //		BEHAT(k,j)+=ValOfJacobian * W[i]*(GDER(1,k)*GDER(1,j)+GDER(2,k)*GDER(2,j)+GDER(0,k)*GDER(0,j));
        // 	}
        // should be here - IP values (properties)*Ni*GDER

        // Step 8. Multiply resulting matrix by W[i] and ValOfJacobianInverse and add to the previous
        //         contributions to the element matrix;
        //	IntNN(1,i)+= BEHAT(1,i)* ValOfJacobian * W[i];
        // Extrapolate to the nodes?
         IntNN=BEHAT;

      }


   }
*/


/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be uptodate
void  IsoparametricQuadraticHexahedron::IntegrationPoint( size_t ip,
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
  void
  IsoparametricQuadraticHexahedron::ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const
  {
    matCoords.Resize(npe, dim);
    matCoords.Fill(0.);
    matCoords = NXYZ;
  }


} // namespace csmp


