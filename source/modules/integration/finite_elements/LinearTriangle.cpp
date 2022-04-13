#include "LinearTriangle.h"
#include "Exception.h"
#include "MJL_Edge.h"

using namespace std;

namespace csmp {

LinearTriangle::LinearTriangle()
 // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
 : FiniteElement( LINEAR_TRIANGLE, false, false, 1U )
 {
   dim = 2;
   itp = 1;
   npf = 2;
   npe = 3;
   fpe = 3;
   spe = 3;
   epe = 3;
   nne = 6;
   cne = 0;
   gpe = 0;
   XY.Resize(npe,dim);
   M.Resize(npe,npe);
   
   UsesLocalCoordinates(false);
   Isoparametric( false );
   SurfaceElement();
   ElementType(LINEAR_TRIANGLE);
 }



LinearTriangle::~LinearTriangle() 
 {
 }



double LinearTriangle::Volume()
{
   return std::fabs( 0.5 * ( XY(1,0)*XY(2,1) + XY(0,0)*XY(1,1) + 
                        XY(0,1)*XY(2,0) - XY(2,1)*XY(0,0) -
                        XY(2,0)*XY(1,1) - XY(1,0)*XY(0,1) ) );
}





void  LinearTriangle::EdgeLengths( std::vector<double>& len )
{
    double  sum;
    len.resize(npe);

    // segment 1
    sum    = (XY(1,0)-XY(0,0)) * (XY(1,0)-XY(0,0));
    sum   += (XY(1,1)-XY(0,1)) * (XY(1,1)-XY(0,1));
    len[0] = std::sqrt(sum);

    // segment 2
    sum    = (XY(2,0)-XY(1,0)) * (XY(2,0)-XY(1,0));
    sum   += (XY(2,1)-XY(1,1)) * (XY(2,1)-XY(1,1));
    len[1] = std::sqrt(sum);

    // segment 3
    sum    = (XY(2,0)-XY(0,0)) * (XY(2,0)-XY(0,0));
    sum   += (XY(2,1)-XY(0,1)) * (XY(2,1)-XY(0,1));
    len[2] = std::sqrt(sum); 

 } // end EdgeLengths
 




void LinearTriangle::TestFunctionCoefficients( DenseMatrix<DM_MIN>& M )
 {
   M.Resize(npe,npe);
   // ---------------------------------------------------------
   // Calculate the value of the test function coefficients
   // in global coordinates
   // result is stored in the 3 x 3 matrix as follows
   //
   //      | ai aj ak |
   //      | bi bj bk |  = coeffs
   //      | ci cj ck |
   //
   // ---------------------------------------------------------
   // calculation the exponents of the interpolation functions
   M(0,0) = XY(1,0) * XY(2,1) - XY(2,0) * XY(1,1);
   M(0,1) = XY(2,0) * XY(0,1) - XY(0,0) * XY(2,1);
   M(0,2) = XY(0,0) * XY(1,1) - XY(1,0) * XY(0,1); 
   
   M(1,0) = XY(1,1) - XY(2,1);
   M(1,1) = XY(2,1) - XY(0,1);
   M(1,2) = XY(0,1) - XY(1,1); 
   
   M(2,0) = XY(2,0) - XY(1,0);
   M(2,1) = XY(0,0) - XY(2,0);
   M(2,2) = XY(1,0) - XY(0,0); 
 }



void LinearTriangle::dN( DenseMatrix<DM_MIN>& B )
{
   B.Resize(dim,npe);
   // call to volume updates coordinate matrix
   double ae2 = 2. * Volume();   
   B(0,0) = (XY(1,1) - XY(2,1)) / ae2;
   B(0,1) = (XY(2,1) - XY(0,1)) / ae2;
   B(0,2) = (XY(0,1) - XY(1,1)) / ae2;
   B(1,0) = (XY(2,0) - XY(1,0)) / ae2;
   B(1,1) = (XY(0,0) - XY(2,0)) / ae2;
   B(1,2) = (XY(1,0) - XY(0,0)) / ae2;
}


/**
    The interpolation function derivates are constant inside of this linear element.
    The method returns the element area.
*/
double  LinearTriangle::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
   B.Resize(dim,npe);
   // call to volume updates coordinate matrix
   double ae2 = 2. * Volume();   
   B(0,0) = (XY(1,1) - XY(2,1)) / ae2;
   B(0,1) = (XY(2,1) - XY(0,1)) / ae2;
   B(0,2) = (XY(0,1) - XY(1,1)) / ae2;
   B(1,0) = (XY(2,0) - XY(1,0)) / ae2;
   B(1,1) = (XY(0,0) - XY(2,0)) / ae2;
   B(1,2) = (XY(1,0) - XY(0,0)) / ae2;
   
   // return the element area
   return ae2 / 2.;
 }


/**
    The interpolation function derivates are constant inside of this linear element.
    Thus, the location of the point itself does not have to be considered.
    The method returns the element area.
*/
double  LinearTriangle::dN_At( DenseMatrix<DM_MIN>& B, const std::vector<double>& )
 {
   return dN_AtBarycenter( B );
 }
 



void LinearTriangle::IntegralNN( DenseMatrix<DM_MIN>& CE )
{
   CE.Resize( npe, npe );
   
   CE(0,0) = Volume() /  6.0;
   CE(0,1) = CE(0,0) / 2.0;
   CE(0,2) = CE(0,1);
   CE(1,0) = CE(0,1);
   CE(1,1) = CE(0,0);
   CE(1,2) = CE(0,1);
   CE(2,0) = CE(0,1);
   CE(2,1) = CE(0,1);
   CE(2,2) = CE(0,0); 
}


/// segments are numbered like faces
void 
LinearTriangle::NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const
 {
    snids.resize(2);
    if ( segm_id == 0 ) {
         snids[0] = 1;
         snids[1] = 2;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 2;
         snids[1] = 0;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 0;
         snids[1] = 1;
      }
 }


/**
For this element, the faces are numbered such that face 0 lies opposite of 
node 0, face 1 node 1 etc.
*/
void LinearTriangle::NodesOfFace( uint32_t face_id, std::vector<uint32_t>& fnids ) const
 {
    fnids.resize(2);
    
    if      ( face_id == 0 )
      {
         fnids[0] = 1;
         fnids[1] = 2;
      }
    else if ( face_id == 1 )
      {
         fnids[0] = 2;
         fnids[1] = 0;
      }
    else if ( face_id == 2 )
      {
         fnids[0] = 0;
         fnids[1] = 1;
      }
    else
    std::cerr <<"\nLinearTriangle::NodesOfFace: Invalid Face ID requested: "<< face_id << std::endl;
 }



vector<uint32_t>  LinearTriangle::CornerNodesOfFace( uint32_t face_id ) const
 {
		switch (face_id) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nLinearTriangle::CornerNodesOfFace: face "<< face_id <<" does not exist.";
    return vector<uint32_t>{};
 }



/// returns the local  numbers of the nodes at the other end of the sgment that the argument node is on
std::vector<uint32_t>  LinearTriangle::NodesConnectedTo( uint32_t node_id ) const
  {
		switch ( node_id ) {
        case 0: return vector<uint32_t>{1,2};
        case 1: return vector<uint32_t>{2,0};
        case 2: return vector<uint32_t>{0,1};
      }
    cerr <<"\nLinearTriangle::NodesConnectedTo: node "<< node_id <<" does not exist.";
    return vector<uint32_t>{};
  }



/// The linear triangle is numbered counter-clockwise by default.
void LinearTriangle::CounterClockwiseNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }

/// clearly all nodes are corner nodes
void LinearTriangle::CornerNodes( std::vector<uint32_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
 }


CSMP_FEM_TYPE LinearTriangle::ElementTypeOfFace( uint32_t ) const
 {
    return LINEAR_BAR;
 }




double  LinearTriangle::AspectRatio()
{
   vector<double>  vec(spe);

   EdgeLengths( vec );

   // order segment
   set<double> segms;

   for ( uint32_t i{0U}; i<spe; i++ ) segms.insert( vec[i] );

   double segm1 = (*segms.begin()), 
            segm2 = (*segms.rbegin());
   
   return segm2 / segm1;
}




double  LinearTriangle::InnerRadius()
{
   vector<double> segms(npe);
   double         sum(0.0), vol;

   EdgeLengths( segms );
   for ( uint32_t i{0U}; i<segms.size(); i++ ) sum += segms[i];
   sum /= 2.0;
   vol  = Volume();
   vol /= sum;
   return vol;
}




void LinearTriangle::N( vector<double>& N, const vector<double>& xyz ) 
{
   // no sign is taken, thus method works with cw and ccw node numbering
   double ae2 = 1.0 / ( ( XY(1,0)*XY(2,1) + XY(0,0)*XY(1,1) + 
                          XY(0,1)*XY(2,0) - XY(2,1)*XY(0,0) -
                          XY(2,0)*XY(1,1) - XY(1,0)*XY(0,1) ) );

   // calculation the coefficients of the element interpolation functions
   a[0] = XY(1,0) * XY(2,1) - XY(2,0) * XY(1,1);
   a[1] = XY(2,0) * XY(0,1) - XY(0,0) * XY(2,1);
   a[2] = XY(0,0) * XY(1,1) - XY(1,0) * XY(0,1); 
   
   b[0] = XY(1,1) - XY(2,1);
   b[1] = XY(2,1) - XY(0,1);
   b[2] = XY(0,1) - XY(1,1); 
   
   c[0] = XY(2,0) - XY(1,0);
   c[1] = XY(0,0) - XY(2,0);
   c[2] = XY(1,0) - XY(0,0); 

   // summing the interpolation functions to get their value at (x,y)   
   N.resize(npe);
   for ( uint32_t i{0U}; i<npe; i++ )
     N[i] = ae2 * (a[i] + b[i] * xyz[0] + c[i] * xyz[1]);
     
} // end N




void LinearTriangle::N_AtBaryCenter( std::vector<double>& IPOL )
 {
    vector<double>  xyz(2);
    xyz[0] = (XY(0,0) + XY(1,0) + XY(2,0)) / 3.;
    xyz[1] = (XY(0,1) + XY(1,1) + XY(2,1)) / 3.;
    
    N( IPOL, xyz ); 
    
 } // N_AtBaryCenter




void  LinearTriangle::ConsecutiveNodesAtBoundary( const vector<uint32_t>& bnodes, 
                                                  vector<uint32_t>& fnids )
 {
     if ( bnodes.size() < 2 ) {
           cerr <<"\n\tnodes at boundary: ";
           for ( uint32_t j{0U}; j<bnodes.size(); j++ ) cerr << bnodes[j] <<" ";
           cerr << endl;
           throw csmp::Exception( ERROR, "LinearTriangle::ConsecutiveNodesAtBoundary",
                                      "Two nodes should be located at a boundary ! -",
                                      "correct input to meet this criterion." );
       }
             
     fnids.resize(2);        
              
     if ( (bnodes[0] == 0 && bnodes[1] == 1) || (bnodes[0] == 1 && bnodes[1] == 0) ) {
           fnids[0] = 0;
           fnids[1] = 1;
           return;
       }
     
     if ( (bnodes[0] == 1 && bnodes[1] == 2) || (bnodes[0] == 2 && bnodes[1] == 1) ) {
           fnids[0] = 1;
           fnids[1] = 2;
           return;
       }
     
     if ( (bnodes[0] == 2 && bnodes[1] == 0) || (bnodes[0] == 0 && bnodes[1] == 2) ) {
           fnids[0] = 2;
           fnids[1] = 0;
           return;
       }

 } // end ConsecutiveNodesAtBoundary




/**
     @author SKM 18/2/2016
*/
void  LinearTriangle::UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const
 {
     assert( face < Faces() );
     unrml.resize(2);
     // nodes 1 and 2
     if ( face == 0 ) {
          mjl::Edge  normal( mjl::Point(XY(1,0),XY(1,1)), mjl::Point(XY(2,0),XY(2,1)) );
          // rotating edge clockwise to find outward pointing normal to face
          normal.Rot();
          normal.NormalizeTo( 1. );
          unrml[0] = normal.Destination()[0];
          unrml[1] = normal.Destination()[1];
          return;
       }
     // nodes 2 and 0
     if ( face == 1 ) {
          mjl::Edge  normal( mjl::Point(XY(2,0),XY(2,1)), mjl::Point(XY(0,0),XY(0,1)) );
          normal.Rot();
          normal.NormalizeTo( 1. );
          unrml[0] = normal.Destination()[0];
          unrml[1] = normal.Destination()[1];
          return;
       }
     // nodes 0 and 1
     if ( face == 2 ) {
          mjl::Edge  normal( mjl::Point(XY(0,0),XY(0,1)), mjl::Point(XY(1,0),XY(1,1)) );
          normal.Rot();
          normal.NormalizeTo( 1. );
          unrml[0] = normal.Destination()[0];
          unrml[1] = normal.Destination()[1];
       }
 }




void LinearTriangle::OutputNodeDataToVTK( const char* file_name, 
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
           cerr <<"\nLinearTriangle::OutputNodeDataToVTK "; 
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
     for ( uint32_t i{0U}; i<npe; i++ )
       {
          // x, y
          for ( uint32_t j{0U}; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          // z
          ofs << 0.0 <<" ";
          ofs << endl;
       }
     ofs << endl;  
       
     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< npe+1 << endl;
     // number of points per cell, point1, ... point n
     ofs << npe <<" 0 1 2";
     ofs << endl;
     ofs << endl;
     
     // 4. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 5 << endl; // VTK_TRIANGLE
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );
     
     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" float"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1x6
           for ( uint32_t i{0U}; i<npe; i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          assert( DATA.Rows() == dim );
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 6
          for ( uint32_t i{0U}; i<DATA.Cols(); i++ ) {
               for ( uint32_t j{0U}; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               // the missing 3rd dimension
               ofs << 0.0 <<"  ";
               ofs << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nLinearTriangle3D::OutputNodeDataToVTK: file '"<< outfile;
     cout <<"' written successfully."<< endl;

} // end OutputNodeDataToVTK

} // end namespace csmp

