#include "mechanics.h"

#ifdef CSMP_WITH_MESCHACH
#include "meschach.h"
#include "matrix.h"
#include "matrix2.h"
#endif

using namespace std;

namespace csmp {

/**
@file mechanics.cpp
@addtogroup CSMPglobalFunctions
@{
*/

/** Computes the 2D material property matrix for plane strain
computations, see Cook et al. p. 21.  

@section arguments Input Arguments 

Young's modulus and  Poisson's ratio as double64 variables and a 
Meschach matrix into which the result is returned.  
*/

void planeStressMatrix( double64 E,  // Young's modulus, 
                        double64 nu, // Poisson's ratio 
                        DenseMatrix<DM_MIN>& D )
 {
    D.Resize(3,3);
    D.Zero();
    assert( nu < 0.5 );
    assert( E > 0. );
    D(1,0) = D(0,1) = nu;
    D(0,0) = D(1,1) = 1.0;
    D(2,2) = (1.0-nu) / 2.0;
    D *= (E / (1.0 - nu*nu));
 }




/** Cook et al. p. 21; Young's modulus,  Poisson's ratio

    stiffness form: to compute stresses from strains

*/
void planeStrainMatrix( double64 E, double64 nu, DenseMatrix<DM_MIN>& D )
 {
    D.Resize(3,3);
    D.Zero();
    assert( nu < 0.5 );
    assert( E > 0. );
    D(0,0) = D(1,1) = (1.0 - nu);
    D(1,0) = D(0,1) = nu;
    D(2,2) = (1. - 2*nu) / 2.;
    D *= (E / ((1. + nu)*(1. - 2*nu)));
 }




/** 3D version: Smith & Griffith (1998), p. 40

is this in compliance form? - i.e.  e = D * signma

void stiffnessMatrix( double64 E, double64 nu, DenseMatrix<DM_MIN>& D )
 {
    D.Resize(6,6);
    assert( nu < 0.5 );
    assert( E > 0. );
    const double64 frac  = nu / (1. - nu);
    const double64 coeff = (E*(1. - nu)) / ((1. + nu)*(1. - 2.*nu));
    
    D.Zero();
    D(0,0) = D(1,1) = D(2,2) = 1.;
    D(3,3) = D(4,4) = D(5,5) = (1. - 2.*nu) / (2.*(1. - nu));
    D(0,1) = D(0,2) = D(1,2) = frac;
    D(1,0) = D(2,0) = D(2,1) = frac;
    D *= coeff;
 }
*/



/** Stresses from strains (Wikipedia), isotropic material
    
    @attention inverse relationship to get from strains to stresses
    
    sigma = D * e
    
    @attention for multiplication, the last 3 terms of the strain vector
    must be multiplied by 2.
    
    see also AFEM.Ch09, p. 9-14
*/
void stiffnessMatrix( double64 E, double64 nu, DenseMatrix<DM_MIN>& D )
 {
    D.Resize(6,6);
    assert( nu < 0.5 );
    assert( E > 0. );
    const double64 oneMinNu      = (1. - nu);
    const double64 oneMin2NuDiv2 = (1./2. - nu);
    const double64 multiplier    =  E / ((1. + nu)*(1. - 2.*nu));
    
    D.Zero();
    D(0,0) = D(1,1) = D(2,2) = oneMinNu;
    D(3,3) = D(4,4) = D(5,5) = oneMin2NuDiv2;
    D(0,1) = D(0,2) = D(1,2) = nu;
    D(1,0) = D(2,0) = D(2,1) = nu;
    D *= multiplier;
 }



/// 1D linear bar
void stiffnessMatrix( double64 E, DenseMatrix<DM_MIN>& D, double64 length )
  {
    D.Resize(2,2);
    D.Zero(); 

    D(0,0) = D(1,1) = 1.;
    D(0,1) = D(1,0) = -1.;
    D *= E/length;
  }





/**
    Deformation only leads to in-plane stresses,
    but material gets extruded perpendicular to the plane of the model.
*/
void planeStressMatrix( const vector<ScalarVariable >& E, 
                        const vector<ScalarVariable >& nu, 
                        vector<DenseMatrix<DM_MIN> >& D )
 {
    if ( E.empty() or E.size() != nu.size() ) {
         cout <<"\nplaneStressMatrix: input data vectors are empty or do not match in size. ";
         cout <<"nothing was computed."<< endl;
         return;
      }
    D.resize(E.size());
    vector<ScalarVariable>::const_iterator  Eit(E.begin()), nit(nu.begin());
    vector<DenseMatrix<DM_MIN> >::iterator  Dit(D.begin());
    while( Eit != E.end() ) {
         planeStressMatrix( (*Eit).Value(), (*nit).Value(), (*Dit) );
         Eit++;
         nit++;
         Dit++;
      }
      
 } // end






/**
    All strain is confined to the plane of the model.
*/
void planeStrainMatrix( const vector<ScalarVariable >& E,
                        const vector<ScalarVariable >& nu, 
                        vector<DenseMatrix<DM_MIN> >& D )
 {
    if ( E.empty() or E.size() != nu.size() ) {
         cout <<"\nplaneStrainMatrix: input data vectors are empty or do not match in size. ";
         cout <<"nothing was computed."<< endl;
         return;
      }
    D.resize(E.size());
    vector<ScalarVariable >::const_iterator  Eit(E.begin()), nit(nu.begin());
    vector<DenseMatrix<DM_MIN> >::iterator   Dit(D.begin());
    while( Eit != E.end() ) {
         planeStrainMatrix( (*Eit).Value(), (*nit).Value(), (*Dit) );
         Eit++;
         nit++;
         Dit++;
      }

 } // end








/**
    Creates array of 3D stiffness matrices, one for each element integration point.
*/
void stiffnessMatrix( const vector<ScalarVariable >& E, 
                      const vector<ScalarVariable >& nu, 
                      vector<DenseMatrix<DM_MIN> >& D )
 {
    if ( E.empty() or E.size() != nu.size() ) {
         cerr <<"\nstiffnessMatrix: input data vectors are empty or do not match in size. ";
         cerr <<"nothing was computed."<< endl;
         return;
      }
    // stiffness matrix
    D.resize(E.size());
    vector<ScalarVariable >::const_iterator  Eit(E.begin()), nit(nu.begin());
    vector<DenseMatrix<DM_MIN> >::iterator   Dit(D.begin());
    while( Eit != E.end() ) {
         stiffnessMatrix( (*Eit).Value(), (*nit).Value(), (*Dit) );
         Eit++;
         nit++;
         Dit++;
      }

 } // end





  
  
  

/// in the order of the biggest positive to the smallest eigenvalue
void sortEigenVectorsAndValues( VectorVariable<2U>& vc, TensorVariable<2U>& ts )
 {
    //  if they are already in the right order nothing has to be done
    if ( vc[0U] >= vc[1U] ) return;
    
    // the eigenvalues and vectors must be flipped
    double64 swap1 = vc[0U];
    vc(0U) = vc[1U];
    vc(1U) = swap1;
    // now the vectors
    swap1           = ts(0,0);
    double64 swap2 = ts(0,1); 
    ts(0,0) = ts(1,0);
    ts(0,1) = ts(1,1);
    ts(1,0) = swap1;
    ts(1,1) = swap2;
 }



void sortEigenVectorsAndValues( VectorVariable<3U>& vc, TensorVariable<3U>& ts )
 {
    //  if they are already in the right order nothing has to be done
    if ( vc[0U] >= vc[1U] and vc[1U] >= vc[2U] ) return;
    
    // the eigenvalues and vectors are ordered using less
    multimap<double64,size_t>  eorder;
    eorder.insert( make_pair(vc[0U],0) );
    eorder.insert( make_pair(vc[1U],1) );
    eorder.insert( make_pair(vc[2U],2) );
    
    multimap<double64,size_t>::reverse_iterator  it(eorder.rbegin());
    double64 v00 = ts( (*it).second, 0 );
    double64 v01 = ts( (*it).second, 1 );
    double64 v02 = ts( (*it).second, 2 );
    vc(0) = (*it++).first;
    double64 v10 = ts( (*it).second, 0 );
    double64 v11 = ts( (*it).second, 1 );
    double64 v12 = ts( (*it).second, 2 );
    vc(1) = (*it++).first;
    double64 v20 = ts( (*it).second, 0 );
    double64 v21 = ts( (*it).second, 1 );
    double64 v22 = ts( (*it).second, 2 );
    vc(2) = (*it).first; 

    // assigning the sorted values to the tensor of eigenvectors
    ts(0,0) = v00;
    ts(0,1) = v01;
    ts(0,2) = v02;
    ts(1,0) = v10;
    ts(1,1) = v11;
    ts(1,2) = v12;
    ts(2,0) = v20;
    ts(2,1) = v21;
    ts(2,2) = v22;
 }





/**

Solves for the surface normal and tangential stresses, i.e. the
first and second term in the equation,  

 t = (t . n) n + n x (t x n)  
 
 respectively, using Cauchy's formula. For description see 
 Pollard & Fletcher, p. 215.  

@section arguments Input Arguments 

The first argument is a tensor of the Cartesian stresses (as opposed to
the principal stresses). The second argument is the outward-pointing unit
normal to the surface of interest.  

@return The computed stress components are returned into the 3rd and 4th 
function arguments.  
*/
void normalAndShearStressOnPlane( const TensorVariable<3U>& ts,
                                  const Point<3U>& un, 
                                  double64& sigma_n, double64& sigma_s )
 {
   // Cauchy's formula applied to find traction vector components, P&F, p. 213
   // VectorVariable<3U> t = ts * un.Coordinates();
   double64  tx = ts(0,0) * un[0] + ts(1,0) * un[1] + ts(2,0) * un[2]; 
   double64  ty = ts(0,1) * un[0] + ts(1,1) * un[1] + ts(2,1) * un[2]; 
   double64  tz = ts(0,2) * un[0] + ts(1,2) * un[1] + ts(2,2) * un[2]; 
 
    // (t . n) n (eqn. 6.49, P&F, p.216)
    sigma_n  = ts(0,0) * un[0] * un[0] + ts(1,1) * un[1] * un[1] + ts(2,2) * un[2] * un[2];
    sigma_n += 2. * ts(0,1) * un[0] * un[1] + 2. * ts(1,2) * un[1] * un[2] + 2. * ts(2,0) * un[2] * un[0];
    
    // n x (t x n) (eqn. 6.52, P&F, p. 216) -> vector product, vp
    double64  vpx = ((1. - un[0] * un[0]) * tx - un[0] * un[1] * ty - un[0] * un[2] * tz);  // * ex;
    double64  vpy = (-un[0] * un[1] * tx + (1. - un[1] * un[1]) * ty - un[1] * un[2] * tz); // * ey;
    double64  vpz = (-un[2] * un[0] * tx - un[2] * un[1] * ty + (1. - un[2] * un[2]) * tz); // * ez;
    
    // the shear stress is the magnitude of the vector product
    sigma_s = sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
    
 } // end normalAndShearStressOnPlane



/**
  @}
  */


} // end namespace csmp














