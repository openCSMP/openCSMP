#ifndef MECHANICS_H
#define MECHANICS_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"
#include "TensorVariable.h"

namespace csmp {

void planeStressMatrix( double E, double nu, DenseMatrix<DM3>& D );

void planeStrainMatrix( double E, double nu, DenseMatrix<DM3>& D );

void stiffnessMatrix( double E, double nu, DenseMatrix<DM6>& D );

void stiffnessMatrix( double E, DenseMatrix<DM6>& D, double length );


void planeStressMatrix( const std::vector<ScalarVariable >& E, 
                        const std::vector<ScalarVariable >& nu, 
                        std::vector<DenseMatrix<DM3> >& D );

void planeStrainMatrix( const std::vector<ScalarVariable >& E, 
                        const std::vector<ScalarVariable >& nu, 
                        std::vector<DenseMatrix<DM3> >& D );

void stiffnessMatrix( const std::vector<ScalarVariable >& E, 
                      const std::vector<ScalarVariable >& nu, 
                      std::vector<DenseMatrix<DM6> >& D );


// by the magnitude of the Eigenvalues
void sortEigenVectorsAndValues( VectorVariable<2U>& vc, TensorVariable<2U>& ts );
void sortEigenVectorsAndValues( VectorVariable<3U>& vc, TensorVariable<3U>& ts );

template<uint32_t DM_SIZE>
void convertColumnTo( const DenseMatrix<DM_SIZE>& INP, uint32_t column, TensorVariable<2U>& ts );

template<uint32_t DM_SIZE>
void convertColumnTo( const DenseMatrix<DM_SIZE>& INP, uint32_t column, TensorVariable<3U>& ts );

void convertTo( const std::vector<double>& INP, uint32_t entry, TensorVariable<2U>& ts );
void convertTo( const std::vector<double>& INP, uint32_t entry, TensorVariable<3U>& ts );

void extractRowTo( const DenseMatrix<DM3>& INP, uint32_t row, VectorVariable<2U>& vc );
void extractRowTo( const DenseMatrix<DM3>& INP, uint32_t row, VectorVariable<3U>& vc );


// various utile functions

/// uses Cauchy's formula to resolve normal and magnitude of shear stress on a surface defined by its normal
void normalAndShearStressOnPlane( const TensorVariable<3U>& cartesian_stress,
                                  const Point<3U>& plane_normal, 
                                  double& sigma_n, double& sigma_s );

/// uses Cauchy's formula to resolve normal and shear stress on a surface defined by its normal
void normalAndShearStressOnPlane( const TensorVariable<3U>& cartesian_stress,
                                  const Point<3U>& plane_normal, 
                                  double& sigma_n, VectorVariable<3U>& sigma_s );


// inline function definitions


/** Extracts row of DenseMatrix into VectorVariable.
*/
inline void extractRowTo( const DenseMatrix<DM_MIN>& INP, 
                          uint32_t row, VectorVariable<2U>& vc )
 {
    vc(0) = INP(row,0); 
    vc(1) = INP(row,1); 
 }

inline void extractRowTo( const DenseMatrix<DM_MIN>& INP, 
                          uint32_t row, VectorVariable<3U>& vc )
 {
    vc(0) = INP(row,0); 
    vc(1) = INP(row,1); 
    vc(2) = INP(row,1); 
 }


/** Extracts column of DenseMatrix into TensorVariable.
*/
template<uint32_t DM_SIZE>
inline void convertColumnTo( const DenseMatrix<DM_SIZE>& INP, 
                             uint32_t column, TensorVariable<2U>& ts )
 {
     static_assert( DM_SIZE >= 2U, "ERROR: convertColumnTo: INP matrix too small");
     // building symmetric 2D tensor
     ts(0,0) = INP(0,column);
     ts(1,1) = INP(1,column);
     ts(0,1) = ts(1,0) = INP(2,column);
 }

template<uint32_t DM_SIZE>
inline void convertColumnTo( const DenseMatrix<DM_SIZE>& INP, 
                             uint32_t column, TensorVariable<3U>& ts )
 {
     static_assert( DM_SIZE >= 3U, "ERROR: convertColumnTo: INP matrix too small");
     // building symmetric 3D tensor
     ts(0,0) = INP(0,column);
     ts(1,1) = INP(1,column);
     ts(2,2) = INP(2,column);
     ts(0,1) = ts(1,0) = INP(3,column);
     ts(0,2) = ts(2,0) = INP(5,column);
     ts(1,2) = ts(2,1) = INP(4,column);
 }



/** Extracts linearily stored tensor from vector into TensorVariable.
*/
inline void convertTo( const std::vector<double>& INP, 
                       uint32_t ip, TensorVariable<2U>& ts )
 {
    const uint32_t components(3U);
    // diagonal elements
    ts(0,0) = INP[ ip * components ];
    ts(1,1) = INP[ ip * components + 1U ];
    // offdiagonal elements of symmetric tensor
    ts(0,1) = ts(1,0) = INP[ ip * components + 2U ];
 }


inline void convertTo( const std::vector<double>& INP, 
                       uint32_t ip, TensorVariable<3U>& ts )
 {
    const uint32_t components(6U);
    // diagonal elements
    ts(0,0) = INP[ ip * components ];
    ts(1,1) = INP[ ip * components + 1U ];
    ts(2,2) = INP[ ip * components + 2U ];
    // off-diagonal elements of symmetric tensor
    ts(0,1) = ts(1,0) = INP[ ip * components + 3U ];
    ts(0,2) = ts(2,0) = INP[ ip * components + 4U ];
    ts(1,2) = ts(2,1) = INP[ ip * components + 5U ];
 }


} // end namespace csmp

#endif
