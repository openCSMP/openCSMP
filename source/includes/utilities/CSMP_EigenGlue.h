// Copyright © 2025 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_EIGEN_GLUE_H
#define CSMP_EIGEN_GLUE_H

#include <Eigen/Dense>
#include "DenseMatrix.h"
#include "Point.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

/** USAGE EXAMPLE
 *
 * @file example_csmp_eigen_usage.cpp
 * @brief Demonstrates conversions and arithmetic operations between CSMP++ types and Eigen types
 *        using the CSMP_EigenGlue.h utilities.

@code
#include <iostream>
#include <Eigen/Dense>
#include "CSMP_EigenGlue.h" 

// Assume CSMP++ namespaces and types exist:
// csmp::Point<dim>, csmp::VectorVariable<dim>, csmp::TensorVariable<dim>, csmp::DenseMatrix<dim,dim>

int main() {
    constexpr uint32_t dim = 3;

    std::cout << "=== CSMP++ and Eigen Interoperability Example ===" << std::endl;

    // -------------------------------------------------------------------------
    // 1. Create CSMP++ types
    // -------------------------------------------------------------------------
    csmp::Point<dim> csmp_point(1.0, 2.0, 3.0);
    csmp::VectorVariable<dim> csmp_vector;
    csmp_vector(0) = 4.0; csmp_vector(1) = 5.0; csmp_vector(2) = 6.0;

    csmp::TensorVariable<dim> csmp_tensor;
    csmp_tensor(0,0) = 1.0; csmp_tensor(0,1) = 2.0; csmp_tensor(0,2) = 3.0;
    csmp_tensor(1,0) = 4.0; csmp_tensor(1,1) = 5.0; csmp_tensor(1,2) = 6.0;
    csmp_tensor(2,0) = 7.0; csmp_tensor(2,1) = 8.0; csmp_tensor(2,2) = 9.0;

    std::cout << "\nCSMP++ Point: " << csmp_point << std::endl;
    std::cout << "CSMP++ Vector: " << csmp_vector << std::endl;

    // -------------------------------------------------------------------------
    // 2. Convert CSMP++ → Eigen
    // -------------------------------------------------------------------------
    Eigen::Vector3d eigen_vec = toEigen(csmp_vector);
    Eigen::Vector3d eigen_point = toEigen(csmp_point);
    Eigen::Matrix3d eigen_tensor = toEigen(csmp_tensor);

    std::cout << "\nConverted to Eigen:\n";
    std::cout << "Eigen Vector: " << eigen_vec.transpose() << std::endl;
    std::cout << "Eigen Point: " << eigen_point.transpose() << std::endl;
    std::cout << "Eigen Tensor:\n" << eigen_tensor << std::endl;

    // -------------------------------------------------------------------------
    // 3. Convert Eigen → CSMP++
    // -------------------------------------------------------------------------
    Eigen::Vector3d eigen_new(10.0, 11.0, 12.0);
    csmp::VectorVariable<dim> csmp_from_eigen = toCSMP<dim>(eigen_new);

    std::cout << "\nConverted Eigen → CSMP Vector: " << csmp_from_eigen << std::endl;

    // -------------------------------------------------------------------------
    // 4. Perform Mixed Arithmetic Operations
    // -------------------------------------------------------------------------
    // Vector addition (CSMP + Eigen)
    auto csmp_plus_eigen = csmp_vector + eigen_vec; // Returns Eigen::VectorXd
    std::cout << "\nCSMP Vector + Eigen Vector = " << csmp_plus_eigen.transpose() << std::endl;

    // Eigen + CSMP
    auto eigen_plus_csmp = eigen_vec + csmp_vector;
    std::cout << "Eigen Vector + CSMP Vector = " << eigen_plus_csmp.transpose() << std::endl;

    // Tensor-vector multiplication (CSMP Tensor * Eigen Vector)
    auto result_tensor_vec = csmp_tensor * eigen_vec;
    std::cout << "\nCSMP Tensor * Eigen Vector = " << result_tensor_vec.transpose() << std::endl;

    // Tensor-tensor addition (CSMP Tensor + Eigen Matrix)
    auto tensor_sum = csmp_tensor + eigen_tensor;
    std::cout << "CSMP Tensor + Eigen Matrix:\n" << tensor_sum << std::endl;

    // -------------------------------------------------------------------------
    // 5. Cross Product Example
    // -------------------------------------------------------------------------
    Eigen::Vector3d eigen_cross = eigen_vec.cross(toEigen(csmp_vector));
    std::cout << "\nEigen Cross Product: " << eigen_cross.transpose() << std::endl;

    return 0;
}
 @endcode
 *
*/

namespace csmp {
namespace eigen_glue {

// =============================================================
// Conversion: csmp::Point<dim> ↔ Eigen::Matrix<double, dim, 1>
// =============================================================

template<uint32_t dim>
inline Eigen::Matrix<double, dim, 1> ToEigen(const csmp::Point<dim>& p) {
    Eigen::Matrix<double, dim, 1> v;
    for (uint32_t i{0u}; i < dim; ++i) v[i] = p[i];
    return v;
}

template<uint32_t dim>
inline csmp::Point<dim> FromEigen(const Eigen::Matrix<double, dim, 1>& v) {
    csmp::Point<dim> p;
    for (uint32_t i{0u}; i < dim; ++i) p[i] = v[i];
    return p;
}

// =============================================================
// Conversion: csmp::VectorVariable<dim> ↔ Eigen::Matrix<double, dim, 1>
// =============================================================

template<uint32_t dim>
inline Eigen::Matrix<double, dim, 1> ToEigen(const csmp::VectorVariable<dim>& vv) {
    Eigen::Matrix<double, dim, 1> v;
    for (uint32_t i{0u}; i < dim; ++i) v[i] = vv[i];
    return v;
}

template<uint32_t dim>
inline csmp::VectorVariable<dim> FromEigenToVectorVar(const Eigen::Matrix<double, dim, 1>& v) {
    csmp::VectorVariable<dim> vv;
    for (uint32_t i{0u}; i < dim; ++i) vv[i] = v[i];
    return vv;
}

// =============================================================
// Conversion: csmp::TensorVariable<dim> ↔ Eigen::Matrix<double, dim, dim>
// =============================================================

template<uint32_t dim>
inline Eigen::Matrix<double, dim, dim> ToEigen(const csmp::TensorVariable<dim>& ts ) {
    Eigen::Matrix<double, dim, dim> m;
    for (uint32_t i{0u}; i < dim; ++i)
        for (uint32_t j{0u}; j < dim; ++j) m(i, j) = ts(i, j);
    return m;
}

template<uint32_t dim>
inline csmp::TensorVariable<dim> FromEigenToTensorVar(const Eigen::Matrix<double, dim, dim>& m) {
    csmp::TensorVariable<dim> tv;
    for (uint32_t i{0u}; i < dim; ++i)
        for (uint32_t j{0u}; j < dim; ++j) tv(i, j) = m(i, j);
    return tv;
}

// =============================================================
// Conversion: csmp::DenseMatrix<Rows, Cols> ↔ Eigen::Matrix<double, Rows, Cols>
// =============================================================

template<uint32_t Rows, uint32_t Cols>
inline Eigen::Matrix<double, Rows, Cols> ToEigen( const csmp::DenseMatrix<DM_MIN>& dm ) {
    Eigen::Matrix<double, Rows, Cols> m;
    for (uint32_t i{0u}; i < Rows; ++i)
        for (uint32_t j{0u}; j < Cols; ++j) m(i, j) = dm(i, j);
    return m;
}

template<uint32_t Rows, uint32_t Cols>
inline csmp::DenseMatrix<DM_MIN> FromEigenToDenseMatrix( const Eigen::Matrix<double,Rows,Cols>& m ) {
    csmp::DenseMatrix<DM_MIN> dm;
    dm.Resize(Rows,Cols);
    for (uint32_t i{0u}; i < Rows; ++i)
        for (uint32_t j{0u}; j < Cols; ++j) dm(i, j) = m(i, j);
    return dm;
}


// ========================
// Arithmetic Operator Overloads
// ========================

// VectorVariable<dim> + Eigen::Vector
template <uint32_t dim>
inline Eigen::Matrix<double, dim, 1> operator+(const csmp::VectorVariable<dim>& lhs, const Eigen::Matrix<double, dim, 1>& rhs) {
    return ToEigen(lhs) + rhs;
}

template <uint32_t dim>
inline Eigen::Matrix<double, dim, 1> operator+(const Eigen::Matrix<double, dim, 1>& lhs, const csmp::VectorVariable<dim>& rhs) {
    return lhs + ToEigen(rhs);
}

// VectorVariable<dim> - Eigen::Vector
template <uint32_t dim>
inline Eigen::Matrix<double, dim, 1> operator-(const csmp::VectorVariable<dim>& lhs, const Eigen::Matrix<double, dim, 1>& rhs) {
    return ToEigen(lhs) - rhs;
}

// Scalar * VectorVariable<dim>
template <uint32_t dim>
inline Eigen::Matrix<double, dim, 1> operator*(double scalar, const csmp::VectorVariable<dim>& v) {
    return scalar * ToEigen(v);
}

// TensorVariable<dim> * Eigen::Vector
template <uint32_t dim>
inline Eigen::Matrix<double, dim, 1> operator*(const csmp::TensorVariable<dim>& T, const Eigen::Matrix<double, dim, 1>& v) {
    return ToEigen(T) * v;
}

// Eigen::Matrix * TensorVariable<dim>
template <uint32_t dim>
inline Eigen::Matrix<double, dim, dim> operator*(const Eigen::Matrix<double, dim, dim>& M, const csmp::TensorVariable<dim>& T) {
    return M * ToEigen(T);
}

// TensorVariable<dim> + Eigen::Matrix
template <uint32_t dim>
inline Eigen::Matrix<double, dim, dim> operator+(const csmp::TensorVariable<dim>& T, const Eigen::Matrix<double, dim, dim>& M) {
    return ToEigen(T) + M;
}



} // namespace eigen_glue
} // namespace csmp

#endif // CSMP_EIGEN_GLUE_H
