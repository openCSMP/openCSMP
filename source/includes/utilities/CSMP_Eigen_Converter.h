// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#pragma once

#include <Eigen/Dense>
#include <type_traits>
#include <cstdint>

/**
   Template that converts csmp variables into Eigen::Dense fixed-size matrices.
   The model dimension plays a role: a ScalarVariable becomes a dim x dim Eigen::DenseMatrix dim by dim
   with the scalar values placed into the diagonal of the zero'ed matrix.
   
   A VectorVariable<dim> goes into the diagonal of the dense matrix.
   A TensorVariable<dim> fully populates the denseMatrix.
 */

namespace csmp {

// ------------------------------------------------------------------
// 1. Forward Declarations
// (Omit or replace these if you include the actual CSMP variable headers)
// ------------------------------------------------------------------
class ScalarVariable;
template<uint32_t dim> class VectorVariable;
template<uint32_t dim> class TensorVariable;

// ------------------------------------------------------------------
// 2. Type Traits for Compile-Time Identification
// ------------------------------------------------------------------
template<typename T> struct is_scalar_variable : std::false_type {};
template<> struct is_scalar_variable<ScalarVariable> : std::true_type {};

template<typename T> struct is_vector_variable : std::false_type {};
template<uint32_t dim> struct is_vector_variable<VectorVariable<dim>> : std::true_type {};

template<typename T> struct is_tensor_variable : std::false_type {};
template<uint32_t dim> struct is_tensor_variable<TensorVariable<dim>> : std::true_type {};

// ------------------------------------------------------------------
// 3. The Converter Template Function
// ------------------------------------------------------------------
/**
 * Converts a CSMP variable (Scalar, Vector, or Tensor) into a fixed-size 
 * dim x dim Eigen::Matrix.
 */
template<uint32_t dim, typename VarType>
Eigen::Matrix<double, dim, dim> ToEigenMatrix(const VarType& var)
{
    using MatrixType = Eigen::Matrix<double, dim, dim>;

    if constexpr (is_scalar_variable<VarType>::value)
    {
        // Scalar -> Diagonal matrix with the scalar value on the main diagonal
        // Note: Adjust the extraction method if ScalarVariable doesn't support cast-to-double
        const double val = static_cast<double>(var); 
        return MatrixType::Identity() * val;
    }
    else if constexpr (is_vector_variable<VarType>::value)
    {
        // Vector -> Diagonal matrix with vector components on the main diagonal
        MatrixType mat = MatrixType::Zero();
        for (uint32_t i{0U}; i < dim; ++i)
        {
            // Note: Adjust to your specific getter if it's not operator[]
            mat(i, i) = var[i]; 
        }
        return mat;
    }
    else if constexpr (is_tensor_variable<VarType>::value)
    {
        // Tensor -> Fully populated dim x dim matrix
        MatrixType mat;
        for (uint32_t i{0U}; i < dim; ++i)
        {
            for (uint32_t j{0U}; j < dim; ++j)
            {
                // Note: Adjust to your specific getter if it's not operator()(i,j)
                mat(i, j) = var(i, j);
            }
        }
        return mat;
    }
    else
    {
        // In C++23, a bare static_assert(false) in an uninstantiated constexpr branch is fully legal
        static_assert(false, "Unsupported CSMP variable type provided to ToEigenMatrix.");
    }
}

template<uint32_t mn_max> class DenseMatrix;

/**
 * Converts a csmp::DenseMatrix<mn_max> into a fixed-size dim x dim Eigen::Matrix.
 * Uses Eigen::Map for zero-overhead direct initialization, accounting for memory strides.
 * * @tparam dim     The exact spatial dimension needed for the Eigen matrix (e.g., 2 or 3).
 * @tparam mn_max  The maximum capacity of the DenseMatrix (auto-deduced by the compiler).
 */
template<uint32_t dim, uint32_t mn_max>
Eigen::Matrix<double, dim, dim> ToEigenMatrix(const DenseMatrix<mn_max>& mat)
{
    // 1. Compile-time and runtime dimension checks
    static_assert(dim <= mn_max, "Requested Eigen matrix dimension exceeds DenseMatrix static capacity (mn_max).");
    assert(mat.Rows() == dim && "DenseMatrix active rows do not match the requested target dimension.");
    assert(mat.Cols() == dim && "DenseMatrix active cols do not match the requested target dimension.");

    // 2. Get a pointer to the contiguous start of the 2D std::array
    const double* raw_ptr = &mat.Data()[0][0];

    // 3. Define an Eigen::Map that understands the memory layout.
    // - std::array is Row-Major.
    // - The distance in memory from one row to the next is exactly mn_max.
    using Stride     = Eigen::OuterStride<mn_max>;
    using MatrixView = Eigen::Map<const Eigen::Matrix<double, dim, dim, Eigen::RowMajor>, 0, Stride>;

    // 4. Construct and return. 
    // Eigen evaluates the map and directly constructs the new matrix. 
    // Note: The return type defaults to Column-Major, so Eigen will automatically 
    // and optimally transpose the layout in memory during this assignment.
    return MatrixView(raw_ptr);
}

} // end namespace csmp
